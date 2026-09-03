#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'ShortRead', 'parallel', 'data.table', 'BiocParallel', 'stringi')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",               type = "character",     required = TRUE,          help = "Directory for output files")
parser$add_argument("--inputData",               type = "character",     required = TRUE,          help = "Path to prepReads module's rds output file.")
parser$add_argument("--softwareRoot",            type = "character",     required = TRUE,          help = "Path to INSPIIRED2 installation.")
parser$add_argument("--threads",                 type = "integer",       default = 50,             help = "Number of threads to use.")
parser$add_argument("--fileTag",                 type = "character",     default = "alignReads",   help = "String appended to output files in the outpt directory.")
parser$add_argument("--ramDiskPath",             type = "character",     default = "/dev/shm",     help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
parser$add_argument("--minPercentID",            type = "double",        default = 95,             help = "Min BLAT alignment percent ID score (0 .. 100).")
parser$add_argument("--minAlignmentCoverage",    type = "double",        default = 95,             help = "Min BLAT alignment percent query coverage (0 .. 100).")
parser$add_argument("--blatStepSize",            type = "integer",       default = 5,              help = "BLAT step size.")
parser$add_argument("--blatTileSize",            type = "integer",       default = 11,             help = "BLAT tile size.")
parser$add_argument("--blatRepMatch",            type = "integer",       default = 3000,           help = "BLAT repMatch value.")
parser$add_argument("--blatMaxtNumInsert",       type = "integer",       default = 1,              help = "BLAT max number of target inserts.")
parser$add_argument("--blatMaxqNumInsert",       type = "integer",       default = 1,              help = "BLAT max number of query inserts.")
parser$add_argument("--blatMaxtBaseInsert",      type = "integer",       default = 1,              help = "BLAT max number of target insert NTs.")
parser$add_argument("--blatMaxqBaseInsert",      type = "integer",       default = 1,              help = "BLAT max number of query insert NTs.")
parser$add_argument("--dataRowChunkSize",        type = "integer",       default = 2500,           help = "Numbers of data rows to process per alignment worker.")

runModule <- function(){
  doneFile <- file.path(args$outputDir, paste0(args$fileTag, '.done'))
  if(file.exists(doneFile) && unlink(doneFile) != 0) stop('Error - could not remove stale completion marker: ', doneFile)
  
  if(is.na(args$threads) || args$threads < 1) stop('Error - threads must be at least 1.')
  if(is.na(args$dataRowChunkSize) || args$dataRowChunkSize < 1) stop('Error - dataRowChunkSize must be at least 1.')
  if(is.na(args$blatStepSize) || args$blatStepSize < 1) stop('Error - blatStepSize must be at least 1.')
  if(is.na(args$blatTileSize) || args$blatTileSize < 1) stop('Error - blatTileSize must be at least 1.')
  if(is.na(args$blatRepMatch) || args$blatRepMatch < 1) stop('Error - blatRepMatch must be at least 1.')
  if(any(is.na(c(args$blatMaxtNumInsert, args$blatMaxqNumInsert, args$blatMaxtBaseInsert, args$blatMaxqBaseInsert))) ||
     any(c(args$blatMaxtNumInsert, args$blatMaxqNumInsert, args$blatMaxtBaseInsert, args$blatMaxqBaseInsert) < 0)){
    stop('Error - BLAT insertion thresholds must be non-negative.')
  }
  if(is.na(args$minPercentID) || !is.finite(args$minPercentID) || args$minPercentID < 0 || args$minPercentID > 100){
    stop('Error - minPercentID must be between 0 and 100.')
  }
  if(is.na(args$minAlignmentCoverage) || !is.finite(args$minAlignmentCoverage) || args$minAlignmentCoverage < 0 || args$minAlignmentCoverage > 100){
    stop('Error - minAlignmentCoverage must be between 0 and 100.')
  }
  
  startModule()
  
  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  moduleSucceeded <- FALSE
  
  on.exit({
    unlink(args$tmpDir,  recursive = TRUE, force = TRUE)
    if(moduleSucceeded) unlink(args$logDir, recursive = TRUE, force = TRUE)
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  updateLog('Starting alignReads module.')
  
  resource_overlay()
  
  if(! file.exists(args$inputData))  stop(paste0('Error - the input data file (', args$inputData, ') does not exist.'))
  if(file.size(args$inputData) == 0) stop(paste0('Error - the input data file (', args$inputData, ') is empty.'))
  
  blatExecutable <- Sys.which('blat')
  if(!nzchar(blatExecutable)) stop('Error - could not find the BLAT executable in PATH.')
  
  pslScoreExecutable <- file.path(args$softwareRoot, 'bin', 'pslScore.pl')
  if(!file.exists(pslScoreExecutable) || file.access(pslScoreExecutable, mode = 1) != 0){
    stop('Error - pslScore.pl does not exist or is not executable: ', pslScoreExecutable)
  }
  
  inputData <- readRDS(args$inputData)
  if(!is.data.frame(inputData)) stop('Error - inputData must contain a data.frame or data.table.')
  if(nrow(inputData) == 0) stop('Error - inputData contains no reads.')
  
  requiredColumns <- c('trial', 'subject', 'sample', 'replicate', 'UMI', 'anchorReadSeq',
                       'adriftReadSeq', 'nReads', 'readID', 'mode', 'refGenome',
                       'vectorFastaFile', 'leaderSeqHMM', 'leaderSeq')
  missingColumns <- setdiff(requiredColumns, names(inputData))
  if(length(missingColumns) > 0){
    stop('Error - inputData is missing required column(s): ', paste(missingColumns, collapse = ', '))
  }
  
  if(anyNA(inputData$readID) || any(!nzchar(as.character(inputData$readID)))){
    stop('Error - readID values must be non-missing and non-empty.')
  }
  if(anyDuplicated(inputData$readID)) stop('Error - readID values must be unique.')
  
  refGenomeIsFactor <- is.factor(inputData$refGenome)
  refGenomeLevels <- if(refGenomeIsFactor) levels(inputData$refGenome) else NULL
  inputData$refGenome <- as.character(inputData$refGenome)
  if(anyNA(inputData$refGenome) || any(!nzchar(inputData$refGenome))){
    stop('Error - refGenome values must be non-missing and non-empty.')
  }
  if(any(basename(inputData$refGenome) != inputData$refGenome)){
    stop('Error - refGenome values must be reference names, not paths.')
  }
  
  for(sequenceColumn in c('anchorReadSeq', 'adriftReadSeq')){
    sequenceValues <- as.character(inputData[[sequenceColumn]])
    if(anyNA(sequenceValues) || any(!nzchar(sequenceValues)) || any(grepl('[\\r\\n>]', sequenceValues))){
      stop('Error - ', sequenceColumn, ' values must be non-missing, non-empty FASTA sequences without header or newline characters.')
    }
  }
  
  refGenomeFiles <- file.path(args$softwareRoot, 'data', 'referenceGenomes', paste0(unique(inputData$refGenome), '.2bit'))
  missingRefGenomeFiles <- refGenomeFiles[!file.exists(refGenomeFiles)]
  if(length(missingRefGenomeFiles) > 0){
    stop('Error - reference genome file(s) do not exist: ', paste(missingRefGenomeFiles, collapse = ', '))
  }
  
  anchorReads <- inputData %>% dplyr::select(-adriftReadSeq) %>% dplyr::rename(seq = anchorReadSeq)
  adriftReads <- inputData %>% dplyr::select(-anchorReadSeq, -leaderSeq) %>% dplyr::rename(seq = adriftReadSeq)
  rm(inputData)
  
  add_sequence_ids <- function(reads){
    sequenceKeys <- dplyr::distinct(reads, refGenome, seq)
    sequenceKeys$seqNum <- paste0('s', seq_len(nrow(sequenceKeys)))
    dplyr::left_join(reads, sequenceKeys, by = c('refGenome', 'seq'), relationship = 'many-to-one')
  }
  
  anchorReads <- add_sequence_ids(anchorReads)
  adriftReads <- add_sequence_ids(adriftReads)
  
  empty_alignment_table <- function(){
    data.table(qName = character(), refGenome = character(), tName = character(),
               strand = character(), tStart = integer(), tEnd = integer())
  }
  
  # Yield chunks from one reference genome at a time while allowing chunks from
  # different genomes to be processed concurrently by MulticoreParam.
  make_reference_iterator <- function(dt, chunk_size, chunk_num_start = 0L){
    dt <- dplyr::arrange(dt, refGenome, seqNum)
    current_row <- 1L
    total_rows <- nrow(dt)
    chunk_num <- as.integer(chunk_num_start)
    reference_end_rows <- if(total_rows > 0) cumsum(rle(as.character(dt$refGenome))$lengths) else integer()
    reference_num <- 1L
    
    function(){
      if(current_row > total_rows) return(NULL)
      while(current_row > reference_end_rows[reference_num]) reference_num <<- reference_num + 1L
      
      end_row <- min(current_row + chunk_size - 1L, reference_end_rows[reference_num])
      chunk_data <- dt[current_row:end_row, ]
      chunk_num <<- chunk_num + 1L
      current_row <<- end_row + 1L
      
      list(data = chunk_data, chunk_num = chunk_num, is_last = (current_row > total_rows))
    }
  }
  
  alignment_worker <- function(chunk, ...) {
    chunkLogDir <- file.path(args$logDir, paste0('chunk_', chunk$chunk_num))
    if(!dir.exists(chunkLogDir)) dir.create(chunkLogDir, recursive = TRUE)
    if(!dir.exists(chunkLogDir)) stop('Error - could not create worker log directory: ', chunkLogDir)
    logFile <- file.path(chunkLogDir, 'log')
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tStarting alignment chunk with ', ppNum(nrow(chunk$data)), ' data rows.'), logFile = logFile)
    
    chunkRefGenome <- unique(as.character(chunk$data$refGenome))
    if(length(chunkRefGenome) != 1 || is.na(chunkRefGenome) || !nzchar(chunkRefGenome)){
      stop('Error - alignment chunk ', chunk$chunk_num, ' does not contain exactly one reference genome.')
    }
    
    ts <- tmpString()
    fastaFile <- file.path(args$ramDisk, paste0(ts, '.fasta'))
    pslFile <- file.path(args$ramDisk, paste0(ts, '.psl'))
    on.exit(unlink(c(fastaFile, pslFile), force = TRUE), add = TRUE)
    
    writeLines(paste0('>', chunk$data$seqNum, '\n', chunk$data$seq), con = fastaFile)
    
    referenceFile <- file.path(args$softwareRoot, 'data', 'referenceGenomes', paste0(chunkRefGenome, '.2bit'))
    
    blatArgs <- c(shQuote(referenceFile), shQuote(fastaFile), shQuote(pslFile),
                  paste0('-tileSize=', args$blatTileSize),
                  paste0('-stepSize=', args$blatStepSize),
                  paste0('-repMatch=', args$blatRepMatch),
                  '-out=psl', '-t=dna', '-q=dna', '-minScore=0', '-minIdentity=0', '-noHead', '-noTrimA')
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\t', shQuote(blatExecutable), ' ', paste(blatArgs, collapse = ' ')), logFile = logFile)
    
    blatOutput <- tryCatch(
      suppressWarnings(system2(blatExecutable, args = blatArgs, stdout = TRUE, stderr = TRUE)),
      error = function(e) stop('Error - BLAT could not be executed for chunk ', chunk$chunk_num, ': ', conditionMessage(e))
    )
    blatStatus <- attr(blatOutput, 'status')
    if(is.null(blatStatus)) blatStatus <- 0L
    if(length(blatOutput) > 0){
      updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tBLAT: ', paste(blatOutput, collapse = '\n')), logFile = logFile)
    }
    if(blatStatus != 0L){
      failureMessage <- paste0('Error - BLAT failed for chunk ', chunk$chunk_num, ' (reference ', chunkRefGenome,
                               ') with exit status ', blatStatus, '. See ', logFile, '.')
      updateLog(paste0('<data chunk #', chunk$chunk_num, '>\t', failureMessage), logFile = logFile)
      stop(failureMessage)
    }
    if(!file.exists(pslFile)) stop('Error - BLAT did not create an output PSL file for chunk ', chunk$chunk_num, '.')
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tParsing BLAT output.'), logFile = logFile)
    
    b <- parseBLAToutput(pslFile)
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\t', ppNum(nrow(b)), ' data rows returned.'), logFile = logFile)
    
    if(nrow(b) > 0){
      if(any(!b$qName %in% chunk$data$seqNum)) stop('Error - BLAT returned an unexpected query ID for chunk ', chunk$chunk_num, '.')
      b$qPercentCoverage <- (b$qWidth / b$qSize)*100
      b <- dplyr::filter(b, queryPercentID >= args$minPercentID, 
                         qPercentCoverage >= args$minAlignmentCoverage,
                         tNumInsert <= args$blatMaxtNumInsert,
                         qNumInsert <= args$blatMaxqNumInsert,
                         tBaseInsert <= args$blatMaxtBaseInsert,
                         qBaseInsert <= args$blatMaxqBaseInsert)
      
      updateLog(paste0('<data chunk #', chunk$chunk_num, '>\t', ppNum(nrow(b)), ' data rows remain after filtering.'), logFile = logFile)
      if(nrow(b) > 0){
        b$refGenome <- chunkRefGenome
        return(data.table::as.data.table(dplyr::select(b, qName, refGenome, tName, strand, tStart, tEnd)))
      }
    }
    
    empty_alignment_table()
  }
  
  run_alignments <- function(reads, chunk_num_start){
    alignmentInputs <- dplyr::distinct(reads, refGenome, seqNum, seq)
    referenceCounts <- table(as.character(alignmentInputs$refGenome))
    nChunks <- sum(ceiling(as.numeric(referenceCounts) / args$dataRowChunkSize))
    my_iter <- make_reference_iterator(alignmentInputs, chunk_size = args$dataRowChunkSize, chunk_num_start = chunk_num_start)
    
    # Do not create more forked workers than there are chunks to process.
    param <- MulticoreParam(workers = min(args$threads, nChunks), stop.on.error = TRUE)
    results <- tryCatch(
      bpiterate(ITER = my_iter, FUN = alignment_worker, BPPARAM = param),
      finally = {
        try(bpstop(param), silent = TRUE)
        closeAllConnections()
      }
    )
    
    list(alignments = rbindlist(results, use.names = TRUE), nChunks = nChunks)
  }
  
  updateLog(paste0('Starting asynchronous BLAT runs. Data chunk logs can be found in ', args$logDir, '/'))
  updateLog('Collated data chunk logs will appear below when done.')
  
  chunk_start_num <- 0
  
  f <- suppressWarnings(as.integer(str_extract(list.files(args$logDir, full.names = FALSE), '\\d+$')))
  if(any(!is.na(f))) chunk_start_num <- max(f, na.rm = TRUE)
  
  o <- list()
  
  updateLog('Starting anchor read alignments.')
  
  alignmentResult <- run_alignments(anchorReads, chunk_start_num)
  anchorReadsAlignments <- alignmentResult$alignments
  chunk_start_num <- chunk_start_num + alignmentResult$nChunks
  rm(alignmentResult)
  
  updateLog('Anchor read alignments completed.')
  updateLog(paste0(ppNum(n_distinct(anchorReadsAlignments$qName)), ' anchor reads returned one or more alignments.'))
  
  o$anchorReads <- left_join(anchorReads, anchorReadsAlignments,
                             by = c('seqNum' = 'qName', 'refGenome' = 'refGenome'),
                             relationship = "many-to-many")
  
  o$anchorReads <- o$anchorReads[!is.na(o$anchorReads$tName), ]
  
  if(nrow(o$anchorReads) == 0) stop('Error - no anchor reads aligned to the reference.')
  
  rm(anchorReads, anchorReadsAlignments)
  gc()
  
  updateLog('Limiting adrift reads to those with anchor read alignments.')
  
  # Limit adrift reads to those with anchor read alignments.
  adriftReads <- adriftReads[adriftReads$readID %in% o$anchorReads$readID, ]
  
  if(nrow(adriftReads) == 0) stop('Error -- no adrift reads remain after limiting reads to those with anchor read mates that aligned to the reference.')
  
  # Align unique adrift read sequences.
  
  updateLog('Starting adrift read alignments.')
  
  alignmentResult <- run_alignments(adriftReads, chunk_start_num)
  adriftReadsAlignments <- alignmentResult$alignments
  rm(alignmentResult)
  
  updateLog('Adrift read alignments completed.')
  updateLog(paste0(ppNum(n_distinct(adriftReadsAlignments$qName)), ' adrift reads returned one or more alignments.'))
  
  o$adriftReads <- left_join(adriftReads, adriftReadsAlignments,
                             by = c('seqNum' = 'qName', 'refGenome' = 'refGenome'),
                             relationship = "many-to-many")
  
  o$adriftReads <- o$adriftReads[!is.na(o$adriftReads$tName), ]
  
  if(nrow(o$adriftReads) == 0) stop('Error - no adrift reads aligned to the reference.')
  
  rm(adriftReads, adriftReadsAlignments)
  gc()
  
  collated_logs <- unlist(lapply(list.files(args$logDir, full.names = TRUE, recursive = TRUE, pattern = '^log$'), readLines))
  unlink(list.files(args$logDir, full.names = TRUE), recursive = TRUE)
  updateLog(collated_logs)
  
  o$anchorReads$seqNum <- NULL
  o$adriftReads$seqNum <- NULL
  
  o$anchorReads$trial     <- as.factor(o$anchorReads$trial);      o$adriftReads$trial     <- as.factor(o$adriftReads$trial)
  o$anchorReads$subject   <- as.factor(o$anchorReads$subject);    o$adriftReads$subject   <- as.factor(o$adriftReads$subject)
  o$anchorReads$sample    <- as.factor(o$anchorReads$sample);     o$adriftReads$sample    <- as.factor(o$adriftReads$sample)
  o$anchorReads$replicate <- as.factor(o$anchorReads$replicate);  o$adriftReads$replicate <- as.factor(o$adriftReads$replicate)
  if(refGenomeIsFactor){
    o$anchorReads$refGenome <- factor(o$anchorReads$refGenome, levels = refGenomeLevels)
    o$adriftReads$refGenome <- factor(o$adriftReads$refGenome, levels = refGenomeLevels)
  }
  
  outputFile <- file.path(args$outputDir, paste0(args$fileTag, '.rds'))
  tmpOutputFile <- file.path(args$outputDir, paste0('.', args$fileTag, '.rds.', Sys.getpid(), '.tmp'))
  on.exit(if(file.exists(tmpOutputFile)) unlink(tmpOutputFile), add = TRUE)
  saveRDS(o, tmpOutputFile)
  if(!file.rename(tmpOutputFile, outputFile)) stop('Error - could not atomically replace output file: ', outputFile)
  updateLog('alignReads module completed.')
  write(date(), doneFile)
  moduleSucceeded <- TRUE
}

args <- parser$parse_args()
source(file.path(args$softwareRoot, 'lib', 'common.R'))

tryCatch({
  runModule()
}, error = function(e) {
  cat("ERROR: ", conditionMessage(e), "\n", sep = "", file = stderr())
  flush(stderr())
  quit(save = "no", status = 1, runLast = FALSE)
})

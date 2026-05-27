#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'ShortRead', 'parallel', 'data.table', 'BiocParallel', 'stringi')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",               type = "character",     required = TRUE,          help = "Directory for output files")
parser$add_argument("--inputData",               type = "character",     required = TRUE,          help = "Path to demultiplex module's rds output file.")
parser$add_argument("--softwareRoot",            type = "character",     required = TRUE,          help = "Path to INSPIIRED2 installation.")
parser$add_argument("--threads",                 type = "integer",       default = 50,             help = "Number of threads to use.")
parser$add_argument("--fileTag",                 type = "character",     default = "alignReads",   help = "String appended to output files in the outpt directory.")
parser$add_argument("--ramDiskPath",             type = "character",     default = "/dev/shm",     help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
parser$add_argument("--minPercentID",            type = "double",        default = 95,             help = "Min BLAT alignment percent ID score (0 .. 100).")
parser$add_argument("--minAlignmentCoverage",    type = "double",        default = 95,             help = "Min BLAT alignment percent query coverage (0 .. 100).")
parser$add_argument("--blatStepSize",            type = "integer",       default = 5,              help = "BLAT step size.")
parser$add_argument("--blatTileSize",            type = "integer",       default = 11,             help = "BLAT tile size.")
parser$add_argument("--blatRepMatch",            type = "integer",       default = 3000,           help = "BLAT repMatch value.")
parser$add_argument("--blatMaxtNumInsert",       type = "integer",       default = 2,              help = "BLAT max number of target inserts.")
parser$add_argument("--blatMaxqNumInsert",       type = "integer",       default = 2,              help = "BLAT max number of query inserts.")
parser$add_argument("--blatMaxtBaseInsert",      type = "integer",       default = 3,              help = "BLAT max number of target insert NTs.")
parser$add_argument("--blatMaxqBaseInsert",      type = "integer",       default = 3,              help = "BLAT max number of target insert NTs.")
parser$add_argument("--dataRowChunkSize",        type = "integer",       default = 1000,           help = "Numbers of data rows to process per alignment worker.")

runModule <- function(){
  startModule()
  
  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  on.exit({
    unlink(args$tmpDir,  recursive = TRUE, force = TRUE)
    unlink(args$logDir,  recursive = TRUE, force = TRUE) 
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  updateLog('Starting alignReads module.')
  
  if(! file.exists(args$inputData))  stop(paste0('Error - the input data file (', args$inputData, ') does not exist.'))
  if(file.size(args$inputData) == 0) stop(paste0('Error - the input data file (', args$inputData, ') is empty.'))
  
  anchorReads <- readRDS(args$inputData) %>% dplyr::select(-adriftReadSeq) %>% dplyr::rename(seq = anchorReadSeq)
  adriftReads <- readRDS(args$inputData) %>% dplyr::select(-anchorReadSeq, -leaderSeq) %>% dplyr::rename(seq = adriftReadSeq)
  
  anchorReads$seqNum <- paste0('s', as.integer(as.factor(anchorReads$seq)))
  adriftReads$seqNum <- paste0('s', as.integer(as.factor(adriftReads$seq)))
  
  alignment_worker <- function(chunk, ...) {
    if(! dir.exists(file.path(args$logDir, paste0('chunk_', chunk$chunk_num)))) dir.create(file.path(args$logDir, paste0('chunk_', chunk$chunk_num)))
    logFile <- file.path(args$logDir, paste0('chunk_', chunk$chunk_num), 'log')
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tStarting alignment chunk with ', ppNum(nrow(chunk$data)), ' data rows.'), logFile = logFile)
  
    ts <- tmpString()
    
    write(paste0('>', chunk$data$seqNum, '\n', chunk$data$seq), file = file.path(args$ramDisk, paste0(ts, '.fasta')))
    
    cmd <- paste0('blat ', file.path(args$softwareRoot, 'data', 'referenceGenomes', paste0(as.character(chunk$data$refGenome[1]), '.2bit')), ' ', 
                  file.path(args$ramDisk, paste0(ts, '.fasta')), ' ', 
                  file.path(args$ramDisk, paste0(ts, '.psl')), 
                  ' -tileSize=', args$blatTileSize, 
                  ' -stepSize=', args$blatStepSize, 
                  ' -repMatch=', args$blatRepMatch, 
                  ' -out=psl -t=dna -q=dna -minScore=0 -minIdentity=0 -noHead -noTrimA')
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\t', cmd), logFile = logFile)
    
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tParsing BLAT output.'), logFile = logFile)
    
    b <- parseBLAToutput(file.path(args$ramDisk, paste0(ts, '.psl')))
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\t', ppNum(nrow(b)), ' data rows returned.'), logFile = logFile)
    
    invisible(file.remove(list.files(args$ramDisk, pattern = ts, full.names = TRUE)))

    if(nrow(b) > 0){
      b$qPercentCoverage <- ((b$qEnd - b$qStart) / b$qSize)*100
      b <- dplyr::filter(b, queryPercentID >= args$minPercentID, 
                              qPercentCoverage >= args$minAlignmentCoverage,
                              tNumInsert <= args$blatMaxtNumInsert,
                              qNumInsert <= args$blatMaxqNumInsert,
                              tBaseInsert <= args$blatMaxtBaseInsert,
                              qBaseInsert <= args$blatMaxqBaseInsert)
      
      updateLog(paste0('<data chunk #', chunk$chunk_num, '>\t', ppNum(nrow(b)), ' data rows remain after filtering.'), logFile = logFile)
      return(b)
    } else {
      return(data.table())
    }
  }
  
  updateLog(paste0('Starting asynchronous BLAT runs. Data chunk logs can be found in ', args$logDir, '/'))
  updateLog('Collated data chunk logs will appear below when done.')

  chunk_start_num <- 0
  
  f <- list.files(args$logDir, full.names = FALSE)
  if(length(f) > 0) chunk_start_num <- max(as.integer(str_extract(f, '\\d+$')))
  
  o <- list()
  
  updateLog('Starting anchor read alignments.')

  my_iter <- make_dt_iterator(anchorReads[! duplicated(anchorReads$seqNum)], chunk_size = args$dataRowChunkSize, chunk_num_start = chunk_start_num)
  
  #param <- SerialParam(stop.on.error = TRUE)
  param <- MulticoreParam(workers = args$threads)
  
  anchorReadsAlignments <- rbindlist(bpiterate(ITER = my_iter, FUN = alignment_worker, BPPARAM = param)) %>% dplyr::select(qName, tName, strand, tStart, tEnd)
  bpstop(param)
  closeAllConnections()
  
  updateLog('Anchor read alignments completed.')
  updateLog(paste0(ppNum(n_distinct(anchorReadsAlignments$qName)), ' anchor reads returned one or more alignments.'))
  
  o$anchorReads <- left_join(anchorReads, anchorReadsAlignments, by = c('seqNum' = 'qName'), relationship = "many-to-many")
  
  o$anchorReads <- o$anchorReads[! is.na(o$anchorReads$tName)]
  
  if(nrow(o$anchorReads) == 0) stop('Error - no anchor reads aligned to the reference.')
  
  rm(param, anchorReads, anchorReadsAlignments)
  gc()
  
  updateLog('Limiting adrift reads to those with anchor read alignments.')
  
  # Limit adrift reads to those with anchor read alignments.
  adriftReads <- adriftReads[adriftReads$readID %in% o$anchorReads$readID]
  
  if(nrow(adriftReads) == 0) stop('Error -- no adrift reads remain after limiting reads to those with anchor read mates that aligned to the reference.')
  
  # Align unique adrift read sequences.
  
  updateLog('Starting adrift read alignments.')
  
  my_iter <- make_dt_iterator(adriftReads[! duplicated(adriftReads$seqNum)], chunk_size = args$dataRowChunkSize, chunk_num_start = chunk_start_num)
  param <- MulticoreParam(workers = args$threads)
  adriftReadsAlignments <- rbindlist(bpiterate(ITER = my_iter, FUN = alignment_worker, BPPARAM = param)) %>% dplyr::select(qName, tName, strand, tStart, tEnd)
  bpstop(param)
  closeAllConnections()
  
  updateLog('Adrift read alignments completed.')
  updateLog(paste0(ppNum(n_distinct(adriftReadsAlignments$qName)), ' adrift reads returned one or more alignments.'))
 
  o$adriftReads <- left_join(adriftReads, adriftReadsAlignments, by = c('seqNum' = 'qName'), relationship = "many-to-many")
  
  o$adriftReads <- o$adriftReads[! is.na(o$adriftReads$tName)]
  
  if(nrow(o$adriftReads) == 0) stop('Error - no adrift reads aligned to the reference.')
  
  rm(param, adriftReads, adriftReadsAlignments)
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
  
  saveRDS(o, file.path(args$outputDir, paste0(args$fileTag, '.rds')))
  updateLog('alignReads module completed.')
  write(date(), file.path(args$outputDir, paste0(args$fileTag, '.done')))
}

#-------------------------------------------------------------------------------

args <- parser$parse_args()
source(file.path(args$softwareRoot, 'lib.R'))

tryCatch({
  runModule()
  q(status = 0)
}, error = function(e) {
  message("Caught error: ", e$message)
  q(status = 1)
})

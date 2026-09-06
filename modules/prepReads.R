#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'ShortRead', 'parallel', 'data.table', 'BiocParallel', 'stringi', 'fst')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",               type = "character",     required = TRUE,          help = "Directory for output files")
parser$add_argument("--inputData",               type = "character",     required = TRUE,          help = "Path to demultiplex module's rds output file.")
parser$add_argument("--softwareRoot",            type = "character",     required = TRUE,          help = "Path to INSPIIRED2 installation.")
parser$add_argument("--threads",                 type = "integer",       default = 50,             help = "Number of threads to use.")
parser$add_argument("--fileTag",                 type  = "character",    default = "prepReads",    help = "String appended to output files in the outpt directory.")
parser$add_argument("--ramDiskPath",             type = "character",     default = "/dev/shm",     help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
parser$add_argument("--disableOverReadTrimming", action = "store_true",  default = FALSE,          help = "Disable over-read trimming.")
parser$add_argument("--disableVectorFilter",     action = "store_true",  default = FALSE,          help = "Disable removing reads due to similarity to vector sequence.")
parser$add_argument("--ORtrimPatternWidth",      type = "integer",       default = 8,              help = "Number of NTs used to build over-reading patterns.")
parser$add_argument("--ORseqMaxMismatch",        type = "double",        default = 0.10,           help = "Max mismatch percentage (0 .. 1) allowed to match over-reading patterns.")
parser$add_argument("--minReadLength",           type = "integer",       default = 30,             help = "Minial read length allowed.")
parser$add_argument("--vectorTestWidth",         type = "integer",       default = 25,             help = "Number of NTs at the end of reads to use to test for vector homology.")
parser$add_argument("--vectorTestMinPercentID",  type  = "double",       default = 90,             help = "Min. perecent ID (0 .. 100) to accept a vector alignment.")
parser$add_argument("--vectorTestMinCoverage",   type = "double",        default = 90,             help = "Min. test sequence converage (0 .. 100) to accept a vector alignment.")
parser$add_argument("--HMMparams",               type = "character",     default = 'none',         help = "Comma delimited shorthand containing HMM parmaters.")

runModule <- function(){
  startModule()
  
  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  on.exit({
    unlink(args$tmpDir,  recursive = TRUE, force = TRUE)
    unlink(args$logDir,  recursive = TRUE, force = TRUE) 
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  updateLog('Starting prepReads module.')

  resource_overlay()
  
  if(! file.exists(args$inputData))  stop(paste0('Error - the input data file (', args$inputData, ') does not exist.'))
  if(file.size(args$inputData) == 0) stop(paste0('Error - the input data file (', args$inputData, ') is empty.'))
  
  d <- readRDS(args$inputData)
  if(nrow(d) == 0) stop('Error -- input data had zero rows of data.')
  
  d$leaderSeqHMM <- as.character(d$leaderSeqHMM)
  d$vectorFastaFile <- as.character(d$vectorFastaFile)
  
  hmm_worker <- function(chunk, ...) {
    if(! dir.exists(file.path(args$logDir, paste0('chunk_', chunk$chunk_num)))) dir.create(file.path(args$logDir, paste0('chunk_', chunk$chunk_num)))
    logFile <- file.path(args$logDir, paste0('chunk_', chunk$chunk_num), 'log')
    
    hmmName <- as.character(chunk$data$leaderSeqHMM[1])
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tStarting HMM chunk with ', ppNum(nrow(chunk$data)), ' data rows.'), logFile = logFile)
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tHMM file: ', hmmName), logFile = logFile)
  
    ts <- tmpString()
    
    write(paste0('>', chunk$data$readID, '\n', chunk$data$anchorReadSeq), file = file.path(args$ramDisk, ts))
    
    if(args$HMMparams != 'none'){
      # HIV1_1-100_U5.hmm,1,5,10,30,TRUE,CA,2|HIV1_1-100_U3_RC.hmm,1,5,30,60,TRUE,CA,2
      updateLog(paste0('Processing --HMMparams ', args$HMMparams), logFile = logFile)
      
      a <- unlist(str_split(args$HMMparams, '\\|'))
   
      hmmList <- lapply(a, function(x){
        p <- unlist(str_split(x, '\\s*,\\s*'))
        p[2:length(p)]
      })
      names(hmmList) <- unlist(lapply(str_split(a, ','), '[[', 1))
      
      if(! hmmName %in% names(hmmList)) stop(paste0('Error - ', hmmName, ' was not defined in --HMMparams.'))
      
      p <- hmmList[[hmmName]]
      if(length(p) != 7) stop(paste0('Error - an argument string passed by --HMMparams did not contain 7 elements as expected for hmm: ', hmmName))
      
      args$HMMminStartPos      <- as.integer(p[1])
      args$HMMmaxStartPos      <- as.integer(p[2])
      args$HMMminFullBitScore  <- as.numeric(p[3])
      args$HMMmaxFullBitScore  <- as.numeric(p[4])
      args$HMMmatchEnd         <- ifelse(grepl('TRUE', p[5], ignore.case = TRUE), TRUE, FALSE)
      args$HMMmatchTerminalSeq <- p[6]
      args$HMMmatchEndRadius   <- as.integer(p[7])
    } else if (file.exists(file.path(args$softwareRoot, 'data', 'hmms', sub('\\.hmm$', '.cfg', hmmName)))){
      updateLog(paste0('--HMMparams was not provided, processing default HMM settings from file ', hmmName), logFile = logFile)
      p <- readr::read_tsv(file.path(args$softwareRoot, 'data', 'hmms', sub('\\.hmm$', '.cfg', hmmName)), col_names = FALSE, show_col_types = FALSE)
      if(nrow(p) != 7 | ncol(p) != 2) stop(paste0('Error - the hmm cfg file for hmm: ', hmmName, ' did not have the expected dimensions.'))
      if(any(is.na(p))) stop(paste0('Error - the hmm cfg file for hmm: ', hmmName, ' contained one or more NA values when read.'))
      
      args$HMMminStartPos      <- as.integer(p[which(p$X1 == 'HMMminStartPos'),]$X2)
      args$HMMmaxStartPos      <- as.integer(p[which(p$X1 == 'HMMmaxStartPos'),]$X2)
      args$HMMminFullBitScore  <- as.numeric(p[which(p$X1 == 'HMMminFullBitScore'),]$X2)
      args$HMMmaxFullBitScore  <- as.numeric(p[which(p$X1 == 'HMMmaxFullBitScore'),]$X2)
      args$HMMmatchEnd         <- ifelse(grepl('TRUE', p[which(p$X1 == 'HMMmatchEnd'),]$X2, ignore.case = TRUE), TRUE, FALSE)
      args$HMMmatchTerminalSeq <- p[which(p$X1 == 'HMMmatchTerminalSeq'),]$X2
      args$HMMmatchEndRadius   <- as.integer(p[which(p$X1 == 'HMMmatchEndRadius'),]$X2)
    } else {
      stop(paste0('Error - could not determine processing parameters for hmm: ', hmmName))
    }
    
    args$HMMparams <- NULL
    HMMparams <- tibble(timeStamp = paste0(base::format(Sys.time(), "%m.%d.%Y"), ' [', timeElapsedString(), "]"),  param = names(args)[grepl('^HMM', names(args))])
    HMMparams$value <- unlist(args[names(args) %in% HMMparams$param])
    names(HMMparams) <- c('', 'parameter', 'value')
    updateLog(paste0('HMM paramters for ', hmmName, ':'), logFile = logFile)
    write_tsv(HMMparams, file = logFile, append = TRUE)
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tCalling nhmmer.'), logFile = logFile)
    
    comm <- paste0('nhmmer --dna --F1 1 --F2 1 --F3 1 -T -5 --incT -5 --nobias --popen 0.15 --pextend 0.05 --tblout ', 
                   file.path(args$ramDisk, paste0(ts, '.tbl')), ' ', file.path(args$softwareRoot, 'data', 'hmms', hmmName), ' ', 
                   file.path(args$ramDisk, ts), ' > ', file.path(args$ramDisk, paste0(ts, '.hmmSearch')))
    
    system(comm)
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tnhmmer completed'), logFile = logFile)
    
    o <- readr::read_table(file.path(args$ramDisk, paste0(ts, '.tbl')), col_names = FALSE, col_types = NULL, comment = "#", show_col_types = FALSE)
    
    invisible(file.remove(list.files(args$ramDisk, pattern = ts, full.names = TRUE)))
    
    if(nrow(o) == 0){
      updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tNo hits were returned by nhmmer.'), logFile = logFile)
      return(data.table())
    }
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tParsing nhmmer output.'), logFile = logFile)
    
    names(o) <- c('targetName', 'targetAcc', 'queryName', 'queryAcc', 'hmmStart', 'hmmEnd', 'targetStart', 'targetEnd', 'envStart', 'envEnd', 'seqLength', 'strand', 'fullEval', 'fullScore', 'bias', 'desc')
    
    # Handle neg strand flipping coords. 
    o1  <- o[o$strand == '+',]
    o2  <- o[o$strand == '-',]
    
    if(nrow(o2) > 0){
      o2x <- o2
      o2x$targetStart <- o2$targetEnd;   o2x$envStart <- o2$envEnd
      o2x$targetEnd   <- o2$targetStart; o2x$envEnd   <- o2$envStart
      o <- bind_rows(o1, o2x)
      invisible(rm(o2x))
    }
    
    invisible(rm(o1, o2))
    
    # Collapse duplicate hits.
    o <- group_by(o, targetName) %>% dplyr::slice_max(fullScore, n = 1, with_ties = FALSE) %>% ungroup()
    
    # Subset the data based on user scoring thresholds.
    o <- subset(o, targetStart >= args$HMMminStartPos     & 
                   targetStart <= args$HMMmaxStartPos     & 
                   fullScore   >= args$HMMminFullBitScore &
                   fullScore   <= args$HMMmaxFullBitScore)
    
    if(nrow(o) == 0){
      updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tNo nhmmer hits reamin after filtering on targetStart and fullScore.'), logFile = logFile)
      return(data.table())
    }
    
    # Read in the HMM so that we can test if we aligned to its end.
    h <- readLines(file.path(args$softwareRoot, 'data', 'hmms', as.character(chunk$data$leaderSeqHMM[1])))
    hmmLength <- as.integer(unlist(strsplit(h[grepl('^LENG', h)], '\\s+'))[2])
    hmmName <- unlist(strsplit(h[grepl('^NAME', h)], '\\s+'))[2]
    
    # If requested, limit HMM hits to those with alignments near the end of the HMM.
    if(args$HMMmatchEnd) o <- o[abs(hmmLength - o$hmmEnd) <= args$HMMmatchEndRadius,]
    
    if(nrow(o) == 0){
      updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tNo nhmmer hits reamin after requiring a match to the full HMM.'), logFile = logFile)
      return(data.table())
    }
    
    # Limit reads to those with HMM matches
    chunk$data <- chunk$data[chunk$data$readID %in% o$targetName]
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\t', ppNum(nrow(chunk$data)), ' data rows remain after removing those with a significant HMM hit.'), logFile = logFile)
    
    chunk$data <- left_join(chunk$data, dplyr::select(o, targetName, targetStart, targetEnd), by = c('readID' = 'targetName'))
    
    if(! grepl('none', args$HMMmatchTerminalSeq, ignore.case = TRUE)){
       updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tSearching for requested terminal sequence match for "', args$HMMmatchTerminalSeq, '"'), logFile = logFile)
       chunk$data$anchorReadSeq <- toupper(chunk$data$anchorReadSeq)
       args$HMMmatchTerminalSeq <- toupper(args$HMMmatchTerminalSeq)
      
       terminal_matchSeq <- substr(chunk$data$anchorReadSeq,  (chunk$data$targetEnd - (nchar(args$HMMmatchTerminalSeq) - 1) - args$HMMmatchEndRadius), (chunk$data$targetEnd + args$HMMmatchEndRadius))
       ends <- stringr::str_locate(terminal_matchSeq, args$HMMmatchTerminalSeq)[, 2]
       
       i <- ! is.na(ends)
       chunk$data <- chunk$data[i]
       ends <- ends[i]
       
       chunk$data$targetEnd <- chunk$data$targetEnd - (nchar(args$HMMmatchTerminalSeq) + args$HMMmatchEndRadius) + ends
       updateLog(paste0('<data chunk #', chunk$chunk_num, '>\t', ppNum(nrow(chunk$data)), ' data rows remain after requiring a terminal sequence match.'), logFile = logFile)
    }
    
    return(chunk$data)
  }
  
  updateLog(paste0('Starting asynchronous HMM calculations. Data chunk logs can be found in ', args$logDir, '/'))
  updateLog('Collated data chunk logs will appear below when done.')
  
  d <- rbindlist(lapply(split(d, d$leaderSeqHMM), function(x){
         chunk_start_num <- 0
         f <- list.files(args$logDir, full.names = FALSE)
         if(length(f) > 0) chunk_start_num <- max(as.integer(str_extract(f, '\\d+$')))
    
         my_iter <- make_dt_iterator(x, chunk_size = ceiling(nrow(x)/args$threads), chunk_num_start = chunk_start_num)
    
         ### param <- SerialParam(stop.on.error = TRUE) # Use SerialParam() for browser() statements.
         param <- MulticoreParam(workers = args$threads)
    
         results <- bpiterate(ITER = my_iter, FUN = hmm_worker, BPPARAM = param)
    
         bpstop(param)
         closeAllConnections()
    
        rbindlist(results)
    }))
  
  if(nrow(d) == 0){
    msg <- 'Error -- no reads remain after selecting for reads with significant HMM signatures.'
    updateLog(msg)
    stop(msg)
  }
  
  collated_logs <- unlist(lapply(list.files(args$logDir, full.names = TRUE, recursive = TRUE, pattern = '^log$'), readLines))
  unlink(list.files(args$logDir, full.names = TRUE), recursive = TRUE)
  updateLog(collated_logs)
  
  d$leaderSeq <- substr(d$anchorReadSeq, 1, d$targetEnd)
  d$anchorReadSeq <- substr(d$anchorReadSeq, d$targetEnd + 1, nchar(d$anchorReadSeq))
  
  d$targetStart <- NULL
  d$targetEnd <- NULL

  if(! args$disableOverReadTrimming){
    updateLog('Trimming over reading.')
  
    d$anchorReadTrimSeq <- as.character(subseq(reverseComplement(DNAStringSet(d$linker2)), 1, args$ORtrimPatternWidth))
    d$adriftReadTrimSeq <- as.character(subseq(reverseComplement(DNAStringSet(d$leaderSeq)), 1, args$ORtrimPatternWidth))
  
    d <- rbindlist(lapply(split(d, d$anchorReadTrimSeq), function(x){
           maxMisMatch <- ceiling(args$ORtrimPatternWidth * args$ORseqMaxMismatch)
           matches <- vmatchPattern(x$anchorReadTrimSeq[1], DNAStringSet(x$anchorReadSeq), max.mismatch = maxMisMatch, fixed = TRUE)
           match_starts <- unlist(lapply(startIndex(matches), function(m) if(length(m) > 0) tail(m, 1) else NA))
           toTrimIndex <- ! is.na(match_starts) & match_starts > 1
    
          if(any(toTrimIndex)){
            x[toTrimIndex]$anchorReadSeq <- substr(x[toTrimIndex]$anchorReadSeq, 1, match_starts[toTrimIndex] - 1)
          }
    
          x
       }))
    
    if(nrow(d) == 0){
      msg <- 'Error - no reads remain after anchorRead overTrimming filter.'
      updateLog(msg)
      stop(msg)
    }
  
    d$anchorReadTrimSeq <- NULL
  
    d <- rbindlist(lapply(split(d, d$adriftReadTrimSeq), function(x){
           maxMisMatch <- ceiling(args$ORtrimPatternWidth * args$ORseqMaxMismatch)
           matches <- vmatchPattern(x$adriftReadTrimSeq[1], DNAStringSet(x$adriftReadSeq), max.mismatch = maxMisMatch, fixed = TRUE)
           match_starts <- unlist(lapply(startIndex(matches), function(m) if(length(m) > 0) tail(m, 1) else NA))
           toTrimIndex <- ! is.na(match_starts) & match_starts > 1
    
          if(any(toTrimIndex)){
            x[toTrimIndex]$adriftReadSeq <- substr(x[toTrimIndex]$adriftReadSeq, 1, match_starts[toTrimIndex] - 1)
          }
    
          x
        }))
  
    d$adriftReadTrimSeq <- NULL
  
    keep_idx <- which(nchar(d$anchorReadSeq) >= args$minReadLength & nchar(d$adriftReadSeq) >= args$minReadLength)
  
    d <- d[keep_idx]
    
    if(nrow(d) == 0){
      msg <- 'Error - no reads remain after adriftRead overTrimming filter.'
      updateLog(msg)
      stop(msg)
    }
  }
  
 
  if(! args$disableVectorFilter){
    d <- rbindlist(lapply(split(d, d$vectorFastaFile), function(x){
      ts <- tmpString()
      system2("makeblastdb", args = c("-in",  file.path(args$softwareRoot, 'data', 'vectors', x$vectorFastaFile[1]), "-dbtype", "nucl", "-out", file.path(args$ramDisk, ts)), stdout = FALSE, stderr = FALSE)
      
      x$testSeq <- substr(x$anchorReadSeq, (nchar(x$anchorReadSeq) - args$vectorTestWidth + 1), nchar(x$anchorReadSeq))
      
      x2 <- dplyr::select(x, readID, testSeq)
      x2 <- x2[! duplicated(x2$testSeq),]
      write(paste0('>', x2$readID, '\n', x2$testSeq), file = file.path(args$ramDisk, paste0(ts, '.fasta')))
      
      blastn_out <- run_blastn_parallel( file.path(args$ramDisk, paste0(ts, '.fasta')), file.path(args$ramDisk, ts), paste0("-word_size 8 -perc_identity ", args$vectorTestMinPercentID, " -gapopen 10 -gapextend 6 -evalue 10 -dust no -soft_masking false"), threads = args$threads)
      
      if(nrow(blastn_out) > 0){
        blastn_out$coverage <- (blastn_out$len / args$vectorTestWidth) * 100         # Calculate alignment coverage 
        blastn_out <- blastn_out[blastn_out$coverage >= args$vectorTestMinCoverage]  # Filter for alignments >= args$vectorTestMinCoverage
        
        if(nrow(blastn_out) > 0){
          x2 <- x2[x2$readID %in% blastn_out$qName]                                  # Limit original test sequences to those with significant hits
          x2$readID <- NULL                                                          # Remove read ID and add vectorHit to create a two column table that can be joined to original
          x2$vectorHit <- TRUE
          x <- left_join(x, x2, by = 'testSeq')                                      # Join table by test sequence.
          x$vectorHit <- ifelse(is.na(x$vectorHit), FALSE, TRUE)                     # Create a boolean to show if a read is a likely internal read
        } else {
          x$vectorHit <- FALSE
        }
      } else {
        x$vectorHit <- FALSE
      }
      
      invisible(file.remove(list.files(args$ramDisk, pattern = ts, full.names = TRUE)))
      x$testSeq <- NULL
      x
    }))
    
    if(nrow(d) == 0){
      msg <- 'Error - no reads remain after vector filter.'
      updateLog(msg)
      stop(msg)
    }
    
    updateLog(paste0(sprintf("%.1f%%", (sum(d$vectorHit) / nrow(d))*100), ' anchorRead ends matched the vector sequences.'))
    updateLog('Writing output.')
    
    write_tsv(d[d$vectorHit == TRUE], file.path(args$outputDir, paste0(args$fileTag, '_vectorHitReads.tsv.gz')))
    
    d <- d[d$vectorHit == FALSE]
    
    d$vectorHit <- NULL
    d$linker1   <- NULL
    d$linker2   <- NULL
  }
  
  d$trial     <- as.factor(d$trial)
  d$subject   <- as.factor(d$subject)
  d$sample    <- as.factor(d$sample)
  d$replicate <- as.factor(d$replicate)

  saveRDS(d, file.path(args$outputDir, paste0(args$fileTag, '.rds')))
  updateLog('prepReads module completed.')
  write(date(), file.path(args$outputDir, paste0(args$fileTag, '.done')))
}

#-------------------------------------------------------------------------------

args <- parser$parse_args()
source(file.path(args$softwareRoot, 'lib', 'common.R'))

tryCatch({
  runModule()
}, error = function(e) {
  cat("ERROR: ", conditionMessage(e), "\n", sep = "", file = stderr())
  flush(stderr())
  quit(save = "no", status = 1, runLast = FALSE)
})
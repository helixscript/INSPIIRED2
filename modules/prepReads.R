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
  doneFile <- file.path(args$outputDir, paste0(args$fileTag, '.done'))
  if(file.exists(doneFile) && unlink(doneFile) != 0) stop('Error - could not remove stale completion marker: ', doneFile)
  
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
  requiredHMMcols <- c("leaderSeqHMM", "vectorFastaFile")
  missingHMMcols <- setdiff(requiredHMMcols, names(d))
  if(length(missingHMMcols))
    stop("Error - input data is missing required column(s): ", paste(missingHMMcols, collapse = ", "))
  
  d$leaderSeqHMM <- as.character(d$leaderSeqHMM)
  d$vectorFastaFile <- as.character(d$vectorFastaFile)
  
  buildHMMParameterTable <- function(hmmNames){
    parseOverrides <- function(x){
      if(is.null(x) || length(x) != 1L || !nzchar(trimws(x)) || tolower(trimws(x)) == "none") return(list())
      entries <- strsplit(x, "|", fixed = TRUE)[[1]]
      parsed <- lapply(entries, function(entry){
        p <- trimws(strsplit(entry, ",", fixed = TRUE)[[1]])
        if(length(p) != 8L)
          stop("Error - --HMMparams entry must contain an HMM name followed by 7 parameters: ", entry)
        p
      })
      hmmNames <- vapply(parsed, `[`, character(1), 1L)
      if(any(!nzchar(hmmNames))) stop("Error - --HMMparams contains an empty HMM name.")
      if(anyDuplicated(hmmNames)) stop("Error - the same HMM was defined more than once in --HMMparams.")
      parsed <- lapply(parsed, function(x) x[-1L])
      names(parsed) <- hmmNames
      parsed
    }
    
    readHMMLength <- function(hmmName){
      hmmFile <- file.path(args$softwareRoot, "data", "hmms", hmmName)
      if(!file.exists(hmmFile)) stop("Error - HMM file does not exist: ", hmmFile)
      h <- readLines(hmmFile, warn = FALSE)
      x <- grep("^LENG\\s+", h, value = TRUE)
      if(length(x) != 1L) stop("Error - could not uniquely determine LENG for hmm: ", hmmName)
      n <- suppressWarnings(as.integer(strsplit(trimws(x), "\\s+")[[1]][2]))
      if(is.na(n) || n < 1L) stop("Error - could not parse a positive HMM length for hmm: ", hmmName)
      n
    }
    
    makeRow <- function(hmmName, p, source){
      if(length(p) != 7L) stop("Error - expected 7 HMM parameters for hmm: ", hmmName)
      p <- trimws(as.character(p))
      if(!grepl("^(TRUE|FALSE)$", p[5], ignore.case = TRUE))
        stop("Error - HMMmatchEnd must be TRUE or FALSE for hmm: ", hmmName)
      
      integerValues <- suppressWarnings(as.numeric(p[c(1L, 2L, 7L)]))
      if(anyNA(integerValues) || any(!is.finite(integerValues)) || any(integerValues != trunc(integerValues)))
        stop("Error - HMM start positions and HMMmatchEndRadius must be finite whole numbers for hmm: ", hmmName)
      
      terminalSeq <- toupper(p[6])
      if(!identical(tolower(terminalSeq), "none") && !grepl("^[ACGTN]+$", terminalSeq))
        stop("Error - HMMmatchTerminalSeq must be 'none' or a sequence containing only A, C, G, T, or N for hmm: ", hmmName)
      
      z <- data.table(
        leaderSeqHMM = hmmName,
        HMMminStartPos = suppressWarnings(as.integer(integerValues[1L])),
        HMMmaxStartPos = suppressWarnings(as.integer(integerValues[2L])),
        HMMminFullBitScore = suppressWarnings(as.numeric(p[3])),
        HMMmaxFullBitScore = suppressWarnings(as.numeric(p[4])),
        HMMmatchEnd = grepl("^TRUE$", p[5], ignore.case = TRUE),
        HMMmatchTerminalSeq = terminalSeq,
        HMMmatchEndRadius = suppressWarnings(as.integer(integerValues[3L])),
        hmmLength = readHMMLength(hmmName),
        parameterSource = source
      )
      
      if(any(is.na(z[, .(HMMminStartPos, HMMmaxStartPos, HMMminFullBitScore,
                         HMMmaxFullBitScore, HMMmatchEndRadius, hmmLength)])) ||
         any(!is.finite(c(z$HMMminFullBitScore, z$HMMmaxFullBitScore))))
        stop("Error - one or more HMM parameters could not be parsed for hmm: ", hmmName)
      if(z$HMMminStartPos < 1L || z$HMMmaxStartPos < z$HMMminStartPos)
        stop("Error - invalid HMM start-position range for hmm: ", hmmName)
      if(z$HMMmaxFullBitScore < z$HMMminFullBitScore)
        stop("Error - invalid HMM score range for hmm: ", hmmName)
      if(z$HMMmatchEndRadius < 0L)
        stop("Error - HMMmatchEndRadius must be >= 0 for hmm: ", hmmName)
      if(!nzchar(z$HMMmatchTerminalSeq))
        stop("Error - HMMmatchTerminalSeq was empty for hmm: ", hmmName)
      z
    }
    
    hmmNames <- unique(as.character(hmmNames))
    if(anyNA(hmmNames) || any(!nzchar(hmmNames)))
      stop("Error - leaderSeqHMM values must be non-missing and non-empty.")
    overrides <- parseOverrides(args$HMMparams)
    missingOverrides <- setdiff(hmmNames, names(overrides))
    if(length(overrides) && length(missingOverrides))
      stop("Error - --HMMparams did not define parameters for input HMM(s): ",
           paste(missingOverrides, collapse = ", "), ".")
    unusedOverrides <- setdiff(names(overrides), hmmNames)
    if(length(unusedOverrides))
      updateLog(paste0("Ignoring --HMMparams entries not used by this input: ",
                       paste(unusedOverrides, collapse = ", "), "."))
    
    expected <- c("HMMminStartPos", "HMMmaxStartPos", "HMMminFullBitScore",
                  "HMMmaxFullBitScore", "HMMmatchEnd", "HMMmatchTerminalSeq",
                  "HMMmatchEndRadius")
    
    parameterRows <- lapply(hmmNames, function(hmmName){
      hmmFile <- file.path(args$softwareRoot, "data", "hmms", hmmName)
      if(!file.exists(hmmFile)) stop("Error - HMM file does not exist: ", hmmFile)
      if(hmmName %in% names(overrides))
        return(makeRow(hmmName, overrides[[hmmName]], "--HMMparams"))
      
      cfgFile <- file.path(args$softwareRoot, "data", "hmms", sub("\\.hmm$", ".cfg", hmmName))
      if(!file.exists(cfgFile))
        stop("Error - could not determine processing parameters for hmm: ", hmmName,
             ". No matching --HMMparams entry or cfg file was found.")
      
      p <- readr::read_tsv(cfgFile, col_names = FALSE,
                           col_types = readr::cols(.default = readr::col_character()),
                           show_col_types = FALSE, progress = FALSE)
      if(nrow(p) != 7L || ncol(p) != 2L)
        stop("Error - the hmm cfg file for hmm: ", hmmName, " did not have the expected dimensions.")
      if(any(is.na(p)))
        stop("Error - the hmm cfg file for hmm: ", hmmName, " contained one or more NA values.")
      names(p) <- c("name", "value")
      if(anyDuplicated(p$name) || !setequal(p$name, expected))
        stop("Error - the hmm cfg file for hmm: ", hmmName, " did not contain the expected parameter names.")
      makeRow(hmmName, p$value[match(expected, p$name)], "cfg")
    })
    
    rbindlist(parameterRows, use.names = TRUE, fill = FALSE)
  }
  
  hmmParameters <- buildHMMParameterTable(d$leaderSeqHMM)
  
  hmm_worker <- function(chunk, hmmParameters, ...) {
    if(! dir.exists(file.path(args$logDir, paste0('chunk_', chunk$chunk_num)))) dir.create(file.path(args$logDir, paste0('chunk_', chunk$chunk_num)))
    logFile <- file.path(args$logDir, paste0('chunk_', chunk$chunk_num), 'log')
    
    hmmName <- as.character(chunk$data$leaderSeqHMM[1])
    hp <- hmmParameters[leaderSeqHMM == hmmName]
    if(nrow(hp) != 1L) stop("Error - expected exactly one parameter record for HMM: ", hmmName)
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tStarting HMM chunk with ', ppNum(nrow(chunk$data)), ' data rows.'), logFile = logFile)
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tHMM file: ', hmmName), logFile = logFile)
    
    ts <- tmpString()
    
    write(paste0('>', chunk$data$readID, '\n', chunk$data$anchorReadSeq), file = file.path(args$ramDisk, ts))
    
    parameterNames <- c("HMMminStartPos", "HMMmaxStartPos", "HMMminFullBitScore",
                        "HMMmaxFullBitScore", "HMMmatchEnd", "HMMmatchTerminalSeq",
                        "HMMmatchEndRadius")
    HMMparameterLog <- tibble(
      timeStamp = paste0(base::format(Sys.time(), "%m.%d.%Y"), ' [', timeElapsedString(), "]"),
      parameter = parameterNames,
      value = vapply(parameterNames, function(nm) as.character(hp[[nm]][1L]), character(1))
    )
    names(HMMparameterLog)[1] <- ''
    updateLog(paste0('HMM parameters for ', hmmName, ' resolved from ', hp$parameterSource[1L], ':'), logFile = logFile)
    write_tsv(HMMparameterLog, file = logFile, append = TRUE)
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tCalling nhmmer.'), logFile = logFile)
    
    comm <- paste0('nhmmer --cpu 1 --dna --F1 1 --F2 1 --F3 1 -T -5 --incT -5 --nobias --popen 0.15 --pextend 0.05 --tblout ', 
                   file.path(args$ramDisk, paste0(ts, '.tbl')), ' ', file.path(args$softwareRoot, 'data', 'hmms', hmmName), ' ', 
                   file.path(args$ramDisk, ts), ' > ', file.path(args$ramDisk, paste0(ts, '.hmmSearch')))
    
    status <- system(comm)
    requireCommandSuccess(status, "nhmmer")
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tnhmmer completed'), logFile = logFile)
    
    o <- readr::read_table(file.path(args$ramDisk, paste0(ts, '.tbl')), col_names = FALSE, col_types = NULL, comment = "#", show_col_types = FALSE)
    
    invisible(file.remove(list.files(args$ramDisk, pattern = ts, full.names = TRUE)))
    
    if(nrow(o) == 0){
      updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tNo hits were returned by nhmmer.'), logFile = logFile)
      return(data.table())
    }
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tParsing nhmmer output.'), logFile = logFile)
    
    names(o) <- c('targetName', 'targetAcc', 'queryName', 'queryAcc', 'hmmStart', 'hmmEnd', 'targetStart', 'targetEnd', 'envStart', 'envEnd', 'seqLength', 'strand', 'fullEval', 'fullScore', 'bias', 'desc')
    
    nMinus <- sum(o$strand == '-')
    if(nMinus > 0) updateLog(paste0('Ignoring ', ppNum(nMinus), ' reverse-strand HMM hits.'), logFile = logFile)
    
    o <- o[o$strand == '+', ]
    if(nrow(o) == 0) return(data.table())
    
    o <- group_by(o, targetName) %>%
      slice_max(fullScore, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    # Subset the data based on user scoring thresholds.
    o <- subset(o, targetStart >= hp$HMMminStartPos[1L]     &
                  targetStart <= hp$HMMmaxStartPos[1L]     &
                  fullScore   >= hp$HMMminFullBitScore[1L] &
                  fullScore   <= hp$HMMmaxFullBitScore[1L])
    
    if(nrow(o) == 0){
      updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tNo nhmmer hits reamin after filtering on targetStart and fullScore.'), logFile = logFile)
      return(data.table())
    }
    
    # If requested, limit HMM hits to those with alignments near the end of the HMM.
    if(isTRUE(hp$HMMmatchEnd[1L]))
      o <- o[abs(hp$hmmLength[1L] - o$hmmEnd) <= hp$HMMmatchEndRadius[1L],]
    
    if(nrow(o) == 0){
      updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tNo nhmmer hits reamin after requiring a match to the full HMM.'), logFile = logFile)
      return(data.table())
    }
    
    # Limit reads to those with HMM matches
    chunk$data <- chunk$data[chunk$data$readID %in% o$targetName]
    
    updateLog(paste0('<data chunk #', chunk$chunk_num, '>\t', ppNum(nrow(chunk$data)), ' data rows remain after removing those with a significant HMM hit.'), logFile = logFile)
    
    chunk$data <- left_join(chunk$data, dplyr::select(o, targetName, targetStart, targetEnd), by = c('readID' = 'targetName'))
    
    terminalSeq <- as.character(hp$HMMmatchTerminalSeq[1L])
    radius <- as.integer(hp$HMMmatchEndRadius[1L])
    if(!identical(tolower(terminalSeq), 'none')){
      terminalSeq <- toupper(terminalSeq)
      updateLog(paste0('<data chunk #', chunk$chunk_num, '>\tSearching for requested terminal sequence match for "', terminalSeq, '"'), logFile = logFile)
      chunk$data$anchorReadSeq <- toupper(chunk$data$anchorReadSeq)
      
      terminal_matchSeq <- substr(chunk$data$anchorReadSeq,
                                  chunk$data$targetEnd - (nchar(terminalSeq) - 1L) - radius,
                                  chunk$data$targetEnd + radius)
      ends <- stringr::str_locate(terminal_matchSeq, stringr::fixed(terminalSeq))[, 2]
      
      i <- ! is.na(ends)
      chunk$data <- chunk$data[i]
      ends <- ends[i]
      
      chunk$data$targetEnd <- chunk$data$targetEnd - (nchar(terminalSeq) + radius) + ends
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
    param <- MulticoreParam(workers = args$threads, stop.on.error = TRUE)
    
    results <- tryCatch(
      bpiterate(ITER = my_iter, FUN = hmm_worker, hmmParameters = hmmParameters, BPPARAM = param),
      finally = {
        try(bpstop(param), silent = TRUE)
        closeAllConnections()
      }
    )
    
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
      status <- system2("makeblastdb", args = c("-in",  file.path(args$softwareRoot, 'data', 'vectors', x$vectorFastaFile[1]), "-dbtype", "nucl", "-out", file.path(args$ramDisk, ts)), stdout = FALSE, stderr = FALSE)
      requireCommandSuccess(status, "makeblastdb")
      
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
    
    updateLog(paste0(sprintf("%.1f%%", (sum(d$vectorHit) / nrow(d))*100), ' anchorRead ends matched the vector sequences.'))
    updateLog('Writing output.')
    
    write_tsv(d[d$vectorHit == TRUE], file.path(args$outputDir, paste0(args$fileTag, '_vectorHitReads.tsv.gz')))
    
    d <- d[d$vectorHit == FALSE]
    
    if(nrow(d) == 0){
      msg <- 'Error - no reads remain after vector filter.'
      updateLog(msg)
      stop(msg)
    }
    
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

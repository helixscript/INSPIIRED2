#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'ShortRead', 'parallel', 'data.table', 'BiocParallel', 'stringi', 'fst')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",                    type = "character",     required = TRUE,                  help = "Directory for output files")
parser$add_argument("--sampleData",                   type = "character",     required = TRUE,                  help = "Sample definition file")
parser$add_argument("--I1",                           type = "character",     required = TRUE,                  help = "Path to the Index1 read FASTQ file")
parser$add_argument("--R1",                           type = "character",     required = TRUE,                  help = "Path to the Forward read FASTQ file")
parser$add_argument("--R2",                           type = "character",     required = TRUE,                  help = "Path to the Reverse read FASTQ file")
parser$add_argument("--softwareRoot",                 type = "character",     required = TRUE,                  help = "Path to AAVengeR installation.")
parser$add_argument("--threads",                      type = "integer",       default = 50,                     help = "Number of threads to use.")
parser$add_argument("--fileTag",                      type = "character",     default = "demultiplex",          help = "String appended to output files in the outpt directory.")
parser$add_argument("--index1ReadMaxMismatch",        type = "integer",       default = 1,                      help = "Number of allowed mismatches to the I1 barcode sequence.")
parser$add_argument("--disableAutoBarcodeOrt",        action = "store_true",  default = FALSE,                  help = "Subsample the data an automatically determine if I1 barcodes need to be reverse complimented.")
parser$add_argument("--disablePostUmiLinker",         action = "store_true",  default = FALSE,                  help = "Disable the requirement to match the post-UMI linker sequence.")
parser$add_argument("--postUmiLinkerMaxMismatch",     type = "integer",       default = 1,                      help = "Number of allowed mismatches to the linker sequence following the UMI sequence.")
parser$add_argument("--qualTrimHalfWidth",            type = "integer",       default = 2,                      help = "Half width of NT window slid along sequence during quality trimming.")
parser$add_argument("--qualTrimEvents",               type = "integer",       default = 3,                      help = "Number of failing events within a window to trigger trimming.")
parser$add_argument("--qualTrimScore",                type = "integer",       default = 10,                     help = "Qual code afterwhich NTs are trimmed.")
parser$add_argument("--polyGfilterPattern",           type = "character",     default = "G{5,}[ATCN]?G{5,}.*$", help = "Pattern to recognize poly-G tail NTs.")
parser$add_argument("--disablePolyGfilter",           action = "store_true",  default = FALSE,                  help = "Disable poly-G filter.")
parser$add_argument("--correctGolayIndexReads",       action = "store_true",  default = FALSE,                  help = "Use a Golay correction algorithm to correct barcode sequences.")
parser$add_argument("--disableAdriftReadLinkers",     action = "store_true",  default = FALSE,                  help = "Use the unique linker sequences on adrift reads with barcode sequences.")
parser$add_argument("--disableSequenceCollapse",      action = "store_true",  default = FALSE,                  help = "Disable the collapse of duplicate sequences.")
parser$add_argument("--adriftReadLinkerMaxMismatch",  type = "integer",       default = 1,                      help = "Number of allowed mismatches to the linker sequence.")
parser$add_argument("--ramDiskPath",                  type = "character",     default = "/dev/shm",             help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")

runModule <- function(){
  startModule()
  
  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  on.exit({
    unlink(args$tmpDir,  recursive = TRUE, force = TRUE)
    unlink(args$logDir,  recursive = TRUE, force = TRUE) 
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  updateLog('Starting demultiplex module.')
  
  # Sanity checks...
  
  if(! file.exists(args$sampleData)){
    updateLog(paste0('Error - sampleData file does not exits. Provided path: ', args$sampleData)) 
    quit(status = 1)
  }
  
  sampleData <- read_tsv(args$sampleData, show_col_types = FALSE)
  
  knownVectors <- list.files(file.path(args$softwareRoot, 'data', 'vectors'))
  knownHMMs <- list.files(file.path(args$softwareRoot, 'data', 'hmms'))
  
  if(! file.exists(args$I1)){
    updateLog(paste0('Error - I1 sequencing file does not exits. Provided path: ', args$I1)) 
    quit(status = 1)
  }
  
  if(! file.exists(args$R1)){
    updateLog(paste0('Error - R1 sequencing file does not exits. Provided path: ', args$R1)) 
    quit(status = 1)
  }
  
  if(! file.exists(args$R2)){
    updateLog(paste0('Error - R2 sequencing file does not exits. Provided path: ', args$R2)) 
    quit(status = 1)
  }
  
  requiredFields <- c('trial', 'subject', 'sample', 'replicate', 'adriftReadLinkerSeq', 'index1Seq', 'refGenome', 'vectorFastaFile', 'leaderSeqHMM', 'mode')
  
  if(! all(requiredFields %in% names(sampleData))){
    updateLog(paste0('Error - these required fields are missing from the sample data file: ', paste0(sQuote(requiredFields[! requiredFields %in% names(sampleData)]), collapse = ', ')))
    quit(status = 1)
  }
 
  if(! all(sampleData$vectorFastaFile %in% knownVectors)){
    missingVectors <- paste0(unique(sampleData$vectorFastaFile)[! unique(sampleData$vectorFastaFile) %in% knownVectors], collapse = ', ')
    updateLog(paste0('Error - These vector file names in the sample data file were not found in ', file.path(args$softwareRoot, 'data', 'vectors'), ': ', missingVectors)) 
    quit(status = 1)
  }
  
  if(! all(sampleData$leaderSeqHMM %in% knownHMMs)){
    missingHMMs <- paste0(unique(sampleData$leaderSeqHMM)[! unique(sampleData$leaderSeqHMM) %in% knownHMMs], collapse = ', ')
    updateLog(paste0('Error - These hmm file names in the sample data file were not found in ', file.path(args$softwareRoot, 'data', 'hmms'), ': ', missingHMMs)) 
    quit(status = 1)
  }
  
  args$reverseComplementI1 <- FALSE
  dataStreamChunkSize <- 1e6
  
  if(! args$disableAutoBarcodeOrt){
    updateLog('Determining if I1 barcodes need to be reverse-complimented.')
    
    stream_I1 <- FastqStreamer(args$I1, n = dataStreamChunkSize)
    I1_test <- yield(stream_I1)
    close(stream_I1)
    
    a <- sum(as.character(I1_test@sread) %in% sampleData$index1Seq)
    b <- sum(as.character(reverseComplement(I1_test)@sread) %in% sampleData$index1Seq)
    
    if(b > a){
      updateLog('I1 barcodes set to their reverse-compliment.')
      args$reverseComplementI1 <- TRUE
    }
    
    suppressWarnings(rm(a, b, I1_test))
  }
  
  stream_I1 <- FastqStreamer(args$I1, n = dataStreamChunkSize)
  stream_R1 <- FastqStreamer(args$R1, n = dataStreamChunkSize)
  stream_R2 <- FastqStreamer(args$R2, n = dataStreamChunkSize)
  
  chunk_num <- 0
  total_reads <- 0
  
  aav_iterator <- function() {
    chunk_I1 <- yield(stream_I1)
    if (length(chunk_I1) == 0) return(NULL)
    chunk_num <<- chunk_num + 1
    total_reads <<- total_reads + length(chunk_I1)
    
    chunk_R1 <- yield(stream_R1)
    chunk_R2 <- yield(stream_R2)
    
    if (length(chunk_I1) != length(chunk_R1) || length(chunk_I1) != length(chunk_R2)) {
      stop("FASTQ files are out of sync! Check for file truncation.")
    }
    
    list(I1 = chunk_I1, R1 = chunk_R1, R2 = chunk_R2, chunk_num = chunk_num)
  }
  
  aav_worker <- function(chunk_data, sampleData, args) {
    for (p in c('dplyr', 'ShortRead', 'data.table', 'stringi', 'fst')) suppressPackageStartupMessages(library(p, character.only = TRUE))
    
    ppNum <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
    
    cI1 <- chunk_data$I1
    cR1 <- chunk_data$R1
    cR2 <- chunk_data$R2
    chunk_num <- chunk_data$chunk_num
    
    if(! dir.exists(file.path(args$logDir, paste0('chunk_', chunk_num)))) dir.create(file.path(args$logDir, paste0('chunk_', chunk_num)))
    logFile <- file.path(args$logDir, paste0('chunk_', chunk_num), 'log')
      
    if(args$reverseComplementI1) cI1 <- reverseComplement(cI1)
    
    updateLog(paste0('<data chunk #', chunk_num, '>\tSeparated reads into groups (', ppNum(length(cI1)), ' reads).'), logFile = logFile)
    
    suppressWarnings(rm(chunk_data))
    
    if(length(cI1) == 0) break
    
    clean_ids <- BStringSet(gsub('^[^:]+:[^:]+:[^:]+:|\\s.+$', '', as.character(cI1@id)))
    cI1@id <- clean_ids
    cR1@id <- clean_ids
    cR2@id <- clean_ids
    
    updateLog(paste0('<data chunk #', chunk_num, '>\tQuality trimming reads. Qual code: ', args$qualTrimCode, ' half width: ', args$qualTrimHalfWidth, ' events: ', args$qualTrimEvents), logFile = logFile)
    
    rng1 <- trimTailw(cR1, args$qualTrimEvents, args$qualTrimCode, args$qualTrimHalfWidth, ranges = TRUE)
    rng2 <- trimTailw(cR2, args$qualTrimEvents, args$qualTrimCode, args$qualTrimHalfWidth, ranges = TRUE)
    
    # Read must be long enough to extract linker sequences.
    minReadLength <- max(nchar(sampleData$adriftReadLinkerSeq)) + 1
    
    keep_idx <- which(width(rng1) >= minReadLength & width(rng2) >= minReadLength)
    updateLog(paste0('<data chunk #', chunk_num, '>\tRanges kept after trimming: ', length(keep_idx)), logFile = logFile)
    
    if(length(keep_idx) == 0) return(data.table())
    
    cI1 <- cI1[keep_idx]
    cR1 <- narrow(cR1[keep_idx], start = start(rng1[keep_idx]), end = end(rng1[keep_idx]))
    cR2 <- narrow(cR2[keep_idx], start = start(rng2[keep_idx]), end = end(rng2[keep_idx]))
    
    suppressWarnings(rm(rng1, rng2, keep_idx))
    
    updateLog(paste0('<data chunk #', chunk_num, '>\tProcessing sample table.'), logFile = logFile)
    
    sampleData$nowNum <- 1:nrow(sampleData)
    
    d <- rbindlist(lapply(split(sampleData, 1:nrow(sampleData)), function(x){
      # Determine if the linker sequence in the sample data table has an UMI defined by Ns.
      # eg. GCGGTGCATCTCTTATAGCGNNNNNNNNNNNNCTCCGCTTAAGGGACT
      
      x$adriftReadLinkerSeq <- toupper(x$adriftReadLinkerSeq)
      n_loc <- str_locate(x$adriftReadLinkerSeq, "N+")
      
      # If an UMI is located, parse the linker sequence.
      if(! is.na(n_loc[1])){
        coords <- list(pre_n  = c(start = 1,            end = n_loc[1] - 1),
                       n_run  = c(start = n_loc[1],     end = n_loc[2]),
                       post_n = c(start = n_loc[2] + 1, end = str_length(x$adriftReadLinkerSeq)))
      }
      
      # I1 barcode test.
      i1 <- vcountPattern(x$index1Seq, cI1@sread, max.mismatch = args$index1ReadMaxMismatch) > 0
      
      # Adrift read linker filter.
      if(! args$disableAdriftReadLinkers){
        if(! is.na(n_loc[1])){
          i2 <- vcountPattern(substr(x$adriftReadLinkerSeq, 1, coords$pre_n[2]), subseq(cR1@sread, 1, coords$pre_n[2]), max.mismatch = args$adriftReadLinkerMaxMismatch) > 0
        } else {
          # Use the entire linker if an UMI was not found.
          i2 <- vcountPattern(x$adriftReadLinkerSeq, subseq(cR1@sread, 1, nchar(x$adriftReadLinkerSeq)), max.mismatch = args$adriftReadLinkerMaxMismatch) > 0
        }
      } else {
        i2 <- rep(TRUE, length(cI1))
      }
      
      # Post UMI linker test.
      if(! args$disablePostUmiLinker){
        pattern <- substr(x$adriftReadLinkerSeq, coords$post_n[1], coords$post_n[2])
        i3 <- vcountPattern(pattern, subseq(cR1@sread, coords$post_n[1], coords$post_n[2]), max.mismatch = args$postUmiLinkerMaxMismatch) > 0
      } else {
        i3 <- rep(TRUE, length(cI1))
      }
      
      d <- data.table()
      
      # Identify reads that pass all tests, i1 through i3.
      i <- which(i1 & i2 & i3)
      
      if(length(i) > 0){
        updateLog(paste0('<data chunk #', chunk_num, '>\t(row #', x$nowNum, ') Subsetting reads for sample row using ', ppNum(length(i)), ' indices.'), logFile = logFile)
        
        sI1 <- cI1[i]
        sR1 <- cR1[i]
        sR2 <- cR2[i]
        
        if(is.na(n_loc[1])) {
          linker1 <- rep(substr(x$adriftReadLinkerSeq, 1, coords$pre_n[2]), length(sI1))
          UMI     <- rep(NA_character_, length(sI1))
          linker2 <- rep(NA_character_, length(sI1))
        } else {
          linker1 <- as.character(subseq(sR1@sread, 1, coords$pre_n[2]))
          UMI     <- as.character(subseq(sR1@sread, coords$n_run[1], coords$n_run[2]))
          linker2 <- as.character(subseq(sR1@sread, coords$post_n[1], coords$post_n[2]))
        }
        
        sR1_qual <- as.character(sR1@quality@quality)
        sR2_qual <- as.character(sR2@quality@quality)
        
        # Poly-G trim.
        if(! args$disablePolyGfilter){
          sR1 <- stri_replace_first_regex(as.character(subseq(sR1@sread, coords$post_n[2]+1, width(sR1))), args$polyGfilterPattern, "")

          start_pos <- coords$post_n[2] + 1
          sR1_qual <- substr(sR1_qual, start_pos, start_pos + nchar(sR1) - 1)
          
          sR2 <- stri_replace_first_regex(as.character(sR2@sread), args$polyGfilterPattern, "")
          sR2_qual <- substr(sR2_qual, 1, nchar(sR2))
          
        } else {
          sR1 <- as.character(subseq(sR1@sread, coords$post_n[2]+1, width(sR1)))
          
          start_pos <- coords$post_n[2] + 1
          sR1_qual <- substr(sR1_qual, start_pos, start_pos + nchar(sR1) - 1)
        
          sR2 <- as.character(sR2@sread)
        }
        
        # Keep these in future version to export FASTQ.
        rm(sR1_qual, sR2_qual)
        
        # Build return table.
        d$readID          <- as.character(sI1@id)
        d$linker1         <- linker1
        d$UMI             <- UMI
        d$linker2         <- linker2
        d$adriftReadSeq   <- sR1
        d$anchorReadSeq   <- sR2
        d$vectorFastaFile <- as.factor(x$vectorFastaFile)
        d$leaderSeqHMM    <- as.factor(x$leaderSeqHMM)
        d$mode            <- as.factor(x$mode)
        d$trial           <- as.factor(x$trial)
        d$subject         <- as.factor(x$subject)
        d$sample          <- as.factor(x$sample)
        d$replicate       <- as.factor(x$replicate)
        d$refGenome       <- as.factor(x$refGenome)
        
        suppressWarnings(rm(sI1, sR1, toTrimIndex))
      }
      
      return(d)
    }))
    
    updateLog(paste0('<data chunk #', chunk_num, '>\t', ppNum(nrow(d)), ' reads demultiplexed.'), logFile = logFile)
    
    suppressWarnings(rm(sampleData, args, cI1, cR1, cR2))
    
    write_fst(d, file.path(args$tmpDir, paste0(chunk_num, '.fst')), compress = 0)
    return(TRUE)
  }
  
  param <- MulticoreParam(workers = args$threads, stop.on.error = TRUE)
  ### param <- SerialParam(stop.on.error = TRUE)
  
  updateLog(paste0('Demultiplexing data across ', args$threads, ' CPUs.'))
  updateLog(paste0('Starting asynchronous calculations. Data chunk logs can be found in ', args$logDir, '/'))
  
  x <- bpiterate(ITER = aav_iterator, 
                 FUN = aav_worker, 
                 BPPARAM = param,
                 sampleData = sampleData, 
                 args = args)
  
  updateLog('Demultiplexing calculations completed.')
  bpstop(param)
  closeAllConnections()
  
  logs <- unlist(lapply(list.files(args$logDir, pattern = '^log$', recursive = TRUE, full.names = TRUE), readLines))
  write(logs, args$logFile, append = TRUE)
  rm(logs)
  
  files <- list.files(args$tmpDir, pattern = '*.fst$', full.names = TRUE)
  
  updateLog(paste0('Collating ', length(files), ' data files from ', args$tmpDir, '/'))
  
  o <- rbindlist(lapply(files, function(x){
    updateLog(paste0('   Loading ', x, ' ...'))
    read_fst(x, as.data.table = TRUE)
  }), use.names = TRUE, fill = TRUE)
  
  
  if(! args$disableSequenceCollapse){
    updateLog('Collapsing duplicate reads.')
    o <- o[, .(
      nReads    = .N, 
      readID    = readID[1],
      linker1   = linker1[1],
      linker2   = linker2[1],
      mode      = mode[1],
      refGenome = refGenome[1],
      vectorFastaFile = vectorFastaFile[1],
      leaderSeqHMM    = leaderSeqHMM[1]
    ), by = .(trial, subject, sample, replicate, UMI, anchorReadSeq, adriftReadSeq)]
  } else {
    o <- o[, .(
      nReads    = .N, 
      readID    = readID[1],
      linker1   = linker1[1],
      linker2   = linker2[1],
      mode      = mode[1],
      refGenome = refGenome[1],
      vectorFastaFile = vectorFastaFile[1],
      leaderSeqHMM    = leaderSeqHMM[1]
    ), by = .(trial, subject, sample, replicate, readID)]
  }
  
  demux_summary <- o[, .(
    demultiplexedReads = sum(as.numeric(nReads), na.rm = TRUE)
  ), by = .(trial, subject, sample, replicate)]
  
  invisible(gc(verbose = FALSE))
  
  updateLog('Writing final data objects.')
  
  updateLog(paste0('Writing ', ppNum(nrow(data)), ' reads.'))
  saveRDS(o, file.path(args$outputDir, paste0(args$fileTag, '.rds')), compress = FALSE)
  
  # Summary table
  setDT(sampleData)
  sampleData[, original_order := .I]
  
  demux_summary[, replicate := as.numeric(as.character(replicate))]
  
  sampleData <- merge(
    sampleData, 
    demux_summary, 
    by = c("trial", "subject", "sample", "replicate"), 
    all.x = TRUE
  )
  
  setorder(sampleData, original_order)
  sampleData[, original_order := NULL]
  
  sampleData[is.na(demultiplexedReads), demultiplexedReads := 0]
  
  write_tsv(sampleData, file.path(args$outputDir, paste0(args$fileTag, '.tbl')))
  
  updateLog('Demultiplex module completed.')
  write(date(), file.path(args$outputDir, paste0(args$fileTag, '.done')))
}

#-------------------------------------------------------------------------------

args <- parser$parse_args()
args$qualTrimCode <- rawToChar(as.raw(args$qualTrimScore + 33))
source(file.path(args$softwareRoot, 'lib.R'))

tryCatch({
  runModule()
  q(status = 0)
}, error = function(e) {
  message("Caught error: ", e$message)
  q(status = 1)
})

#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'ShortRead', 'parallel', 'data.table', 'BiocParallel', 'stringi', 'fst')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",                    type = "character",     required = TRUE,                  help = "Directory for output files")
parser$add_argument("--sampleData",                   type = "character",     required = TRUE,                  help = "Sample definition file")
parser$add_argument("--anchorReads",                  type = "character",     required = TRUE,                  help = "Path to the Reverse read FASTQ file")
parser$add_argument("--softwareRoot",                 type = "character",     required = TRUE,                  help = "Path to AAVengeR installation.")
parser$add_argument("--threads",                      type = "integer",       default = 50,                     help = "Number of threads to use.")
parser$add_argument("--fileTag",                      type = "character",     default = "testHMMs",             help = "String appended to output files in the outpt directory.")
parser$add_argument("--ramDiskPath",                  type = "character",     default = "/dev/shm",             help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
parser$add_argument("--vectorDir",                    type = "character",     default = 'none',                 help = "Path to custom vector files.")
parser$add_argument("--hmmDir",                       type = "character",     default = 'none',                 help = "Path to custom hmm files.")
parser$add_argument("--HMMmatchEnd",                  action = "store_true",  default = FALSE,                  help = "Require a match to the end of the HMM (True / FALSE).")
parser$add_argument("--HMMmatchTerminalSeq",          type = "character",     default = 'none',                 help = "Path to custom hmm files.")
parser$add_argument("--HMMmatchEndRadius",            type = "integer",       default = 2,                      help = "Search radisuf for HMMmatchTerminalSeq")

runModule <- function(){
  startModule()
  
  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  on.exit({
    unlink(args$tmpDir,  recursive = TRUE, force = TRUE)
    unlink(args$logDir,  recursive = TRUE, force = TRUE) 
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  updateLog('Starting testHMMs module.')
  
  # Sanity checks...
  
  if(! file.exists(args$sampleData)){
    updateLog(paste0('Error - sampleData file does not exits. Provided path: ', args$sampleData)) 
    quit(status = 1)
  }
  
  vector_hmm_copy()
  
  sampleData <- read_tsv(args$sampleData, show_col_types = FALSE)
  
  knownHMMs <- list.files(file.path(args$softwareRoot, 'data', 'hmms'))
  
  if(! file.exists(args$anchorReads)){
    updateLog(paste0('Error - R2 sequencing file does not exits. Provided path: ', args$anchorReads)) 
    quit(status = 1)
  }
  
  requiredFields <- c('trial', 'subject', 'sample', 'replicate', 'adriftReadLinkerSeq', 'index1Seq', 'refGenome', 'vectorFastaFile', 'leaderSeqHMM', 'mode')
  
  if(! all(requiredFields %in% names(sampleData))){
    updateLog(paste0('Error - these required fields are missing from the sample data file: ', paste0(sQuote(requiredFields[! requiredFields %in% names(sampleData)]), collapse = ', ')))
    quit(status = 1)
  }
 
  if(! all(sampleData$leaderSeqHMM %in% knownHMMs)){
    missingHMMs <- paste0(unique(sampleData$leaderSeqHMM)[! unique(sampleData$leaderSeqHMM) %in% knownHMMs], collapse = ', ')
    updateLog(paste0('Error - These hmm file names in the sample data file were not found in ', file.path(args$softwareRoot, 'data', 'hmms'), ': ', missingHMMs)) 
    quit(status = 1)
  }
  
  R2 <- readFastq(args$anchorReads)
  
  browser()
  
  clean_ids <- BStringSet(gsub('^[^:]+:[^:]+:[^:]+:|\\s.+$', '', as.character(R2@id)))
  R2@id <- clean_ids
  
  writeFasta(R2, file.path(args$tmpDir, 'R2.fasta'))
  
  plots <- lapply(unique(sampleData$leaderSeqHMM), function(hmm){
    comm <- paste0('nhmmer --dna --F1 1 --F2 1 --F3 1 -T -5 --incT -5 --nobias --popen 0.15 --pextend 0.05 --cpu ', args$threads, ' --tblout ', 
                   file.path(args$tmpDir, 'out.tbl'), ' ', file.path(args$softwareRoot, 'data', 'hmms', hmm), ' ', 
                   file.path(args$tmpDir, 'R2.fasta'), ' > /dev/null')
    
    system(comm)
    o <- readr::read_table(file.path(args$tmpDir, 'out.tbl'), col_names = FALSE, col_types = NULL, comment = "#", show_col_types = FALSE)
    names(o) <- c('targetName', 'targetAcc', 'queryName', 'queryAcc', 'hmmStart', 'hmmEnd', 'targetStart', 'targetEnd', 'envStart', 'envEnd', 'seqLength', 'strand', 'fullEval', 'fullScore', 'bias', 'desc')
    
    invisible(file.remove(file.path(args$tmpDir, 'out.tbl')))
    
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
    
    o2 <- group_by(o, targetName) %>% slice_max(fullScore, with_ties = FALSE) %>% ungroup()
    
    # Read in the HMM so that we can test if we aligned to its end.
    h <- readLines(file.path(args$softwareRoot, 'data', 'hmms', hmm))
    hmmLength <- as.integer(unlist(strsplit(h[grepl('^LENG', h)], '\\s+'))[2])
    hmmName <- unlist(strsplit(h[grepl('^NAME', h)], '\\s+'))[2]
    
    if(args$HMMmatchEnd) o2 <- o2[abs(hmmLength - o2$hmmEnd) <= args$HMMmatchEndRadius,]
    
    data <- data.table(readID = sub('\\s.+$', '', as.character(R2@id)), anchorReadSeq = as.character(R2@sread))
    data <- data[data$readID %in% o2$targetName]
    
    data <- left_join(data, dplyr::select(o2, targetName, targetStart, targetEnd, fullScore), by = c('readID' = 'targetName'))
    
    if(! grepl('none', args$HMMmatchTerminalSeq, ignore.case = TRUE)){
      data$anchorReadSeq <- toupper(data$anchorReadSeq)
      args$HMMmatchTerminalSeq <- toupper(args$HMMmatchTerminalSeq)
      
      terminal_matchSeq <- substr(data$anchorReadSeq,  (data$targetEnd - (nchar(args$HMMmatchTerminalSeq) - 1) - args$HMMmatchEndRadius), (data$targetEnd + args$HMMmatchEndRadius))
      ends <- stringr::str_locate(terminal_matchSeq, args$HMMmatchTerminalSeq)[, 2]
      
      i <- ! is.na(ends)
      data <- data[i]
      ends <- ends[i]
      
      data$targetEnd <- data$targetEnd - (nchar(args$HMMmatchTerminalSeq) + args$HMMmatchEndRadius) + ends
    }
    
    lower_bound <- quantile(data$fullScore, 0.001)
    upper_bound <- quantile(data$fullScore, 0.999)
    data_clean <- data[fullScore >= lower_bound & fullScore <= upper_bound]

    p1 <- ggplot(data_clean, aes(x = "All Reads", y = fullScore)) + 
          geom_violin(fill = "blue4", color = NA, alpha = 0.8) +
          geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", alpha = 0.6) +
          stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "darkred", color = 'black') +
          scale_y_continuous(breaks = seq(0, max(data_clean$fullScore, na.rm = TRUE) + 5, by = 3)) +
          theme_bw() +
          labs(x = NULL, y = "HMM Score") +
          ggtitle(paste0(hmm, ' - HMM max scores per read')) +
          theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
          
    lower_bound <- quantile(data$targetStart, 0.001)
    upper_bound <- quantile(data$targetStart, 0.999)
    data_clean <- data[targetStart >= lower_bound & targetStart <= upper_bound]
    
    p2 <- ggplot(data_clean, aes(x = "All Reads", y = targetStart)) + 
      geom_violin(fill = "blue4", color = NA, alpha = 0.8) +
      geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", alpha = 0.6) +
      stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "darkred", color = 'black') +
      scale_y_continuous(breaks = seq(0, max(data_clean$targetStart, na.rm = TRUE) + 3, by = 3)) +
      theme_bw() +
      labs(x = NULL, y = "HMM read start position") +
      ggtitle(paste0(hmm, ' - HMM read alignment start position')) +
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
    
    list(p1, p2)
  })

  p <- unlist(plots, recursive = FALSE)
  pdf(file.path(args$outputDir, paste0(args$fileTag, '.pdf')), width = 8, height = 6)
  invisible(lapply(p, print))
  dev.off()
  
  updateLog('testHMMs module completed.')
  write(date(), file.path(args$outputDir, paste0(args$fileTag, '.done')))
}

args <- parser$parse_args()
source(file.path(args$softwareRoot, 'lib.R'))

tryCatch({
  runModule()
  q(status = 0)
}, error = function(e) {
  message("Caught error: ", e$message)
  q(status = 1)
})

#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'ShortRead', 'parallel', 'data.table')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",        type = "character", required = TRUE,       help = "Directory for output files.")
parser$add_argument("--inputDir",        type = "character", required = TRUE,       help = "Directory for output files.")
parser$add_argument("--softwareRoot",     type = "character", required = TRUE,       help = "Path to INSPIIRED2 installation.")
parser$add_argument("--indexReads",                   type = "character",     required = TRUE,                  help = "Path to the Index1 read FASTQ file")
parser$add_argument("--adriftReads",                  type = "character",     required = TRUE,                  help = "Path to the Forward read FASTQ file")
parser$add_argument("--anchorReads",                  type = "character",     required = TRUE,                  help = "Path to the Reverse read FASTQ file")
parser$add_argument("--fileTag",          type = "character", default  = "testHMMs", help = "String appended to output files in the output directory.")
parser$add_argument("--ramDiskPath",      type = "character", default  = "/dev/shm", help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
parser$add_argument("--topIndexRows", type = "integer", default = 50, help = "Number of top index rows to display")
parser$add_argument("--sampleLinkerLen", type = "integer", default = 20, help = "Length of adrift linker sequence to consider in reports.")

runModule <- function() {
  startModule()
  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  on.exit({
    unlink(args$tmpDir, recursive = TRUE, force = TRUE)
    unlink(args$logDir, recursive = TRUE, force = TRUE)
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  a <- bind_rows(lapply(list.files(args$inputDir, pattern = 'demultiplex_[^\\.]+\\.rds', full.names = TRUE), readRDS))
  tab1 <- group_by(a, subject, sample, refGenome) %>% summarise(reads = n_distinct(readID)) %>% ungroup()
  
  a <- bind_rows(lapply(list.files(args$inputDir, pattern = 'prepReads_[^\\.]+\\.rds', full.names = TRUE), readRDS))
  tab2 <- group_by(a, subject, sample, refGenome) %>% summarise(prepedReads = n_distinct(readID)) %>% ungroup()
  
  a <- bind_rows(lapply(list.files(args$inputDir, pattern = 'alignReads_[^\\.]+\\.rds', full.names = TRUE), function(x){ d <- readRDS(x); d[[1]]} ))
  tab3 <- group_by(a, subject, sample, refGenome) %>% summarise(alignments = n_distinct(readID)) %>% ungroup()
  
  a <- bind_rows(lapply(list.files(args$inputDir, pattern = 'buildFragments_[^\\.]+\\.rds', full.names = TRUE), readRDS))
  a$id <- paste0(a$fragChromosome, a$fragStrand, a$fragStart, '-', a$fragEnd)
  tab4 <- group_by(a, subject, sample, refGenome) %>% summarise(fragments = n_distinct(id)) %>% ungroup()
  
  f  <- list.files(args$inputDir, pattern = 'buildStdFragments_[^\\.]+\\.rds', full.names = TRUE) 
  f <- f[! grepl('cluster', f, ignore.case = TRUE)]
  f <- f[! grepl('multiHit', f, ignore.case = TRUE)]
  a <- bind_rows(lapply(f, readRDS))
  tab5 <- group_by(a, subject, sample, refGenome) %>% summarise(stdFragments = n()) %>% ungroup()
  
  f  <- list.files(args$inputDir, pattern = 'buildStdFragments_[^\\.]+\\.rds', full.names = TRUE) 
  f <- f[grepl('multiHitClusters', f, ignore.case = TRUE)] 
  a <- bind_rows(lapply(f, function(x){
         o <- readRDS(x)
         o$refGenome <- rev(unlist(str_split(x, '_')))[2]
         o
       }))
  tab6  <- group_by(a, subject, sample, refGenome) %>% summarise(multiHitClusts= n_distinct(clusterID)) %>% ungroup()
  
  a <- bind_rows(lapply(list.files(args$inputDir, pattern = 'buildSites_[^\\.]+\\.rds', full.names = TRUE), readRDS))
  tab7 <- group_by(a, subject, sample, refGenome) %>% summarise(sites = n_distinct(posid)) %>% ungroup()
  
  rm(a)
  
  attTab <- purrr::reduce(list(tab1, tab2, tab3, tab4, tab5, tab6, tab7), full_join, by = c("subject", "sample", "refGenome")) %>% arrange(refGenome)
  
  if(any(grepl('Yeast', attTab$subject, ignore.case = TRUE))) attTab[grepl('Yeast', attTab$subject, ignore.case = TRUE),]$subject <- 'YeastPosCtrl'
  
  if(any(is.na(attTab$subject))) attTab <- attTab[! is.na(attTab$subject),]
  
  system(paste0('Rscript ', file.path(args$softwareRoot, 'modules', 'buildSeqDataMap.R'),
         ' --outputDir ', args$outputDir,
         ' --inputData ', args$adriftReads,
         ' --softwareRoot ', args$softwareRoot,
         ' --fileTag adriftReadsMap'))

  system(paste0('Rscript ', file.path(args$softwareRoot, 'modules', 'buildSeqDataMap.R'),
         ' --outputDir ', args$outputDir,
         ' --inputData ', args$anchorReads,
         ' --softwareRoot ', args$softwareRoot,
         ' --fileTag anchorReadsMap'))

  
  indexReads  <- readFastq(args$indexReads)@sread
  adriftReads <- subseq(readFastq(args$adriftReads)@sread, 1, args$sampleLinkerLen)
  
  indexTab <- data.frame(sort(table(as.character(indexReads)), decreasing = TRUE)[1:args$topIndexRows])
  indexLinkerTab <- data.frame(sort(table(paste0(as.character(indexReads), ':', as.character(adriftReads))), decreasing = TRUE)[1:args$topIndexRows])
  
  seqsTab <- bind_cols(indexTab, indexLinkerTab)
  names(seqsTab) <- c('Barcode', '', 'Barcode:sample linker', '')
  
  sacCer_reads <- bind_rows(lapply(list.files(args$inputDir, pattern = 'demultiplex_[^\\.]+\\.rds', full.names = TRUE), readRDS)) %>% filter(grepl('sacCer', refGenome))
  sacCer_sites <- bind_rows(lapply(list.files(args$inputDir, pattern = 'buildSites_[^\\.]+\\.rds', full.names = TRUE), readRDS)) %>% filter(grepl('sacCer', refGenome)) %>% select(posid, reads, sonicLengths)

  rmarkdown::render(file.path(args$softwareRoot, 'modules', 'buildAnalysisReport.Rmd'),
                    output_dir = args$outputDir,
                    output_file = paste0(args$fileTag, '.pdf'),
                    params = list('date'  = format(Sys.Date(), format="%B %d, %Y"),
                                  'title' = '220124_MN01490_0059_A000H3M2VG'))
  
}

args <- parser$parse_args()

# Testing params...
# args <- list()
# args$outputDir <- 'INSPIIRED2_U5_only_with_yeast_output/'
# args$inputDir <- 'INSPIIRED2_U5_only_with_yeast_output'
# args$softwareRoot <- '/scratch/super1/everett/INSPIIRED2_dev/INSPIIRED2'
# args$indexReads <- 'Undetermined_S0_I1_001.fastq.gz'
# args$adriftReads <- 'Undetermined_S0_R1_001.fastq.gz'
# args$anchorReads <- 'Undetermined_S0_R2_001.fastq.gz'
# args$fileTag <- 'buildAnalysisReport'
# args$ramDiskPath <- '/dev/shm'
# args$topIndexRows <- 25
# args$sampleLinkerLen <- 20


args$outputDir <- normalizePath(args$outputDir, mustWork = TRUE)

source(file.path(args$softwareRoot, 'lib', 'common.R'))

tryCatch({
  runModule()
}, error = function(e) {
  cat("ERROR: ", conditionMessage(e), "\n", sep = "", file = stderr())
  flush(stderr())
  quit(save = "no", status = 1, runLast = FALSE)
})
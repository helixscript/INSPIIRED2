#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'ShortRead', 'parallel', 'data.table', 'BiocParallel', 'stringi', 'fst')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",               type = "character",     required = TRUE,            help = "Directory for output files")
parser$add_argument("--inputData",               type = "character",     required = TRUE,            help = "Path to demultiplex module's rds output file.")
parser$add_argument("--softwareRoot",            type = "character",     required = TRUE,            help = "Path to AAVengeR installation.")
parser$add_argument("--threads",                 type = "integer",       default = 50,               help = "Number of threads to use.")
parser$add_argument("--fileTag",                 type = "character",     default = "buildFragments", help = "String appended to output files in the outpt directory.")
parser$add_argument("--ramDiskPath",             type = "character",     default = "/dev/shm",       help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
parser$add_argument("--dataRowChunkSize",        type = "integer",       default = 5000,             help = "Numbers of data rows to process per alignment worker.")
parser$add_argument("--minFrgamentLength",       type = "integer",       default = 50,               help = "Min. Fragment length.")
parser$add_argument("--maxFrgamentLength",       type = "integer",       default = 100000L,           help = "Max. Fragment length.")

runModule <- function(){
  startModule()
  
  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  on.exit({
    unlink(args$tmpDir,  recursive = TRUE, force = TRUE)
    unlink(args$logDir,  recursive = TRUE, force = TRUE) 
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  updateLog('Starting buildFragment module.')
  
  o <- readRDS(args$inputData)
  setkey(o$anchorReads, readID)
  setkey(o$adriftReads, readID)
  
  readIDs <- unique(o$anchorReads$readID)
  readIDs <- split(readIDs, ceiling(seq_along(readIDs) / args$dataRowChunkSize))
  
  frags <- rbindlist(lapply(readIDs, function(x){
    a <- o$anchorReads[.(x), on = .(readID), nomatch = NULL]
    b <- o$adriftReads[.(x), on = .(readID), nomatch = NULL] 
    b <- b[, .(readID = readID, adrift_seq = seq, adrift_tName = tName, adrift_strand = strand, adrift_tStart = tStart, adrift_tEnd = tEnd)]
    r <- b[a, on = .(readID), nomatch = NULL, allow.cartesian = TRUE]
    rm(a, b)
    gc()
    
    r <- r[r$tName == r$adrift_tName]
    if(nrow(r) == 0) return(data.table())
    
    r <- r[r$strand != r$adrift_strand]
    if(nrow(r) == 0) return(data.table())
    
    r$fragStart  <- ifelse(r$strand == '+', r$tStart, r$adrift_tStart)
    r$fragEnd    <- ifelse(r$strand == '+', r$adrift_tEnd, r$tEnd)
    r$fragStrand <- ifelse(r$strand == '+', '+', '-')
    r$fragChromosome <- r$tName
    r$fragWidth = (r$fragEnd - r$fragStart) + 1
  
    r[, c('tName', 'tStart', 'tEnd', 'adrift_tName', 'adrift_strand', 'adrift_tStart', 'adrift_tEnd', 'strand') := NULL]
    
    r[fragWidth >= args$minFrgamentLength & fragWidth <= args$maxFrgamentLength]
  }))
  
  setnames(frags, "seq", "anchor_seq")
  
  frags$trial     <- as.character(frags$trial)
  frags$subject   <- as.character(frags$subject)
  frags$sample    <- as.character(frags$sample)
  frags$replicate <- as.integer(frags$replicate)
  frags$fragStart <- as.integer(frags$fragStart)
  frags$fragEnd   <- as.integer(frags$fragEnd)
  
  saveRDS(frags, file.path(args$outputDir, paste0(args$fileTag, '.rds')))
  updateLog('buildFragments module completed.')
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

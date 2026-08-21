#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'data.table', 'GenomicRanges')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",     type = "character",     required = TRUE,             help = "Directory for output files")
parser$add_argument("--inputData",     type = "character",     required = TRUE,             help = "Path to demultiplex module's rds output file.")
parser$add_argument("--softwareRoot",  type = "character",     required = TRUE,             help = "Path to INSPIIRED2 installation.")
parser$add_argument("--threads",       type = "integer",       default = 50,                help = "Number of threads to use.")
parser$add_argument("--fileTag",       type = "character",     default = "annotateRepeats", help = "String appended to output files in the outpt directory.")
parser$add_argument("--ramDiskPath",   type = "character",     default = "/dev/shm",        help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")

runModule <- function(){
  startModule()
  
  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  on.exit({
    unlink(args$tmpDir,  recursive = TRUE, force = TRUE)
    unlink(args$logDir,  recursive = TRUE, force = TRUE) 
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  updateLog('Starting annotateRepeats module.')
  
  if(! file.exists(args$inputData))  stop(paste0('Error - the input data file (', args$inputData, ') does not exist.'))
  if(file.size(args$inputData) == 0) stop(paste0('Error - the input data file (', args$inputData, ') is empty.'))
  
  d <- readRDS(args$inputData)

  o <- rbindlist(lapply(split(d, d$refGenome), function(x){ 
    tab <- read_tsv(file.path(args$softwareRoot, 'data', 'genomeAnnotations', paste0(x$refGenome[1], '.repeatTable.gz')), show_col_types = FALSE)
    
    tab$strand <- sub('C', '-', tab$strand)
    tab <- subset(tab, strand %in% c('+', '-'))
    
    tab <- makeGRangesFromDataFrame(dplyr::select(tab, query_seq, query_start, query_end, strand, repeat_name, repeat_class), 
                                    seqnames.field = 'query_seq', start.field = 'query_start', end.field = 'query_end', 
                                    strand.field = 'strand', keep.extra.columns = TRUE)
    
    i <- sub('\\.\\d+', '', unique(x$posid))
    
    g <- makeGRangesFromDataFrame(tibble(seqnames = unlist(lapply(str_split(i, '[\\+\\-]'), '[', 1)),
                                         start = unlist(lapply(str_split(i, '[\\+\\-]'), '[', 2)),
                                         end = start,
                                         strand = str_extract(i, '[\\+\\-]')))
    
    g$posid <- paste0(seqnames(g), strand(g), start(g))
    
    o <- data.frame(suppressWarnings(GenomicRanges::distanceToNearest(g, tab, select='all', ignore.strand=TRUE)))
    
    r <- unlist(GenomicRanges::GRangesList(lapply(1:length(g), function(xx){
           gg <- g[xx]
           oo <- unique(o[o$queryHits == xx,]) %>% dplyr::filter(distance == 0)
           
           gg$repeat_name  <- NA
           gg$repeat_class <- NA
           
           hits <- distinct(data.frame(tab[oo$subjectHits,])[, c('repeat_name', 'repeat_class')])
           
           gg$repeat_name <- paste0(hits$repeat_name, collapse = ',')
           gg$repeat_class <- paste0(hits$repeat_class, collapse = ',')
           gg
         })))
    
    as.data.table(data.frame(r))
  }))
  
  d <- dplyr::left_join(d, o[,6:length(o)], by = 'posid')
  d <- dplyr::relocate(d, repeat_name, .after = posid)
  d <- dplyr::relocate(d, repeat_class, .after = repeat_name)
  
  d[nchar(d$repeat_name) == 0,]$repeat_name   <- NA
  d[nchar(d$repeat_class) == 0,]$repeat_class <- NA

  saveRDS(d, file.path(args$outputDir, paste0(args$fileTag, '.rds')))
  write(date(), file.path(args$outputDir, paste0(args$fileTag, '.done')))
  updateLog('Completed annotateRepeats module.')
}

args <- parser$parse_args()
source(file.path(args$softwareRoot, 'lib', 'common.R'))

tryCatch({
  runModule()
  q(status = 0)
}, error = function(e) {
  message("Caught error: ", e$message)
  q(status = 1)
})

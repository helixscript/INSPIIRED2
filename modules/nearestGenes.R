#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'data.table', 'GenomicRanges')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",     type = "character",     required = TRUE,          help = "Directory for output files")
parser$add_argument("--inputData",     type = "character",     required = TRUE,          help = "Path to demultiplex module's rds output file.")
parser$add_argument("--softwareRoot",  type = "character",     required = TRUE,          help = "Path to INSPIIRED2 installation.")
parser$add_argument("--threads",       type = "integer",       default = 50,             help = "Number of threads to use.")
parser$add_argument("--fileTag",       type = "character",     default = "nearestGenes", help = "String appended to output files in the outpt directory.")
parser$add_argument("--ramDiskPath",   type = "character",     default = "/dev/shm",     help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")

runModule <- function(){
  startModule()
  
  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  on.exit({
    unlink(args$tmpDir,  recursive = TRUE, force = TRUE)
    unlink(args$logDir,  recursive = TRUE, force = TRUE) 
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  updateLog('Starting nearestGenes module.')
  
  if(! file.exists(args$inputData))  stop(paste0('Error - the input data file (', args$inputData, ') does not exist.'))
  if(file.size(args$inputData) == 0) stop(paste0('Error - the input data file (', args$inputData, ') is empty.'))
  
  d <- readRDS(args$inputData)

  o <- rbindlist(lapply(split(d, d$refGenome), function(x){
    tu <- readRDS(file.path(args$softwareRoot, 'data', 'genomeAnnotations', paste0(x$refGenome[1], '.TUs.rds')))
    ex <- readRDS(file.path(args$softwareRoot, 'data', 'genomeAnnotations', paste0(x$refGenome[1], '.exons.rds')))
    
    i <- sub('\\.\\d+', '', unique(x$posid))
    
    g <- makeGRangesFromDataFrame(tibble(seqnames = unlist(lapply(str_split(i, '[\\+\\-]'), '[', 1)),
                                         start = unlist(lapply(str_split(i, '[\\+\\-]'), '[', 2)),
                                         end = start,
                                         strand = str_extract(i, '[\\+\\-]')))
    
    g$posid <- paste0(seqnames(g), strand(g), start(g))
    
    o <- data.frame(suppressWarnings(GenomicRanges::distanceToNearest(g, tu, select='all', ignore.strand=TRUE)))
    e <- data.frame(suppressWarnings(GenomicRanges::distanceToNearest(g, ex, select='all', ignore.strand=TRUE)))
    
    r <- unlist(GenomicRanges::GRangesList(lapply(1:length(g), function(xx){
           gg <- g[xx]
           oo <- unique(o[o$queryHits == xx,])
           ee <- unique(e[e$queryHits == xx,])
           
           hits <- distinct(data.frame(tu[oo$subjectHits,])[, c('name2', 'strand')])
           
           gg$nearestGene <- paste0(hits$name2, collapse = ',')
           gg$nearestGeneDist <- min(oo$distance)
           gg$nearestGeneStrand <- paste0(hits$strand, collapse = ',')
           gg$beforeNearestGene <- start(gg) < min(start(tu[oo$subjectHits,]))
           gg$inExon <- min(ee$distance) == 0
           gg$inGene <- min(oo$distance) == 0
           gg
         })))
    
    as.data.table(data.frame(r))
  }))
  
  d <- dplyr::left_join(d, o[,6:length(o)], by = 'posid')
  d <- dplyr::relocate(d, inGene, .after = posid)
  d <- dplyr::relocate(d, inExon, .after = inGene)
  d <- dplyr::relocate(d, nearestGene, .after = inExon)
  d <- dplyr::relocate(d, nearestGeneDist, .after = nearestGene)
  d <- dplyr::relocate(d, nearestGeneStrand, .after = nearestGeneDist)

  saveRDS(d, file.path(args$outputDir, paste0(args$fileTag, '.rds')))
  write(date(), file.path(args$outputDir, paste0(args$fileTag, '.done')))
  updateLog('Completed nearestGenes module.')
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

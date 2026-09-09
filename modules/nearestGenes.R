#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'data.table', 'GenomicRanges')) suppressPackageStartupMessages(library(p, character.only = TRUE))
if(!requireNamespace("GenomeInfoDb", quietly = TRUE)) stop("Error - required R package GenomeInfoDb is not installed.")

parser <- ArgumentParser()
parser$add_argument("--outputDir",     type = "character", required = TRUE,          help = "Directory for output files")
parser$add_argument("--inputData",     type = "character", required = TRUE,          help = "Path to demultiplex module's rds output file.")
parser$add_argument("--softwareRoot",  type = "character", required = TRUE,          help = "Path to INSPIIRED2 installation.")
parser$add_argument("--threads",       type = "integer",   default = 50,             help = "Number of threads to use.")
parser$add_argument("--fileTag",       type = "character", default = "nearestGenes", help = "String appended to output files in the outpt directory.")
parser$add_argument("--ramDiskPath",   type = "character", default = "/dev/shm",     help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")

runModule <- function(){
  doneFile <- file.path(args$outputDir, paste0(args$fileTag, '.done'))
  if(file.exists(doneFile) && unlink(doneFile) != 0) stop('Error - could not remove stale completion marker: ', doneFile)
  
  startModule()

  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  on.exit({
    unlink(args$tmpDir, recursive = TRUE, force = TRUE)
    unlink(args$logDir, recursive = TRUE, force = TRUE)
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  updateLog('Starting nearestGenes module.')
  resource_overlay()
  
  if(!file.exists(args$inputData)) stop(paste0('Error - the input data file (', args$inputData, ') does not exist.'))
  if(file.size(args$inputData) == 0) stop(paste0('Error - the input data file (', args$inputData, ') is empty.'))
  
  d <- readRDS(args$inputData)
  
  o <- rbindlist(lapply(split(d, d$refGenome), function(x){
    tu <- readRDS(file.path(args$softwareRoot, 'data', 'genomeAnnotations', paste0(x$refGenome[1], '.TUs.rds')))
    ex <- readRDS(file.path(args$softwareRoot, 'data', 'genomeAnnotations', paste0(x$refGenome[1], '.exons.rds')))
    
    # i <- sub('\\.\\d+', '', unique(x$posid))
    i <- unique(x$posid)
    
    g <- makeGRangesFromDataFrame(tibble(
      seqnames = unlist(lapply(str_split(i, '[\\+\\-]'), '[', 1)),
      start = unlist(lapply(str_split(i, '[\\+\\-]'), '[', 2)),
      end = start,
      strand = str_extract(i, '[\\+\\-]')
    ))
    
    g$posid <- paste0(seqnames(g), strand(g), start(g))
    g$refGenome <- x$refGenome[1]
    
    o <- data.frame(suppressWarnings(GenomicRanges::distanceToNearest(g, tu, select = 'all', ignore.strand = TRUE)))
    e <- data.frame(suppressWarnings(GenomicRanges::distanceToNearest(g, ex, select = 'all', ignore.strand = TRUE)))
    
    r <- unlist(GenomicRanges::GRangesList(lapply(seq_along(g), function(xx){
      gg <- g[xx]
      oo <- unique(o[o$queryHits == xx, , drop = FALSE])
      
      gg$inGene <- IRanges::overlapsAny(gg, tu, ignore.strand = TRUE)
      gg$inExon <- IRanges::overlapsAny(gg, ex, ignore.strand = TRUE)
      gg$nearestGene <- NA_character_
      gg$nearestGeneStrand <- NA_character_
      gg$nearestGeneDist <- NA_integer_
      gg$beforeNearestGene <- NA
      
      if(nrow(oo) > 0L){
        hits <- distinct(data.frame(tu[oo$subjectHits, ])[, c("name2", "strand")])
        
        gg$nearestGene <- paste(hits$name2, collapse = ",")
        gg$nearestGeneStrand <- paste(hits$strand, collapse = ",")
        gg$nearestGeneDist <- as.integer(min(oo$distance) + as.integer(!gg$inGene))
        gg$beforeNearestGene <- start(gg) < min(start(tu[oo$subjectHits]))
      }
      
      gg
    })))
    
    as.data.table(data.frame(r))
  }))
  
  o <- o[,6:length(o)]
  if(anyDuplicated(o[, .(refGenome, posid)])) stop('Error - nearestGenes produced duplicate annotations for the same refGenome and posid.')
  
  inputRows <- nrow(d)
  d <- dplyr::left_join(d, o, by = c('refGenome', 'posid'))
  if(nrow(d) != inputRows) stop('Error - nearestGenes annotation join changed the input row count.')
  
  d <- dplyr::relocate(d, inGene, .after = posid)
  d <- dplyr::relocate(d, inExon, .after = inGene)
  d <- dplyr::relocate(d, nearestGene, .after = inExon)
  d <- dplyr::relocate(d, nearestGeneDist, .after = nearestGene)
  d <- dplyr::relocate(d, nearestGeneStrand, .after = nearestGeneDist)
  d <- dplyr::relocate(d, beforeNearestGene, .after = nearestGeneStrand)
  
  saveRDS(d, file.path(args$outputDir, paste0(args$fileTag, '.rds')))
  updateLog('Completed nearestGenes module.')
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
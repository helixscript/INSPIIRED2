#!/usr/bin/env Rscript                                                                                                                                                                                                                                                                                                                                                  

for (p in c('argparse', 'dplyr', 'tidyr', 'readr', 'GenomicRanges', 'rtracklayer', 'stringdist', 'stringr', 'ggplot2')) suppressPackageStartupMessages(library(p, character.only = TRUE))                                                                                                                                                                               
options(warn=-1)                                                                                                                                                                                                                                                                                                                                                        

# Parse command line arguments.                                                                                                                                                                                                                                                                                                                                         
parser <- ArgumentParser()                                                                                                                                                                                                                                                                                                                                              
parser$add_argument("-s", "--sitesFile",    type="character",  default='sites.rds',            help="Path to sites output file (rds).", metavar="")                                                                                                                                                                                                                     
parser$add_argument("-m", "--multiHitFile", type="character",  default='multiHitClusters.rds', help="Path to multiHitCluster output file (rds).", metavar="")                                                                                                                                                                                                           
parser$add_argument("-f", "--filterFile",   type="character",  default='anchorReadClustAttrition.rds', help="Path to INSPIIRED2 anchor read cluster attrition file (rds).", metavar="")                                                                                                                                                                                 
parser$add_argument("-t", "--truthFile",    type="character",  default='truth.tsv',            help="Path to synthetic data truth file (tsv).", metavar="")                                                                                                                                                                                                             
parser$add_argument("-b", "--twoBitPath",   type="character",  default='../../data/referenceGenomes/hg38.2bit', help="Path to local 2bit reference genome.", metavar="")                                                                                                                                                                                                
parser$add_argument("-o", "--outputDir",    type="character",  default='out',                  help="Path to output directory.", metavar="")                                                                                                                                                                                                                            
parser$add_argument("-w", "--siteWidth",    type="integer",    default=5,                      help="Number of NTs to expand truth positions during evaluation.", metavar="")                                                                                                                                                                                           
parser$add_argument("-v", "--version",      type="integer",    default=2,                      help="INSPIIRED pipeline version (1 or 2). Default: 2", metavar="")                                                                                                                                                                                                      
args <- parser$parse_args()             

# args$sitesFile <- '/home/everett/scratch/processing_runs/INSPIIRED2/U5_sites5000_seed1_error0.00/out/buildSites.rds'
# args$multiHitFile <- '/home/everett/scratch/processing_runs/INSPIIRED2/U5_sites5000_seed1_error0.00/out/buildStdFragments_multiHitClusters.rds'
# args$filterFile <- '/home/everett/scratch/processing_runs/INSPIIRED2/U5_sites5000_seed1_error0.00/out/buildStdFragments_anchorReadClusters.rds'
# args$truthFile <- '/home/everett/scratch/processing_runs/INSPIIRED2/U5_sites5000_seed1_error0.00/truth.tsv'
# args$outputDir <- '/home/everett/scratch/processing_runs/INSPIIRED2/U5_sites5000_seed1_error0.00/eval'

if(! file.exists(args$sitesFile))    { message('Error - the sites file could not be found.'); q(save = "no", status = 1, runLast = FALSE) }                                                                                                                                                                                                                             
if(! file.exists(args$truthFile))    { message('Error - the truth file could not be found.'); q(save = "no", status = 1, runLast = FALSE) }                                                                                                                                                                                                                             
if(! dir.exists(args$outputDir))     { dir.create(args$outputDir, recursive = TRUE) }                                                                                                                                                                                                                                                                                   

parse_posid <- function(df, posid_col = "posid") {                                                                                                                                                                                                                                                                                                                      
  df %>%                                                                                                                                                                                                                                                                                                                                                                
    dplyr::mutate(                                                                                                                                                                                                                                                                                                                                                      
      chr = stringr::str_extract(!!rlang::sym(posid_col), "^[^\\+\\-]+"),                                                                                                                                                                                                                                                                                               
      strand = stringr::str_extract(!!rlang::sym(posid_col), "[\\+\\-]"),                                                                                                                                                                                                                                                                                               
      pos = as.numeric(stringr::str_extract(!!rlang::sym(posid_col), "\\d+$"))                                                                                                                                                                                                                                                                                          
    )                                                                                                                                                                                                                                                                                                                                                                   
} 

t <- readr::read_tsv(args$truthFile, show_col_types = FALSE) %>% parse_posid()                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                                    
s <- readRDS(args$sitesFile) %>% parse_posid()
s$found <- FALSE

if(file.exists(args$filterFile)){
  f <- readRDS(args$filterFile) %>% filter(remove == TRUE) %>% parse_posid()
} else {
  f <- tibble()
}

m <- readRDS(args$multiHitFile) %>% tidyr::unnest(posids) %>% parse_posid(posid_col = 'posids')

r <- bind_rows(lapply(split(t, 1:nrow(t)), function(x){
  u <- s[s$trial == x$trial     & 
         s$subject == x$subject & 
         s$sample == x$sample   & 
         s$chr == x$chr         & 
         s$strand == x$strand   & 
         s$pos >= x$pos - args$siteWidth & 
         s$pos <= x$pos + args$siteWidth,]
  
  s[s$posid %in% u$posid,]$found <<- TRUE
  
  k <- m[m$trial == x$trial     & 
         m$subject == x$subject & 
         m$sample == x$sample   & 
         m$chr == x$chr         & 
         m$strand == x$strand   & 
         m$pos >= x$pos - args$siteWidth & 
         m$pos <= x$pos + args$siteWidth,]
  
  o <- f[f$trial == x$trial     & 
         f$subject == x$subject & 
         f$sample == x$sample   & 
         f$chr == x$chr         & 
         f$strand == x$strand   & 
         f$pos >= x$pos - args$siteWidth & 
         f$pos <= x$pos + args$siteWidth,]
  
  x$foundUnique   <- 0
  x$foundMulti    <- 0
  x$foundFiltered <- 0
  x$missing       <- 0
  
  if(nrow(u) > 0){
    x$foundUnique <- 1
  } else if (nrow(k) > 0){
    x$foundMulti <- 1
  } else if (nrow(o) > 0){
    x$foundFiltered <- 1
  } else {
    x$missing <- 1
  }
  
  x
})) %>% dplyr::select(-chr, -strand, -pos)

tab1 <- tibble(uniqueSites   = sprintf("%.2f%%", (sum(r$foundUnique)/nrow(r))*100),
               multiHitSites = sprintf("%.2f%%", (sum(r$foundMulti)/nrow(r))*100),
               filteredSites = sprintf("%.2f%%", (sum(r$foundFiltered)/nrow(r))*100),
               missingSites  = sprintf("%.2f%%", (sum(r$missing)/nrow(r))*100),
               totalRecovery = sprintf("%.2f%%", ((sum(r$foundUnique) + sum(r$foundMulti) + sum(r$foundFiltered))/nrow(r))*100),
               imaginedSites = sum(! s$found))

readr::write_tsv(tab1, file.path(args$outputDir, 'table1.tsv'))
readr::write_tsv(r, file.path(args$outputDir, 'table2.tsv'))
readr::write_tsv(s[s$found == FALSE,], file.path(args$outputDir, 'table3.tsv'))

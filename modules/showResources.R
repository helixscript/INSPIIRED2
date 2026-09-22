#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'data.table')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--softwareRoot", type = "character", required = TRUE,  help = "Path to INSPIIRED2 installation.")

runModule <- function(){
  resource_overlay(logging = FALSE)
  #tree_output <- system2(
    fs::dir_tree(file.path(args$softwareRoot, 'data'), recurse = TRUE)
  #)
  
  #writeLines(tree_output)
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
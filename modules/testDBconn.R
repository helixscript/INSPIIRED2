#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'RMariaDB')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--dbConfigFile", type = "character", default = 'none', help = "Path to db credential file.")
parser$add_argument("--dbConfigID",   type = "character", default = 'none', help = "DB credential block identifier in db credential file.")
parser$add_argument("--softwareRoot", type = "character", required = TRUE,  help = "Path to INSPIIRED2 installation.")

runModule <- function(){
  if(args$dbConfigFile == 'none'){
    stop('Error - the db config file was not provided with the --dbConfigFile flag.')  
  }
  
  if(args$dbConfigID == 'none'){
    stop('Error - the db config block id was not provided with the --dbConfigID flag.')  
  }
  
  if(! file.exists(args$dbConfigFile)){
    stop(paste0('Error - the db config file "', args$dbConfigFile, '" could not be found.'))
  }
  
  conn <-createDBconnection()
  message('Connection successful.')
  dbDisconnect(conn)
}

args <- parser$parse_args()
source(file.path(args$softwareRoot, 'lib', 'common.R'))

tryCatch({
  runModule()
}, error = function(e) {
  cat("ERROR: ", conditionMessage(e), "\n", sep = "", file = stderr())
  flush(stderr())
  quit(save = "no", status = 1, runLast = FALSE)
})
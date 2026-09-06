#!/usr/bin/env -S Rscript --vanilla
options(useFancyQuotes = FALSE)

for (p in c('argparse', 'data.table', 'RMariaDB')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--dbConfigFile",     type = "character", default = 'none', help = "Path to db credential file.")
parser$add_argument("--dbConfigID",       type = "character", default = 'none', help = "DB credential block identifier in db credential file.")
parser$add_argument("--outputFilePath",   type = "character", default = 'buildFragments.rds', help = "Path to export final rds file.")
parser$add_argument("--dataPath",         type = "character", default = 'none', help = "Path to INSPIIRED2 parquet file collection.")
parser$add_argument("--trials",           type = "character", default = 'none', help = "Comma delimited list of trial identifiers.")
parser$add_argument("--subjects",         type = "character", default = 'none', help = "Comma delimited list of subject identifiers.")
parser$add_argument("--samples",          type = "character", default = 'none', help = "Comma delimited list of sample identifiers.")
parser$add_argument("--refGenomes",       type = "character", default = 'none', help = "Comma delimited list of reference genome identifiers.")
parser$add_argument("--modes",            type = "character", default = 'none', help = "Comma delimited list of mode identifiers.")
parser$add_argument("--softwareRoot",     type = "character", required = TRUE,  help = "Path to INSPIIRED2 installation.")

runModule <- function(){
  if(args$dbConfigFile == 'none'){
    stop('Error - the db config file was not provided with the --dbConfigFile flag.')  
  }
  
  if(args$dbConfigID == 'none'){
    stop('Error - the db config block id was not provided with the --dbConfigID flag.')  
  }
  
  if(args$dataPath == 'none'){
    stop('Error - the data path containing fragment parquet files must be defined with the --dataPath flag.')  
  }
  
  if(! dir.exists(args$dataPath)){
    stop('Error - the provided data path does not exist.')  
  }
  
  if(args$trials == 'none'){
    stop('Error - one or more trial identifiers were not provided with the --trials flag.')  
  }
  
  if(! file.exists(args$dbConfigFile)){
    stop(paste0('Error - the db config file "', args$dbConfigFile, '" could not be found.'))
  }
  
  conn <- createDBconnection()
  on.exit(DBI::dbDisconnect(conn), add = TRUE)
  
  q1 <- 'x'
  if(args$trials != 'none'){
    q1 <- paste0('trial in (', paste(sQuote(unlist(strsplit(args$trials, '\\s*,\\s*'))), collapse = ', '), ')')
  }
  
  q2 <- 'x'
  if(args$subjects != 'none'){
    q2 <- paste0('subject in (', paste(sQuote(unlist(strsplit(args$subjects, '\\s*,\\s*'))), collapse = ', '), ')')
  }
  
  q3 <- 'x'
  if(args$samples != 'none'){
    q3 <- paste0('sample in (', paste(sQuote(unlist(strsplit(args$samples, '\\s*,\\s*'))), collapse = ', '), ')')
  }
  
  q4 <- 'x'
  if(args$refGenomes != 'none'){
    q4 <- paste0('ref_genome in (', paste(sQuote(unlist(strsplit(args$refGenomes, '\\s*,\\s*'))), collapse = ', '), ')')
  }
  
  q5 <- 'x'
  if(args$modes != 'none'){
    q5 <- paste0('mode in (', paste(sQuote(unlist(strsplit(args$modes, '\\s*,\\s*'))), collapse = ', '), ')')
  }
  
  q <- paste('select data_file_name from fragments where', gsub('AND x', '', paste(q1, 'AND', q2, 'AND', q3, 'AND', q4, 'AND', q5)))
  q <- sub('\\s*$', '', q)
  q <- gsub('\\s*AND\\s*', ' AND ', q)
  
  q <- paste(q, 'ORDER BY trial, subject, sample, replicate, ref_genome, mode')
  
  message(paste('Constructed query:', q))
  
  o <- dbGetQuery(conn, q)
  
  if (nrow(o) == 0L) {
    stop("No database records matched the requested filters.")
  }
  
  f <- file.path(args$dataPath, o$data_file_name)
  
  if(! all(file.exists(f))) stop('Error -- all of the expected parquest files could not be found.')
  
  d <- as.data.table(rbindlist(lapply(f, arrow::read_parquet)))
  
  # Rework the factor levels. 
  cols <- names(d)[vapply(d, is.factor, logical(1))]
  d[, (cols) := lapply(.SD, as.character), .SDcols = cols]
  d[, (cols) := lapply(.SD, as.factor), .SDcols = cols]
  
  saveRDS(d, args$outputFilePath, compress = FALSE)
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
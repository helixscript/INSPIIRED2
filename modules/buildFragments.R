#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'ShortRead', 'parallel', 'data.table', 'BiocParallel', 'stringi', 'RMariaDB')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",               type = "character",     required = TRUE,            help = "Directory for output files")
parser$add_argument("--inputData",               type = "character",     required = TRUE,            help = "Path to demultiplex module's rds output file.")
parser$add_argument("--softwareRoot",            type = "character",     required = TRUE,            help = "Path to INSPIIRED2 installation.")
parser$add_argument("--threads",                 type = "integer",       default = 50,               help = "Number of threads to use.")
parser$add_argument("--fileTag",                 type = "character",     default = "buildFragments", help = "String appended to output files in the outpt directory.")
parser$add_argument("--ramDiskPath",             type = "character",     default = "/dev/shm",       help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
parser$add_argument("--dataRowChunkSize",        type = "integer",       default = 5000L,            help = "Numbers of data rows to process per alignment worker.")
parser$add_argument("--minFrgamentLength",       type = "integer",       default = 50,               help = "Min. Fragment length.")
parser$add_argument("--maxFrgamentLength",       type = "integer",       default = 100000L,          help = "Max. Fragment length.")
parser$add_argument("--dbConfigFile",            type = "character",     default = 'none',           help = "Path to db credential file.")
parser$add_argument("--dbConfigID",              type = "character",     default = 'none',           help = "DB credential block identifier in db credential file.")

runModule <- function(){
  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  startModule()
  
  on.exit({
    unlink(args$tmpDir,  recursive = TRUE, force = TRUE)
    unlink(args$logDir,  recursive = TRUE, force = TRUE) 
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
    if(! is.null(args$dbConn)) dbDisconnect(args$dbConn)
  }, add = TRUE)
  
  updateLog('Starting buildFragment module.')
  
  if(! file.exists(args$inputData))  stop(paste0('Error - the input data file (', args$inputData, ') does not exist.'))
  if(file.size(args$inputData) == 0) stop(paste0('Error - the input data file (', args$inputData, ') is empty.'))
    
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
  
  frags$trial     <- as.factor(frags$trial)
  frags$subject   <- as.factor(frags$subject)
  frags$sample    <- as.factor(frags$sample)
  frags$replicate <- as.factor(frags$replicate)
  frags$fragStart <- as.integer(frags$fragStart)
  frags$fragEnd   <- as.integer(frags$fragEnd)
  
  saveRDS(frags, file.path(args$outputDir, paste0(args$fileTag, '.rds')))
  
  if(! is.null(args$dbConn)){
    frags[, fragID := paste(trial, subject, sample, replicate, fragChromosome, fragStrand, fragStart, fragEnd, sep = ":")]
    frags[, `:=`(g = .GRP, totalFrags = uniqueN(fragID)), by = .(trial, subject, sample, replicate, mode, refGenome)]
    
    for(x in split(frags, frags$g)){
      
      common_params <- list(
        trial   = as.character(x$trial[1]), 
        subject = as.character(x$subject[1]), 
        sample  = as.character(x$sample[1]), 
        rep     = as.integer(x$replicate[1]),
        genome  = as.character(x$refGenome[1]),
        mode    = as.character(x$mode[1])
      )
      
      record_tag <- paste(common_params, collapse = '|')
      updateLog(paste0('Processing data entry: ', record_tag))
      
      check_query <- "SELECT 1 FROM fragments WHERE trial = ?trial AND subject = ?subject AND sample = ?sample AND replicate = ?rep AND ref_genome = ?genome AND mode = ?mode LIMIT 1;"
      record_exists <- dbGetQuery(args$dbConn, DBI::sqlInterpolate(args$dbConn, check_query, .dots = common_params))
      
      if(nrow(record_exists) > 0){
        msg <- "Error - this entry is already in the database. In order to run the buildFragments module with databasing enabled, either remove the previously uploaded entry from the database or remove it from the module's input data object."
        updateLog(msg)
        stop(msg)
      }
      
      ts <- tmpString()
      arrow::write_parquet(x, file.path(args$outputDir, paste0(ts, '.parquet')))
      md5sum <- unname(tools::md5sum(file.path(args$outputDir, paste0(ts, '.parquet'))))
      file.rename(file.path(args$outputDir, paste0(ts, '.parquet')), file.path(args$outputDir, paste0(md5sum, '.parquet')))
      
      updateLog(paste0('Copying parquet file to data lake (', paste0(md5sum, '.parquet'), ').'))
      
      copyResult <- file.copy(file.path(args$outputDir, paste0(md5sum, '.parquet')),  file.path('/data', paste0(md5sum, '.parquet')), overwrite = TRUE)
      if(! copyResult) stop(paste0('Error - failed to copy ', file.path(args$outputDir, paste0(md5sum, '.parquet')), ' to ',  file.path('/data', paste0(md5sum, '.parquet'))))
      
      invisible(file.remove(file.path(args$outputDir, paste0(md5sum, '.parquet'))))
      
      updateLog('Inserting record into database.')
      
      insert_query <- "INSERT INTO fragments (trial, subject, sample, replicate, ref_genome, mode, total_fragments, data_file_name) 
                       VALUES (?trial, ?subject, ?sample, ?rep, ?genome, ?mode, ?total, ?file);"
      
      insert_params <- c(common_params, list(
        total = as.integer(x$totalFrags[1]),
        file  = paste0(md5sum, '.parquet')
      ))
      
      database_error <- NA
      insert_success <- tryCatch({
        rows_affected <- dbExecute(args$dbConn, DBI::sqlInterpolate(args$dbConn, insert_query, .dots = insert_params))
        rows_affected == 1
      }, error = function(cond) {
        message(paste("Database Error:", cond$message))
        database_error <<- cond$message
        return(FALSE)
      })
      
      if (! insert_success) stop(paste0('Database error caught: ', database_error, ' for record tag: ', record_tag))
      
      updateLog('Entry successfully processed.')
    }
  }
  
  updateLog('buildFragments module completed.')
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

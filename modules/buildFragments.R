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
parser$add_argument("--overwriteDBrecords",      action = "store_true",  default  = FALSE,           help = "Allow existing database records to be overwritten.") 

runModule <- function(){
  doneFile <- file.path(args$outputDir, paste0(args$fileTag, '.done'))
  if(file.exists(doneFile) && unlink(doneFile) != 0) stop('Error - could not remove stale completion marker: ', doneFile)
  
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
  
  updateLog(paste0(ppNum(length(readIDs)), ' aligned reads will be chunked into ', ppNum(args$dataRowChunkSize), ' read ID chunks.'))
  
  readIDs <- split(readIDs, ceiling(seq_along(readIDs) / args$dataRowChunkSize))
  
  updateLog(paste0('Aligned readIDs broken into ', ppNum(length(readIDs)), ' chunks for fragment generation.'))
  chunkNum <- 1
  
  frags <- rbindlist(lapply(readIDs, function(x){
    updateLog(paste0('Processing readID data chunk ',   ppNum(chunkNum), '/',  ppNum(length(readIDs))))
    chunkNum <<- chunkNum + 1
    
    a <- o$anchorReads[.(x), on = .(readID), nomatch = NULL]
    b <- o$adriftReads[.(x), on = .(readID), nomatch = NULL]
    
    b <- b[, .(readID = readID, adrift_seq = seq, adrift_tName = tName, adrift_strand = strand, adrift_tStart = tStart, adrift_tEnd = tEnd)]
    
    r <- b[a, on = .(readID), nomatch = NULL, allow.cartesian = TRUE]
    
    rm(a, b)
    ### gc()
    
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
  
  updateLog('All readID alignment data chunks processed.')
  
  if(nrow(frags) == 0){
    msg <- 'Error - no fragments could be built from the alignment data.'
    updateLog(msg)
    stop(msg)
  }
  
  setnames(frags, "seq", "anchor_seq")
  
  frags$trial     <- as.factor(frags$trial)
  frags$subject   <- as.factor(frags$subject)
  frags$sample    <- as.factor(frags$sample)
  frags$replicate <- as.factor(frags$replicate)
  frags$fragStart <- as.integer(frags$fragStart)
  frags$fragEnd   <- as.integer(frags$fragEnd)
  
  updateLog('Saving outputs.')
  
  saveRDS(frags, file.path(args$outputDir, paste0(args$fileTag, '.rds')))
  
  if(!is.null(args$dbConn)){
    updateLog("Database upload beginning.")
    
    data_lake <- "/data"
    if(!DBI::dbIsValid(args$dbConn)) stop("Error - database connection is not valid.")
    if(!dir.exists(data_lake) || file.access(data_lake, 2L) != 0L)
      stop("Error - data lake directory does not exist or is not writable: ", data_lake)
    
    engine_query <- paste(
      "SELECT ENGINE FROM information_schema.TABLES",
      "WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'fragments';"
    )
    engine <- DBI::dbGetQuery(args$dbConn, engine_query)
    engine_name <- if(nrow(engine) == 1L && ncol(engine) == 1L)
      as.character(engine[[1L]][1L]) else NA_character_
    
    if(is.na(engine_name) || toupper(engine_name) != "INNODB")
      stop("Error - the fragments table must use InnoDB for transactional database updates.")
    
    key_where <- paste(
      "trial = ?trial AND subject = ?subject AND sample = ?sample",
      "AND replicate = ?rep AND ref_genome = ?genome AND mode = ?mode"
    )
    check_query <- paste("SELECT data_file_name FROM fragments WHERE", key_where, ";")
    lock_query <- paste("SELECT data_file_name FROM fragments WHERE", key_where, "FOR UPDATE;")
    verify_query <- paste("SELECT data_file_name, total_fragments FROM fragments WHERE", key_where, ";")
    insert_query <- paste0(
      "INSERT INTO fragments (trial, subject, sample, replicate, ref_genome, mode, total_fragments, data_file_name) ",
      "VALUES (?trial, ?subject, ?sample, ?rep, ?genome, ?mode, ?total, ?file);"
    )
    update_query <- paste(
      "UPDATE fragments SET total_fragments = ?total, data_file_name = ?file,",
      "processed_date = CURRENT_TIMESTAMP WHERE", key_where, ";"
    )
    
    stage_parquet <- function(x, record_tag){
      local_path <- tempfile(pattern = "buildFragments_", tmpdir = args$tmpDir,
                             fileext = ".parquet")
      staged_path <- tempfile(pattern = ".inspiired2_", tmpdir = data_lake,
                              fileext = ".pending")
      
      on.exit(invisible(unlink(c(local_path, staged_path), force = TRUE)), add = TRUE)
      
      arrow::write_parquet(x, local_path)
      local_size <- if(file.exists(local_path)) file.size(local_path) else NA_real_
      
      if(is.na(local_size) || local_size < 1)
        stop("Error - failed to create a non-empty parquet file for record tag: ",
             record_tag)
      
      md5 <- unname(tools::md5sum(local_path))
      if(length(md5) != 1L || is.na(md5) ||
         !grepl("^[[:xdigit:]]{32}$", md5))
        stop("Error - failed to calculate the parquet checksum for record tag: ",
             record_tag)
      
      md5 <- tolower(md5)
      new_file <- paste0(md5, ".parquet")
      new_path <- file.path(data_lake, new_file)
      
      # Content-addressed file already exists. Validate rather than overwrite it.
      if(file.exists(new_path)){
        existing_md5 <- unname(tools::md5sum(new_path))
        if(length(existing_md5) != 1L || is.na(existing_md5) ||
           tolower(existing_md5) != md5)
          stop("Error - existing data-lake file does not match its checksum-derived name: ",
               new_path)
        
        return(new_file)
      }
      
      # Copy to a unique temporary path within /data and verify it before rename.
      copied <- file.copy(local_path, staged_path, overwrite = FALSE)
      staged_md5 <- if(isTRUE(copied) && file.exists(staged_path))
        unname(tools::md5sum(staged_path)) else NA_character_
      
      if(!isTRUE(copied) || length(staged_md5) != 1L ||
         is.na(staged_md5) || tolower(staged_md5) != md5)
        stop("Error - failed to copy and verify the staged parquet file for record tag: ",
             record_tag)
      
      # Rename inside /data so the final filename never exposes a partial copy.
      if(!file.rename(staged_path, new_path)){
        final_md5 <- if(file.exists(new_path))
          unname(tools::md5sum(new_path)) else NA_character_
        
        # Another process may have installed the identical content concurrently.
        if(length(final_md5) != 1L || is.na(final_md5) ||
           tolower(final_md5) != md5)
          stop("Error - failed to install the parquet file in the data lake: ",
               new_path)
      }
      
      final_md5 <- if(file.exists(new_path))
        unname(tools::md5sum(new_path)) else NA_character_
      
      if(length(final_md5) != 1L || is.na(final_md5) ||
         tolower(final_md5) != md5)
        stop("Error - final parquet verification failed for record tag: ",
             record_tag)
      
      new_file
    }
    
    split_cols <- c("trial", "subject", "sample", "replicate", "mode",
                    "refGenome")
    frag_groups <- split(frags, by = split_cols, keep.by = TRUE,
                         flatten = TRUE, sorted = TRUE, drop = TRUE)
    uploads <- vector("list", length(frag_groups))
    
    # Check the complete batch before writing replacement files.
    for(i in seq_along(frag_groups)){
      x <- frag_groups[[i]]
      common_params <- list(
        trial = as.character(x$trial[1]),
        subject = as.character(x$subject[1]),
        sample = as.character(x$sample[1]),
        rep = as.integer(as.character(x$replicate[1])),
        genome = as.character(x$refGenome[1]),
        mode = as.character(x$mode[1])
      )
      
      record_tag <- paste(unlist(common_params, use.names = FALSE),
                          collapse = "|")
      existing <- DBI::dbGetQuery(
        args$dbConn,
        DBI::sqlInterpolate(args$dbConn, check_query,
                            .dots = common_params)
      )
      
      if(nrow(existing) > 1L)
        stop("Error - multiple database records found for record tag: ",
             record_tag)
      
      if(nrow(existing) == 1L && !isTRUE(args$overwriteDBrecords))
        stop(
          "Error - database entry already exists; remove it, exclude it ",
          "from the input, or enable overwrite: ", record_tag
        )
      
      # Equivalent to the former unique fragID count because the remaining
      # fragID fields are constant within this database group.
      total_frags <- data.table::uniqueN(
        x, by = c("fragChromosome", "fragStrand", "fragStart", "fragEnd")
      )
      
      uploads[[i]] <- list(
        common = common_params,
        tag = record_tag,
        total = as.integer(total_frags),
        new_file = NA_character_
      )
    }
    
    # Make every replacement parquet durable before changing the database.
    for(i in seq_along(frag_groups)){
      updateLog(paste0("Preparing data entry: ", uploads[[i]]$tag))
      uploads[[i]]$new_file <- stage_parquet(
        frag_groups[[i]], uploads[[i]]$tag
      )
      updateLog(paste0(
        "Parquet ready in data lake (", uploads[[i]]$new_file, ")."
      ))
    }
    
    rm(frag_groups)
    
    prepared_files <- vapply(uploads, `[[`, character(1), "new_file")
    if(anyNA(prepared_files) ||
       !all(file.exists(file.path(data_lake, prepared_files))))
      stop(
        "Error - one or more prepared parquet files disappeared before ",
        "the database transaction."
      )
    
    # Commit the batch together so ordinary SQL errors roll back every row.
    old_files <- tryCatch({
      DBI::dbWithTransaction(args$dbConn, {
        previous_files <- rep(NA_character_, length(uploads))
        
        for(i in seq_along(uploads)){
          u <- uploads[[i]]
          current <- DBI::dbGetQuery(
            args$dbConn,
            DBI::sqlInterpolate(args$dbConn, lock_query,
                                .dots = u$common)
          )
          
          if(nrow(current) > 1L)
            stop("Error - multiple database records found for record tag: ",
                 u$tag)
          
          if(nrow(current) == 1L && !isTRUE(args$overwriteDBrecords))
            stop(
              "Error - entry was added concurrently and overwrite is disabled: ",
              u$tag
            )
          
          previous_files[i] <- if(nrow(current) == 1L)
            as.character(current$data_file_name[1]) else NA_character_
          
          params <- c(
            u$common,
            list(total = u$total, file = u$new_file)
          )
          
          if(nrow(current) == 1L){
            rows <- DBI::dbExecute(
              args$dbConn,
              DBI::sqlInterpolate(args$dbConn, update_query,
                                  .dots = params)
            )
            
            # MariaDB can report zero when an UPDATE does not change values.
            if(length(rows) != 1L || is.na(rows) ||
               !(rows %in% c(0L, 1L)))
              stop("Error - unexpected number of rows updated for record tag: ",
                   u$tag)
          } else {
            rows <- DBI::dbExecute(
              args$dbConn,
              DBI::sqlInterpolate(args$dbConn, insert_query,
                                  .dots = params)
            )
            
            if(length(rows) != 1L || is.na(rows) || rows != 1L)
              stop("Error - failed to insert database record for record tag: ",
                   u$tag)
          }
          
          stored <- DBI::dbGetQuery(
            args$dbConn,
            DBI::sqlInterpolate(args$dbConn, verify_query,
                                .dots = u$common)
          )
          
          if(nrow(stored) != 1L ||
             is.na(stored$data_file_name[1]) ||
             as.character(stored$data_file_name[1]) != u$new_file ||
             is.na(stored$total_fragments[1]) ||
             as.integer(stored$total_fragments[1]) != u$total)
            stop("Error - database verification failed for record tag: ",
                 u$tag)
        }
        
        previous_files
      })
    }, error = function(e){
      updateLog(
        "Database transaction failed; database state should be verified. ",
        "Prepared parquet files were retained."
      )
      stop(e)
    })
    
    for(i in seq_along(uploads)){
      updateLog(paste0(
        "Entry successfully processed: ", uploads[[i]]$tag
      ))
      
      if(!is.na(old_files[i]) && nzchar(old_files[i]) &&
         old_files[i] != uploads[[i]]$new_file)
        updateLog(paste0(
          "Previous parquet retained for safe garbage collection: ",
          old_files[i]
        ))
    }
  }
  
  updateLog('buildFragments module completed.')
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
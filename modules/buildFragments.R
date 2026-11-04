#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'ShortRead', 'parallel', 'data.table', 'BiocParallel', 'stringi', 'RMariaDB')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",               type = "character",     required = TRUE,            help = "Directory for output files.")
parser$add_argument("--inputData",               type = "character",     required = TRUE,            help = "Path to alignReads module's rds output file.")
parser$add_argument("--softwareRoot",            type = "character",     required = TRUE,            help = "Path to INSPIIRED2 installation.")
parser$add_argument("--threads",                 type = "integer",       default = 50,               help = "Number of threads to use.")
parser$add_argument("--fileTag",                 type = "character",     default = "buildFragments", help = "String appended to output files in the outpt directory.")
parser$add_argument("--ramDiskPath",             type = "character",     default = "/dev/shm",       help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
parser$add_argument("--dataRowChunkSize",        type = "integer",       default = 5000L,            help = "Numbers of data rows to process per alignment worker.")
parser$add_argument("--minFrgamentLength",       type = "integer",       default = 40,               help = "Min. Fragment length.")
parser$add_argument("--maxFrgamentLength",       type = "integer",       default = 100000L,          help = "Max. Fragment length.")
parser$add_argument("--dbConfigFile",            type = "character",     default = 'none',           help = "Path to db credential file.")
parser$add_argument("--dbConfigID",              type = "character",     default = 'none',           help = "DB credential block identifier in db credential file.")
parser$add_argument("--overwriteDBrecords",      action = "store_true",  default  = FALSE,           help = "Allow existing database records to be overwritten.") 

runModule <- function(){
  doneFile <- file.path(args$outputDir, paste0(args$fileTag, '.done'))
  if(file.exists(doneFile) && unlink(doneFile) != 0) stop('Error - could not remove stale completion marker: ', doneFile)
  
  startModule()
  
  yaml::write_yaml(args[names(args) != "dbConn"], file.path(args$outputDir, paste0(args$fileTag, ".yml")))
  
  on.exit({
    unlink(args$tmpDir,  recursive = TRUE, force = TRUE)
    unlink(args$logDir,  recursive = TRUE, force = TRUE) 
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
    if(! is.null(args$dbConn)) dbDisconnect(args$dbConn)
  }, add = TRUE)
  
  if(xor(args$dbConfigFile != "none", args$dbConfigID != "none")){
    msg <- "Error - supply both --dbConfigFile and --dbConfigID, or neither."
    updateLog(msg)
    stop(msg, call. = FALSE)
  }
  
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
  
  # Keep this module's database implementation local to the template runModule.
  uploadFragmentsToDB <- function(frags){
    key_columns <- c("trial", "subject", "sample", "replicate", "refGenome", "mode")
    db_columns <- c("trial", "subject", "sample", "replicate", "ref_genome", "mode")
    coordinate_columns <- c("fragChromosome", "fragStrand", "fragStart", "fragEnd")
    if(!is.data.frame(frags) || !nrow(frags) || anyDuplicated(names(frags)) ||
       !all(c(key_columns, coordinate_columns) %in% names(frags)))
      stop("Error - fragment uploads require a non-empty table with complete fragment keys and coordinates.", call. = FALSE)
    data <- data.table::as.data.table(data.table::copy(frags))
    keys <- as.data.frame(data)[key_columns]
    for(column in key_columns){
      values <- keys[[column]]
      if(!is.atomic(values) || !is.null(dim(values)) || is.complex(values) || anyNA(values))
        stop("Error - invalid fragment identifier column: ", column, call. = FALSE)
      values <- as.character(values)
      if(any(!nzchar(trimws(values))) || any(values != trimws(values)))
        stop("Error - blank or whitespace-padded fragment identifiers: ", column, call. = FALSE)
      keys[[column]] <- values
    }
    reps <- suppressWarnings(as.numeric(keys$replicate))
    if(any(!is.finite(reps) | reps != trunc(reps) |
           reps > .Machine$integer.max | reps < -.Machine$integer.max))
      stop("Error - fragment replicate identifiers must be valid integers.", call. = FALSE)
    if(!is.logical(args$overwriteDBrecords) || length(args$overwriteDBrecords) != 1L ||
       is.na(args$overwriteDBrecords))
      stop("Error - overwriteDBrecords must be TRUE or FALSE.", call. = FALSE)
    
    data_lake_root <- "/data"
    if(!dir.exists(data_lake_root) || !file.exists(file.path(data_lake_root, ".inspiired")) ||
       file.access(data_lake_root, 2L) != 0L)
      stop("Error - data lake must be writable and contain its .inspiired marker: ", data_lake_root, call. = FALSE)
    data_lake <- normalizePath(data_lake_root, mustWork = TRUE)
    tmp_dir <- args$tmpDir
    if(is.null(tmp_dir) || !dir.exists(tmp_dir) || file.access(tmp_dir, 2L) != 0L)
      stop("Error - upload temporary directory is missing or not writable.", call. = FALSE)
    if(!requireNamespace("arrow", quietly = TRUE))
      stop("Error - the arrow R package is required to upload fragments.", call. = FALSE)
    conn <- args$dbConn
    if(is.null(conn) || !DBI::dbIsValid(conn))
      stop("Error - database connection is not valid.", call. = FALSE)
    
    locked <- transaction_open <- verified <- FALSE
    published_files <- character()
    on.exit({
      if(transaction_open) tryCatch(DBI::dbRollback(conn), error = function(e) {
        warning("Fragment upload rollback failed: ", conditionMessage(e), call. = FALSE)
      })
      # Finalized checksum files may already be referenced by another upload.
      # Keep them on failure; staging removes only its local/pending files.
      if(!verified && length(published_files))
        updateLog(paste0("Prepared Parquet files retained after an unsuccessful or unverified fragment upload: ",
                         paste(published_files, collapse = ", ")))
      if(locked) tryCatch(DBI::dbGetQuery(conn,
                                          "SELECT RELEASE_LOCK(CONCAT('INSPIIRED2.uploadSites:', MD5(DATABASE()))) AS released"
      ), error = function(e) warning("Upload lock release failed: ", conditionMessage(e), call. = FALSE))
    }, add = TRUE)
    
    # Share the existing sites lock: both uploaders may reference the same lake
    # file. Hold it from preflight through staging, commit, readback and cleanup.
    lock <- DBI::dbGetQuery(conn,
                            "SELECT GET_LOCK(CONCAT('INSPIIRED2.uploadSites:', MD5(DATABASE())), 30) AS acquired")
    if(nrow(lock) != 1L || !isTRUE(lock$acquired[[1L]] == 1L))
      stop("Error - could not acquire the fragment/sites upload lock; retry after other uploads finish.", call. = FALSE)
    locked <- TRUE
    engine <- DBI::dbGetQuery(conn, paste(
      "SELECT ENGINE AS engine FROM information_schema.TABLES",
      "WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'fragments'"))
    if(nrow(engine) != 1L || !isTRUE(toupper(engine$engine[[1L]]) == "INNODB"))
      stop("Error - the fragments table must use InnoDB for transactional database updates.", call. = FALSE)
    schema <- DBI::dbGetQuery(conn, paste(
      "SELECT COLUMN_NAME AS column_name, CHARACTER_MAXIMUM_LENGTH AS max_length",
      "FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'fragments'"))
    if(!all(c(db_columns, "total_fragments", "processed_date", "data_file_name") %in% schema$column_name))
      stop("Error - the fragments table is missing required columns.", call. = FALSE)
    primary <- DBI::dbGetQuery(conn, paste(
      "SELECT COLUMN_NAME AS column_name FROM information_schema.KEY_COLUMN_USAGE",
      "WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'fragments'",
      "AND CONSTRAINT_NAME = 'PRIMARY' ORDER BY ORDINAL_POSITION"))
    if(!identical(as.character(primary$column_name), db_columns))
      stop("Error - fragments must have the six-column trial/subject/sample/replicate/ref_genome/mode primary key.", call. = FALSE)
    lengths <- setNames(as.numeric(schema$max_length), schema$column_name)
    for(i in which(key_columns != "replicate")){
      limit <- lengths[[db_columns[[i]]]]
      if(is.na(limit) || any(nchar(keys[[i]], type = "chars") > limit))
        stop("Error - fragments.", db_columns[[i]], " is too short for the supplied identifiers.", call. = FALSE)
    }
    if(is.na(lengths[["data_file_name"]]) || lengths[["data_file_name"]] < 40L)
      stop("Error - fragments.data_file_name must allow at least 40 characters.", call. = FALSE)
    
    key_where <- paste(paste(db_columns, "= ?"), collapse = " AND ")
    select_sql <- paste("SELECT data_file_name, total_fragments FROM fragments WHERE", key_where)
    delete_sql <- paste("DELETE FROM fragments WHERE", key_where)
    insert_sql <- paste(
      "INSERT INTO fragments (trial, subject, sample, replicate, ref_genome, mode, total_fragments, processed_date, data_file_name)",
      "VALUES (?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP, ?)"
    )
    is_symlink <- function(path){
      target <- Sys.readlink(path)
      !is.na(target) && nzchar(target)
    }
    old_path <- function(name){
      if(is.na(name) || !nzchar(name)) return(NULL)
      path <- if(startsWith(name, "/")) name else file.path(data_lake, name)
      if(is_symlink(path))
        stop("Error - refusing to remove an unsafe fragment file path: ", name, call. = FALSE)
      path <- normalizePath(path, mustWork = FALSE)
      if(dirname(path) != data_lake || dir.exists(path) ||
         (file.exists(path) && (!utils::file_test("-f", path) || !grepl("[.]parquet$", path))))
        stop("Error - refusing to remove an unsafe fragment file path: ", name, call. = FALSE)
      path
    }
    same_bytes <- function(left, right){
      if(!isTRUE(file.info(left)$size == file.info(right)$size)) return(FALSE)
      a <- file(left, open = "rb")
      on.exit(close(a), add = TRUE)
      b <- file(right, open = "rb")
      on.exit(close(b), add = TRUE)
      repeat {
        x <- readBin(a, what = "raw", n = 1048576L)
        y <- readBin(b, what = "raw", n = 1048576L)
        if(!identical(x, y)) return(FALSE)
        if(!length(x)) return(TRUE)
      }
    }
    stage_file <- function(x){
      local <- tempfile("buildFragments_", tmpdir = tmp_dir, fileext = ".parquet")
      pending <- tempfile(".fragments_", tmpdir = data_lake, fileext = ".pending")
      on.exit(unlink(c(local, pending)), add = TRUE)
      arrow::write_parquet(x, local)
      checksum <- if(file.exists(local) && file.size(local) > 0) unname(tools::md5sum(local)) else NA_character_
      if(length(checksum) != 1L || is.na(checksum) || !grepl("^[[:xdigit:]]{32}$", checksum))
        stop("Error - failed to write and checksum fragments Parquet file.", call. = FALSE)
      checksum <- tolower(checksum)
      filename <- paste0(checksum, ".parquet")
      path <- file.path(data_lake, filename)
      verify_file <- function(candidate){
        if(is_symlink(candidate) || !utils::file_test("-f", candidate))
          stop("Error - fragments Parquet path is not a regular, non-symlink file: ", candidate, call. = FALSE)
        if(!identical(unname(tools::md5sum(candidate)), checksum))
          stop("Error - fragments Parquet file does not match its checksum-derived name: ", candidate, call. = FALSE)
        if(!same_bytes(local, candidate))
          stop("Error - MD5 collision: different Parquet bytes have the same filename: ", candidate, call. = FALSE)
      }
      if(file.exists(path) || is_symlink(path)){
        verify_file(path)
        return(filename)
      }
      if(!isTRUE(file.copy(local, pending, overwrite = FALSE)) ||
         !identical(unname(tools::md5sum(pending)), checksum))
        stop("Error - failed to copy and verify fragments Parquet in the data lake.", call. = FALSE)
      # As in the sites uploader, publish atomically without overwriting a name.
      if(isTRUE(suppressWarnings(file.link(pending, path)))){
        published_files <<- c(published_files, path)
      } else if(!file.exists(path) && !is_symlink(path)){
        stop("Error - could not publish fragments Parquet without overwriting; the data-lake filesystem must support hard links: ",
             path, call. = FALSE)
      }
      verify_file(path)
      filename
    }
    
    # Preserve the existing grouping and the Parquet table's columns/types/order.
    split_cols <- c("trial", "subject", "sample", "replicate", "mode", "refGenome")
    frag_groups <- split(data, by = split_cols, keep.by = TRUE,
                         flatten = TRUE, sorted = TRUE, drop = TRUE)
    uploads <- vector("list", length(frag_groups))
    # Preflight the complete batch before publishing any replacement files.
    for(i in seq_along(frag_groups)){
      x <- frag_groups[[i]]
      params <- unname(lapply(key_columns, function(column) as.character(x[[column]][1L])))
      params[[4L]] <- as.integer(params[[4L]])
      tag <- paste(unlist(params, use.names = FALSE), collapse = "|")
      existing <- DBI::dbGetQuery(conn, select_sql, params = params)
      if(nrow(existing) > 1L)
        stop("Error - multiple database records found for record tag: ", tag, call. = FALSE)
      if(nrow(existing) == 1L && !args$overwriteDBrecords)
        stop("Error - database entry already exists; exclude it from the input or enable --overwriteDBrecords: ",
             tag, call. = FALSE)
      if(nrow(existing)) old_path(as.character(existing$data_file_name[[1L]]))
      total <- data.table::uniqueN(x, by = coordinate_columns)
      if(total > .Machine$integer.max)
        stop("Error - fragment count exceeds the database INT range for record tag: ", tag, call. = FALSE)
      uploads[[i]] <- list(params = params, tag = tag, total = as.integer(total), file = NA_character_)
    }
    for(i in seq_along(frag_groups)){
      uploads[[i]]$file <- stage_file(frag_groups[[i]])
      updateLog(paste0("Prepared fragments Parquet: ", uploads[[i]]$tag, " -> ", uploads[[i]]$file))
    }
    rm(frag_groups)
    new_names <- vapply(uploads, `[[`, character(1), "file")
    if(anyNA(new_names) || !all(file.exists(file.path(data_lake, new_names))))
      stop("Error - prepared fragments Parquet files disappeared before the database transaction.", call. = FALSE)
    verify_records <- function(connection){
      for(upload in uploads){
        stored <- DBI::dbGetQuery(connection, select_sql, params = upload$params)
        if(nrow(stored) != 1L ||
           !isTRUE(as.character(stored$data_file_name[[1L]]) == upload$file) ||
           !isTRUE(as.numeric(stored$total_fragments[[1L]]) == upload$total))
          stop("Error - fragments database verification failed for ", upload$tag, call. = FALSE)
      }
      invisible(TRUE)
    }
    
    DBI::dbBegin(conn)
    transaction_open <- TRUE
    old_files <- character()
    old_records <- lapply(uploads, function(upload){
      old <- DBI::dbGetQuery(conn, paste(select_sql, "FOR UPDATE"), params = upload$params)
      if(nrow(old) > 1L)
        stop("Error - multiple database records found for record tag: ", upload$tag, call. = FALSE)
      if(nrow(old) && !args$overwriteDBrecords)
        stop("Error - entry was added concurrently and overwrite is disabled: ", upload$tag, call. = FALSE)
      if(nrow(old)) old_path(as.character(old$data_file_name[[1L]]))
      old
    })
    # Delete ALL exact-key replacements before inserting ANY new rows. If input
    # keys compare equal in SQL, the primary key rejects the second INSERT and
    # rolls back the entire batch instead of silently replacing the first group.
    for(i in seq_along(uploads)){
      upload <- uploads[[i]]
      old <- old_records[[i]]
      if(nrow(old)) old_files <- c(old_files, as.character(old$data_file_name))
      deleted <- DBI::dbExecute(conn, delete_sql, params = upload$params)
      # Colliding input keys can already have deleted this same old record.
      # Requiring the exact count here also fails safely before any INSERT.
      if(!isTRUE(deleted == nrow(old)))
        stop("Error - unexpected fragments deletion count (possible duplicate database keys): ", upload$tag, call. = FALSE)
    }
    for(upload in uploads){
      inserted <- DBI::dbExecute(conn, insert_sql,
                                 params = c(upload$params, list(upload$total, upload$file)))
      if(!isTRUE(inserted == 1L))
        stop("Error - failed to insert a fragments record: ", upload$tag, call. = FALSE)
    }
    verify_records(conn)
    tryCatch({
      # Use the same checked COMMIT and DBI bookkeeping as the sites uploader.
      DBI::dbExecute(conn, "COMMIT", immediate = TRUE)
      DBI::dbCommit(conn)
    }, error = function(e) {
      stop("Error - fragments COMMIT failed; its outcome may be uncertain. Old and new Parquet files were retained. ",
           conditionMessage(e), call. = FALSE)
    })
    transaction_open <- FALSE
    tryCatch({
      verify <- createDBconnection()
      tryCatch({
        if(!DBI::dbIsValid(verify)) stop("Verification connection is not valid.")
        verify_records(verify)
      }, finally = DBI::dbDisconnect(verify))
    }, error = function(e) {
      stop("Error - fragments COMMIT completed but independent verification failed. Old and new Parquet files were retained. ",
           conditionMessage(e), call. = FALSE)
    })
    verified <- TRUE
    
    # Only the displaced files from this batch are eligible; never sweep /data.
    for(name in unique(old_files[!is.na(old_files) & nzchar(old_files)])){
      path <- old_path(name)
      if(!file.exists(path)) next
      if(basename(path) %in% new_names) next
      aliases <- unique(c(name, basename(path), path, file.path(data_lake_root, basename(path))))
      placeholders <- paste(rep("?", length(aliases)), collapse = ", ")
      referenced <- FALSE
      for(table in c("sites", "fragments")){
        refs <- tryCatch(DBI::dbGetQuery(conn, paste0("SELECT data_file_name FROM ", table,
                                                      " WHERE data_file_name IN (", placeholders, ") LIMIT 1"), params = as.list(aliases)),
                         error = function(e) stop("Error - fragments upload committed, but old-file reference checking failed. Retained ",
                                                  path, ": ", conditionMessage(e), call. = FALSE))
        referenced <- referenced || nrow(refs) > 0L
      }
      if(referenced){
        updateLog(paste0("Retained Parquet still referenced by another record: ", path))
      } else if(unlink(path) != 0L || file.exists(path)){
        stop("Error - fragments upload committed, but could not remove obsolete Parquet file: ", path, call. = FALSE)
      } else {
        updateLog(paste0("Removed obsolete fragments Parquet: ", path))
      }
    }
    updateLog(paste0("Fragments upload committed and verified: ", length(uploads), " record(s)."))
    invisible(TRUE)
  }
  
  if(!is.null(args$dbConn)){
    updateLog("Database upload beginning.")
    uploadFragmentsToDB(frags)
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
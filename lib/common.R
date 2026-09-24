tmpString <- function() paste0(Sys.getpid(), '___', stringi::stri_rand_strings(1, 30, '[A-Za-z0-9]'))


ppNum <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)


parse_bool <- function(x) {
  if (tolower(x) %in% c("true", "t", "1", "yes")) return(TRUE)
  if (tolower(x) %in% c("false", "f", "0", "no")) return(FALSE)
  stop("Must be a valid boolean value.")
}


timeElapsedString <- function(){
  elapsed_period <- lubridate::as.period(lubridate::now() - args$startTime)
  
  sprintf("%dh %02dmin %02.0fsec", 
          lubridate::hour(elapsed_period),
          lubridate::minute(elapsed_period), 
          lubridate::second(elapsed_period))
}


updateLog <- function(msg, logFile = NULL){
  if(is.null(logFile)) logFile <- args$logFile
  msg <- paste0(base::format(Sys.time(), "%m.%d.%Y"), ' [', timeElapsedString(), "]\t", msg)
  write(msg, file = logFile, append = TRUE)
}


createDBconnection <- function(){
  tryCatch({
    message(paste0('Connecting to the INSPIIRED database using cnf file "',  args$dbConfigFile, '" and group ID "', args$dbConfigID, '".'))
    dbConnect(RMariaDB::MariaDB(), group = args$dbConfigID, default.file = args$dbConfigFile)
  },
  error=function(cond) {
    stop(paste0('Error - could not connect to the database. Caught error: ', cond$message))
  })
}


commandStatus <- function(x){
  status <- if(is.character(x)) attr(x, "status", exact = TRUE) else x
  if(is.null(status)) 0L else suppressWarnings(as.integer(status))
}


requireCommandSuccess <- function(x, label){
  status <- commandStatus(x)
  if(length(status) != 1L || is.na(status) || status != 0L)
    stop("Error - ", label, " failed with exit status ",
         if(length(status) == 1L && !is.na(status)) status else "unknown", ".",
         call. = FALSE)
  invisible(x)
}


startModule <- function(){
  options(useFancyQuotes = FALSE)
  if(! file.access(args$ramDiskPath, mode = 2) == 0) args$ramDiskPath <- args$outputDir
  
  args$command      <<- paste(commandArgs(trailingOnly = FALSE), collapse = " ")
  args$version      <<- paste('INSPIIRED2', readLines(file.path(args$softwareRoot, 'VERSION')))
  args$logFile      <<- file.path(args$outputDir, paste0(args$fileTag, '.log'))
  args$tmpDir       <<- file.path(args$outputDir, paste0(args$fileTag, '_tmp'))
  args$logDir       <<- file.path(args$outputDir, paste0(args$fileTag, '_log'))
  args$defaultDelim <<- '___'
  args$ramDisk      <<- file.path(args$ramDiskPath,  "INSPIIRED2", gsub("\\.", "", paste0(format(Sys.time(), "%Y%m%d_%H%M%OS6"), "_", Sys.getpid())))
  args$startTime    <<- lubridate::now()
  
  if (isNamespaceLoaded("data.table")) data.table::setDTthreads(args$threads)
  
  if(! dir.exists(args$outputDir)) dir.create(args$outputDir, recursive = TRUE)
  if(! dir.exists(args$outputDir)) stop('Error -- output directory could not be created.')
  
  if(! dir.exists(args$logDir)) dir.create(args$logDir, recursive = TRUE)
  if(! dir.exists(args$logDir)) stop('Error -- log directory could not be created.')
  
  if(! dir.exists(args$tmpDir)) dir.create(args$tmpDir, recursive = TRUE)
  if(! dir.exists(args$tmpDir)) stop('Error -- tmp directory could not be created.')
  
  if(! dir.exists(args$ramDisk)) dir.create(args$ramDisk, recursive = TRUE)
  if(! dir.exists(args$ramDisk)) stop('Error -- could not create ram disk.')
  
  if(all(c('dbConfigFile', 'dbConfigID') %in% names(args))){
    if(args$dbConfigFile != 'none' & args$dbConfigID != 'none'){
      args$dbConn <<- createDBconnection()
      updateLog('Established conection with database.')
      
      # Database function is tied to a mounted external file system.
      # Test that this file system is mounted and if writable. 
      if(! file.exists('/data/.inspiired')) 
        stop('Error - an external file system with a file name .inspiired was expected to be mounted to /data in the Docker container.')
      
      probe <- '/data/.inspiired_writeTest'
      writable <- tryCatch({
        writeBin(as.raw(1L), probe)
        isTRUE(file.size(probe) == 1L)
      }, error = function(e) FALSE)
      
      removed <- !file.exists(probe) || unlink(probe) == 0L
      if(!writable || !removed)
        stop('Error - external file system is not writable')
    } else {
      args$dbConn <<- NULL
    }
  }
}



resource_overlay <- function(logging = TRUE){
  resourceRoot <- '/resources'
  dataRoot <- file.path(args$softwareRoot, 'data')
  
  if(!dir.exists(resourceRoot)) return(invisible(NULL))
  
  files <- list.files(resourceRoot, recursive=TRUE, full.names=TRUE, all.files=TRUE, no..=TRUE)
  files <- files[file.exists(files) & !file.info(files)$isdir]
  
  if(!length(files)) return(invisible(NULL))
  
  resourceRoot <- normalizePath(resourceRoot)
  rel <- substring(files, nchar(resourceRoot) + 2L)
  destinations <- file.path(dataRoot, rel)
  
  nLinked <- 0L
  
  for(i in seq_along(files)){
    source <- normalizePath(files[i])
    destination <- destinations[i]
    
    dir.create(dirname(destination), recursive=TRUE, showWarnings=FALSE)
    
    currentLink <- Sys.readlink(destination)
    
    # Already correctly overlaid.
    if(nzchar(currentLink)){
      currentTarget <- tryCatch(normalizePath(currentLink, mustWork=FALSE), error=function(e) currentLink)
      if(identical(currentTarget, source)) next
    }
    
    # Remove bundled file, old symlink, or other existing destination.
    if(file.exists(destination) || nzchar(currentLink)) unlink(destination)
    
    if(!file.symlink(source, destination)){
      stop('Could not overlay resource: ', source, ' -> ', destination)
    } else {
      if(logging) updateLog(paste0('Link created: ', source, ' -> ', destination))
    }
    
    nLinked <- nLinked + 1L
  }
  
  if(nLinked)
    if(logging) updateLog(paste0('Applied ', nLinked, ' user resource overlay(s) from ', resourceRoot))
  
  invisible(NULL)
}


make_dt_iterator <- function(dt, chunk_size, chunk_num_start = 0) {
  current_row <- 1
  total_rows <- nrow(dt)
  chunk_num <- chunk_num_start
  
  function() {
    if (current_row > total_rows) return(NULL)
    
    end_row <- min(current_row + chunk_size - 1, total_rows)
    chunk_data <- dt[current_row:end_row, ]
    chunk_num <<- chunk_num + 1
    current_row <<- end_row + 1
    
    list(
      data = chunk_data,
      chunk_num = chunk_num,
      is_last = (current_row > total_rows)
    )
  }
}


parse_cdhit_clstr <- function(file_path) {
  if(! file.exists(file_path)) stop()
  data.frame(raw_text = readLines(file_path), stringsAsFactors = FALSE) %>%
    dplyr::mutate(is_header = stringr::str_starts(raw_text, ">"),
                  cluster_id = ifelse(is_header, str_remove(raw_text, ">"), NA)) %>%
    tidyr::fill(cluster_id, .direction = "down") %>%
    dplyr::filter(!is_header) %>%
    dplyr::mutate(is_rep = str_detect(raw_text, "\\*\\s*$"),
                  readID = str_extract(raw_text, "(?<=>).+?(?=\\.\\.\\.)")) %>%
    dplyr::add_count(cluster_id, name = "cluster_size") %>%
    dplyr::select(readID, cluster_id, is_rep, cluster_size) %>%
    dplyr::arrange(desc(cluster_size), cluster_id, desc(is_rep))
}


run_blastn_parallel <- function(fastaFile, dbPath, params, threads = 60){
  if(length(threads) != 1L || !is.numeric(threads) || is.na(threads) ||
     !is.finite(threads) || threads < 1 || threads %% 1 != 0)
    stop("Error - threads must be a positive integer.")
  
  threads <- as.integer(threads)
  seqs <- Biostrings::readDNAStringSet(fastaFile)
  if(length(seqs) == 0L) return(data.table())
  
  num_chunks <- min(length(seqs), threads)
  
  # Avoid cut(..., breaks = 1) and unnecessary parallel overhead.
  if(num_chunks == 1L)
    return(run_blastn(fastaFile, dbPath, params, threads = 1L))
  
  chunks <- split(seqs, cut(seq_along(seqs), breaks = num_chunks, labels = FALSE))
  param <- BiocParallel::MulticoreParam(workers = num_chunks, tasks = num_chunks, stop.on.error = TRUE)
  
  results_list <- tryCatch(
    BiocParallel::bplapply(seq_along(chunks), function(i){
      tmp_chunk <- paste0(fastaFile, ".chunk_", i)
      on.exit(unlink(tmp_chunk, force = TRUE), add = TRUE)
      Biostrings::writeXStringSet(chunks[[i]], tmp_chunk)
      run_blastn(tmp_chunk, dbPath, params, threads = 1L)
    }, BPPARAM = param),
    finally = try(BiocParallel::bpstop(param), silent = TRUE)
  )
  
  data.table::rbindlist(results_list, use.names = TRUE, fill = TRUE)
}


run_blastn <- function(fastaFile, dbPath, params, threads = 1){
  
  outFile <- paste0(fastaFile, '.blastn')
  comm <- paste0("blastn ",  params,
                " -query ", fastaFile, 
                " -db ",    dbPath,
                " -out ",   outFile,
                " -num_threads ", threads,
                " -outfmt '6 qseqid qstart qend sstart send sstrand length pident gaps gapopen bitscore'")
  
  status <- system(comm)
  requireCommandSuccess(status, paste0("blastn query ", basename(fastaFile)))
  if(!file.exists(outFile)) stop("Error - blastn returned success without creating ", outFile, ".")
  
  hits <- data.table()
  if(file.info(outFile)$size > 0) hits <- fread(outFile, col.names = c("qName", "qS", "qE", "vS", "vE", "strand", "len", "pident", "gaps", "gapsopen", "bitscore"))
  invisible(file.remove(outFile))
  hits
}


parseBLAToutput <- function(f){
  if(!file.exists(f) || file.info(f)$size == 0) return(tibble::tibble())
  
  b <- readr::read_delim(f, delim='\t', col_names=FALSE, col_types='iiiiiiiicciiiciiiiccc')
  names(b) <- c('matches','misMatches','repMatches','nCount','qNumInsert','qBaseInsert','tNumInsert','tBaseInsert',
                'strand','qName','qSize','qStart','qEnd','tName','tSize','tStart','tEnd','blockCount',
                'blockSizes','qStarts','tStarts')
  
  scoreOutput <- suppressWarnings(system(paste(shQuote(file.path(args$softwareRoot, "bin", "pslScore.pl")), shQuote(f)), intern = TRUE))
  requireCommandSuccess(scoreOutput, "pslScore.pl")
  x <- read.table(textConnection(scoreOutput), sep = "\t")
  
  names(x) <- c('tName','tStart','tEnd','hit','pslScore','percentIdentity')
  
  if(nrow(x) != nrow(b)) stop('pslScore.pl output does not match PSL record count')
  
  b$queryPercentID <- as.numeric(x$percentIdentity)
  b$pslScore <- as.numeric(x$pslScore)
  
  b$qStart <- as.integer(b$qStart + 1L)
  b$qEnd   <- as.integer(b$qEnd)
  b$tStart <- as.integer(b$tStart + 1L)
  b$tEnd   <- as.integer(b$tEnd)
  
  b$qWidth <- as.integer(b$qEnd - b$qStart + 1L)
  b$tWidth <- as.integer(b$tEnd - b$tStart + 1L)
  
  dplyr::select(b, qName, matches, strand, qSize, qStart, qEnd, tName, tNumInsert, qNumInsert,
                tBaseInsert, qBaseInsert, tStart, tEnd, queryPercentID, pslScore, qWidth, tWidth)
}



# Same mounted data lake used by buildFragments. Kept separate from DB paths,
# which are relative filenames so the lake can be mounted elsewhere for reading.
.sitesDBDataLake <- "/data"

uploadSitesToDB <- function(sites){
  if(!is.list(args) || !all(c("dbConfigFile", "dbConfigID") %in% names(args)) ||
     any(vapply(args[c("dbConfigFile", "dbConfigID")], function(x)
       !is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x) || x == "none", logical(1))) ||
     !file.exists(args$dbConfigFile))
    stop("Error - uploadSitesToDB requires args$dbConfigID and an existing args$dbConfigFile.", call. = FALSE)
  key_columns <- c("trial", "subject", "sample", "refGenome", "mode")
  db_columns <- c("trial", "subject", "sample", "ref_genome", "mode")
  if(!is.data.frame(sites) || !nrow(sites))
    stop("Error - uploadSitesToDB requires a non-empty sites data frame.", call. = FALSE)
  if(anyDuplicated(names(sites)) || !all(c(key_columns, "posid") %in% names(sites)))
    stop("Error - sites must have unique column names and include trial, subject, sample, refGenome, mode and posid.", call. = FALSE)
  data <- as.data.frame(sites)
  for(column in c(key_columns, "posid")){
    values <- data[[column]]
    # read_tsv can infer numeric/date sample metadata. Convert SQL keys only;
    # preserve the original column types in the uploaded frame.
    valid_type <- is.atomic(values) && is.null(dim(values))
    if(column == "posid") valid_type <- valid_type && (is.character(values) || is.factor(values))
    if(!valid_type || anyNA(values) || any(!nzchar(trimws(as.character(values)))))
      stop("Error - sites column ", column, " must contain non-missing, non-empty scalar identifiers (strings for posid).", call. = FALSE)
  }
  keys <- as.data.frame(lapply(data[key_columns], as.character), stringsAsFactors = FALSE)
  groups <- dplyr::group_rows(dplyr::group_by(keys, dplyr::across(dplyr::all_of(key_columns))))
  if(!requireNamespace("arrow", quietly = TRUE))
    stop("Error - the arrow R package is required to upload sites.", call. = FALSE)
  if(!dir.exists(.sitesDBDataLake) || !file.exists(file.path(.sitesDBDataLake, ".inspiired")) ||
     file.access(.sitesDBDataLake, 2L) != 0L)
    stop("Error - data lake must be writable and contain its .inspiired marker: ", .sitesDBDataLake, call. = FALSE)
  data_lake <- normalizePath(.sitesDBDataLake, mustWork = TRUE)
  tmp_dir <- if(is.null(args$tmpDir)) tempdir() else args$tmpDir
  if(!dir.exists(tmp_dir) || file.access(tmp_dir, 2L) != 0L)
    stop("Error - upload temporary directory is missing or not writable: ", tmp_dir, call. = FALSE)
  log_message <- function(text){
    if(is.null(args$logFile)) message(text) else updateLog(text)
  }
  
  own_connection <- is.null(args$dbConn)
  conn <- if(own_connection) createDBconnection() else args$dbConn
  if(own_connection) on.exit(tryCatch(DBI::dbDisconnect(conn), error = function(e) warning(conditionMessage(e))), add = TRUE)
  if(!DBI::dbIsValid(conn)) stop("Error - database connection is not valid.", call. = FALSE)
  locked <- transaction_open <- commit_attempted <- FALSE
  created_files <- character()
  on.exit({
    if(transaction_open) tryCatch(DBI::dbRollback(conn), error = function(e) {
      warning("Upload rollback failed: ", conditionMessage(e), call. = FALSE)
    })
    # Before COMMIT, these UUID-named files cannot be referenced by a committed
    # upload. Once COMMIT is attempted its outcome may be uncertain; retain them.
    if(!commit_attempted && length(created_files)){
      if(unlink(created_files) != 0L || any(file.exists(created_files)))
        warning("Uncommitted sites files could not all be removed: ", paste(created_files, collapse = ", "), call. = FALSE)
    }
    if(locked) tryCatch(DBI::dbGetQuery(conn,
                                        "SELECT RELEASE_LOCK(CONCAT('INSPIIRED2.uploadSites:', MD5(DATABASE()))) AS released"
    ), error = function(e) warning("Upload lock release failed: ", conditionMessage(e), call. = FALSE))
  }, add = TRUE, after = FALSE)
  
  lock <- DBI::dbGetQuery(conn,
                          "SELECT GET_LOCK(CONCAT('INSPIIRED2.uploadSites:', MD5(DATABASE())), 30) AS acquired")
  if(nrow(lock) != 1L || !isTRUE(lock$acquired[[1L]] == 1L))
    stop("Error - could not acquire the sites upload lock; check the selected database and retry after other uploads finish.", call. = FALSE)
  locked <- TRUE
  engine <- DBI::dbGetQuery(conn, paste(
    "SELECT ENGINE AS engine FROM information_schema.TABLES",
    "WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'sites'"))
  if(nrow(engine) != 1L || !isTRUE(toupper(engine$engine[[1L]]) == "INNODB"))
    stop("Error - sites must be an InnoDB table.", call. = FALSE)
  schema <- DBI::dbGetQuery(conn, paste(
    "SELECT COLUMN_NAME AS column_name, CHARACTER_MAXIMUM_LENGTH AS max_length",
    "FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'sites'"))
  required <- c(db_columns, "total_sites", "processed_date", "data_file_name")
  if(!all(required %in% schema$column_name))
    stop("Error - the sites table is missing required columns.", call. = FALSE)
  lengths <- setNames(as.numeric(schema$max_length), schema$column_name)
  for(i in seq_along(key_columns)){
    limit <- lengths[[db_columns[[i]]]]
    if(is.na(limit) || any(nchar(keys[[i]], type = "chars") > limit)){
      hint <- if(db_columns[[i]] == "mode")
        " Run migrations/sites_mode_width.sql to support 'dual detect'." else ""
      stop("Error - sites.", db_columns[[i]], " is too short for the supplied identifiers.", hint, call. = FALSE)
    }
  }
  if(is.na(lengths[["data_file_name"]]) || lengths[["data_file_name"]] < 46L)
    stop("Error - sites.data_file_name must allow at least 46 characters.", call. = FALSE)
  
  key_where <- paste(paste(db_columns, "= ?"), collapse = " AND ")
  select_sql <- paste("SELECT data_file_name, total_sites FROM sites WHERE", key_where)
  insert_sql <- paste(
    "INSERT INTO sites (trial, subject, sample, ref_genome, mode, total_sites, processed_date, data_file_name)",
    "VALUES (?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP, ?)"
  )
  uploads <- vector("list", length(groups))
  is_symlink <- function(path){
    target <- Sys.readlink(path)
    !is.na(target) && nzchar(target)
  }
  stage_file <- function(x, path){
    local <- tempfile("sites_", tmpdir = tmp_dir, fileext = ".parquet")
    pending <- tempfile(".sites_", tmpdir = data_lake, fileext = ".pending")
    on.exit(unlink(c(local, pending)), add = TRUE)
    arrow::write_parquet(x, local)
    checksum <- if(file.exists(local) && file.size(local) > 0) unname(tools::md5sum(local)) else NA_character_
    if(length(checksum) != 1L || is.na(checksum)) stop("Error - failed to write sites Parquet file.", call. = FALSE)
    if(!isTRUE(file.copy(local, pending, overwrite = FALSE)) ||
       !identical(unname(tools::md5sum(pending)), checksum))
      stop("Error - failed to copy and verify sites Parquet file in the data lake.", call. = FALSE)
    if(!file.rename(pending, path) || !identical(unname(tools::md5sum(path)), checksum))
      stop("Error - failed to install and verify sites Parquet file: ", path, call. = FALSE)
  }
  
  # Prepare every file before deleting any database records. UUID names prevent
  # a failed replacement from overwriting the old file or another upload's file.
  for(i in seq_along(groups)){
    indices <- groups[[i]]
    token <- DBI::dbGetQuery(conn, "SELECT REPLACE(UUID(), '-', '') AS id")$id
    if(length(token) != 1L || is.na(token) || !grepl("^[[:xdigit:]]{32}$", token))
      stop("Error - could not allocate a unique sites Parquet filename.", call. = FALSE)
    filename <- paste0("sites_", tolower(token), ".parquet")
    path <- file.path(data_lake, filename)
    if(file.exists(path) || is_symlink(path))
      stop("Error - new sites Parquet path already exists: ", path, call. = FALSE)
    created_files <- c(created_files, path)
    stage_file(data[indices, , drop = FALSE], path)
    uploads[[i]] <- list(params = unname(as.list(keys[indices[[1L]], , drop = FALSE])),
                         total = length(unique(data$posid[indices])), file = filename)
    log_message(paste0("Prepared sites Parquet: ", filename))
  }
  new_names <- vapply(uploads, `[[`, character(1), "file")
  old_files <- character()
  old_path <- function(name){
    if(is.na(name) || !nzchar(name)) return(NULL)
    path <- if(startsWith(name, "/")) name else file.path(data_lake, name)
    if(is_symlink(path)) stop("Error - refusing to remove an unsafe sites file path: ", name, call. = FALSE)
    path <- normalizePath(path, mustWork = FALSE)
    if(dirname(path) != data_lake || dir.exists(path))
      stop("Error - refusing to remove an unsafe sites file path: ", name, call. = FALSE)
    if(file.exists(path) && !grepl("[.]parquet$", path))
      stop("Error - existing sites file is not a Parquet file: ", name, call. = FALSE)
    path
  }
  verify_record <- function(connection, upload){
    stored <- DBI::dbGetQuery(connection, select_sql, params = upload$params)
    if(nrow(stored) != 1L || !isTRUE(as.character(stored$data_file_name[[1L]]) == upload$file) ||
       !isTRUE(as.numeric(stored$total_sites[[1L]]) == upload$total))
      stop("Error - sites database verification failed for ", paste(upload$params, collapse = "/"), call. = FALSE)
  }
  
  DBI::dbBegin(conn)
  transaction_open <- TRUE
  for(upload in uploads){
    old <- DBI::dbGetQuery(conn, paste(select_sql, "FOR UPDATE"), params = upload$params)
    if(nrow(old) > 1L)
      stop("Error - multiple sites records found for one primary key.", call. = FALSE)
    if(nrow(old)){
      # Detect input keys that are distinct in R but equal under the DB collation.
      if(!is.na(old$data_file_name[[1L]]) && old$data_file_name[[1L]] %in% new_names)
        stop("Error - supplied site groups collide under the database key collation.", call. = FALSE)
      old_path(as.character(old$data_file_name[[1L]]))
      old_files <- c(old_files, as.character(old$data_file_name[[1L]]))
    }
    deleted <- DBI::dbExecute(conn, paste("DELETE FROM sites WHERE", key_where), params = upload$params)
    if(!isTRUE(deleted == nrow(old))) stop("Error - unexpected sites deletion count.", call. = FALSE)
    inserted <- DBI::dbExecute(conn, insert_sql, params = c(upload$params, list(upload$total, upload$file)))
    if(!isTRUE(inserted == 1L)) stop("Error - failed to insert a sites record.", call. = FALSE)
    verify_record(conn, upload)
  }
  commit_attempted <- TRUE
  tryCatch({
    # Older RMariaDB dbCommit() methods ignore C API errors. The SQL query path
    # checks the server's response; dbCommit() then clears driver bookkeeping.
    DBI::dbExecute(conn, "COMMIT", immediate = TRUE)
    DBI::dbCommit(conn)
  }, error = function(e) {
    stop("Error - sites COMMIT failed; its outcome may be uncertain. Old and new Parquet files were retained. ",
         conditionMessage(e), call. = FALSE)
  })
  transaction_open <- FALSE
  
  # A new connection must see the replacements before old files can be removed.
  tryCatch({
    verify <- createDBconnection()
    tryCatch({
      for(upload in uploads) verify_record(verify, upload)
    }, finally = DBI::dbDisconnect(verify))
  }, error = function(e) {
    stop("Error - sites COMMIT completed but independent verification failed. Old and new Parquet files were retained. ",
         conditionMessage(e), call. = FALSE)
  })
  
  for(name in unique(old_files[!is.na(old_files) & nzchar(old_files)])){
    path <- old_path(name)
    if(!file.exists(path)) next
    aliases <- unique(c(name, basename(path), path, file.path(.sitesDBDataLake, basename(path))))
    placeholders <- paste(rep("?", length(aliases)), collapse = ", ")
    referenced <- FALSE
    for(table in c("sites", "fragments")){
      refs <- tryCatch(DBI::dbGetQuery(conn, paste0("SELECT data_file_name FROM ", table,
                                                    " WHERE data_file_name IN (", placeholders, ") LIMIT 1"), params = as.list(aliases)),
                       error = function(e) stop("Error - sites upload committed, but old-file reference checking failed. Retained ",
                                                path, ": ", conditionMessage(e), call. = FALSE))
      referenced <- referenced || nrow(refs) > 0L
    }
    if(referenced){
      log_message(paste0("Retained Parquet still referenced by another record: ", path))
    } else if(unlink(path) != 0L || file.exists(path)){
      stop("Error - sites upload committed, but could not remove obsolete Parquet file: ", path, call. = FALSE)
    }
  }
  log_message(paste0("Sites upload committed and verified: ", length(uploads), " record(s)."))
  invisible(TRUE)
}



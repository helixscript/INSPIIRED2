tmpString <- function() paste0(Sys.getpid(), '___', stringi::stri_rand_strings(1, 30, '[A-Za-z0-9]'))


ppNum <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)


# Canonical sample-sheet identifiers. Keep IDs as character at TSV import so
# leading zeros are not lost before normalization. Formatting variants are aliases.
normalizeSampleIDs <- function(x, field = "identifier", strict = TRUE){
  if(!(is.character(x) || is.factor(x)) || !is.null(dim(x)))
    stop("Error - ", field, " identifiers must be a character or factor vector; read sample-sheet IDs as text.", call. = FALSE)
  if(!is.logical(strict) || length(strict) != 1L || is.na(strict))
    stop("Error - strict must be TRUE or FALSE.", call. = FALSE)
  original_names <- names(x)
  values <- enc2utf8(as.character(x))
  # PCRE horizontal/vertical whitespace includes tabs and non-breaking spaces.
  values <- gsub("[\\h\\v]+", " ", values, perl = TRUE)
  values <- trimws(values, whitespace = " ")
  values <- gsub(" ", "_", values, fixed = TRUE)
  # ASCII case conversion is independent of the process locale (e.g. Turkish i).
  values <- chartr("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", values)
  names(values) <- original_names
  if(strict){
    invalid <- is.na(values) | !nzchar(values) | values == "NA" |
      !grepl("^[A-Z0-9_.-]+$", values, perl = TRUE)
    if(any(invalid))
      stop("Error - invalid ", field, " identifier(s) at row(s) ",
           paste(which(invalid), collapse = ", "),
           ". After normalization, IDs must be non-empty and contain only A-Z, 0-9, underscore, hyphen or period; NA is reserved for missing values.",
           call. = FALSE)
  }
  values
}


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
    ### message(paste0('Connecting to the INSPIIRED database using cnf file "',  args$dbConfigFile, '" and group ID "', args$dbConfigID, '".'))
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


startModule <- function(connectDB = TRUE){
  options(useFancyQuotes = FALSE, scipen = 999)
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
  
  # Readers open their own connection and do not need a writable data lake.
  if(connectDB && all(c('dbConfigFile', 'dbConfigID') %in% names(args))){
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

# Pipeline modules supply the validated four-key scope to replace all modes in
# each processed sample/genome. Every upload requires complete sites provenance.
# Omitted scope retains the legacy exact-mode helper behavior, with the same gate.
uploadSitesToDB <- function(sites, processedGroups = NULL){
  validateSiteRebuild <- function(x, provenance){
    instruction <- paste0(' Re-run buildSites using the complete output of ',
                          'buildStdFragments --pullSubjectFragments with database credentials, ',
                          'then pass the complete sites RDS through the annotation chain.')
    if(!is.data.frame(x) || !nrow(x) || anyDuplicated(names(x)))
      stop('Error - database site updates require a non-empty sites table with unique column names.',
           instruction, call. = FALSE)
    flag <- x[['pullSubjectFragments']]
    if(!(is.factor(flag) || is.logical(flag)) || !is.null(dim(flag)) ||
       length(flag) != nrow(x) || anyNA(flag) ||
       !all(as.character(flag) == 'TRUE'))
      stop('Error - database site updates require pullSubjectFragments=TRUE for every input site; missing, FALSE, NA, or mixed values are not accepted.',
           instruction, call. = FALSE)
    required <- c('version', 'pullSubjectFragments', 'siteRows', 'groups')
    if(!is.list(provenance) || anyDuplicated(names(provenance)) ||
       !setequal(names(provenance), required) ||
       !identical(provenance$version, 2L) ||
       !identical(provenance$pullSubjectFragments, TRUE))
      stop('Error - database site updates require valid version-2 inspiired2.siteRebuild provenance.',
           instruction, call. = FALSE)
    rows <- provenance$siteRows
    if(!is.numeric(rows) || !is.null(dim(rows)) || length(rows) != 1L ||
       is.na(rows) || !is.finite(rows) || rows != nrow(x))
      stop('Error - the site row count does not match its site rebuild provenance.',
           instruction, call. = FALSE)
    
    keyColumns <- c('trial', 'subject', 'sample', 'refGenome')
    normalizeKeys <- function(tab, label){
      if(!is.data.frame(tab) || anyDuplicated(names(tab)) ||
         !all(keyColumns %in% names(tab)))
        stop('Error - ', label, ' must contain trial, subject, sample, and refGenome.',
             instruction, call. = FALSE)
      keys <- as.data.frame(tab)[, keyColumns, drop = FALSE]
      for(column in keyColumns){
        values <- keys[[column]]
        if(!is.atomic(values) || !is.null(dim(values)) || is.complex(values) || anyNA(values))
          stop('Error - ', label, ' contains invalid ', column, ' identifiers.',
               instruction, call. = FALSE)
        values <- as.character(values)
        if(any(!nzchar(trimws(values))) || any(values != trimws(values)))
          stop('Error - ', label, ' contains blank or whitespace-padded ', column, ' identifiers.',
               instruction, call. = FALSE)
        keys[[column]] <- values
      }
      keys
    }
    if(!is.data.frame(provenance$groups) ||
       !setequal(names(provenance$groups), keyColumns))
      stop('Error - site rebuild groups must contain exactly trial, subject, sample, and refGenome.',
           instruction, call. = FALSE)
    groups <- normalizeKeys(provenance$groups, 'Site rebuild groups')
    if(!nrow(groups) || anyDuplicated(groups))
      stop('Error - site rebuild groups must be non-empty and unique.',
           instruction, call. = FALSE)
    observed <- unique(normalizeKeys(x, 'Input sites'))
    if(nrow(dplyr::anti_join(observed, groups, by = keyColumns)))
      stop('Error - input sites include a sample/genome group absent from the site rebuild provenance.',
           instruction, call. = FALSE)
    groups
  }
  provenance_groups <- validateSiteRebuild(
    sites, attr(sites, 'inspiired2.siteRebuild', exact = TRUE))
  if(!is.list(args) || !all(c("dbConfigFile", "dbConfigID") %in% names(args)) ||
     any(vapply(args[c("dbConfigFile", "dbConfigID")], function(x)
       !is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x) || x == "none", logical(1))) ||
     !file.exists(args$dbConfigFile))
    stop("Error - uploadSitesToDB requires args$dbConfigID and an existing args$dbConfigFile.", call. = FALSE)
  key_columns <- c("trial", "subject", "sample", "refGenome", "mode")
  db_columns <- c("trial", "subject", "sample", "ref_genome", "mode")
  scope_columns <- key_columns[seq_len(4L)]
  full_scope <- !is.null(processedGroups)
  if(anyDuplicated(names(sites)) || !all(c(key_columns, "posid") %in% names(sites)))
    stop("Error - sites must have unique column names and include trial, subject, sample, refGenome, mode and posid.", call. = FALSE)
  processingProvenance <- provenanceFromObject(sites)
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
  if(full_scope){
    # Only an explicit four-key scope authorizes replacement across every mode.
    # Do not infer this scope from surviving sites or from object attributes.
    if(!is.data.frame(processedGroups) || anyDuplicated(names(processedGroups)) ||
       !setequal(names(processedGroups), scope_columns) || !nrow(processedGroups))
      stop("Error - processedGroups must be a non-empty data frame containing exactly trial, subject, sample and refGenome.", call. = FALSE)
    scope_keys <- as.data.frame(processedGroups)[scope_columns]
    for(column in scope_columns){
      values <- scope_keys[[column]]
      if(!is.atomic(values) || !is.null(dim(values)) || is.complex(values) || anyNA(values))
        stop("Error - processedGroups.", column, " must contain non-missing scalar identifiers.", call. = FALSE)
      values <- as.character(values)
      if(any(!nzchar(trimws(values))) || any(values != trimws(values)))
        stop("Error - processedGroups.", column, " contains blank or whitespace-padded identifiers.", call. = FALSE)
      scope_keys[[column]] <- values
    }
    if(anyDuplicated(scope_keys))
      stop("Error - processedGroups must contain unique sample/genome groups.", call. = FALSE)
    if(nrow(dplyr::anti_join(scope_keys, provenance_groups, by = scope_columns)) ||
       nrow(dplyr::anti_join(provenance_groups, scope_keys, by = scope_columns)))
      stop("Error - processedGroups must match the complete sample/genome scope in the sites provenance.", call. = FALSE)
    key_ids <- function(x){
      if(!nrow(x)) return(character())
      do.call(paste0, lapply(x, function(values)
        paste0(nchar(values, type = "bytes"), ":", values)))
    }
    site_scope <- match(key_ids(keys[scope_columns]), key_ids(scope_keys))
    if(anyNA(site_scope))
      stop("Error - sites include a sample/genome group absent from processedGroups.", call. = FALSE)
  }
  groups <- if(nrow(keys))
    dplyr::group_rows(dplyr::group_by(keys, dplyr::across(dplyr::all_of(key_columns)))) else list()
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
  locked <- transaction_open <- verified <- FALSE
  published_files <- character()
  on.exit({
    if(transaction_open) tryCatch(DBI::dbRollback(conn), error = function(e) {
      warning("Upload rollback failed: ", conditionMessage(e), call. = FALSE)
    })
    # A checksum-named file can be reused by another upload. As in buildFragments,
    # retain finalized files on failure; only temporary/pending files are removed.
    if(!verified && length(published_files))
      message("Prepared Parquet files retained after an unsuccessful or unverified sites upload: ",
              paste(published_files, collapse = ", "))
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
  primary <- DBI::dbGetQuery(conn, paste(
    "SELECT COLUMN_NAME AS column_name FROM information_schema.KEY_COLUMN_USAGE",
    "WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'sites'",
    "AND CONSTRAINT_NAME = 'PRIMARY' ORDER BY ORDINAL_POSITION"))
  if(!identical(as.character(primary$column_name), db_columns))
    stop("Error - sites must have the five-column sample/genome/mode primary key.", call. = FALSE)
  lengths <- setNames(as.numeric(schema$max_length), schema$column_name)
  for(i in seq_along(key_columns)){
    limit <- lengths[[db_columns[[i]]]]
    values <- keys[[i]]
    if(full_scope && i <= length(scope_columns)) values <- c(values, scope_keys[[i]])
    if(is.na(limit) || any(nchar(values, type = "chars") > limit)){
      hint <- if(db_columns[[i]] == "mode")
        " Run migrations/sites_mode_width.sql to support 'dual detect'." else ""
      stop("Error - sites.", db_columns[[i]], " is too short for the supplied identifiers.", hint, call. = FALSE)
    }
  }
  if(is.na(lengths[["data_file_name"]]) || lengths[["data_file_name"]] < 40L)
    stop("Error - sites.data_file_name must allow at least 40 characters.", call. = FALSE)
  
  key_where <- paste(paste(db_columns, "= ?"), collapse = " AND ")
  select_sql <- paste("SELECT data_file_name, total_sites FROM sites WHERE", key_where)
  replacement_where <- if(full_scope)
    paste(paste(db_columns[seq_along(scope_columns)], "= ?"), collapse = " AND ") else key_where
  replacement_select <- paste("SELECT data_file_name, total_sites FROM sites WHERE", replacement_where)
  insert_sql <- paste(
    "INSERT INTO sites (trial, subject, sample, ref_genome, mode, total_sites, processed_date, data_file_name)",
    "VALUES (?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP, ?)"
  )
  uploads <- vector("list", length(groups))
  is_symlink <- function(path){
    target <- Sys.readlink(path)
    !is.na(target) && nzchar(target)
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
    local <- tempfile("sites_", tmpdir = tmp_dir, fileext = ".parquet")
    pending <- tempfile(".sites_", tmpdir = data_lake, fileext = ".pending")
    on.exit(unlink(c(local, pending)), add = TRUE)
    provenanceWriteParquet(x, local, provenance = processingProvenance)
    checksum <- if(file.exists(local) && file.size(local) > 0) unname(tools::md5sum(local)) else NA_character_
    if(length(checksum) != 1L || is.na(checksum) || !grepl("^[[:xdigit:]]{32}$", checksum))
      stop("Error - failed to write and checksum sites Parquet file.", call. = FALSE)
    checksum <- tolower(checksum)
    filename <- paste0(checksum, ".parquet")
    path <- file.path(data_lake, filename)
    verify_file <- function(candidate){
      if(is_symlink(candidate) || !utils::file_test("-f", candidate))
        stop("Error - sites Parquet path is not a regular, non-symlink file: ", candidate, call. = FALSE)
      if(!identical(unname(tools::md5sum(candidate)), checksum))
        stop("Error - sites Parquet file does not match its checksum-derived name: ", path, call. = FALSE)
      # Equal MD5 values alone cannot rule out an actual hash collision.
      if(!same_bytes(local, candidate))
        stop("Error - MD5 collision: different Parquet bytes have the same filename: ", path, call. = FALSE)
    }
    if(file.exists(path) || is_symlink(path)){
      verify_file(path)
      return(filename)
    }
    if(!isTRUE(file.copy(local, pending, overwrite = FALSE)) ||
       !identical(unname(tools::md5sum(pending)), checksum))
      stop("Error - failed to copy and verify sites Parquet file in the data lake.", call. = FALSE)
    # Both paths are in the lake. A hard link publishes the complete file
    # atomically and cannot overwrite an existing name (unlike POSIX rename).
    # If another uploader won the race, verify its file before reusing it.
    if(isTRUE(suppressWarnings(file.link(pending, path)))){
      published_files <<- c(published_files, path)
    } else if(!file.exists(path) && !is_symlink(path)){
      stop("Error - could not publish sites Parquet without overwriting; the data-lake filesystem must support hard links: ",
           path, call. = FALSE)
    }
    verify_file(path)
    filename
  }
  
  # Prepare every file before deleting any database records. Existing content
  # is verified and reused; it is never overwritten by this uploader.
  for(i in seq_along(groups)){
    indices <- groups[[i]]
    filename <- stage_file(data[indices, , drop = FALSE])
    uploads[[i]] <- list(params = unname(as.list(keys[indices[[1L]], , drop = FALSE])),
                         total = length(unique(data$posid[indices])), file = filename)
    log_message(paste0("Prepared sites Parquet: ", filename))
  }
  if(full_scope){
    group_scope <- vapply(groups, function(indices) site_scope[indices[[1L]]], integer(1))
    expected_counts <- tabulate(group_scope, nbins = nrow(scope_keys))
    replacements <- lapply(seq_len(nrow(scope_keys)), function(i)
      list(params = unname(as.list(scope_keys[i, , drop = FALSE])), expected = expected_counts[[i]]))
  } else {
    replacements <- lapply(uploads, function(upload)
      list(params = upload$params, expected = 1L))
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
  verify_scope <- function(connection){
    # Checking inserted records alone cannot prove that obsolete modes or an
    # entirely eliminated sample were cleared. Verify absence within every scope.
    for(replacement in replacements){
      stored <- DBI::dbGetQuery(connection, replacement_select, params = replacement$params)
      if(nrow(stored) != replacement$expected)
        stop("Error - sites replacement-scope verification failed for ",
             paste(replacement$params, collapse = "/"), ": expected ", replacement$expected,
             " record(s), found ", nrow(stored), ".", call. = FALSE)
    }
    for(upload in uploads) verify_record(connection, upload)
    invisible(TRUE)
  }
  
  DBI::dbBegin(conn)
  transaction_open <- TRUE
  # Delete all old rows before inserting any replacements. If two input groups
  # compare equal under the DB collation, the primary key rejects the second
  # INSERT and the entire batch rolls back. Filenames cannot detect this: an
  # unchanged re-upload legitimately has the same filename as the old record.
  for(replacement in replacements){
    old <- DBI::dbGetQuery(conn, paste(replacement_select, "FOR UPDATE"), params = replacement$params)
    if(!full_scope && nrow(old) > 1L)
      stop("Error - multiple sites records found for one primary key.", call. = FALSE)
    if(nrow(old)){
      filenames <- as.character(old$data_file_name)
      for(name in filenames) old_path(name)
      old_files <- c(old_files, filenames)
    }
    deleted <- DBI::dbExecute(conn, paste("DELETE FROM sites WHERE", replacement_where), params = replacement$params)
    if(!isTRUE(deleted == nrow(old))) stop("Error - unexpected sites deletion count.", call. = FALSE)
  }
  for(upload in uploads){
    inserted <- DBI::dbExecute(conn, insert_sql, params = c(upload$params, list(upload$total, upload$file)))
    if(!isTRUE(inserted == 1L)) stop("Error - failed to insert a sites record.", call. = FALSE)
  }
  verify_scope(conn)
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
      verify_scope(verify)
    }, finally = DBI::dbDisconnect(verify))
  }, error = function(e) {
    stop("Error - sites COMMIT completed but independent verification failed. Old and new Parquet files were retained. ",
         conditionMessage(e), call. = FALSE)
  })
  verified <- TRUE
  
  for(name in unique(old_files[!is.na(old_files) & nzchar(old_files)])){
    path <- old_path(name)
    if(!file.exists(path)) next
    if(basename(path) %in% new_names) next
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
  log_message(paste0("Sites upload committed and verified: ", length(uploads), " record(s) across ",
                     length(replacements), if(full_scope) " sample/genome replacement group(s)." else " mode replacement group(s)."))
  invisible(TRUE)
}

# Processing provenance is independent of inspiired2.siteRebuild, which remains
# the authorization gate for complete site replacements. Version 1 stores the
# complete ancestry inside each artifact; it is not a database success ledger.
.provenanceAttribute <- "inspiired2.provenance"
.provenanceMetadataKey <- "inspiired2.provenance"

.provenancePlain <- function(x){
  if(is.null(x)) return(NULL)
  if(inherits(x, "POSIXt")) return(format(x, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"))
  if(inherits(x, "Date")) return(as.character(x))
  if(is.factor(x)) x <- as.character(x)
  if(is.data.frame(x)){
    return(list(columns = as.list(names(x)), rows = lapply(seq_len(nrow(x)), function(i)
      lapply(as.list(x[i, , drop = FALSE]), function(v) .provenancePlain(v[[1L]])))))
  }
  if(is.list(x)){
    if(!is.null(names(x)) && (anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))))
      stop("Error - provenance lists must have unique non-empty names or no names.", call. = FALSE)
    return(lapply(x, .provenancePlain))
  }
  if(!is.atomic(x) || !is.null(dim(x)) || is.complex(x) || is.raw(x))
    stop("Error - unsupported object in provenance; use plain scalar values, lists, or data frames.", call. = FALSE)
  if(length(x) != 1L || !is.null(names(x))) return(lapply(as.list(x), .provenancePlain))
  if(is.na(x)) return(NULL)
  if(is.numeric(x)){
    if(!is.finite(x)) return(as.character(x))
    return(as.double(x))
  }
  if(is.logical(x)) return(as.logical(x))
  enc2utf8(as.character(x))
}

.provenanceYAML <- function(x){
  if(!requireNamespace("yaml", quietly = TRUE))
    stop("Error - the yaml package is required for processing provenance.", call. = FALSE)
  yaml::as.yaml(.provenancePlain(x), precision = 17, indent.mapping.sequence = TRUE)
}

.provenanceHash <- function(x){
  path <- tempfile("inspiired2_provenance_")
  on.exit(unlink(path), add = TRUE)
  writeBin(charToRaw(enc2utf8(x)), path)
  value <- unname(tools::md5sum(path))
  if(length(value) != 1L || is.na(value)) stop("Error - could not checksum provenance.", call. = FALSE)
  value
}

.provenanceIDs <- function(x, label){
  if(is.list(x) && !length(x)) return(character())
  if(is.list(x) && all(vapply(x, function(v) is.character(v) && length(v) == 1L && !is.na(v), logical(1))))
    x <- unlist(x, use.names = FALSE)
  if(!is.character(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x))
    stop("Error - invalid ", label, " in processing provenance.", call. = FALSE)
  unname(x)
}

provenanceValidate <- function(bundle){
  if(!is.list(bundle) || is.null(names(bundle)) || anyDuplicated(names(bundle)) ||
     !setequal(names(bundle), c("version", "heads", "records", "artifact")) ||
     length(bundle$version) != 1L || !is.numeric(bundle$version) ||
     !isTRUE(bundle$version == 1))
    stop("Error - unsupported or malformed INSPIIRED2 processing provenance.", call. = FALSE)
  bundle$version <- 1L
  bundle$heads <- .provenanceIDs(bundle$heads, "head IDs")
  records <- bundle$records
  if(!is.list(records) || (length(records) && (is.null(names(records)) || anyNA(names(records)) ||
                                               any(!nzchar(names(records))) || anyDuplicated(names(records)))))
    stop("Error - invalid provenance record map.", call. = FALSE)
  for(id in names(records)){
    r <- records[[id]]
    if(!is.list(r) || is.null(names(r)) || anyDuplicated(names(r)) ||
       !identical(r$id, id) || length(r$kind) != 1L ||
       !is.character(r$kind) || !r$kind %in% c("execution", "source") ||
       !"parents" %in% names(r))
      stop("Error - malformed provenance record: ", id, call. = FALSE)
    r$parents <- .provenanceIDs(r$parents, paste0("parent IDs for ", id))
    if(r$kind == "execution"){
      if(!is.character(r$module) || length(r$module) != 1L || is.na(r$module) || !nzchar(r$module) ||
         !is.character(r$recorded_at) || length(r$recorded_at) != 1L || is.na(r$recorded_at) ||
         !is.character(r$arguments_yaml) || length(r$arguments_yaml) != 1L || is.na(r$arguments_yaml) ||
         !is.list(r$software) || !is.list(r$details) || !is.list(r$resources) || !is.list(r$input_files))
        stop("Error - incomplete provenance execution: ", id, call. = FALSE)
    } else {
      if(!is.character(r$history) || length(r$history) != 1L ||
         !r$history %in% c("embedded", "unavailable") || !is.list(r$source))
        stop("Error - incomplete provenance source: ", id, call. = FALSE)
    }
    if(any(!r$parents %in% names(records)))
      stop("Error - missing parent records in processing provenance: ", id, call. = FALSE)
    records[[id]] <- r
  }
  if(any(!bundle$heads %in% names(records)) || (length(records) && !length(bundle$heads)))
    stop("Error - missing head records in processing provenance.", call. = FALSE)
  # Iterative topological validation avoids recursion limits on imported graphs.
  remaining <- names(records)
  resolved <- character()
  while(length(remaining)){
    ready <- remaining[vapply(records[remaining], function(r) all(r$parents %in% resolved), logical(1))]
    if(!length(ready)) stop("Error - cycle in processing provenance.", call. = FALSE)
    resolved <- c(resolved, ready)
    remaining <- setdiff(remaining, ready)
  }
  if(!is.list(bundle$artifact)) stop("Error - invalid provenance artifact descriptor.", call. = FALSE)
  # Keep the portable representation identical before and after YAML loading.
  records <- lapply(records, function(r){
    p <- r$parents
    r <- .provenancePlain(r)
    r$parents <- p
    r
  })
  bundle$records <- records
  bundle$artifact <- .provenancePlain(bundle$artifact)
  bundle
}

provenanceToYAML <- function(bundle){
  bundle <- provenanceValidate(bundle)
  bundle$heads <- as.list(bundle$heads)
  bundle$records <- lapply(bundle$records, function(r){ r$parents <- as.list(r$parents); r })
  .provenanceYAML(bundle)
}

provenanceFromYAML <- function(text){
  if(!is.character(text) || length(text) != 1L || is.na(text) || !nzchar(text) ||
     nchar(text, type = "bytes") > 64 * 1024^2)
    stop("Error - invalid or excessively large provenance YAML metadata.", call. = FALSE)
  bundle <- tryCatch(yaml::yaml.load(text, eval.expr = FALSE), error = function(e)
    stop("Error - could not decode provenance YAML: ", conditionMessage(e), call. = FALSE))
  provenanceValidate(bundle)
}

.provenanceEmpty <- function(){
  list(version = 1L, heads = character(), records = list(), artifact = list())
}

.provenanceFile <- function(path, checksum = FALSE){
  if(!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
    stop("Error - provenance file paths must be non-empty scalar strings.", call. = FALSE)
  exists <- file.exists(path) && !dir.exists(path)
  info <- if(exists) file.info(path) else NULL
  out <- list(path = path, exists = exists, identity = "path_size_mtime_only",
              size_bytes = if(exists) as.double(info$size) else NULL,
              modified_at = if(exists) format(info$mtime, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC") else NULL)
  if(checksum){
    if(!exists || file.access(path, 4L) != 0L)
      stop("Error - cannot fingerprint provenance resource: ", path, call. = FALSE)
    out$md5 <- unname(tools::md5sum(path))
    if(length(out$md5) != 1L || is.na(out$md5))
      stop("Error - cannot checksum provenance resource: ", path, call. = FALSE)
    out$identity <- "md5"
  }
  .provenancePlain(out)
}

.provenanceScope <- function(x){
  if(is.data.frame(x)){
    columns <- intersect(c("trial", "subject", "sample", "replicate", "refGenome", "mode"), names(x))
    # Scope depends only on identifiers. Do not convert/copy read sequences or
    # coordinates, and collapse repeated keys before expanding factor labels.
    if(!length(columns)){
      keys <- data.frame()
    } else {
      if(inherits(x, "data.table") && requireNamespace("data.table", quietly = TRUE)){
        keys <- unique(x[, columns, with = FALSE])
      } else {
        keys <- unique(x[, columns, drop = FALSE])
      }
      keys <- unique(as.data.frame(lapply(keys, as.character), stringsAsFactors = FALSE))
    }
    return(list(type = "table", row_count = as.double(nrow(x)),
                scope = .provenancePlain(keys)))
  }
  if(is.list(x)){
    members <- lapply(x, function(y) if(is.data.frame(y)) .provenanceScope(y) else list(type = class(y)[1L]))
    return(list(type = "list", members = members))
  }
  list(type = class(x)[1L])
}

provenanceScope <- function(x, bundle){
  bundle <- provenanceValidate(bundle)
  bundle$artifact <- .provenanceScope(x)
  bundle
}

.provenanceSource <- function(bundle, source, history, artifact){
  if(is.null(source)) source <- list(description = "Unattributed input object")
  if(is.character(source) && length(source) == 1L) source <- list(path = source)
  if(!is.list(source)) stop("Error - provenance source must be a list or path.", call. = FALSE)
  if(!is.null(source$path) && is.null(source$file)) source$file <- .provenanceFile(source$path)
  source$artifact <- artifact
  record <- list(kind = "source", parents = bundle$heads, history = history,
                 source = .provenancePlain(source))
  id <- paste0("source_", .provenanceHash(.provenanceYAML(record)))
  record <- c(list(id = id), record)
  if(!is.null(bundle$records[[id]]) && !identical(.provenanceYAML(bundle$records[[id]]), .provenanceYAML(record)))
    stop("Error - conflicting provenance source ID.", call. = FALSE)
  bundle$records[[id]] <- record
  bundle$heads <- id
  bundle$artifact <- artifact
  provenanceValidate(bundle)
}

provenanceFromObject <- function(x, source = NULL){
  bundle <- attr(x, .provenanceAttribute, exact = TRUE)
  if(is.null(bundle))
    return(.provenanceSource(.provenanceEmpty(), source, "unavailable", .provenanceScope(x)))
  bundle <- provenanceValidate(bundle)
  if(!is.null(source)) bundle <- .provenanceSource(bundle, source, "embedded", bundle$artifact)
  bundle
}

provenanceMerge <- function(bundles){
  if(!is.list(bundles)) stop("Error - provenanceMerge expects a list of bundles.", call. = FALSE)
  out <- .provenanceEmpty()
  for(bundle in bundles){
    bundle <- provenanceValidate(bundle)
    for(id in names(bundle$records)){
      record <- bundle$records[[id]]
      if(!is.null(out$records[[id]]) &&
         !identical(.provenanceYAML(out$records[[id]]), .provenanceYAML(record)))
        stop("Error - conflicting immutable provenance record: ", id, call. = FALSE)
      out$records[[id]] <- record
    }
    out$heads <- union(out$heads, bundle$heads)
  }
  provenanceValidate(out)
}

.provenanceArguments <- function(x){
  if(!is.list(x)) stop("Error - provenance arguments must be a list.", call. = FALSE)
  # Config filenames/groups and the raw command are unnecessary for scientific
  # reproduction and may reveal credentials. Never read credential-file contents.
  redact <- function(value){
    if(is.environment(value) || is.function(value) || isS4(value) || typeof(value) == "externalptr") return(NULL)
    if(is.list(value)){
      if(!is.null(names(value))){
        secret <- grepl("password|passwd|secret|token|credential|dbConn|dbConfig|^command$", names(value), ignore.case = TRUE)
        value <- value[!secret]
      }
      value <- lapply(value, redact)
    }
    value
  }
  .provenancePlain(redact(x))
}

provenanceRecord <- function(module, parents = list(), details = list(), resources = character(),
                             inputFiles = character(), argsSnapshot = args){
  if(!is.character(module) || length(module) != 1L || is.na(module) || !nzchar(module))
    stop("Error - provenance module must be a non-empty string.", call. = FALSE)
  bundle <- provenanceMerge(parents)
  recorded <- format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")
  # tempfile uses the process's temporary-name generator, not R's .Random.seed.
  id <- paste0("run_", .provenanceHash(paste(recorded, Sys.getpid(), tempfile("run_"), sep = "|")))
  root <- argsSnapshot$softwareRoot
  code <- character()
  if(is.character(root) && length(root) == 1L && !is.na(root) && nzchar(root)){
    code <- c(file.path(root, c("VERSION", "inspiired2.R", paste0("modules/", module, ".R"))),
              list.files(file.path(root, "lib"), pattern = "[.]R$", full.names = TRUE))
    code <- unique(code[file.exists(code) & !dir.exists(code)])
  }
  version <- argsSnapshot$version
  if(is.null(version) && length(code) && file.exists(file.path(root, "VERSION")))
    version <- paste("INSPIIRED2", readLines(file.path(root, "VERSION"), n = 1L, warn = FALSE))
  if(is.null(version)) version <- "unavailable" else version <- as.character(version)
  software <- list(version = version, R_version = R.version.string,
                   code = lapply(code, .provenanceFile, checksum = TRUE))
  resources <- unique(as.character(resources))
  inputFiles <- unique(as.character(inputFiles))
  if(anyNA(resources) || any(!nzchar(resources)) || anyNA(inputFiles) || any(!nzchar(inputFiles)))
    stop("Error - missing provenance input/resource paths.", call. = FALSE)
  record <- list(id = id, kind = "execution", parents = bundle$heads, module = module,
                 recorded_at = recorded, software = software,
                 arguments_yaml = .provenanceYAML(.provenanceArguments(argsSnapshot)),
                 details = .provenancePlain(details),
                 resources = lapply(resources, .provenanceFile, checksum = TRUE),
                 input_files = lapply(inputFiles, .provenanceFile, checksum = FALSE))
  bundle$records[[id]] <- record
  bundle$heads <- id
  provenanceValidate(bundle)
}

provenanceAttach <- function(x, bundle){
  attr(x, .provenanceAttribute) <- provenanceScope(x, bundle)
  x
}

provenanceStrip <- function(x){
  attr(x, .provenanceAttribute) <- NULL
  if(is.list(x) && !is.data.frame(x)){
    for(i in seq_along(x)) x[i] <- list(provenanceStrip(x[[i]]))
  }
  x
}

provenanceWriteParquet <- function(x, path, provenance = NULL){
  if(!requireNamespace("arrow", quietly = TRUE))
    stop("Error - the arrow package is required to archive processing provenance.", call. = FALSE)
  if(is.null(provenance)) provenance <- provenanceFromObject(x)
  bundle <- provenanceScope(x, provenance)
  metadata <- provenanceToYAML(bundle)
  # Store the explicit portable metadata once. Other Arrow/R attributes, notably
  # the independent site-rebuild authorization attribute, retain their behavior.
  table <- arrow::Table$create(provenanceStrip(x))
  table$metadata[[.provenanceMetadataKey]] <- metadata
  arrow::write_parquet(table, path)
  reader <- arrow::ParquetFileReader$create(path)
  actual <- reader$GetSchema()$metadata[[.provenanceMetadataKey]]
  if(is.null(actual) || !identical(provenanceToYAML(provenanceFromYAML(actual)), metadata))
    stop("Error - written Parquet processing provenance did not round-trip correctly: ", path, call. = FALSE)
  invisible(path)
}

provenanceReadParquet <- function(path, source = NULL){
  if(!requireNamespace("arrow", quietly = TRUE))
    stop("Error - the arrow package is required to read processing provenance.", call. = FALSE)
  table <- arrow::read_parquet(path, as_data_frame = FALSE)
  text <- table$metadata[[.provenanceMetadataKey]]
  data <- as.data.frame(table)
  if(is.null(source)) source <- list()
  if(is.character(source) && length(source) == 1L) source <- list(description = source)
  if(!is.list(source)) stop("Error - provenance archive source must be a list.", call. = FALSE)
  source$file <- .provenanceFile(path, checksum = TRUE)
  if(is.null(text)){
    # Explicitly missing even if an unrelated/older R-specific attribute exists.
    bundle <- .provenanceSource(.provenanceEmpty(), source, "unavailable", .provenanceScope(data))
  } else {
    bundle <- provenanceFromYAML(text)
    if(!identical(.provenanceYAML(bundle$artifact), .provenanceYAML(.provenanceScope(data))))
      stop("Error - Parquet rows or sample keys do not match their provenance artifact scope: ", path, call. = FALSE)
    bundle <- .provenanceSource(bundle, source, "embedded", bundle$artifact)
  }
  provenanceAttach(data, bundle)
}

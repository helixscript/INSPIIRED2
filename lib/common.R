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
    } else {
      args$dbConn <<- NULL
    }
  }
}



resource_overlay <- function(){
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
      updateLog(paste0('Link created: ', source, ' -> ', destination))
    }
    
    nLinked <- nLinked + 1L
  }
  
  if(nLinked)
    updateLog(paste0('Applied ', nLinked, ' user resource overlay(s) from ', resourceRoot))
  
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


run_blastn_parallel <- function(fastaFile, dbPath, params, threads = 60) {
  
  seqs <- Biostrings::readDNAStringSet(fastaFile)
  if(length(seqs) == 0) return(data.table())
  
  num_chunks <- min(length(seqs), threads)
  chunks <- split(seqs, cut(seq_along(seqs), num_chunks, labels = FALSE))
  
  results_list <- parallel::mclapply(seq_along(chunks), function(i) {
    tmp_chunk <- paste0(fastaFile, ".chunk_", i)
    
    Biostrings::writeXStringSet(chunks[[i]], tmp_chunk)
    
    res <- run_blastn(tmp_chunk, dbPath, params, threads = 1)
    
    if(file.exists(tmp_chunk)) unlink(tmp_chunk)
    return(res)
  }, mc.cores = num_chunks)
  
  return(data.table::rbindlist(results_list))
}


run_blastn <- function(fastaFile, dbPath, params, threads = 1){
  system(paste0("blastn ",  params,
                " -query ", fastaFile, 
                " -db ",    dbPath,
                " -out ",   paste0(fastaFile, '.blastn'),
                " -num_threads ", threads,
                " -outfmt '6 qseqid qstart qend sstart send sstrand length pident gaps gapopen bitscore'"))
  
  hits <- data.table()
  if(file.info(paste0(fastaFile, '.blastn'))$size > 0) hits <- fread(paste0(fastaFile, '.blastn'), col.names = c("qName", "qS", "qE", "vS", "vE", "strand", "len", "pident", "gaps", "gapsopen", "bitscore"))
  invisible(file.remove(paste0(fastaFile, '.blastn')))
  hits
}


parseBLAToutput <- function(f){
  if(!file.exists(f) || file.info(f)$size == 0) return(tibble::tibble())
  
  b <- readr::read_delim(f, delim='\t', col_names=FALSE, col_types='iiiiiiiicciiiciiiiccc')
  names(b) <- c('matches','misMatches','repMatches','nCount','qNumInsert','qBaseInsert','tNumInsert','tBaseInsert',
                'strand','qName','qSize','qStart','qEnd','tName','tSize','tStart','tEnd','blockCount',
                'blockSizes','qStarts','tStarts')
  
  x <- read.table(textConnection(system(paste(file.path(args$softwareRoot, 'bin', 'pslScore.pl'), f), intern=TRUE)), sep='\t')
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

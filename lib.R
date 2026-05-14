tmpString <- function() paste0(Sys.getpid(), '___', stringi::stri_rand_strings(1, 30, '[A-Za-z0-9]'))


ppNum <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)


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


startModule <- function(){
  options(useFancyQuotes = FALSE)
  if(! file.access(args$ramDiskPath, mode = 2) == 0) args$ramDiskPath <- args$outputDir
  
  args$command      <<- paste(commandArgs(trailingOnly = FALSE), collapse = " ")
  args$version      <<- paste('AVIA', readLines(file.path(args$softwareRoot, 'VERSION')))
  args$logFile      <<- file.path(args$outputDir, paste0(args$fileTag, '.log'))
  args$tmpDir       <<- file.path(args$outputDir, paste0(args$fileTag, '_tmp'))
  args$logDir       <<- file.path(args$outputDir, paste0(args$fileTag, '_log'))
  args$defaultDelim <<- '___'
  args$ramDisk      <<- file.path(args$ramDiskPath,  "AVIA", gsub("\\.", "", paste0(format(Sys.time(), "%Y%m%d_%H%M%OS6"), "_", Sys.getpid())))
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
}


parseBLAToutput <- function(f){
  if(! file.exists(f) | file.info(f)$size == 0) return(tibble::tibble())
  b <- readr::read_delim(f, delim = '\t', col_names = FALSE, col_types = 'iiiiiiiicciiiciiiiccc')
  
  x <- read.table(textConnection(system(paste(file.path(args$softwareRoot, 'bin', 'pslScore.pl') ,  f), intern = TRUE)), sep = '\t')
  names(x) <- c('tName', 'tStart', 'tEnd', 'hit', 'pslScore', 'percentIdentity')
  
  names(b) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand',
                'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')
  
  b$queryPercentID <- x$percentIdentity
  b$pslScore <- x$pslScore
  b$tStart   <- x$tStart + 1
  b$tEnd     <- x$tEnd
  b$qWidth   <- b$qEnd - b$qStart + 1
  b$tWidth   <- b$tEnd - b$tStart + 1
  
  dplyr::select(b, qName, matches, strand, qSize, qStart, qEnd, tName, tNumInsert, qNumInsert, tBaseInsert, qBaseInsert, tStart, tEnd, queryPercentID, pslScore, qWidth, tWidth)
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


standardize_positions <- function(df, side = "left", window = 10, local_radius = 2, sd_shrink = 4) {
  if (nrow(df) == 0) return(df)
  
  # 1. Force data.table and stabilize join columns
  dt <- as.data.table(copy(df))
  dt[, `:=`(seqnames = as.character(seqnames), 
            strand = as.character(strand))]
  
  sigma <- window / sd_shrink
  
  # 2. Target the active coordinate with strict rounding for precision
  if (side == "left") {
    dt[, coord := round(as.numeric(start))]
  } else {
    dt[, coord := round(as.numeric(end))]
  }
  
  # 3. Aggregate Strength (Sums Rows 16 & 17 to create a 117-read anchor)
  counts <- dt[, .(reads = sum(as.numeric(reads), na.rm = TRUE)), 
               by = .(seqnames, strand, coord)]
  setorder(counts, seqnames, strand, coord)
  
  # 4. Identify Anchors via Strict Genomic Distance Join
  counts[, `:=`(rad_min = coord - local_radius, rad_max = coord + local_radius)]
  neighbors <- counts[counts, 
                      on = .(seqnames, strand, 
                             coord >= rad_min, 
                             coord <= rad_max), 
                      allow.cartesian = TRUE]
  
  anchor_check <- neighbors[, .(is_anchor = all(i.reads >= reads)), 
                            by = .(seqnames, strand, i.coord)]
  
  counts <- merge(counts, anchor_check, 
                  by.x = c("seqnames", "strand", "coord"), 
                  by.y = c("seqnames", "strand", "i.coord"))
  
  anchors <- counts[is_anchor == TRUE]
  anchors[, anchor_pos := coord] # Explicitly preserve anchor position
  
  # 5. Competitive Mapping via Search Window Join
  counts[, `:=`(win_min = coord - window, win_max = coord + window)]
  mapping <- anchors[counts, 
                     on = .(seqnames, strand, 
                            coord >= win_min, 
                            coord <= win_max), 
                     allow.cartesian = TRUE]
  
  setnames(mapping, "i.coord", "orig_pos")
  mapping[, pull := reads * exp(-((anchor_pos - orig_pos)^2) / (2 * sigma^2))]
  
  best_mapping <- mapping[, .(corrected_coord = anchor_pos[which.max(pull)]), 
                          by = .(seqnames, strand, orig_pos)]
  
  # 6. Final Join and Fallback Update with Type Safety
  best_mapping[, orig_pos := as.numeric(orig_pos)]
  
  res <- merge(dt, best_mapping, 
               by.x = c("seqnames", "strand", "coord"), 
               by.y = c("seqnames", "strand", "orig_pos"), 
               all.x = TRUE)
  
  # Use fcoalesce with numeric coercion to avoid type mismatch errors
  if (side == "left") {
    res[, start := data.table::fcoalesce(as.numeric(corrected_coord), as.numeric(start))]
  } else {
    res[, end := data.table::fcoalesce(as.numeric(corrected_coord), as.numeric(end))]
  }
  
  # Silent Cleanup: Only remove columns that actually exist in the final table
  temp_cols <- c("coord", "corrected_coord", "rad_min", "rad_max", "win_min", "win_max")
  cols_to_remove <- intersect(names(res), temp_cols)
  if (length(cols_to_remove) > 0) res[, (cols_to_remove) := NULL]
  
  return(res[])
}

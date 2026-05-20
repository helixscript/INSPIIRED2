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
}


parseBLAToutput <- function(f){
  if(! file.exists(f) | file.info(f)$size == 0) return(tibble::tibble())
  b <- readr::read_delim(f, delim = '\t', col_names = FALSE, col_types = 'iiiiiiiicciiiciiiiccc')
  
  x <- read.table(textConnection(system(paste(file.path(args$softwareRoot, 'bin', 'pslScore.pl') ,  f), intern = TRUE)), sep = '\t')
  names(x) <- c('tName', 'tStart', 'tEnd', 'hit', 'pslScore', 'percentIdentity')
  
  names(b) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand',
                'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')
  
  b$queryPercentID <- as.numeric(x$percentIdentity)
  b$pslScore <- as.numeric(x$pslScore)
  b$tStart   <- as.integer(x$tStart + 1)
  b$tEnd     <- as.integer(x$tEnd)
  b$qWidth   <- as.integer(b$qEnd - b$qStart + 1)
  b$tWidth   <- as.integer(b$tEnd - b$tStart + 1)
  
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




build_multiHit_clusters <- function(frags_multPosIDs){
  # Ensure incoming data is treated as a data.table for maximum optimization
  setDT(frags_multPosIDs)
  
  # 2. Vectorized Network Builder & Clonal Abundance Processor
  # ------------------------------------------------------------------------------
  multiHitClusters <- data.table()
    
    # Grouping columns (trial, subject, sample) are processed natively at the C-level
    multiHitClusters <- frags_multPosIDs[, {
      
      # Filter out reads that only hit a single position (breakpoint-only variation)
      valid_reads_dt <- .SD[, .(pos_count = uniqueN(posid)), by = readID][pos_count > 1]
      
      if (nrow(valid_reads_dt) == 0) {
        # Return empty schema mirroring your target output layout if no true multi-hits exist
        data.table(
          clusterID = character(), nodes = integer(), reads = integer(), UMIs = integer(),
          posids = list(), readIDs = list(), clusterSonicLengths = numeric(), nodeSonicLengths = list()
        )
      } else {
        # Isolate only the valid multi-hit data for processing
        sub_sd <- .SD[readID %in% valid_reads_dt$readID]
        
        # Step A: Batch Execute CD-HIT (Once per sample grouping to kill disk I/O loops)
        # ------------------------------------------------------------------------
        unique_seqs <- unique(sub_sd[, .(readID, adrift_seq)])
        
        # Establish clean path and safe prefix configuration
        ts_id <- paste0("mhc_", sample[1], "_", data.table::frank(unique_seqs)[1])
        tmp_dir <- args$ramDisk
        if(!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
        
        fasta_path <- file.path(tmp_dir, paste0(ts_id, '.fasta'))
        out_prefix <- file.path(tmp_dir, paste0(ts_id, "_cdhit"))
        
        # Fast vectorized FASTA assembly (avoids memory-hogging string aggregations)
        fasta_lines <- character(nrow(unique_seqs) * 2)
        fasta_lines[c(TRUE, FALSE)] <- paste0(">", unique_seqs$readID)
        fasta_lines[c(FALSE, TRUE)] <- unique_seqs$adrift_seq
        writeLines(fasta_lines, con = fasta_path)
        
        # Execute system call using your specified cluster parameters
        cmd <- paste0("cd-hit-est ", args$anchorReadClusterParams, 
                      " -T ", args$threads, 
                      " -i ", fasta_path, " -o ", out_prefix)
        system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
        
        clstr_path <- paste0(out_prefix, ".clstr")
        if(!file.exists(clstr_path)) stop('Error - cd-hit-est failed to return a clstr file.')
        
        # Parse the file directly back using your custom text reader function
        cdhit_lookup <- as.data.table(parse_cdhit_clstr(clstr_path))
        setkey(cdhit_lookup, readID)
        
        # Instantly wipe temporary scratch disk files
        unlink(c(fasta_path, out_prefix, clstr_path, paste0(out_prefix, ".bak")))
        
        # Step B: Build Bipartite Graph Projection
        # ------------------------------------------------------------------------
        # Establish edge mappings exclusively from readID <-> genomic posid
        bipartite_edges <- unique(sub_sd[, .(from = readID, to = posid)])
        
        # Extract flat connected components instantly 
        g <- graph_from_data_frame(bipartite_edges, directed = FALSE)
        comp <- components(g)
        
        # Build mapping table of nodes to cluster configurations
        mem_dt <- data.table(node_name = names(comp$membership), comp_num = comp$membership)
        mem_dt[, clusterID := paste0("MHC.", comp_num)]
        
        # Isolate positions and map back their cluster identifiers
        pos_mem <- mem_dt[node_name %in% sub_sd$posid]
        setnames(pos_mem, "node_name", "posid")
        
        # Merge network clusters, original inputs, and CD-HIT group identifiers
        dt_joined <- merge(sub_sd, pos_mem, by = "posid")
        dt_joined <- merge(dt_joined, cdhit_lookup[, .(readID, cluster_id)], by = "readID", all.x = TRUE)
        
        # Step C: Aggregate Features Simultaneously
        # ------------------------------------------------------------------------
        ans <- dt_joined[, {
          u_posids   <- unique(posid)
          u_reads    <- unique(readID)
          u_umis     <- unique(UMI)
          
          # Calculate cluster abundance via unique CD-HIT groups inside the network
          tot_sonic  <- uniqueN(cluster_id)
          
          # Calculate node-specific abundance breakdowns matching your layout format
          node_table <- .SD[, .(sonicLengths = uniqueN(cluster_id)), by = .(posid = posid)][, .(posid, sonicLengths)]
          
          .(
            nodes               = length(u_posids),
            reads               = length(u_reads),
            UMIs                = length(u_umis),
            posids              = list(u_posids),
            readIDs             = list(u_reads),
            clusterSonicLengths = tot_sonic,
            nodeSonicLengths    = list(node_table)
          )
        }, by = .(clusterID)]
        
        ans
      }
    }, by = .(trial, subject, sample)]
    
    # Step 3: Global Serialization and Metadata Ingestion
    # ------------------------------------------------------------------------------
    saveRDS(multiHitClusters, file.path(args$outputDir, paste0(args$fileTag, '_multiHitClusters.rds')))
    
    multiHitClusters
}


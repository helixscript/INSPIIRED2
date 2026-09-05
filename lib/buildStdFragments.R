# buildStdFragments

#' Standardize Genomic Jitter using Gaussian Weighting
#'
#' This function standardizes high-throughput genomic fragment boundaries by 
#' executing an optimized, two-tier "Competitive Mapping" algorithm engineered 
#' for performance at scale using data.table.
#'
#' @param df A data.table or data.frame containing genomic fragment 
#'   records. Must include seqnames, strand, reads, 
#'   start, and end columns.
#' @param side Character string, either "left" or "right". Controls 
#'   which side of the genomic fragments is targeted as the active track for 
#'   stabilization. Setting it to "left" focuses processing and 
#'   adjustments entirely on the fragment start positions, whereas "right" 
#'   targets the fragment end positions. Changing this parameter 
#'   shifts the directional focus, allowing back-to-back runs to cleanly "box in" 
#'   both edges. Default is "left".
#' @param window Numeric scalar indicating the maximum search boundary fence (in 
#'   nucleotides) for identifying candidate anchors. It creates an inclusive 
#'   window of plus/minus window nucleotides around each observed coordinate. 
#'   Increasing the value widens the search fence to capture widely dispersed noise, 
#'   while decreasing it (e.g., to 3 or 5) restricts corrections to a localized range, 
#'   preventing distinct biological peaks from blending together. Default is 10.
#' @param local_radius Numeric scalar defining the strict genomic distance threshold 
#'   (in nucleotides) used to identify true local maxima (anchors). A 
#'   coordinate must have an aggregated read count greater than or equal to all neighboring 
#'   sites within a strict plus/minus local_radius span to qualify as an anchor. 
#'   Increasing this parameter minimizes false anchors by ignoring small background spikes across 
#'   a wider area, whereas decreasing it to 1 preserves finer resolution 
#'   by letting closely spaced peak shoulders form separate clusters. Default is 2.
#' @param sd_shrink Numeric scalar acting as a mathematical divider to calculate the 
#'   standard deviation (sigma = window / sd_shrink) of the Gaussian probability 
#'   curve. This controls the "tightness" of the gravitational pull decay. 
#'   Increasing this parameter (e.g., to 6 or 8) sharpens the curve into a narrow spike, 
#'   punishing distance aggressively so only fragments very close to an anchor can snap. 
#'   Decreasing it flattens and widens the curve, broadening its reach so prominent anchors can 
#'   easily grab far-flung reads from the outer tails of the jitter distribution. Default is 4.
#'
#' @return A data.table with updated and standardized start or end coordinates, 
#'   depending on the side evaluated. Temporary math and range columns are silently 
#'   cleaned up prior to return.
#' @export
standardize_positions <- function(df, side = "left", window = 8, local_radius = 2, sd_shrink = 4) {
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
  setDT(frags_multPosIDs)
  nt_len <- args$multiHitclusteringNTlen
  save_details <- isTRUE(args$saveMultiHitClusteringDetails)
  
  if(length(nt_len) != 1L || is.na(nt_len) || nt_len < 1L)
    stop("Error - multiHitclusteringNTlen must be a positive integer.")
  if(!"adrift_seq" %in% names(frags_multPosIDs))
    stop("Error - adrift_seq is missing from the multi-hit fragment table.")
  
  if(!dir.exists(args$ramDisk)) dir.create(args$ramDisk, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(args$ramDisk))
    stop("Error - unable to create the multi-hit clustering temporary directory.")
  
  multiHitClusters <- frags_multPosIDs[, {
    valid_reads_dt <- .SD[, .(pos_count = uniqueN(posid)), by = readID][pos_count > 1]
    
    if(nrow(valid_reads_dt) == 0){
      empty <- data.table(clusterID = character(), nodes = integer(), reads = integer(), UMIs = integer(),
                          posids = list(), readIDs = list(), clusterSonicLengths = integer(),
                          nodeSonicLengths = list())
      if(save_details) empty[, cdhitAssignments := vector("list", .N)]
      empty
    } else {
      sub_sd <- .SD[readID %in% valid_reads_dt$readID]
      adrift_seqs <- as.character(sub_sd$adrift_seq)
      
      if(anyNA(adrift_seqs) || any(nchar(adrift_seqs) < nt_len))
        stop("Error - one or more adrift reads are shorter than ", nt_len,
             " nt or contain missing sequences.")
      
      # Identify sample-level connected-component networks first.
      edge_map <- unique(sub_sd[, .(readID, posid,
                                    from = paste0("read:", readID),
                                    to = paste0("pos:", posid))])
      setorder(edge_map, from, to)
      
      graph_membership <- components(
        graph_from_data_frame(edge_map[, .(from, to)], directed = FALSE)
      )$membership
      
      read_mem <- unique(edge_map[, .(readID, node_name = from)])
      read_mem[, clusterID := paste0("MHC.", unname(graph_membership[node_name]))]
      
      dt_joined <- merge(sub_sd, read_mem[, .(readID, clusterID)],
                         by = "readID", all.x = TRUE, sort = FALSE)
      if(anyNA(dt_joined$clusterID))
        stop("Error - one or more multi-hit reads were not assigned to a network.")
      
      # Run CD-HIT separately within every connected-component network.
      rbindlist(lapply(
        split(dt_joined, by = "clusterID", keep.by = TRUE, sorted = TRUE),
        function(net){
          unique_seqs <- unique(net[, .(
            readID,
            testSeq = substr(as.character(adrift_seq), 1L, nt_len)
          )])
          
          if(nrow(unique_seqs[, .N, by = readID][N != 1L]) > 0)
            stop("Error - a readID has more than one adrift-read sequence within a multi-hit network.")
          
          # Deterministic FASTA order is important because CD-HIT -g 0 is greedy.
          setorder(unique_seqs, readID)
          
          ts <- file.path(args$ramDisk, paste0("mhc_", tmpString()))
          fasta_path <- paste0(ts, ".fasta")
          out_prefix <- paste0(ts, "_cdhit")
          clstr_path <- paste0(out_prefix, ".clstr")
          on.exit(unlink(c(fasta_path, out_prefix, clstr_path,
                           paste0(out_prefix, ".bak"))), add = TRUE)
          
          fasta_lines <- character(nrow(unique_seqs) * 2L)
          fasta_lines[c(TRUE, FALSE)] <- paste0(">", unique_seqs$readID)
          fasta_lines[c(FALSE, TRUE)] <- unique_seqs$testSeq
          writeLines(fasta_lines, fasta_path)
          
          cmd <- paste("cd-hit-est", args$multiHitclusteringParams,
                       "-T", args$threads, "-i", shQuote(fasta_path),
                       "-o", shQuote(out_prefix))
          status <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
          
          if(status != 0L || !file.exists(clstr_path))
            stop("Error - cd-hit-est failed for multi-hit network ", net$clusterID[1], ".")
          
          cdhit_lookup <- as.data.table(parse_cdhit_clstr(clstr_path))
          if(anyDuplicated(cdhit_lookup$readID) ||
             !setequal(unique_seqs$readID, cdhit_lookup$readID))
            stop("Error - incomplete or duplicated CD-HIT assignments for multi-hit network ",
                 net$clusterID[1], ".")
          
          if(save_details){
            assignment_table <- merge(unique_seqs, cdhit_lookup,
                                      by = "readID", all.x = TRUE, sort = FALSE)
            setnames(assignment_table,
                     c("testSeq", "cluster_id", "is_rep", "cluster_size"),
                     c("adriftSeqSegment", "cdhitClusterID", "isRep", "clusterSize"))
            setorder(assignment_table, cdhitClusterID, readID)
          }
          
          net <- merge(net, cdhit_lookup[, .(readID, cluster_id)],
                       by = "readID", all.x = TRUE, sort = FALSE)
          if(anyNA(net$cluster_id))
            stop("Error - missing CD-HIT assignments in multi-hit network ",
                 net$clusterID[1], ".")
          
          ans <- net[, {
            u_posids <- unique(posid)
            u_reads <- unique(readID)
            u_umis <- unique(UMI)
            node_table <- .SD[, .(sonicLengths = uniqueN(cluster_id)), by = posid]
            
            .(nodes = length(u_posids), reads = length(u_reads), UMIs = length(u_umis),
              posids = list(u_posids), readIDs = list(u_reads),
              clusterSonicLengths = uniqueN(cluster_id),
              nodeSonicLengths = list(node_table))
          }, by = clusterID]
          
          if(save_details) ans[, cdhitAssignments := list(assignment_table)]
          ans
        }
      ), use.names = TRUE, fill = FALSE)
    }
  }, by = .(trial, subject, sample, refGenome)]
  
  multiHitClusters
}
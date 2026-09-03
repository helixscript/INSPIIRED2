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
  # Ensure incoming data is treated as a data.table for maximum optimization
  setDT(frags_multPosIDs)
  
  if(length(args$multiHitclusteringNTlen) != 1L ||
     is.na(args$multiHitclusteringNTlen) ||
     args$multiHitclusteringNTlen < 1L){
    stop("Error - multiHitclusteringNTlen must be a positive integer.")
  }
  
  # Vectorized Network Builder & Clonal Abundance Processor
  # ------------------------------------------------------------------------------
  multiHitClusters <- data.table()
  
  # Including refGenome prevents reads aligned to different references from
  # entering the same network.
  multiHitClusters <- frags_multPosIDs[, {
    
    # Filter out reads that only hit a single position
    # (breakpoint-only variation).
    valid_reads_dt <- .SD[
      ,
      .(pos_count = uniqueN(posid)),
      by = readID
    ][pos_count > 1]
    
    if(nrow(valid_reads_dt) == 0){
      # Return an empty schema if no true multi-hits exist.
      data.table(
        clusterID = character(),
        nodes = integer(),
        reads = integer(),
        UMIs = integer(),
        posids = list(),
        readIDs = list(),
        clusterSonicLengths = numeric(),
        nodeSonicLengths = list()
      )
    } else {
      # Isolate valid multi-hit reads.
      sub_sd <- .SD[readID %in% valid_reads_dt$readID]
      
      # Cluster linker-adjacent adrift-read sequences.
      # ------------------------------------------------------------------------
      if(anyNA(sub_sd$adrift_seq) ||
         any(nchar(as.character(sub_sd$adrift_seq)) <
             args$multiHitclusteringNTlen)){
        stop(
          "Error - one or more adrift reads are shorter than ",
          args$multiHitclusteringNTlen,
          " nt or contain missing sequences."
        )
      }
      
      # adrift_seq begins immediately after the linker removed during
      # demultiplexing. Its first nucleotides therefore represent sequence
      # adjacent to the sonic-shearing boundary.
      unique_seqs <- unique(sub_sd[, .(
        readID,
        adrift_end_seq = substr(
          as.character(adrift_seq),
          1L,
          args$multiHitclusteringNTlen
        )
      )])
      
      # Establish temporary file paths.
      ts_id <- paste0(
        "mhc_",
        sample[1], "_",
        refGenome[1], "_",
        data.table::frank(unique_seqs)[1]
      )
      
      tmp_dir <- args$ramDisk
      
      if(!dir.exists(tmp_dir)){
        dir.create(
          tmp_dir,
          recursive = TRUE,
          showWarnings = FALSE
        )
      }
      
      fasta_path <- file.path(
        tmp_dir,
        paste0(ts_id, ".fasta")
      )
      
      out_prefix <- file.path(
        tmp_dir,
        paste0(ts_id, "_cdhit")
      )
      
      # Write the extracted adrift-read ends to FASTA.
      fasta_lines <- character(nrow(unique_seqs) * 2)
      fasta_lines[c(TRUE, FALSE)] <- paste0(
        ">",
        unique_seqs$readID
      )
      fasta_lines[c(FALSE, TRUE)] <- unique_seqs$adrift_end_seq
      
      writeLines(
        fasta_lines,
        con = fasta_path
      )
      
      # Execute CD-HIT using the dedicated multi-hit clustering parameters.
      cmd <- paste0(
        "cd-hit-est ",
        args$multiHitclusteringParams,
        " -T ", args$threads,
        " -i ", fasta_path,
        " -o ", out_prefix
      )
      
      system(
        cmd,
        ignore.stdout = TRUE,
        ignore.stderr = TRUE
      )
      
      clstr_path <- paste0(
        out_prefix,
        ".clstr"
      )
      
      if(!file.exists(clstr_path)){
        stop(
          "Error - cd-hit-est failed to return a clstr file."
        )
      }
      
      # Parse CD-HIT output.
      cdhit_lookup <- as.data.table(
        parse_cdhit_clstr(clstr_path)
      )
      
      setkey(
        cdhit_lookup,
        readID
      )
      
      # Remove temporary files.
      unlink(c(
        fasta_path,
        out_prefix,
        clstr_path,
        paste0(out_prefix, ".bak")
      ))
      
      # Build Bipartite Graph Projection
      # ------------------------------------------------------------------------
      bipartite_edges <- unique(sub_sd[, .(
        from = readID,
        to = posid
      )])
      
      g <- graph_from_data_frame(
        bipartite_edges,
        directed = FALSE
      )
      
      comp <- components(g)
      
      # Map graph vertices to network identifiers.
      mem_dt <- data.table(
        node_name = names(comp$membership),
        comp_num = comp$membership
      )
      
      mem_dt[, clusterID := paste0(
        "MHC.",
        comp_num
      )]
      
      # Isolate position vertices.
      pos_mem <- mem_dt[
        node_name %in% sub_sd$posid
      ]
      
      setnames(
        pos_mem,
        "node_name",
        "posid"
      )
      
      # Add network and sequence-cluster identifiers to each read-position
      # association.
      dt_joined <- merge(
        sub_sd,
        pos_mem,
        by = "posid"
      )
      
      dt_joined <- merge(
        dt_joined,
        cdhit_lookup[, .(
          readID,
          cluster_id
        )],
        by = "readID",
        all.x = TRUE
      )
      
      # Aggregate Network Features
      # ------------------------------------------------------------------------
      ans <- dt_joined[, {
        u_posids <- unique(posid)
        u_reads  <- unique(readID)
        u_umis   <- unique(UMI)
        
        # Estimate network abundance using unique clusters of linker-adjacent
        # adrift-read sequences.
        tot_sonic <- uniqueN(cluster_id)
        
        # Estimate abundance separately for each candidate position.
        node_table <- .SD[
          ,
          .(sonicLengths = uniqueN(cluster_id)),
          by = .(posid)
        ][, .(
          posid,
          sonicLengths
        )]
        
        .(
          nodes = length(u_posids),
          reads = length(u_reads),
          UMIs = length(u_umis),
          posids = list(u_posids),
          readIDs = list(u_reads),
          clusterSonicLengths = tot_sonic,
          nodeSonicLengths = list(node_table)
        )
      }, by = .(clusterID)]
      
      ans
    }
  }, by = .(
    trial,
    subject,
    sample,
    refGenome
  )]
  
  saveRDS(
    multiHitClusters,
    file.path(
      args$outputDir,
      paste0(
        args$fileTag,
        "_multiHitClusters.rds"
      )
    )
  )
  
  multiHitClusters
}
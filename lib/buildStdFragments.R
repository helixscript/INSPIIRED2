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




# Splitting is an optional additional report. The original connected components,
# CD-HIT assignments, summary schema and database payload remain unchanged.
multiHitSplitSettings <- function(arguments){
  positive_integer <- function(name, default){
    value <- arguments[[name]]
    if(is.null(value)) value <- default
    if(!is.numeric(value) || length(value) != 1L || is.na(value) ||
       !is.finite(value) || value < 1 || value != floor(value) ||
       value > .Machine$integer.max)
      stop("Error - --", name, " must be a positive integer.", call. = FALSE)
    as.integer(value)
  }
  list(minConnectionReads = positive_integer("multiHitMinConnectionReads", 1L),
       maxCandidates = positive_integer("multiHitMaxSplitCandidates", 100000L),
       minPosids = 2L, minReads = 2L)
}

empty_multiHit_split_tables <- function(){
  list(
    subclusters = data.table(
      clusterID = character(), status = character(), totalReads = integer(),
      readIDs = list(), UMIs = integer(), posids = list(), totalPosids = integer(),
      totalAdriftReadClusters = integer(), totalFragments = integer(),
      adriftReadClusters = list(), fragments = list(), reads = list()),
    readAssignments = data.table(
      readID = character(), subclusterID = character(), role = character(),
      cutID = character(), cdhitClusterID = character(), posids = list(),
      candidateSubclusterIDs = list()),
    positionAssignments = data.table(
      posid = character(), subclusterID = character(), role = character()),
    cuts = data.table(
      cutID = character(), parentCutID = character(), depth = integer(),
      connectingReads = integer(), readIDs = list(), supportedComponents = integer(),
      residualComponents = integer(), positionGroups = list()),
    statistics = data.table(
      splitCount = integer(), bridgeReads = integer(), candidateCutsTested = integer())
  )
}

# Exact bounded search over read-only separators. Identical read neighbourhoods
# can be bundled: removing only some of those reads cannot disconnect positions.
# Position twins can also be bundled without changing any supported component.
# Neither optimisation projects reads into pairwise position edges.
find_multiHit_read_splits <- function(read_pos, settings){
  pairs <- unique(copy(as.data.table(read_pos))[, .(
    readID = as.character(readID), posid = as.character(posid))])
  if(!nrow(pairs) || anyNA(pairs) || any(!nzchar(pairs$readID)) || any(!nzchar(pairs$posid)))
    stop("Error - multi-hit splitting requires nonempty read IDs and position IDs.", call. = FALSE)
  setorder(pairs, readID, posid)
  posids <- sort(unique(pairs$posid), method = "radix")
  pairs[, p := match(posid, posids)]
  read_sets <- pairs[, .(signature = paste(p, collapse = ",")), by = readID]
  read_sets[, bundle := match(signature, unique(signature))]
  read_bundles <- read_sets[, .(readIDs = list(readID), nReads = .N), by = bundle]
  read_bundles[, name := paste0("r", bundle)]
  pairs[read_sets, on = "readID", bundle := i.bundle]
  incidence <- unique(pairs[, .(p, bundle)])
  setorder(incidence, p, bundle)
  pos_sets <- incidence[, .(signature = paste(bundle, collapse = ",")), by = p]
  pos_sets[, posBundle := match(signature, unique(signature))]
  pos_bundles <- pos_sets[, .(posids = list(posids[p]), nPosids = .N), by = posBundle]
  pos_bundles[, name := paste0("p", posBundle)]
  incidence[pos_sets, on = "p", posBundle := i.posBundle]
  edges <- unique(incidence[, .(from = paste0("r", bundle), to = paste0("p", posBundle))])
  vertices <- rbindlist(list(
    read_bundles[, .(name, nReads, nPosids = 0L)],
    pos_bundles[, .(name, nReads = 0L, nPosids)]))
  graph <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = vertices)
  if(igraph::components(graph)$no != 1L)
    stop("Error - the split helper expects one connected parent network.", call. = FALSE)
  
  tests <- 0L
  cuts <- list()
  leaves <- list()
  read_names <- read_bundles$name
  read_weights <- setNames(read_bundles$nReads, read_names)
  expand_reads <- function(nodes) sort(as.character(unlist(read_bundles$readIDs[
    match(intersect(nodes, read_names), read_names)], use.names = FALSE)), method = "radix")
  expand_positions <- function(nodes) sort(as.character(unlist(pos_bundles$posids[
    match(intersect(nodes, pos_bundles$name), pos_bundles$name)], use.names = FALSE)), method = "radix")
  
  evaluate_cut <- function(g, removed){
    if(tests >= settings$maxCandidates)
      stop("Error - multi-hit split search exceeded --multiHitMaxSplitCandidates (",
           settings$maxCandidates, "). Increase the limit or lower ",
           "--multiHitMinConnectionReads; no incomplete split report will be saved.", call. = FALSE)
    tests <<- tests + 1L
    reduced <- igraph::delete_vertices(g, removed)
    membership <- igraph::components(reduced)$membership
    if(length(unique(membership)) < 2L) return(NULL)
    nodes <- split(names(membership), membership)
    nr <- setNames(igraph::vertex_attr(reduced, "nReads"), igraph::V(reduced)$name)
    np <- setNames(igraph::vertex_attr(reduced, "nPosids"), igraph::V(reduced)$name)
    supported <- vapply(nodes, function(x)
      sum(nr[x]) >= settings$minReads && sum(np[x]) >= settings$minPosids, logical(1))
    if(sum(supported) < 2L) return(NULL)
    list(removed = removed, nodes = nodes, supported = supported)
  }
  
  find_cut <- function(g){
    if(settings$minConnectionReads <= 1L ||
       sum(igraph::V(g)$nReads) < 2L * settings$minReads + 1L ||
       sum(igraph::V(g)$nPosids) < 2L * settings$minPosids) return(NULL)
    available <- intersect(read_names, igraph::V(g)$name)
    # Positions sharing >= threshold distinct reads cannot be separated by an
    # allowed cut. Greedily join such position blocks to eliminate impossible
    # candidates, without making a weighted position-edge projection. This is
    # a conservative acceleration: unmerged blocks are still searched exactly.
    # Recompute for each child, since ancestor cuts have removed some reads.
    current_positions <- igraph::V(g)$name[igraph::V(g)$nPosids > 0L]
    current_edges <- edges[from %chin% available & to %chin% current_positions]
    current_edges[, p_index := match(to, current_positions)]
    parents <- seq_along(current_positions)
    covered <- rep(FALSE, length(parents))
    root_of <- function(i){
      while(parents[i] != i) i <- parents[i]
      i
    }
    for(i in seq_along(parents)){
      if(covered[i]) next
      neighbours <- current_edges[p_index == i, from]
      shared <- current_edges[from %chin% neighbours,
                              .(support = sum(as.numeric(read_weights[from]))), by = p_index]
      joined <- shared[support >= settings$minConnectionReads, p_index]
      for(j in joined){
        a <- root_of(i)
        b <- root_of(j)
        if(a != b) parents[max(a, b)] <- min(a, b)
      }
      covered[c(i, joined)] <- TRUE
    }
    blocks <- vapply(seq_along(parents), root_of, integer(1))
    if(uniqueN(blocks) == 1L) return(NULL)
    crossing_bundles <- current_edges[, .(blocks = uniqueN(blocks[p_index])), by = from
    ][blocks > 1L, from]
    available <- available[available %chin% crossing_bundles]
    max_cost <- min(settings$minConnectionReads - 1,
                    sum(igraph::V(g)$nReads) - 2L * settings$minReads)
    candidates <- available[read_weights[available] <= max_cost]
    if(!length(candidates)) return(NULL)
    
    # Every one-read separator is an articulation vertex. Restrict to reads;
    # removing a position would have a different scientific meaning.
    articulation <- igraph::as_ids(igraph::articulation_points(g))
    singles <- candidates[read_weights[candidates] == 1L & candidates %in% articulation]
    for(candidate in singles){
      answer <- evaluate_cut(g, candidate)
      if(!is.null(answer)) return(answer)
    }
    if(max_cost < 2L) return(NULL)
    
    # Enumerate by total number of distinct read IDs, then stable bundle order.
    # This is exhaustive up to the threshold unless the explicit budget errors.
    # Generate combinations lazily instead of materialising choose(N, k) sets.
    weights <- unname(read_weights[candidates])
    suffix_weight <- rev(cumsum(rev(as.numeric(weights))))
    answer <- NULL
    enumerate <- function(start, remaining, chosen){
      if(remaining == 0L){
        answer <<- evaluate_cut(g, chosen)
        return(!is.null(answer))
      }
      if(start > length(candidates) || suffix_weight[start] < remaining) return(FALSE)
      for(i in seq.int(start, length(candidates))){
        if(weights[i] <= remaining &&
           enumerate(i + 1L, remaining - weights[i], c(chosen, candidates[i]))) return(TRUE)
      }
      FALSE
    }
    for(cost in seq.int(2L, max_cost)){
      if(enumerate(1L, cost, character())) return(answer)
    }
    NULL
  }
  
  # An explicit work stack avoids recursion through a long chain of networks.
  work <- list(list(nodes = igraph::V(graph)$name, parentCut = NA_integer_, depth = 0L))
  while(length(work)){
    item <- work[[length(work)]]
    work[[length(work)]] <- NULL
    g <- igraph::induced_subgraph(graph, item$nodes)
    cut <- find_cut(g)
    if(is.null(cut)){
      leaves[[length(leaves) + 1L]] <- list(
        readIDs = expand_reads(item$nodes), posids = expand_positions(item$nodes))
    } else {
      cut_number <- length(cuts) + 1L
      ordering <- order(vapply(cut$nodes, function(x) expand_positions(x)[1L], character(1)), method = "radix")
      groups <- cut$nodes[ordering]
      supported <- cut$supported[ordering]
      cuts[[cut_number]] <- list(
        parentCut = item$parentCut, depth = item$depth,
        readIDs = expand_reads(cut$removed),
        supportedComponents = as.integer(sum(supported)),
        residualComponents = as.integer(sum(!supported)),
        positionGroups = lapply(groups, expand_positions))
      for(i in rev(seq_along(groups))){
        if(supported[i]){
          work[[length(work) + 1L]] <- list(nodes = groups[[i]], parentCut = cut_number,
                                            depth = item$depth + 1L)
        } else {
          leaves[[length(leaves) + 1L]] <- list(
            readIDs = expand_reads(groups[[i]]), posids = expand_positions(groups[[i]]))
        }
      }
    }
  }
  ordering <- order(vapply(leaves, function(x) x$posids[1L], character(1)), method = "radix")
  list(leaves = leaves[ordering], cuts = cuts, candidateCutsTested = tests)
}

build_multiHit_split_report <- function(net, parent_summary, settings){
  result <- empty_multiHit_split_tables()
  parent_id <- parent_summary$clusterID[1L]
  partition <- find_multiHit_read_splits(net[, .(readID, posid)], settings)
  did_split <- length(partition$cuts) > 0L
  members <- list()
  position_maps <- list()
  summaries <- list()
  child_number <- 0L
  for(leaf in partition$leaves){
    if(!length(leaf$readIDs)){
      position_maps[[length(position_maps) + 1L]] <- data.table(
        posid = leaf$posids, subclusterID = NA_character_, role = "bridge_only")
      next
    }
    child_number <- child_number + 1L
    child_id <- if(did_split) paste0(parent_id, ".S", child_number) else parent_id
    status <- if(!did_split) "unsplit" else if(
      length(leaf$readIDs) >= settings$minReads && length(leaf$posids) >= settings$minPosids
    ) "supported" else "residual"
    x <- net[readID %chin% leaf$readIDs]
    positions <- sort(unique(x$posid), method = "radix")
    node_counts <- x[, .(
      adriftReadClusters = uniqueN(cluster_id),
      fragments = uniqueN(as.numeric(fragEnd) - as.numeric(fragStart) + 1),
      reads = uniqueN(readID)), by = posid]
    node_order <- match(positions, node_counts$posid)
    summaries[[length(summaries) + 1L]] <- data.table(
      clusterID = child_id, status = status, totalReads = as.integer(length(leaf$readIDs)),
      readIDs = list(leaf$readIDs), UMIs = as.integer(uniqueN(x$real_UMI)),
      posids = list(positions), totalPosids = as.integer(length(positions)),
      totalAdriftReadClusters = as.integer(uniqueN(x$cluster_id)),
      totalFragments = as.integer(uniqueN(as.numeric(x$fragEnd) - as.numeric(x$fragStart) + 1)),
      adriftReadClusters = list(node_counts$adriftReadClusters[node_order]),
      fragments = list(node_counts$fragments[node_order]), reads = list(node_counts$reads[node_order]))
    members[[length(members) + 1L]] <- data.table(readID = leaf$readIDs, subclusterID = child_id)
    position_maps[[length(position_maps) + 1L]] <- data.table(
      posid = leaf$posids, subclusterID = child_id, role = status)
  }
  result$subclusters <- rbindlist(summaries)
  result$positionAssignments <- rbindlist(position_maps)
  setorder(result$positionAssignments, posid)
  membership <- rbindlist(members)
  assignments <- unique(net[, .(readID, posid, cdhitClusterID = as.character(cluster_id))])
  setorder(assignments, readID, posid)
  result$readAssignments <- assignments[, .(
    subclusterID = membership$subclusterID[match(readID[1L], membership$readID)],
    role = if(readID[1L] %chin% membership$readID) "member" else "bridge",
    cutID = NA_character_, cdhitClusterID = cdhitClusterID[1L], posids = list(posid),
    candidateSubclusterIDs = list(sort(unique(na.omit(result$positionAssignments$subclusterID[
      match(posid, result$positionAssignments$posid)])), method = "radix"))
  ), by = readID]
  if(did_split){
    result$cuts <- rbindlist(lapply(seq_along(partition$cuts), function(i){
      cut <- partition$cuts[[i]]
      cut_id <- paste0(parent_id, ".C", i)
      result$readAssignments[readID %chin% cut$readIDs, cutID := cut_id]
      data.table(cutID = cut_id,
                 parentCutID = if(is.na(cut$parentCut)) NA_character_ else paste0(parent_id, ".C", cut$parentCut),
                 depth = cut$depth, connectingReads = as.integer(length(cut$readIDs)),
                 readIDs = list(cut$readIDs), supportedComponents = cut$supportedComponents,
                 residualComponents = cut$residualComponents, positionGroups = list(cut$positionGroups))
    }))
  }
  # Enforce read and position conservation before the optional report can escape.
  counted_reads <- unlist(result$subclusters$readIDs, use.names = FALSE)
  bridge_reads <- result$readAssignments[role == "bridge", readID]
  if(anyDuplicated(counted_reads) || length(intersect(counted_reads, bridge_reads)) ||
     !setequal(c(counted_reads, bridge_reads), unique(net$readID)) ||
     anyDuplicated(result$positionAssignments$posid) ||
     !setequal(result$positionAssignments$posid, unique(net$posid)) ||
     anyNA(result$readAssignments[role == "bridge", cutID]))
    stop("Error - multi-hit split report failed evidence-conservation checks.", call. = FALSE)
  result$statistics <- data.table(splitCount = as.integer(length(partition$cuts)),
                                  bridgeReads = as.integer(length(bridge_reads)), candidateCutsTested = partition$candidateCutsTested)
  result
}

collect_multiHit_split_reports <- function(clusters, settings){
  templates <- empty_multiHit_split_tables()
  meta_columns <- c("trial", "subject", "sample", "mode", "refGenome")
  metadata <- copy(clusters[, ..meta_columns])
  metadata[, parentClusterID := clusters$clusterID]
  out <- lapply(names(templates), function(name){
    if(!nrow(clusters)) return(cbind(metadata, templates[[name]]))
    rbindlist(lapply(seq_len(nrow(clusters)), function(i){
      tab <- clusters$splitReport[[i]][[name]]
      cbind(metadata[rep(i, nrow(tab))], tab)
    }), use.names = TRUE)
  })
  names(out) <- names(templates)
  c(list(settings = c(settings, list(
    readSupport = "distinct retained readID; nReads is not a cut weight",
    adriftClusters = "inherited parent CD-HIT assignments; counted on subgroup member reads only"))), out)
}

build_multiHit_clusters <- function(frags_multPosIDs){
  setDT(frags_multPosIDs)
  nt_len <- args$multiHitclusteringNTlen
  save_details <- isTRUE(args$saveMultiHitClusteringDetails)
  split_settings <- multiHitSplitSettings(args)
  
  if(length(nt_len) != 1L || is.na(nt_len) || nt_len < 1L)
    stop("Error - multiHitclusteringNTlen must be a positive integer.")
  if(!"adrift_seq" %in% names(frags_multPosIDs))
    stop("Error - adrift_seq is missing from the multi-hit fragment table.")
  if(!all(c("fragStart", "fragEnd") %in% names(frags_multPosIDs)))
    stop("Error - fragStart and fragEnd are required for multi-hit fragment lengths.")
  valid_coordinates <- vapply(frags_multPosIDs[, .(fragStart, fragEnd)], function(x){
    is.numeric(x) && all(is.finite(x) & x >= 1 & x == floor(x))
  }, logical(1))
  if(!all(valid_coordinates) || any(frags_multPosIDs$fragEnd < frags_multPosIDs$fragStart))
    stop("Error - multi-hit fragment coordinates must be finite positive integers with fragEnd >= fragStart.")
  
  empty <- data.table(clusterID = character(), totalReads = integer(), readIDs = list(),
                      UMIs = integer(), posids = list(), totalPosids = integer(),
                      totalAdriftReadClusters = integer(), totalFragments = integer(),
                      adriftReadClusters = list(), fragments = list(), reads = list())
  if(save_details) empty[, cdhitAssignments := vector("list", .N)]
  if(split_settings$minConnectionReads > 1L) empty[, splitReport := vector("list", .N)]
  if(nrow(frags_multPosIDs) == 0L)
    return(cbind(frags_multPosIDs[, .(trial, subject, sample, mode, refGenome)], empty))
  
  if(!dir.exists(args$ramDisk)) dir.create(args$ramDisk, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(args$ramDisk))
    stop("Error - unable to create the multi-hit clustering temporary directory.")
  
  multiHitClusters <- frags_multPosIDs[, {
    valid_reads_dt <- .SD[, .(pos_count = uniqueN(posid)), by = readID][pos_count > 1]
    
    if(nrow(valid_reads_dt) == 0){
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
            u_umis <- unique(real_UMI)
            
            # Coordinates have already been standardized by buildStdFragments.
            # Count distinct lengths independently of the CD-HIT assignments.
            node_counts <- .SD[, .(
              adriftReadClusters = uniqueN(cluster_id),
              fragments = uniqueN(as.numeric(fragEnd) - as.numeric(fragStart) + 1),
              reads = uniqueN(readID)
            ), by = posid]
            node_order <- match(u_posids, node_counts$posid)
            
            .(totalReads = length(u_reads), readIDs = list(u_reads), UMIs = length(u_umis),
              posids = list(u_posids), totalPosids = length(u_posids),
              totalAdriftReadClusters = uniqueN(cluster_id),
              totalFragments = uniqueN(as.numeric(fragEnd) - as.numeric(fragStart) + 1),
              adriftReadClusters = list(node_counts$adriftReadClusters[node_order]),
              fragments = list(node_counts$fragments[node_order]),
              reads = list(node_counts$reads[node_order]))
          }, by = clusterID]
          
          if(save_details) ans[, cdhitAssignments := list(assignment_table)]
          if(split_settings$minConnectionReads > 1L)
            ans[, splitReport := list(tryCatch(
              build_multiHit_split_report(net, ans, split_settings),
              error = function(e) stop("Multi-hit splitting failed for ", net$clusterID[1L],
                                       " (", paste(net[1L, .(trial, subject, sample, mode, refGenome)], collapse = "/"),
                                       "): ", conditionMessage(e), call. = FALSE)))]
          ans
        }
      ), use.names = TRUE, fill = FALSE)
    }
  }, by = .(trial, subject, sample, mode, refGenome)]
  
  multiHitClusters
}





# Return additional, unstandardized fragments for the incoming trial/subject
# pairs. Incoming full database keys take precedence over stored versions.
pullDBfragments <- function(frags){
  if(!is.data.frame(frags)) stop('Error - pullDBfragments requires a fragment data frame.', call. = FALSE)
  input <- as.data.frame(frags)
  empty <- input[FALSE, , drop = FALSE]
  if(!nrow(input)){
    updateLog('pullDBfragments: no incoming trial/subject groups; pulled 0 fragment rows.')
    return(empty)
  }
  key_columns <- c('trial', 'subject', 'sample', 'replicate', 'refGenome', 'mode')
  db_columns <- c('trial', 'subject', 'sample', 'replicate', 'ref_genome', 'mode')
  get_keys <- function(x, label){
    if(anyDuplicated(names(x)) || !all(key_columns %in% names(x)))
      stop('Error - ', label, ' must have unique column names and include ', paste(key_columns, collapse = ', '), '.', call. = FALSE)
    if(any(!vapply(x[key_columns], function(v) is.atomic(v) && is.null(dim(v)), logical(1))))
      stop('Error - ', label, ' contains non-scalar fragment identifiers.', call. = FALSE)
    keys <- as.data.frame(lapply(x[key_columns], as.character), stringsAsFactors = FALSE)
    if(anyNA(keys) || any(vapply(keys, function(v) any(!nzchar(trimws(v))), logical(1))))
      stop('Error - ', label, ' contains missing or blank fragment identifiers.', call. = FALSE)
    reps <- suppressWarnings(as.numeric(keys$replicate))
    if(any(!is.finite(reps) | reps != trunc(reps) | reps > .Machine$integer.max | reps < -.Machine$integer.max))
      stop('Error - ', label, ' contains invalid replicate identifiers.', call. = FALSE)
    keys$replicate <- as.integer(reps)
    keys
  }
  input_keys <- unique(get_keys(input, 'Incoming fragments'))
  pairs <- unique(input_keys[c('trial', 'subject')])
  if(!is.list(args) || !all(c('dbConfigFile', 'dbConfigID') %in% names(args)) ||
     any(vapply(args[c('dbConfigFile', 'dbConfigID')], function(x)
       !is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x) || x == 'none', logical(1))) ||
     !file.exists(args$dbConfigFile))
    stop('Error - pullDBfragments requires --dbConfigID and an existing --dbConfigFile.', call. = FALSE)
  
  conn <- createDBconnection()
  on.exit(tryCatch(DBI::dbDisconnect(conn), error = function(e)
    warning('Database disconnect failed: ', conditionMessage(e), call. = FALSE)), add = TRUE)
  row_params <- function(x) unlist(lapply(seq_len(nrow(x)), function(i)
    unname(as.list(x[i, , drop = FALSE]))), recursive = FALSE, use.names = FALSE)
  input_match <- paste(rep(paste0('(', paste(paste(db_columns, '= ?'), collapse = ' AND '), ')'),
                           nrow(input_keys)), collapse = ' OR ')
  pair_match <- paste(rep('(trial = ? AND subject = ?)', nrow(pairs)), collapse = ' OR ')
  # SQL marks existing input keys using the database's own comparison rules.
  # One SELECT captures all requested records, without a trial/subject cross product.
  query <- paste0('SELECT trial, subject, sample, replicate, ref_genome AS refGenome, mode, data_file_name, ',
                  'CASE WHEN ', input_match, ' THEN 1 ELSE 0 END AS in_input ',
                  'FROM fragments WHERE ', pair_match,
                  ' ORDER BY trial, subject, sample, replicate, ref_genome, mode')
  records <- as.data.frame(DBI::dbGetQuery(conn, query, params = c(row_params(input_keys), row_params(pairs))))
  record_keys <- get_keys(records, 'Database records')
  if(anyDuplicated(record_keys)) stop('Error - duplicate fragment keys returned by the database.', call. = FALSE)
  if(nrow(dplyr::anti_join(unique(record_keys[c('trial', 'subject')]), pairs, by = c('trial', 'subject'))))
    stop('Error - database trial/subject identifiers differ from the incoming identifiers; check spelling and capitalization.', call. = FALSE)
  records[key_columns] <- record_keys
  
  data_path <- if(is.null(args$dataPath) || identical(args$dataPath, 'none')) '/data' else args$dataPath
  pulled <- list()
  for(i in seq_len(nrow(pairs))){
    group <- records[records$trial == pairs$trial[i] & records$subject == pairs$subject[i], , drop = FALSE]
    extra <- group[group$in_input == 0L, , drop = FALSE]
    label <- paste0('trial="', pairs$trial[i], '", subject="', pairs$subject[i], '"')
    updateLog(paste0('pullDBfragments: ', label, ': ', nrow(group), ' DB record(s) available; ',
                     nrow(group) - nrow(extra), ' already represented in the incoming data.'))
    group_rows <- 0L
    if(nrow(extra)){
      if(!is.character(data_path) || length(data_path) != 1L || is.na(data_path) || !dir.exists(data_path))
        stop('Error - fragment data lake directory does not exist: ', data_path, call. = FALSE)
      if(!requireNamespace('arrow', quietly = TRUE)) stop('Error - the arrow R package is required to read fragment Parquet files.', call. = FALSE)
      for(j in seq_len(nrow(extra))){
        name <- extra$data_file_name[j]
        if(is.na(name) || !nzchar(name) || basename(name) != name || !grepl('[.]parquet$', name))
          stop('Error - invalid fragment Parquet filename for ', label, '.', call. = FALSE)
        path <- file.path(data_path, name)
        if(!file.exists(path) || dir.exists(path) || file.access(path, 4L) != 0L)
          stop('Error - fragment Parquet file is missing or unreadable: ', path, call. = FALSE)
        x <- as.data.frame(arrow::read_parquet(path))
        missing_columns <- setdiff(names(input), names(x))
        if(length(missing_columns)) stop('Error - fragment Parquet file ', name, ' lacks input columns: ',
                                         paste(missing_columns, collapse = ', '), call. = FALSE)
        keys <- get_keys(x, paste0('Parquet file ', name))
        if(any(vapply(key_columns, function(column) any(keys[[column]] != extra[[column]][j]), logical(1))))
          stop('Error - fragment Parquet contents do not match their database key: ', path, call. = FALSE)
        # Match the normalization at the start of buildStdFragments.
        x[] <- lapply(x, function(column) if(is.factor(column)) as.character(column) else column)
        x[key_columns] <- keys
        pulled[[length(pulled) + 1L]] <- x
        group_rows <- group_rows + nrow(x)
      }
    }
    updateLog(paste0('pullDBfragments: ', label, ': pulled ', ppNum(group_rows),
                     ' fragment rows from ', nrow(extra), ' DB record(s) into the analysis.'))
  }
  result <- if(length(pulled)) dplyr::bind_rows(pulled) else empty
  updateLog(paste0('pullDBfragments: total ', ppNum(nrow(result)), ' additional fragment rows pulled for ',
                   nrow(pairs), ' trial/subject group(s).'))
  result
}




# Database-only archiving; the analysis helpers and source objects are unchanged.
.multiHitDBDataLake <- "/data"

.multiHitKeys <- function(x, label){
  columns <- c("trial", "subject", "sample", "refGenome", "mode")
  if(!is.data.frame(x) || anyDuplicated(names(x)) || !all(columns %in% names(x)))
    stop("Error - ", label, " must include trial, subject, sample, refGenome and mode.", call. = FALSE)
  keys <- as.data.frame(x)[columns]
  for(column in columns){
    v <- keys[[column]]
    if(!is.atomic(v) || !is.null(dim(v)) || is.complex(v) || anyNA(v))
      stop("Error - ", label, ".", column, " contains invalid identifiers.", call. = FALSE)
    v <- as.character(v)
    if(any(!nzchar(trimws(v))) || any(v != trimws(v)))
      stop("Error - ", label, ".", column, " contains blank or whitespace-padded identifiers.", call. = FALSE)
    keys[[column]] <- v
  }
  keys
}

# Only used on the existing no-unique-position error path, with DB flags enabled.
# No rescue is possible without a unique position. Run the existing cluster
# builder on a copy, archive its summary, then let the module keep its error.
archiveMultiHitOnlyToDB <- function(frags, processedGroups){
  clusters <- data.table::rbindlist(lapply(
    split(data.table::copy(frags), by = c("trial", "subject"),
          flatten = TRUE, sorted = TRUE), function(x){
            result <- build_multiHit_clusters(x)
            if("cdhitAssignments" %in% names(result))
              result[, cdhitAssignments := NULL]
            if("splitReport" %in% names(result)) result[, splitReport := NULL]
            result
          }),
    use.names = TRUE, fill = TRUE)
  uploadMultiHitClustersToDB(clusters, processedGroups)
}

uploadMultiHitClustersToDB <- function(clusters, processedGroups){
  if(!is.list(args) || !all(c("dbConfigFile", "dbConfigID") %in% names(args)) ||
     any(vapply(args[c("dbConfigFile", "dbConfigID")], function(x)
       !is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x) || x == "none", logical(1))) ||
     !file.exists(args$dbConfigFile))
    stop("Error - uploadMultiHitClustersToDB requires args$dbConfigID and an existing args$dbConfigFile.", call. = FALSE)
  key_columns <- c("trial", "subject", "sample", "refGenome", "mode")
  db_columns <- c("trial", "subject", "sample", "ref_genome", "mode")
  keys <- .multiHitKeys(processedGroups, "Processed groups")
  if(!nrow(keys) || anyDuplicated(keys))
    stop("Error - processed groups must be non-empty and unique.", call. = FALSE)
  if(!is.data.frame(clusters) || anyDuplicated(names(clusters)))
    stop("Error - clusters must be a data frame with unique column names.", call. = FALSE)
  data <- data.table::copy(clusters)
  if(nrow(data)){
    cluster_keys <- .multiHitKeys(data, "Clusters")
    if(!"clusterID" %in% names(data) ||
       !(is.character(data$clusterID) || is.factor(data$clusterID)) ||
       anyNA(data$clusterID) || any(!nzchar(trimws(as.character(data$clusterID)))))
      stop("Error - clusters must contain non-empty clusterID values.", call. = FALSE)
    if(anyDuplicated(cbind(cluster_keys, clusterID = as.character(data$clusterID))))
      stop("Error - duplicate clusterID within a sample/genome/mode group.", call. = FALSE)
    key_ids <- function(x) do.call(paste0, lapply(x, function(v)
      paste0(nchar(v, type = "bytes"), ":", v)))
    group_id <- match(key_ids(cluster_keys), key_ids(keys))
    if(anyNA(group_id))
      stop("Error - clusters include a group absent from processedGroups.", call. = FALSE)
  } else {
    group_id <- integer()
  }
  if(!dir.exists(.multiHitDBDataLake) || !file.exists(file.path(.multiHitDBDataLake, ".inspiired")) ||
     file.access(.multiHitDBDataLake, 2L) != 0L)
    stop("Error - data lake must be writable and contain its .inspiired marker: ", .multiHitDBDataLake, call. = FALSE)
  data_lake <- normalizePath(.multiHitDBDataLake, mustWork = TRUE)
  tmp_dir <- if(is.null(args$tmpDir)) tempdir() else args$tmpDir
  if(!dir.exists(tmp_dir) || file.access(tmp_dir, 2L) != 0L)
    stop("Error - upload temporary directory is missing or not writable: ", tmp_dir, call. = FALSE)
  log_message <- function(text){
    if(is.null(args$logFile)) message(text) else updateLog(text)
  }
  
  conn <- createDBconnection()
  on.exit(tryCatch(DBI::dbDisconnect(conn), error = function(e) warning(conditionMessage(e))), add = TRUE)
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
      message("Prepared RDS files retained after an unsuccessful or unverified multihit_clusters upload: ",
              paste(published_files, collapse = ", "))
    if(locked) tryCatch(DBI::dbGetQuery(conn,
                                        "SELECT RELEASE_LOCK(CONCAT('INSPIIRED2.uploadMHC:', MD5(DATABASE()))) AS released"
    ), error = function(e) warning("Upload lock release failed: ", conditionMessage(e), call. = FALSE))
  }, add = TRUE, after = FALSE)
  
  lock <- DBI::dbGetQuery(conn,
                          "SELECT GET_LOCK(CONCAT('INSPIIRED2.uploadMHC:', MD5(DATABASE())), 30) AS acquired")
  if(nrow(lock) != 1L || !isTRUE(lock$acquired[[1L]] == 1L))
    stop("Error - could not acquire the multihit_clusters upload lock; check the selected database and retry after other uploads finish.", call. = FALSE)
  locked <- TRUE
  engine <- DBI::dbGetQuery(conn, paste(
    "SELECT ENGINE AS engine FROM information_schema.TABLES",
    "WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'multihit_clusters'"))
  if(nrow(engine) != 1L || !isTRUE(toupper(engine$engine[[1L]]) == "INNODB"))
    stop("Error - multihit_clusters must be an InnoDB table.", call. = FALSE)
  schema <- DBI::dbGetQuery(conn, paste(
    "SELECT COLUMN_NAME AS column_name, CHARACTER_MAXIMUM_LENGTH AS max_length, IS_NULLABLE AS nullable",
    "FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'multihit_clusters'"))
  required <- c(db_columns, "total_clusters", "processed_date", "data_file_name")
  if(!all(required %in% schema$column_name))
    stop("Error - the multihit_clusters table is missing required columns.", call. = FALSE)
  primary <- DBI::dbGetQuery(conn, paste(
    "SELECT COLUMN_NAME AS column_name FROM information_schema.KEY_COLUMN_USAGE",
    "WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'multihit_clusters'",
    "AND CONSTRAINT_NAME = 'PRIMARY' ORDER BY ORDINAL_POSITION"))
  if(!identical(as.character(primary$column_name), db_columns))
    stop("Error - multihit_clusters must have the five-column sample-level primary key.", call. = FALSE)
  if(!isTRUE(schema$nullable[match("data_file_name", schema$column_name)] == "YES"))
    stop("Error - multihit_clusters.data_file_name must allow NULL for zero-cluster results.", call. = FALSE)
  lengths <- setNames(as.numeric(schema$max_length), schema$column_name)
  for(i in seq_along(key_columns)){
    limit <- lengths[[db_columns[[i]]]]
    if(is.na(limit) || any(nchar(keys[[i]], type = "chars") > limit)){
      stop("Error - multihit_clusters.", db_columns[[i]], " is too short for the supplied identifiers.", call. = FALSE)
    }
  }
  if(is.na(lengths[["data_file_name"]]) || lengths[["data_file_name"]] < 36L)
    stop("Error - multihit_clusters.data_file_name must allow at least 36 characters.", call. = FALSE)
  
  key_where <- paste(paste(db_columns, "= ?"), collapse = " AND ")
  select_sql <- paste("SELECT data_file_name, total_clusters FROM multihit_clusters WHERE", key_where)
  insert_sql <- paste(
    "INSERT INTO multihit_clusters (trial, subject, sample, ref_genome, mode, total_clusters, processed_date, data_file_name)",
    "VALUES (?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP, ?)"
  )
  uploads <- vector("list", nrow(keys))
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
    local <- tempfile("multihit_clusters_", tmpdir = tmp_dir, fileext = ".rds")
    pending <- tempfile(".multihit_clusters_", tmpdir = data_lake, fileext = ".pending")
    on.exit(unlink(c(local, pending)), add = TRUE)
    saveRDS(x, local, version = 3)
    checksum <- if(file.exists(local) && file.size(local) > 0) unname(tools::md5sum(local)) else NA_character_
    if(length(checksum) != 1L || is.na(checksum) || !grepl("^[[:xdigit:]]{32}$", checksum))
      stop("Error - failed to write and checksum multihit_clusters RDS file.", call. = FALSE)
    checksum <- tolower(checksum)
    filename <- paste0(checksum, ".rds")
    path <- file.path(data_lake, filename)
    verify_file <- function(candidate){
      if(is_symlink(candidate) || !utils::file_test("-f", candidate))
        stop("Error - multihit_clusters RDS path is not a regular, non-symlink file: ", candidate, call. = FALSE)
      if(!identical(unname(tools::md5sum(candidate)), checksum))
        stop("Error - multihit_clusters RDS file does not match its checksum-derived name: ", path, call. = FALSE)
      # Equal MD5 values alone cannot rule out an actual hash collision.
      if(!same_bytes(local, candidate))
        stop("Error - MD5 collision: different RDS bytes have the same filename: ", path, call. = FALSE)
    }
    if(file.exists(path) || is_symlink(path)){
      verify_file(path)
      return(filename)
    }
    if(!isTRUE(file.copy(local, pending, overwrite = FALSE)) ||
       !identical(unname(tools::md5sum(pending)), checksum))
      stop("Error - failed to copy and verify multihit_clusters RDS file in the data lake.", call. = FALSE)
    # Both paths are in the lake. A hard link publishes the complete file
    # atomically and cannot overwrite an existing name (unlike POSIX rename).
    # If another uploader won the race, verify its file before reusing it.
    if(isTRUE(suppressWarnings(file.link(pending, path)))){
      published_files <<- c(published_files, path)
    } else if(!file.exists(path) && !is_symlink(path)){
      stop("Error - could not publish multihit_clusters RDS without overwriting; the data-lake filesystem must support hard links: ",
           path, call. = FALSE)
    }
    verify_file(path)
    filename
  }
  
  # Prepare every file before deleting any database records. Existing content
  # is verified and reused; it is never overwritten by this uploader.
  for(i in seq_len(nrow(keys))){
    indices <- which(group_id == i)
    filename <- NA_character_
    if(length(indices)){
      subset <- if(data.table::is.data.table(data)) data[indices] else data[indices, , drop = FALSE]
      filename <- stage_file(subset)
      log_message(paste0("Prepared multi-hit RDS: ", filename))
    }
    uploads[[i]] <- list(params = unname(as.list(keys[i, , drop = FALSE])),
                         total = length(indices), file = filename)
  }
  new_names <- vapply(uploads, `[[`, character(1), "file")
  old_files <- character()
  old_path <- function(name){
    if(is.na(name) || !nzchar(name)) return(NULL)
    path <- if(startsWith(name, "/")) name else file.path(data_lake, name)
    if(is_symlink(path)) stop("Error - refusing to remove an unsafe multihit_clusters file path: ", name, call. = FALSE)
    path <- normalizePath(path, mustWork = FALSE)
    if(dirname(path) != data_lake || dir.exists(path))
      stop("Error - refusing to remove an unsafe multihit_clusters file path: ", name, call. = FALSE)
    if(file.exists(path) && !grepl("[.]rds$", path))
      stop("Error - existing multihit_clusters file is not a RDS file: ", name, call. = FALSE)
    path
  }
  verify_record <- function(connection, upload){
    stored <- DBI::dbGetQuery(connection, select_sql, params = upload$params)
    valid_file <- nrow(stored) == 1L && if(is.na(upload$file))
      is.na(stored$data_file_name[[1L]]) else
        isTRUE(as.character(stored$data_file_name[[1L]]) == upload$file)
    if(!valid_file || !isTRUE(as.numeric(stored$total_clusters[[1L]]) == upload$total))
      stop("Error - multihit_clusters database verification failed for ", paste(upload$params, collapse = "/"), call. = FALSE)
  }
  
  DBI::dbBegin(conn)
  transaction_open <- TRUE
  # Delete all old rows before inserting any replacements. If two input groups
  # compare equal under the DB collation, the primary key rejects the second
  # INSERT and the entire batch rolls back. Filenames cannot detect this: an
  # unchanged re-upload legitimately has the same filename as the old record.
  for(upload in uploads){
    old <- DBI::dbGetQuery(conn, paste(select_sql, "FOR UPDATE"), params = upload$params)
    if(nrow(old) > 1L)
      stop("Error - multiple multihit_clusters records found for one primary key.", call. = FALSE)
    if(nrow(old)){
      old_path(as.character(old$data_file_name[[1L]]))
      old_files <- c(old_files, as.character(old$data_file_name[[1L]]))
    }
    deleted <- DBI::dbExecute(conn, paste("DELETE FROM multihit_clusters WHERE", key_where), params = upload$params)
    if(!isTRUE(deleted == nrow(old))) stop("Error - unexpected multihit_clusters deletion count.", call. = FALSE)
  }
  for(upload in uploads){
    inserted <- DBI::dbExecute(conn, insert_sql, params = c(upload$params, list(upload$total, upload$file)))
    if(!isTRUE(inserted == 1L)) stop("Error - failed to insert a multihit_clusters record.", call. = FALSE)
    verify_record(conn, upload)
  }
  tryCatch({
    # Older RMariaDB dbCommit() methods ignore C API errors. The SQL query path
    # checks the server's response; dbCommit() then clears driver bookkeeping.
    DBI::dbExecute(conn, "COMMIT", immediate = TRUE)
    DBI::dbCommit(conn)
  }, error = function(e) {
    stop("Error - multihit_clusters COMMIT failed; its outcome may be uncertain. Old and new RDS files were retained. ",
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
    stop("Error - multihit_clusters COMMIT completed but independent verification failed. Old and new RDS files were retained. ",
         conditionMessage(e), call. = FALSE)
  })
  verified <- TRUE
  
  for(name in unique(old_files[!is.na(old_files) & nzchar(old_files)])){
    path <- old_path(name)
    if(!file.exists(path)) next
    if(basename(path) %in% new_names) next
    aliases <- unique(c(name, basename(path), path, file.path(.multiHitDBDataLake, basename(path))))
    placeholders <- paste(rep("?", length(aliases)), collapse = ", ")
    referenced <- FALSE
    for(table in c("multihit_clusters", "sites", "fragments")){
      refs <- tryCatch(DBI::dbGetQuery(conn, paste0("SELECT data_file_name FROM ", table,
                                                    " WHERE data_file_name IN (", placeholders, ") LIMIT 1"), params = as.list(aliases)),
                       error = function(e) stop("Error - multihit_clusters upload committed, but old-file reference checking failed. Retained ",
                                                path, ": ", conditionMessage(e), call. = FALSE))
      referenced <- referenced || nrow(refs) > 0L
    }
    if(referenced){
      log_message(paste0("Retained RDS still referenced by another record: ", path))
    } else if(unlink(path) != 0L || file.exists(path)){
      stop("Error - multihit_clusters upload committed, but could not remove obsolete RDS file: ", path, call. = FALSE)
    }
  }
  log_message(paste0("Multi-hit clusters upload committed and verified: ", length(uploads), " record(s)."))
  invisible(TRUE)
}


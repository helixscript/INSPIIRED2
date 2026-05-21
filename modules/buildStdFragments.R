#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'parallel', 'data.table', 'BiocParallel', 'stringi', 'igraph')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",                     type = "character",   required = TRUE,                help = "Directory for output files")
parser$add_argument("--inputData",                     type = "character",   required = TRUE,                help = "Path to demultiplex module's rds output file.")
parser$add_argument("--softwareRoot",                  type = "character",   required = TRUE,                help = "Path to INSPIIRED2 installation.")
parser$add_argument("--threads",                       type = "integer",     default  = 30,                  help = "Number of threads to use.")
parser$add_argument("--fileTag",                       type = "character",   default  = "buildStdFragments", help = "String appended to output files in the outpt directory.")
parser$add_argument("--ramDiskPath",                   type = "character",   default  = "/dev/shm",          help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
parser$add_argument("--clusterLeaderSeqs",           action = "store_true",  default  = FALSE,               help = 'Cluster leader sequences and consider when building fragments.')
parser$add_argument("--disableBreakPointPosStd",     action = "store_true",  default  = FALSE,               help = 'Disable break point standardization.')
parser$add_argument("--disableIntSitePosStd",        action = "store_true",  default  = FALSE,               help = 'Disable intSite position standardization.')
parser$add_argument("--anchorReadClusterLen",          type = "integer",     default  = 30,                  help = 'Length of anchor read sequences to test for rearrangments.')
parser$add_argument("--anchorReadClusterMinAbundDiff", type = "integer",     default  = 5,                   help = 'When clustering anchor read sequences, min. difference between 1st and 2nd most abundant sequence clusters to pick a winner.')
parser$add_argument("--anchorReadClusterMinReadMult",  type = "integer",     default  = 10,                  help = 'When clustering anchor read sequences, multiplier for 1st and 2nd most read sequence clusters to pick a winner.')
parser$add_argument("--UMIclusterMinReadMult",         type = "integer",     default  = 5,                   help = 'When clustering UMI sequences, multiplier for 1st and 2nd most read sequence clusters to pick a winning target.')
parser$add_argument("--minReadsPerFrag",               type = "integer",     default  = 1,                   help = 'Min. number of reads to accept a fragment.')
parser$add_argument("--UMIclusteringParams",           type = "character",   default  =  "-c 0.80 -d 0 -M 0 -g 1 -r 0 -n 4 -G 1",                                       help = 'Clustering params for clustering UMIs.')
parser$add_argument("--leaderSeqClusteringParams",     type = "character",   default  =  "-c 0.87 -d 0 -M 0 -g 0 -r 0 -n 5 -G 1 -aS 0.80",                              help = 'Clustering params for clustering leader sequences.')
parser$add_argument("--multiHitclusteringParams",      type = "character",   default  =  "-c 0.87 -d 0 -M 0 -g 0 -r 0 -n 5 -G 1 -gap -5 -gap-ext -1 -aS 0.93",          help = 'Clustering params for clustering building multi-hit clusters.')
parser$add_argument("--anchorReadClusterParams",       type = "character",   default  =  "-c 0.87 -d 0 -M 0 -g 0 -r 0 -n 5 -G 1 -gap -5 -gap-ext -2 -aS 0.93 -aL 0.93", help = 'Clustering params for clustering the start of anchor read sequences.')

runModule <- function(){
  startModule()
  
  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  on.exit({
    unlink(args$tmpDir,  recursive = TRUE, force = TRUE)
    unlink(args$logDir,  recursive = TRUE, force = TRUE) 
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  updateLog('Starting buildStdFragment module.')
  
  if(! file.exists(args$inputData))  stop(paste0('Error - the input data file (', args$inputData, ') does not exist.'))
  if(file.size(args$inputData) == 0) stop(paste0('Error - the input data file (', args$inputData, ') is empty.'))
  
  frags <- setDT(readRDS(args$inputData))
 
  
  # leaderSeq clustering (optional)
  #-----------------------------------------------------------------------------
  if(args$clusterLeaderSeqs){
    orgFragRowCount <- nrow(frags)
    o <- dplyr::arrange(data.frame(table(frags$leaderSeq)), desc(Freq))
    o$n <- 1:nrow(o)
    o$readID <- paste0('s',  o$n)
    
    # Cluster unique leader sequences.
    ts <- tmpString()
    write(paste0('>', o$readID, '\n', o$Var1), file = file.path(args$ramDisk, paste0(ts, '.fasta')))
    out_prefix <- file.path(args$ramDisk, paste0(ts, "_cdhit"))
    cmd <- paste0("cd-hit-est ", args$leaderSeqClusteringParams, " -T ", args$threads, " -i ", file.path(args$ramDisk, paste0(ts, '.fasta')), " -o ", out_prefix)
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    clstr_path <- paste0(out_prefix, ".clstr")
    if(!file.exists(clstr_path)) stop(paste0('Error - cd-hit-est failed to return a clstr file.'))
    
    # Rename the clusters using table `o` so that the clusters with the highest number of reads are numbered the lowest.
    r <- parse_cdhit_clstr(clstr_path)
    
    browser()
    
    r <- left_join(r, o[, c('n', 'readID')], by = 'readID')
    k <- group_by(r, cluster_id) %>% summarise(newClusterID = paste('Cluster', min(n))) %>% ungroup()
    r <- left_join(r, k, by = 'cluster_id')
    o <- left_join(o, r[, c('readID', 'newClusterID')], by = 'readID')
    frags <- left_join(frags, o[, c('Var1', 'newClusterID')], by = c('leaderSeq' = 'Var1'))

    # Use re-named cluster ids to determine leader sequencing group numbers.
    frags$leaderSeqGroupNum <-  as.integer(str_extract(frags$newClusterID, '\\d+'))
    frags$newClusterID <- NULL
    if(nrow(frags) != orgFragRowCount) stop('Error sorting and rennanming leader sequence clusters.')
  } else {
    frags$leaderSeqGroupNum <- 1
  }
  
  
  # Build fragment ids and separate reads for position standardization.
  #-----------------------------------------------------------------------------
  frags[, fragID := paste(trial, subject, sample, replicate, fragChromosome, 
                          fragStrand, fragStart, fragEnd, leaderSeqGroupNum, 
                          UMI, sep = ":")]
  
  posFrags <- frags[frags$fragStrand == '+']
  negFrags <- frags[frags$fragStrand == '-']
  rm(frags)
  
  posSubjectFrags <- list()
  negSubjectFrags <- list()
  
  if(nrow(posFrags) > 0) posSubjectFrags <- split(posFrags, by = c('trial', 'subject', 'fragChromosome', 'leaderSeqGroupNum'), flatten = TRUE, sorted = TRUE)
  if(nrow(negFrags) > 0) negSubjectFrags <- split(negFrags, by = c('trial', 'subject', 'fragChromosome', 'leaderSeqGroupNum'), flatten = TRUE, sorted = TRUE)

  
  # Standardize intSite positions.
  #-----------------------------------------------------------------------------
  if(! args$disableIntSitePosStd){
    if(length(posSubjectFrags) > 0){
      posSubjectFrags <- lapply(posSubjectFrags, function(x){
               tab <- group_by(x, seqnames = fragChromosome, 
                                  strand = fragStrand, 
                                  start = fragStart, 
                                  end = fragEnd) %>% summarise(reads = sum(nReads), .groups = "drop") %>% arrange(start, end, reads)

               tab2 <- standardize_positions(tab, side = 'left', window = 10)
               update <- unique(data.table(fragStart = tab$start, newFragStart = tab2$start))
               
               preJoinRows <- nrow(x)
               x <- left_join(x, update, by = 'fragStart')
               if(nrow(x) != preJoinRows) stop('intSite posRepFrags join Error')
               if(any(is.na(x$newFragStart))) stop('intSite posRepFrags NA Error')
               x
      })
    }
    
    if(length(negSubjectFrags) > 0){
      negSubjectFrags <- lapply(negSubjectFrags, function(x){
        tab <- group_by(x, seqnames = fragChromosome, 
                           strand = fragStrand, 
                           start = fragStart, 
                           end = fragEnd) %>% summarise(reads = sum(nReads), .groups = "drop") %>% arrange(end, start, reads)
        tab2 <- standardize_positions(tab, side = 'right', window = 10)

        update <- unique(data.table(fragEnd = tab$end, newFragEnd = tab2$end))
        
        preJoinRows <- nrow(x)
        x <- left_join(x, update, by = 'fragEnd')
        if(nrow(x) != preJoinRows) stop('intSite negRepFrags join Error')
        if(any(is.na(x$newFragEnd))) stop('intSite posRepFrags NA Error')
        x
      })
    }
  }
  
  posSubjectFrags <- rbindlist(posSubjectFrags)
  posMaxUpdatedDist <- abs(max(posSubjectFrags$newFragStart - posSubjectFrags$fragStart))
  posPercentUpdated <- sprintf("%.2f%%", (sum(posSubjectFrags$fragStart != posSubjectFrags$newFragStart) / nrow(posSubjectFrags))*100)
  
  updateLog(paste0('Intsite positions updated for positive strand fragments. Max position shift: ', posMaxUpdatedDist, 
                   ', percent fragments updated: ', posPercentUpdated))
  
  negSubjectFrags <- rbindlist(negSubjectFrags)
  negMaxUpdatedDist <- abs(max(negSubjectFrags$newFragEnd - negSubjectFrags$fragEnd))
  negPercentUpdated <- sprintf("%.2f%%", (sum(negSubjectFrags$fragEnd != negSubjectFrags$newFragEnd) / nrow(negSubjectFrags))*100)
  
  updateLog(paste0('Intsite positions updated for negative strand fragments. Max position shift: ', negMaxUpdatedDist, 
                   ', percent fragments updated: ', negPercentUpdated))
  
  
  # Assign updated coordinates.
  posSubjectFrags$fragStart <- posSubjectFrags$newFragStart
  negSubjectFrags$fragEnd   <- negSubjectFrags$newFragEnd
  
  posSubjectFrags$newFragStart <- NULL
  negSubjectFrags$newFragEnd   <- NULL
  
  if(nrow(posFrags) > 0) posRepFrags <- split(posSubjectFrags, by = c('trial', 'subject', 'sample', 'replicate', 'fragChromosome', 'leaderSeqGroupNum'), flatten = TRUE, sorted = TRUE)
  if(nrow(negFrags) > 0) negRepFrags <- split(negSubjectFrags, by = c('trial', 'subject', 'sample', 'replicate', 'fragChromosome', 'leaderSeqGroupNum'), flatten = TRUE, sorted = TRUE)
  
  
  # Standardize break point positions.
  #-----------------------------------------------------------------------------
  if(! args$disableBreakPointPosStd){
    if(length(posRepFrags) > 0){
      posRepFrags <- lapply(posRepFrags, function(x){
        tab <- group_by(x, seqnames = fragChromosome, 
                        strand = fragStrand, 
                        start = fragStart, 
                        end = fragEnd) %>% summarise(reads = sum(nReads), .groups = "drop") %>% arrange(end, start, reads)
        tab2 <- standardize_positions(tab, side = 'right', window = 5)
        update <- unique(data.table(fragEnd = tab$end, newFragEnd = tab2$end))
        preJoinRows <- nrow(x)
        x <- left_join(x, update, by = 'fragEnd')
        if(nrow(x) != preJoinRows)   stop('breakPoint posRepFrags join Error')
        if(any(is.na(x$newFragEnd))) stop('breakPoint posRepFrags NA Error')
        x
      })
    }
    
    if(length(negRepFrags) > 0){
      negRepFrags <- lapply(negRepFrags, function(x){
        tab <- group_by(x, seqnames = fragChromosome, 
                        strand = fragStrand, 
                        start = fragStart, 
                        end = fragEnd) %>% summarise(reads = sum(nReads), .groups = "drop") %>% arrange(start, end, reads)
        tab2 <- standardize_positions(tab, side = 'left', window = 5)
        update <- unique(data.table(fragStart = tab$start, newFragStart = tab2$start))
        preJoinRows <- nrow(x)
        x <- left_join(x, update, by = 'fragStart')
        if(nrow(x) != preJoinRows) stop('negRepFrags join Error')
        if(any(is.na(x$newFragStart))) stop('breakPoint negRepFrags NA Error')
        x
      })
    }
  }
  
  posRepFrags <- rbindlist(posRepFrags)
  posMaxUpdatedDist <- abs(max(posRepFrags$newFragEnd - posRepFrags$fragEnd))
  posPercentUpdated <- sprintf("%.2f%%", (sum(posRepFrags$fragEnd != posRepFrags$newFragEnd) / nrow(posRepFrags))*100)
  
  updateLog(paste0('Intsite break points updated for positive strand fragments. Max position shift: ', posMaxUpdatedDist, 
                   ', percent fragments updated: ', posPercentUpdated))
  
  negRepFrags <- rbindlist(negRepFrags)
  negMaxUpdatedDist <- abs(max(negRepFrags$newFragStart - negRepFrags$fragStart))
  negPercentUpdated <- sprintf("%.2f%%", (sum(negRepFrags$fragStart != negRepFrags$newFragStart) / nrow(negRepFrags))*100)
  
  updateLog(paste0('Intsite break points updated for negative strand fragments. Max position shift: ', negMaxUpdatedDist, 
                   ', percent fragments updated: ', negPercentUpdated))
  
  # Assign updated coordinates.
  posRepFrags$fragEnd <- posRepFrags$newFragEnd
  negRepFrags$fragStart <- negRepFrags$newFragStart
  
  posRepFrags$newFragEnd <- NULL
  negRepFrags$newFragStart <- NULL
  
  frags <- rbindlist(list(posRepFrags, negRepFrags))
  
  
  # UMI cleanup.
  #-----------------------------------------------------------------------------
  
  updateLog('Clustering UMIs within fragment records.')

  frags[, fragID := paste(trial, subject, sample, replicate, fragChromosome, 
                          fragStrand, fragStart, fragEnd, leaderSeqGroupNum, sep = ":")]
  
  UMI_worker <- function(f) {
    results <- lapply(split(f$data, f$data$fragID), function(x) {
      
      unique_umis <- unique(x$UMI) # Ensure we only cluster unique sequences
      
      if (length(unique_umis) > 1) {
        ts <- tmpString()
        fasta_path <- file.path(args$ramDisk, paste0(ts, '.fasta'))
        out_prefix <- file.path(args$ramDisk, paste0(ts, "_cdhit"))
        
        # 1. Write ONLY unique UMIs to the fasta to prevent row explosion
        write(paste0('>', unique_umis, '\n', unique_umis), file = fasta_path)
        
        cmd <- paste0("cd-hit-est ", args$UMIclusteringParams, " -T 1 -i ", fasta_path, " -o ", out_prefix)
        system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
        
        clstr_path <- paste0(out_prefix, ".clstr")
        if (!file.exists(clstr_path)) stop("cd-hit-est failed.")
        
        # Parse and create a 1:1 mapping table: original_UMI -> newUMI
        r <- parse_cdhit_clstr(clstr_path)
        
        # Identify the most frequent UMI in each cluster to act as the 'newUMI'
        mapping <- as.data.table(r)[, .(newUMI = .SD[1, readID]), by = cluster_id]
        r_mapped <- merge(as.data.table(r), mapping, by = "cluster_id")
        
        # Final 1:1 lookup table
        lookup <- unique(r_mapped[, .(readID, newUMI)])
        
        # Join back to fragment data 
        x_dt <- as.data.table(x)
        x_dt[lookup, newUMI := i.newUMI, on = .(UMI = readID)]
        return(x_dt)
        
      } else {
        # If only one UMI, newUMI is just the original UMI
        x_dt <- as.data.table(x)
        x_dt[, newUMI := UMI]
        return(x_dt)
      }
    })
    
    rbindlist(results)
  }
  

  make_fragment_iterator <- function(dt, chunk_num_start = 0) {
    chunks <- sort(unique(dt$chunkNum))
    current_idx <- 1
    total_chunks <- length(chunks)
    chunk_count <- chunk_num_start
    
    function() {
      if (current_idx > total_chunks) return(NULL)
      this_chunk_id <- chunks[current_idx]
      chunk_data <- dt[chunkNum == this_chunk_id]
      chunk_count <<- chunk_count + 1
      current_idx <<- current_idx + 1
      
      list(
        data = chunk_data,
        chunk_num = chunk_count,
        is_last = (current_idx > total_chunks)
      )
    }
  }
  
  unique_frags <- frags[, .(fragID = unique(fragID))]
  unique_frags[, chunkNum := cut(seq_len(.N), breaks = args$threads, labels = FALSE)]
  frags[unique_frags, chunkNum := i.chunkNum, on = "fragID"]
  frag_iterator <- make_fragment_iterator(frags)
  
  param <- MulticoreParam(workers = args$threads)
  frags <- rbindlist(bpiterate(ITER = frag_iterator, FUN = UMI_worker, BPPARAM = param))
  
  bpstop(param)
  closeAllConnections()
  
  updateLog(paste0(sprintf("%.2f%%", (sum(frags$UMI != frags$newUMI) / nrow(frags))*100), ' UMI sequences updated.'))
  frags$UMI <- frags$newUMI
  frags$newUMI <- NULL
  
  
  # Add UMIs back to fragments records.
  frags[, fragID := paste(trial, subject, sample, replicate, fragChromosome, 
                          fragStrand, fragStart, fragEnd, leaderSeqGroupNum, 
                          UMI, sep = ":")]
  
  # Create position ids.
  frags$posid <- paste0(frags$fragChromosome, frags$fragStrand, ifelse(frags$fragStrand == '+', frags$fragStart, frags$fragEnd))

  
  # Identify uniquely called reads.
  #-----------------------------------------------------------------------------
  
  updateLog('Identifying uniquely called positions.')
  
  u <- group_by(frags, readID) %>% 
       mutate(nPosIDs = n_distinct(posid)) %>% 
       ungroup() %>% filter(nPosIDs == 1)  %>% 
       pull(readID)
  
  frags_uniqPosIDs <- frags[frags$readID %in% u,] 
  frags_multPosIDs <- frags[! frags$readID %in% u,]
  
  updateLog(paste0(sprintf("%.2f%%", (n_distinct(frags_uniqPosIDs$readID) / n_distinct(frags$readID))*100), ' of fragment reads mapped uniquely to the genome.'))
  
  # Do not continue unless we have at least one uniquely called fragment.
  if(nrow(frags_uniqPosIDs) == 0) quitOnErorr('Error - No unique position remain after filtering.')
  
  
  # Correct for instances where a read maps to more than fragment but all fragments 
  # have the same integration position. These are instances of fuzzy break points 
  # and here we select the shortest fragments lengths.
  #-----------------------------------------------------------------------------
  z <- frags_uniqPosIDs$readID[duplicated(frags_uniqPosIDs$readID)]
  
  if(length(z) > 0){
    updateLog('Correcting fuzzy break points.')
    
    a <- subset(frags_uniqPosIDs, readID %in% z)
    b <- subset(frags_uniqPosIDs, ! readID %in% z)
    
    a2 <- rbindlist(lapply(split(a, a$readID), function(x){
      if(n_distinct(x$fragID) > 1 & n_distinct(x$posid) == 1){
        if(x$fragStrand[1] == '+'){
          i <- which(x$fragEnd == min(x$fragEnd))[1]
          x <- x[i,]
        } else {
          i <- which(x$fragStart == max(x$fragStart))[1]
          x <- x[i,]
        }
      } else {
        # Only one frag.
        x <- x[1,]
      }
      x
    }))
    
    frags_uniqPosIDs <- bind_rows(a2, b)
    
    invisible(rm(a, b, a2))
    invisible(gc())
  }
  
  
  # Multi-hit read rescue.
  #-----------------------------------------------------------------------------
  updateLog('Rescuing multihit reads using list on uniquely called positions.')
  
  if(nrow(frags_multPosIDs) > 0){
    frags_multPosIDs <- rbindlist(lapply(split(frags_multPosIDs, paste0(frags_multPosIDs$trial, frags_multPosIDs$subject)), function(x){
      unique_subject_posids <- unique(subset(frags_uniqPosIDs, trial == x$trial[1] & subject == x$subject[1])$posid)
      
      rbindlist(lapply(split(x, x$readID), function(xx){
        xx$rescue <- FALSE
        i <- which(xx$posid %in% unique_subject_posids)
        if(length(i) == 1) xx[i,]$rescue <- TRUE
        xx
      }))
    }))
    
    if(any(frags_multPosIDs$rescue == TRUE)){
      r <- frags_multPosIDs[frags_multPosIDs$rescue == TRUE]
      updateLog(paste0(ppNum(nrow(r)), ' reads rescued from multihit read table.'))
      frags_uniqPosIDs <- rbindlist(list(frags_uniqPosIDs, r[, 'rescue' := NULL]))
      frags_multPosIDs <- frags_multPosIDs[! frags_multPosIDs$readID %in% frags_uniqPosIDs$readID]
    }
  }
  
  saveRDS(frags_multPosIDs, file = file.path(args$outputDir, paste0(args$fileTag, '_multHitFrags.rds')))
  
  
  # Build multihit clusters
  #-----------------------------------------------------------------------------
  updateLog('Building multi-hit clusters.')
  multiHit_clusters <- rbindlist(lapply(split(frags_multPosIDs, by = c('trial', 'subject'), flatten = TRUE, sorted = TRUE), build_multiHit_clusters))
  saveRDS(multiHit_clusters, file = file.path(args$outputDir, paste0(args$fileTag, '_multHitClusters.rds')))


  
  # Anchor read cluster filter
  #-----------------------------------------------------------------------------
  updateLog('Cluster the beginnings of anchor read sequences.')
  anchorReadClusterDecisionTable <- tibble()
  

  frags_uniqPosIDs <- bind_rows(lapply(split(frags_uniqPosIDs, by = c('trial', 'subject', 'sample'), flatten = TRUE, sorted = TRUE), function(s){ 
    # Sort fragment records so that fragments likely to contribute to high abund / high read count fragments apart first.
    s <- group_by(s, posid) %>%
         mutate(potentialAbund = n_distinct(abs(fragEnd - fragStart))) %>%
         ungroup() %>%
         group_by(anchor_seq) %>%
         mutate(anchorReadSeqReads = n()) %>%
         ungroup() %>%
         mutate(anchorReadSeqLen = nchar(anchor_seq)) %>%
         arrange(desc(potentialAbund), desc(anchorReadSeqReads), desc(anchorReadSeqLen)) %>%
         select(-potentialAbund, -anchorReadSeqReads, -anchorReadSeqLen)
    
    # Determine anchor read test sequences, write to disk, cluster, and parse.
    s$testSeq <- substr(s$anchor_seq, 1, args$anchorReadClusterLen)
    ts <- tmpString()
    write(paste0('>', s$readID, '\n', s$testSeq), file = file.path(args$ramDisk, paste0(ts, '.fasta')))
    out_prefix <- file.path(args$ramDisk, paste0(ts, "_cdhit"))
    cmd <- paste0("cd-hit-est ", args$anchorReadClusterParams, " -T ", args$threads, " -i ", file.path(args$ramDisk, paste0(ts, '.fasta')), " -o ", out_prefix)
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    clstr_path <- paste0(out_prefix, ".clstr")
    if(! file.exists(clstr_path)) stop(paste0('Error - cd-hit-est failed to return a clstr file: ', file.exists(clstr_path)))
    
    r <- parse_cdhit_clstr(clstr_path)
    s <- left_join(s, dplyr::select(r, readID, cluster_id), by = 'readID')
    s$remove <- FALSE
    clusterNum <- 1
    
    s <- bind_rows(lapply(split(s, s$cluster_id), function(x){
      seqs <- s[s$readID %in% x$readID,]$testSeq
      repSeq <- names(sort(table(as.character(seqs)), decreasing = TRUE)[1])
      clusterNum <<- clusterNum + 1
      
      if(n_distinct(x$posid) > 1){
        
        x$anchorReadCluster <- TRUE
        
        z  <- group_by(x, trial, subject, sample, posid) %>% 
              summarise(frags = n_distinct(fragEnd - fragStart + 1), reads = n(), readIDs = list(readID), .groups = "drop") %>% 
              ungroup() %>% 
              arrange(desc(frags), desc(reads)) %>%
              mutate(clusterNum = clusterNum - 1,
                     clusterRepSeq = repSeq,
                     selected = FALSE,
                     remove = TRUE,
                     criteria = NA)
        
        z2 <- arrange(z, desc(reads), desc(frags)) # Re-sort table for read based decisions.

        if((z[1,]$frags - z[2,]$frags) >= args$anchorReadClusterMinAbundDiff){
          x$remove <- ifelse(x$posid == z[1,]$posid, FALSE, TRUE)
          z$remove <- ifelse(z$posid == z[1,]$posid, FALSE, TRUE)
          z$criteria <- ifelse(z$posid == z[1,]$posid, 'fragment counts', NA)
          z$selected <- z$posid == z[1,]$posid
          anchorReadClusterDecisionTable <<- bind_rows(anchorReadClusterDecisionTable, z)
        } else if(z2[1,]$reads >= (z2[2,]$reads * args$anchorReadClusterMinReadMult)){
          x$remove  <- ifelse(x$posid == z2[1,]$posid, FALSE, TRUE)
          z2$remove <- ifelse(z2$posid == z2[1,]$posid, FALSE, TRUE)
          z2$criteria <- ifelse(z2$posid == z2[1,]$posid, 'read counts', NA)
          z2$selected <- z2$posid == z2[1,]$posid
          anchorReadClusterDecisionTable <<- bind_rows(anchorReadClusterDecisionTable, z2)
        } else {
          x$remove <- TRUE
          z$remove <- TRUE
          z$selected <- FALSE
          anchorReadClusterDecisionTable <<- bind_rows(anchorReadClusterDecisionTable, z)
        }
      } else {
        x$anchorReadCluster <- FALSE
        x$remove <- FALSE
      }
      
      x
    }))
    
    s
  }))

  saveRDS(anchorReadClusterDecisionTable, file.path(args$outputDir, paste0(args$fileTag, '_anchorReadClusters.rds')))
  updateLog(paste0('Removing ', n_distinct(subset(frags_uniqPosIDs, remove == TRUE)$readID), ' reads due to not being the clear choice.'))
  frags_uniqPosIDs <- subset(frags_uniqPosIDs, remove == FALSE)
  frags_uniqPosIDs <- setDT(dplyr::select(frags_uniqPosIDs, -testSeq, -cluster_id, -remove))
  
  
  # UMI filter
  #-----------------------------------------------------------------------------
  
  UMIclusterDecisionTable <- tibble()
  updateLog('Ensuring that each fragment record is associated with a single UMI sequence.')

  frags_uniqPosIDs <- rbindlist(lapply(split(frags_uniqPosIDs, by = c('trial', 'subject', 'sample', 'UMI'), flatten = TRUE, sorted = TRUE), function(x){
    if(n_distinct(x$posid) > 1){
      
      z  <- group_by(x, trial, subject, sample, UMI, posid) %>% 
            summarise(frags = n_distinct(fragEnd - fragStart + 1), reads = n(), readIDs = list(readID), .groups = "drop") %>% 
            ungroup() %>% 
            arrange(desc(reads)) %>%
            mutate(selected = FALSE,
                   remove = TRUE,
                   criteria = NA)
      
      if(any(z$frags) > 1){
        x$remove <- TRUE
        z$remove <- TRUE
        z$selected <- FALSE
      } else if(z[1,]$reads >= (z[2,]$reads * args$UMIclusterMinReadMult)){
        x$remove <- ifelse(x$posid == z[1,]$posid, FALSE, TRUE)
        z$remove <- ifelse(z$posid == z[1,]$posid, FALSE, TRUE)
        z$criteria <- ifelse(z$posid == z[1,]$posid, 'read counts', NA)
        z$selected <- z$posid == z[1,]$posid
      } else {
        x$remove <- TRUE
        z$remove <- TRUE
        z$selected <- FALSE
      }
      
      UMIclusterDecisionTable <<- bind_rows(UMIclusterDecisionTable, z)
    } else {
      x$remove <- FALSE
    }
    
    x
  }))
  
  saveRDS(UMIclusterDecisionTable, file.path(args$outputDir, paste0(args$fileTag, '_UMIclusters.rds')))
  updateLog(paste0('UMI filter - removing ', n_distinct(subset(frags_uniqPosIDs, remove == TRUE)$readID), ' reads due to not being the clear choice.'))
  frags_uniqPosIDs <- subset(frags_uniqPosIDs, remove == FALSE)
  frags_uniqPosIDs <- dplyr::select(frags_uniqPosIDs, -remove)
  
  
  # Count the number of reads associated with each fragment.
  # fragments with more than one read, i > 1, need additional processing.
  frags <- group_by(data.frame(frags_uniqPosIDs), fragID) %>% mutate(i = n()) %>% ungroup()
  
  a <- subset(frags, i == 1)   
  b <- subset(frags, i > 1)    
  
  commonLeaderSeq <- function(x){
    dplyr::group_by(x, leaderSeq) %>% 
    dplyr::summarise(nReads = n()) %>% 
    dplyr::ungroup() %>%
    dplyr::slice_max(nReads, with_ties = FALSE) %>%
    dplyr::pull(leaderSeq)
  }
  
  if(nrow(b) > 0){
    b2 <- rbindlist(lapply(split(b, b$fragID), function(x){
      totalReads <- sum(x$nReads)
      
      if(totalReads < args$minReadsPerFrag) return(data.table())
      
      x$repLeaderSeq <- commonLeaderSeq(x)
      
      readList <- x$readID
      UMIs <- unique(x$UMI)
      
      x <- x[1,]
      
      x$reads <- totalReads
      x$readIDs   <- list(sort(readList))
      x$UMIs <- list(sort(UMIs))
      x
    }))
  } else {
    b2 <- data.table()
  }
  
  if(nrow(a) > 0){
    a2 <- dplyr::group_by(a, fragID) %>%
          dplyr::mutate(repLeaderSeq = leaderSeq[1],
                        reads = nReads, 
                        readIDs = list(sort(readID)),
                        UMIs = list(sort(UMI))) %>%
         dplyr::slice(1) %>%
         dplyr::ungroup() %>%
         dplyr::filter(reads >= args$minReadsPerFrag) %>%
         as.data.table()
  } else {
    a2 <- data.table()
  } 
  
  frags <- rbindlist(list(a2, b2))
  
  frags$replicate <- as.integer(frags$replicate)
  frags$fragStart <- as.integer(frags$fragStart)
  frags$fragEnd   <- as.integer(frags$fragEnd)
  
  frags <- frags[, .(mode, refGenome, trial, subject, sample, replicate, UMI, posid, reads, repLeaderSeq, fragChromosome, fragStrand, fragStart, fragEnd, anchorReadCluster, readIDs, UMIs, leaderSeqGroupNum)]
  
  frags$clusterLeaderSeqs <- as.factor(args$clusterLeaderSeqs)
    
  saveRDS(frags, file.path(args$outputDir, paste0(args$fileTag, '.rds')))
  updateLog('buildFragments module completed.')
  write(date(), file.path(args$outputDir, paste0(args$fileTag, '.done')))
}

args <- parser$parse_args()
source(file.path(args$softwareRoot, 'lib.R'))

tryCatch({
  runModule()
  q(status = 0)
}, error = function(e) {
  message("Caught error: ", e$message)
  q(status = 1)
})

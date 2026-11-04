#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'data.table', 'stringi', 'RMariaDB')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",                    type = "character",     required = TRUE,         help = "Directory for output files.")
parser$add_argument("--inputData",                    type = "character",     required = TRUE,         help = "Path to buildStdFragments output rds file.")
parser$add_argument("--softwareRoot",                 type = "character",     required = TRUE,         help = "Path to INSPIIRED2 installation.")
parser$add_argument("--threads",                      type = "integer",       default = 50,            help = "Number of threads to use.")
parser$add_argument("--fileTag",                      type = "character",     default = "buildSites",  help = "String appended to output files in the outpt directory.")
parser$add_argument("--ramDiskPath",                  type = "character",     default = "/dev/shm",    help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
parser$add_argument("--disableDualDetect",            action = "store_true",  default = FALSE,         help = "Diable the merging of U5 and U3 samples into dual-detection sites.")
parser$add_argument("--disableOrientationCorrection", action = "store_true",  default = FALSE,         help = "Disable the changing of fragment strands to reflect integrated vector orientation.")
parser$add_argument("--dualDetectWidth",              type = "integer",       default = 6,             help = "Radius for searching for dual-detections.")
parser$add_argument("--integraseCorrectionDist",      type = "integer",       default = 2,             help = "Integrase correction value (NT) to account for gDNA duplication caused by integration.")
parser$add_argument("--sumSonicBreaksWithin",         type = "character",     default = "replicates",  help = "Sum sonic breaks within either 'replicates' (default) or within sample 'samples'.") 
parser$add_argument("--leadSeqClusteringParms",       type = "character",     default = "-c 0.90 -n 5 -G 0 -aS 0.95 -gap -2 -gap-ext -1 -d 0 -M 0", help = "CLustering parameters used to determine representative leaders sequence.")
parser$add_argument("--dbConfigFile",                 type = "character",     default = 'none',        help = "Path to db credential file.")
parser$add_argument("--dbConfigID",                   type = "character",     default = 'none',        help = "DB credential block identifier in db credential file.")

# Dev notes. 
# --threads not implemented yet.

runModule <- function(){
  doneFile <- file.path(args$outputDir, paste0(args$fileTag, '.done'))
  if(file.exists(doneFile) && unlink(doneFile) != 0) stop('Error - could not remove stale completion marker: ', doneFile)
  
  on.exit({
    if(!is.null(args$dbConn)) tryCatch({
      if(DBI::dbIsValid(args$dbConn))
        DBI::dbDisconnect(args$dbConn)
    }, error = function(e) {
      warning("Database disconnect failed: ",
              conditionMessage(e), call. = FALSE)
    })
  }, add = TRUE)
  
  # Validate the input's rebuild provenance before opening any DB connection.
  args$dbConn <<- NULL
  startModule(connectDB = FALSE)
  provenanceArgs <- args[names(args) != 'dbConn']
  
  yaml::write_yaml(args[names(args) != "dbConn"], file.path(args$outputDir, paste0(args$fileTag, ".yml")))
  
  on.exit({
    unlink(args$tmpDir,  recursive = TRUE, force = TRUE)
    unlink(args$logDir,  recursive = TRUE, force = TRUE) 
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  updateLog('Starting buildSites module.')
  
  if(xor(args$dbConfigFile != "none", args$dbConfigID != "none")){
    msg <- "Error - supply both --dbConfigFile and --dbConfigID, or neither."
    updateLog(msg)
    stop(msg, call. = FALSE)
  }
  useDB <- args$dbConfigFile != 'none' && args$dbConfigID != 'none'
  
  validateSiteRebuild <- function(x, provenance){
    instruction <- paste0(' Re-run buildStdFragments with --pullSubjectFragments ',
                          'and database credentials, then use its complete RDS output.')
    if(!is.data.frame(x) || !nrow(x) || anyDuplicated(names(x)))
      stop('Error - database site rebuilding requires a non-empty fragment table with unique column names.',
           instruction, call. = FALSE)
    flag <- x[['pullSubjectFragments']]
    if(!(is.factor(flag) || is.logical(flag)) || !is.null(dim(flag)) ||
       length(flag) != nrow(x) || anyNA(flag) ||
       !all(as.character(flag) == 'TRUE'))
      stop('Error - database site rebuilding requires pullSubjectFragments=TRUE for every input fragment; missing, FALSE, NA, or mixed values are not accepted.',
           instruction, call. = FALSE)
    required <- c('version', 'pullSubjectFragments', 'fragmentRows', 'groups')
    if(!is.list(provenance) || anyDuplicated(names(provenance)) ||
       !all(required %in% names(provenance)) ||
       !identical(provenance$version, 1L) ||
       !identical(provenance$pullSubjectFragments, TRUE))
      stop('Error - database site rebuilding requires valid, completed inspiired2.siteRebuild provenance.',
           instruction, call. = FALSE)
    rows <- provenance$fragmentRows
    if(!is.numeric(rows) || !is.null(dim(rows)) || length(rows) != 1L ||
       is.na(rows) || !is.finite(rows) || rows != nrow(x))
      stop('Error - the fragment row count does not match its site rebuild provenance.',
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
    observed <- unique(normalizeKeys(x, 'Input fragments'))
    if(nrow(dplyr::anti_join(observed, groups, by = keyColumns)))
      stop('Error - input fragments include a sample/genome group absent from the site rebuild provenance.',
           instruction, call. = FALSE)
    groups
  }
  
  addSiteRebuildProvenance <- function(x, groups){
    # Never pass fragment-level or stale site metadata into a new site result.
    attr(x, 'inspiired2.siteRebuild') <- NULL
    if(is.null(groups)) return(x)
    flag <- x[['pullSubjectFragments']]
    if(!is.data.frame(x) || !nrow(x) ||
       !(is.factor(flag) || is.logical(flag)) || !is.null(dim(flag)) ||
       length(flag) != nrow(x) || anyNA(flag) ||
       !all(as.character(flag) == 'TRUE'))
      stop('Error - cannot attach completed rebuild provenance to an empty or unverified site result.',
           call. = FALSE)
    attr(x, 'inspiired2.siteRebuild') <- list(
      version = 2L,
      pullSubjectFragments = TRUE,
      siteRows = nrow(x),
      groups = groups
    )
    x
  }
  
  if(! args$sumSonicBreaksWithin %in% c('samples', 'replicates')) stop("Error - the flag --sumSonicBreaksWithin must be set to with 'samples' or 'replicates'.")
  if(! file.exists(args$inputData))  stop(paste0('Error - the input data file (', args$inputData, ') does not exist.'))
  if(file.size(args$inputData) == 0) stop(paste0('Error - the input data file (', args$inputData, ') is empty.'))
  
  # Read in standardized fragments.
  frags <- readRDS(args$inputData)
  inputProvenance <- provenanceFromObject(frags,
                                          source = list(path = args$inputData, format = 'rds'))
  
  # Validate even for file-only runs so a later annotation module can perform
  # the first database upload using the original complete replacement scope.
  siteRebuildGroups <- tryCatch({
    validateSiteRebuild(frags, attr(frags, 'inspiired2.siteRebuild', exact = TRUE))
  }, error = function(e){
    if(useDB){
      updateLog(conditionMessage(e))
      stop(conditionMessage(e), call. = FALSE)
    }
    updateLog(paste0('Continuing without database rebuild provenance; this site output cannot be uploaded to the database. ',
                     sub('^Error - ', '', conditionMessage(e))))
    NULL
  })
  if(!is.null(siteRebuildGroups)){
    updateLog(paste0('Validated completed archive retrieval for ',
                     ppNum(nrow(siteRebuildGroups)), ' sample/genome group(s).'))
  }
  
  if(nrow(frags) == 0){
    msg <- 'Error - fragment input data file has zero rows.'
    updateLog(msg)
    stop(msg)
  }
  
  # Preserve a uniform valid flag when present; older non-DB inputs need none.
  sitePullStatus <- NULL
  flag <- frags[['pullSubjectFragments']]
  if((is.factor(flag) || is.logical(flag)) && is.null(dim(flag)) &&
     length(flag) == nrow(frags) && !anyNA(flag)){
    values <- unique(as.character(flag))
    if(length(values) == 1L && values %in% c('FALSE', 'TRUE'))
      sitePullStatus <- values
  }
  
  frags$mode <- as.character(frags$mode)
  
  # Define fragment widths.
  frags$fragWidths <- frags$fragEnd - frags$fragStart + 1
  
  if(!args$disableDualDetect && all(c('U3', 'U5') %in% frags$mode)){
    updateLog('Searching for dual detections.')
    if(length(args$dualDetectWidth) != 1L || is.na(args$dualDetectWidth) || args$dualDetectWidth < 0L)
      stop('Error - dualDetectWidth must be a non-negative integer.')
    if(anyNA(frags$posid)) stop('Error - missing posids cannot be evaluated for dual detection.')
    if('.dualRowID' %in% names(frags)) stop('Error - reserved column .dualRowID is already present.')
    
    # Build all candidate relationships before changing any fragments. Only
    # reciprocal one-to-one U3/U5 relationships are accepted.
    originalRowCount <- nrow(frags)
    frags[, .dualRowID := .I]
    dualGroupCols <- c('trial', 'subject', 'sample', 'refGenome')
    
    frags <- rbindlist(lapply(split(frags, by = dualGroupCols, keep.by = TRUE,
                                    flatten = TRUE, sorted = TRUE, drop = TRUE), function(x){
                                      if(!all(c('U3', 'U5') %in% x$mode)) return(x)
                                      
                                      u3 <- x[mode == 'U3']
                                      u5 <- x[mode == 'U5']
                                      
                                      candidatePairs <- rbindlist(lapply(sort(unique(as.character(u3$posid))), function(u3Posid){
                                        parts <- unlist(strsplit(u3Posid, '[\\+\\-]'))
                                        if(length(parts) != 2L || is.na(suppressWarnings(as.integer(parts[2]))))
                                          stop('Error - unable to parse U3 posid during dual-detection search: ', u3Posid)
                                        
                                        chrom <- parts[1]
                                        pos <- as.integer(parts[2])
                                        strand <- str_extract(u3Posid, '[\\+\\-]')
                                        alts <- paste0(chrom, ifelse(strand == '+', '-', '+'),
                                                       (pos - args$dualDetectWidth):(pos + args$dualDetectWidth))
                                        candidateU5 <- sort(unique(as.character(u5$posid[u5$posid %in% alts])))
                                        
                                        if(!length(candidateU5)) return(NULL)
                                        data.table(u3Posid = u3Posid, u5Posid = candidateU5,
                                                   chromosome = chrom, u3Strand = strand)
                                      }), use.names = TRUE)
                                      
                                      if(!nrow(candidatePairs)) return(x)
                                      
                                      # For each U3 site, this counts how many distinct U5 sites it could match.
                                      candidatePairs <- unique(candidatePairs, by = c('u3Posid', 'u5Posid'))
                                      candidatePairs[, nU5forU3 := uniqueN(u5Posid), by = u3Posid]
                                      candidatePairs[, nU3forU5 := uniqueN(u3Posid), by = u5Posid]
                                      setorder(candidatePairs, u3Posid, u5Posid)
                                      
                                      # Remove ambiguous pairings.
                                      rejected <- candidatePairs[nU5forU3 != 1L | nU3forU5 != 1L]
                                      
                                      
                                      # Warn about ambiguous pairings.
                                      if(nrow(rejected)){
                                        groupLabel <- paste(x$trial[1], x$subject[1], x$sample[1], x$refGenome[1], sep = '/')
                                        updateLog(paste0('Warning - ambiguous dual-detection candidates in ', groupLabel,
                                                         ': U3 sites [', paste(sort(unique(rejected$u3Posid)), collapse = ', '),
                                                         ']; U5 sites [', paste(sort(unique(rejected$u5Posid)), collapse = ', '),
                                                         ']. Leaving all involved sites unmerged.'))
                                      }
                                      
                                      accepted <- candidatePairs[nU5forU3 == 1L & nU3forU5 == 1L]
                                      
                                      if(!nrow(accepted)) return(x)
                                      
                                      if(anyDuplicated(accepted$u3Posid) || anyDuplicated(accepted$u5Posid))
                                        stop('Error - a dual-detection site was assigned to more than one accepted pair.')
                                      
                                      # Accepted data table stores unambiguous rationale dual detection pairs 
                                      # u3Posid           u5Posid chromosome u3Strand nU5forU3 nU3forU5
                                      # 1:  chr17+7683122   chr17-7683126       chr17        +         1         1
                                      # 2: chr3-170572672  chr3+170572668        chr3        -         1         1
                                      
                                      for(j in seq_len(nrow(accepted))){
                                        targetU3 <- accepted$u3Posid[j]
                                        targetU5 <- accepted$u5Posid[j]
                                        u3Rows <- u3[posid == targetU3, .dualRowID]  # Return row ids of targeted U3 candidates
                                        u5Rows <- u5[posid == targetU5, .dualRowID]  # Return row ids of targeted U5 candidates
                                        pairRows <- c(u3Rows, u5Rows)
                                        
                                        if(!length(u3Rows) || !length(u5Rows))
                                          stop('Error - accepted dual-detection pair has no associated fragments.')
                                        
                                        f1 <- x[.dualRowID %in% u3Rows] # Retrieve U3 records as f1
                                        f2 <- x[.dualRowID %in% u5Rows] # Retrieve U5 records as f2
                                        
                                        updateLog(paste0('   Processing U3 posid ', targetU3, ' as a dual detection with ',
                                                         nrow(f2), ' U5 fragments at ', targetU5, '.'))
                                        
                                        # Apply the integrase correction factor to positive and negative strand positions in the pairing.
                                        i <- which(x$.dualRowID %in% pairRows & x$fragStrand == '+')
                                        if(length(i)) x[i, posid := unlist(lapply(strsplit(as.character(posid), '[\\+\\-]', perl = TRUE),
                                                                                  function(z) paste0(z[1], '+', as.integer(z[2]) + args$integraseCorrectionDist)))]
                                        
                                        i <- which(x$.dualRowID %in% pairRows & x$fragStrand == '-')
                                        if(length(i)) x[i, posid := unlist(lapply(strsplit(as.character(posid), '[\\+\\-]', perl = TRUE),
                                                                                  function(z) paste0(z[1], '-', as.integer(z[2]) - args$integraseCorrectionDist)))]
                                        
                                        # Record the most common leader sequences and report both in repLeaseSeq.
                                        # Set mode and assign a reserved leaderSeqGroupNum, 0, reserved for dual detections.
                                        i <- which(x$.dualRowID %in% pairRows)
                                        u3Leader <- names(sort(table(f1$repLeaderSeq), decreasing = TRUE))[1]
                                        u5Leader <- names(sort(table(f2$repLeaderSeq), decreasing = TRUE))[1]
                                        x[i, repLeaderSeq := paste0(u3Leader, '/', u5Leader)]
                                        x[i, mode := 'dual detect']
                                        x[i, leaderSeqGroupNum := 0]
                                        
                                        # Identify th mos common, corrected position, and correct all fragments in the pair to that position.
                                        # Assign the the orientation strand used in the posid by U3 alignment.
                                        pos <- names(sort(table(sub('[\\+\\-]', '', stringr::str_extract(as.character(x[i, posid]), '[\\+\\-]\\d+'))), decreasing = TRUE))[1]
                                        
                                        if(length(pos) != 1L || is.na(pos))
                                          stop('Error - unable to determine a common position for accepted dual-detection pair.')
                                        
                                        if(accepted$u3Strand[j] == '-'){
                                          x[i, fragStrand := '+']
                                          x[i, posid := paste0(accepted$chromosome[j], '+', pos)]
                                        } else {
                                          x[i, fragStrand := '-']
                                          x[i, posid := paste0(accepted$chromosome[j], '-', pos)]
                                        }
                                      }
                                      
                                      x
                                    }), use.names = TRUE, fill = FALSE)
    
    setorder(frags, .dualRowID)
    if(nrow(frags) != originalRowCount || anyDuplicated(frags$.dualRowID) ||
       !identical(frags$.dualRowID, seq_len(originalRowCount)))
      stop('Error - dual-detection processing changed the fragment row set.')
    frags[, .dualRowID := NULL]
  }
  
  if(! args$disableOrientationCorrection & ('U5' %in% frags$mode | 'U3' %in% frags$mode)){
    updateLog('Updating strandedness of U5 and U3 intSite calls.')
    
    frags <- bind_rows(lapply(split(frags, paste(frags$trial, frags$subject, frags$sample, frags$refGenome)), function(x){
      a <- x[x$mode == "dual detect",]
      b <- x[x$mode != "dual detect",]
      
      if(nrow(b)){
        
        # Shift positions to reflect duplication caused by integrase.
        b1 <- subset(b, fragStrand == '+')
        if(nrow(b1) > 0) b1$posid <- unlist(lapply(strsplit(b1$posid, '[\\+\\-]', perl = TRUE), function(x) paste0(x[1], '+', as.integer(x[2]) + args$integraseCorrectionDist)))
        
        b2 <- subset(b, fragStrand == '-')
        if(nrow(b2) > 0) b2$posid <- unlist(lapply(strsplit(b2$posid, '[\\+\\-]', perl = TRUE), function(x) paste0(x[1], '-', as.integer(x[2]) - args$integraseCorrectionDist)))
        
        b <- bind_rows(b1, b2)
        rm(b1, b2)
        
        updatePosIdStrand <- function(x, s){
          o <- unlist(strsplit(x, '[\\+\\-]'))
          paste0(o[1], s, o[2])
        }
        
        # Change strand to reflect orientation. 
        b1 <- subset(b, fragStrand == '+' & grepl('U3', b$mode))
        b2 <- subset(b, fragStrand == '-' & grepl('U3', b$mode))
        b3 <- subset(b, fragStrand == '+' & grepl('U5', b$mode))
        b4 <- subset(b, fragStrand == '-' & grepl('U5', b$mode))
        
        if(nrow(b1) > 0) b1$posid <- sapply(b1$posid, updatePosIdStrand, '-')
        if(nrow(b2) > 0) b2$posid <- sapply(b2$posid, updatePosIdStrand, '+')
        if(nrow(b3) > 0) b3$posid <- sapply(b3$posid, updatePosIdStrand, '+')
        if(nrow(b4) > 0) b4$posid <- sapply(b4$posid, updatePosIdStrand, '-')
        
        b <- bind_rows(b1, b2, b3, b4)
      }
      
      bind_rows(a, b)
    }))
  }
  
  # At this point, now that we're done parsing position ids, we can add leaderSeq
  # identifiers if more than one leaderSeqGroupNum is present. 
  # 
  #  This feature is disabled in this version of the software.
  #  if(n_distinct(frags$leaderSeqGroupNum) > 1) frags$posid <- paste0(frags$posid, '.', frags$leaderSeqGroupNum)
  
  consensusLeaderSeq <- function(x){
    tab <- dplyr::group_by(x, repLeaderSeq) %>% 
      dplyr::summarise(nWidths = n_distinct(fragWidths), nReads = sum(reads)) %>% 
      dplyr::ungroup() %>%
      dplyr::arrange(desc(nWidths), desc(nReads))
    as.character(tab[1, 'repLeaderSeq'])
  }
  
  clusterSeqs <- function(seqs){
    if(length(unique(seqs)) > 1){
      ts <- tmpString()
      write(paste0('>', paste0('s', 1:length(seqs)), '\n', seqs), file = file.path(args$ramDisk, paste0(ts, '.fasta')))
      out_prefix <- file.path(args$ramDisk, paste0(ts, "_cdhit"))
      cmd <- paste0("cd-hit-est ", args$leadSeqClusteringParms, " -T 1 -i ", file.path(args$ramDisk, paste0(ts, '.fasta')), " -o ", out_prefix)
      status <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
      requireCommandSuccess(status, "cd-hit-est")
      clstr_path <- paste0(out_prefix, ".clstr")
      if(!file.exists(clstr_path)) stop(paste0('Error - cd-hit-est failed to return a clstr file.'))
      return(parse_cdhit_clstr(clstr_path))
    } else {
      return(data.table(readID = 's1', cluster_id = 'Cluster 0', is_rep = TRUE, cluster_size = 1))
    }
  }
  
  collapsePrepMetadata <- function(x){
    x <- as.character(x)
    x <- sort(unique(x[!is.na(x) & nzchar(x)]))
    if(length(x)) paste(x, collapse = ";") else NA_character_
  }
  
  frags <- group_by(frags, trial, subject, sample, mode, refGenome, posid) %>%
    mutate(g = cur_group_id()) %>%
    ungroup() %>%
    data.table()
  
  
  updateLog('Gather fragments into intSite events.')
  sites <- bind_rows(lapply(split(frags, frags$g), function(x){
    # Loop through replicates for this site defined by 'g'
    r <- bind_cols(lapply(min(frags$replicate):max(frags$replicate), function(r){
      
      b <- tibble(UMIs = NA, sonicLengths = NA, reads = NA, repLeaderSeq = NA)
      o <- x[x$replicate == r,]
      
      if(nrow(o) >= 1){
        b$UMIs <- n_distinct(unlist(o$UMIs))
        b$sonicLengths <- n_distinct(o$fragWidths)
        b$reads <- sum(o$reads)
        b$repLeaderSeq <- consensusLeaderSeq(o)
      } 
      
      names(b) <- paste0('rep', r, '-', names(b))
      b
    }))
    
    bind_cols(tibble(trial = x$trial[1], 
                     subject = x$subject[1], 
                     sample = x$sample[1],
                     refGenome = x$refGenome[1],
                     mode = x$mode[1],
                     leaderSeqHMM = collapsePrepMetadata(x$leaderSeqHMM),
                     vectorFastaFile = collapsePrepMetadata(x$vectorFastaFile),
                     posid = x$posid[1],
                     UMIs = n_distinct(unlist(x$UMIs)),
                     sonicLengths = ifelse(args$sumSonicBreaksWithin == 'replicates',
                                           sum(r[, grepl('sonicLengths', names(r))], na.rm = TRUE),  
                                           n_distinct(x$fragWidths)),
                     reads = sum(x$reads),
                     repLeaderSeq = consensusLeaderSeq(x),
                     repLeaderSeqClusters = n_distinct(clusterSeqs(unique(x$repLeaderSeq))$cluster_id),
                     nRepsObs = sum(! is.na(unlist(r[, which(grepl('reads', names(r)))])))), r)
  })) %>% arrange(desc(sonicLengths))
  
  # Set nRepsObs to NA for dual detections since these have values of 1 after moving dual detection to rep-0.
  sites[sites$mode == 'dual detect',]$nRepsObs <- NA
  
  updateLog('Collapsing replicate level sites into sample level records.')
  sites <- group_by(sites, trial, subject, sample, refGenome) %>%
    mutate(
      sampleAbund = sum(sonicLengths),
      percentSampleRelAbund =
        round((sonicLengths / sampleAbund) * 100, 2),
      .after = nRepsObs
    ) %>%
    ungroup() %>%
    select(-sampleAbund)
  
  # Change UMI counts of 1 to NA if UMIs were not captured.
  umiValues <- unlist(frags$UMIs, use.names = FALSE)
  if(length(umiValues) > 0L && isTRUE(all(umiValues == "AAAAAAAAAAAA"))){
    umiCols <- grep("^(UMIs|rep[0-9]+-UMIs)$", names(sites), value = TRUE)
    for(nm in umiCols) sites[[nm]] <- NA_integer_
  }
  
  if(!is.null(sitePullStatus))
    sites$pullSubjectFragments <- factor(rep(sitePullStatus, nrow(sites)),
                                         levels = c('FALSE', 'TRUE'))
  sites <- addSiteRebuildProvenance(sites, siteRebuildGroups)
  if(!is.null(siteRebuildGroups))
    updateLog(paste0('Preserved complete database replacement scope for ',
                     ppNum(nrow(siteRebuildGroups)), ' sample/genome group(s) in the site output.'))
  
  updateLog('Sample level site summary:')
  ts <- paste0(base::format(Sys.time(), "%m.%d.%Y"), ' [', timeElapsedString(), "]")
  siteSummary <- group_by(sites, trial, subject, sample, refGenome) %>% 
    summarise(nSites = n_distinct(posid), .groups = 'drop') %>% 
    ungroup() %>%
    mutate(timeStamp = ts, .before = trial) %>%
    mutate(across(everything(), as.character))
  siteSummary <- rbind(names(siteSummary), siteSummary)
  siteSummary[1,1] <- ts
  write.table(siteSummary, file =  args$logFile, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE,  append = TRUE)
  
  sites <- provenanceAttach(sites, provenanceRecord('buildSites',
                                                    parents = list(inputProvenance), inputFiles = args$inputData,
                                                    details = list(siteRows = nrow(sites),
                                                                   sampleGenomeScope = siteRebuildGroups),
                                                    argsSnapshot = provenanceArgs))
  updateLog('Preserved input processing history and recorded buildSites parameters.')
  
  if(useDB){
    updateLog('Database site rebuild beginning; replacing all modes for the validated sample/genome groups.')
    uploadSitesToDB(sites, processedGroups = siteRebuildGroups)
  }
  
  saveRDS(sites, file.path(args$outputDir, paste0(args$fileTag, '.rds')))
  updateLog('buildSites module completed.')
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

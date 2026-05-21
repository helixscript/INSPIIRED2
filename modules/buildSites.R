#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'data.table', 'stringi')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",               type = "character",     required = TRUE,         help = "Directory for output files")
parser$add_argument("--inputData",               type = "character",     required = TRUE,         help = "Path to demultiplex module's rds output file.")
parser$add_argument("--softwareRoot",            type = "character",     required = TRUE,         help = "Path to INSPIIRED2 installation.")
parser$add_argument("--threads",                 type = "integer",       default = 50,            help = "Number of threads to use.")
parser$add_argument("--fileTag",                 type = "character",     default = "buildSites",  help = "String appended to output files in the outpt directory.")
parser$add_argument("--ramDiskPath",             type = "character",     default = "/dev/shm",    help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
parser$add_argument("--disableDualDetect",       action = "store_true",  default = FALSE,         help = "Diable the merging of U5 and U3 samples into dual-detection sites.")
parser$add_argument("--dualDetectWidth",         type = "integer",       default = 6,             help = "Radius for searching for dual-detections.")
parser$add_argument("--integraseCorrectionDist", type = "integer",       default = 2,             help = "Integrase correction value (NT) to account for gDNA duplication caused by integration.")
parser$add_argument("--sumSonicBreaksWithin",    type = "character",     default = "samples",     help = "Sum sonic breaks within either 'samples' (default) or within sample 'replicates'.") 
parser$add_argument("--leadSeqClusteringParms",  type = "character",     default = "-c 0.90 -n 5 -G 0 -aS 0.95 -gap -2 -gap-ext -1 -d 0 -M 0", help = "CLustering parameters used to determine representative leaders sequence.")

runModule <- function(){
  startModule()
  
  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  on.exit({
    unlink(args$tmpDir,  recursive = TRUE, force = TRUE)
    unlink(args$logDir,  recursive = TRUE, force = TRUE) 
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  
  updateLog('Starting buildSites module.')
  if(! args$sumSonicBreaksWithin %in% c('samples', 'replicates')) stop("Error - the flag --sumSonicBreaksWithin must be set to with 'samples' or 'replicates'.")
  
  if(! file.exists(args$inputData))  stop(paste0('Error - the input data file (', args$inputData, ') does not exist.'))
  if(file.size(args$inputData) == 0) stop(paste0('Error - the input data file (', args$inputData, ') is empty.'))
  
  frags <- readRDS(args$inputData)
  
  # Define fragment widths.
  frags$fragWidths <- frags$fragEnd - frags$fragStart + 1
  
  # Define fragment IDs.
  frags[, fragID := paste(trial, subject, sample, replicate, fragChromosome, 
                          fragStrand, fragStart, fragEnd, sep = ":")]
  
  if(! args$disableDualDetect & 'IN_u5' %in% frags$mode & 'IN_u3' %in% frags$mode){
    updateLog('Searching for dual detections.')
    
    # Loop through IN_u3 fragments and search for close by IN_u5 fragments oriented in the opposite direction. 
    # For each identified pair, switch mode to 'dual detect', assign all fragments the U3 strand,
    # and assign a common posid which will be centered between the U3 and U5 posids.
    # Track which fragments have been merged so that closely spaced fragments are not merged more than once.
    
    invisible(lapply(split(frags, paste(frags$trial, frags$subject, frags$sample)), function(x){
      processed_fragments <- data.table()
      
      if('IN_u5' %in% x$mode & 'IN_u3' %in% x$mode){
        u3 <- x[x$mode == 'IN_u3']
        u5 <- x[x$mode == 'IN_u5']
      
        invisible(lapply(unique(u3$posid), function(u3_posid){
          parts  <- unlist(strsplit(u3_posid, '[\\+\\-]'))
          chrom  <- parts[1]
          pos    <- as.integer(parts[2])
          strand <- str_extract(u3_posid, '[\\+\\-]')
          
          alts <- paste0(chrom, ifelse(strand == '+', '-', '+'), (pos - args$dualDetectWidth):(pos + args$dualDetectWidth))
          
          z <- subset(u5, posid %in% alts)
          
          if(nrow(z) > 0){
            # Retrieve fragments for both sites.
            f1 <- subset(u3, posid == u3_posid & ! fragID %in% processed_fragments$fragID)
            f2 <- subset(u5, posid %in% alts & ! fragID %in% processed_fragments$fragID)
            
            if(nrow(f1) == 0 | nrow(f2) == 0) return()
            
            updateLog(paste0('   Processing U3 posid ', u3_posid, ' as a dual detection with ', nrow(f2), ' U5 fragments.'))
            
            # Records processed u5 fragments 
            i <- which(frags$fragID %in% c(f1$fragID, f2$fragID))
            processed_fragments <<- bind_rows(processed_fragments, frags[i,])
            
            # Apply positive strand position corrections.
            i <- which(frags$fragID %in% c(f1$fragID, f2$fragID) & frags$fragStrand == '+')
            frags[i,]$posid <<- unlist(lapply(strsplit(frags[i,]$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '+', as.integer(x[2]) + args$integraseCorrectionDist)))
            
            # Apply negative strand position corrections.
            i <- which(frags$fragID %in% c(f1$fragID, f2$fragID) & frags$fragStrand == '-')
            frags[i,]$posid <<- unlist(lapply(strsplit(frags[i,]$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '-', as.integer(x[2]) - args$integraseCorrectionDist)))
            
            # Create combined repLeaderSeq string.
            i <- which(frags$fragID %in% c(f1$fragID, f2$fragID))
            frags[i,]$repLeaderSeq <<- paste0(names(sort(table(f1$repLeaderSeq), decreasing = TRUE))[1], '/', names(sort(table(f2$repLeaderSeq), decreasing = TRUE))[1])
            
            # Set new mode.
            frags[i,]$mode <<- 'dual detect'
            
            # Find the most common fragment intSite position from combined fragments.
            pos <- names(sort(table(sub('[\\+\\-]', '', stringr::str_extract(frags[i,]$posid, '[\\+\\-]\\d+'))), decreasing = TRUE))[1]
            
            # Let a leaderSeqGroupNum value of zero represent dual-detections.
            frags[i,]$leaderSeqGroupNum <<- 0

            if(strand == '-'){
              frags[i,]$fragStrand <<- '+'                  # Set the u3 frag strands to positive to reflect correct orientation. U5 posid already '+'.
              frags[i,]$posid  <<- paste0(chrom, '+', pos)  # Set the u3 frag posids to the u5 posid which is positive causing its fragments to merge with U3 fragments.
            } else {
              frags[i,]$fragStrand <<- '-'                  # Set the u3 frag strands to negative to reflect reverse orientation. U5 posid already '-'.
              frags[i,]$posid  <<- paste0(chrom, '-', pos)  # Set the u3 frag posids to the u5 posid which is negative causing its fragments to merge with U3 fragments.
            }
          }
        }))
      }
    }))
  }
    
  if('IN_u5' %in% frags$mode | 'IN_u3' %in% frags$mode){
    updateLog('Updating strandedness of U5 and U3 intSite calls.')
    
    frags <- bind_rows(lapply(split(frags, paste(frags$trial, frags$subject, frags$sample)), function(x){
      a <- subset(frags, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1] & mode == 'dual detect')
      b <- subset(frags, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1] & mode != 'dual detect')
      
      if(nrow(b)){
        
        # Shift positions to reflect duplication caused by integrase.
        b1 <- subset(b, fragStrand == '+')
        if(nrow(b1) > 0) b1$posid <- unlist(lapply(strsplit(b1$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '+', as.integer(x[2]) + args$integraseCorrectionDist)))
        
        b2 <- subset(b, fragStrand == '-')
        if(nrow(b2) > 0) b2$posid <- unlist(lapply(strsplit(b2$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '-', as.integer(x[2]) - args$integraseCorrectionDist)))
        
        b <- bind_rows(b1, b2)
        rm(b1, b2)
        
        updatePosIdStrand <- function(x, s){
          o <- unlist(strsplit(x, '[\\+\\-]'))
          paste0(o[1], s, o[2])
        }
        
        # Change strand to reflect orientation. 
        b1 <- subset(b, fragStrand == '+' & grepl('IN_u3', b$mode))
        b2 <- subset(b, fragStrand == '-' & grepl('IN_u3', b$mode))
        b3 <- subset(b, fragStrand == '+' & grepl('IN_u5', b$mode))
        b4 <- subset(b, fragStrand == '-' & grepl('IN_u5', b$mode))
        
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
  
  if(as.logical(frags$clusterLeaderSeqs[1]) == TRUE) frags$posid <- paste0(frags$posid, '.', frags$leaderSeqGroupNum)
  
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
      system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
      clstr_path <- paste0(out_prefix, ".clstr")
      if(!file.exists(clstr_path)) stop(paste0('Error - cd-hit-est failed to return a clstr file.'))
      return(parse_cdhit_clstr(clstr_path))
    } else {
      return(data.table(readID = 's1', cluster_id = 'Cluster 0', is_rep = TRUE, cluster_size = 1))
    }
  }
  
  # Create a sample + posid grouping vector.
  frags <- group_by(frags, trial, subject, sample, posid) %>% mutate(g = cur_group_id()) %>% ungroup() %>% data.table()
  
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
                    b$repLeaderSeq <- consensusLeaderSeq(x)
                } 
      
                names(b) <- paste0('rep', r, '-', names(b))
                b
             }))
      
             bind_cols(tibble(trial = x$trial[1], 
                              subject = x$subject[1], 
                              sample = x$sample[1],
                              refGenome = x$refGenome[1],
                              mode = x$mode[1],
                              posid = x$posid[1],
                              UMIs = n_distinct(unlist(x$UMIs)),
                              sonicLengths = ifelse(args$sumSonicBreaksWithin == 'replicates',
                                                    sum(r[, grepl('sonicLengths', names(r))], na.rm = TRUE),  
                                                    n_distinct(x$fragWidths)),
                              reads = sum(x$reads),
                              repLeaderSeq = consensusLeaderSeq(x),
                              repLeaderSeqClusters = n_distinct(clusterSeqs(unique(x$repLeaderSeq))$cluster_id),
                              nRepsObs = sum(! is.na(unlist(r[, which(grepl('reads', names(r)))]))),
                              vector = x$vectorFastaFile[1]), r)
           })) %>% arrange(desc(sonicLengths))
  
  # Set nRepsObs to NA for dual detections since these have values of 1 after moving dual detection to rep-0.
  sites[sites$mode == 'dual detect',]$nRepsObs <- NA
  
  updateLog('Collapsing replicate level sites into sample level records.')
  sites <- group_by(sites, trial, subject, sample) %>%
           mutate(sampleAbund = sum(sonicLengths)) %>%
           ungroup() %>%
           group_by(posid) %>%
           mutate(percentSampleRelAbund = round((sonicLengths/sampleAbund[1]) * 100, 2), .after = 'nRepsObs') %>%
           ungroup() %>%
           select(-sampleAbund)
  
  updateLog('Sample level site summary:')
  ts <- paste0(base::format(Sys.time(), "%m.%d.%Y"), ' [', timeElapsedString(), "]")
  siteSummary <- group_by(sites, trial, subject, sample) %>% 
                 summarise(nSites = n_distinct(posid), .groups = 'drop') %>% 
                 ungroup() %>%
                 mutate(timeStamp = ts, .before = trial) %>%
                 mutate(across(everything(), as.character))
  siteSummary <- rbind(names(siteSummary), siteSummary)
  siteSummary[1,1] <- ts
  write.table(siteSummary, file =  args$logFile, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE,  append = TRUE)
  
  saveRDS(sites, file.path(args$outputDir, paste0(args$fileTag, '.rds')))
  updateLog('buildSites module completed.')
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

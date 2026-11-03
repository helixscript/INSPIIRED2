#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'data.table', 'GenomicRanges', 'RMariaDB')) suppressPackageStartupMessages(library(p, character.only = TRUE))
if(!requireNamespace("GenomeInfoDb", quietly = TRUE)) stop("Error - required R package GenomeInfoDb is not installed.")

parser <- ArgumentParser()
parser$add_argument("--outputDir",     type = "character",  required = TRUE,          help = "Directory for output files.")
parser$add_argument("--inputData",     type = "character",  required = TRUE,          help = "Path to rds output containing posid and refGenome columns")
parser$add_argument("--softwareRoot",  type = "character",  required = TRUE,          help = "Path to INSPIIRED2 installation.")
parser$add_argument("--threads",       type = "integer",    default = 50,             help = "Number of threads to use.")
parser$add_argument("--fileTag",       type = "character",  default = "nearestGenes", help = "String appended to output files in the outpt directory.")
parser$add_argument("--ramDiskPath",   type = "character",  default = "/dev/shm",     help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
parser$add_argument("--dbConfigFile",  type = "character",  default = 'none',         help = "Path to db credential file.")
parser$add_argument("--dbConfigID",    type = "character",  default = 'none',         help = "DB credential block identifier in db credential file.")

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
  
  # Validate rebuild provenance before opening any database connection.
  args$dbConn <<- NULL
  startModule(connectDB = FALSE)
  
  yaml::write_yaml(args[names(args) != "dbConn"], file.path(args$outputDir, paste0(args$fileTag, ".yml")))
  
  on.exit({
    unlink(args$tmpDir, recursive = TRUE, force = TRUE)
    unlink(args$logDir, recursive = TRUE, force = TRUE)
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  updateLog('Starting nearestGenes module.')
  
  if(xor(args$dbConfigFile != "none", args$dbConfigID != "none")){
    msg <- "Error - supply both --dbConfigFile and --dbConfigID, or neither."
    updateLog(msg)
    stop(msg, call. = FALSE)
  }
  useDB <- args$dbConfigFile != 'none' && args$dbConfigID != 'none'
  
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
  
  resource_overlay()
  
  if(!file.exists(args$inputData)) stop(paste0('Error - the input data file (', args$inputData, ') does not exist.'))
  if(file.size(args$inputData) == 0) stop(paste0('Error - the input data file (', args$inputData, ') is empty.'))
  
  d <- readRDS(args$inputData)
  
  siteRebuildProvenance <- attr(d, 'inspiired2.siteRebuild', exact = TRUE)
  siteRebuildGroups <- tryCatch(
    validateSiteRebuild(d, siteRebuildProvenance),
    error = function(e){
      if(useDB){
        updateLog(conditionMessage(e))
        stop(e)
      }
      updateLog('No valid site-rebuild provenance; annotating without database update eligibility.')
      NULL
    }
  )
  if(is.null(siteRebuildGroups)){
    siteRebuildProvenance <- NULL
    attr(d, 'inspiired2.siteRebuild') <- NULL
  } else {
    updateLog(paste0('Validated completed archive retrieval for ',
                     ppNum(nrow(siteRebuildGroups)), ' sample/genome group(s).'))
  }
  
  if(!is.data.frame(d))
    stop('Error - input data must be a table.', call. = FALSE)
  if(!nrow(d))
    stop('Error - input contains no integration sites to annotate.', call. = FALSE)
  if(!'refGenome' %in% names(d) || anyNA(d$refGenome) || any(!nzchar(trimws(as.character(d$refGenome)))))
    stop('Error - input contains missing or blank refGenome values.', call. = FALSE)
  
  o <- rbindlist(lapply(split(d, d$refGenome, drop = TRUE), function(x){
    tu <- readRDS(file.path(args$softwareRoot, 'data', 'genomeAnnotations', paste0(x$refGenome[1], '.TUs.rds')))
    ex <- readRDS(file.path(args$softwareRoot, 'data', 'genomeAnnotations', paste0(x$refGenome[1], '.exons.rds')))
    
    # i <- sub('\\.\\d+', '', unique(x$posid))
    i <- unique(x$posid)
    
    g <- makeGRangesFromDataFrame(tibble(
      seqnames = unlist(lapply(str_split(i, '[\\+\\-]'), '[', 1)),
      start = unlist(lapply(str_split(i, '[\\+\\-]'), '[', 2)),
      end = start,
      strand = str_extract(i, '[\\+\\-]')
    ))
    
    g$posid <- paste0(seqnames(g), strand(g), start(g))
    g$refGenome <- x$refGenome[1]
    
    o <- data.frame(suppressWarnings(GenomicRanges::distanceToNearest(g, tu, select = 'all', ignore.strand = TRUE)))
    e <- data.frame(suppressWarnings(GenomicRanges::distanceToNearest(g, ex, select = 'all', ignore.strand = TRUE)))
    
    r <- unlist(GenomicRanges::GRangesList(lapply(seq_along(g), function(xx){
      gg <- g[xx]
      oo <- unique(o[o$queryHits == xx, , drop = FALSE])
      
      gg$inGene <- IRanges::overlapsAny(gg, tu, ignore.strand = TRUE)
      gg$inExon <- IRanges::overlapsAny(gg, ex, ignore.strand = TRUE)
      gg$nearestGene <- NA_character_
      gg$nearestGeneStrand <- NA_character_
      gg$nearestGeneDist <- NA_integer_
      gg$beforeNearestGene <- NA
      
      if(nrow(oo) > 0L){
        hits <- distinct(data.frame(tu[oo$subjectHits, ])[, c("name2", "strand")])
        
        gg$nearestGene <- paste(hits$name2, collapse = ",")
        gg$nearestGeneStrand <- paste(hits$strand, collapse = ",")
        gg$nearestGeneDist <- as.integer(min(oo$distance) + as.integer(!gg$inGene))
        gg$beforeNearestGene <- start(gg) < min(start(tu[oo$subjectHits]))
      }
      
      gg
    })))
    
    as.data.table(data.frame(r))
  }))
  
  o <- o[,6:length(o)]
  if(anyDuplicated(o[, .(refGenome, posid)])) stop('Error - nearestGenes produced duplicate annotations for the same refGenome and posid.')
  
  inputRows <- nrow(d)
  d <- dplyr::left_join(d, o, by = c('refGenome', 'posid'))
  if(nrow(d) != inputRows) stop('Error - nearestGenes annotation join changed the input row count.')
  
  d <- dplyr::relocate(d, inGene, .after = posid)
  d <- dplyr::relocate(d, inExon, .after = inGene)
  d <- dplyr::relocate(d, nearestGene, .after = inExon)
  d <- dplyr::relocate(d, nearestGeneDist, .after = nearestGene)
  d <- dplyr::relocate(d, nearestGeneStrand, .after = nearestGeneDist)
  d <- dplyr::relocate(d, beforeNearestGene, .after = nearestGeneStrand)
  
  # Keep the original scope, including processed groups with no surviving sites.
  # Revalidate after annotation so a changed flag or row count cannot reach upload.
  if(!is.null(siteRebuildProvenance))
    validateSiteRebuild(d, siteRebuildProvenance)
  attr(d, 'inspiired2.siteRebuild') <- siteRebuildProvenance
  
  if(useDB){
    updateLog('Database site rebuild beginning; replacing all modes for the validated sample/genome groups.')
    uploadSitesToDB(d, processedGroups = siteRebuildGroups)
  }
  
  saveRDS(d, file.path(args$outputDir, paste0(args$fileTag, '.rds')))
  updateLog('Completed nearestGenes module.')
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
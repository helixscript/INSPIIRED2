#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'data.table', 'GenomicRanges', 'readr', 'DBI')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",     type = "character", required = TRUE,                     help = "Directory for output files.")
parser$add_argument("--inputData",     type = "character", required = TRUE,                     help = "Path to the multihit cluster RDS file produced by buildStdFragments.")
parser$add_argument("--softwareRoot",  type = "character", required = TRUE,                     help = "Path to INSPIIRED2 installation.")
parser$add_argument("--threads",       type = "integer",   default = 50,                        help = "Number of threads to use.")
parser$add_argument("--fileTag",       type = "character", default = "annotateMultiHitRepeats", help = "Prefix for output files in the output directory.")
parser$add_argument("--ramDiskPath",   type = "character", default = "/dev/shm",                help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
parser$add_argument("--dbConfigFile",  type = "character", default = "none",                    help = "Path to database credentials. With --dbConfigID, replaces stored multihit clusters for the input trial/subject/sample/refGenome/mode groups.")
parser$add_argument("--dbConfigID",    type = "character", default = "none",                    help = "Database credential block identifier; requires --dbConfigFile.")

runModule <- function(){
  doneFile <- file.path(args$outputDir, paste0(args$fileTag, '.done'))
  if(file.exists(doneFile) && unlink(doneFile) != 0) stop('Error - could not remove stale completion marker: ', doneFile)
  
  if(xor(args$dbConfigFile != 'none', args$dbConfigID != 'none'))
    stop('Error - supply both --dbConfigFile and --dbConfigID, or neither.', call. = FALSE)
  
  startModule(connectDB = FALSE)
  # As in buildStdFragments, the uploader manages its own short-lived connections.
  useDB <- args$dbConfigFile != 'none' && args$dbConfigID != 'none'
  options(scipen = 999)
  
  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  on.exit({
    unlink(args$tmpDir,  recursive = TRUE, force = TRUE)
    unlink(args$logDir,  recursive = TRUE, force = TRUE)
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  updateLog('Starting annotateMultiHitRepeats module.')
  
  # Expand scientific notation without changing vector order, names, or length.
  # Other malformed IDs are rejected before constructing genomic ranges.
  expandPosids <- function(posids){
    if(!is.character(posids)) stop('Error - posids must be character vectors.', call. = FALSE)
    
    pattern <- '^[^+-]+[+-]([0-9]+([.][0-9]*)?|[.][0-9]+)[eE][+-]?[0-9]+$'
    i <- grep(pattern, posids)
    if(length(i)){
      prefix <- sub('^([^+-]+[+-]).*$', '\\1', posids[i])
      position <- as.numeric(sub('^[^+-]+[+-]', '', posids[i]))
      if(any(!is.finite(position) | position != trunc(position)))
        stop('Error - scientific-notation coordinates must represent finite integers.', call. = FALSE)
      posids[i] <- paste0(prefix, sprintf('%.0f', position))
    }
    
    posids
  }
  
  annotatePositions <- function(posids, genome, softwareRoot){
    invalid <- is.na(posids) | !grepl('^[^+-]+[+-][0-9]+$', posids)
    if(any(invalid))
      stop('Error - invalid posid(s): ', paste(head(posids[invalid], 5L), collapse = ', '), call. = FALSE)
    
    positions <- as.numeric(sub('^[^+-]+[+-]', '', posids))
    if(any(!is.finite(positions) | positions < 1 | positions > .Machine$integer.max))
      stop('Error - posid coordinates must be positive integers within the supported genomic range.', call. = FALSE)
    
    repeatFile <- file.path(softwareRoot, 'data', 'genomeAnnotations', paste0(genome, '.repeatTable.gz'))
    if(!file.exists(repeatFile)) stop('Error - repeat annotation file does not exist: ', repeatFile, call. = FALSE)
    
    tab <- as.data.table(read_tsv(repeatFile, col_types = cols(.default = col_character()),
                                  show_col_types = FALSE, progress = FALSE))
    required <- c('query_seq', 'query_start', 'query_end', 'strand', 'repeat_name', 'repeat_class')
    if(!all(required %in% names(tab)))
      stop('Error - repeat annotation file is missing required columns: ', repeatFile, call. = FALSE)
    tab <- tab[, ..required]
    
    # Match annotateRepeats.R: RepeatMasker coordinates are 1-based and inclusive;
    # C denotes the negative strand, and overlaps ignore strand.
    tab$strand <- sub('C', '-', tab$strand, fixed = TRUE)
    tab <- tab[strand %in% c('+', '-')]
    tab$query_start <- suppressWarnings(as.numeric(tab$query_start))
    tab$query_end <- suppressWarnings(as.numeric(tab$query_end))
    if(any(!is.finite(tab$query_start) | !is.finite(tab$query_end) |
           tab$query_start < 1 | tab$query_end < tab$query_start |
           tab$query_end > .Machine$integer.max |
           tab$query_start != trunc(tab$query_start) | tab$query_end != trunc(tab$query_end)))
      stop('Error - repeat annotation file contains invalid coordinates: ', repeatFile, call. = FALSE)
    for(column in c('query_seq', 'repeat_name', 'repeat_class')){
      if(anyNA(tab[[column]]) || any(!nzchar(trimws(tab[[column]]))))
        stop('Error - repeat annotation file contains missing or blank ', column, ' values: ', repeatFile, call. = FALSE)
    }
    
    lookup <- data.table(posid = posids,
                         repeatName = rep('none', length(posids)),
                         repeatClass = rep('none', length(posids)),
                         classes = rep(list('none'), length(posids)))
    if(!nrow(tab)) return(lookup)
    
    repeats <- makeGRangesFromDataFrame(tab, seqnames.field = 'query_seq',
                                        start.field = 'query_start', end.field = 'query_end',
                                        strand.field = 'strand', keep.extra.columns = TRUE)
    sites <- GRanges(seqnames = sub('[+-][0-9]+$', '', posids),
                     ranges = IRanges::IRanges(start = as.integer(positions), width = 1L),
                     strand = sub('^[^+-]+([+-])[0-9]+$', '\\1', posids))
    overlaps <- findOverlaps(sites, repeats, ignore.strand = TRUE)
    
    if(length(overlaps)){
      hits <- unique(data.table(index = queryHits(overlaps),
                                repeatName = tab$repeat_name[subjectHits(overlaps)],
                                repeatClass = tab$repeat_class[subjectHits(overlaps)]))
      annotations <- hits[, .(repeatName = paste(repeatName, collapse = ','),
                              repeatClass = paste(repeatClass, collapse = ','),
                              classes = list(unique(repeatClass))), by = index]
      lookup$repeatName[annotations$index] <- annotations$repeatName
      lookup$repeatClass[annotations$index] <- annotations$repeatClass
      lookup$classes[annotations$index] <- annotations$classes
    }
    
    lookup
  }
  
  summarizeClasses <- function(classes){
    # Each position contributes once to each distinct overlapping class.
    # "none" participates in the same way as a repeat class.
    counts <- table(unlist(lapply(classes, unique), use.names = FALSE))
    maximum <- max(counts)
    leaders <- sort(names(counts)[counts == maximum], method = 'radix')
    
    list(dominant = if(length(leaders) == 1L) leaders else 'tie',
         percentage = as.numeric(maximum) / length(classes),
         tied = if(length(leaders) > 1L) leaders else character())
  }
  
  resource_overlay()
  
  if(!file.exists(args$inputData)) stop('Error - the input data file does not exist: ', args$inputData, call. = FALSE)
  if(file.size(args$inputData) == 0) stop('Error - the input data file is empty: ', args$inputData, call. = FALSE)
  
  o <- readRDS(args$inputData)
  if(!is.data.frame(o)) stop('Error - input data must be a multihit cluster table.', call. = FALSE)
  if(!all(c('refGenome', 'posids') %in% names(o)))
    stop('Error - input data must contain refGenome and posids columns.', call. = FALSE)
  setDT(o)
  
  if(anyNA(o$refGenome) || any(!nzchar(trimws(as.character(o$refGenome)))))
    stop('Error - input contains missing or blank refGenome values.', call. = FALSE)
  if(!is.list(o$posids) || !all(vapply(o$posids, is.character, logical(1))))
    stop('Error - posids must be a list of character vectors.', call. = FALSE)
  if(any(lengths(o$posids) == 0L))
    stop('Error - every multihit cluster must contain at least one posid.', call. = FALSE)
  
  # Replace complete sample/genome/mode groups represented in this input.
  if(useDB) dbMultiHitGroups <- unique(.multiHitKeys(o, 'Input clusters'))
  
  # Correct IDs before either constructing the lookup or matching back to rows.
  # Keep nested node IDs consistent with the corrected cluster IDs.
  correctedPosids <- lapply(o$posids, expandPosids)
  nCorrected <- sum(vapply(seq_len(nrow(o)), function(i){
    sum(o$posids[[i]] != correctedPosids[[i]], na.rm = TRUE)
  }, numeric(1)))
  o$posids <- correctedPosids
  
  nNodeCorrected <- 0L
  if('nodeSonicLengths' %in% names(o)){
    if(!is.list(o$nodeSonicLengths))
      stop('Error - nodeSonicLengths must be a list of tables.', call. = FALSE)
    for(i in seq_len(nrow(o))){
      nodeTable <- o$nodeSonicLengths[[i]]
      if(!is.data.frame(nodeTable) || !'posid' %in% names(nodeTable))
        stop('Error - each nodeSonicLengths table must contain a posid column.', call. = FALSE)
      corrected <- expandPosids(nodeTable$posid)
      if(!identical(nodeTable$posid, corrected)){
        nNodeCorrected <- nNodeCorrected + sum(nodeTable$posid != corrected, na.rm = TRUE)
        nodeTable <- copy(nodeTable)
        nodeTable$posid <- corrected
        o$nodeSonicLengths[[i]] <- nodeTable
      }
    }
  }
  updateLog(paste0('Expanded scientific notation in ', ppNum(nCorrected),
                   ' posids entries and ', ppNum(nNodeCorrected), ' nodeSonicLengths posid entries.'))
  
  repeatNames <- rep(list(character()), nrow(o))
  repeatClasses <- rep(list(character()), nrow(o))
  dominantRepeatClass <- rep(NA_character_, nrow(o))
  # Retain the requested column spelling and the fraction scale (0 to 1).
  dominatRepeatClassPercentage <- rep(NA_real_, nrow(o))
  tiedRepeatClasses <- rep(list(character()), nrow(o))
  
  for(genome in unique(as.character(o$refGenome))){
    rows <- which(as.character(o$refGenome) == genome)
    allPosids <- unlist(o$posids[rows], use.names = FALSE)
    updateLog(paste0('Annotating ', ppNum(length(rows)), ' clusters and ',
                     ppNum(length(allPosids)), ' positions for ', genome, '.'))
    
    lookup <- annotatePositions(unique(allPosids), genome, args$softwareRoot)
    allIndices <- match(allPosids, lookup$posid)
    if(anyNA(allIndices)) stop('Error - an input posid is missing from the annotation lookup.', call. = FALSE)
    
    offset <- 0L
    for(row in rows){
      size <- length(o$posids[[row]])
      indices <- allIndices[seq.int(offset + 1L, offset + size)]
      offset <- offset + size
      
      repeatNames[[row]] <- lookup$repeatName[indices]
      repeatClasses[[row]] <- lookup$repeatClass[indices]
      summary <- summarizeClasses(lookup$classes[indices])
      dominantRepeatClass[row] <- summary$dominant
      dominatRepeatClassPercentage[row] <- summary$percentage
      tiedRepeatClasses[[row]] <- summary$tied
    }
  }
  
  o$repeatNames <- repeatNames
  o$repeatClasses <- repeatClasses
  o$dominantRepeatClass <- dominantRepeatClass
  o$dominatRepeatClassPercentage <- dominatRepeatClassPercentage
  o$tiedRepeatClasses <- tiedRepeatClasses
  
  saveRDS(o, file.path(args$outputDir, paste0(args$fileTag, '.rds')))
  if(useDB){
    if(nrow(dbMultiHitGroups)){
      updateLog(paste0('Replacing database multihit clusters for ', ppNum(nrow(dbMultiHitGroups)), ' input group(s).'))
      uploadMultiHitClustersToDB(o, dbMultiHitGroups)
      updateLog('Annotated multihit clusters database upload completed.')
    } else {
      updateLog('Input contains no multihit clusters; no database groups to replace.')
    }
  }
  updateLog(paste0('Completed annotateMultiHitRepeats module. Annotated ', ppNum(nrow(o)),
                   ' clusters; ', ppNum(sum(o$dominantRepeatClass == 'tie')), ' have tied leading classes.'))
  write(date(), file.path(args$outputDir, paste0(args$fileTag, '.done')))
}

#-------------------------------------------------------------------------------

args <- parser$parse_args()
source(file.path(args$softwareRoot, 'lib', 'common.R'))
source(file.path(args$softwareRoot, 'lib', 'buildStdFragments.R'))

tryCatch({
  runModule()
}, error = function(e) {
  cat("ERROR: ", conditionMessage(e), "\n", sep = "", file = stderr())
  flush(stderr())
  quit(save = "no", status = 1, runLast = FALSE)
})

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
})

# -------------------------------------------------------------------------
# Options
# -------------------------------------------------------------------------

optionList <- list(
  make_option(c("-g", "--genome"), type="character", dest="genome",
              help="Genome/assembly name, e.g. hg38, hs1, felCat9 [required]"),
  make_option(c("-o", "--output-dir"), type="character", dest="outputDir", default=NA_character_,
              help="Output directory [default: <genome>_output]"),
  make_option(c("-t", "--threads"), type="integer", dest="threads", default=50L,
              help="Maximum requested CPU threads [default: %default]"),
  
  make_option("--two-bit-source", type="character", dest="twoBitSource", default=NA_character_,
              help="Genome 2bit: local file or URL [default: UCSC <genome>.2bit]"),
  make_option("--annotation-source", type="character", dest="annotationSource", default=NA_character_,
              help="genePred/genePredExt annotation: local file or URL [default: auto-detect UCSC RefSeq]"),
  
  make_option("--chrom-regex", type="character", dest="chromRegex", default="^chr([1-9][0-9]*|X|Y)$",
              help="Regex defining chromosomes to retain [default: chr1...chrN, chrX, chrY]"),
  make_option("--chromosome-list", type="character", dest="chromosomeList", default=NA_character_,
              help="Text file containing chromosome names to retain, one per line; overrides --chrom-regex"),
  make_option("--all-chromosomes", action="store_true", dest="allChromosomes", default=FALSE,
              help="Retain all chromosomes, chrM, scaffolds, alt/fix/random sequences, etc."),
  
  make_option("--repeat-masker-species", type="character", dest="repeatMaskerSpecies", default=NA_character_,
              help='RepeatMasker species, e.g. "Homo sapiens"'),
  make_option("--repeat-masker-path", type="character", dest="repeatMaskerPath", default="RepeatMasker",
              help="RepeatMasker executable [default: %default]"),
  make_option("--repeat-masker-sensitive", action="store_true", dest="repeatMaskerSensitive", default=FALSE,
              help="Use slow/sensitive RepeatMasker mode (-s); standard sensitivity is default"),
  make_option("--repeat-masker-jobs", type="integer", dest="repeatMaskerJobs", default=NA_integer_,
              help="Concurrent chromosome-level RepeatMasker jobs [default: floor(threads/4)]"),
  make_option("--repeat-masker-cores-per-job", type="integer", dest="repeatMaskerCoresPerJob", default=4L,
              help="Cores assumed for each RMBlast RepeatMasker -pa 1 job [default: %default]"),
  make_option("--skip-repeat-masker", action="store_true", dest="skipRepeatMasker", default=FALSE,
              help="Skip RepeatMasker"),
  
  make_option("--two-bit-info-path", type="character", dest="twoBitInfoPath", default="twoBitInfo",
              help="UCSC twoBitInfo executable [default: %default]"),
  make_option("--two-bit-to-fa-path", type="character", dest="twoBitToFaPath", default="twoBitToFa",
              help="UCSC twoBitToFa executable [default: %default]"),
  make_option("--fa-to-two-bit-path", type="character", dest="faToTwoBitPath", default="faToTwoBit",
              help="UCSC faToTwoBit executable [default: %default]"),
  
  make_option("--keep-work-dir", action="store_true", dest="keepWorkDir", default=FALSE,
              help="Keep temporary build files and individual RepeatMasker logs"),
  make_option("--force", action="store_true", dest="force", default=FALSE,
              help="Overwrite an existing build")
)

parser <- OptionParser(
  usage="%prog --genome <assembly> --repeat-masker-species <species> [options]",
  description="Build INSPIIRED2 reference genome objects from a 2bit genome and genePred annotation.",
  option_list=optionList
)

opts <- parse_args(parser)

if(is.null(opts$genome) || !nzchar(opts$genome)){
  print_help(parser)
  stop("--genome is required", call.=FALSE)
}

if(opts$threads < 1L) stop("--threads must be >= 1", call.=FALSE)
if(opts$repeatMaskerCoresPerJob < 1L) stop("--repeat-masker-cores-per-job must be >= 1", call.=FALSE)
if(!is.na(opts$repeatMaskerJobs) && opts$repeatMaskerJobs < 1L) stop("--repeat-masker-jobs must be >= 1", call.=FALSE)
if(!is.na(opts$chromosomeList) && opts$allChromosomes)
  stop("--chromosome-list and --all-chromosomes cannot be used together", call.=FALSE)

if(!is.na(opts$chromosomeList)){
  opts$chromosomeList <- path.expand(opts$chromosomeList)
  if(!file.exists(opts$chromosomeList))
    stop("Chromosome list file does not exist: ", opts$chromosomeList, call.=FALSE)
}

if(is.na(opts$outputDir)) opts$outputDir <- paste0(opts$genome, "_output")

ucscBase <- paste0("https://hgdownload.soe.ucsc.edu/goldenPath/", opts$genome)
if(is.na(opts$twoBitSource)) opts$twoBitSource <- paste0(ucscBase, "/bigZips/", opts$genome, ".2bit")

if(!opts$skipRepeatMasker && (is.na(opts$repeatMaskerSpecies) || !nzchar(opts$repeatMaskerSpecies)))
  stop("--repeat-masker-species is required unless --skip-repeat-masker is used", call.=FALSE)

# -------------------------------------------------------------------------
# Directories / output
# -------------------------------------------------------------------------

dir.create(opts$outputDir, recursive=TRUE, showWarnings=FALSE)
opts$outputDir <- normalizePath(opts$outputDir)

workDir <- file.path(opts$outputDir, paste0(".", opts$genome, "_build"))
logFile <- file.path(opts$outputDir, paste0(opts$genome, "_buildGenomeObjects.log"))

twoBitFile <- file.path(opts$outputDir, paste0(opts$genome, ".2bit"))
tuFile <- file.path(opts$outputDir, paste0(opts$genome, ".TUs.rds"))
exonFile <- file.path(opts$outputDir, paste0(opts$genome, ".exons.rds"))
repeatFile <- file.path(opts$outputDir, paste0(opts$genome, ".repeatTable.gz"))
repeatInfoFile <- file.path(opts$outputDir, paste0(opts$genome, ".repeatMaskerInfo.txt"))

knownOutputs <- c(twoBitFile, tuFile, exonFile, repeatFile, repeatInfoFile, logFile)

if((any(file.exists(knownOutputs)) || dir.exists(workDir)) && !opts$force)
  stop("Existing build files found. Use --force to overwrite.", call.=FALSE)

if(opts$force){
  unlink(knownOutputs)
  unlink(workDir, recursive=TRUE)
}

dir.create(workDir, recursive=TRUE)
cat("", file=logFile)

logMessage <- function(...){
  x <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", paste0(...))
  message(x)
  cat(x, "\n", file=logFile, append=TRUE)
}

# -------------------------------------------------------------------------
# General helpers
# -------------------------------------------------------------------------

resolveExecutable <- function(x){
  if(file.exists(x)) return(normalizePath(x))
  p <- Sys.which(x)
  if(nzchar(p)) return(unname(p))
  stop("Executable not found: ", x, call.=FALSE)
}

commandStatus <- function(x){
  s <- attr(x, "status")
  if(is.null(s)) 0L else as.integer(s)
}

captureCommand <- function(command, args=character()){
  x <- suppressWarnings(system2(command, args=args, stdout=TRUE, stderr=TRUE))
  list(status=commandStatus(x), output=x)
}

runCommand <- function(command, args=character()){
  logMessage("Running: ", basename(command), if(length(args)) paste0(" ", paste(args, collapse=" ")) else "")
  status <- system2(command, args=args)
  if(status != 0L) stop("Command failed with exit status ", status, ": ", command, call.=FALSE)
  invisible(status)
}

isRemoteSource <- function(x) grepl("^(https?|ftp)://", x, ignore.case=TRUE)

tryDownload <- function(url, destination, quiet=TRUE){
  unlink(destination)
  err <- ""
  
  status <- tryCatch(
    suppressWarnings(download.file(url, destination, mode="wb", quiet=quiet)),
    error=function(e){ err <<- conditionMessage(e); NA_integer_ }
  )
  
  ok <- !is.na(status) && status == 0L && file.exists(destination) &&
    !is.na(file.info(destination)$size) && file.info(destination)$size > 0
  
  if(!ok) unlink(destination)
  
  list(ok=ok, status=status, error=err, path=if(ok) destination else NA_character_)
}

getSource <- function(source, label, fileext=""){
  if(isRemoteSource(source)){
    destination <- file.path(workDir, paste0(label, fileext))
    logMessage("Downloading: ", source)
    
    x <- tryDownload(source, destination, quiet=FALSE)
    
    if(!x$ok)
      stop("Unable to download required source: ", source,
           if(nzchar(x$error)) paste0("\n", x$error) else "", call.=FALSE)
    
    return(destination)
  }
  
  if(grepl("^file://", source, ignore.case=TRUE))
    source <- sub("^file://", "", source, ignore.case=TRUE)
  
  source <- path.expand(source)
  
  if(!file.exists(source)) stop("Local source does not exist: ", source, call.=FALSE)
  source <- normalizePath(source)
  if(is.na(file.info(source)$size) || file.info(source)$size == 0)
    stop("Local source is empty: ", source, call.=FALSE)
  
  logMessage("Using local source: ", source)
  source
}

resolveAnnotationSource <- function(){
  if(!is.na(opts$annotationSource)){
    path <- getSource(opts$annotationSource, "geneAnnotation")
    return(list(path=path, source=opts$annotationSource, label="user-supplied annotation"))
  }
  
  candidates <- c(
    "UCSC ncbiRefSeqCurated"=paste0(ucscBase, "/database/ncbiRefSeqCurated.txt.gz"),
    "UCSC ncbiRefSeq genePred"=paste0(ucscBase, "/bigZips/genes/", opts$genome, ".ncbiRefSeq.gp.gz")
  )
  
  destination <- file.path(workDir, "geneAnnotation.gz")
  
  for(i in seq_along(candidates)){
    logMessage("Checking annotation source: ", candidates[i])
    x <- tryDownload(candidates[i], destination)
    
    if(x$ok){
      if(i > 1L)
        logMessage("Preferred ncbiRefSeqCurated annotation unavailable; using UCSC ncbiRefSeq genePred fallback.")
      
      logMessage("Annotation source selected: ", candidates[i])
      return(list(path=destination, source=unname(candidates[i]), label=names(candidates)[i]))
    }
    
    logMessage("Annotation source not available: ", candidates[i])
  }
  
  stop(
    "No default UCSC RefSeq annotation was found for genome '", opts$genome, "'.\n",
    "Supply a local file or URL with --annotation-source.", call.=FALSE
  )
}

isGzip <- function(path){
  con <- base::file(path, "rb")
  on.exit(close(con))
  x <- readBin(con, what="raw", n=2L)
  length(x) == 2L && identical(as.integer(x), c(31L, 139L))
}

writeGzTable <- function(d, path){
  con <- gzfile(path, "wt")
  on.exit(close(con))
  write.table(d, con, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

publishFile <- function(source, destination){
  if(file.exists(destination)) unlink(destination)
  
  if(!file.rename(source, destination)){
    if(!file.copy(source, destination, overwrite=TRUE))
      stop("Unable to publish output file: ", destination, call.=FALSE)
    unlink(source)
  }
  
  invisible(destination)
}

readChromosomeList <- function(path){
  path <- normalizePath(path)
  x <- trimws(readLines(path, warn=FALSE))
  
  if(length(x)) x[1L] <- sub("^\ufeff", "", x[1L])
  
  x <- x[nzchar(x)]
  x <- x[!grepl("^#", x)]
  x <- unique(x)
  
  if(!length(x))
    stop("Chromosome list contains no chromosome names: ", path, call.=FALSE)
  
  logMessage("Chromosome list: ", path)
  x
}

# -------------------------------------------------------------------------
# 2bit helpers
# -------------------------------------------------------------------------

readTwoBitInfo <- function(twoBit, output){
  runCommand(twoBitInfoPath, c(shQuote(twoBit), shQuote(output)))
  d <- read.delim(output, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  
  if(ncol(d) != 2L) stop("Unexpected twoBitInfo output.", call.=FALSE)
  
  names(d) <- c("chrom", "length")
  d$length <- as.numeric(d$length)
  d
}

# -------------------------------------------------------------------------
# genePred / genePredExt
# -------------------------------------------------------------------------

genePredColumns <- c(
  "bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
  "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat",
  "cdsEndStat", "exonFrames"
)

readGenePred <- function(path){
  con <- if(isGzip(path)) gzfile(path, "rt") else base::file(path, "rt")
  on.exit(close(con))
  
  d <- read.delim(con, header=FALSE, sep="\t", quote="", comment.char="",
                  stringsAsFactors=FALSE, check.names=FALSE)
  
  if(ncol(d) == 16L){
    names(d) <- genePredColumns
    format <- "UCSC bin + genePredExt (16 columns)"
    
  } else if(ncol(d) == 15L){
    names(d) <- genePredColumns[-1L]
    d$bin <- NA_integer_
    d <- d[, genePredColumns]
    format <- "genePredExt (15 columns)"
    
  } else if(ncol(d) == 10L){
    names(d) <- c("name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
                  "exonCount", "exonStarts", "exonEnds")
    
    d$bin <- NA_integer_
    d$score <- 0L
    d$name2 <- NA_character_
    d$cdsStartStat <- NA_character_
    d$cdsEndStat <- NA_character_
    d$exonFrames <- NA_character_
    d <- d[, genePredColumns]
    format <- "genePred (10 columns)"
    
  } else {
    stop("Expected 10-column genePred, 15-column genePredExt, or 16-column UCSC bin+genePredExt; found ",
         ncol(d), " columns.", call.=FALSE)
  }
  
  intCols <- c("bin", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "score")
  d[intCols] <- lapply(d[intCols], as.integer)
  
  if(any(!d$strand %in% c("+", "-")))
    stop("Annotation contains invalid strand values.", call.=FALSE)
  
  if(any(d$txStart < 0L | d$txEnd <= d$txStart))
    stop("Annotation contains invalid transcript coordinates.", call.=FALSE)
  
  if(any(d$exonCount < 1L))
    stop("Annotation contains records with exonCount < 1.", call.=FALSE)
  
  list(data=d, format=format)
}

validateAnnotation <- function(d, chromInfo){
  sourceRows <- nrow(d)
  sourceGenes <- length(unique(d$name2[!is.na(d$name2) & nzchar(d$name2)]))
  
  d <- d[d$chrom %in% chromInfo$chrom, , drop=FALSE]
  
  if(!nrow(d))
    stop("No annotation records correspond to retained reference chromosomes.", call.=FALSE)
  
  chrLength <- chromInfo$length[match(d$chrom, chromInfo$chrom)]
  
  badTx <- which(d$txStart < 0L | d$txEnd > chrLength | d$txStart >= d$txEnd)
  if(length(badTx))
    stop("Annotation contains ", length(badTx),
         " transcript(s) outside reference chromosome bounds.", call.=FALSE)
  
  badCds <- which(d$cdsStart < d$txStart | d$cdsEnd > d$txEnd | d$cdsStart > d$cdsEnd)
  if(length(badCds))
    stop("Annotation contains ", length(badCds),
         " invalid CDS coordinate record(s).", call.=FALSE)
  
  starts <- strsplit(sub(",+$", "", d$exonStarts), ",", fixed=TRUE)
  ends <- strsplit(sub(",+$", "", d$exonEnds), ",", fixed=TRUE)
  
  if(any(lengths(starts) != lengths(ends)))
    stop("exonStarts/exonEnds length mismatch.", call.=FALSE)
  
  if(any(lengths(starts) != d$exonCount))
    stop("exonCount disagrees with exonStarts/exonEnds.", call.=FALSE)
  
  badExon <- vapply(seq_len(nrow(d)), function(i){
    s <- suppressWarnings(as.integer(starts[[i]]))
    e <- suppressWarnings(as.integer(ends[[i]]))
    
    anyNA(s) || anyNA(e) || any(s < d$txStart[i]) || any(e > d$txEnd[i]) ||
      any(s < 0L) || any(e > chrLength[i]) || any(s >= e)
  }, logical(1))
  
  if(any(badExon))
    stop("Annotation contains ", sum(badExon),
         " transcript(s) with invalid exon coordinates.", call.=FALSE)
  
  retainedGenes <- length(unique(d$name2[!is.na(d$name2) & nzchar(d$name2)]))
  
  list(data=d, sourceRows=sourceRows, retainedRows=nrow(d),
       sourceGenes=sourceGenes, retainedGenes=retainedGenes)
}

applySeqinfo <- function(g, chromInfo, genome){
  seqs <- GenomeInfoDb::seqlevels(g)
  i <- match(seqs, chromInfo$chrom)
  
  if(anyNA(i))
    stop("GRanges contains chromosomes absent from reference genome.", call.=FALSE)
  
  GenomeInfoDb::seqinfo(g) <- GenomeInfoDb::Seqinfo(
    seqnames=seqs, seqlengths=chromInfo$length[i], genome=genome
  )
  
  g
}

createTranscriptGRanges <- function(d, chromInfo, genome){
  g <- makeGRangesFromDataFrame(
    d, seqnames.field="chrom", start.field="txStart", end.field="txEnd",
    strand.field="strand", starts.in.df.are.0based=TRUE, keep.extra.columns=TRUE
  )
  
  applySeqinfo(g, chromInfo, genome)
}

createExonGRanges <- function(d, chromInfo, genome){
  starts <- strsplit(sub(",+$", "", d$exonStarts), ",", fixed=TRUE)
  ends <- strsplit(sub(",+$", "", d$exonEnds), ",", fixed=TRUE)
  nExons <- lengths(starts)
  rowIndex <- rep.int(seq_len(nrow(d)), nExons)
  
  e <- d[rowIndex, , drop=FALSE]
  e$name <- paste("exon", sequence(nExons))
  e$exonStarts <- as.integer(unlist(starts, use.names=FALSE))
  e$exonEnds <- as.integer(unlist(ends, use.names=FALSE))
  
  g <- makeGRangesFromDataFrame(
    e, seqnames.field="chrom", start.field="exonStarts", end.field="exonEnds",
    strand.field="strand", starts.in.df.are.0based=TRUE, keep.extra.columns=TRUE
  )
  
  applySeqinfo(g, chromInfo, genome)
}

# -------------------------------------------------------------------------
# RepeatMasker
# -------------------------------------------------------------------------

repeatMaskerNames <- c(
  "SW_score", "percent_div", "percent_del", "percent_ins", "query_seq",
  "query_start", "query_end", "query_after", "strand", "repeat_name",
  "repeat_class", "repeat_start", "repeat_end", "repeat_after", "ID", "alt"
)

emptyRepeatMasker <- function(){
  data.frame(
    SW_score=integer(), percent_div=numeric(), percent_del=numeric(), percent_ins=numeric(),
    query_seq=character(), query_start=integer(), query_end=integer(), query_after=character(),
    strand=character(), repeat_name=character(), repeat_class=character(),
    repeat_start=character(), repeat_end=character(), repeat_after=character(),
    ID=character(), alt=character(), stringsAsFactors=FALSE
  )
}

readRepeatMasker <- function(path){
  if(!file.exists(path))
    stop("RepeatMasker output not found: ", path, call.=FALSE)
  
  d <- tryCatch(
    read.table(path, skip=3L, fill=TRUE, header=FALSE, quote="", comment.char="",
               stringsAsFactors=FALSE),
    error=function(e) NULL
  )
  
  if(is.null(d) || !nrow(d)) return(emptyRepeatMasker())
  
  d[[1L]] <- suppressWarnings(as.integer(d[[1L]]))
  d <- d[!is.na(d[[1L]]), , drop=FALSE]
  
  if(!nrow(d)) return(emptyRepeatMasker())
  
  if(ncol(d) == 15L){
    names(d) <- repeatMaskerNames[-16L]
    d$alt <- NA_character_
  } else if(ncol(d) == 16L){
    names(d) <- repeatMaskerNames
  } else {
    stop("Unexpected RepeatMasker output: ", ncol(d),
         " columns in ", path, call.=FALSE)
  }
  
  d[, repeatMaskerNames]
}

repeatMaskerSensitivityArgs <- function(){
  if(opts$repeatMaskerSensitive) "-s" else character(0)
}

getRepeatMaskerVersion <- function(){
  x <- captureCommand(repeatMaskerPath, "-v")
  hit <- grep("RepeatMasker.*version", x$output, value=TRUE, ignore.case=TRUE)
  
  if(!length(hit)){
    x <- captureCommand(repeatMaskerPath, "-h")
    hit <- grep("RepeatMasker.*version", x$output, value=TRUE, ignore.case=TRUE)
  }
  
  if(length(hit)) trimws(hit[1L]) else "RepeatMasker version unavailable"
}

repeatMaskerPreflight <- function(){
  preflightDir <- file.path(workDir, "repeatMaskerPreflight")
  unlink(preflightDir, recursive=TRUE)
  dir.create(preflightDir)
  
  fasta <- file.path(preflightDir, "test.fa")
  seq <- paste(rep("ACGTGCGCATAT", 500L), collapse="")
  writeLines(c(">RepeatMasker_preflight", seq), fasta)
  
  args <- c(
    repeatMaskerSensitivityArgs(), "-pa", "1", "-e", "rmblast",
    "-species", shQuote(opts$repeatMaskerSpecies),
    "-dir", shQuote(preflightDir), shQuote(fasta)
  )
  
  logMessage("Running RepeatMasker preflight.")
  result <- captureCommand(repeatMaskerPath, args)
  
  writeLines(
    c(paste0("Command: ", repeatMaskerPath, " ", paste(args, collapse=" ")),
      "", result$output),
    file.path(preflightDir, "RepeatMasker.preflight.log")
  )
  
  outFile <- file.path(preflightDir, "test.fa.out")
  
  if(result$status != 0L || !file.exists(outFile)){
    logMessage("RepeatMasker preflight FAILED.")
    logMessage("Preflight directory retained: ", preflightDir)
    
    x <- tail(result$output, 30L)
    if(length(x)) message(paste(x, collapse="\n"))
    
    stop("RepeatMasker preflight failed. See RepeatMasker.preflight.log.", call.=FALSE)
  }
  
  info <- grep(
    "RepeatMasker version|Search Engine:|Using FamDB:|Families\\s*:|families in ancestor|lineage-specific",
    result$output, value=TRUE
  )
  
  if(length(info)) invisible(lapply(unique(trimws(info)), logMessage))
  
  logMessage("RepeatMasker preflight passed.")
  unlink(preflightDir, recursive=TRUE)
  
  result$output
}

runRepeatMaskerJob <- function(chrom, index, total){
  safeChrom <- gsub("[^A-Za-z0-9_.-]", "_", chrom)
  jobDir <- file.path(repeatWorkDir, sprintf("%05d_%s", index, safeChrom))
  fasta <- file.path(jobDir, paste0(safeChrom, ".fa"))
  jobLog <- file.path(jobDir, "RepeatMasker.log")
  
  dir.create(jobDir, recursive=TRUE)
  
  message(sprintf("RepeatMasker [%d/%d] starting %s", index, total, chrom))
  startTime <- Sys.time()
  
  extractArgs <- c(shQuote(stagedTwoBit), shQuote(fasta), shQuote(paste0("-seq=", chrom)))
  extract <- captureCommand(twoBitToFaPath, extractArgs)
  
  cat(
    paste0("Sequence: ", chrom, "\n"),
    paste0("Start: ", format(startTime), "\n"),
    paste0("twoBitToFa command: ", twoBitToFaPath, " ", paste(extractArgs, collapse=" "), "\n\n"),
    paste(extract$output, collapse="\n"), "\n\n", file=jobLog
  )
  
  if(extract$status != 0L){
    elapsed <- as.numeric(difftime(Sys.time(), startTime, units="secs"))
    return(list(chrom=chrom, status=extract$status, outFile=NA_character_,
                logFile=jobLog, seconds=elapsed, error="twoBitToFa failed"))
  }
  
  rmArgs <- c(
    repeatMaskerSensitivityArgs(), "-pa", "1", "-e", "rmblast",
    "-species", shQuote(opts$repeatMaskerSpecies),
    "-dir", shQuote(jobDir), shQuote(fasta)
  )
  
  rm <- captureCommand(repeatMaskerPath, rmArgs)
  
  cat(
    paste0("RepeatMasker command: ", repeatMaskerPath, " ", paste(rmArgs, collapse=" "), "\n\n"),
    paste(rm$output, collapse="\n"), "\n", file=jobLog, append=TRUE
  )
  
  outFile <- file.path(jobDir, paste0(basename(fasta), ".out"))
  status <- rm$status
  
  if(status == 0L && !file.exists(outFile)){
    candidates <- list.files(jobDir, pattern="\\.out$", full.names=TRUE)
    if(length(candidates) == 1L) outFile <- candidates else status <- 1L
  }
  
  elapsed <- as.numeric(difftime(Sys.time(), startTime, units="secs"))
  
  message(sprintf(
    "RepeatMasker [%d/%d] finished %s (%s, %.1f min)",
    index, total, chrom, if(status == 0L) "OK" else "FAILED", elapsed / 60
  ))
  
  list(
    chrom=chrom, status=status,
    outFile=if(status == 0L) outFile else NA_character_,
    logFile=jobLog, seconds=elapsed,
    error=if(status == 0L) NA_character_ else "RepeatMasker failed"
  )
}

# -------------------------------------------------------------------------
# Resolve executables
# -------------------------------------------------------------------------

twoBitInfoPath <- resolveExecutable(opts$twoBitInfoPath)
twoBitToFaPath <- resolveExecutable(opts$twoBitToFaPath)

if(!opts$allChromosomes) faToTwoBitPath <- resolveExecutable(opts$faToTwoBitPath)
if(!opts$skipRepeatMasker) repeatMaskerPath <- resolveExecutable(opts$repeatMaskerPath)

# -------------------------------------------------------------------------
# INPUT PREFLIGHT
# All required inputs are acquired and validated before RepeatMasker starts.
# -------------------------------------------------------------------------

logMessage("Genome: ", opts$genome)
logMessage("2bit source: ", opts$twoBitSource)

# Genome

sourceTwoBit <- getSource(opts$twoBitSource, "sourceGenome", ".2bit")
sourceChromInfoFile <- file.path(workDir, "source.chromInfo.txt")
sourceChromInfo <- readTwoBitInfo(sourceTwoBit, sourceChromInfoFile)

logMessage("Sequences in source genome: ", format(nrow(sourceChromInfo), big.mark=","))

# Chromosome selection

chromInfo <- sourceChromInfo

if(!is.na(opts$chromosomeList)){
  allowedChromosomes <- readChromosomeList(opts$chromosomeList)
  missingChromosomes <- setdiff(allowedChromosomes, sourceChromInfo$chrom)
  
  if(length(missingChromosomes))
    stop(
      "Chromosome list contains sequence(s) not present in the reference genome:\n",
      paste(missingChromosomes, collapse=", "),
      call.=FALSE
    )
  
  chromInfo <- sourceChromInfo[sourceChromInfo$chrom %in% allowedChromosomes, , drop=FALSE]
  
  logMessage("Chromosome selection: explicit chromosome list")
  logMessage("Chromosomes requested: ", length(allowedChromosomes))
  
} else if(opts$allChromosomes){
  logMessage("Chromosome selection: all sequences (--all-chromosomes)")
  
} else {
  chromInfo <- sourceChromInfo[
    grepl(opts$chromRegex, sourceChromInfo$chrom, perl=TRUE),
    ,
    drop=FALSE
  ]
  
  logMessage("Chromosome selection: regex")
  logMessage("Chromosome filter: ", opts$chromRegex)
}

if(!nrow(chromInfo))
  stop("No chromosomes remain after chromosome filtering.", call.=FALSE)

logMessage("Sequences retained: ", format(nrow(chromInfo), big.mark=","))
logMessage("Sequences excluded: ", format(nrow(sourceChromInfo) - nrow(chromInfo), big.mark=","))
logMessage("Reference size: ", format(sum(chromInfo$length), big.mark=","), " bp")

if(nrow(chromInfo) <= 100L)
  logMessage("Sequences: ", paste(chromInfo$chrom, collapse=", "))

# Annotation

annotation <- resolveAnnotationSource()

logMessage("Reading annotation.")
annotationData <- readGenePred(annotation$path)
d <- annotationData$data

logMessage("Annotation source: ", annotation$source)
logMessage("Annotation source type: ", annotation$label)
logMessage("Annotation format: ", annotationData$format)

validated <- validateAnnotation(d, chromInfo)
d <- validated$data

logMessage("Source transcript records: ", format(validated$sourceRows, big.mark=","))
logMessage("Retained transcript records: ", format(validated$retainedRows, big.mark=","))

if(validated$sourceGenes > 0L)
  logMessage("Source gene symbols: ", format(validated$sourceGenes, big.mark=","))

if(validated$retainedGenes > 0L)
  logMessage("Retained gene symbols: ", format(validated$retainedGenes, big.mark=","))

if(grepl("genePred \\(10 columns\\)", annotationData$format))
  logMessage("WARNING: 10-column genePred has no gene-symbol/name2 metadata.")

logMessage("Genome and annotation preflight passed.")

# RepeatMasker installation/database preflight

rmVersion <- NA_character_
rmPreflightOutput <- character()

if(!opts$skipRepeatMasker){
  rmVersion <- getRepeatMaskerVersion()
  logMessage("RepeatMasker: ", rmVersion)
  rmPreflightOutput <- repeatMaskerPreflight()
}

# -------------------------------------------------------------------------
# Build staged reference 2bit
# -------------------------------------------------------------------------

stagedTwoBit <- file.path(workDir, paste0(opts$genome, ".2bit"))

if(opts$allChromosomes){
  logMessage("Copying complete source 2bit.")
  
  if(!file.copy(sourceTwoBit, stagedTwoBit, overwrite=TRUE))
    stop("Unable to copy source 2bit.", call.=FALSE)
  
} else {
  chromListFile <- file.path(workDir, "chromosomes.txt")
  genomeFasta <- file.path(workDir, paste0(opts$genome, ".fa"))
  
  writeLines(chromInfo$chrom, chromListFile)
  
  logMessage("Extracting retained chromosomes.")
  runCommand(
    twoBitToFaPath,
    c(shQuote(sourceTwoBit), shQuote(genomeFasta), shQuote(paste0("-seqList=", chromListFile)))
  )
  
  logMessage("Building filtered 2bit.")
  runCommand(faToTwoBitPath, c(shQuote(genomeFasta), shQuote(stagedTwoBit)))
  
  unlink(genomeFasta)
}

verifyInfoFile <- file.path(workDir, "final.chromInfo.txt")
verifyInfo <- readTwoBitInfo(stagedTwoBit, verifyInfoFile)

if(!identical(verifyInfo$chrom, chromInfo$chrom) ||
   length(verifyInfo$length) != length(chromInfo$length) ||
   any(verifyInfo$length != chromInfo$length))
  stop("Final 2bit does not match expected retained chromosomes.", call.=FALSE)

logMessage("Final 2bit verified.")

# -------------------------------------------------------------------------
# RepeatMasker
# -------------------------------------------------------------------------

if(!opts$skipRepeatMasker){
  repeatStart <- Sys.time()
  repeatWorkDir <- file.path(workDir, "repeatMasker")
  dir.create(repeatWorkDir)
  
  if(is.na(opts$repeatMaskerJobs))
    repeatMaskerJobs <- max(1L, floor(opts$threads / opts$repeatMaskerCoresPerJob))
  else
    repeatMaskerJobs <- opts$repeatMaskerJobs
  
  repeatMaskerJobs <- min(repeatMaskerJobs, nrow(chromInfo))
  estimatedCores <- repeatMaskerJobs * opts$repeatMaskerCoresPerJob
  sensitivityLabel <- if(opts$repeatMaskerSensitive) "slow/sensitive (-s)" else "standard"
  
  logMessage("RepeatMasker species: ", opts$repeatMaskerSpecies)
  logMessage("RepeatMasker sensitivity: ", sensitivityLabel)
  logMessage("RepeatMasker engine: rmblast")
  logMessage("Requested threads: ", opts$threads)
  logMessage("Assumed cores per RepeatMasker -pa 1 job: ", opts$repeatMaskerCoresPerJob)
  logMessage("Concurrent RepeatMasker jobs: ", repeatMaskerJobs)
  logMessage("Estimated peak RMBlast cores: ", estimatedCores)
  
  if(estimatedCores > opts$threads)
    logMessage("WARNING: RepeatMasker configuration may exceed requested thread count.")
  
  rmChromInfo <- chromInfo[order(chromInfo$length, decreasing=TRUE), , drop=FALSE]
  
  logMessage(
    "Starting ", nrow(rmChromInfo),
    " chromosome/scaffold-level RepeatMasker jobs, largest sequences first."
  )
  
  rmResults <- parallel::mclapply(
    seq_len(nrow(rmChromInfo)),
    function(i){
      tryCatch(
        runRepeatMaskerJob(rmChromInfo$chrom[i], i, nrow(rmChromInfo)),
        error=function(e){
          list(
            chrom=rmChromInfo$chrom[i], status=1L, outFile=NA_character_,
            logFile=NA_character_, seconds=NA_real_, error=conditionMessage(e)
          )
        }
      )
    },
    mc.cores=repeatMaskerJobs,
    mc.preschedule=FALSE
  )
  
  failed <- vapply(
    rmResults,
    function(x) is.null(x$status) || x$status != 0L,
    logical(1)
  )
  
  if(any(failed)){
    failedChroms <- vapply(rmResults[failed], `[[`, character(1), "chrom")
    logMessage("RepeatMasker failures: ", paste(failedChroms, collapse=", "))
    logMessage("Work directory retained: ", repeatWorkDir)
    stop("RepeatMasker failed for ", length(failedChroms), " sequence(s).", call.=FALSE)
  }
  
  resultChroms <- vapply(rmResults, `[[`, character(1), "chrom")
  rmResults <- rmResults[match(chromInfo$chrom, resultChroms)]
  repeatFiles <- vapply(rmResults, `[[`, character(1), "outFile")
  
  logMessage("Reading RepeatMasker outputs.")
  
  repeatTables <- lapply(repeatFiles, readRepeatMasker)
  repeats <- do.call(rbind, repeatTables)
  
  logMessage("RepeatMasker records: ", format(nrow(repeats), big.mark=","))
  
  stagedRepeatFile <- file.path(workDir, paste0(opts$genome, ".repeatTable.gz"))
  writeGzTable(repeats, stagedRepeatFile)
  
  repeatEnd <- Sys.time()
  
  commandTemplate <- paste0(
    repeatMaskerPath, " ",
    if(opts$repeatMaskerSensitive) "-s " else "",
    '-pa 1 -e rmblast -species "', opts$repeatMaskerSpecies,
    '" -dir <jobDir> <sequence.fa>'
  )
  
  rmEnvironment <- unique(trimws(grep(
    "RepeatMasker version|Search Engine:|Using FamDB:|Families\\s*:|families in ancestor|lineage-specific",
    rmPreflightOutput, value=TRUE
  )))
  
  chromosomeMode <- if(!is.na(opts$chromosomeList)){
    paste0("chromosome list: ", normalizePath(opts$chromosomeList))
  } else if(opts$allChromosomes){
    "all reference sequences"
  } else {
    paste0("regex: ", opts$chromRegex)
  }
  
  provenance <- c(
    paste0("Genome: ", opts$genome),
    paste0("Build date: ", format(Sys.time())),
    paste0("Chromosome selection: ", chromosomeMode),
    paste0("Sequences processed: ", nrow(chromInfo)),
    paste0("Reference size: ", sum(chromInfo$length), " bp"),
    paste0("RepeatMasker executable: ", repeatMaskerPath),
    paste0("RepeatMasker version: ", rmVersion),
    paste0("Search engine: rmblast"),
    paste0("Sensitivity: ", sensitivityLabel),
    paste0("Species: ", opts$repeatMaskerSpecies),
    paste0("Requested threads: ", opts$threads),
    paste0("Assumed RMBlast cores per -pa 1 job: ", opts$repeatMaskerCoresPerJob),
    paste0("Concurrent RepeatMasker jobs: ", repeatMaskerJobs),
    paste0("Estimated peak RMBlast cores: ", estimatedCores),
    paste0("Repeat records: ", nrow(repeats)),
    paste0("Elapsed hours: ", round(as.numeric(difftime(repeatEnd, repeatStart, units="hours")), 3)),
    paste0("Command template: ", commandTemplate)
  )
  
  if(length(rmEnvironment))
    provenance <- c(provenance, "", "RepeatMasker environment:", rmEnvironment)
  
  stagedRepeatInfoFile <- file.path(workDir, paste0(opts$genome, ".repeatMaskerInfo.txt"))
  writeLines(provenance, stagedRepeatInfoFile)
  
  rm(repeats, repeatTables)
  gc()
}

# -------------------------------------------------------------------------
# GRanges
# -------------------------------------------------------------------------

logMessage("Building transcript GRanges.")

transcripts <- createTranscriptGRanges(d, chromInfo, opts$genome)
stagedTuFile <- file.path(workDir, paste0(opts$genome, ".TUs.rds"))
saveRDS(transcripts, stagedTuFile)

logMessage("Transcript GRanges: ", format(length(transcripts), big.mark=","), " ranges")

logMessage("Building exon GRanges.")

exons <- createExonGRanges(d, chromInfo, opts$genome)
stagedExonFile <- file.path(workDir, paste0(opts$genome, ".exons.rds"))
saveRDS(exons, stagedExonFile)

logMessage("Exon GRanges: ", format(length(exons), big.mark=","), " ranges")

# -------------------------------------------------------------------------
# Publish completed build
# -------------------------------------------------------------------------

logMessage("Publishing completed genome build.")

publishFile(stagedTwoBit, twoBitFile)
publishFile(stagedTuFile, tuFile)
publishFile(stagedExonFile, exonFile)

if(!opts$skipRepeatMasker){
  publishFile(stagedRepeatFile, repeatFile)
  publishFile(stagedRepeatInfoFile, repeatInfoFile)
}

if(!opts$keepWorkDir) unlink(workDir, recursive=TRUE)

logMessage("Done.")
logMessage("Output directory: ", opts$outputDir)
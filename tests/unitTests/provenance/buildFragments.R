#!/usr/bin/env Rscript
# Execute the actual buildFragments runModule against paired-alignment fixtures.
# Run: Rscript tests/unitTests/provenance/buildFragments.R [softwareRoot] [baselineRoot]
# No aligner, FASTQs or database are used; only startup/logging are substituted.
required <- c('data.table', 'dplyr', 'yaml')
for(package in required){
  if(!requireNamespace(package, quietly = TRUE)) stop('Install test dependency: ', package)
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}
argv <- commandArgs(trailingOnly = TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1L])
softwareRoot <- normalizePath(if(length(argv)) argv[[1L]] else
  file.path(dirname(script), '..', '..', '..'), mustWork = TRUE)
baselineRoot <- if(length(argv) > 1L) normalizePath(argv[[2L]], mustWork = TRUE) else NULL
source(file.path(softwareRoot, 'lib', 'common.R'))

run_tests <- function(){
  private <- tempfile('provenance_buildFragments_')
  dir.create(private)
  on.exit(unlink(private, recursive = TRUE), add = TRUE)
  anchors <- data.table(
    readID = c('A_U5', 'B_U3', 'C_cross_chromosome', 'D_same_strand',
               'E_too_short', 'F_too_long', 'G_multihit', 'G_multihit', 'H_unpaired'),
    tName = c('chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr2', 'chr3', 'chr1'),
    strand = c('+', '-', '+', '+', '+', '+', '+', '+', '+'),
    tStart = c(100L, 300L, 100L, 100L, 400L, 500L, 1000L, 2000L, 100L),
    tEnd = c(129L, 329L, 129L, 129L, 409L, 529L, 1029L, 2029L, 129L),
    seq = 'ACGTACGTACGT', trial = factor('LOD'), subject = factor('CLONE'),
    sample = factor(c('DILUTION_1', 'DILUTION_1', rep('DILUTION_2', 7L))),
    replicate = factor(c(1L, 5L, rep(1L, 7L))),
    refGenome = factor('hs1'), mode = factor(c('U5', 'U3', rep('U5', 7L))),
    nReads = c(8L, 4L, 1L, 1L, 1L, 1L, 2L, 2L, 1L),
    UMI = 'AAAAAAAAAAAA', leaderSeq = 'GAAAATCTCTAGCA',
    leaderSeqHMM = c('U5.hmm', 'U3.hmm', rep('U5.hmm', 7L)),
    vectorFastaFile = 'vector.fasta')
  adrift <- data.table(
    readID = c('A_U5', 'B_U3', 'C_cross_chromosome', 'D_same_strand',
               'E_too_short', 'F_too_long', 'G_multihit', 'G_multihit'),
    seq = 'TGCATGCATGCA',
    tName = c('chr1', 'chr1', 'chr2', 'chr1', 'chr1', 'chr1', 'chr2', 'chr3'),
    strand = c('-', '+', '-', '+', '-', '-', '-', '-'),
    tStart = c(170L, 250L, 170L, 170L, 410L, 200000L, 1070L, 2070L),
    tEnd = c(199L, 279L, 199L, 199L, 420L, 200029L, 1099L, 2099L))
  input <- list(anchorReads = anchors, adriftReads = adrift)
  ancestry <- lapply(c('U5', 'U3'), function(mode){
    arguments <- list(softwareRoot = softwareRoot, version = 'INSPIIRED2 test',
      HMMparams = paste0(mode, ',1,5,', if(mode == 'U5') '55,75' else '40,65', ',TRUE,CA,2'))
    demux <- provenanceRecord('demultiplex', argsSnapshot = arguments)
    prep <- provenanceRecord('prepReads', parents = list(demux), argsSnapshot = arguments,
      details = list(hmmParameters = data.frame(mode = mode,
        HMMminFullBitScore = if(mode == 'U5') 55 else 40)))
    provenanceRecord('alignReads', parents = list(prep), argsSnapshot = arguments)
  })
  inputWithHistory <- provenanceAttach(input, provenanceMerge(ancestry))
  extract_module <- function(path){
    for(expression in parse(path, keep.source = FALSE)){
      if(is.call(expression) && identical(expression[[1L]], as.name('<-')) &&
         identical(expression[[2L]], as.name('runModule'))) return(expression[[3L]])
    }
    stop('runModule not found: ', path)
  }
  execute <- function(root, object, label){
    out <- file.path(private, label)
    dir.create(out)
    inputFile <- file.path(out, 'alignReads.rds')
    saveRDS(object, inputFile)
    env <- new.env(parent = globalenv())
    env$args <- list(outputDir = out, inputData = inputFile,
      softwareRoot = root, threads = 1L, fileTag = 'buildFragments',
      ramDiskPath = out, dataRowChunkSize = 2L,
      minFrgamentLength = 40L, maxFrgamentLength = 100000L,
      dbConfigFile = 'none', dbConfigID = 'none', overwriteDBrecords = FALSE)
    env$startModule <- function(){
      for(name in c('tmpDir', 'logDir', 'ramDisk')){
        env$args[[name]] <- file.path(out, name)
        dir.create(env$args[[name]])
      }
      env$args$version <- 'INSPIIRED2 test'
      env$args$startTime <- Sys.time()
      env$args$logFile <- file.path(out, 'buildFragments.log')
    }
    env$updateLog <- function(...) invisible(NULL)
    env$runModule <- eval(extract_module(file.path(root, 'modules', 'buildFragments.R')), env)
    env$runModule()
    stopifnot(file.exists(file.path(out, 'buildFragments.done')))
    readRDS(file.path(out, 'buildFragments.rds'))
  }
  result <- execute(softwareRoot, inputWithHistory, 'current')
  stopifnot(nrow(result) == 4L,
    identical(as.character(result$readID), c('A_U5', 'B_U3', 'G_multihit', 'G_multihit')),
    identical(result$fragStart, c(100L, 250L, 1000L, 2000L)),
    identical(result$fragEnd, c(199L, 329L, 1099L, 2099L)),
    identical(as.integer(result$fragWidth), c(100L, 80L, 100L, 100L)),
    identical(result$fragStrand, c('+', '-', '+', '+')),
    identical(result$fragChromosome, c('chr1', 'chr1', 'chr2', 'chr3')),
    identical(as.character(result$mode), c('U5', 'U3', 'U5', 'U5')),
    identical(result$nReads, c(8L, 4L, 2L, 2L)),
    identical(result$anchor_seq, rep('ACGTACGTACGT', 4L)),
    identical(result$adrift_seq, rep('TGCATGCATGCA', 4L)))
  cat('PASS: actual buildFragments pairing, filtering, strand, coordinates, counts and multihits\n')
  p <- provenanceFromObject(result)
  executions <- Filter(function(record) identical(record$kind, 'execution'), p$records)
  modules <- vapply(executions, `[[`, '', 'module')
  stopifnot(identical(p$records[[p$heads]]$module, 'buildFragments'),
    sum(modules == 'prepReads') == 2L, sum(modules == 'alignReads') == 2L,
    all(vapply(ancestry, function(parent) all(names(parent$records) %in% names(p$records)), logical(1))))
  arguments <- yaml::yaml.load(p$records[[p$heads]]$arguments_yaml, eval.expr = FALSE)
  stopifnot(arguments$minFrgamentLength == 40L, arguments$maxFrgamentLength == 100000L,
    arguments$dataRowChunkSize == 2L)
  cat('PASS: U5/U3 histories and actual buildFragments settings survive processing\n')
  if(!is.null(baselineRoot)){
    baseline <- execute(baselineRoot, input, 'baseline')
    stopifnot(identical(provenanceStrip(result), baseline))
    cat('PASS: stripped current buildFragments output is identical to unmodified baseline\n')
  } else {
    cat('Baseline-code comparison not requested; explicit scientific expectations were checked.\n')
  }
}
run_tests()

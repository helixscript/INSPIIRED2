#!/usr/bin/env Rscript
# Run real annotation module bodies against small GRanges/repeat fixtures.
# Optional second argument: a v1.0.8 source directory for before/after equality.
argv <- commandArgs(trailingOnly = TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1L])
softwareRoot <- normalizePath(if(length(argv)) argv[1L] else
  file.path(dirname(script), '..', '..', '..'), mustWork = TRUE)
baseRoot <- if(length(argv) > 1L) normalizePath(argv[2L], mustWork = TRUE) else NULL
required <- c('data.table', 'dplyr', 'tibble', 'stringr', 'readr',
              'GenomicRanges', 'GenomeInfoDb', 'IRanges', 'yaml', 'DBI', 'arrow')
for(p in required){
  if(!requireNamespace(p, quietly = TRUE)) stop('Install test dependency: ', p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
source(file.path(softwareRoot, 'lib', 'common.R'))

run_tests <- function(){
  previousOptions <- options(useFancyQuotes = FALSE)
  on.exit(options(previousOptions), add = TRUE)
  testRoot <- tempfile('provenance_downstream_')
  dir.create(testRoot)
  on.exit(unlink(testRoot, recursive = TRUE), add = TRUE)
  fixtures <- file.path(testRoot, 'software')
  annotations <- file.path(fixtures, 'data', 'genomeAnnotations')
  dir.create(annotations, recursive = TRUE)
  file.copy(file.path(softwareRoot, 'VERSION'), file.path(fixtures, 'VERSION'))
  # The execution fingerprint identifies the actual source under test.
  for(part in c('modules', 'lib')){
    dir.create(file.path(fixtures, part))
    paths <- list.files(file.path(softwareRoot, part), full.names = TRUE)
    file.copy(paths[file.info(paths)$isdir == FALSE], file.path(fixtures, part))
  }
  file.copy(file.path(softwareRoot, 'inspiired2.R'), fixtures)

  tuA <- GRanges('chr1', IRanges(c(100L, 400L), c(200L, 450L)), strand = c('+', '-'))
  tuA$name2 <- c('GENE_A', 'GENE_B')
  exA <- GRanges('chr1', IRanges(c(110L, 180L, 410L), c(120L, 190L, 420L)),
                 strand = c('+', '+', '-'))
  exA$name2 <- c('GENE_A', 'GENE_A', 'GENE_B')
  tuB <- GRanges('chr1', IRanges(600L, 700L), strand = '+')
  tuB$name2 <- 'GENE_C'
  exB <- tuB
  for(item in c('tuA', 'exA', 'tuB', 'exB')){
    genome <- if(endsWith(item, 'A')) 'gA' else 'gB'
    suffix <- if(startsWith(item, 'tu')) '.TUs.rds' else '.exons.rds'
    saveRDS(get(item), file.path(annotations, paste0(genome, suffix)))
  }
  repeatA <- data.frame(query_seq = 'chr1', query_start = c(105L, 105L, 400L),
    query_end = c(115L, 115L, 430L), strand = c('+', 'C', '+'),
    repeat_name = c('AluY', 'AT_rich', 'L1HS'),
    repeat_class = c('SINE/Alu', 'Simple_repeat', 'LINE/L1'))
  repeatB <- data.frame(query_seq = 'chr1', query_start = 490L, query_end = 510L,
    strand = '+', repeat_name = 'ERVK', repeat_class = 'LTR/ERVK')
  for(genome in c('gA', 'gB')){
    con <- gzfile(file.path(annotations, paste0(genome, '.repeatTable.gz')), 'wt')
    write.table(if(genome == 'gA') repeatA else repeatB, con,
                sep = '\t', row.names = FALSE, quote = FALSE)
    close(con)
  }

  provenanceArgs <- list(softwareRoot = fixtures, version = 'INSPIIRED2 test',
                         HMMparams = 'U5_PARAMETERS')
  u5 <- provenanceRecord('prepReads', argsSnapshot = provenanceArgs,
                         details = list(mode = 'U5', effectiveScore = 55))
  provenanceArgs$HMMparams <- 'U3_PARAMETERS'
  u3 <- provenanceRecord('prepReads', argsSnapshot = provenanceArgs,
                         details = list(mode = 'U3', effectiveScore = 40))
  parent <- provenanceRecord('buildSites', parents = list(u5, u3),
                             argsSnapshot = provenanceArgs)
  sites <- data.table(trial = 'TRIAL', subject = 'SUBJECT',
    sample = c('DILUTION_1', 'DILUTION_1', 'DILUTION_2', 'DILUTION_3'),
    refGenome = c('gA', 'gA', 'gA', 'gB'),
    mode = c('dual detect', 'U5', 'U3', 'U5'),
    posid = c('chr1+110', 'chr1-250', 'chr1+425', 'chr1+500'),
    sonicLengths = c(8L, 3L, 1L, 1L),
    pullSubjectFragments = factor(rep('TRUE', 4L), levels = c('FALSE', 'TRUE')))
  rebuild <- list(version = 2L, pullSubjectFragments = TRUE, siteRows = 4L,
    groups = unique(as.data.frame(sites)[c('trial', 'subject', 'sample', 'refGenome')]))
  attr(sites, 'inspiired2.siteRebuild') <- rebuild
  sites <- provenanceAttach(sites, parent)

  extract_module <- function(path){
    for(expr in parse(path, keep.source = FALSE)){
      if(is.call(expr) && identical(expr[[1L]], as.name('<-')) &&
         identical(expr[[2L]], as.name('runModule'))) return(expr[[3L]])
    }
    stop('runModule not found: ', path)
  }
  execute <- function(module, input, label, root = softwareRoot){
    out <- file.path(testRoot, label)
    dir.create(out)
    inputFile <- file.path(out, 'input.rds')
    saveRDS(input, inputFile)
    env <- new.env(parent = globalenv())
    env$args <- list(outputDir = out, inputData = inputFile, fileTag = module,
      softwareRoot = fixtures, threads = 1L, ramDiskPath = out,
      dbConfigFile = 'none', dbConfigID = 'none')
    env$startModule <- function(...){
      for(nm in c('tmpDir', 'logDir', 'ramDisk')){
        env$args[[nm]] <- file.path(out, nm)
        dir.create(env$args[[nm]])
      }
      env$args$logFile <- file.path(out, paste0(module, '.log'))
      env$args$version <- 'INSPIIRED2 test'
      env$args$startTime <- Sys.time()
    }
    env$resource_overlay <- function(...) invisible(TRUE)
    env$updateLog <- function(...) invisible(NULL)
    env$runModule <- eval(extract_module(file.path(root, 'modules', paste0(module, '.R'))), env)
    env$runModule()
    stopifnot(file.exists(file.path(out, paste0(module, '.done'))))
    readRDS(file.path(out, paste0(module, '.rds')))
  }
  check_history <- function(x, module, resourceCount){
    p <- provenanceFromObject(x)
    execution <- p$records[[p$heads]]
    stopifnot(identical(execution$module, module),
      length(execution$resources) == resourceCount,
      all(names(u5$records) %in% names(p$records)),
      all(names(u3$records) %in% names(p$records)))
    for(resource in execution$resources)
      stopifnot(identical(resource$md5, unname(tools::md5sum(resource$path))))
    invisible(TRUE)
  }

  genes <- execute('nearestGenes', sites, 'genes')
  stopifnot(identical(genes$nearestGene, c('GENE_A', 'GENE_A', 'GENE_B', 'GENE_C')),
    identical(genes$nearestGeneDist, c(0L, 50L, 0L, 100L)),
    identical(genes$inGene, c(TRUE, FALSE, TRUE, FALSE)),
    identical(genes$inExon, c(TRUE, FALSE, FALSE, FALSE)),
    identical(attr(genes, 'inspiired2.siteRebuild'), rebuild),
    identical(genes$sonicLengths, sites$sonicLengths))
  check_history(genes, 'nearestGenes', 4L)

  repeats <- execute('annotateRepeats', genes, 'repeats')
  stopifnot(identical(repeats$repeat_name, c('AluY,AT_rich', NA_character_, 'L1HS', 'ERVK')),
    identical(repeats$repeat_class, c('SINE/Alu,Simple_repeat', NA_character_, 'LINE/L1', 'LTR/ERVK')),
    identical(attr(repeats, 'inspiired2.siteRebuild'), rebuild),
    identical(repeats$sonicLengths, sites$sonicLengths))
  check_history(repeats, 'annotateRepeats', 2L)

  clusters <- data.table(trial = 'TRIAL', subject = 'SUBJECT', sample = 'DILUTION_1',
    mode = 'U5', refGenome = c('gA', 'gA', 'gB'), clusterID = c('MHC.1', 'MHC.2', 'MHC.3'),
    posids = list(c('chr1+1.1e+02', 'chr1-250', 'chr1+425'),
                  c('chr1+105', 'chr1+110'), c('chr1+500', 'chr1+700')))
  clusters$nodeSonicLengths <- lapply(clusters$posids, function(p)
    data.table(posid = p, sonicLengths = rep(1L, length(p))))
  clusters <- provenanceAttach(clusters, parent)
  annotated <- execute('annotateMultiHitRepeats', clusters, 'multihit')
  stopifnot(identical(annotated$posids[[1L]][1L], 'chr1+110'),
    identical(annotated$nodeSonicLengths[[1L]]$posid[1L], 'chr1+110'),
    identical(annotated$repeatNames[[1L]], c('AluY,AT_rich', 'none', 'L1HS')),
    identical(annotated$dominantRepeatClass, rep('tie', 3L)),
    identical(annotated$dominatRepeatClassPercentage, c(1/3, 1, 0.5)),
    identical(annotated$tiedRepeatClasses[[2L]], c('SINE/Alu', 'Simple_repeat')))
  check_history(annotated, 'annotateMultiHitRepeats', 2L)
  empty <- execute('annotateMultiHitRepeats', clusters[0L], 'empty_multihit')
  stopifnot(nrow(empty) == 0L)
  check_history(empty, 'annotateMultiHitRepeats', 0L)

  legacy <- provenanceStrip(sites)
  legacyGenes <- execute('nearestGenes', legacy, 'legacy_genes')
  legacyRecords <- provenanceFromObject(legacyGenes)$records
  stopifnot(any(vapply(legacyRecords, function(r)
    identical(r$kind, 'source') && identical(r$history, 'unavailable'), logical(1))))

  # Real Parquet retrieval through the complete pullDBrecords body. Only the
  # database connection/query are fixtures; no MariaDB server is contacted.
  lake <- file.path(testRoot, 'lake')
  dir.create(lake)
  sourceFragments <- data.table(trial = 'TRIAL', subject = 'SUBJECT',
    sample = factor('DILUTION_1'), replicate = c(1L, 5L), refGenome = 'gA',
    mode = factor(c('U5', 'U3')), fragChromosome = 'chr1',
    fragStrand = c('+', '-'), fragStart = c(100L, 102L),
    fragEnd = c(200L, 202L), readID = c('FLOWCELL_U5_READ1', 'FLOWCELL_U3_READ1'))
  archiveNames <- vapply(seq_len(2L), function(i){
    temporary <- file.path(lake, paste0('stage_', i, '.parquet'))
    provenanceWriteParquet(sourceFragments[i], temporary,
      provenance = if(i == 1L) u5 else u3)
    name <- paste0(unname(tools::md5sum(temporary)), '.parquet')
    stopifnot(file.rename(temporary, file.path(lake, name)))
    name
  }, character(1))
  records <- as.data.frame(sourceFragments)[c('trial', 'subject', 'sample',
                                             'replicate', 'refGenome', 'mode')]
  names(records)[names(records) == 'refGenome'] <- 'ref_genome'
  records$data_file_name <- archiveNames
  config <- file.path(testRoot, 'fixture.cnf')
  writeLines('# Fixture only; never opened as a database connection.', config)
  methods::setClass('INSPIIRED2ProvenanceFixtureConnection', contains = 'DBIConnection')
  methods::setMethod('dbDisconnect', 'INSPIIRED2ProvenanceFixtureConnection',
                     function(conn, ...) invisible(TRUE))
  execute_pull <- function(root, label){
    env <- new.env(parent = globalenv())
    env$args <- list(dbConfigFile = config, dbConfigID = 'fixture', dataPath = lake,
      trials = 'TRIAL', subjects = 'none', samples = 'none', refGenomes = 'none',
      modes = 'none', softwareRoot = fixtures, outputFile = file.path(testRoot, label))
    env$createDBconnection <- function(...)
      methods::new('INSPIIRED2ProvenanceFixtureConnection')
    env$dbGetQuery <- function(conn, query){
      stopifnot(grepl('from fragments where', query, fixed = TRUE))
      if(startsWith(query, 'select data_file_name')) records['data_file_name'] else records
    }
    env$runModule <- eval(extract_module(file.path(root, 'modules', 'pullDBrecords.R')), env)
    env$runModule()
    readRDS(env$args$outputFile)
  }
  pulled <- execute_pull(softwareRoot, 'pulled.rds')
  pullHistory <- provenanceFromObject(pulled)
  stopifnot(nrow(pulled) == 2L,
    identical(pullHistory$records[[pullHistory$heads]]$module, 'pullDBrecords'),
    all(names(u5$records) %in% names(pullHistory$records)),
    all(names(u3$records) %in% names(pullHistory$records)))
  archivedSources <- Filter(function(r) identical(r$kind, 'source') &&
    identical(r$source$table, 'fragments'), pullHistory$records)
  stopifnot(length(archivedSources) == 2L,
    identical(unname(sort(vapply(archivedSources, function(r)
      r$source$database_key$mode, character(1)))), c('U3', 'U5')))

  if(!is.null(baseRoot)){
    oldGenes <- execute('nearestGenes', sites, 'old_genes', baseRoot)
    oldRepeats <- execute('annotateRepeats', oldGenes, 'old_repeats', baseRoot)
    oldClusters <- execute('annotateMultiHitRepeats', clusters, 'old_multihit', baseRoot)
    oldEmpty <- execute('annotateMultiHitRepeats', clusters[0L], 'old_empty', baseRoot)
    oldPulled <- execute_pull(baseRoot, 'old_pulled.rds')
    for(pair in list(list(genes, oldGenes), list(repeats, oldRepeats),
                     list(annotated, oldClusters), list(empty, oldEmpty),
                     list(pulled, oldPulled))){
      result <- all.equal(provenanceStrip(pair[[1L]]), provenanceStrip(pair[[2L]]))
      if(!isTRUE(result)) stop(paste(result, collapse = '\n'))
    }
    cat('PASS: all five outputs equal v1.0.8 after removing only processing provenance.\n')
  }
  cat('PASS: actual nearestGenes, annotateRepeats, and annotateMultiHitRepeats module bodies; two-genome resources; U5/U3 history; legacy and empty MHC inputs.\n')
  cat('PASS: actual pullDBrecords body with real Parquet files, per-key inherited histories, and a simulated read-only database.\n')
}
run_tests()

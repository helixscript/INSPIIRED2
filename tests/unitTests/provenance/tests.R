#!/usr/bin/env Rscript
# Run: Rscript tests/unitTests/provenance/tests.R [softwareRoot] [baselineRoot]
# Uses real R/Arrow/filesystem operations. Only database access is substituted.
# An optional unmodified v1.0.8 checkout enables an additional exact comparison
# of the buildSites scientific result. No sequencing data or live DB is needed.

required <- c("yaml", "arrow", "data.table", "dplyr", "tibble", "stringr", "lubridate")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if(length(missing)) stop("Missing test packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tibble)
  library(stringr); library(lubridate)
})
argv <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(if(length(argv)) argv[[1L]] else
  file.path(dirname(script), "..", "..", ".."), mustWork = TRUE)
baseline <- if(length(argv) > 1L) normalizePath(argv[[2L]], mustWork = TRUE) else NULL

main <- function(){
  private <- tempfile("inspiired2_provenance_tests_")
  dir.create(private)
  on.exit(unlink(private, recursive = TRUE), add = TRUE)
  e <- new.env(parent = globalenv())
  sys.source(file.path(root, "lib", "common.R"), envir = e)
  sys.source(file.path(root, "lib", "buildStdFragments.R"), envir = e)
  e$args <- list(softwareRoot = root, version = "INSPIIRED2 provenance tests",
                 outputDir = private, fileTag = "fixture", HMMparams = "U5,1,5,55,75,TRUE,CA,2")
  e$updateLog <- function(...) invisible(NULL)
  n <- 0L
  expect <- function(ok, message = "assertion failed") if(!isTRUE(ok)) stop(message, call. = FALSE)
  expect_error <- function(expr, pattern = NULL){
    error <- tryCatch({ force(expr); NULL }, error = function(err) conditionMessage(err))
    expect(!is.null(error), "Expected an error, but the operation succeeded")
    if(!is.null(pattern)) expect(grepl(pattern, error, ignore.case = TRUE), paste("Unexpected error:", error))
    invisible(error)
  }
  test <- function(name, code){
    force(code)
    n <<- n + 1L
    cat(sprintf("PASS %02d: %s\n", n, name))
  }
  record <- function(module, parents = list(), details = list(), args = e$args,
                     resources = character(), inputs = character())
    e$provenanceRecord(module, parents, details, resources, inputs, argsSnapshot = args)
  modules <- function(p) vapply(Filter(function(r) identical(r$kind, "execution"), p$records), `[[`, "", "module")
  records_of <- function(p, kind) Filter(function(r) identical(r$kind, kind), p$records)
  text_of <- function(p) e$provenanceToYAML(p)
  key <- c("trial", "subject", "sample", "replicate", "refGenome", "mode")
  fragment <- function(mode = "U5", replicate = 1L, sample = "DILUTION_1"){
    data.frame(trial = "LOD", subject = "CLONE", sample = sample,
      replicate = replicate, refGenome = "hs1", mode = mode, readID = paste(mode, sample),
      leaderSeqHMM = paste0(mode, ".hmm"), vectorFastaFile = "vector.fasta",
      fragChromosome = "chr1", fragStrand = if(mode == "U5") "-" else "+",
      fragStart = if(mode == "U5") 4L else 100L,
      fragEnd = if(mode == "U5") 104L else 200L,
      UMI = "AAAAAAAAAAAA", nReads = 4L, stringsAsFactors = FALSE)
  }
  make_history <- function(mode){
    a <- e$args; a$HMMparams <- paste0(mode, ",1,5,", if(mode == "U5") "55,75" else "40,65", ",TRUE,CA,2")
    demux <- record("demultiplex", details = list(mode = mode), args = a)
    prep <- record("prepReads", list(demux), details = list(hmmParameters =
      data.frame(leaderSeqHMM = paste0(mode, ".hmm"), HMMminFullBitScore = if(mode == "U5") 55 else 40)), args = a)
    aligned <- record("alignReads", list(prep), args = a)
    record("buildFragments", list(aligned), args = a)
  }

  test("capture preserves effective settings, resources, and excludes credentials", {
    model <- file.path(private, "U5.hmm"); writeLines("test model", model)
    input <- file.path(private, "reads.fastq"); writeLines("test input", input)
    a <- e$args; a$dbConfigFile <- "DO_NOT_STORE.cnf"; a$dbConfigID <- "PRIVATE_ADMIN"
    a$dbConn <- new.env(); a$password <- "DO_NOT_STORE_PASSWORD"; a$command <- "--password DO_NOT_STORE_COMMAND"
    p <- record("prepReads", details = list(effective = data.frame(score = 55, end = "CA")),
                args = a, resources = model, inputs = input)
    r <- p$records[[p$heads[[1L]]]]
    expect(grepl("HMMparams", r$arguments_yaml))
    expect(!grepl("DO_NOT_STORE|PRIVATE_ADMIN", text_of(p)))
    expect(identical(r$resources[[1L]]$md5, unname(tools::md5sum(model))))
    expect(identical(r$input_files[[1L]]$identity, "path_size_mtime_only"))
    expect(grepl("55", text_of(p)) && grepl("CA", text_of(p)))
    expect(length(r$software$code) >= 2L, "Code identity was not captured")
  })

  test("execution capture neither creates nor advances R's RNG state", {
    if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    invisible(record("buildSites"))
    expect(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    set.seed(579); before <- get(".Random.seed", .GlobalEnv)
    invisible(record("buildSites")); expect(identical(before, get(".Random.seed", .GlobalEnv)))
  })

  u5 <- make_history("U5"); u3 <- make_history("U3")
  test("distinct executions and shared ancestors merge without duplication", {
    merged <- e$provenanceMerge(list(u5, u5, u3))
    expect(length(merged$records) == length(u5$records) + length(u3$records))
    expect(setequal(merged$heads, c(u5$heads, u3$heads)))
    repeated <- record("buildFragments", list(u5))
    expect(!identical(repeated$heads, u5$heads))
    expect(sum(modules(repeated) == "buildFragments") == 2L)
    conflicted <- u5; conflicted$records[[u5$heads[[1L]]]]$module <- "differentModule"
    expect_error(e$provenanceMerge(list(u5, conflicted)), "conflict")
  })

  test("version, missing parents, cycles, and malformed YAML are rejected", {
    p <- u5; p$version <- 99L; expect_error(e$provenanceValidate(p), "unsupported|malformed")
    p <- u5; p$records[[p$heads[[1L]]]]$parents <- "missing"; expect_error(e$provenanceValidate(p), "missing")
    p <- u5; id <- p$heads[[1L]]; p$records[[id]]$parents <- id; expect_error(e$provenanceValidate(p), "cycle")
    expect_error(e$provenanceFromYAML("not: [valid"), "decode")
    expect_error(e$provenanceFromYAML("version: 99"), "unsupported|malformed")
  })

  test("portable YAML roundtrip preserves immutable record content", {
    p <- record("buildSites", list(u5, u3), details = list(unicode = "αβ", flags = c(TRUE, FALSE), empty = list()))
    expect(identical(text_of(p), text_of(e$provenanceFromYAML(text_of(p)))))
  })

  test("legacy input is explicitly unavailable, never assigned current HMM settings", {
    p <- e$provenanceFromObject(fragment(), source = list(path = "legacy.rds"))
    expect(length(records_of(p, "execution")) == 0L)
    expect(all(vapply(records_of(p, "source"), function(r) identical(r$history, "unavailable"), logical(1))))
    expect(!grepl("HMMparams", text_of(p)))
    x <- fragment(); attr(x, "inspiired2.provenance") <- list(version = 999)
    expect_error(e$provenanceFromObject(x), "malformed|unsupported")
  })

  test("attach/strip preserve data, rebuild authorization, and NULL list members", {
    d <- as.data.table(fragment()); d$mode <- factor(d$mode, levels = c("U5", "U3"))
    attr(d, "inspiired2.siteRebuild") <- list(version = 2L, test = TRUE)
    attr(d, "unrelated") <- "keep"
    nested <- list(nothing = NULL, data = d, another = NULL)
    with_metadata <- e$provenanceAttach(nested, u5)
    with_metadata$data <- e$provenanceAttach(with_metadata$data, u5)
    stripped <- e$provenanceStrip(with_metadata)
    expect(identical(names(stripped), names(nested)) && length(stripped) == 3L)
    expect(identical(stripped, nested))
  })

  test("strip/reserialize reproduces original compressed RDS bytes", {
    d <- as.data.table(fragment()); d$mode <- factor(d$mode)
    d$pullSubjectFragments <- factor(TRUE, levels = c(FALSE, TRUE))
    attr(d, "inspiired2.siteRebuild") <- list(version = 2L, pullSubjectFragments = TRUE,
      siteRows = nrow(d), groups = as.data.frame(d)[c("trial", "subject", "sample", "refGenome")])
    paths <- file.path(private, c("baseline.rds", "with_history.rds", "stripped.rds"))
    saveRDS(d, paths[[1L]])
    saveRDS(e$provenanceAttach(d, u5), paths[[2L]])
    saveRDS(e$provenanceStrip(readRDS(paths[[2L]])), paths[[3L]])
    hashes <- unname(tools::md5sum(paths))
    expect(identical(hashes[[1L]], hashes[[3L]]), "Scientific serialization changed after stripping")
    expect(hashes[[1L]] != hashes[[2L]])
  })

  archive <- function(d, p, label){
    path <- file.path(private, paste0(label, ".parquet"))
    e$provenanceWriteParquet(e$provenanceAttach(d, p), path)
    path
  }
  u5_file <- archive(fragment(), u5, "u5")
  u3_file <- archive(fragment("U3", 5L), u3, "u3")

  test("real Parquet preserves values, factors, rebuild attributes and metadata", {
    d <- fragment(); d$mode <- factor(d$mode, levels = c("U5", "U3"))
    attr(d, "inspiired2.siteRebuild") <- list(version = 2L, pullSubjectFragments = TRUE)
    path <- archive(d, u5, "roundtrip")
    metadata <- arrow::ParquetFileReader$create(path)$GetSchema()$metadata[["inspiired2.provenance"]]
    expect(is.character(metadata) && length(metadata) == 1L)
    expect(identical(e$provenanceStrip(e$provenanceReadParquet(path)), d))
    decoded <- e$provenanceFromYAML(metadata)
    expect(identical(decoded$heads, u5$heads))
    expect(decoded$artifact$row_count == nrow(d))
    expect(all(key %in% unlist(decoded$artifact$scope$columns)))
  })

  test("archive imports retain exact file/checksum/key and both HMM histories", {
    a <- e$provenanceReadParquet(u5_file, source = list(database_key = as.list(fragment()[1L, key])))
    b <- e$provenanceReadParquet(u3_file, source = list(database_key = as.list(fragment("U3", 5L)[1L, key])))
    pa <- e$provenanceFromObject(a); pb <- e$provenanceFromObject(b)
    imported <- pa$records[[pa$heads[[1L]]]]
    expect(identical(imported$source$file$md5, unname(tools::md5sum(u5_file))))
    expect(identical(imported$source$database_key$mode, "U5"))
    expect(identical(imported$parents, u5$heads))
    rebuilt <- record("buildStdFragments", list(pa, pb), details = list(pullSubjectFragments = TRUE))
    expect(sum(modules(rebuilt) == "prepReads") == 2L)
    expect(all(names(u5$records) %in% names(rebuilt$records)))
    expect(all(names(u3$records) %in% names(rebuilt$records)))
    expect(grepl("55", text_of(rebuilt)) && grepl("40", text_of(rebuilt)))
  })

  test("legacy Parquet remains readable and malformed present metadata fails", {
    path <- file.path(private, "legacy.parquet"); arrow::write_parquet(fragment(), path)
    legacy <- e$provenanceFromObject(e$provenanceReadParquet(path))
    expect(length(records_of(legacy, "execution")) == 0L)
    expect(identical(legacy$records[[legacy$heads[[1L]]]]$history, "unavailable"))
    table <- arrow::Table$create(fragment()); table$metadata[["inspiired2.provenance"]] <- "version: 999"
    bad <- file.path(private, "malformed.parquet"); arrow::write_parquet(table, bad)
    expect_error(e$provenanceReadParquet(bad), "unsupported|malformed")
  })

  test("embedded artifact scope must match the actual Parquet rows and keys", {
    original <- e$provenanceScope(fragment(), u5)
    wrong <- rbind(fragment(), fragment(sample = "OTHER_SAMPLE"))
    table <- arrow::Table$create(wrong)
    table$metadata[["inspiired2.provenance"]] <- text_of(original)
    path <- file.path(private, "wrong_scope.parquet"); arrow::write_parquet(table, path)
    expect_error(e$provenanceReadParquet(path), "scope|rows|keys")
  })

  test("derived archive remains self-contained after source archives are deleted", {
    a <- e$provenanceFromObject(e$provenanceReadParquet(u5_file))
    b <- e$provenanceFromObject(e$provenanceReadParquet(u3_file))
    std <- record("buildStdFragments", list(a, b))
    sites <- record("buildSites", list(std))
    genes <- record("nearestGenes", list(sites))
    repeats <- record("annotateRepeats", list(genes))
    path <- archive(fragment(), repeats, "final_sites")
    unlink(c(u5_file, u3_file))
    restored <- e$provenanceFromObject(e$provenanceReadParquet(path))
    expect(all(names(repeats$records) %in% names(restored$records)))
    expect(sum(modules(restored) == "prepReads") == 2L)
    expect(all(c("buildSites", "nearestGenes", "annotateRepeats") %in% modules(restored)))
  })

  test("partition metadata scopes each file while retaining all computational parents", {
    p <- record("buildSites", list(u5, u3))
    all <- rbind(fragment(sample = "DILUTION_1"), fragment(sample = "DILUTION_10"))
    attached <- e$provenanceAttach(all, p)
    # Explicit parent bundle is passed after a transformation that drops attrs.
    partition <- as.data.frame(lapply(attached[1L, , drop = FALSE], identity))
    path <- file.path(private, "partition.parquet")
    e$provenanceWriteParquet(partition, path, provenance = e$provenanceFromObject(attached))
    encoded <- arrow::ParquetFileReader$create(path)$GetSchema()$metadata[["inspiired2.provenance"]]
    decoded <- e$provenanceFromYAML(encoded)
    expect(length(decoded$artifact$scope$rows) == 1L)
    expect(identical(decoded$artifact$scope$rows[[1L]]$sample, "DILUTION_1"))
    expect(all(names(p$records) %in% names(decoded$records)))
  })

  test("actual pullDBfragments merges per-key archive histories and skips incoming key", {
    incoming <- e$provenanceAttach(fragment("U3", 5L), u3)
    u5_path <- archive(fragment(), u5, "pull_u5")
    more_path <- file.path(private, "pull_legacy.parquet")
    arrow::write_parquet(fragment("U5", 2L), more_path)
    make_row <- function(d, file, present = 0L){
      r <- d[1L, key, drop = FALSE]; r$data_file_name <- basename(file); r$in_input <- present; r
    }
    rows <- rbind(make_row(fragment(), u5_path), make_row(fragment("U5", 2L), more_path),
                  make_row(fragment("U3", 5L), "does_not_exist.parquet", 1L))
    mock <- new.env(parent = e)
    cfg <- file.path(private, "test.cnf"); file.create(cfg)
    mock$args <- c(e$args, list(dbConfigFile = cfg, dbConfigID = "fixture", dataPath = private))
    mock$createDBconnection <- function() list(test = TRUE)
    mock$mock <- list(dbGetQuery = function(conn, query, params){
      expect(length(params) == 8L, "Pull query lost exact input-key/subject parameters")
      rows
    }, dbDisconnect = function(conn) TRUE)
    rewrite <- function(node){
      if(missing(node)) return(quote(expr = ))
      if(is.call(node)){
        if(identical(node[[1L]], as.name("::")) && identical(node[[2L]], as.name("DBI")))
          return(call("$", as.name("mock"), node[[3L]]))
        return(as.call(lapply(as.list(node), rewrite)))
      }
      if(is.pairlist(node)) return(as.pairlist(lapply(node, rewrite)))
      node
    }
    pull <- e$pullDBfragments; body(pull) <- rewrite(body(pull)); environment(pull) <- mock
    out <- pull(incoming)
    expect(nrow(out) == 2L && all(out$mode == "U5"))
    p <- e$provenanceFromObject(out)
    expect(all(names(u5$records) %in% names(p$records)))
    expect(!any(names(u3$records) %in% names(p$records)), "Skipped incoming data was falsely re-imported")
    sources <- records_of(p, "source")
    expect(any(vapply(sources, function(r) identical(r$history, "unavailable"), logical(1))))
    expect(grepl("pull_u5.parquet", text_of(p)) && grepl("pull_legacy.parquet", text_of(p)))
    rows <- rows[FALSE, , drop = FALSE]
    empty <- pull(incoming)
    expect(nrow(empty) == 0L)
    expect(length(e$provenanceFromObject(empty)$records) == 0L,
           "A zero-row pull inherited history from the incoming batch")
  })

  test("export utility CLI reads real RDS/Parquet provenance without modifying inputs", {
    if(!requireNamespace("argparse", quietly = TRUE)) stop("Install argparse for the export CLI test")
    tool <- file.path(root, "tools", "exportProvenance.R")
    rds <- file.path(private, "input with spaces.rds")
    saveRDS(e$provenanceAttach(fragment(), u5), rds)
    parquet <- archive(fragment(), u5, "export")
    for(path in c(rds, parquet)){
      before <- unname(tools::md5sum(path))
      output <- paste0(path, ".yml")
      command <- c("--vanilla", shQuote(tool), "--inputData", shQuote(path),
                   "--outputFile", shQuote(output), "--softwareRoot", shQuote(root))
      result <- system2(file.path(R.home("bin"), "Rscript"), command, stdout = TRUE, stderr = TRUE)
      expect(is.null(attr(result, "status")), paste(result, collapse = "\n"))
      exported <- e$provenanceFromYAML(paste(readLines(output), collapse = "\n"))
      expect(identical(text_of(exported), text_of(e$provenanceScope(fragment(), u5))))
      expect(identical(unname(tools::md5sum(path)), before))
    }
    result <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      c("--vanilla", shQuote(tool), "--inputData", shQuote(rds), "--outputFile", shQuote(rds)),
      stdout = TRUE, stderr = TRUE))
    expect(!is.null(attr(result, "status")) && attr(result, "status") != 0L)
    expect(grepl("must not replace", paste(result, collapse = "\n")))
  })

  # Actual buildSites execution: no DB, no external sequence clustering needed.
  module_expression <- function(path){
    expressions <- parse(path, keep.source = FALSE)
    matches <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
      identical(x[[2L]], as.name("runModule")), as.list(expressions))
    expect(length(matches) == 1L, "Could not locate runModule")
    matches[[1L]][[3L]]
  }
  run_sites <- function(moduleRoot, frags, disabled = FALSE){
    outdir <- tempfile("sites_", tmpdir = private); dir.create(outdir)
    input <- file.path(outdir, "input.rds"); saveRDS(frags, input)
    run <- new.env(parent = e)
    run$args <- list(inputData = input, outputDir = outdir, softwareRoot = root,
      fileTag = "buildSites", dbConfigFile = "none", dbConfigID = "none",
      disableDualDetect = disabled, disableOrientationCorrection = disabled,
      dualDetectWidth = 6L, integraseCorrectionDist = 2L, sumSonicBreaksWithin = "replicates",
      leadSeqClusteringParms = "unused", threads = 1L, ramDiskPath = outdir)
    run$startModule <- function(...){
      run$args$version <- "INSPIIRED2 test"; run$args$command <- "test fixture"
      run$args$startTime <- Sys.time()
      for(nm in c("tmpDir", "logDir", "ramDisk")){
        run$args[[nm]] <- file.path(outdir, nm); dir.create(run$args[[nm]])
      }
      run$args$logFile <- file.path(outdir, "buildSites.log")
    }
    run$updateLog <- function(...) invisible(NULL)
    run$timeElapsedString <- function() "0h 00min 00sec"
    fn <- eval(module_expression(file.path(moduleRoot, "modules", "buildSites.R")), run)
    invisible(capture.output(suppressMessages(fn())))
    expect(file.exists(file.path(outdir, "buildSites.done")))
    readRDS(file.path(outdir, "buildSites.rds"))
  }
  std_fixture <- function(){
    d <- rbindlist(list(fragment("U5", 1L), fragment("U3", 5L),
                       fragment("U5", 1L, "DILUTION_10"), fragment("U3", 5L, "DILUTION_10")))
    d[, posid := ifelse(mode == "U5", "chr1-104", "chr1+100")]
    d[, reads := c(8L, 4L, 2L, 1L)]
    d[, repLeaderSeq := "ACGTACGT"]
    d[, leaderSeqGroupNum := 1L]
    d$UMIs <- rep(list("AAAAAAAAAAAA"), nrow(d))
    d$pullSubjectFragments <- factor(rep(TRUE, nrow(d)), levels = c(FALSE, TRUE))
    attr(d, "inspiired2.siteRebuild") <- list(version = 1L, pullSubjectFragments = TRUE,
      fragmentRows = nrow(d), groups = unique(as.data.frame(d)[c("trial", "subject", "sample", "refGenome")]))
    e$provenanceAttach(d, record("buildStdFragments", list(u5, u3)))
  }
  test("actual buildSites preserves U5/U3 histories and yields expected dual results", {
    d <- std_fixture(); sites <- run_sites(root, d)
    expect(nrow(sites) == 2L && all(sites$mode == "dual detect"))
    expect(all(sites$posid == "chr1-102") && identical(as.integer(sites$reads), c(12L, 3L)))
    expect(all(sites$sonicLengths == 2L))
    expect(identical(attr(sites, "inspiired2.siteRebuild")$version, 2L))
    p <- e$provenanceFromObject(sites)
    expect(sum(modules(p) == "prepReads") == 2L && tail(unname(modules(p)), 1L) == "buildSites")
    if(!is.null(baseline)){
      old <- run_sites(baseline, e$provenanceStrip(d))
      expect(identical(e$provenanceStrip(sites), old), "Default buildSites result differs from unmodified baseline")
    }
  })
  test("actual alternative buildSites run records its flags without losing ancestry", {
    d <- std_fixture(); sites <- run_sites(root, d, disabled = TRUE)
    expect(nrow(sites) == 4L && setequal(as.character(sites$mode), c("U5", "U3")))
    expect(setequal(sites$posid, c("chr1-104", "chr1+100")))
    p <- e$provenanceFromObject(sites); r <- p$records[[p$heads[[1L]]]]
    a <- yaml::yaml.load(r$arguments_yaml, eval.expr = FALSE)
    expect(isTRUE(a$disableDualDetect) && isTRUE(a$disableOrientationCorrection))
    expect(sum(modules(p) == "prepReads") == 2L)
    if(!is.null(baseline)){
      old <- run_sites(baseline, e$provenanceStrip(d), disabled = TRUE)
      expect(identical(e$provenanceStrip(sites), old), "Alternative buildSites result differs from unmodified baseline")
    }
  })
  cat(sprintf("\nAll %d provenance test groups passed with real Arrow %s.\n", n, as.character(packageVersion("arrow"))))
  if(is.null(baseline)) cat("Baseline-code comparison not requested; fixture biological expectations were checked.\n")
  cat("Database queries were simulated; no live MariaDB or FASTQ pipeline was exercised.\n")
}
main()

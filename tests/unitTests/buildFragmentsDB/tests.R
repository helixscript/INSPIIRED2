#!/usr/bin/env Rscript
# Exercises the real nested uploader using an in-memory transactional DB double.
# Real temporary files, checksums, hard links, and unlink operations are used.
# Arrow serialization is replaced with deterministic R serialization; these are
# fault-injection unit tests, not live MariaDB or Arrow compatibility tests.
# Provenance transport is substituted here; tests/unitTests/provenance/tests.R
# separately exercises the real metadata helpers with actual Arrow Parquet.

if(!requireNamespace("data.table", quietly = TRUE))
  stop("Install data.table to run these tests.", call. = FALSE)
suppressPackageStartupMessages(library(data.table))
argv <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- sub("^--file=", "", script_arg[[1L]])
module <- if(length(argv)) argv[[1L]] else
  file.path(dirname(script), "..", "..", "..", "modules", "buildFragments.R")
module <- normalizePath(module, mustWork = TRUE)

find_function <- function(node, name){
  if(missing(node)) return(NULL)
  if(is.call(node) && as.character(node[[1L]])[1L] %in% c("<-", "=") &&
     is.symbol(node[[2L]]) && identical(as.character(node[[2L]]), name) &&
     is.call(node[[3L]]) && identical(node[[3L]][[1L]], as.name("function")))
    return(node[[3L]])
  if(is.call(node) || is.expression(node) || is.pairlist(node)){
    for(part in as.list(node)){
      found <- find_function(part, name)
      if(!is.null(found)) return(found)
    }
  }
  NULL
}
parsed <- parse(module, keep.source = FALSE)
run_module <- find_function(parsed, "runModule")
uploader_expression <- find_function(run_module, "uploadFragmentsToDB")
if(is.null(uploader_expression))
  stop("Expected uploadFragmentsToDB() nested inside runModule().", call. = FALSE)

key_names <- c("trial", "subject", "sample", "replicate", "ref_genome", "mode")
empty_records <- function(){
  data.frame(trial = character(), subject = character(), sample = character(),
             replicate = integer(), ref_genome = character(), mode = character(),
             total_fragments = integer(), data_file_name = character(),
             stringsAsFactors = FALSE)
}
fixture <- function(sample = "SAMPLE", replicate = 1L, mode = "U5", trial = "TRIAL"){
  data.table(trial = trial, subject = "SUBJECT", sample = sample,
             replicate = replicate, refGenome = "hs1", mode = mode,
             fragChromosome = "chr1", fragStrand = "+",
             fragStart = 100L, fragEnd = 199L, readID = paste0(sample, "_read"))
}
sql_key <- function(values){
  values <- as.character(unlist(values, use.names = FALSE))
  values[4L] <- as.character(as.integer(values[4L]))
  # Model common case-insensitive, trailing-space-insensitive SQL comparisons.
  values <- tolower(sub(" +$", "", values))
  paste0(nchar(values), ":", values, collapse = "|")
}
record_for <- function(x, file, total = 1L){
  data.frame(trial = as.character(x$trial[1L]), subject = as.character(x$subject[1L]),
             sample = as.character(x$sample[1L]), replicate = as.integer(as.character(x$replicate[1L])),
             ref_genome = as.character(x$refGenome[1L]), mode = as.character(x$mode[1L]),
             total_fragments = as.integer(total), data_file_name = file,
             stringsAsFactors = FALSE)
}
expect <- function(ok, message = "assertion failed"){
  if(!isTRUE(ok)) stop(message, call. = FALSE)
}

make_harness <- function(input = fixture(), overwrite = TRUE, faults = list()){
  s <- new.env(parent = emptyenv())
  s$root <- tempfile("buildFragmentsDB_tests_")
  dir.create(s$root)
  s$lake <- file.path(s$root, "data")
  s$tmp <- file.path(s$root, "tmp")
  dir.create(s$lake); dir.create(s$tmp)
  file.create(file.path(s$lake, ".inspiired"))
  s$input <- input
  s$committed <- empty_records()
  s$working <- NULL
  s$references <- list(sites = character(), fragments = character())
  s$events <- character()
  s$logs <- character()
  s$old_files <- character()
  s$in_transaction <- FALSE
  s$locked <- FALSE
  s$insert_count <- 0L
  s$reader_queries <- 0L
  s$writer <- structure(list(id = "writer"), class = "fake_connection")
  mark <- function(x) s$events <- c(s$events, x)
  fault <- function(name, default = FALSE) if(is.null(faults[[name]])) default else faults[[name]]
  write_payload <- function(x, path){
    # Deterministic bytes; the product still computes and verifies real MD5s.
    # A data.table split can reorder data.frame attributes without changing
    # the table. Serialize named columns so that irrelevant attribute ordering
    # does not manufacture a different mock Parquet payload.
    saveRDS(as.list(as.data.frame(x)), path, compress = FALSE, version = 2)
  }
  s$write_payload <- write_payload
  s$add_old <- function(x = input, content = "old fragment data", file = NULL){
    if(is.null(file)){
      p <- tempfile(tmpdir = s$lake, fileext = ".parquet")
      writeBin(charToRaw(content), p)
      file <- paste0(unname(tools::md5sum(p)), ".parquet")
      file.rename(p, file.path(s$lake, file))
    }
    s$committed <- rbind(s$committed, record_for(x, file))
    s$old_files <- unique(c(s$old_files, file))
    file
  }
  s$add_identical <- function(x = input){
    p <- tempfile(tmpdir = s$lake, fileext = ".parquet")
    write_payload(x, p)
    name <- paste0(unname(tools::md5sum(p)), ".parquet")
    file.rename(p, file.path(s$lake, name))
    s$add_old(x, file = name)
  }
  current <- function(conn){
    if(identical(conn$id, "writer") && s$in_transaction) s$working else s$committed
  }
  matching <- function(tab, params){
    wanted <- sql_key(params[seq_len(6L)])
    if(!nrow(tab)) return(integer())
    which(vapply(seq_len(nrow(tab)), function(i) sql_key(tab[i, key_names, drop = FALSE]), character(1)) == wanted)
  }
  mock <- list()
  mock$dbIsValid <- function(conn) TRUE
  mock$dbBegin <- function(conn){
    expect(!s$in_transaction, "nested transaction")
    mark("begin"); s$working <- s$committed; s$in_transaction <- TRUE
    if(isTRUE(fault("concurrent_duplicate"))){
      s$working <- rbind(s$working, record_for(input, "concurrent.parquet"))
      s$committed <- s$working
    }
    invisible(TRUE)
  }
  mock$dbRollback <- function(conn){
    mark("rollback"); s$working <- NULL; s$in_transaction <- FALSE
    invisible(TRUE)
  }
  mock$dbCommit <- function(conn){
    mark("commit_driver")
    if(isTRUE(fault("driver_commit_error"))) stop("injected driver commit failure")
    invisible(TRUE)
  }
  mock$dbDisconnect <- function(conn){ mark(paste0("disconnect_", conn$id)); invisible(TRUE) }
  mock$dbGetQuery <- function(conn, statement, params = list(), ...){
    sql <- as.character(statement)
    if(grepl("GET_LOCK", sql, fixed = TRUE)){
      mark("lock_acquire")
      if(isTRUE(fault("lock_denied"))) return(data.frame(acquired = 0L))
      s$locked <- TRUE
      return(data.frame(acquired = 1L))
    }
    if(grepl("RELEASE_LOCK", sql, fixed = TRUE)){
      mark("lock_release"); s$locked <- FALSE
      return(data.frame(released = 1L))
    }
    if(grepl("information_schema.TABLES", sql, fixed = TRUE))
      return(data.frame(engine = if(isTRUE(fault("wrong_engine"))) "MyISAM" else "InnoDB"))
    if(grepl("information_schema.COLUMNS", sql, fixed = TRUE)){
      cols <- c(key_names, "total_fragments", "processed_date", "data_file_name")
      return(data.frame(column_name = cols, max_length = c(100,100,100,NA,10,20,NA,NA,50)))
    }
    if(grepl("information_schema.KEY_COLUMN_USAGE", sql, fixed = TRUE))
      return(data.frame(column_name = key_names))
    if(grepl("data_file_name IN", sql, fixed = TRUE)){
      table <- if(grepl("FROM sites", sql, fixed = TRUE)) "sites" else "fragments"
      mark(paste0("refs_", table))
      if(identical(fault("reference_error", NULL), table)) stop("injected reference query failure")
      files <- s$references[[table]]
      if(table == "fragments") files <- c(files, s$committed$data_file_name)
      return(data.frame(data_file_name = intersect(files, as.character(unlist(params)))))
    }
    if(grepl("FROM fragments WHERE", sql, fixed = TRUE)){
      tab <- current(conn)
      found <- tab[matching(tab, params), , drop = FALSE]
      if(identical(conn$id, "reader")){
        mark("reader_select"); s$reader_queries <- s$reader_queries + 1L
        if(isTRUE(fault("readback_error"))) stop("injected readback failure")
        if(isTRUE(fault("readback_wrong_count")) && nrow(found))
          found$total_fragments <- found$total_fragments + 1L
      } else mark(if(grepl("FOR UPDATE", sql, fixed = TRUE)) "locked_select" else "writer_select")
      return(found[, c("data_file_name", "total_fragments"), drop = FALSE])
    }
    stop("Unhandled mock SQL query: ", sql)
  }
  mock$dbExecute <- function(conn, statement, params = list(), ...){
    sql <- as.character(statement)
    if(identical(trimws(sql), "COMMIT")){
      mark("commit_sql")
      if(identical(fault("commit", "ok"), "throw_before")) stop("injected COMMIT failure")
      if(!identical(fault("commit", "ok"), "silent_drop")) s$committed <- s$working
      s$in_transaction <- FALSE
      if(identical(fault("commit", "ok"), "throw_after")) stop("injected lost COMMIT reply")
      return(0L)
    }
    expect(s$in_transaction, "mutation outside a transaction")
    if(grepl("^DELETE FROM fragments", sql)){
      mark("delete_sql")
      indices <- matching(s$working, params)
      if(length(indices)) s$working <- s$working[-indices, , drop = FALSE]
      return(as.integer(length(indices)))
    }
    if(grepl("^INSERT INTO fragments", sql)){
      mark("insert_sql"); s$insert_count <- s$insert_count + 1L
      if(identical(as.integer(fault("insert_error_at", 0L)), s$insert_count))
        stop("injected INSERT failure")
      if(length(matching(s$working, params))) stop("duplicate SQL primary key")
      values <- params[seq_len(6L)]
      names(values) <- key_names
      row <- as.data.frame(values, stringsAsFactors = FALSE)
      row$replicate <- as.integer(row$replicate)
      row$total_fragments <- as.integer(params[[7L]])
      row$data_file_name <- as.character(params[[8L]])
      s$working <- rbind(s$working, row)
      return(1L)
    }
    stop("Unhandled mock SQL mutation: ", sql)
  }
  mock$write_parquet <- function(x, path, ...){
    mark("write_parquet")
    if(isTRUE(fault("write_error"))) stop("injected serialization failure")
    write_payload(x, path)
  }
  e <- new.env(parent = globalenv())
  e$mock <- mock
  provenance_fixture <- list(test = "immutable parent history")
  e$provenanceFromObject <- function(x) provenance_fixture
  e$provenanceWriteParquet <- function(x, path, provenance = NULL){
    expect(identical(provenance, provenance_fixture),
           "uploader did not explicitly preserve its input provenance")
    mock$write_parquet(x, path)
  }
  e$args <- list(dbConn = s$writer, dbConfigFile = file.path(s$root, "test.cnf"),
                 dbConfigID = "test", tmpDir = s$tmp, overwriteDBrecords = overwrite,
                 logFile = file.path(s$root, "test.log"))
  file.create(e$args$dbConfigFile)
  e$updateLog <- function(text, ...) s$logs <- c(s$logs, text)
  e$createDBconnection <- function(){
    mark("open_reader")
    if(isTRUE(fault("new_connection_error"))) stop("injected new connection failure")
    structure(list(id = "reader"), class = "fake_connection")
  }
  e$mock_unlink <- function(x, ...){
    obsolete <- x[dirname(x) == s$lake & basename(x) %in% basename(s$old_files)]
    if(length(obsolete)){
      mark("delete_old")
      expect("commit_sql" %in% s$events && s$reader_queries > 0L,
             "old file deleted before committed independent verification")
      if(isTRUE(fault("unlink_error"))) return(1L)
    }
    base::unlink(x, ...)
  }
  e$mock_file_link <- function(from, to){
    mark("publish_file")
    if(isTRUE(fault("link_error"))) return(FALSE)
    base::file.link(from, to)
  }
  rewrite <- function(node){
    if(missing(node)) return(quote(expr = ))
    if(is.character(node)){
      node[node == "/data"] <- s$lake
      return(node)
    }
    if(is.call(node)){
      if(identical(node[[1L]], as.name("::")) &&
         as.character(node[[2L]]) %in% c("DBI", "arrow"))
        return(call("$", as.name("mock"), node[[3L]]))
      if(identical(node[[1L]], as.name("requireNamespace")) &&
         length(node) >= 2L && identical(node[[2L]], "arrow")) return(TRUE)
      if(identical(node[[1L]], as.name("unlink"))) node[[1L]] <- as.name("mock_unlink")
      if(identical(node[[1L]], as.name("file.link"))) node[[1L]] <- as.name("mock_file_link")
      return(as.call(lapply(as.list(node), rewrite)))
    }
    if(is.pairlist(node)) return(as.pairlist(lapply(node, rewrite)))
    node
  }
  fn <- eval(rewrite(uploader_expression), envir = e)
  s$run <- function(){
    s$error <- tryCatch({
      suppressWarnings(suppressMessages(fn(input)))
      NULL
    }, error = function(err) conditionMessage(err))
    invisible(s$error)
  }
  s$exists <- function(name) file.exists(file.path(s$lake, name))
  s$files <- function() list.files(s$lake, pattern = "[.]parquet$")
  s$clean <- function() unlink(s$root, recursive = TRUE)
  s
}

tests_passed <- 0L
test <- function(name, code){
  force(name)
  tryCatch({force(code); tests_passed <<- tests_passed + 1L; cat("PASS:", name, "\n")},
           error = function(e) stop("FAIL: ", name, ": ", conditionMessage(e), call. = FALSE))
}
with_case <- function(code, ...){
  s <- make_harness(...)
  on.exit(s$clean(), add = TRUE)
  code(s)
}
expect_ok <- function(s){expect(is.null(s$error), paste("Unexpected error:", s$error)); expect(!s$locked, "lock leaked")}
expect_failed <- function(s, pattern = NULL){
  expect(is.character(s$error) && nzchar(s$error), "expected failure")
  if(!is.null(pattern)) expect(grepl(pattern, s$error, ignore.case = TRUE), paste("Unexpected error:", s$error))
  expect(!s$locked, "lock leaked after failure")
}

test("new keys insert without overwrite and preserve modes", with_case(function(s){
  s$run(); expect_ok(s)
  expect(nrow(s$committed) == 2L)
  expect(setequal(s$committed$mode, c("U3", "U5")))
  expect(all(s$committed$data_file_name %in% s$files()))
}, input = rbind(fixture(mode = "U3"), fixture(mode = "U5")), overwrite = FALSE))

test("default duplicate aborts the whole batch before any files are staged", with_case(function(s){
  old <- s$add_old(fixture("B")); before <- s$committed
  s$run(); expect_failed(s, "already exists|overwrite")
  expect(identical(before, s$committed)); expect(s$exists(old))
  expect(!any(c("write_parquet", "begin", "insert_sql", "delete_sql") %in% s$events))
}, input = rbind(fixture("A"), fixture("B")), overwrite = FALSE))

test("overwrite deletes only displaced unreferenced files after fresh readback", with_case(function(s){
  old <- s$add_old(); unrelated <- s$add_old(fixture("UNRELATED"), "unrelated")
  s$run(); expect_ok(s)
  expect(!s$exists(old)); expect(s$exists(unrelated)); expect(nrow(s$committed) == 2L)
  expect(which(s$events == "delete_old")[1L] > max(which(s$events == "reader_select")))
  expect(which(s$events == "lock_release")[1L] > which(s$events == "delete_old")[1L])
}))

test("byte-identical replacement reuses its file", with_case(function(s){
  old <- s$add_identical(); s$run(); expect_ok(s)
  expect(s$exists(old)); expect(identical(s$committed$data_file_name, old))
  expect(!"delete_old" %in% s$events); expect(!"publish_file" %in% s$events)
}))

for(table in c("sites", "fragments")) local({
  target <- table
  test(paste("references from", target, "retain the old file"), with_case(function(s){
    old <- s$add_old(); s$references[[target]] <- file.path(s$lake, old)
    s$run(); expect_ok(s); expect(s$exists(old)); expect(!"delete_old" %in% s$events)
  }))
})

for(outcome in c("throw_before", "throw_after", "silent_drop")) local({
  action <- outcome
  test(paste("uncertain COMMIT retains old and prepared files:", action), with_case(function(s){
    old <- s$add_old(); s$run(); expect_failed(s)
    expect(s$exists(old)); expect(length(s$files()) == 2L); expect(!"delete_old" %in% s$events)
    if(action != "throw_after") expect(identical(s$committed$data_file_name, old))
  }, faults = list(commit = action)))
})

for(failure in c("new_connection_error", "readback_error", "readback_wrong_count", "driver_commit_error")) local({
  selected <- failure
  test(paste("verification failure retains files:", selected), with_case(function(s){
    old <- s$add_old(); s$run(); expect_failed(s)
    expect(s$exists(old)); expect(length(s$files()) == 2L); expect(!"delete_old" %in% s$events)
  }, faults = setNames(list(TRUE), selected)))
})

test("SQL-equivalent input keys roll back without losing the original record", with_case(function(s){
  old <- s$add_old(fixture(trial = "TRIAL")); before <- s$committed
  s$run(); expect_failed(s, "duplicate")
  expect(identical(before, s$committed)); expect(s$exists(old)); expect("rollback" %in% s$events)
  expect(!"commit_sql" %in% s$events); expect(!"delete_old" %in% s$events)
}, input = rbind(fixture(trial = "Trial"), fixture(trial = "TRIAL"))))

test("SQL-equivalent new keys fail INSERT and roll back the entire batch", with_case(function(s){
  s$run(); expect_failed(s, "duplicate")
  expect(nrow(s$committed) == 0L); expect("rollback" %in% s$events)
  expect(sum(s$events == "insert_sql") == 2L)
  expect(!"commit_sql" %in% s$events); expect(!"delete_old" %in% s$events)
}, input = rbind(fixture(trial = "Trial"), fixture(trial = "TRIAL"))))

test("overwrite preserves other modes and replicates of the same sample", with_case(function(s){
  old <- s$add_old(fixture(mode = "U5", replicate = 1L), "old U5 replicate1")
  other_mode <- s$add_old(fixture(mode = "U3", replicate = 1L), "old U3 replicate1")
  other_rep <- s$add_old(fixture(mode = "U5", replicate = 2L), "old U5 replicate2")
  before <- s$committed[2:3, , drop = FALSE]
  s$run(); expect_ok(s)
  after <- s$committed[s$committed$mode == "U3" | s$committed$replicate == 2L, , drop = FALSE]
  rownames(before) <- rownames(after) <- NULL
  expect(identical(before, after)); expect(nrow(s$committed) == 3L)
  expect(!s$exists(old)); expect(s$exists(other_mode)); expect(s$exists(other_rep))
}))

test("an INSERT failure rolls back every replacement", with_case(function(s){
  old_a <- s$add_old(fixture("A"), "old A"); old_b <- s$add_old(fixture("B"), "old B")
  before <- s$committed; s$run(); expect_failed(s, "INSERT")
  expect(identical(before, s$committed)); expect(s$exists(old_a) && s$exists(old_b))
  expect(!"delete_old" %in% s$events)
}, input = rbind(fixture("A"), fixture("B")), faults = list(insert_error_at = 2L)))

test("an unsafe old path cannot remove a file outside the lake", with_case(function(s){
  outside <- file.path(s$root, "outside.parquet"); writeLines("keep me", outside)
  s$add_old(file = "../outside.parquet"); before <- s$committed
  s$run(); expect_failed(s, "unsafe|path|Parquet")
  expect(file.exists(outside)); expect(identical(before, s$committed)); expect(!"commit_sql" %in% s$events)
}))

test("an old symlink cannot remove its target", with_case(function(s){
  outside <- file.path(s$root, "outside.parquet"); writeLines("keep me", outside)
  link <- file.path(s$lake, "linked.parquet")
  expect(file.symlink(outside, link), "symlink fixture unavailable")
  s$add_old(file = basename(link)); before <- s$committed
  s$run(); expect_failed(s, "unsafe|path|symlink")
  expect(file.exists(outside)); expect(nzchar(Sys.readlink(link)))
  expect(identical(before, s$committed))
}))

test("reference-check failure retains old files after committing the new result", with_case(function(s){
  old <- s$add_old(); s$run(); expect_failed(s, "reference")
  expect(s$exists(old)); expect(s$committed$data_file_name != old)
  expect(!"delete_old" %in% s$events)
}, faults = list(reference_error = "sites")))

test("cleanup failure reports failure without undoing the committed replacement", with_case(function(s){
  old <- s$add_old(); s$run(); expect_failed(s, "remove|cleanup|delet")
  expect(s$exists(old)); expect(s$committed$data_file_name != old)
  expect(!"rollback" %in% s$events)
}, faults = list(unlink_error = TRUE)))

test("denied upload lock prevents publication and mutation", with_case(function(s){
  s$run(); expect_failed(s, "lock")
  expect(!any(c("write_parquet", "begin", "insert_sql") %in% s$events))
}, faults = list(lock_denied = TRUE)))

test("serialization failure leaves existing database and files untouched", with_case(function(s){
  old <- s$add_old(); before <- s$committed; s$run(); expect_failed(s, "serialization")
  expect(identical(before, s$committed)); expect(s$exists(old)); expect(!"begin" %in% s$events)
}, faults = list(write_error = TRUE)))

test("failed atomic publication does not change the database", with_case(function(s){
  old <- s$add_old(); before <- s$committed; s$run(); expect_failed(s, "publish|hard link")
  expect(identical(before, s$committed)); expect(s$exists(old)); expect(!"begin" %in% s$events)
  expect(length(list.files(s$lake, pattern = "[.]pending$", all.files = TRUE)) == 0L)
}, faults = list(link_error = TRUE)))

test("corrupt existing checksum-named content is never overwritten", with_case(function(s){
  old <- s$add_old(); before <- s$committed
  p <- tempfile(tmpdir = s$tmp)
  s$write_payload(s$input, p)
  corrupt <- paste0(unname(tools::md5sum(p)), ".parquet")
  writeLines("wrong bytes", file.path(s$lake, corrupt))
  old_hash <- unname(tools::md5sum(file.path(s$lake, corrupt)))
  s$run(); expect_failed(s, "checksum")
  expect(identical(before, s$committed)); expect(s$exists(old))
  expect(identical(unname(tools::md5sum(file.path(s$lake, corrupt))), old_hash))
  expect(!"begin" %in% s$events)
}))

test("concurrent duplicate is rejected without overwrite", with_case(function(s){
  s$run(); expect_failed(s, "concurrent|overwrite|exists")
  expect(identical(s$committed$data_file_name, "concurrent.parquet"))
  expect(!"delete_old" %in% s$events)
}, overwrite = FALSE, faults = list(concurrent_duplicate = TRUE)))

cat("\n", tests_passed, " fault-injection test groups passed.\n", sep = "")

#!/usr/bin/env -S Rscript --vanilla
# Revision: fixed-rows-reviewed-v4-2026-09-23
reportDBTestError <- function(e){
  cat(conditionMessage(e), "\n", sep = "", file = stderr())
  flush(stderr())
  quit(save = "no", status = 1, runLast = FALSE)
}

# The parser itself may warn while loading. Defer those diagnostics until the
# command line tells us whether verbose output was requested.
startupWarnings <- character()
tryCatch(withCallingHandlers(
  suppressPackageStartupMessages(library('argparse')),
  warning = function(w) {
    startupWarnings <<- c(startupWarnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
), error = reportDBTestError)

parser <- ArgumentParser()
parser$add_argument("--dbConfigFile", type = "character", default = 'none', help = "Path to db credential file.")
parser$add_argument("--dbConfigID",   type = "character", default = 'none', help = "DB credential block identifier in db credential file.")
parser$add_argument("--readData",    action = "store_true", default = FALSE, help = "Require an actual row to be read from both fragments and sites; empty tables fail.")
parser$add_argument("--writeData",   action = "store_true", default = FALSE, help = "Insert or refresh and commit one reusable test row per table; requires INSERT and UPDATE.")
parser$add_argument("--deleteData",  action = "store_true", default = FALSE, help = "Delete and commit marked test rows from both tables; fails if no test rows exist.")
parser$add_argument("--verbose",     action = "store_true", default = FALSE, help = "Print diagnostic messages; by default only errors are printed.")
parser$add_argument("--softwareRoot", type = "character", required = TRUE, help = "Path to INSPIIRED2 installation.")

# Include diagnostics from shared connection helpers and DBI cleanup. Errors
# still propagate to the handler below, which prints their text and exits 1.
withDBTestDiagnostics <- function(code){
  if(isTRUE(args$verbose)) return(force(code))
  suppressWarnings(suppressMessages(force(code)))
}

# This trial/subject/genome/mode combination is reserved for database tests.
# Permanent 'fauxtrial' read seeds from inspiired.sql use different keys and are
# never changed by these write/delete tests.
# Fixed keys bound the test footprint to one row per table. A fresh marker in
# data_file_name proves each write changes data, even within the same second.
# It is a diagnostic token, not a Parquet filename.
dbTestRows <- function(sample = "__dbtest__", token = NA_character_){
  common <- list(trial = "__INSPIIRED2_DBTEST__", subject = "__dbtest__",
                 sample = sample, ref_genome = "test", mode = "test")
  list(
    fragments = c(common, list(replicate = 1L, total_fragments = 123L,
                               data_file_name = token)),
    sites = c(common, list(total_sites = 45L, data_file_name = token))
  )
}

dbTestKey <- function(conn, row, includeSample = TRUE, includeToken = FALSE){
  keys <- intersect(c("trial", "subject", "sample", "replicate", "ref_genome", "mode"), names(row))
  if(!includeSample) keys <- setdiff(keys, "sample")
  if(includeToken) keys <- c(keys, "data_file_name")
  list(where = paste(paste(DBI::dbQuoteIdentifier(conn, keys), "= ?"), collapse = " AND "),
       params = unname(row[keys]))
}

withDBTestConnection <- function(code){
  conn <- createDBconnection()
  open <- TRUE
  on.exit({
    if(open) tryCatch(DBI::dbDisconnect(conn), error = function(e) {
      message("WARNING: disconnect failed: ", conditionMessage(e))
    })
  }, add = TRUE)
  result <- code(conn)
  DBI::dbDisconnect(conn)
  open <- FALSE
  invisible(result)
}

withDBTestLock <- function(code){
  # Keep a dedicated connection open across the stage connections and COMMITs.
  # Disconnect releases this advisory lock on both successful and failed tests.
  withDBTestConnection(function(conn){
    lock <- DBI::dbGetQuery(conn,
                            "SELECT GET_LOCK(CONCAT('INSPIIRED2.testDBconn:', MD5(DATABASE())), 10) AS acquired"
    )
    if(nrow(lock) != 1L || !isTRUE(lock$acquired[[1L]] == 1L)){
      stop("Could not acquire the database test lock within 10 seconds. ",
           "Check that a database is selected and retry after other tests finish.", call. = FALSE)
    }
    code()
  })
}

requireDBTestInnoDB <- function(conn){
  for(table in c("fragments", "sites")){
    engine <- DBI::dbGetQuery(conn, paste(
      "SELECT ENGINE AS engine FROM information_schema.TABLES",
      "WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?"
    ), params = list(table))
    if(nrow(engine) != 1L) stop(table, ": table is missing or not visible.", call. = FALSE)
    if(!isTRUE(toupper(engine$engine[[1L]]) == "INNODB")){
      stop(table, ": InnoDB is required for transactional tests.", call. = FALSE)
    }
  }
}

commitDBTestTransaction <- function(conn, code){
  DBI::dbBegin(conn)
  open <- TRUE
  on.exit({
    if(open) tryCatch(DBI::dbRollback(conn), error = function(e) {
      message("WARNING: rollback failed: ", conditionMessage(e))
    })
  }, add = TRUE)
  result <- code()
  # RMariaDB 1.3.5's dbCommit() ignores mysql_commit()'s error return. Execute
  # COMMIT through the checked query path first, then finish DBI bookkeeping.
  # This detects server/connection errors even without --readData privileges.
  tryCatch({
    DBI::dbExecute(conn, "COMMIT", immediate = TRUE)
    DBI::dbCommit(conn)
  }, error = function(e) {
    stop("COMMIT failed; its outcome may be uncertain: ", conditionMessage(e), call. = FALSE)
  })
  open <- FALSE
  invisible(result)
}

testWriteData <- function(conn, rows){
  requireDBTestInnoDB(conn)
  message("Test records: trial=", rows$sites$trial, "; sample=", rows$sites$sample)
  commitDBTestTransaction(conn, function(){
    for(table in names(rows)){
      tryCatch({
        row <- rows[[table]]
        columns <- DBI::dbQuoteIdentifier(conn, c(if(table == "fragments") "total_fragments" else "total_sites",
                                                  "data_file_name"))
        updates <- c(paste0(columns, " = VALUES(", columns, ")"), "processed_date = CURRENT_TIMESTAMP")
        query <- sprintf("INSERT INTO %s (%s) VALUES (%s) ON DUPLICATE KEY UPDATE %s",
                         DBI::dbQuoteIdentifier(conn, table),
                         paste(DBI::dbQuoteIdentifier(conn, names(row)), collapse = ", "),
                         paste(rep("?", length(row)), collapse = ", "),
                         paste(updates, collapse = ", "))
        affected <- DBI::dbExecute(conn, query, params = unname(row))
        if(length(affected) != 1L || is.na(affected) || !(affected %in% c(1L, 2L))){
          stop("Expected one inserted or updated test row (affected count 1 or 2).")
        }
      }, error = function(e) {
        stop(table, " INSERT/UPDATE: ", conditionMessage(e), call. = FALSE)
      })
    }
  })
  for(table in names(rows)) message("writeData passed for ", table, ": reusable test row committed.")
  invisible(TRUE)
}

testReadData <- function(conn, rows = NULL){
  templates <- if(is.null(rows)) dbTestRows() else rows
  failures <- character()
  for(table in names(templates)){
    tryCatch({
      row <- templates[[table]]
      query <- sprintf("SELECT %s FROM %s",
                       paste(DBI::dbQuoteIdentifier(conn, c(names(row), "processed_date")), collapse = ", "),
                       DBI::dbQuoteIdentifier(conn, table))
      if(is.null(rows)){
        actual <- DBI::dbGetQuery(conn, paste(query, "LIMIT 1"))
      } else {
        key <- dbTestKey(conn, row)
        actual <- DBI::dbGetQuery(conn, paste(query, "WHERE", key$where), params = key$params)
      }
      if(nrow(actual) == 0L){
        if(!is.null(rows)) stop("The committed test row was not returned by this new connection.")
        stop("No rows available to read. Have an administrator seed both tables using inspiired.sql, ",
             "or first run --writeData with writer credentials, then retry --readData with this account.")
      }
      if(nrow(actual) != 1L) stop("Expected exactly one row.")
      if(!is.null(rows)){
        matches <- vapply(names(row), function(column) {
          if(is.na(row[[column]])) is.na(actual[[column]][[1L]]) else
            isTRUE(actual[[column]][[1L]] == row[[column]])
        }, logical(1))
        if(!all(matches)) stop("Read-back mismatch in: ", paste(names(row)[!matches], collapse = ", "))
        if(is.na(actual$processed_date[[1L]])) stop("processed_date was not populated.")
      }
      message("readData passed for ", table, ": retrieved one row",
              if(is.null(rows)) "." else " and verified its committed values on a new connection.")
    }, error = function(e) {
      failures <<- c(failures, paste0(table, " SELECT: ", conditionMessage(e)))
    })
  }
  if(length(failures)) stop(paste(failures, collapse = "\n"), call. = FALSE)
  invisible(TRUE)
}

testDeleteData <- function(conn, rows = NULL){
  requireDBTestInnoDB(conn)
  templates <- if(is.null(rows)) dbTestRows() else rows
  counts <- commitDBTestTransaction(conn, function(){
    vapply(names(templates), function(table) {
      tryCatch({
        key <- dbTestKey(conn, templates[[table]], includeSample = !is.null(rows), includeToken = !is.null(rows))
        query <- paste("DELETE FROM", DBI::dbQuoteIdentifier(conn, table), "WHERE", key$where)
        affected <- DBI::dbExecute(conn, query, params = key$params)
        if(length(affected) != 1L || is.na(affected) || affected < 1L){
          stop("No marked test row was deleted. Run --writeData first, or include it in this command.")
        }
        if(!is.null(rows) && !isTRUE(affected == 1L)) stop("Expected exactly one deleted test row.")
        as.numeric(affected)
      }, error = function(e) {
        stop(table, " DELETE: ", conditionMessage(e), call. = FALSE)
      })
    }, numeric(1))
  })
  for(table in names(counts)){
    message("deleteData passed for ", table, ": committed deletion of ", counts[[table]], " test row(s).")
  }
  invisible(TRUE)
}

verifyDBTestDeletion <- function(conn, rows = NULL){
  templates <- if(is.null(rows)) dbTestRows() else rows
  for(table in names(templates)){
    key <- dbTestKey(conn, templates[[table]], includeSample = !is.null(rows), includeToken = !is.null(rows))
    query <- paste("SELECT trial FROM", DBI::dbQuoteIdentifier(conn, table), "WHERE", key$where, "LIMIT 1")
    remaining <- DBI::dbGetQuery(conn, query, params = key$params)
    if(nrow(remaining) != 0L) stop(table, ": test rows remain after committed deletion.", call. = FALSE)
  }
  message("Committed deletions verified on a new connection.")
  invisible(TRUE)
}

runDBTestStages <- function(){
  failures <- character()
  rows <- NULL
  writeCommitted <- FALSE
  attempt <- function(label, code){
    tryCatch(force(code), error = function(e) {
      failures <<- c(failures, paste0(label, ": ", conditionMessage(e)))
      invisible(FALSE)
    })
  }
  
  # Selected stages always run in write -> read -> delete order, each on a new
  # connection. Requested deletion is still attempted if read verification fails.
  if(isTRUE(args$writeData)){
    attempt("--writeData", withDBTestConnection(function(conn){
      token <- paste0("__dbtest__", DBI::dbGetQuery(conn, "SELECT UUID() AS id")$id[[1L]])
      rows <<- dbTestRows(token = token)
      testWriteData(conn, rows)
      writeCommitted <<- TRUE
    }))
  }
  if(isTRUE(args$readData)){
    attempt("--readData", withDBTestConnection(function(conn){
      # Verify this run's committed values after a successful write. If writing
      # failed, independently check SELECT access using any existing row instead
      # of reporting a missing write fixture as a read-access failure.
      testReadData(conn, if(writeCommitted) rows else NULL)
    }))
  }
  if(isTRUE(args$deleteData)){
    attempt("--deleteData", {
      # Match the new token as well as the key, so a failed refresh cannot delete
      # a previously committed fixture. No identifier means no safe cleanup scope.
      if(isTRUE(args$writeData) && is.null(rows)) stop("The write test did not establish a test identifier.")
      withDBTestConnection(function(conn) testDeleteData(conn, rows))
      if(isTRUE(args$readData)) withDBTestConnection(function(conn) verifyDBTestDeletion(conn, rows))
    })
  }
  if(length(failures)) stop("Database test(s) failed:\n", paste(failures, collapse = "\n"), call. = FALSE)
  message('All requested database tests passed.')
  invisible(TRUE)
}

runModule <- function(){
  if(args$dbConfigFile == 'none') stop('Error - provide --dbConfigFile.')
  if(args$dbConfigID == 'none') stop('Error - provide --dbConfigID.')
  if(!file.exists(args$dbConfigFile)) stop('Error - db config file not found: ', args$dbConfigFile)
  
  message('testDBconn revision: fixed-rows-reviewed-v4-2026-09-23')
  if(!any(c(isTRUE(args$readData), isTRUE(args$writeData), isTRUE(args$deleteData)))){
    withDBTestConnection(function(conn) message('Connection successful.'))
    return(invisible(TRUE))
  }
  withDBTestLock(runDBTestStages)
}

#-------------------------------------------------------------------------------

args <- parser$parse_args()

tryCatch({
  withDBTestDiagnostics({
    for(text in startupWarnings) warning(text, call. = FALSE, immediate. = TRUE)
    suppressPackageStartupMessages(library('RMariaDB'))
    source(file.path(args$softwareRoot, 'lib', 'common.R'))
    runModule()
  })
}, error = reportDBTestError)

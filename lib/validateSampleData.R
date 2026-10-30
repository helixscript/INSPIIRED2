


}
}

.vsd_connect <- function(args){
  if(!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("RMariaDB", quietly = TRUE))
    stop("Database checking requires the DBI and RMariaDB R packages.", call. = FALSE)
  DBI::dbConnect(RMariaDB::MariaDB(), group = args$dbConfigID, default.file = args$dbConfigFile)
}
.vsd_schema <- function(report, conn, limits){
  schema <- tryCatch(DBI::dbGetQuery(conn, paste(
    "SELECT TABLE_NAME AS table_name, COLUMN_NAME AS column_name, CHARACTER_MAXIMUM_LENGTH AS max_length",
    "FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()",
    "AND TABLE_NAME IN ('fragments', 'sites') AND COLUMN_NAME IN ('trial', 'subject', 'sample')")),
    error = function(e){.vsd_issue(report, NA, "database schema", conditionMessage(e)); NULL})
  if(is.null(schema)) return(limits)
  expected <- as.vector(outer(c("fragments", "sites"), names(limits), paste, sep = "."))
  actual <- paste(schema$table_name, schema$column_name, sep = ".")
  if(anyDuplicated(actual) || !setequal(expected, actual) || any(!is.finite(as.numeric(schema$max_length)) | as.numeric(schema$max_length) < 1)){
    .vsd_issue(report, NA, "database schema", "Could not read all trial/subject/sample character limits for fragments and sites. Check the selected database, table schema, and SELECT access to both tables.")
    return(limits)
  }
  for(column in names(limits)) limits[[column]] <- min(as.numeric(schema$max_length[schema$column_name == column]))
  report$limit_source <- "database schema (the smaller limit across fragments and sites)"
  limits
}
.vsd_duplicates <- function(report, sheet, integer_ok, conn){
  if(is.null(sheet$typed)) return(invisible(NULL))
  keys <- as.data.frame(lapply(sheet$typed[.vsd_keys], as.character), stringsAsFactors = FALSE)
  eligible <- integer_ok & Reduce(`&`, lapply(keys, .vsd_present))
  keys <- keys[eligible, , drop = FALSE]
  lines <- sheet$lines[eligible]
  keys$replicate <- as.integer(suppressWarnings(as.numeric(keys$replicate)))
  unique_rows <- which(!duplicated(keys))
  report$db_checked <- 0L
  query <- paste("SELECT trial, subject, sample, replicate, ref_genome AS refGenome, mode FROM fragments",
                 "WHERE trial = ? AND subject = ? AND sample = ? AND replicate = ? AND ref_genome = ? AND mode = ? LIMIT 1")
  for(i in unique_rows){
    found <- tryCatch(DBI::dbGetQuery(conn, query, params = unname(as.list(keys[i, , drop = FALSE]))),
                      error = function(e){.vsd_issue(report, lines[i], "database", paste0("Could not check fragment records: ", conditionMessage(e))); NULL})
    if(is.null(found)) return(invisible(NULL))
    report$db_checked <- report$db_checked + 1L
    if(nrow(found)){
      same <- Reduce(`&`, lapply(names(keys), function(column) keys[[column]] == keys[[column]][i]))
      description <- paste(paste(names(keys), vapply(keys[i, , drop = FALSE], as.character, character(1)), sep = "="), collapse = ", ")
      for(j in which(same)) .vsd_issue(report, lines[j], "database", paste0("Fragment record already exists: ", description,
                                                                            ". Exclude this sample replicate or resolve the existing record before processing."))
    }
  }
  report$db_complete <- TRUE
}
.vsd_database <- function(report, sheet, integer_ok, args, limits, connect){
  supplied <- c(args$dbConfigFile != "none", args$dbConfigID != "none")
  if(!any(supplied)) return(limits)
  report$db_requested <- TRUE
  if(!all(supplied) || any(!nzchar(trimws(c(args$dbConfigFile, args$dbConfigID))))){
    .vsd_issue(report, NA, "database flags", "Supply both --dbConfigFile and --dbConfigID, or omit both for validation without a database.")
    return(limits)
  }
  if(!file.exists(args$dbConfigFile) || dir.exists(args$dbConfigFile) || file.access(args$dbConfigFile, 4L) != 0L){
    .vsd_issue(report, NA, "dbConfigFile", "The credential file is missing or unreadable. Correct --dbConfigFile.")
    return(limits)
  }
  conn <- tryCatch(connect(args), error = function(e){
    .vsd_issue(report, NA, "database connection", paste0(conditionMessage(e), " Check the credential group, database name, host, and port.")); NULL
  })
  if(is.null(conn)) return(limits)
  open <- TRUE
  on.exit(if(open) try(DBI::dbDisconnect(conn), silent = TRUE), add = TRUE)
  limits <- .vsd_schema(report, conn, limits)
  if(!is.null(sheet)) .vsd_duplicates(report, sheet, integer_ok, conn)
  tryCatch({DBI::dbDisconnect(conn); open <- FALSE}, error = function(e)
    .vsd_issue(report, NA, "database connection", paste0("Could not close the connection cleanly: ", conditionMessage(e))))
  limits
}
.vsd_lengths <- function(report, sheet, limits){
  for(column in names(limits)){
    lengths <- nchar(sheet$raw[[column]], type = "chars")
    for(i in which(lengths > limits[[column]]))
      .vsd_issue(report, sheet$lines[i], column, paste0(.vsd_value(sheet$raw[[column]][i]), " has ", lengths[i],
                                                        " characters; the limit is ", limits[[column]], ". Shorten the name consistently across the study's sample files."))
  }
}

validateSampleData <- function(args, connect = .vsd_connect){
  report <- new.env(parent = emptyenv())
  report$issues <- list(); report$notes <- character(); report$n_rows <- 0L
  report$db_requested <- FALSE; report$db_complete <- FALSE; report$db_checked <- 0L
  report$limit_source <- "bundled schema (used when live limits are unavailable)"
  sheet <- .vsd_read(report, args$sampleData)
  integer_ok <- logical()
  if(!is.null(sheet)){
    integer_ok <- .vsd_cells(report, sheet)
    .vsd_sequence_pairs(report, sheet)
    .vsd_resources(report, sheet, args$softwareRoot, args$resourceDir)
  }
  limits <- .vsd_database(report, sheet, integer_ok, args,
                          c(trial = 100L, subject = 100L, sample = 100L), connect)
  if(!is.null(sheet)) .vsd_lengths(report, sheet, limits)
  report$limits <- limits
  as.list(report)
}

printSampleDataValidation <- function(report, path, verbose = FALSE, output = NULL){
  failed <- length(report$issues) > 0L
  if(is.null(output)) output <- if(failed) stderr() else stdout()
  writeLines(paste0(if(failed) "FAILED: " else "PASS: ", path, " (", report$n_rows, " sample rows; ",
                    length(report$issues), " issues)."), output)
  issues <- report$issues[order(vapply(report$issues, `[[`, numeric(1), "line"), na.last = FALSE)]
  for(issue in issues){
    location <- if(is.na(issue$line)) "File" else paste("Line", issue$line)
    writeLines(paste0("  ", location, " [", issue$column, "]: ", issue$message), output)
  }
  writeLines(paste0("Name limits: ", paste(paste(names(report$limits), report$limits, sep = "="), collapse = ", "),
                    " characters; ", report$limit_source, "."), output)
  database <- if(!report$db_requested) "skipped (no credentials supplied)" else if(report$db_complete)
    paste0("checked ", report$db_checked, " unique sample replicate keys") else "incomplete; resolve the reported input/connection/query issues"
  writeLines(paste0("Database: ", database, "."), output)
  if(verbose) writeLines(report$notes, output)
  if(failed) writeLines("Correct the reported issues, then rerun validateSampleData before starting the pipeline.", output)
  flush(output)
  invisible(!failed)
}


# Supporting library for modules/validateSampleData.R. Helpers keep parsing, resource checks,
# database checks, and the technician-facing report separate.
.vsd_columns <- c("trial", "subject", "sample", "replicate", "index1Seq",
                  "adriftReadLinkerSeq", "refGenome", "leaderSeqHMM", "vectorFastaFile", "mode")
.vsd_keys <- c("trial", "subject", "sample", "replicate", "refGenome", "mode")

.vsd_issue <- function(report, line, column, message){
  report$issues[[length(report$issues) + 1L]] <- list(line = line, column = column, message = message)
}
.vsd_value <- function(x){
  if(is.na(x)) return("<missing>")
  if(nchar(x) > 90L) x <- paste0(substr(x, 1L, 87L), "...")
  encodeString(x, quote = '"')
}
.vsd_present <- function(x) !is.na(x) & nzchar(trimws(x)) & x != "NA"

.vsd_record <- function(line){
  warnings <- character()
  value <- tryCatch(withCallingHandlers(utils::read.table(
    text = line, header = FALSE, sep = "\t", quote = '"', comment.char = "",
    colClasses = "character", na.strings = character(), strip.white = FALSE,
    blank.lines.skip = FALSE, fill = FALSE, check.names = FALSE, row.names = NULL),
    warning = function(w){warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
    error = identity)
  if(inherits(value, "error") || length(warnings) || nrow(value) != 1L)
    return(NULL)
  unname(unlist(value[1L, ], use.names = FALSE))
}

.vsd_read <- function(report, path){
  if(!file.exists(path) || dir.exists(path) || file.access(path, 4L) != 0L){
    .vsd_issue(report, NA, "file", "Cannot read the sample file. Check --sampleData and file permissions.")
    return(NULL)
  }
  warnings <- character()
  lines <- tryCatch(withCallingHandlers(readLines(path, warn = FALSE, encoding = "UTF-8"),
    warning = function(w){warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
    error = identity)
  if(inherits(lines, "error") || length(warnings) || any(!validUTF8(lines))){
    .vsd_issue(report, NA, "file", "Cannot read clean UTF-8 text. Re-export as UTF-8 tab-delimited text, not an Excel workbook or UTF-16 file.")
    return(NULL)
  }
  if(!length(lines) || !nzchar(lines[1L])){
    .vsd_issue(report, 1L, "header", "The file is empty or has no header row.")
    return(NULL)
  }
  lines[1L] <- sub("^\ufeff", "", lines[1L])
  header <- .vsd_record(lines[1L])
  if(is.null(header) || length(header) < 2L){
    .vsd_issue(report, 1L, "header", "Expected a tab-delimited header. Check the separator and quotation marks; export TSV rather than comma- or space-delimited text.")
    return(NULL)
  }
  if(any(!nzchar(header)) || anyDuplicated(header) || any(header != trimws(header))){
    .vsd_issue(report, 1L, "header", "Column names must be non-empty, unique, and have no surrounding spaces.")
    return(NULL)
  }
  missing <- setdiff(.vsd_columns, header)
  if(length(missing)){
    .vsd_issue(report, 1L, "header", paste0("Missing required columns: ", paste(missing, collapse = ", "), ". Names are case-sensitive; column order may vary."))
    return(NULL)
  }
  report$n_rows <- length(lines) - 1L
  if(!report$n_rows){
    .vsd_issue(report, NA, "file", "The header is present, but there are no sample rows.")
    return(NULL)
  }
  rows <- list(); line_numbers <- integer()
  for(i in seq.int(2L, length(lines))){
    cells <- if(nzchar(trimws(lines[i]))) .vsd_record(lines[i]) else NULL
    if(is.null(cells)){
      .vsd_issue(report, i, "TSV", "Blank row or malformed quotation marks. Use one sample per line and remove blank rows/unclosed quotes.")
    } else if(length(cells) != length(header)){
      .vsd_issue(report, i, "TSV", paste0("Found ", length(cells), " fields; the header has ", length(header),
        ". Check for missing/extra tabs or shifted cells. Keep empty cells between their tabs."))
    } else {
      rows[[length(rows) + 1L]] <- cells
      line_numbers <- c(line_numbers, i)
    }
  }
  if(!length(rows)) return(NULL)
  raw <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE, check.names = FALSE)
  names(raw) <- header
  # Match demultiplex's readr type inference for SQL identifiers (e.g. numeric
  # sample names), while retaining raw cells for corrections and length checks.
  typed <- tryCatch(suppressWarnings(readr::read_tsv(I(paste(c(lines[1L], lines[line_numbers]), collapse = "\n")),
    show_col_types = FALSE, name_repair = "minimal", progress = FALSE, num_threads = 1)), error = identity)
  if(inherits(typed, "error")){
    .vsd_issue(report, NA, "TSV", paste0("The pipeline's TSV reader failed: ", conditionMessage(typed)))
    typed <- NULL
  } else {
    problems <- readr::problems(typed)
    for(i in seq_len(nrow(problems))){
      row <- problems$row[i] - 1L # readr counts the header as row 1.
      line <- if(row %in% seq_along(line_numbers)) line_numbers[row] else NA
      column <- if(problems$col[i] %in% seq_along(header)) header[problems$col[i]] else "TSV"
      .vsd_issue(report, line, column, paste0("The pipeline's TSV reader expected ", problems$expected[i],
        " but found ", .vsd_value(problems$actual[i]), ". Check the column and use consistent value types."))
    }
    if(nrow(typed) != nrow(raw))
      .vsd_issue(report, NA, "TSV", "The pipeline's TSV reader found a different row count. Check tabs and quotation marks.")
    if(nrow(problems) || nrow(typed) != nrow(raw)) typed <- NULL
  }
  list(raw = raw, typed = typed, lines = line_numbers)
}

.vsd_cells <- function(report, sheet){
  d <- sheet$raw
  for(column in .vsd_columns){
    for(i in which(!.vsd_present(d[[column]])))
      .vsd_issue(report, sheet$lines[i], column, "A value is required. Blank cells and literal NA are read as missing by the pipeline.")
    for(i in which(.vsd_present(d[[column]]) & (d[[column]] != trimws(d[[column]]) | grepl("[\t\r\n]", d[[column]]))))
      .vsd_issue(report, sheet$lines[i], column, "Remove surrounding whitespace or embedded tabs/newlines; check that this value is in the correct column.")
  }
  reps <- suppressWarnings(as.numeric(d$replicate))
  integer_ok <- is.finite(reps) & reps == trunc(reps) & abs(reps) <= .Machine$integer.max
  for(i in which(.vsd_present(d$replicate) & !integer_ok))
    .vsd_issue(report, sheet$lines[i], "replicate", paste0(.vsd_value(d$replicate[i]),
      " must be a whole number between -2147483647 and 2147483647 (for example, 1 or 2)."))
  for(i in which(.vsd_present(d$adriftReadLinkerSeq) & !grepl("^[ACGT]{3,}N{5,}[ACGT]{3,}$", d$adriftReadLinkerSeq, ignore.case = TRUE)))
    .vsd_issue(report, sheet$lines[i], "adriftReadLinkerSeq", paste0(.vsd_value(d$adriftReadLinkerSeq[i]),
      " must have one run of at least 5 Ns, with at least 3 A/C/G/T bases on each side (for example ACGNNNNNTGA). Check the linker/UMI sequence."))
  for(i in which(.vsd_present(d$index1Seq) & !grepl("^[ACGTN]+$", d$index1Seq, ignore.case = TRUE)))
    .vsd_issue(report, sheet$lines[i], "index1Seq", paste0(.vsd_value(d$index1Seq[i]), " is not a DNA barcode (A/C/G/T/N). Check for a shifted column."))
  for(i in which(.vsd_present(d$mode) & !d$mode %in% c("U3", "U5")))
    .vsd_issue(report, sheet$lines[i], "mode", paste0(.vsd_value(d$mode[i]), " must be U3 or U5 for the current demultiplex reader. Check for a shifted column."))
  integer_ok
}

.vsd_sequence_pairs <- function(report, sheet){
  columns <- c("index1Seq", "adriftReadLinkerSeq")
  present <- .vsd_present(sheet$raw$index1Seq) & .vsd_present(sheet$raw$adriftReadLinkerSeq)
  pairs <- as.data.frame(lapply(sheet$raw[present, columns, drop = FALSE], toupper), stringsAsFactors = FALSE)
  lines <- sheet$lines[present]
  repeated <- unique(pairs[duplicated(pairs), , drop = FALSE])
  for(i in seq_len(nrow(repeated))){
    matches <- which(pairs$index1Seq == repeated$index1Seq[i] &
                     pairs$adriftReadLinkerSeq == repeated$adriftReadLinkerSeq[i])
    message <- paste0("Duplicate barcode/linker combination on TSV lines ", paste(lines[matches], collapse = ", "),
      ": index1Seq=", .vsd_value(repeated$index1Seq[i]),
      "; adriftReadLinkerSeq=", .vsd_value(repeated$adriftReadLinkerSeq[i]),
      ". Remove the repeated row or correct the barcode/linker assignments so each pair occurs once. Letter case is ignored.")
    for(j in matches) .vsd_issue(report, lines[j], "index1Seq + adriftReadLinkerSeq", message)
  }
}

.vsd_resources <- function(report, sheet, softwareRoot, resourceDir){
  bundled <- file.path(softwareRoot, "data")
  if(resourceDir != "/resources" && !dir.exists(resourceDir))
    .vsd_issue(report, NA, "resources", paste0("Resource overlay directory does not exist: ", resourceDir, ". Correct --resourceDir."))
  definitions <- list(leaderSeqHMM = c("hmms", "[.]hmm$"),
                      refGenome = c("referenceGenomes", "[.]2bit$"), vectorFastaFile = c("vectors", ""))
  for(column in names(definitions)){
    definition <- definitions[[column]]
    roots <- unique(file.path(c(resourceDir, bundled), definition[1L]))
    # A read-only view of resource_overlay(): external files take precedence.
    files <- unlist(lapply(roots, function(path) list.files(path, pattern = definition[2L], full.names = TRUE)), use.names = FALSE)
    files <- files[file.exists(files) & !dir.exists(files)]
    files <- files[!duplicated(basename(files))]
    names(files) <- if(column == "refGenome") sub("[.]2bit$", "", basename(files)) else basename(files)
    report$notes <- c(report$notes, paste0(column, ": searched ", paste(roots, collapse = " and ")))
    if(!length(files)){
      .vsd_issue(report, NA, column, paste0("No resource files found in ", paste(roots, collapse = " or "),
        ". Install/mount the resources or correct --resourceDir."))
      next
    }
    values <- sheet$raw[[column]]
    for(value in unique(values[.vsd_present(values)])){
      if(!value %in% names(files)){
        choices <- names(files)[order(utils::adist(value, names(files))[1L, ])]
        hint <- paste(head(choices, 3L), collapse = ", ")
        description <- if(column == "refGenome") "Use the genome name without .2bit" else "Use the exact filename, including its extension"
        msg <- paste0(.vsd_value(value), " was not found. ", description, ". Available names to check: ", hint, ".")
      } else if(file.access(files[[value]], 4L) != 0L || !isTRUE(file.size(files[[value]]) > 0)){
        msg <- paste0("Resource ", files[[value]], " is empty or unreadable. Ask the informatician to check it.")
      } else next
      for(i in which(values == value)) .vsd_issue(report, sheet$lines[i], column, msg)
    }
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

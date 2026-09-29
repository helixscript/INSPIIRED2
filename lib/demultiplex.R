## Golay correction

## Nucleotide->bit mappings (decode allows N; encode does not)
DECODE_GOLAY_NT_TO_BITS <- c(A="11", C="00", T="10", G="01", N="11")
ENCODE_GOLAY_NT_TO_BITS <- c(A="11", C="00", T="10", G="01")

## Parity submatrix P (12x12), as in the original code
DEFAULT_P <- matrix(
  c(
    0,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,0,1,1,1,0,0,0,1,0,
    1,1,0,1,1,1,0,0,0,1,0,1,
    1,0,1,1,1,0,0,0,1,0,1,1,
    1,1,1,1,0,0,0,1,0,1,1,0,
    1,1,1,0,0,0,1,0,1,1,0,1,
    1,1,0,0,0,1,0,1,1,0,1,1,
    1,0,0,0,1,0,1,1,0,1,1,1,
    1,0,0,1,0,1,1,0,1,1,1,0,
    1,0,1,0,1,1,0,1,1,1,0,0,
    1,1,0,1,1,0,1,1,1,0,0,0,
    1,0,1,1,0,1,1,1,0,0,0,1
  ),
  nrow = 12, byrow = TRUE
)

## Generator (G) and parity-check (H) matrices
DEFAULT_G <- cbind(DEFAULT_P, diag(12))
DEFAULT_H <- cbind(diag(12), t(DEFAULT_P))

## Build all error vectors with Hamming weight <= 3 (length 24)
.make_3bit_errors <- function(veclen = 24L) {
  errs <- list()
  idx <- 1L
  z <- rep(0L, veclen)
  errs[[idx]] <- z; idx <- idx + 1L
  
  for (i in 1L:veclen) {
    v <- z; v[i] <- 1L; errs[[idx]] <- v; idx <- idx + 1L
  }
  for (i in 1L:(veclen-1L)) for (j in (i+1L):veclen) {
    v <- z; v[i] <- 1L; v[j] <- 1L; errs[[idx]] <- v; idx <- idx + 1L
  }
  for (i in 1L:(veclen-2L)) for (j in (i+1L):(veclen-1L)) for (k in (j+1L):veclen) {
    v <- z; v[i] <- 1L; v[j] <- 1L; v[k] <- 1L; errs[[idx]] <- v; idx <- idx + 1L
  }
  errs
}

## Syndrome lookup table: key = 12-bit syndrome as a string, value = 24-bit error vec
.build_syndrome_lut <- function(H) {
  lut <- new.env(parent = emptyenv())
  for (err in .make_3bit_errors(24L)) {
    syn <- as.integer((H %*% err) %% 2)
    key <- paste(syn, collapse = "")
    lut[[key]] <- err
  }
  lut
}

.DEFAULT_SYNDROME_LUT <- .build_syndrome_lut(DEFAULT_H)

## --- Helpers: sequence<->bits ------------------------------------------------

.seq_to_bits <- function(seq, nt_to_bits = DECODE_GOLAY_NT_TO_BITS) {
  nts <- strsplit(toupper(seq), "", fixed = TRUE)[[1]]
  bits_each <- vapply(nts, function(n) {
    v <- nt_to_bits[[n]]
    if (is.null(v)) stop(sprintf("Invalid nucleotide '%s' in '%s'", n, seq), call. = FALSE)
    v
  }, character(1))
  bit_chars <- strsplit(paste(bits_each, collapse = ""), "", fixed = TRUE)[[1]]
  as.integer(bit_chars)
}

.bits_to_seq <- function(bits, nt_to_bits = ENCODE_GOLAY_NT_TO_BITS) {
  if (length(bits) %% 2L != 0L) stop("Bit vector length must be even.", call. = FALSE)
  bits_to_nt <- setNames(names(nt_to_bits), nt_to_bits)  # e.g., "11" -> "A"
  pairs <- vapply(seq.int(1L, length(bits), by = 2L),
                  function(i) paste0(bits[i], bits[i+1]), character(1))
  paste0(bits_to_nt[pairs], collapse = "")
}

## --- Core: decode a 24-bit received vector via syndrome decoding -------------

.decode_bits <- function(rec_bits, H = DEFAULT_H, lut = .DEFAULT_SYNDROME_LUT) {
  syn <- as.integer((H %*% rec_bits) %% 2)
  key <- paste(syn, collapse = "")
  err <- lut[[key]]
  if (is.null(err)) {
    return(list(corrected = NULL, num_errors = 4L))  # uncorrectable (likely 4-bit)
  }
  corrected <- as.integer((rec_bits + err) %% 2)
  list(corrected = corrected, num_errors = sum(err))
}

## --- Public API: decode Golay 12-nt sequences --------------------------------

#' Decode 12-nt DNA Golay barcodes (vectorized).
#'
#' @param x Character string, character vector, or Biostrings::DNAString/ DNAStringSet.
#' @param nt_to_bits Optional named character vector mapping nts to 2-bit strings.
#'                   Defaults allow N during decoding ("N" -> "11").
#' @return Data frame with columns:
#'         input, corrected (NA if uncorrectable), num_bit_errors (int),
#'         uncorrectable (logical).
correctGolay12 <- function(x, nt_to_bits = DECODE_GOLAY_NT_TO_BITS) {
  # Accept DNAString / DNAStringSet transparently if Biostrings is present
  if (inherits(x, "DNAString") || inherits(x, "DNAStringSet")) x <- as.character(x)
  if (length(x) == 0L) return(data.frame(input=character(), corrected=character(),
                                         num_bit_errors=integer(), uncorrectable=logical(),
                                         stringsAsFactors = FALSE))
  if (!is.character(x)) stop("Input must be character or DNAString(Set).", call. = FALSE)
  
  decode_one <- function(s) {
    if (is.na(s) || nchar(s) == 0L) {
      return(c(corrected = NA_character_, num_bit_errors = NA_integer_, uncorrectable = TRUE))
    }
    if (nchar(s) != 12L) {
      return(c(corrected = NA_character_, num_bit_errors = NA_integer_, uncorrectable = TRUE))
    }
    rec_bits <- .seq_to_bits(s, nt_to_bits = nt_to_bits)
    out <- .decode_bits(rec_bits, H = DEFAULT_H, lut = .DEFAULT_SYNDROME_LUT)
    if (is.null(out$corrected)) {
      c(corrected = NA_character_, num_bit_errors = 4L, uncorrectable = TRUE)
    } else {
      corrected_seq <- .bits_to_seq(out$corrected, nt_to_bits = ENCODE_GOLAY_NT_TO_BITS)
      c(corrected = corrected_seq, num_bit_errors = out$num_errors, uncorrectable = FALSE)
    }
  }
  
  res <- t(vapply(x, decode_one,
                  FUN.VALUE = c(corrected = "", num_bit_errors = 0, uncorrectable = FALSE)))
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  res$num_bit_errors <- as.integer(res$num_bit_errors)
  res$uncorrectable <- as.logical(res$uncorrectable)
  data.frame(input = x, res, row.names = NULL, check.names = FALSE)
}


# Sample keys and reports do not require database packages. Connections below are
# local to each phase and are closed before this module starts parallel workers.
.demultiplexKeyColumns <- c("trial", "subject", "sample", "replicate", "refGenome", "mode")
.demultiplexDBColumns <- c("trial", "subject", "sample", "replicate", "ref_genome", "mode")

.demultiplexKeyLabel <- function(keys, i){
  paste(paste0(names(keys), "=", vapply(keys, function(x)
    encodeString(as.character(x[[i]]), quote = '"'), character(1))), collapse = ", ")
}

validateDemultiplexSampleKeys <- function(sampleData){
  if(!is.data.frame(sampleData) || !nrow(sampleData) || anyDuplicated(names(sampleData)) ||
     !all(.demultiplexKeyColumns %in% names(sampleData)))
    stop("Error - sample data must contain rows and unique columns including trial, subject, sample, replicate, refGenome and mode.", call. = FALSE)
  keys <- as.data.frame(sampleData)[.demultiplexKeyColumns]
  limits <- c(trial = 100L, subject = 100L, sample = 100L, refGenome = 10L, mode = 20L)
  for(column in names(limits)){
    values <- keys[[column]]
    if(!is.atomic(values) || !is.null(dim(values)) || is.complex(values) || anyNA(values))
      stop("Error - sample key ", column, " must contain non-missing scalar identifiers.", call. = FALSE)
    values <- as.character(values)
    if(any(!nzchar(values)) || any(values != trimws(values)) || any(nchar(values, type = "chars") > limits[[column]]))
      stop("Error - sample key ", column, " must be non-empty, have no surrounding whitespace and be at most ", limits[[column]], " characters.", call. = FALSE)
    keys[[column]] <- values
  }
  replicate <- keys$replicate
  if(!is.atomic(replicate) || !is.null(dim(replicate)) || is.logical(replicate) || is.complex(replicate))
    stop("Error - replicate must contain whole numbers from -2147483647 to 2147483647.", call. = FALSE)
  replicate <- suppressWarnings(as.numeric(as.character(replicate)))
  if(anyNA(replicate) || any(!is.finite(replicate)) || any(replicate != trunc(replicate)) ||
     any(replicate < -2147483647 | replicate > 2147483647))
    stop("Error - replicate must contain whole numbers from -2147483647 to 2147483647.", call. = FALSE)
  keys$replicate <- as.integer(replicate)
  duplicate <- which(duplicated(keys))
  if(length(duplicate))
    stop("Error - duplicate sample key: ", .demultiplexKeyLabel(keys, duplicate[[1L]]), call. = FALSE)
  keys
}

.demultiplexCounts <- function(values){
  if(!is.atomic(values) || !is.null(dim(values)) || is.logical(values) || is.complex(values))
    stop("Error - demultiplexedReads must contain non-negative whole counts.", call. = FALSE)
  counts <- if(is.numeric(values) && !inherits(values, "integer64")) as.numeric(values) else
    suppressWarnings(as.numeric(as.character(values)))
  if(anyNA(counts) || any(!is.finite(counts)) || any(counts < 0 | counts > 2^53 - 1) ||
     any(counts != trunc(counts)))
    stop("Error - demultiplexedReads must contain non-negative whole counts no greater than 2^53 - 1.", call. = FALSE)
  counts
}

.demultiplexKeyIDs <- function(keys){
  # Length prefixes keep arbitrary identifier punctuation unambiguous.
  do.call(paste0, lapply(keys, function(x){
    x <- as.character(x)
    paste0(nchar(x, type = "bytes"), ":", x)
  }))
}

buildDemultiplexSampleReport <- function(sampleData, counts = NULL){
  keys <- validateDemultiplexSampleKeys(sampleData)
  report <- as.data.frame(sampleData)
  report[.demultiplexKeyColumns] <- keys
  report$demultiplexedReads <- rep(0, nrow(report))
  if(is.null(counts)) return(report)
  if(!is.data.frame(counts)) stop("Error - demultiplex counts must be a data frame.", call. = FALSE)
  if(!nrow(counts)) return(report)
  count_keys <- validateDemultiplexSampleKeys(counts)
  if(!"demultiplexedReads" %in% names(counts))
    stop("Error - demultiplex counts are missing demultiplexedReads.", call. = FALSE)
  values <- .demultiplexCounts(counts$demultiplexedReads)
  positions <- match(.demultiplexKeyIDs(count_keys), .demultiplexKeyIDs(keys))
  if(anyNA(positions))
    stop("Error - counts contain an unregistered sample key: ",
         .demultiplexKeyLabel(count_keys, which(is.na(positions))[[1L]]), call. = FALSE)
  report$demultiplexedReads[positions] <- values
  report
}

.withDemultiplexDB <- function(code){
  if(!is.list(args) || !all(c("dbConfigFile", "dbConfigID") %in% names(args)) ||
     any(vapply(args[c("dbConfigFile", "dbConfigID")], function(x)
       !is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x) || x == "none", logical(1))) ||
     !file.exists(args$dbConfigFile))
    stop("Error - sample database operations require args$dbConfigID and an existing args$dbConfigFile.", call. = FALSE)
  conn <- DBI::dbConnect(RMariaDB::MariaDB(), group = args$dbConfigID, default.file = args$dbConfigFile)
  on.exit(DBI::dbDisconnect(conn), add = TRUE)
  code(conn)
}

.checkDemultiplexDBSchema <- function(conn){
  engine <- DBI::dbGetQuery(conn, paste(
    "SELECT ENGINE AS engine FROM information_schema.TABLES",
    "WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'samples'"))
  if(nrow(engine) != 1L || !isTRUE(toupper(engine$engine[[1L]]) == "INNODB"))
    stop("Error - samples must exist as an InnoDB table; check its definition in inspiired2_dbSetup.sql.", call. = FALSE)
  schema <- DBI::dbGetQuery(conn, paste(
    "SELECT COLUMN_NAME AS column_name, DATA_TYPE AS data_type,",
    "CHARACTER_MAXIMUM_LENGTH AS max_length, IS_NULLABLE AS nullable",
    "FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'samples'"))
  required <- c(.demultiplexDBColumns, "demultiplexed_reads", "processed_date")
  if(!all(required %in% schema$column_name))
    stop("Error - samples is missing required columns.", call. = FALSE)
  schema <- schema[match(required, schema$column_name), , drop = FALSE]
  if(!identical(as.character(schema$data_type), c("varchar", "varchar", "varchar", "int", "varchar", "varchar", "bigint", "timestamp")) ||
     any(as.numeric(schema$max_length[c(1L, 2L, 3L, 5L, 6L)]) < c(100, 100, 100, 10, 20)) ||
     any(schema$nullable[seq_along(.demultiplexDBColumns)] != "NO") || schema$nullable[[7L]] != "YES")
    stop("Error - samples column types, key widths or count nullability do not match the required schema.", call. = FALSE)
  primary <- DBI::dbGetQuery(conn, paste(
    "SELECT COLUMN_NAME AS column_name FROM information_schema.KEY_COLUMN_USAGE",
    "WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'samples' AND CONSTRAINT_NAME = 'PRIMARY'",
    "ORDER BY ORDINAL_POSITION"))
  if(!identical(as.character(primary$column_name), .demultiplexDBColumns))
    stop("Error - samples must have the six-column primary key (trial, subject, sample, replicate, ref_genome, mode).", call. = FALSE)
}

.demultiplexDBWhere <- paste(paste0("`", .demultiplexDBColumns, "` = ?"), collapse = " AND ")
.demultiplexDBSelect <- paste("SELECT CAST(demultiplexed_reads AS CHAR) AS demultiplexed_reads, processed_date FROM samples WHERE", .demultiplexDBWhere)

.verifyDemultiplexDBRows <- function(conn, keys, counts = NULL){
  for(i in seq_len(nrow(keys))){
    actual <- DBI::dbGetQuery(conn, .demultiplexDBSelect, params = unname(as.list(keys[i, ])))
    valid <- nrow(actual) == 1L && !is.na(actual$processed_date[[1L]])
    if(valid) valid <- if(is.null(counts)) is.na(actual$demultiplexed_reads[[1L]]) else
      identical(actual$demultiplexed_reads[[1L]], sprintf("%.0f", counts[[i]]))
    if(!valid) stop("Error - samples verification failed for ", .demultiplexKeyLabel(keys, i), call. = FALSE)
  }
  invisible(TRUE)
}

.writeDemultiplexDBRows <- function(sampleData, completion = FALSE){
  keys <- validateDemultiplexSampleKeys(sampleData)
  counts <- NULL
  if(completion){
    if(!"demultiplexedReads" %in% names(sampleData))
      stop("Error - sample report is missing demultiplexedReads.", call. = FALSE)
    counts <- .demultiplexCounts(sampleData$demultiplexedReads)
  }
  .withDemultiplexDB(function(conn){
    .checkDemultiplexDBSchema(conn)
    DBI::dbBegin(conn)
    transaction_open <- TRUE
    on.exit({
      if(transaction_open) tryCatch(DBI::dbRollback(conn), error = function(e)
        warning("Sample transaction rollback failed: ", conditionMessage(e), call. = FALSE))
    }, add = TRUE)
    for(i in seq_len(nrow(keys))){
      params <- unname(as.list(keys[i, ]))
      tryCatch({
        if(completion){
          affected <- DBI::dbExecute(conn, paste(
            "UPDATE samples SET demultiplexed_reads = ?, processed_date = CURRENT_TIMESTAMP WHERE", .demultiplexDBWhere,
            "AND demultiplexed_reads IS NULL"), params = c(list(counts[[i]]), params))
          if(!isTRUE(affected == 1L))
            stop("Expected one registered sample with a NULL count; the key is missing or already completed.")
        } else {
          if(nrow(DBI::dbGetQuery(conn, .demultiplexDBSelect, params = params)))
            stop("Sample key already exists; existing sample records cannot be reset or overwritten.")
          affected <- DBI::dbExecute(conn, paste(
            "INSERT INTO samples (trial, subject, sample, replicate, ref_genome, mode, demultiplexed_reads)",
            "VALUES (?, ?, ?, ?, ?, ?, NULL)"), params = params)
          if(!isTRUE(affected == 1L)) stop("Expected exactly one inserted sample.")
        }
      }, error = function(e) stop("Error - samples ", if(completion) "completion" else "registration",
                                  " for ", .demultiplexKeyLabel(keys, i), ": ", conditionMessage(e), call. = FALSE))
    }
    .verifyDemultiplexDBRows(conn, keys, counts)
    # Older RMariaDB releases ignore mysql_commit()'s error return. Check the
    # SQL COMMIT first, then finish DBI's transaction bookkeeping.
    tryCatch({
      DBI::dbExecute(conn, "COMMIT", immediate = TRUE)
      DBI::dbCommit(conn)
    }, error = function(e) stop("Error - sample COMMIT failed; its outcome may be uncertain: ", conditionMessage(e), call. = FALSE))
    transaction_open <- FALSE
  })
  # A fresh connection proves registrations are committed before workers start,
  # and completed counts remain visible after the writing connection closes.
  .withDemultiplexDB(function(conn) .verifyDemultiplexDBRows(conn, keys, counts))
  invisible(TRUE)
}

registerSamplesInDB <- function(sampleData) .writeDemultiplexDBRows(sampleData)
updateSampleCountsInDB <- function(sampleData) .writeDemultiplexDBRows(sampleData, completion = TRUE)


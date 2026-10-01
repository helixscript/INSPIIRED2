memoryMonitorArgs <- c("disableMemoryMonitor", "memorySampleSeconds", "memoryPlotUpdateSeconds")


addMemoryMonitorArgs <- function(parser){
  parser$add_argument("--disableMemoryMonitor", action = "store_true", default = FALSE,
                      help = "Disable container memory monitoring and memory-use reports.")
  parser$add_argument("--memorySampleSeconds", type = "double", default = 2,
                      help = "Seconds between container memory measurements.")
  parser$add_argument("--memoryPlotUpdateSeconds", type = "double", default = 30,
                      help = "Seconds between updates of the memory-use PNG and PDF reports.")
  invisible(parser)
}


memoryMonitorID <- function(prefix){
  paste0(prefix, "_", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "_",
         Sys.getpid(), "_", sprintf("%06d", sample.int(999999L, 1L)))
}


atomicWriteLines <- function(text, path){
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  writeLines(text, tmp, useBytes = TRUE)
  if(!file.rename(tmp, path)){
    unlink(tmp)
    stop("could not atomically write ", path, call. = FALSE)
  }
  invisible(path)
}


startMemoryMonitor <- function(args, pipelineRoot){
  if(isTRUE(args$disableMemoryMonitor) || is.null(args$outputDir)) return(NULL)

  sampleSeconds <- args$memorySampleSeconds
  plotSeconds <- args$memoryPlotUpdateSeconds
  if(length(sampleSeconds) != 1L || !is.finite(sampleSeconds) || sampleSeconds <= 0)
    stop("Error - memorySampleSeconds must be one finite number greater than zero.", call. = FALSE)
  if(length(plotSeconds) != 1L || !is.finite(plotSeconds) || plotSeconds <= 0)
    stop("Error - memoryPlotUpdateSeconds must be one finite number greater than zero.", call. = FALSE)

  outputDir <- normalizePath(args$outputDir, mustWork = FALSE)
  monitorDir <- file.path(outputDir, "memory")
  controlDir <- file.path(monitorDir, "control")
  if(!dir.exists(monitorDir)) dir.create(monitorDir, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(controlDir)) dir.create(controlDir, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(monitorDir) || !dir.exists(controlDir)){
    warning("Container memory monitoring disabled: could not create ", controlDir, ".",
            call. = FALSE)
    return(NULL)
  }

  runID <- memoryMonitorID("run")
  invocationID <- memoryMonitorID(args$module)
  controlPrefix <- file.path(controlDir, invocationID)
  paths <- list(
    trace = file.path(monitorDir, "INSPIIRED2_memoryUsage.tsv"),
    png = file.path(monitorDir, "INSPIIRED2_memoryUsage.png"),
    pdf = file.path(monitorDir, "INSPIIRED2_memoryUsage.pdf"),
    summary = file.path(monitorDir, "INSPIIRED2_memoryUsage_summary.tsv"),
    ready = paste0(controlPrefix, ".ready"),
    stop = paste0(controlPrefix, ".stop"),
    status = paste0(controlPrefix, ".status"),
    done = paste0(controlPrefix, ".done"),
    error = paste0(controlPrefix, ".error"),
    log = paste0(controlPrefix, ".log")
  )

  cleanText <- function(x) gsub("[\t\r\n]", "_", as.character(x), perl = TRUE)
  script <- file.path(pipelineRoot, "tools", "memoryMonitor.sh")
  if(!file.exists(script)){
    warning("Container memory monitoring disabled: sampler script was not found.", call. = FALSE)
    return(NULL)
  }

  scriptArgs <- c(
    shQuote(script),
    "--trace", shQuote(paths$trace), "--png", shQuote(paths$png),
    "--pdf", shQuote(paths$pdf), "--summary", shQuote(paths$summary),
    "--run-id", shQuote(runID), "--invocation-id", shQuote(invocationID),
    "--module", shQuote(cleanText(args$module)),
    "--file-tag", shQuote(cleanText(if(is.null(args$fileTag)) args$module else args$fileTag)),
    "--sample-seconds", shQuote(format(sampleSeconds, scientific = FALSE, trim = TRUE)),
    "--plot-seconds", shQuote(as.character(max(1L, ceiling(plotSeconds)))),
    "--parent-pid", shQuote(as.character(Sys.getpid())),
    "--ready-file", shQuote(paths$ready), "--stop-file", shQuote(paths$stop),
    "--status-file", shQuote(paths$status), "--done-file", shQuote(paths$done),
    "--error-file", shQuote(paths$error)
  )

  launched <- tryCatch({
    status <- suppressWarnings(system2("bash", args = scriptArgs, stdout = FALSE,
                                       stderr = paths$log, wait = FALSE))
    is.null(status) || (length(status) == 1L && !is.na(status) && as.integer(status) == 0L)
  }, error = function(e) FALSE)

  if(!launched){
    warning("Container memory monitoring disabled: the sampler could not be launched.",
            call. = FALSE)
    return(NULL)
  }

  deadline <- Sys.time() + 5
  while(!file.exists(paths$ready) && !file.exists(paths$done) && Sys.time() < deadline)
    Sys.sleep(0.05)

  if(!file.exists(paths$ready)){
    file.create(paths$stop)
    stopDeadline <- Sys.time() + max(2, sampleSeconds + 1)
    while(!file.exists(paths$done) && Sys.time() < stopDeadline) Sys.sleep(0.05)
    detail <- if(file.exists(paths$error)) paste(readLines(paths$error, warn = FALSE), collapse = " ") else "sampler did not become ready"
    if(file.exists(paths$done)){
      unlink(c(paths$ready, paths$stop, paths$status, paths$done, paths$error), force = TRUE)
      if(file.exists(paths$log) && isTRUE(file.size(paths$log) == 0L)) unlink(paths$log)
    }
    warning("Container memory monitoring disabled: ", detail, ".", call. = FALSE)
    return(NULL)
  }

  ### runFile <- file.path(monitorDir, "currentRunID.txt")
  ### recordedRunID <- tryCatch(trimws(readLines(runFile, n = 1L, warn = FALSE)), error = function(e) "")

  runFile <- file.path(controlDir, "currentRunID.txt")
  recordedRunID <- if(file.exists(runFile)){
    tryCatch(trimws(readLines(runFile, n = 1L, warn = FALSE)),
             error = function(e) "")
  } else ""


  if(length(recordedRunID) == 1L && nzchar(recordedRunID)) runID <- recordedRunID

  ### message("Container memory monitoring enabled. Live plot: ", paths$png, "; PDF report: ", paths$pdf)
  c(paths, list(active = TRUE, sampleSeconds = sampleSeconds, runID = runID,
                invocationID = invocationID))
}


stopMemoryMonitor <- function(monitor, status){
  if(is.null(monitor) || !isTRUE(monitor$active)) return(invisible(NULL))

  status <- suppressWarnings(as.integer(status))
  if(length(status) != 1L || is.na(status)) status <- 1L

  try(atomicWriteLines(as.character(status), monitor$status), silent = TRUE)
  try(file.create(monitor$stop), silent = TRUE)

  deadline <- Sys.time() + max(10, monitor$sampleSeconds + 5)
  while(!file.exists(monitor$done) && Sys.time() < deadline) Sys.sleep(0.05)
  if(!file.exists(monitor$done)){
    warning("The memory sampler did not confirm shutdown; it will exit when the launcher process ends.",
            call. = FALSE)
    # Preserve stop/status controls so a renderer that is still finishing can
    # record the module's exit status after this launcher exits.
    return(invisible(NULL))
  }
  if(file.exists(monitor$error)){
    detail <- paste(readLines(monitor$error, warn = FALSE), collapse = " ")
    warning("Container memory monitoring warning: ", detail, ".", call. = FALSE)
  }

  unlink(c(monitor$ready, monitor$stop, monitor$status, monitor$done, monitor$error),
         force = TRUE)
  if(file.exists(monitor$log) && isTRUE(file.size(monitor$log) == 0L)) unlink(monitor$log)
  invisible(NULL)
}

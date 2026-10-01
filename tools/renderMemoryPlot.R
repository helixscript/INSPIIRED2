#!/usr/bin/env Rscript

# Live INSPIIRED2 memory plot and PDF report. Plot settings are intentionally
# kept together so the report can be restyled without changing its data logic.
PLOT_WIDTH <- 12
PLOT_HEIGHT <- 7.2
PNG_DPI <- 150
MAX_PLOT_POINTS <- 2400L
TOTAL_MEMORY_COLOR <- "#8f8f8f"
BOUNDARY_COLOR <- "#bdbdbd"
LIMIT_COLOR <- "#b2182b"

MODULE_COLORS <- c(
  demultiplex = "#1b9e77", prepReads = "#d95f02", alignReads = "#7570b3",
  buildFragments = "#e7298a", buildStdFragments = "#66a61e",
  buildSites = "#e6ab02", nearestGenes = "#a6761d",
  annotateRepeats = "#1f78b4", testHMMs = "#6a3d9a",
  buildSeqDataMap = "#b15928", validateSampleData = "#666666"
)
FALLBACK_COLORS <- c("#17becf", "#e41a1c", "#377eb8", "#4daf4a",
                     "#984ea3", "#ff7f00", "#a65628", "#f781bf")

fail <- function(...){
  message("Error - unable to render the memory report: ", paste0(..., collapse = ""))
  quit(save = "no", status = 1L)
}

require_package <- function(package){
  if(!requireNamespace(package, quietly = TRUE))
    fail("required R package '", package, "' is not installed.")
}

fread_strict <- function(...){
  withCallingHandlers(data.table::fread(...),
                      warning = function(w) stop(conditionMessage(w), call. = FALSE))
}

as_finite_numeric <- function(x, name, allow_na = FALSE){
  original_missing <- is.na(x) | trimws(as.character(x)) == ""
  suppressWarnings(value <- as.numeric(x))
  if(any(!original_missing & is.na(value)))
    fail("column '", name, "' contains a non-numeric value.")
  if(!allow_na && anyNA(value)) fail("column '", name, "' contains a missing value.")
  if(any(!is.na(value) & !is.finite(value)))
    fail("column '", name, "' contains a non-finite value.")
  value
}

module_palette <- function(modules){
  modules <- unique(as.character(modules))
  known <- intersect(modules, names(MODULE_COLORS))
  extra <- setdiff(modules, names(MODULE_COLORS))
  colors <- MODULE_COLORS[known]
  if(length(extra)){
    extra_colors <- rep(FALLBACK_COLORS, length.out = length(extra))
    names(extra_colors) <- extra
    colors <- c(colors, extra_colors)
  }
  colors[modules]
}

memory_theme <- function(){
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10, color = "#555555", lineheight = 1.15),
      plot.caption = ggplot2::element_text(size = 9, color = "#555555", hjust = 0),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom", legend.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(12, 16, 8, 12)
    )
}

validate_output <- function(path, label){
  size <- if(file.exists(path)) file.info(path)$size else NA_real_
  if(is.na(size) || size == 0) stop("the ", label, " device did not produce an output file.")
}

render_png <- function(plot, path){
  if(!isTRUE(capabilities("cairo")))
    stop("R was built without Cairo support, which is required for headless PNG output.")
  
  device_open <- FALSE
  on.exit(if(device_open) try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  
  grDevices::png(filename = path, width = PLOT_WIDTH, height = PLOT_HEIGHT,
                 units = "in", res = PNG_DPI, bg = "white", type = "cairo-png")
  device_open <- TRUE
  print(plot)
  grDevices::dev.off()
  device_open <- FALSE
  validate_output(path, "PNG")
}

format_duration <- function(seconds){
  seconds <- max(0, round(seconds))
  if(seconds < 60) return(sprintf("%ds", seconds))
  if(seconds < 3600) return(sprintf("%dm %02ds", seconds %/% 60, seconds %% 60))
  sprintf("%dh %02dm %02ds", seconds %/% 3600, (seconds %% 3600) %/% 60, seconds %% 60)
}

format_gib <- function(bytes){
  if(!length(bytes) || is.na(bytes)) return("not set")
  sprintf("%.2f GiB", bytes / 1024^3)
}

shorten <- function(x, width){
  x <- as.character(x)
  ifelse(nchar(x) <= width, x, paste0(substr(x, 1L, width - 3L), "..."))
}

draw_summary_pages <- function(d, summary_data, run_id, status, limit_gib){
  grid::grid.newpage()
  grid::grid.text("INSPIIRED2 container memory-use summary",
                  x = grid::unit(0.05, "npc"), y = grid::unit(0.95, "npc"),
                  just = c("left", "top"),
                  gp = grid::gpar(fontfamily = "sans", fontsize = 18, fontface = "bold"))

  if(!nrow(d)){
    grid::grid.text(paste("Run", run_id, "contains no memory measurements."),
                    x = grid::unit(0.05, "npc"), y = grid::unit(0.86, "npc"),
                    just = c("left", "top"), gp = grid::gpar(fontfamily = "sans", fontsize = 11))
    return(invisible(NULL))
  }

  duration <- max(d$epoch_seconds) - min(d$epoch_seconds)
  limit_text <- if(is.finite(limit_gib)) sprintf("%.2f GiB", limit_gib) else "not set"
  peak_percent <- if(is.finite(limit_gib)) sprintf("%.2f%%", 100 * max(d$total_gib) / limit_gib) else "not available"
  overview <- c(
    paste0("Run: ", run_id), paste0("Status: ", status),
    paste0("Start: ", d$timestamp_utc[1L],
           "    Last measurement: ", d$timestamp_utc[nrow(d)]),
    paste0("Elapsed: ", format_duration(duration), "    Measurements: ", nrow(d)),
    paste0("Peak working set: ", format_gib(max(d$working_set_bytes)),
           "    Peak total memory: ", format_gib(max(d$memory_bytes))),
    paste0("Peak anonymous memory: ", format_gib(max(d$anon_bytes)),
           "    Peak shared memory: ", format_gib(max(d$shmem_bytes))),
    paste0("Docker memory limit: ", limit_text, "    Peak percent of limit: ", peak_percent),
    paste0("OOM events: ", max(d$oom_events), "    OOM-kill events: ", max(d$oom_kill_events)),
    paste0("Report updated: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  )
  grid::grid.text(paste(overview, collapse = "\n"),
                  x = grid::unit(0.05, "npc"), y = grid::unit(0.87, "npc"),
                  just = c("left", "top"),
                  gp = grid::gpar(fontfamily = "sans", fontsize = 10.5, lineheight = 1.25))

  if(!nrow(summary_data)) return(invisible(NULL))
  labels <- ifelse(summary_data$module == summary_data$file_tag, summary_data$module,
                   paste0(summary_data$module, " [", summary_data$file_tag, "]"))
  exit_status <- trimws(as.character(summary_data$exit_status))
  exit_status[is.na(summary_data$exit_status) | !nzchar(exit_status)] <- "running"
  rows <- sprintf("%-28s %12s %7d %11.2f %11.2f %10.2f %11.2f %8s",
                  shorten(labels, 28L), vapply(summary_data$duration_seconds, format_duration, character(1)),
                  as.integer(summary_data$samples), summary_data$peak_working_set_bytes / 1024^3,
                  summary_data$peak_memory_bytes / 1024^3, summary_data$peak_anon_bytes / 1024^3,
                  summary_data$peak_shmem_bytes / 1024^3, shorten(exit_status, 8L))
  header <- sprintf("%-28s %12s %7s %11s %11s %10s %11s %8s",
                    "Module [file tag]", "Duration", "Samples", "Work GiB", "Total GiB",
                    "Anon GiB", "Shmem GiB", "Exit")

  rows_per_page <- 14L
  groups <- split(seq_along(rows), ceiling(seq_along(rows) / rows_per_page))
  for(page in seq_along(groups)){
    if(page > 1L){
      grid::grid.newpage()
      grid::grid.text("INSPIIRED2 module memory-use summary (continued)",
                      x = grid::unit(0.05, "npc"), y = grid::unit(0.95, "npc"),
                      just = c("left", "top"),
                      gp = grid::gpar(fontfamily = "sans", fontsize = 16, fontface = "bold"))
    }
    top <- if(page == 1L) 0.54 else 0.86
    grid::grid.text(paste(c(header, rows[groups[[page]]]), collapse = "\n"),
                    x = grid::unit(0.05, "npc"), y = grid::unit(top, "npc"),
                    just = c("left", "top"),
                    gp = grid::gpar(fontfamily = "mono", fontsize = 8.3, lineheight = 1.22))
  }
  invisible(NULL)
}

render_pdf <- function(plot, path, d, summary_data, run_id, status, limit_gib){
  device_open <- FALSE
  on.exit(if(device_open) try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  grDevices::pdf(file = path, width = PLOT_WIDTH, height = PLOT_HEIGHT, onefile = TRUE,
                 family = "sans", bg = "white", paper = "special", useDingbats = FALSE)
  device_open <- TRUE
  print(plot)
  draw_summary_pages(d, summary_data, run_id, status, limit_gib)
  grDevices::dev.off()
  device_open <- FALSE
  validate_output(path, "PDF")
}

build_empty_plot <- function(run_id){
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = "Waiting for memory measurements...",
                      hjust = 0, size = 4.5, color = "#333333") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(-0.5, 0.5) +
    ggplot2::labs(title = "INSPIIRED2 container memory use", subtitle = paste("Run", run_id)) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 18, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 10, color = "#555555"),
                   plot.margin = ggplot2::margin(12, 16, 8, 12))
}

read_module_summary <- function(path, target_run){
  if(!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0)
    fail("memory summary does not exist or is empty: ", path)
  required <- c("run_id", "invocation_id", "module", "file_tag", "start_utc", "end_utc",
                "duration_seconds", "samples", "peak_memory_bytes", "peak_working_set_bytes",
                "peak_anon_bytes", "peak_shmem_bytes", "memory_limit_bytes",
                "peak_percent_limit", "exit_status")
  x <- tryCatch(fread_strict(path, sep = "\t", na.strings = c("", "NA"),
                             showProgress = FALSE, nThread = 1L),
                error = function(e) fail(conditionMessage(e)))
  missing_cols <- setdiff(required, names(x))
  if(length(missing_cols))
    fail("memory summary is missing column(s): ", paste(missing_cols, collapse = ", "))
  x <- x[run_id == target_run]
  if(!nrow(x)) return(x)
  for(name in c("duration_seconds", "samples", "peak_memory_bytes", "peak_working_set_bytes",
                "peak_anon_bytes", "peak_shmem_bytes"))
    x[, (name) := as_finite_numeric(get(name), name)]
  x
}

main <- function(){
  cli <- commandArgs(trailingOnly = TRUE)
  if(length(cli) != 5L)
    fail("expected: renderMemoryPlot.R TRACE_TSV SUMMARY_TSV OUTPUT_PNG OUTPUT_PDF RUN_ID")
  trace_file <- cli[1]; summary_file <- cli[2]; output_png <- cli[3]
  output_pdf <- cli[4]; target_run <- cli[5]
  if(!nzchar(target_run)) fail("RUN_ID must not be blank.")

  require_package("data.table")
  require_package("ggplot2")
  if(!file.exists(trace_file)) fail("memory trace does not exist: ", trace_file)
  trace_size <- file.info(trace_file)$size
  if(is.na(trace_size)) fail("memory trace cannot be read: ", trace_file)
  if(trace_size == 0) fail("memory trace is empty: ", trace_file)

  required <- c("run_id", "invocation_id", "timestamp_utc", "epoch_seconds", "module",
                "file_tag", "event", "exit_status", "memory_bytes", "working_set_bytes",
                "anon_bytes", "shmem_bytes", "swap_bytes", "memory_limit_bytes",
                "oom_events", "oom_kill_events")
  header_line <- readLines(trace_file, n = 1L, warn = FALSE)
  if(length(header_line) != 1L || !nzchar(header_line)) fail("memory trace has no header.")
  header_names <- strsplit(sub("\r$", "", header_line), "\t", fixed = TRUE)[[1L]]
  if(anyDuplicated(header_names)) fail("memory trace contains duplicated column names.")
  missing_cols <- setdiff(required, header_names)
  if(length(missing_cols))
    fail("memory trace is missing column(s): ", paste(missing_cols, collapse = ", "))

  d <- tryCatch(fread_strict(trace_file, sep = "\t", select = required,
                             na.strings = c("", "NA"), showProgress = FALSE, nThread = 1L),
                error = function(e) fail(conditionMessage(e)))
  d <- d[run_id == target_run]
  summary_data <- read_module_summary(summary_file, target_run)
  if(!nrow(d)){
    p <- build_empty_plot(target_run)
    render_png(p, output_png)
    render_pdf(p, output_pdf, d, summary_data, target_run, "waiting", NA_real_)
    return(invisible(NULL))
  }

  if(anyNA(d$run_id) || anyNA(d$invocation_id) || anyNA(d$timestamp_utc) ||
     anyNA(d$module) || anyNA(d$file_tag) || anyNA(d$event))
    fail("the selected run contains a missing identifier, timestamp, module, file tag, or event.")
  if(any(!nzchar(trimws(as.character(d$invocation_id)))) ||
     any(!nzchar(trimws(as.character(d$module)))))
    fail("the selected run contains a blank invocation identifier or module name.")
  if(any(!d$event %in% c("start", "sample", "end")))
    fail("column 'event' contains a value other than start, sample, or end.")

  d[, source_order := .I]
  numeric_cols <- c("epoch_seconds", "memory_bytes", "working_set_bytes", "anon_bytes",
                    "shmem_bytes", "swap_bytes", "oom_events", "oom_kill_events")
  for(name in numeric_cols) d[, (name) := as_finite_numeric(get(name), name)]
  d[, memory_limit_bytes := as_finite_numeric(memory_limit_bytes, "memory_limit_bytes",
                                               allow_na = TRUE)]
  check_cols <- c(numeric_cols, "memory_limit_bytes")
  if(any(vapply(d[, ..check_cols], function(x) any(x < 0, na.rm = TRUE), logical(1))))
    fail("memory trace contains a negative time, memory measurement, or event count.")
  if(any(d$working_set_bytes > d$memory_bytes))
    fail("working-set memory exceeds total memory in the selected run.")
  if(any(d[, data.table::uniqueN(module), by = invocation_id]$V1 != 1L))
    fail("an invocation identifier is associated with more than one module.")

  data.table::setorder(d, epoch_seconds, source_order)
  gib <- 1024^3
  d[, `:=`(total_gib = memory_bytes / gib, working_gib = working_set_bytes / gib)]
  start_epoch <- min(d$epoch_seconds)
  duration_seconds <- max(d$epoch_seconds) - start_epoch
  use_hours <- duration_seconds >= 3600
  time_divisor <- if(use_hours) 3600 else 60
  time_unit <- if(use_hours) "hours" else "minutes"
  d[, elapsed_time := (epoch_seconds - start_epoch) / time_divisor]

  starts <- d[, .SD[1L], by = invocation_id]
  data.table::setorder(starts, epoch_seconds, source_order)
  module_levels <- unique(as.character(starts$module))
  d[, module := factor(module, levels = module_levels)]
  colors <- module_palette(module_levels)

  d[, row_in_invocation := seq_len(.N), by = invocation_id]
  d[, invocation_rows := .N, by = invocation_id]
  stride <- max(1L, ceiling(nrow(d) / MAX_PLOT_POINTS))
  d[, plot_order := seq_len(.N)]
  plot_data <- d[row_in_invocation == 1L | row_in_invocation == invocation_rows |
                   (plot_order - 1L) %% stride == 0L]
  line_data <- plot_data[, if(.N > 1L) .SD else NULL, by = invocation_id]
  boundaries <- if(nrow(starts) > 1L) starts[-1L] else starts[0L]

  peak_working <- max(d$working_gib)
  peak_total <- max(d$total_gib)
  peak <- d[which.max(working_gib)]
  peak[, label := sprintf("%.2f GiB", working_gib)]
  peak[, label_hjust := if(elapsed_time > max(d$elapsed_time) * 0.82) 1.08 else -0.08]

  limits <- d[!is.na(memory_limit_bytes) & memory_limit_bytes > 0, memory_limit_bytes / gib]
  limit_gib <- if(length(limits)) tail(limits, 1L) else NA_real_
  reference_max <- max(peak_total, 0.25)
  draw_limit <- is.finite(limit_gib) && limit_gib <= reference_max * 1.2

  last <- d[.N]
  last_status <- if(is.na(last$exit_status)) "" else trimws(as.character(last$exit_status))
  status <- if(last$event != "end") "updating" else if(nzchar(last_status))
    paste0("completed (exit ", last_status, ")") else "interrupted/status unavailable"
  subtitle <- sprintf("Run %s | %s | working-set peak %.2f GiB | total peak %.2f GiB",
                      target_run, status, peak_working, peak_total)
  if(is.finite(limit_gib)) subtitle <- paste0(subtitle, sprintf(" | limit %.2f GiB", limit_gib))

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = elapsed_time)) +
    ggplot2::geom_vline(data = boundaries, ggplot2::aes(xintercept = elapsed_time),
                        inherit.aes = FALSE, color = BOUNDARY_COLOR, linewidth = 0.35,
                        linetype = "dashed")
  if(nrow(line_data))
    p <- p + ggplot2::geom_line(data = line_data,
                                ggplot2::aes(y = total_gib, group = invocation_id),
                                color = TOTAL_MEMORY_COLOR, linewidth = 0.55)

  p <- p +
    ggplot2::geom_point(ggplot2::aes(y = working_gib, color = module), size = 1.25, alpha = 0.82) +
    ggplot2::geom_point(data = peak, ggplot2::aes(y = working_gib, color = module),
                        shape = 21, fill = "white", size = 3.4, stroke = 0.8,
                        show.legend = FALSE) +
    ggplot2::geom_text(data = peak, ggplot2::aes(y = working_gib, label = label,
                                                 hjust = label_hjust),
                       vjust = 1.5, size = 3.2, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = colors, breaks = module_levels, drop = FALSE) +
    ggplot2::guides(color = ggplot2::guide_legend(reverse = FALSE, byrow = TRUE)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.05))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.10))) +
    ggplot2::labs(
      title = "INSPIIRED2 container memory use", subtitle = subtitle,
      x = paste0("Elapsed wall-clock time (", time_unit, ")"), y = "Container memory (GiB)",
      caption = paste("Colored points: Docker-like working set.",
                      "Grey line: total cgroup memory, including cache and shared memory.")
    ) + memory_theme()

  if(draw_limit){
    p <- p +
      ggplot2::geom_hline(yintercept = limit_gib, color = LIMIT_COLOR,
                          linewidth = 0.45, linetype = "dashed") +
      ggplot2::annotate("text", x = Inf, y = limit_gib, label = "memory limit",
                        hjust = 1.05, vjust = -0.35, size = 3, color = LIMIT_COLOR)
  }
  render_png(p, output_png)
  render_pdf(p, output_pdf, d, summary_data, target_run, status, limit_gib)
  invisible(NULL)
}

tryCatch(main(), error = function(e) fail(conditionMessage(e)))

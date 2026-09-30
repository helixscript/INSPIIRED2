#!/usr/bin/env Rscript

# Live INSPIIRED2 memory plot. The settings in this block are intentionally
# kept together so the plot can be restyled without touching its data logic.
PLOT_WIDTH <- 12
PLOT_HEIGHT <- 7.2
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
  message("Error - unable to render the memory plot: ", paste0(..., collapse = ""))
  quit(save = "no", status = 1L)
}

require_package <- function(package){
  if(!requireNamespace(package, quietly = TRUE)) fail("required R package '", package, "' is not installed.")
}

fread_strict <- function(...){
  withCallingHandlers(data.table::fread(...),
                      warning = function(w) stop(conditionMessage(w), call. = FALSE))
}

as_finite_numeric <- function(x, name, allow_na = FALSE){
  original_missing <- is.na(x) | trimws(as.character(x)) == ""
  suppressWarnings(value <- as.numeric(x))
  if(any(!original_missing & is.na(value))) fail("column '", name, "' contains a non-numeric value.")
  if(!allow_na && anyNA(value)) fail("column '", name, "' contains a missing value.")
  if(any(!is.finite(value), na.rm = TRUE)) fail("column '", name, "' contains a non-finite value.")
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

render_svg <- function(plot, path){
  device_open <- FALSE
  on.exit(if(device_open) try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  grDevices::svg(filename = path, width = PLOT_WIDTH, height = PLOT_HEIGHT,
                 onefile = FALSE, family = "sans", bg = "white")
  device_open <- TRUE
  print(plot)
  grDevices::dev.off()
  device_open <- FALSE
  if(!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0)
    stop("the SVG device did not produce an output file.")
}

build_empty_plot <- function(run_id){
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = "Waiting for memory measurements...",
                      hjust = 0, size = 4.5, color = "#333333") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(-0.5, 0.5) +
    ggplot2::labs(title = "INSPIIRED2 container memory use",
                  subtitle = paste("Run", run_id)) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 18, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 10, color = "#555555"),
                   plot.margin = ggplot2::margin(12, 16, 8, 12))
}

main <- function(){
  cli <- commandArgs(trailingOnly = TRUE)
  if(length(cli) != 3L) fail("expected: renderMemoryPlot.R TRACE_TSV OUTPUT_SVG RUN_ID")
  trace_file <- cli[1]; output_svg <- cli[2]; target_run <- cli[3]
  if(!nzchar(target_run)) fail("RUN_ID must not be blank.")

  require_package("data.table")
  require_package("ggplot2")
  if(!file.exists(trace_file)) fail("memory trace does not exist: ", trace_file)
  trace_size <- file.info(trace_file)$size
  if(is.na(trace_size)) fail("memory trace cannot be read: ", trace_file)
  if(trace_size == 0) fail("memory trace is empty: ", trace_file)

  required <- c("run_id", "invocation_id", "epoch_seconds", "module", "event", "exit_status",
                "memory_bytes", "working_set_bytes", "memory_limit_bytes")
  header_line <- readLines(trace_file, n = 1L, warn = FALSE)
  if(length(header_line) != 1L || !nzchar(header_line)) fail("memory trace has no header.")
  header_names <- strsplit(sub("\r$", "", header_line), "\t", fixed = TRUE)[[1L]]
  if(anyDuplicated(header_names)) fail("memory trace contains duplicated column names.")
  missing_cols <- setdiff(required, header_names)
  if(length(missing_cols)) fail("memory trace is missing column(s): ", paste(missing_cols, collapse = ", "))

  d <- tryCatch(fread_strict(trace_file, sep = "\t", select = required,
                             na.strings = c("", "NA"), showProgress = FALSE,
                             nThread = 1L),
                error = function(e) fail(conditionMessage(e)))
  d <- d[run_id == target_run]
  if(!nrow(d)){
    render_svg(build_empty_plot(target_run), output_svg)
    return(invisible(NULL))
  }

  if(anyNA(d$run_id) || anyNA(d$invocation_id) || anyNA(d$module) || anyNA(d$event))
    fail("the selected run contains a missing identifier, module, or event.")
  if(any(!nzchar(trimws(as.character(d$invocation_id)))) || any(!nzchar(trimws(as.character(d$module)))))
    fail("the selected run contains a blank invocation identifier or module name.")
  if(any(!d$event %in% c("start", "sample", "end")))
    fail("column 'event' contains a value other than start, sample, or end.")

  d[, source_order := .I]
  d[, epoch_seconds := as_finite_numeric(epoch_seconds, "epoch_seconds")]
  d[, memory_bytes := as_finite_numeric(memory_bytes, "memory_bytes")]
  d[, working_set_bytes := as_finite_numeric(working_set_bytes, "working_set_bytes")]
  d[, memory_limit_bytes := as_finite_numeric(memory_limit_bytes, "memory_limit_bytes", allow_na = TRUE)]
  if(any(d$epoch_seconds < 0) || any(d$memory_bytes < 0) || any(d$working_set_bytes < 0) ||
     any(d$memory_limit_bytes < 0, na.rm = TRUE))
    fail("memory trace contains a negative time or memory measurement.")
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

  module_levels <- unique(as.character(d$module))
  d[, module := factor(module, levels = module_levels)]
  colors <- module_palette(module_levels)

  d[, row_in_invocation := seq_len(.N), by = invocation_id]
  d[, invocation_rows := .N, by = invocation_id]
  stride <- max(1L, ceiling(nrow(d) / MAX_PLOT_POINTS))
  d[, plot_order := seq_len(.N)]
  plot_data <- d[row_in_invocation == 1L | row_in_invocation == invocation_rows |
                   (plot_order - 1L) %% stride == 0L]
  line_data <- plot_data[, if(.N > 1L) .SD else NULL, by = invocation_id]

  starts <- d[, .SD[1L], by = invocation_id]
  data.table::setorder(starts, epoch_seconds, source_order)
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

  # Avoid ggplot's one-observation line warning on the first live update.
  if(nrow(line_data))
    p <- p + ggplot2::geom_line(data = line_data,
                                ggplot2::aes(y = total_gib, group = invocation_id),
                                color = TOTAL_MEMORY_COLOR, linewidth = 0.55)

  p <- p +
    ggplot2::geom_point(ggplot2::aes(y = working_gib, color = module),
                        size = 1.25, alpha = 0.82) +
    ggplot2::geom_point(data = peak, ggplot2::aes(y = working_gib, color = module),
                        shape = 21, fill = "white", size = 3.4, stroke = 0.8,
                        show.legend = FALSE) +
    ggplot2::geom_text(data = peak, ggplot2::aes(y = working_gib, label = label,
                                                 hjust = label_hjust),
                       vjust = 1.5, size = 3.2, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
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
  render_svg(p, output_svg)
  invisible(NULL)
}

tryCatch(main(), error = function(e) fail(conditionMessage(e)))

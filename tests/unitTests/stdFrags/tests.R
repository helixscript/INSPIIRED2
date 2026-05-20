library(data.table)

#' Standardize Genomic Jitter using Gaussian Weighting
#'
#' @param df A data.table with 'seqnames', 'start', 'end', 'strand', and 'reads'.
#' @param side Character, "left" (start) or "right" (end) to determine which edge to fix.
#' @param window The search fence (NT) for pulling in minor peaks. 
#' @param local_radius Radius (NT) for identifying a local maximum (anchor). 
#' @param sd_shrink Factor to divide window by for Gaussian sigma. 
#'
standardize_positions <- function(df, side = "left", window = 10, local_radius = 2, sd_shrink = 4) {
  if (nrow(df) == 0) return(df)
  
  # Ensure data.table and work on a copy to prevent side effects 
  dt <- as.data.table(copy(df))
  sigma <- window / sd_shrink
  
  # 1. Target the active coordinate 
  if (side == "left") {
    dt[, coord := start]
  } else {
    dt[, coord := end]
  }
  
  # 2. Aggregate reads at each unique coordinate 
  counts <- dt[, .(reads = sum(reads)), by = .(seqnames, strand, coord)]
  setorder(counts, seqnames, strand, coord)
  
  # 3. Identify Anchors (Local Maxima) 
  counts[, is_anchor := TRUE]
  if (local_radius > 0) {
    for (i in 1:local_radius) {
      counts[, is_anchor := is_anchor & 
               (reads >= data.table::shift(reads, n = i, fill = 0, type = "lag")) & 
               (reads >= data.table::shift(reads, n = i, fill = 0, type = "lead")), 
             by = .(seqnames, strand)]
    }
  }
  
  anchors <- counts[is_anchor == TRUE]
  anchors[, anchor_pos := coord] # Preserve anchor pos for joining 
  
  # 4. Competitive Mapping via Non-Equi Join 
  unique_coords <- counts[, .(seqnames, strand, coord)]
  unique_coords[, `:=`(win_min = coord - window, win_max = coord + window)] # Pre-calc bounds 
  
  mapping <- anchors[unique_coords, 
                     on = .(seqnames, strand, 
                            coord >= win_min, 
                            coord <= win_max), 
                     allow.cartesian = TRUE]
  
  setnames(mapping, "i.coord", "orig_pos")
  
  # Calculate Gaussian Gravitational Pull 
  mapping[, pull := reads * exp(-((anchor_pos - orig_pos)^2) / (2 * sigma^2))]
  
  best_mapping <- mapping[, .(corrected_coord = anchor_pos[which.max(pull)]), 
                          by = .(seqnames, strand, orig_pos)]
  
  # 5. Update and clean up 
  res <- merge(dt, best_mapping, 
               by.x = c("seqnames", "strand", "coord"), 
               by.y = c("seqnames", "strand", "orig_pos"), 
               all.x = TRUE)
  
  # fcoalesce ensures we keep original positions if no anchor is nearby
  # numeric coercion prevents type mismatch errors 
  if (side == "left") {
    res[, start := data.table::fcoalesce(as.numeric(corrected_coord), as.numeric(start))]
  } else {
    res[, end := data.table::fcoalesce(as.numeric(corrected_coord), as.numeric(end))]
  }
  
  res[, c("coord", "corrected_coord") := NULL] # Efficient in-place removal 
  return(res[])
}




# Simple case.
d1 <- data.table(seqnames = 'chrX', strand = '+',
                 start = c(50, 50, 50, 50, 50, 48, 52, 47, 53, 55),
                 end = c(100, 101, 98, 102, 100, 101, 99, 101, 99, 90),
                 reads = c(10, 20, 10, 10, 12, 5, 5, 3, 2, 1))

d1
standardize_positions(d1, side = "left")




# Two closely spaced peaks.
d2 <- data.table(seqnames = 'chrX', strand = '+',
                 start = c(48, 49, 50, 51, 52, 51, 52, 53, 54, 55, 60, 62),
                 end = c(100, 101, 98, 102, 100, 101, 99, 101, 99, 90, 100, 100),
                 reads = c(1, 3, 10, 3, 1, 1, 3, 10, 3, 1, 5, 2))
              
d2                 
standardize_positions(d2, side = "left")



d2_right<- data.table(seqnames = 'chrX', strand = '+',
                 start = c(50,   52,  48,  51,  52,  51,  52,  53,  54 , 55, 60),
                 end   = c(100, 100, 100, 100, 100, 105, 106, 104, 106, 105, 106),
                 reads = c(1,     3,  10,   8,   1,   3,   1,  3,   10,   3,   2))

d2_right 
standardize_positions(d2_right , side = "right")







# Peak rising out of a tall peaks tail.
d3 <- data.table(seqnames = 'chrX', strand = '+',
                 start = c(48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 68),
                 end = c(100, 101, 98, 102, 100, 101, 99, 101, 99, 90, 100),
                 reads = c(5, 20, 50, 20, 10, 8, 15, 10, 5, 3, 1))
 
d3
standardize_positions(d3, side = "left")


# Tall peak with varability in its tail
d4 <- data.table(seqnames = 'chrX',  strand = '+',
                 start = c(50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60),
                 end = c(100, 101, 98, 102, 100, 101, 99, 101, 99, 90, 100),
                 reads = c(1, 10, 25, 20, 40, 50, 40, 20, 25, 10, 1))

d4
standardize_positions(d4, side = "left")






# Tall peak with varability in its tail
d5 <- data.table(seqnames = 'chrX',  strand = '+',
                 start = c(50,   51, 52,  53,  54,  55, 56,  57, 58, 59,  60, 61,  62, 63, 64,  65),
                 end   = c(100, 101, 98, 102, 100, 101, 99, 101, 99, 90, 100, 99, 100, 99, 100, 99),
                 reads = c(1,     3,  1,   3,   1,   3,  1,   3,  1,  3,   1,  2,   3,  5,   2,  1))

d5
standardize_positions(d5, side = "left")









#------------------------


library(data.table)

# --- SETUP: Hybrid Signal Parameters ---
n_true_peaks   <- 2000      # Strong, high-read biological sites
n_chatter_sites <- 3000     # Closely spaced, low-read noise/chatter regions
reads_per_peak  <- 50       # High read count for strong peaks
reads_per_chat  <- 5        # Low read count for chatter
jitter_sd       <- 2.0      # Narrow biological jitter
test_window     <- 8        # Moderate search fence
test_radius     <- 2        # Standard anchor radius

message("Generating Hybrid Signal: Peaks + Chatter...")

# 1. Create Strong "Golden" Peaks
peaks <- data.table(
  seqnames = sample(paste0("chr", 1:22), n_true_peaks, replace = TRUE),
  true_pos = sample(1:1e6, n_true_peaks),
  type = "Peak",
  strand = sample(c("+", "-"), n_true_peaks, replace = TRUE)
)

# 2. Create "Chatter" Regions (clusters of 3-5 nearby small peaks)
chatter_centers <- sample(1:1e6, n_chatter_sites)
chatter <- data.table(
  seqnames = sample(paste0("chr", 1:22), n_chatter_sites, replace = TRUE),
  true_pos = unlist(lapply(chatter_centers, function(x) x + c(-4, -2, 0, 2, 4))),
  type = "Chatter",
  strand = sample(c("+", "-"), n_chatter_sites, replace = TRUE)
)

true_sites <- rbind(peaks, chatter)
true_sites[, true_start := true_pos][, true_end := true_pos + 150]

# 3. Apply Gaussian Jitter
stress_data <- true_sites[rep(1:.N, times = ifelse(type == "Peak", reads_per_peak, reads_per_chat))]
stress_data[, `:=`(
  start = round(true_start + rnorm(.N, 0, jitter_sd)),
  end = round(true_end + rnorm(.N, 0, jitter_sd)),
  reads = 1
)]

# --- RUN: Pipeline ---
start_time <- Sys.time()
results <- stress_data |>
  standardize_positions(side = "left", window = test_window, local_radius = test_radius) |>
  standardize_positions(side = "right", window = test_window, local_radius = test_radius)
total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

# --- EVALUATION ---
# Calculate recovery by type
results[, is_near := abs(start - true_start) <= 1]
stats <- results[, .(Recovery_1bp = mean(is_near) * 100), by = type]
compression <- (1 - (nrow(results[, .(sum(reads)), by = .(seqnames, start, end, strand)]) / nrow(stress_data))) * 100

cat("\n--- HYBRID STRESS TEST SUMMARY ---\n")
cat(sprintf("Total Rows: %d\n", nrow(stress_data)))
cat(sprintf("Execution Time: %.3f seconds\n", total_time))
cat(sprintf("Peak Recovery (Within 1bp): %.2f%%\n", stats[type == "Peak", Recovery_1bp]))
cat(sprintf("Chatter Stability (Within 1bp): %.2f%%\n", stats[type == "Chatter", Recovery_1bp]))
cat(sprintf("Data Compression Ratio: %.2f%%\n", compression))
cat("----------------------------------\n")





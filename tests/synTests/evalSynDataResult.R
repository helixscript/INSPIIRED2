#!/usr/bin/env Rscript

for (p in c('argparse', 'dplyr', 'tidyr', 'readr', 'GenomicRanges', 'rtracklayer', 'stringdist', 'stringr', 'ggplot2')) suppressPackageStartupMessages(library(p, character.only = TRUE))
options(warn=-1)

# Parse command line arguments.
parser <- ArgumentParser()
parser$add_argument("-s", "--sitesFile",    type="character",  default='sites.rds',            help="Path to sites output file (rds).", metavar="")
parser$add_argument("-m", "--multiHitFile", type="character",  default='multiHitClusters.rds', help="Path to multiHitCluster output file (rds).", metavar="")
parser$add_argument("-f", "--filterFile",   type="character",  default='anchorReadClustAttrition.rds', help="Path to INSPIIRED2 anchor read cluster attrition file (rds).", metavar="")
parser$add_argument("-t", "--truthFile",    type="character",  default='truth.tsv',            help="Path to synthetic data truth file (tsv).", metavar="")
parser$add_argument("-b", "--twoBitPath",   type="character",  default='../../data/referenceGenomes/hg38.2bit', help="Path to local 2bit reference genome.", metavar="")
parser$add_argument("-o", "--outputDir",    type="character",  default='out',                  help="Path to output directory.", metavar="")
parser$add_argument("-w", "--siteWidth",    type="integer",    default=5,                      help="Number of NTs to expand truth positions during evaluation.", metavar="")
parser$add_argument("-v", "--version",      type="integer",    default=2,                      help="INSPIIRED pipeline version (1 or 2). Default: 2", metavar="")
args <- parser$parse_args()

# Test for input errors.
if(! file.exists(args$sitesFile))    { message('Error - the sites file could not be found.'); q(save = "no", status = 1, runLast = FALSE) }
if(! file.exists(args$truthFile))    { message('Error - the truth file could not be found.'); q(save = "no", status = 1, runLast = FALSE) }
if(! dir.exists(args$outputDir))     { dir.create(args$outputDir, recursive = TRUE) }

# ==========================================
# 1. LOAD & PARSE DATA
# ==========================================
message(paste0("Loading and parsing data for INSPIIRED Version ", args$version, "..."))

# Helper to extract V2 coordinates
parse_posid <- function(df, posid_col = "posid") {
  df %>%
    dplyr::mutate(
      chr = stringr::str_extract(!!rlang::sym(posid_col), "^[^\\+\\-]+"),
      strand = stringr::str_extract(!!rlang::sym(posid_col), "[\\+\\-]"),
      pos = as.numeric(stringr::str_extract(!!rlang::sym(posid_col), "\\d+$"))
    )
}

# Truth Data
t <- readr::read_tsv(args$truthFile, show_col_types = FALSE) %>%
  parse_posid() %>%
  dplyr::mutate(
    composite_seqname = paste(trial, subject, sample, chr, strand, sep = "___"),
    t_idx = dplyr::row_number() 
  )

# Sites Data
s_raw <- readRDS(args$sitesFile)
if(inherits(s_raw, "GRanges")) s_raw <- as.data.frame(s_raw)

if(args$version == 1) {
  s <- s_raw %>%
    dplyr::mutate(
      chr = as.character(seqnames),
      pos = start,
      strand = as.character(strand),
      posid = paste0(chr, strand, pos),
      reads = if("reads" %in% names(.)) reads else 1,
      sonicLengths = if("estAbund" %in% names(.)) estAbund else 1,
      repLeaderSeq = NA_character_
    )
  if(!("trial" %in% names(s))) s$trial <- "test"
  if(!("subject" %in% names(s)) & "sampleName" %in% names(s)) s$subject <- s$sampleName
  if(!("sample" %in% names(s)) & "sampleName" %in% names(s)) s$sample <- s$sampleName
} else {
  s <- s_raw %>%
    dplyr::mutate(posid = sub('\\.\\d+$', '', posid)) %>%
    parse_posid()
}

s <- s %>% dplyr::mutate(
  composite_seqname = paste(trial, subject, sample, chr, strand, sep = "___"),
  s_idx = dplyr::row_number()
)

# MultiHits Data
if(file.exists(args$multiHitFile)) {
  m_raw <- readRDS(args$multiHitFile)
  
  if(args$version == 1) {
    if(inherits(m_raw, "GRangesList")) m_raw <- as.data.frame(unlist(m_raw))
    if(inherits(m_raw, "GRanges")) m_raw <- as.data.frame(m_raw)
    
    if(nrow(m_raw) > 0) {
      m <- m_raw %>%
        dplyr::mutate(chr = as.character(seqnames), pos = start, strand = as.character(strand))
      if(!("trial" %in% names(m))) m$trial <- "test"
      if(!("subject" %in% names(m)) & "sampleName" %in% names(m)) m$subject <- m$sampleName
      if(!("sample" %in% names(m)) & "sampleName" %in% names(m)) m$sample <- m$sampleName
      m <- m %>% dplyr::mutate(composite_seqname = paste(trial, subject, sample, chr, strand, sep = "___"))
    } else { m <- data.frame(composite_seqname = character(), pos = numeric()) }
    
  } else {
    if(nrow(m_raw) > 0) {
      m <- m_raw %>%
        dplyr::select(trial, subject, sample, posids) %>%
        tidyr::unnest(posids) %>%
        dplyr::mutate(posid = sub('\\.\\d+$', '', posids)) %>%
        parse_posid() %>%
        dplyr::mutate(composite_seqname = paste(trial, subject, sample, chr, strand, sep = "___"))
    } else { m <- data.frame(composite_seqname = character(), pos = numeric()) }
  }
} else {
  message("Warning - Multi-hit file not found.")
  m <- data.frame(composite_seqname = character(), pos = numeric())
}

# Filter Data (Attrition) - Only applicable to V2
if(args$version == 2 && file.exists(args$filterFile)) {
  f_raw <- readRDS(args$filterFile)
  if(nrow(f_raw) > 0) {
    f <- f_raw %>%
      dplyr::filter(remove == TRUE) %>%
      dplyr::mutate(posid = sub('\\.\\d+$', '', posid)) %>%
      parse_posid() %>%
      dplyr::mutate(composite_seqname = paste(trial, subject, sample, chr, strand, sep = "___"))
  } else { f <- data.frame(composite_seqname = character(), pos = numeric()) }
} else {
  if(args$version == 1) {
    message("INSPIIRED1 selected: Skipping anchor read attrition filter (not applicable).")
  } else {
    message("Warning - Filter file not found. Skipping attrition evaluation.")
  }
  f <- data.frame(composite_seqname = character(), pos = numeric())
}

# ==========================================
# 2. CALCULATE SPACING STATISTICS
# ==========================================
message("Calculating spatial statistics...")
spacing_stats <- t %>%
  dplyr::group_by(trial, subject, sample, chr) %>%
  dplyr::arrange(pos) %>%
  dplyr::mutate(dist_to_next = abs(pos - dplyr::lead(pos))) %>%
  dplyr::filter(!is.na(dist_to_next)) %>%
  dplyr::summarise(
    n_sites = dplyr::n() + 1,
    min_dist = min(dist_to_next),
    max_dist = max(dist_to_next),
    mean_dist = round(mean(dist_to_next), 1),
    median_dist = median(dist_to_next),
    .groups = "drop"
  )

# ==========================================
# 3. VECTORIZED OVERLAP SEARCH
# ==========================================
message("Executing vectorized overlap search...")

g_truth <- GRanges(seqnames = t$composite_seqname, ranges = IRanges(start = t$pos, end = t$pos)) + args$siteWidth
g_sites <- GRanges(seqnames = s$composite_seqname, ranges = IRanges(start = s$pos, end = s$pos)) + args$siteWidth

if(nrow(m) > 0) { g_multi <- GRanges(seqnames = m$composite_seqname, ranges = IRanges(start = m$pos, end = m$pos)) + args$siteWidth } else { g_multi <- GRanges() }
if(nrow(f) > 0) { g_filt <- GRanges(seqnames = f$composite_seqname, ranges = IRanges(start = f$pos, end = f$pos)) + args$siteWidth } else { g_filt <- GRanges() }

hits_sites <- findOverlaps(g_truth, g_sites)
hits_multi <- findOverlaps(g_truth, g_multi)
hits_filt <- findOverlaps(g_truth, g_filt)

# ==========================================
# 4. RESOLVE HITS AND SPLIT SITES
# ==========================================
message("Resolving multi-matches and calculating metrics...")

hit_df <- dplyr::tibble(
  t_idx = as.integer(queryHits(hits_sites)),
  s_idx = as.integer(subjectHits(hits_sites))
)

hit_calc <- hit_df %>%
  dplyr::left_join(dplyr::select(t, t_idx, t_pos = pos, t_nReads = nReads, t_nFrags = nFrags, t_leader = leaderSeq), by = "t_idx") %>%
  dplyr::left_join(dplyr::select(s, s_idx, s_pos = pos, s_reads = reads, s_frags = sonicLengths, s_leader = repLeaderSeq), by = "s_idx") %>%
  dplyr::mutate(
    absDiff = abs(t_pos - s_pos),
    leaderSeqDist = if(args$version == 2) stringdist::stringdist(t_leader, s_leader, method = "osa") else NA_real_
  )

hit_resolved <- hit_calc %>%
  dplyr::group_by(t_idx) %>%
  dplyr::arrange(absDiff, dplyr::desc(s_reads)) %>%
  dplyr::mutate(
    splitSite = dplyr::n() > 1,
    posDiff = t_pos - s_pos,
    readDiff = t_nReads - s_reads,
    fragDiff = t_nFrags - s_frags
  ) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(t_idx, posDiff, readDiff, fragDiff, leaderSeqDist, splitSite)

# ==========================================
# 5. ASSEMBLE FINAL TABLES
# ==========================================
message("Formatting final outputs...")

found_in_multi <- unique(as.integer(queryHits(hits_multi)))
found_in_filt <- unique(as.integer(queryHits(hits_filt)))

tab <- t %>%
  dplyr::left_join(hit_resolved, by = "t_idx") %>%
  dplyr::mutate(
    splitSite = tidyr::replace_na(splitSite, FALSE),
    found_in_m = t_idx %in% found_in_multi,
    found_in_f = t_idx %in% found_in_filt,
    found = !is.na(posDiff) | found_in_m | found_in_f
  )

sumTab <- dplyr::tibble(
  nSitesExpected = dplyr::n_distinct(t$posid), 
  nSitesRecovered = dplyr::n_distinct(s$posid), 
  percentUniqueRecovery = sprintf("%.2f%%", (sum(!is.na(tab$posDiff)) / dplyr::n_distinct(t$posid)) * 100),
  percentTotalRecovery = sprintf("%.2f%%",  (sum(tab$found == TRUE) / dplyr::n_distinct(t$posid)) * 100),
  percentFiltered = if(args$version == 2 && file.exists(args$filterFile)) sprintf("%.2f%%", (sum(tab$found_in_f == TRUE) / dplyr::n_distinct(t$posid)) * 100) else "NA",
  leaderSeqDistMean = sprintf("%.2f", mean(tab$leaderSeqDist, na.rm = TRUE)),
  leaderSeqDistSD = sprintf("%.2f", sd(tab$leaderSeqDist, na.rm = TRUE))
)

# Clean up internal boolean columns before export
tab <- tab %>%
  dplyr::select(-chr, -strand, -pos, -composite_seqname, -t_idx, -found_in_m, -found_in_f)

readr::write_tsv(sumTab, file.path(args$outputDir, 'table1_summary.tsv'))
readr::write_tsv(tab, file.path(args$outputDir, 'table2_site_details.tsv'))
readr::write_tsv(spacing_stats, file.path(args$outputDir, 'table3_spacing_stats.tsv'))

# ==========================================
# 6. GENERATE GENOMIC MAP
# ==========================================
if(file.exists(args$twoBitPath)) {
  message("Extracting chromosome lengths from ", args$twoBitPath, "...")
  seq_info <- seqinfo(TwoBitFile(args$twoBitPath))
  
  lengths_df <- dplyr::tibble(
    chr = seqnames(seq_info),
    length = seqlengths(seq_info)
  ) %>%
    dplyr::filter(chr %in% paste0("chr", c(1:22, "X", "Y"))) %>%
    dplyr::mutate(chr_factor = factor(chr, levels = paste0("chr", c(1:22, "X", "Y"))))
  
  message("Generating genomic map visualization...")
  
  plot_data <- tab %>%
    dplyr::mutate(
      chr = stringr::str_extract(posid, "^[^\\+\\-]+"),
      pos = as.numeric(stringr::str_extract(posid, "\\d+$")),
      chr_factor = factor(chr, levels = levels(lengths_df$chr_factor)),
      Status = dplyr::case_when(
        !is.na(posDiff) ~ "1. Unique Recovery",
        is.na(posDiff) & found == TRUE ~ "2. Multi-hit / Filtered",
        found == FALSE ~ "3. Lost (Unmapped)"
      )
    ) %>%
    dplyr::filter(!is.na(chr_factor))
  
  p <- ggplot() +
    geom_segment(data = lengths_df, 
                 aes(x = 0, xend = length / 1e6, y = chr_factor, yend = chr_factor),
                 color = "gray80", linewidth = 3, lineend = "round") +
    geom_jitter(data = plot_data, 
                aes(x = pos / 1e6, y = chr_factor, color = Status),
                height = 0.3, size = 1.2, alpha = 0.7) +
    scale_color_manual(values = c(
      "1. Unique Recovery" = "#2c7bb6",
      "2. Multi-hit / Filtered" = "#fdae61",
      "3. Lost (Unmapped)" = "#d7191c"
    )) +
    labs(
      title = "Genomic Distribution of Simulated Integration Sites",
      subtitle = paste0("Synthetic Data Evaluation (INSPIIRED Version ", args$version, ")"),
      x = "Genomic Position (Mb)",
      y = "",
      color = "Site Status"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.y = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    ) +
    scale_y_discrete(limits = rev(levels(lengths_df$chr_factor)))
  
  plot_path <- file.path(args$outputDir, 'figure1_genome_map.png')
  ggsave(plot_path, plot = p, width = 12, height = 8, dpi = 300, bg = "white")
  message("Map saved to ", plot_path)
} else {
  message("Warning: 2bit file not found at ", args$twoBitPath, ". Skipping genomic map generation.")
}

message("Done! Results saved to ", args$outputDir)
q(save = 'no', status = 0, runLast = FALSE)
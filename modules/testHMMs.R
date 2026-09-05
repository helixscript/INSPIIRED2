#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'ShortRead', 'parallel', 'data.table', 'BiocParallel', 'stringi')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir", type = "character", required = TRUE, help = "Directory for output files")
parser$add_argument("--inputData", type = "character", required = TRUE, help = "Path to demultiplex module output")
parser$add_argument("--softwareRoot", type = "character", required = TRUE, help = "Path to INSPIIRED2 installation.")
parser$add_argument("--threads", type = "integer", default = 50, help = "Number of threads to use.")
parser$add_argument("--fileTag", type = "character", default = "testHMMs", help = "String appended to output files in the output directory.")
parser$add_argument("--ramDiskPath", type = "character", default = "/dev/shm", help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
parser$add_argument("--maxReadStartPos", type = "integer", default = 50, help = "Max. read position to show on plot before overflow bin.")
parser$add_argument("--startPosBinWidth", type = "integer", default = 3, help = "Bin size of alignment start positions before overflow bin.")
parser$add_argument("--scoreBinWidth", type = "double", default = 3, help = "Bin size for HMM scores.")
parser$add_argument("--minScoreBinPct", type = "double", default = 1, help = "Min. percent of total reads falling into a grid square needed to print square.")
parser$add_argument("--facetCols", type = "integer", default = 4, help = "Number of facet columns in the output plot.")
parser$add_argument("--disableHorizontalGuides", action = "store_true",  default = FALSE,          help = "Disable horizontal grid lines.")
parser$add_argument("--horizontalGuideEvery", type = "double", default = 6, help = "Number of HMM score points between horizontal grid lines.")
parser$add_argument("--HMMparams", type = "character", required = FALSE, default = 'none', help = "HMM parameter string.")

runModule <- function() {
  startModule()
  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  on.exit({
    unlink(args$tmpDir, recursive = TRUE, force = TRUE)
    unlink(args$logDir, recursive = TRUE, force = TRUE)
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  
  updateLog('Starting testHMMs module.')
  if (!file.exists(args$inputData)) stop(paste0('Error - input file "', args$inputData, '" does not exist.'))
  
  o <- readRDS(args$inputData)
  if (nrow(o) == 0) stop(paste0('Error - input file "', args$inputData, '" has zero rows of data.'))
  
  buildHMMParameterTable <- function(hmmNames) {
    parseOverrides <- function(x) {
      if (is.null(x) || length(x) != 1L || !nzchar(trimws(x)) || tolower(trimws(x)) == "none") return(list())
      
      a <- strsplit(x, "|", fixed = TRUE)[[1]]
      
      z <- lapply(a, function(x) {
        p <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
        if (length(p) != 8L) stop("Error - --HMMparams entry must contain an HMM name followed by 7 parameters: ", x)
        p
      })
      
      n <- vapply(z, `[`, character(1), 1L)
      if (anyDuplicated(n)) stop("Error - the same HMM was defined more than once in --HMMparams.")
      
      z <- lapply(z, function(x) x[-1L])
      names(z) <- n
      z
    }
    
    readHMMLength <- function(hmmName) {
      f <- file.path(args$softwareRoot, "data", "hmms", hmmName)
      if (!file.exists(f)) stop("Error - HMM file does not exist: ", f)
      
      h <- readLines(f, warn = FALSE)
      x <- grep("^LENG\\s+", h, value = TRUE)
      
      if (length(x) != 1L) stop("Error - could not uniquely determine LENG for hmm: ", hmmName)
      
      n <- suppressWarnings(as.integer(strsplit(trimws(x), "\\s+")[[1]][2]))
      
      if (is.na(n)) stop("Error - could not parse HMM length for hmm: ", hmmName)
      n
    }
    
    makeRow <- function(hmmName, p, source) {
      if (length(p) != 7L) stop("Error - expected 7 HMM parameters for hmm: ", hmmName)
      if (!grepl("^(TRUE|FALSE)$", p[5], ignore.case = TRUE)) stop("Error - HMMmatchEnd must be TRUE or FALSE for hmm: ", hmmName)
      
      z <- data.table(
        leaderSeqHMM = hmmName,
        HMMminStartPos = suppressWarnings(as.integer(p[1])),
        HMMmaxStartPos = suppressWarnings(as.integer(p[2])),
        HMMminFullBitScore = suppressWarnings(as.numeric(p[3])),
        HMMmaxFullBitScore = suppressWarnings(as.numeric(p[4])),
        HMMmatchEnd = grepl("^TRUE$", p[5], ignore.case = TRUE),
        HMMmatchTerminalSeq = as.character(p[6]),
        HMMmatchEndRadius = suppressWarnings(as.integer(p[7])),
        hmmLength = readHMMLength(hmmName),
        parameterSource = source
      )
      
      if (any(is.na(z[, .(HMMminStartPos, HMMmaxStartPos, HMMminFullBitScore, HMMmaxFullBitScore, HMMmatchEndRadius, hmmLength)])))
        stop("Error - one or more HMM parameters could not be parsed for hmm: ", hmmName)
      
      if (z$HMMminStartPos < 1L || z$HMMmaxStartPos < z$HMMminStartPos)
        stop("Error - invalid HMM start-position range for hmm: ", hmmName)
      
      if (z$HMMmaxFullBitScore < z$HMMminFullBitScore)
        stop("Error - invalid HMM score range for hmm: ", hmmName)
      
      if (z$HMMmatchEndRadius < 0L)
        stop("Error - HMMmatchEndRadius must be >= 0 for hmm: ", hmmName)
      
      if (!nzchar(z$HMMmatchTerminalSeq))
        stop("Error - HMMmatchTerminalSeq was empty for hmm: ", hmmName)
      
      z
    }
    
    overrides <- parseOverrides(args$HMMparams)
    
    expected <- c(
      "HMMminStartPos",
      "HMMmaxStartPos",
      "HMMminFullBitScore",
      "HMMmaxFullBitScore",
      "HMMmatchEnd",
      "HMMmatchTerminalSeq",
      "HMMmatchEndRadius"
    )
    
    rbindlist(lapply(unique(as.character(hmmNames)), function(hmmName) {
      if (hmmName %in% names(overrides))
        return(makeRow(hmmName, overrides[[hmmName]], "--HMMparams"))
      
      cfgFile <- file.path(args$softwareRoot, "data", "hmms", sub("\\.hmm$", ".cfg", hmmName))
      
      if (!file.exists(cfgFile))
        stop("Error - could not determine processing parameters for hmm: ", hmmName)
      
      p <- readr::read_tsv(
        cfgFile,
        col_names = FALSE,
        col_types = readr::cols(.default = readr::col_character()),
        show_col_types = FALSE,
        progress = FALSE
      )
      
      if (nrow(p) != 7L || ncol(p) != 2L)
        stop("Error - the hmm cfg file for hmm: ", hmmName, " did not have the expected dimensions.")
      
      if (any(is.na(p)))
        stop("Error - the hmm cfg file for hmm: ", hmmName, " contained one or more NA values.")
      
      names(p) <- c("name", "value")
      
      if (anyDuplicated(p$name) || !setequal(p$name, expected))
        stop("Error - the hmm cfg file for hmm: ", hmmName, " did not contain the expected parameter names.")
      
      makeRow(hmmName, p$value[match(expected, p$name)], "cfg")
    }), use.names = TRUE, fill = TRUE)
  }
  
  testHMM <- function(x, hmmParameters) {
    emptyResult <- function() data.table(
      trial = character(),
      subject = character(),
      sample = character(),
      leaderSeqHMM = character(),
      targetStart = integer(),
      targetEnd = integer(),
      fullScore = numeric()
    )
    
    hmmName <- as.character(x$leaderSeqHMM[1L])
    hp <- hmmParameters[leaderSeqHMM == hmmName]
    
    if (nrow(hp) != 1L) stop("Error - could not uniquely resolve HMM parameters for: ", hmmName)
    
    prefix <- tempfile(pattern = "nhmmer_", tmpdir = args$ramDisk)
    fastaFile <- paste0(prefix, ".fa")
    tblFile <- paste0(prefix, ".tbl")
    searchFile <- paste0(prefix, ".hmmSearch")
    
    on.exit(unlink(c(fastaFile, tblFile, searchFile), force = TRUE), add = TRUE)
    
    s <- DNAStringSet(x$anchorReadSeq)
    names(s) <- paste0("s", seq_along(s))
    writeXStringSet(s, fastaFile)
    
    hmmFile <- file.path(args$softwareRoot, "data", "hmms", hmmName)
    
    status <- system2(
      "nhmmer",
      args = c(
        "--dna", "--cpu", "1",
        "--F1", "1", "--F2", "1", "--F3", "1",
        "-T", "-5", "--incT", "-5", "--nobias",
        "--popen", "0.15", "--pextend", "0.05",
        "--tblout", shQuote(tblFile),
        shQuote(hmmFile),
        shQuote(fastaFile)
      ),
      stdout = searchFile,
      stderr = searchFile
    )
    
    if (status != 0L) stop("nhmmer failed with exit status ", status)
    
    lines <- readLines(tblFile, warn = FALSE)
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    
    if (!length(lines)) return(emptyResult())
    
    hits <- fread(
      text = paste(lines, collapse = "\n"),
      header = FALSE,
      select = c(1L, 6L, 7L, 8L, 12L, 14L),
      fill = TRUE,
      nThread = 1L,
      showProgress = FALSE
    )
    
    setnames(
      hits,
      c(
        "targetName",
        "hmmEnd",
        "targetStart",
        "targetEnd",
        "strand",
        "fullScore"
      )
    )
    
    neg <- hits$strand == "-"
    
    if (any(neg)) {
      tmp <- hits$targetStart[neg]
      hits$targetStart[neg] <- hits$targetEnd[neg]
      hits$targetEnd[neg] <- tmp
    }
    
    hits <- hits[hits[, .I[which.max(fullScore)], by = targetName]$V1]
    
    if (isTRUE(hp$HMMmatchEnd[1L])) {
      hits <- hits[
        abs(
          hp$hmmLength[1L] -
            hmmEnd
        ) <= hp$HMMmatchEndRadius[1L]
      ]
    }
    
    if (!nrow(hits)) return(emptyResult())
    
    terminalSeq <- as.character(hp$HMMmatchTerminalSeq[1L])
    radius <- as.integer(hp$HMMmatchEndRadius[1L])
    
    if (!grepl("none", terminalSeq, ignore.case = TRUE)) {
      terminalSeq <- toupper(terminalSeq)
      
      readIndex <- match(hits$targetName, names(s))
      
      if (anyNA(readIndex))
        stop("Error - could not map one or more nhmmer hits back to input reads for hmm: ", hmmName)
      
      anchorReadSeq <- toupper(as.character(x$anchorReadSeq[readIndex]))
      
      terminalMatchSeq <- substr(
        anchorReadSeq,
        hits$targetEnd - (nchar(terminalSeq) - 1L) - radius,
        hits$targetEnd + radius
      )
      
      ends <- stringr::str_locate(
        terminalMatchSeq,
        terminalSeq
      )[, 2]
      
      keep <- !is.na(ends)
      
      hits <- hits[keep]
      ends <- ends[keep]
      
      if (!nrow(hits)) return(emptyResult())
      
      hits[
        ,
        targetEnd :=
          targetEnd -
          (nchar(terminalSeq) + radius) +
          ends
      ]
    }
    
    hits[
      ,
      `:=`(
        trial = as.character(x$trial[1L]),
        subject = as.character(x$subject[1L]),
        sample = as.character(x$sample[1L]),
        leaderSeqHMM = hmmName
      )
    ]
    
    hits[
      ,
      .(
        trial,
        subject,
        sample,
        leaderSeqHMM,
        targetStart,
        targetEnd,
        fullScore
      )
    ]
  }
  
  build_HMM_tests <- function(o, hmmParameters) {
    hmmInput <- o[, .(trial, subject, sample, leaderSeqHMM, anchorReadSeq)]
    groupCols <- c("trial", "subject", "sample", "leaderSeqHMM")
    
    hmmInput[, (groupCols) := lapply(.SD, as.character), .SDcols = groupCols]
    groupIdx <- hmmInput[, .(idx = list(.I)), by = groupCols]
    
    param <- MulticoreParam(
      workers = min(args$threads, nrow(groupIdx)),
      tasks = nrow(groupIdx),
      stop.on.error = TRUE
    )
    
    rbindlist(
      bplapply(
        groupIdx$idx,
        function(ii) testHMM(hmmInput[ii], hmmParameters),
        BPPARAM = param
      ),
      use.names = TRUE,
      fill = TRUE
    )
  }
  
  plotHMMStartHeatmap <- function(
    hmmResults,
    hmmParameters,
    maxStart = 50L,
    startBinWidth = 3L,
    scoreBinWidth = 3,
    minScoreBinPct = 1,
    facetCols = 4L,
    showHorizontalGuides = TRUE,
    horizontalGuideEvery = 6,
    horizontalGuideLinewidth = 0.25,
    horizontalGuideColour = "grey65",
    hitBorderColour = "grey70",
    hitBorderLinewidth = 0.10,
    parameterBoxColour = "blue3",
    parameterBoxLinewidth = 0.8,
    axisTextSizeX = 8.5,
    axisTextSizeY = 8.5
  ) {
    d <- as.data.table(copy(hmmResults))
    hp <- as.data.table(copy(hmmParameters))
    
    roundPctTo100 <- function(counts) {
      counts <- as.numeric(counts)
      total <- sum(counts)
      
      if (total <= 0) return(rep(0L, length(counts)))
      
      rawPct <- 100 * counts / total
      out <- floor(rawPct)
      remainder <- rawPct - out
      missing <- as.integer(round(100 - sum(out)))
      
      if (missing > 0) {
        idx <- order(remainder, decreasing = TRUE)
        out[idx[seq_len(missing)]] <- out[idx[seq_len(missing)]] + 1L
      } else if (missing < 0) {
        idx <- order(remainder, decreasing = FALSE)
        out[idx[seq_len(abs(missing))]] <- out[idx[seq_len(abs(missing))]] - 1L
      }
      
      as.integer(out)
    }
    
    requiredCols <- c("trial", "subject", "sample", "leaderSeqHMM", "targetStart", "fullScore")
    requiredParamCols <- c(
      "leaderSeqHMM",
      "HMMminStartPos",
      "HMMmaxStartPos",
      "HMMminFullBitScore",
      "HMMmaxFullBitScore"
    )
    
    missingCols <- setdiff(requiredCols, names(d))
    missingParamCols <- setdiff(requiredParamCols, names(hp))
    
    if (length(missingCols))
      stop("Missing required columns: ", paste(missingCols, collapse = ", "))
    
    if (length(missingParamCols))
      stop("Missing required HMM parameter columns: ", paste(missingParamCols, collapse = ", "))
    
    maxStart <- as.integer(maxStart)
    startBinWidth <- as.integer(startBinWidth)
    scoreBinWidth <- as.numeric(scoreBinWidth)
    minScoreBinPct <- as.numeric(minScoreBinPct)
    facetCols <- as.integer(facetCols)
    horizontalGuideEvery <- as.numeric(horizontalGuideEvery)
    horizontalGuideLinewidth <- as.numeric(horizontalGuideLinewidth)
    hitBorderLinewidth <- as.numeric(hitBorderLinewidth)
    parameterBoxLinewidth <- as.numeric(parameterBoxLinewidth)
    axisTextSizeX <- as.numeric(axisTextSizeX)
    axisTextSizeY <- as.numeric(axisTextSizeY)
    
    if (maxStart < 1L) stop("maxStart must be >= 1.")
    if (startBinWidth < 1L) stop("startBinWidth must be >= 1.")
    if (!is.finite(scoreBinWidth) || scoreBinWidth <= 0) stop("scoreBinWidth must be greater than zero.")
    if (!is.finite(minScoreBinPct) || minScoreBinPct <= 0 || minScoreBinPct > 100) stop("minScoreBinPct must be > 0 and <= 100.")
    if (facetCols < 1L) stop("facetCols must be >= 1.")
    if (!is.finite(horizontalGuideEvery) || horizontalGuideEvery <= 0) stop("horizontalGuideEvery must be greater than zero.")
    if (!is.finite(hitBorderLinewidth) || hitBorderLinewidth < 0) stop("hitBorderLinewidth must be >= 0.")
    if (!is.finite(parameterBoxLinewidth) || parameterBoxLinewidth < 0) stop("parameterBoxLinewidth must be >= 0.")
    
    d[, `:=`(
      trial = as.character(trial),
      subject = as.character(subject),
      sample = as.character(sample),
      leaderSeqHMM = as.character(leaderSeqHMM),
      targetStart = as.integer(targetStart),
      fullScore = as.numeric(fullScore)
    )]
    
    hp[, `:=`(
      leaderSeqHMM = as.character(leaderSeqHMM),
      HMMminStartPos = as.integer(HMMminStartPos),
      HMMmaxStartPos = as.integer(HMMmaxStartPos),
      HMMminFullBitScore = as.numeric(HMMminFullBitScore),
      HMMmaxFullBitScore = as.numeric(HMMmaxFullBitScore)
    )]
    
    d <- d[
      !is.na(targetStart) &
        targetStart >= 1L &
        !is.na(fullScore) &
        is.finite(fullScore)
    ]
    
    if (!nrow(d)) stop("No valid HMM results remain.")
    if (anyDuplicated(hp$leaderSeqHMM)) stop("HMM parameter table contains duplicate HMM names.")
    
    missingHMMs <- setdiff(unique(d$leaderSeqHMM), hp$leaderSeqHMM)
    
    if (length(missingHMMs))
      stop("No HMM processing parameters were found for: ", paste(missingHMMs, collapse = ", "))
    
    hp <- hp[leaderSeqHMM %in% unique(d$leaderSeqHMM)]
    
    groupCols <- c("trial", "subject", "sample", "leaderSeqHMM")
    
    sampleOrderDT <- unique(
      d[
        ,
        .(
          trial,
          subject,
          sample,
          leaderSeqHMM
        )
      ]
    )
    
    sampleOrderDT[, sampleOrder := .I]
    
    rawMaxScore <- max(d$fullScore, na.rm = TRUE)
    
    d[, scoreBin := floor(fullScore / scoreBinWidth) * scoreBinWidth]
    
    d[
      ,
      `:=`(
        scoreBinUpper = scoreBin + scoreBinWidth,
        scoreMid = scoreBin + scoreBinWidth / 2
      )
    ]
    
    scoreStats <- d[
      ,
      .N,
      by = c(
        groupCols,
        "scoreBin",
        "scoreBinUpper"
      )
    ]
    
    scoreStats[
      ,
      pctScoreBin :=
        100 *
        N /
        sum(N),
      by = groupCols
    ]
    
    qualifyingScores <- scoreStats[
      scoreBin >= 0 &
        pctScoreBin >= minScoreBinPct
    ]
    
    if (nrow(qualifyingScores)) {
      maxDisplayBin <- max(
        qualifyingScores$scoreBin,
        na.rm = TRUE
      )
    } else {
      nonnegativeBins <- d[
        scoreBin >= 0,
        scoreBin
      ]
      
      if (!length(nonnegativeBins))
        stop("No non-negative HMM scores are available for plotting.")
      
      warning(
        "No non-negative HMM score bin contains >= ",
        minScoreBinPct,
        "% of any sample/HMM group; using the full non-negative score range."
      )
      
      maxDisplayBin <- max(
        nonnegativeBins,
        na.rm = TRUE
      )
    }
    
    parameterMaxBin <- max(
      floor(
        hp$HMMmaxFullBitScore /
          scoreBinWidth
      ) *
        scoreBinWidth,
      na.rm = TRUE
    )
    
    maxDisplayBin <- max(
      maxDisplayBin,
      parameterMaxBin
    )
    
    maxDisplayScore <-
      maxDisplayBin +
      scoreBinWidth
    
    d[
      ,
      startBinLower :=
        (
          (targetStart - 1L) %/%
            startBinWidth
        ) *
        startBinWidth +
        1L
    ]
    
    d[
      ,
      startBinUpper :=
        startBinLower +
        startBinWidth -
        1L
    ]
    
    startBinStarts <- seq(
      1L,
      maxStart,
      by = startBinWidth
    )
    
    startBinEnds <- pmin(
      startBinStarts +
        startBinWidth -
        1L,
      maxStart
    )
    
    startBinLabels <- paste0(
      startBinStarts,
      "-",
      startBinEnds
    )
    
    overflowLabel <- paste0(
      ">",
      maxStart
    )
    
    xLevels <- c(
      startBinLabels,
      overflowLabel
    )
    
    d[, xCategory := overflowLabel]
    
    d[
      targetStart <= maxStart,
      xCategory :=
        paste0(
          startBinLower,
          "-",
          pmin(
            startBinUpper,
            maxStart
          )
        )
    ]
    
    sampleStats <- d[
      ,
      .(
        totalReads = .N,
        withinReads = sum(
          targetStart >= 1L &
            targetStart <= maxStart
        ),
        overflowReads = sum(
          targetStart > maxStart
        )
      ),
      by = groupCols
    ]
    
    sampleStats[
      ,
      c(
        "pctWithin",
        "pctOverflow"
      ) := {
        p <- roundPctTo100(
          c(
            withinReads,
            overflowReads
          )
        )
        
        .(
          p[1],
          p[2]
        )
      },
      by = groupCols
    ]
    
    sampleStats <- merge(
      sampleOrderDT,
      sampleStats,
      by = groupCols,
      all.x = TRUE,
      sort = FALSE
    )
    
    setorder(
      sampleStats,
      sampleOrder
    )
    
    sampleStats[
      ,
      facetID := paste0(
        "HMM: ",
        leaderSeqHMM,
        "\n",
        trial,
        " | ",
        subject,
        " | ",
        sample,
        "\n",
        "Total reads = ",
        format(
          totalReads,
          big.mark = ",",
          scientific = FALSE
        ),
        "\n",
        "positions 1-",
        maxStart,
        " = ",
        pctWithin,
        "% | >",
        maxStart,
        " = ",
        pctOverflow,
        "%"
      )
    ]
    
    facetOrder <- sampleStats$facetID
    
    heatmapAll <- d[
      ,
      .N,
      by = c(
        groupCols,
        "xCategory",
        "scoreBin",
        "scoreBinUpper",
        "scoreMid"
      )
    ]
    
    heatmapAll <- merge(
      heatmapAll,
      sampleStats[
        ,
        .(
          trial,
          subject,
          sample,
          leaderSeqHMM,
          totalReads,
          facetID
        )
      ],
      by = groupCols,
      all.x = TRUE,
      sort = FALSE
    )
    
    heatmapAll[
      ,
      pct :=
        100 *
        N /
        totalReads
    ]
    
    heatmapDT <- heatmapAll[
      scoreBin >= 0 &
        scoreBin <= maxDisplayBin
    ]
    
    heatmapDT[
      ,
      xIndex :=
        match(
          xCategory,
          xLevels
        )
    ]
    
    heatmapDT[
      ,
      facetID :=
        factor(
          facetID,
          levels = facetOrder
        )
    ]
    
    startLeft <- function(p) {
      p <- as.integer(p)
      
      if (p > maxStart)
        return(
          length(startBinLabels) +
            0.5
        )
      
      i <- findInterval(
        p,
        startBinStarts
      )
      
      n <- startBinEnds[i] -
        startBinStarts[i] +
        1
      
      i -
        0.5 +
        (
          p -
            startBinStarts[i]
        ) /
        n
    }
    
    startRight <- function(p) {
      p <- as.integer(p)
      
      if (p > maxStart)
        return(
          length(startBinLabels) +
            1.5
        )
      
      i <- findInterval(
        p,
        startBinStarts
      )
      
      n <- startBinEnds[i] -
        startBinStarts[i] +
        1
      
      i -
        0.5 +
        (
          p -
            startBinStarts[i] +
            1
        ) /
        n
    }
    
    parameterBoxes <- merge(
      sampleStats[
        ,
        .(
          trial,
          subject,
          sample,
          leaderSeqHMM,
          facetID
        )
      ],
      hp[
        ,
        .(
          leaderSeqHMM,
          HMMminStartPos,
          HMMmaxStartPos,
          HMMminFullBitScore,
          HMMmaxFullBitScore,
          HMMmatchEnd,
          HMMmatchTerminalSeq,
          HMMmatchEndRadius,
          parameterSource
        )
      ],
      by = "leaderSeqHMM",
      all.x = TRUE,
      sort = FALSE,
      allow.cartesian = TRUE
    )
    
    parameterBoxes[
      ,
      `:=`(
        xmin = vapply(
          HMMminStartPos,
          startLeft,
          numeric(1)
        ),
        xmax = vapply(
          HMMmaxStartPos,
          startRight,
          numeric(1)
        ),
        ymin = HMMminFullBitScore,
        ymax = HMMmaxFullBitScore,
        facetID = factor(
          facetID,
          levels = facetOrder
        )
      )
    ]
    
    fillScale <- scale_fill_gradientn(
      colours = c(
        "white",
        "#fefffd",
        "#f8fcf4",
        "#eef8e8",
        "#d9f0d3",
        "#31a354",
        "yellow",
        "orange",
        "orangered",
        "red"
      ),
      values = scales::rescale(
        c(
          0,
          0.1,
          1,
          3,
          5,
          10,
          18,
          25,
          40,
          50
        ),
        from = c(
          0,
          50
        )
      ),
      limits = c(
        0,
        50
      ),
      breaks = c(
        1,
        10,
        25,
        50
      ),
      labels = c(
        "1",
        "10",
        "25",
        ">50"
      ),
      oob = scales::squish,
      name = "% reads"
    )
    
    guideScores <- numeric()
    
    if (maxDisplayScore > horizontalGuideEvery) {
      guideScores <- seq(
        horizontalGuideEvery,
        floor(
          (
            maxDisplayScore -
              .Machine$double.eps
          ) /
            horizontalGuideEvery
        ) *
          horizontalGuideEvery,
        by = horizontalGuideEvery
      )
    }
    
    yTopBreak <- floor(
      maxDisplayScore /
        horizontalGuideEvery
    ) *
      horizontalGuideEvery
    
    yBreaks <- if (yTopBreak > 0)
      seq(
        0,
        yTopBreak,
        by = horizontalGuideEvery
      )
    else
      0
    
    nFacets <- nrow(
      sampleStats
    )
    
    nRows <- ceiling(
      nFacets /
        facetCols
    )
    
    p <- ggplot(
      heatmapDT,
      aes(
        x = xIndex,
        y = scoreMid,
        fill = pct
      )
    ) +
      geom_tile(
        width = 1,
        height = scoreBinWidth,
        colour = hitBorderColour,
        linewidth = hitBorderLinewidth
      ) +
      facet_wrap(
        ~facetID,
        ncol = facetCols,
        drop = FALSE
      ) +
      scale_x_continuous(
        breaks = seq_along(xLevels),
        labels = xLevels,
        limits = c(
          0.5,
          length(xLevels) + 0.5
        ),
        expand = c(
          0,
          0
        )
      ) +
      scale_y_continuous(
        breaks = yBreaks,
        limits = c(
          0,
          maxDisplayScore
        ),
        expand = c(
          0,
          0
        )
      ) +
      fillScale +
      theme_bw() +
      theme(
        strip.text = element_text(
          size = 7.2,
          lineheight = 0.95
        ),
        strip.background = element_rect(
          fill = "grey90"
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(
          size = axisTextSizeX,
          angle = 45,
          hjust = 1
        ),
        axis.text.y = element_text(
          size = axisTextSizeY
        ),
        axis.title = element_text(
          size = 10.5
        ),
        plot.title = element_text(
          size = 14
        ),
        plot.subtitle = element_text(
          size = 10
        ),
        legend.position = "right"
      ) +
      labs(
        x = "HMM alignment start position",
        y = "HMM full score",
        title = "HMM alignment start positions vs HMM score heatmap by sample and HMM",
        subtitle = paste0(
          "Read starts are grouped into ",
          startBinWidth,
          "-position bins through position ",
          maxStart,
          "; >",
          maxStart,
          " is the overflow bin; scores are grouped into ",
          scoreBinWidth,
          "-unit bins"
        )
      )
    
    if (
      showHorizontalGuides &&
      length(guideScores)
    ) {
      p <- p +
        geom_hline(
          yintercept = guideScores,
          linewidth = horizontalGuideLinewidth,
          colour = horizontalGuideColour
        )
    }
    
    p <- p +
      geom_vline(
        xintercept =
          length(startBinLabels) +
          0.5,
        linewidth = 0.7
      )
    
    p <- p +
      geom_rect(
        data = parameterBoxes,
        aes(
          xmin = xmin,
          xmax = xmax,
          ymin = ymin,
          ymax = ymax
        ),
        inherit.aes = FALSE,
        fill = NA,
        colour = parameterBoxColour,
        linewidth = parameterBoxLinewidth
      )
    
    list(
      plot = p,
      heatmap = heatmapDT,
      heatmapAllScores = heatmapAll,
      scoreStats = scoreStats,
      sampleStats = sampleStats,
      hmmParameters = hp,
      parameterBoxes = parameterBoxes,
      rawMaxScore = rawMaxScore,
      maxDisplayBin = maxDisplayBin,
      maxDisplayScore = maxDisplayScore,
      startBinWidth = startBinWidth,
      scoreBinWidth = scoreBinWidth,
      minScoreBinPct = minScoreBinPct,
      maxStart = maxStart,
      nFacets = nFacets,
      nRows = nRows,
      facetCols = facetCols,
      xLevels = xLevels
    )
  }
  
  hmmNames <- unique(
    as.character(
      o$leaderSeqHMM
    )
  )
  
  hmmParameters <- buildHMMParameterTable(
    hmmNames
  )
  
  hmmResults <- build_HMM_tests(
    o,
    hmmParameters
  )
  
  if (!nrow(hmmResults))
    stop(
      "No HMM hits remained after applying the configured HMM-end and terminal-sequence requirements."
    )
  
  hmmAnalysis <- plotHMMStartHeatmap(
    hmmResults,
    hmmParameters,
    maxStart = args$maxReadStartPos,
    startBinWidth = args$startPosBinWidth,
    scoreBinWidth = args$scoreBinWidth,
    minScoreBinPct = args$minScoreBinPct,
    facetCols = args$facetCols,
    showHorizontalGuides = args$disableHorizontalGuides,
    horizontalGuideEvery = args$horizontalGuideEvery,
    hitBorderColour = "grey70",
    hitBorderLinewidth = 0.10,
    parameterBoxColour = "blue3",
    parameterBoxLinewidth = 0.8
  )
  
  ggsave(
    filename = file.path(
      args$outputDir,
      paste0(
        args$fileTag,
        '_HMMscores.pdf'
      )
    ),
    plot = hmmAnalysis$plot,
    width = 16,
    height = 3.6 * hmmAnalysis$nRows,
    units = "in"
  )
  
  updateLog(
    'testHMMs module completed.'
  )
  
  write(
    date(),
    file.path(
      args$outputDir,
      paste0(
        args$fileTag,
        '.done'
      )
    )
  )
}

args <- parser$parse_args()
source(file.path(args$softwareRoot, 'lib', 'common.R'))

tryCatch({
  runModule()
}, error = function(e) {
  cat("ERROR: ", conditionMessage(e), "\n", sep = "", file = stderr())
  flush(stderr())
  quit(save = "no", status = 1, runLast = FALSE)
})
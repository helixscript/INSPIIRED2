library(tidyverse)

in_dir  <- 'U5_sites10000_seed1_error0.00'
out_dir <- '101010_M03249_0302_000000000-SYN00'

dir.create(paste0(out_dir, '/Data/Intensities/BaseCalls'), recursive = TRUE)

o <- read_tsv(file.path(in_dir, 'sampleData.tsv'))

t <- tibble(subject = o$subject,
            sample = o$sample,
            replicate = o$replicate,
            alias = 'x',
            linkerSequence = o$adriftReadLinkerSeq,
            bcSeq = o$index1Seq,
            gender = 'u',
            primer = 'TCTGCGC',
            ltrBit = 'GCTCGCTCGCTCA',
            largeLTRFrag = 'TGAGCGAGCGAGCGCGCAGA',
            vectorSeq = 'synDataTest.fasta',
            refGenome = 'hg38') 

# alias,linkerSequence,bcSeq,gender,primer,ltrBit,largeLTRFrag,vectorSeq,refGenome

i <- which(t$subject == 'subjectA' & t$sample == 'sample1')
t[i,]$alias <- paste0('GTSP9990-', t[i,]$replicate)

t <- dplyr::select(t, -subject, -sample, -replicate)

t2 <- t %>% rbind(names(.), .) %>% `colnames<-`(paste0("X", 1:ncol(.)))

h <- read_csv('resources/fauxMiSeqHeader.csv', show_col_types = FALSE, col_names = FALSE)

t3 <- rbind(h, t2)

write_csv(t3, file.path(out_dir, 'SampleSheet.csv'), col_names = FALSE, na = "")

generate_miseq_ids <- function(n, read = 1) {
  if (!read %in% c(1, 2)) {
    stop("Error: Read flag must be either 1 or 2.")
  }
  
  # Set a realistic upper bound for the Y coordinate to mimic MiSeq optics
  max_y <- 20000 
  
  # Use modulo arithmetic to create a unique, repeating grid of X and Y coordinates.
  # This guarantees no two IDs are identical while keeping numbers under 65,535.
  y_coords <- ((seq_len(n) - 1) %% max_y) + 1000
  x_coords <- ((seq_len(n) - 1) %/% max_y) + 1000
  
  # Randomly assign a few standard MiSeq tile numbers for visual realism
  tiles <- sample(c(1101, 1102, 1103, 2101, 2102, 2103), n, replace = TRUE)
  
  # Construct the static portions of the ID
  base_id <- "M03249:43:000000000-M8R8L:1"
  
  # We can inject the read flag (1:N or 2:N) seamlessly here
  barcode_suffix <- paste0(" ", read, ":N:0:GCNNNNNNNNNN")
  
  # Vectorized paste0 is highly optimized in R for multi-million row character vectors
  ids <- paste0(base_id, ":", tiles, ":", x_coords, ":", y_coords, barcode_suffix)
  
  return(ids)
}

I1 <- readFastq(file.path(in_dir, 'I1.fastq.gz'))

ids <- generate_miseq_ids(length(I1), read = 1)

I1@id <- BStringSet(ids)
writeFastq(I1, file.path(out_dir, 'Data', 'Intensities', 'BaseCalls', 'Undetermined_S0_I1_001.fastq.gz'), compress = TRUE)

R1 <- readFastq(file.path(in_dir, 'R1.fastq.gz'))
R1@id <- BStringSet(ids)
writeFastq(R1, file.path(out_dir, 'Data', 'Intensities', 'BaseCalls', 'Undetermined_S0_R1_001.fastq.gz'), compress = TRUE)

ids <- sub(" 1:N:0:", " 2:N:0:", ids)

R2 <- readFastq(file.path(in_dir, 'R2.fastq.gz'))
R2@id <- BStringSet(ids)
writeFastq(R2, file.path(out_dir, 'Data', 'Intensities', 'BaseCalls', 'Undetermined_S0_R2_001.fastq.gz'), compress = TRUE)

rm(list = ls())


library(tidyverse)

in_dir  <- 'U5_sites5000_seed1_error0.00'
out_dir <- '101010_M03249_0302_000000000-SYN00'

dir.create(paste0(out_dir, '/Data/Intensities/BaseCalls'), recursive = TRUE)

o <- read_tsv(file.path(in_dir, 'sampleData.tsv'))

t <- tibble(subject = o$subject,
            sample = o$sample,
            replicate = o$replicate,
            alias = paste0(o$sample, '-', o$replicate),
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

I1 <- readFastq(file.path(in_dir, 'I1.fastq.gz'))
ids <- readLines('resources/realMiSeq_readIDs_I1.txt.gz')
I1@id <- BStringSet(ids[1:length(I1)])
writeFastq(I1, file.path(out_dir, 'Data', 'Intensities', 'BaseCalls', 'Undetermined_S0_I1_001.fastq.gz'), compress = TRUE)

R1 <- readFastq(file.path(in_dir, 'R1.fastq.gz'))
ids <- readLines('resources/realMiSeq_readIDs_R1.txt.gz')
R1@id <- BStringSet(ids[1:length(R1)])
writeFastq(R1, file.path(out_dir, 'Data', 'Intensities', 'BaseCalls', 'Undetermined_S0_R1_001.fastq.gz'), compress = TRUE)

R2 <- readFastq(file.path(in_dir, 'R2.fastq.gz'))
ids <- readLines('resources/realMiSeq_readIDs_R2.txt.gz')
R2@id <- BStringSet(ids[1:length(R2)])
writeFastq(R2, file.path(out_dir, 'Data', 'Intensities', 'BaseCalls', 'Undetermined_S0_R2_001.fastq.gz'), compress = TRUE)

rm(list = ls())


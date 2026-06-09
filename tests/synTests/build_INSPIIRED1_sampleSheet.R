library(readr)
# Change into each syn run directory, set as cwd, and run this script.

o <- read_tsv('sampleData.tsv')

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
            vectorSeq = 'none.fa',
            refGenome = 'hg38')

i <- which(t$subject == 'subjectA' & t$sample == 'sample1')
t[i,]$alias <- paste0('GTSP9999-', t[i,]$replicate)






t$alias <- gsub('sample1', 'GTSP9999', t$alias)
t$alias <- gsub('sample2', 'GTSP9998', t$alias)
t$alias <- gsub('sample3', 'GTSP9997', t$alias)

t2 <- t %>% rbind(names(.), .) %>% `colnames<-`(paste0("X", 1:ncol(.)))

h <- read_csv('../fauxMisSeqHeader.csv', show_col_types = FALSE, col_names = FALSE)

t3 <- rbind(h, t2)

write_csv(t3, 'SampleSheet_INSPIIRED1.csv', col_names = FALSE, na = "")
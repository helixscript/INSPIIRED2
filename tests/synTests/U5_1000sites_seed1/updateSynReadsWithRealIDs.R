library(readr)
library(ShortRead)

o <- read_tsv('sampleData.tsv')
m <- stringdist::stringdistmatrix(o$index1Seq, o$index1Seq)[1,]
message('Max distance between barcodes: ', min(m[2:length(m)]))

I1 <- readFastq('I1.fastq.gz')
ids <- readLines('../realMiSeq_readIDS_I1.txt.gz')
I1@id <- BStringSet(ids[1:length(I1)])
writeFastq(I1, 'I1.fastq.realIDs.gz')

R1 <- readFastq('R1.fastq.gz')
ids <- readLines('../realMiSeq_readIDS_R1.txt.gz')
R1@id <- BStringSet(ids[1:length(R1)])
writeFastq(R1, 'R1.fastq.realIDs.gz')

R2 <- readFastq('R2.fastq.gz')
ids <- readLines('../realMiSeq_readIDS_R2.txt.gz')
R2@id <- BStringSet(ids[1:length(R2)])
writeFastq(R2, 'R2.fastq.realIDs.gz')


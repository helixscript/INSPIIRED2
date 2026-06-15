library(tidyverse)
library(GenomicRanges)

o <- read_tsv('truth.tsv')
a <- str_split(o$posid, '[\\+\\-]')
d <- data.frame(seqnames = unlist(lapply(a, '[[', 1)),
               start = as.integer(unlist(lapply(a, '[[', 2))),
               end = as.integer(unlist(lapply(a, '[[', 2))),
               strand = stringr::str_extract(o$posid, '[\\+\\-]'))

g <- makeGRangesFromDataFrame(d)

d <- bind_rows(lapply(1:length(g), function(i){
       message(i)
       GenomicRanges::distanceToNearest(g[i], g[-i])
       data.frame(posidIndex = i, dist = mcols(GenomicRanges::distanceToNearest(g[i], g[-i]))$distance)
     }))

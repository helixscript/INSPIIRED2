library(RMySQL)
library(dplyr)
library(readr)

# Input file is sites list need to get GTSP ids.

### input_file  <- 'Jones_limit_of_detection2_INSPIIRED1.tsv'
### output_file <- 'Jones_limit_of_detection2_INSPIIRED1_multiHitClusters.rds'

input_file  <- 'U5_sites10000_seed1_error0.00_INSPIIRED1.tsv'
output_file <- 'U5_sites10000_seed1_error0.00_INSPIIRED1_multiHitClusters.rds'

dbConn  <- dbConnect(MySQL(), group = 'intsites_miseq')
o <- read_tsv(input_file, show_col_types = FALSE)
sample_ids <- unique(o$internalSampleID)
in_clause <- paste0("'", sample_ids, "'", collapse = ", ")

sql <- paste0(
  "SELECT mp.multihitID, s.sampleName, mp.position, mp.chr, mp.strand, COUNT(DISTINCT ml.length) AS unique_lengths ",
  "FROM samples s ",
  "JOIN multihitpositions mp ON s.sampleID = mp.sampleID ",
  "LEFT JOIN multihitlengths ml ON mp.multihitID = ml.multihitID ",
  "WHERE SUBSTRING_INDEX(s.sampleName, '-', 1) IN (", in_clause, ") ",
  "GROUP BY mp.multihitID, s.sampleName, mp.position, mp.chr, mp.strand;"
)

message(sql)
multi <- unique(dbGetQuery(dbConn, sql))
saveRDS(multi, output_file)
dbDisconnect(dbConn)

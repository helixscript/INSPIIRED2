#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'dplyr', 'ShortRead', 'data.table')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",               type = "character",     required = TRUE,            help = "Directory for output files")
parser$add_argument("--inputData",               type = "character",     required = TRUE,            help = "Path to demultiplex module's rds output file.")
parser$add_argument("--dbConfigFile",            type = "character",     default = 'none',           help = "Path to db credential file.")
parser$add_argument("--dbConfigID",              type = "character",     default = 'none',           help = "DB credential block identifier in db credential file.")
args <- parser$parse_args()

args <- list()
args$inputData    <- 'SampleSheet.csv'
args$dbConfigFile <- '/home/everett/.my.cnf'
args$dbConfigID   <- 'specimen_management'

createDBconnection <- function(){
  tryCatch({
    dbConnect(RMariaDB::MariaDB(), group = args$dbConfigID, default.file = args$dbConfigFile)
  },
  error=function(cond) {
    stop(paste0('Error - could not connect to the database. Caught error: ', cond$message))
  })
}

dbConn <- createDBconnection()
d <- dbGetQuery(dbConn, "select * from gtsp")
dbDisconnect(dbConn)

o <- readLines(args$inputData)
o <- o[grepl('NNNNNN', o)]

tab <- read.table(textConnection(o), sep = ',', header = FALSE)

# Assume first row contains sample name.

samples <- bind_rows(lapply(tab[,1], function(x){
  if(grepl('GTSP', x)){
    dd <- d[d$SpecimenAccNum == sub('\\-\\d+$', '', x),]
    if(nrow(dd) != 1) stop(paste0('Error -- could not find information for sample: ', sub('\\-\\d+$', '', x), ' in the sample database.'))
    return(tibble(trial = dd$Trial,
                  subject = dd$Patient,
                  sample = sub('\\-\\d+$', '', x),
                  replicate = str_extract(x, '\\d+$')))
  } else {
    return(tibble(trial = 'Control',
                  subject = 'Control',
                  sample = sub('\\-\\d+$', '', x),
                  replicate = str_extract(x, '\\d+$')))
  }
}))


samples$adriftReadLinkerSeq <- unlist(str_extract_all(o, '[ATCG]+[N]+[ATCG]+'))
samples$index1Seq <- gsub(',', '', unlist(str_extract_all(o, ',[ATCG]{12},')))
samples$refGenome <- gsub(',', '', unlist(str_match_all(o, ',[A-Za-z]+\\d+,?')))
samples$vectorFastaFile <- 'Exuma.fasta'
samples$leaderSeqHMM <- 'genericCART19.hmm'
samples$mode <- 'U5'

readr::write_tsv(samples, 'sampleData.tsv')

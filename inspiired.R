#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'dplyr')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser(prog = "inspiired", description = "inspiired: Automated Vector Integration Analysis")

version_file <- file.path(this.path::this.dir(), "VERSION")
version <- if (file.exists(version_file)) readLines(version_file, n = 1) else NA
parser$add_argument("-v", "--version", action = "version", version = paste("inspiired", version))

subparsers <- parser$add_subparsers(dest = "module", help = "inspiired modules")

# Define parameters for each module.
# The parameters must match the parameters defined at the top of each module except that --softwareRoot flags should be excluded here.

demux_parser <- subparsers$add_parser("demultiplex", help = "Separate reads by barcode")
demux_parser$add_argument("--outputDir",                    type = "character",     required = TRUE,                  help = "Directory for output files")
demux_parser$add_argument("--sampleData",                   type = "character",     required = TRUE,                  help = "Sample definition file")
demux_parser$add_argument("--I1",                           type = "character",     required = TRUE,                  help = "Path to the Index1 read FASTQ file")
demux_parser$add_argument("--R1",                           type = "character",     required = TRUE,                  help = "Path to the Forward read FASTQ file")
demux_parser$add_argument("--R2",                           type = "character",     required = TRUE,                  help = "Path to the Reverse read FASTQ file")
demux_parser$add_argument("--threads",                      type = "integer",       default = 50,                     help = "Number of threads to use.")
demux_parser$add_argument("--fileTag",                      type  = "character",    default = "demultiplex",          help = "String appended to output files in the outpt directory.")
demux_parser$add_argument("--index1ReadMaxMismatch",        type = "integer",       default = 1,                      help = "Number of allowed mismatches to the I1 barcode sequence.")
demux_parser$add_argument("--disableAutoBarcodeOrt",        action = "store_true",  default = FALSE,                  help = "Subsample the data an automatically determine if I1 barcodes need to be reverse complimented.")
demux_parser$add_argument("--disablePostUmiLinker",         action = "store_true",  default = FALSE,                  help = "Disable the requirement to match the post-UMI linker sequence.")
demux_parser$add_argument("--postUmiLinkerMaxMismatch",     type = "integer",       default = 1,                      help = "Number of allowed mismatches to the linker sequence following the UMI sequence.")
demux_parser$add_argument("--qualTrimHalfWidth",            type = "integer",       default = 3,                      help = "Half width of NT window slid along sequence during quality trimming.")
demux_parser$add_argument("--qualTrimEvents",               type = "integer",       default = 2,                      help = "Number of failing events within a window to trigger trimming.")
demux_parser$add_argument("--qualTrimScore",                type = "integer",       default = 10,                     help = "Qual code afterwhich NTs are trimmed.")
demux_parser$add_argument("--polyGfilterPattern",           type = "character",     default = "G{5,}[ATCN]?G{5,}.*$", help = "Pattern to recognize poly-G tail NTs.")
demux_parser$add_argument("--disablePolyGfilter",           action = "store_true",  default = FALSE,                  help = "Disable poly-G filter.")
demux_parser$add_argument("--correctGolayIndexReads",       action = "store_true",  default = FALSE,                  help = "Use a Golay correction algorithm to correct barcode sequences.")
demux_parser$add_argument("--disableAdriftReadLinkers",     action = "store_true",  default = FALSE,                  help = "Use the unique linker sequences on adrift reads with barcode sequences.")
demux_parser$add_argument("--adriftReadLinkerMaxMismatch",  type = "integer",       default = 1,                      help = "Number of allowed mismatches to the linker sequence.")
demux_parser$add_argument("--ramDiskPath",                  type = "character",     default = "/dev/shm",             help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
demux_parser$add_argument("--disableSequenceCollapse",      action = "store_true",  default = FALSE,                  help = "Disable the collapse of duplicate sequences.")
demux_parser$add_argument("--vectorDir",                    type = "character",     default = 'none',                 help = "Path to custom vector files.")
demux_parser$add_argument("--hmmDir",                       type = "character",     default = 'none',                 help = "Path to custom hmm files.")


prp_parser <- subparsers$add_parser("prepReads", help = "Prepare demultiplexed reads for alignment to a reference genome.")
prp_parser$add_argument("--outputDir",               type = "character",     required = TRUE,          help = "Directory for output files")
prp_parser$add_argument("--inputData",               type = "character",     required = TRUE,          help = "Path to demultiplex module's rds output file.")
prp_parser$add_argument("--threads",                 type = "integer",       default = 50,             help = "Number of threads to use.")
prp_parser$add_argument("--fileTag",                 type  = "character",    default = "prepReads",    help = "String appended to output files in the outpt directory.")
prp_parser$add_argument("--ramDiskPath",             type = "character",     default = "/dev/shm",     help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
prp_parser$add_argument("--ORtrimPatternWidth",      type = "integer",       default = 8,              help = "Number of NTs used to build over-reading patterns.")
prp_parser$add_argument("--ORseqMaxMismatch",        type = "double",        default = 0.10,           help = "Max mismatch percentage (0 .. 1) allowed to match over-reading patterns.")
prp_parser$add_argument("--minReadLength",           type = "integer",       default = 30,             help = "Minial read length allowed.")
prp_parser$add_argument("--vectorTestWidth",         type = "integer",       default = 25,             help = "Number of NTs at the end of reads to use to test for vector homology.")
prp_parser$add_argument("--vectorTestMinPercentID",  type  = "double",       default = 90,             help = "Min. perecent ID (0 .. 100) to accept a vector alignment.")
prp_parser$add_argument("--vectorTestMinCoverage",   type = "double",        default = 90,             help = "Min. test sequence converage (0 .. 100) to accept a vector alignment.")
prp_parser$add_argument("--vectorDir",               type = "character",     default = 'none',         help = "Path to custom vector files.")
prp_parser$add_argument("--hmmDir",                  type = "character",     default = 'none',         help = "Path to custom hmm files.")


alr_parser <- subparsers$add_parser("alignReads", help = "Align reads to a reference genome.")
alr_parser$add_argument("--outputDir",               type = "character",     required = TRUE,          help = "Directory for output files")
alr_parser$add_argument("--inputData",               type = "character",     required = TRUE,          help = "Path to demultiplex module's rds output file.")
alr_parser$add_argument("--threads",                 type = "integer",       default = 50,             help = "Number of threads to use.")
alr_parser$add_argument("--fileTag",                 type = "character",     default = "alignReads",   help = "String appended to output files in the outpt directory.")
alr_parser$add_argument("--ramDiskPath",             type = "character",     default = "/dev/shm",     help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
alr_parser$add_argument("--minPercentID",            type = "double",        default = 95,             help = "Min BLAT alignment percent ID score (0 .. 100).")
alr_parser$add_argument("--minAlignmentCoverage",    type = "double",        default = 95,             help = "Min BLAT alignment percent query coverage (0 .. 100).")
alr_parser$add_argument("--blatStepSize",            type = "integer",       default = 5,              help = "BLAT step size.")
alr_parser$add_argument("--blatTileSize",            type = "integer",       default = 11,             help = "BLAT tile size.")
alr_parser$add_argument("--blatRepMatch",            type = "integer",       default = 3000,           help = "BLAT repMatch value.")
alr_parser$add_argument("--blatMaxtNumInsert",       type = "integer",       default = 1,              help = "BLAT max number of target inserts.")
alr_parser$add_argument("--blatMaxqNumInsert",       type = "integer",       default = 1,              help = "BLAT max number of query inserts.")
alr_parser$add_argument("--blatMaxtBaseInsert",      type = "integer",       default = 1,              help = "BLAT max number of target insert NTs.")
alr_parser$add_argument("--blatMaxqBaseInsert",      type = "integer",       default = 1,              help = "BLAT max number of target insert NTs.")
alr_parser$add_argument("--dataRowChunkSize",        type = "integer",       default = 1000,           help = "Numbers of data rows to process per alignment worker.")

bdf_parser <- subparsers$add_parser("buildFragments", help = "Build genomic fragments from alignment data.")
bdf_parser$add_argument("--outputDir",               type = "character",     required = TRUE,            help = "Directory for output files")
bdf_parser$add_argument("--inputData",               type = "character",     required = TRUE,            help = "Path to demultiplex module's rds output file.")
bdf_parser$add_argument("--threads",                 type = "integer",       default = 50,               help = "Number of threads to use.")
bdf_parser$add_argument("--fileTag",                 type = "character",     default = "buildFragments", help = "String appended to output files in the outpt directory.")
bdf_parser$add_argument("--ramDiskPath",             type = "character",     default = "/dev/shm",       help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
bdf_parser$add_argument("--dataRowChunkSize",        type = "integer",       default = 5000,             help = "Numbers of data rows to process per alignment worker.")
bdf_parser$add_argument("--minFrgamentLength",       type = "integer",       default = 50,               help = "Min. Fragment length.")
bdf_parser$add_argument("--maxFrgamentLength",       type = "integer",       default = 100000L,          help = "Max. Fragment length.")

bsf_parser <- subparsers$add_parser("buildStdFragments",   help = "Standardize genomic fragments.")
bsf_parser$add_argument("--outputDir",                     type = "character",   required = TRUE,                help = "Directory for output files")
bsf_parser$add_argument("--inputData",                     type = "character",   required = TRUE,                help = "Path to demultiplex module's rds output file.")
bsf_parser$add_argument("--threads",                       type = "integer",     default  = 30,                  help = "Number of threads to use.")
bsf_parser$add_argument("--fileTag",                       type = "character",   default  = "buildStdFragments", help = "String appended to output files in the outpt directory.")
bsf_parser$add_argument("--ramDiskPath",                   type = "character",   default  = "/dev/shm",          help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
bsf_parser$add_argument("--clusterLeaderSeqs",           action = "store_true",  default  = FALSE,               help = 'Cluster leader sequences and consider when building fragments.')
bsf_parser$add_argument("--disableBreakPointPosStd",     action = "store_true",  default  = FALSE,               help = 'Disable break point standardization.')
bsf_parser$add_argument("--disableIntSitePosStd",        action = "store_true",  default  = FALSE,               help = 'Disable intSite position standardization.')
bsf_parser$add_argument("--anchorReadClusterLen",          type = "integer",     default  = 30,                  help = 'Length of anchor read sequences to test for rearrangments.')
bsf_parser$add_argument("--anchorReadClusterMinAbundDiff", type = "integer",     default  = 5,                   help = 'When clustering anchor read sequences, min. difference between 1st and 2nd most abundant sequence clusters to pick a winner.')
bsf_parser$add_argument("--anchorReadClusterMinReadMult",  type = "integer",     default  = 10,                  help = 'When clustering anchor read sequences, multiplier for 1st and 2nd most read sequence clusters to pick a winner.')
bsf_parser$add_argument("--UMIclusterMinReadMult",         type = "integer",     default  = 5,                   help = 'When clustering UMI sequences, multiplier for 1st and 2nd most read sequence clusters to pick a winning target.')
bsf_parser$add_argument("--minReadsPerFrag",               type = "integer",     default  = 1,                   help = 'Min. number of reads to accept a fragment.')
bsf_parser$add_argument("--UMIclusteringParams",           type = "character",   default  =  "-c 0.80 -d 0 -M 0 -g 1 -r 0 -n 4 -G 1",                                       help = 'Clustering params for clustering UMIs.')
bsf_parser$add_argument("--leaderSeqClusteringParams",     type = "character",   default  =  "-c 0.87 -d 0 -M 0 -g 0 -r 0 -n 5 -G 1 -aS 0.80",                              help = 'Clustering params for clustering leader sequences.')
bsf_parser$add_argument("--multiHitclusteringParams",      type = "character",   default  =  "-c 0.87 -d 0 -M 0 -g 0 -r 0 -n 5 -G 1 -gap -5 -gap-ext -1 -aS 0.93",          help = 'Clustering params for clustering building multi-hit clusters.')
bsf_parser$add_argument("--anchorReadClusterParams",       type = "character",   default  =  "-c 0.87 -d 0 -M 0 -g 0 -r 0 -n 5 -G 1 -gap -5 -gap-ext -2 -aS 0.93 -aL 0.93", help = 'Clustering params for clustering the start of anchor read sequences.')




bst_parser <- subparsers$add_parser("buildSites",   help = "Assemble standardized fragments into integration events.")
bst_parser$add_argument("--outputDir",               type = "character",     required = TRUE,         help = "Directory for output files")
bst_parser$add_argument("--inputData",               type = "character",     required = TRUE,         help = "Path to demultiplex module's rds output file.")
bst_parser$add_argument("--threads",                 type = "integer",       default = 50,            help = "Number of threads to use.")
bst_parser$add_argument("--fileTag",                 type = "character",     default = "buildSites",  help = "String appended to output files in the outpt directory.")
bst_parser$add_argument("--ramDiskPath",             type = "character",     default = "/dev/shm",    help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
bst_parser$add_argument("--disableDualDetect",       action = "store_true",  default = FALSE,         help = "Diable the merging of U5 and U3 samples into dual-detection sites.")
bst_parser$add_argument("--dualDetectWidth",         type = "integer",       default = 6,             help = "Radius for searching for dual-detections.")
bst_parser$add_argument("--integraseCorrectionDist", type = "integer",       default = 2,             help = "Integrase correction value (NT) to account for gDNA duplication caused by integration.")
bst_parser$add_argument("--sumSonicBreaksWithin",    type = "character",     default = "samples",     help = "Sum sonic breaks within either 'samples' (default) or within sample 'replicates'.") 
bst_parser$add_argument("--leadSeqClusteringParms",  type = "character",     default = "-c 0.90 -n 5 -G 0 -aS 0.95 -gap -2 -gap-ext -1 -d 0 -M 0", help = "CLustering parameters used to determine representative leaders sequence.")

if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  parser$print_help()
  quit(status = 0)
}

args <- parser$parse_args()

if (is.null(args$module)) {
  parser$print_help()
  quit(status = 1)
}

pipeline_root <- this.path::this.dir()

module_script <- file.path(pipeline_root, "modules", paste0(args$module, ".R"))

clean_args <- args[names(args)!= "module"]

cmd_args <- sapply(names(clean_args), function(n) {
  val <- clean_args[[n]]
  if (is.logical(val)) {
    if (val) return(paste0("--", n)) else return("")
  } else {
    return(paste0("--", n, " ", shQuote(as.character(val))))
  }
})

final_cmd <- paste(cmd_args[cmd_args != ""], collapse = " ")

system2("Rscript", args = c("--vanilla", module_script, final_cmd, "--softwareRoot", shQuote(pipeline_root)))
quit(status = 0)

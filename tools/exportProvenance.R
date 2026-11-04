#!/usr/bin/env -S Rscript --vanilla
suppressPackageStartupMessages(library(argparse))

scriptArgument <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
defaultRoot <- if(length(scriptArgument))
  dirname(dirname(normalizePath(sub('^--file=', '', scriptArgument[[1L]]), mustWork = TRUE))) else getwd()

parser <- ArgumentParser(description = 'Export embedded INSPIIRED2 processing provenance as YAML.')
parser$add_argument('--inputData', type = 'character', required = TRUE,
                    help = 'Input RDS or Parquet file. Parquet data rows are not loaded.')
parser$add_argument('--outputFile', type = 'character', required = TRUE,
                    help = 'Path for the exported provenance YAML.')
parser$add_argument('--softwareRoot', type = 'character', default = defaultRoot,
                    help = 'INSPIIRED2 installation root; defaults to the parent of this tools directory.')

runModule <- function(){
  if(!file.exists(args$inputData) || dir.exists(args$inputData) || file.access(args$inputData, 4L) != 0L)
    stop('Error - input file is missing or unreadable: ', args$inputData, call. = FALSE)
  inputPath <- normalizePath(args$inputData, mustWork = TRUE)
  outputDirectory <- dirname(args$outputFile)
  if(!dir.exists(outputDirectory) && !dir.create(outputDirectory, recursive = TRUE))
    stop('Error - could not create output directory: ', outputDirectory, call. = FALSE)
  outputPath <- file.path(normalizePath(outputDirectory, mustWork = TRUE), basename(args$outputFile))
  if(identical(inputPath, normalizePath(outputPath, mustWork = FALSE)))
    stop('Error - provenance output must not replace the input data file.', call. = FALSE)
  if(dir.exists(outputPath)) stop('Error - outputFile is a directory.', call. = FALSE)

  extension <- tolower(tools::file_ext(inputPath))
  if(extension == 'rds'){
    object <- readRDS(inputPath)
    bundle <- if(is.null(attr(object, 'inspiired2.provenance', exact = TRUE)))
      provenanceFromObject(object, source = list(path = inputPath, format = 'rds')) else
      provenanceFromObject(object)
  } else if(extension == 'parquet'){
    if(!requireNamespace('arrow', quietly = TRUE))
      stop('Error - the arrow package is required to inspect a Parquet file.', call. = FALSE)
    reader <- arrow::ParquetFileReader$create(inputPath)
    text <- reader$GetSchema()$metadata[['inspiired2.provenance']]
    # Inspect only the footer. Existing provenance is exported without appending
    # an execution/source node or changing the archived history.
    bundle <- if(is.null(text))
      provenanceFromObject(NULL, source = list(
        path = inputPath, format = 'parquet',
        note = 'No embedded processing provenance; data rows were not loaded.')) else
      provenanceFromYAML(text)
  } else {
    stop('Error - inputData must have a .rds or .parquet extension.', call. = FALSE)
  }

  temporary <- tempfile('.provenance_', tmpdir = outputDirectory, fileext = '.yml')
  on.exit(unlink(temporary), add = TRUE)
  cat(provenanceToYAML(bundle), file = temporary)
  if(!file.rename(temporary, outputPath))
    stop('Error - could not install the exported provenance file: ', outputPath, call. = FALSE)
  cat('Exported ', length(bundle$records), ' provenance record(s) to ', outputPath, '\n', sep = '')
}

#-------------------------------------------------------------------------------

args <- parser$parse_args()

tryCatch({
  source(file.path(args$softwareRoot, 'lib', 'common.R'))
  runModule()
}, error = function(e){
  cat('ERROR: ', conditionMessage(e), '\n', sep = '', file = stderr())
  flush(stderr())
  quit(save = 'no', status = 1, runLast = FALSE)
})

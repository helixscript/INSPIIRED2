#!/usr/bin/env -S Rscript --vanilla
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = 'Check sample definitions and resources before processing; report all detected problems.')
parser$add_argument('--outputDir', type = 'character', required = TRUE, help = 'Directory for the validation report; created if necessary.')
parser$add_argument('--fileTag', type = 'character', default = 'validateSampleData', help = 'Report filename prefix; writes <outputDir>/<fileTag>.txt.')
parser$add_argument('--sampleData', type = 'character', required = TRUE, help = 'Sample definition TSV used by demultiplex.')
parser$add_argument('--dbConfigFile', type = 'character', default = 'none', help = 'Optional database credential file; requires --dbConfigID.')
parser$add_argument('--dbConfigID', type = 'character', default = 'none', help = 'Optional credential group; requires --dbConfigFile.')
parser$add_argument('--resourceDir', type = 'character', default = '/resources', help = 'Resource overlay containing hmms/, referenceGenomes/, and vectors/; bundled data/ is also checked.')
parser$add_argument('--verbose', action = 'store_true', default = FALSE, help = 'Include resource search locations in the report.')
parser$add_argument('--softwareRoot', type = 'character', required = TRUE, help = 'Path to the INSPIIRED2 installation.')

args <- parser$parse_args()

runModule <- function(){
  if(!nzchar(trimws(args$fileTag)) || grepl('[[:cntrl:]/\\\\]', args$fileTag))
    stop('--fileTag must be a non-empty filename prefix without path separators or control characters.', call. = FALSE)
  if(!dir.exists(args$outputDir)) dir.create(args$outputDir, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(args$outputDir)) stop('Could not create --outputDir: ', args$outputDir, call. = FALSE)
  reportFile <- file.path(args$outputDir, paste0(args$fileTag, '.txt'))
  inputs <- c(args$sampleData, if(args$dbConfigFile != 'none') args$dbConfigFile)
  if(normalizePath(reportFile, mustWork = FALSE) %in% normalizePath(inputs, mustWork = FALSE))
    stop('The report path must not overwrite the sample data or credential file. Change --outputDir or --fileTag.', call. = FALSE)

  # Clear an older result before starting so an interrupted run cannot leave a
  # previous PASS report looking like the result of the current validation.
  writeLines(c(paste0('FAILED: ', args$sampleData), 'Validation did not complete.'), reportFile)
  report <- tryCatch({
    source(file.path(args$softwareRoot, 'lib', 'validateSampleData.R'))
    validateSampleData(args)
  }, error = identity)

  output <- file(reportFile, open = 'wt', encoding = 'UTF-8')
  on.exit(close(output), add = TRUE)
  if(inherits(report, 'error')){
    writeLines(c(paste0('FAILED: ', args$sampleData),
                 paste0('Validation could not complete: ', conditionMessage(report))), output)
    return(1L)
  }
  printSampleDataValidation(report, args$sampleData, args$verbose, output = output)
  as.integer(length(report$issues) > 0L)
}

status <- tryCatch(runModule(), error = function(e){
  cat('ERROR: Could not complete sample validation or save its report: ', conditionMessage(e), '\n', sep = '', file = stderr())
  flush(stderr())
  1L
})
quit(save = 'no', status = status, runLast = FALSE)

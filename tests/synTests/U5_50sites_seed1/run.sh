#!/bin/bash
set -euo pipefail

mkdir -p output
rm -rf output/*

inspiired2 demultiplex --outputDir   output          \
                       --threads 10                  \
                       --sampleData  sampleData.tsv  \
                       --indexReads  I1.fastq.gz     \
                       --adriftReads R1.fastq.gz     \
                       --anchorReads R2.fastq.gz
inspiired2 prepReads         --outputDir output  --threads 10 --inputData output/demultiplex.rds
inspiired2 alignReads        --outputDir output  --threads 10 --inputData output/prepReads.rds
inspiired2 buildFragments    --outputDir output  --threads 10 --inputData output/alignReads.rds
inspiired2 buildStdFragments --outputDir output  --threads 10 --inputData output/buildFragments.rds
inspiired2 buildSites        --outputDir output  --threads 10 --inputData output/buildStdFragments.rds
inspiired2 nearestGenes      --outputDir output  --threads 10 --inputData output/buildSites.rds
inspiired2 annotateRepeats   --outputDir output  --threads 10 --inputData output/nearestGenes.rds

output_file="output/annotateRepeats.rds"
md5_file="expected_md5sum"

# Execution timestamps/IDs are deliberately different on each run. Compare the
# original scientific object after removing ONLY the new provenance attribute.
# Leave the archived output intact, and require its provenance to be valid.
software_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
actual_md5=$(Rscript --vanilla - "$output_file" "$software_root/lib/common.R" <<'RS'
argv <- commandArgs(trailingOnly = TRUE)
source(argv[[2L]])
local({
  output <- readRDS(argv[[1L]])
  if(is.null(attr(output, "inspiired2.provenance", exact = TRUE)))
    stop("The synthetic output is missing processing provenance.", call. = FALSE)
  provenanceFromObject(output)  # Validate; malformed metadata must fail the test.
  canonical <- tempfile(fileext = ".rds")
  on.exit(unlink(canonical), add = TRUE)
  saveRDS(provenanceStrip(output), canonical)
  cat(unname(tools::md5sum(canonical)))
})
RS
)
expected_md5=$(awk '{print $1}' "$md5_file")

if [[ "$actual_md5" == "$expected_md5" ]]; then
    echo "PASS"
    echo "DATA MD5 MATCH (provenance excluded): $actual_md5"
    exit 0
else
    echo "FAILED"
    echo "DATA MD5 MISMATCH (provenance excluded)"
    echo "  Expected: $expected_md5"
    echo "  Actual:   $actual_md5"
    exit 1
fi

#!/bin/bash
set -euo pipefail

mkdir -p output
rm -rf output/*

inspiired2 demultiplex --outputDir   output          \
                       --sampleData  sampleData.tsv  \
                       --indexReads  I1.fastq.gz     \
                       --adriftReads R1.fastq.gz     \
                       --anchorReads R2.fastq.gz
inspiired2 prepReads         --outputDir output  --inputData output/demultiplex.rds
inspiired2 alignReads        --outputDir output  --inputData output/prepReads.rds
inspiired2 buildFragments    --outputDir output  --inputData output/alignReads.rds
inspiired2 buildStdFragments --outputDir output  --inputData output/buildFragments.rds
inspiired2 buildSites        --outputDir output  --inputData output/buildStdFragments.rds

output_file="output/buildSites.rds"
md5_file="expected_md5sum"

actual_md5=$(md5sum "$output_file" | awk '{print $1}')
expected_md5=$(awk '{print $1}' "$md5_file")

if [[ "$actual_md5" == "$expected_md5" ]]; then
    echo "MD5 MATCH: $actual_md5"
else
    echo "MD5 MISMATCH"
    echo "  Expected: $expected_md5"
    echo "  Actual:   $actual_md5"
fi

#!/bin/bash
set -euo pipefail

inspiired demultiplex --outputDir output          \
                      --sampleData sampleData.tsv \
                      --I1 I1.fastq.gz            \
                      --R1 R1.fastq.gz            \
                      --R2 R2.fastq.gz
inspiired prepReads --outputDir output  --inputData output/demultiplex.rds
inspiired alignReads --outputDir output  --inputData output/prepReads.rds
inspiired buildFragments --outputDir output  --inputData output/alignReads.rds
inspiired buildStdFragments --outputDir output  --inputData output/buildFragments.rds
inspiired buildSites --outputDir output  --inputData output/buildStdFragments.rds

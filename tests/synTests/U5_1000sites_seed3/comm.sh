#!/bin/bash
set -euo pipefail

inspiired demultiplex --outputDir out             \
                      --sampleData sampleData.tsv \
                      --I1 I1.fastq.gz            \
                      --R1 R1.fastq.gz            \
                      --R2 R2.fastq.gz
inspiired prepReads --outputDir out  --inputData out/demultiplex.rds
inspiired alignReads --outputDir out  --inputData out/prepReads.rds
inspiired buildFragments --outputDir out  --inputData out/alignReads.rds
inspiired buildStdFragments --outputDir out  --inputData out/buildFragments.rds
inspiired buildSites --outputDir out  --inputData out/buildStdFragments.rds

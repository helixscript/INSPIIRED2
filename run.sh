#!/bin/bash
set -euo pipefail

inspiired2 demultiplex --outputDir   out             \
                       --sampleData  sampleData.tsv  \
                       --indexReads  I1.fastq.gz     \
                       --adriftReads R1.fastq.gz     \
                       --anchorReads R2.fastq.gz
                       
inspiired2 prepReads         --outputDir out  --inputData out/demultiplex.rds
inspiired2 alignReads        --outputDir out  --inputData out/prepReads.rds
inspiired2 buildFragments    --outputDir out  --inputData out/alignReads.rds
inspiired2 buildStdFragments --outputDir out  --inputData out/buildFragments.rds
inspiired2 buildSites        --outputDir out  --inputData out/buildStdFragments.rds
inspiired2 nearestGenes      --outputDir out  --inputData out/buildSites.rds
inspiired2 annotateRepeats   --outputDir out  --inputData out/nearestGenes.rds

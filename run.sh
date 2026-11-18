#!/bin/bash
set -euo pipefail

inspiired2 validateSampleData --outputDir out --sampleData sampleData.tsv

inspiired2 demultiplex --outputDir out                                      \
                       --sampleData  sampleData.tsv                          \
                       --indexReads  Undetermined_S0_I1_001.fastq.gz         \
                       --adriftReads Undetermined_S0_R1_001.fastq.gz        \
                       --anchorReads Undetermined_S0_R2_001.fastq.gz

inspiired2 testHMMs                 --outputDir out  --inputData out/demultiplex.rds
inspiired2 prepReads                --outputDir out  --inputData out/demultiplex.rds
inspiired2 alignReads               --outputDir out  --inputData out/prepReads.rds
inspiired2 buildFragments           --outputDir out  --inputData out/alignReads.rds
inspiired2 buildStdFragments        --outputDir out  --inputData out/buildFragments.rds
inspiired2 annotateMultiHitRepeats  --outputDir out  --inputData out/buildStdFragments_multiHitClusters.rds
inspiired2 buildSites               --outputDir out  --inputData out/buildStdFragments.rds
inspiired2 nearestGenes             --outputDir out  --inputData out/buildSites.rds
inspiired2 annotateRepeats          --outputDir out  --inputData out/nearestGenes.rds

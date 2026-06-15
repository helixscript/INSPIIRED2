#!/bin/bash
set -euo pipefail

inspiired2 demultiplex --outputDir output                       \
                       --sampleData sampleData.tsv              \
                       --I1 I1.fastq.gz     \
                       --R1 R1.fastq.gz     \
                       --R2 R2.fastq.gz

inspiired2 prepReads --outputDir output  --inputData output/demultiplex.rds
inspiired2 alignReads --outputDir output  --inputData output/prepReads.rds
inspiired2 buildFragments --outputDir output  --inputData output/alignReads.rds
inspiired2 buildStdFragments --outputDir output  --inputData output/buildFragments.rds --anchorReadClusterGrouping sample
inspiired2 buildSites --outputDir output  --inputData output/buildStdFragments.rds --disableDualDetect
inspiired2 nearestGenes --outputDir output  --inputData output/buildSites.rds
inspiired2 annotateRepeats --outputDir output --inputData output/nearestGenes.rds

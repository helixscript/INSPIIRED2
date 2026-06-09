
source inspiired2/bin/activate

./buildSynData.py --seed 1 --outputDir U5_1000sites_seed1
./buildSynData.py --seed 2 --outputDir U5_1000sites_seed2
./buildSynData.py --seed 3 --outputDir U5_1000sites_seed3

./evalSynDataResult.R -s U5_1000sites_seed1/out/buildSites.rds -m U5_1000sites_seed1/out/buildStdFragments_multHitClusters.rds -t U5_1000sites_seed1/truth.tsv -o U5_1000sites_seed1
./evalSynDataResult.R -s U5_1000sites_seed2/out/buildSites.rds -m U5_1000sites_seed2/out/buildStdFragments_multHitClusters.rds -t U5_1000sites_seed2/truth.tsv -o U5_1000sites_seed2
./evalSynDataResult.R -s U5_1000sites_seed3/out/buildSites.rds -m U5_1000sites_seed3/out/buildStdFragments_multHitClusters.rds -t U5_1000sites_seed3/truth.tsv -o U5_1000sites_seed3



Faux testing seq directories.
101010_M03249_0302_000000000-ZZ101
101010_M03249_0302_000000000-ZZ102
101010_M03249_0302_000000000-ZZ103


scp SampleSheet_INSPIIRED1.csv microb120:/media/sequencing/Illumina/101010_M03249_0302_000000000-ZZ101/SampleSheet.csv
scp I1.fastq.realIDs.gz  microb120:/media/sequencing/Illumina/101010_M03249_0302_000000000-ZZ101/Data/Intensities/BaseCalls/Undetermined_S0_I1_001.fastq.gz
scp R1.fastq.realIDs.gz  microb120:/media/sequencing/Illumina/101010_M03249_0302_000000000-ZZ101/Data/Intensities/BaseCalls/Undetermined_S0_R1_001.fastq.gz
scp R2.fastq.realIDs.gz  microb120:/media/sequencing/Illumina/101010_M03249_0302_000000000-ZZ101/Data/Intensities/BaseCalls/Undetermined_S0_R2_001.fastq.gz
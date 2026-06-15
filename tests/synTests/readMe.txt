
source inspiired2/bin/activate

./buildSynData.py --seed 1 --nSites 10000 --outputDir U5_sites10000_seed1_error0.00
./buildSynData.py --seed 1 --nSites 10000 --percentGenomicError 0.01 --outputDir U5_sites10000_seed1_error0.01
./buildSynData.py --seed 1 --nSites 10000 --percentGenomicError 0.02 --outputDir U5_sites10000_seed1_error0.02
./buildSynData.py --seed 1 --nSites 10000 --percentGenomicError 0.03 --outputDir U5_sites10000_seed1_error0.03

# Use convert_test_to_INSPIIRED1.R to convert each synthetic data sets to a format INSPIIRED1 can process. 
# Use script to create 101010_M03249_0302_000000000-SYN00, 101010_M03249_0302_000000000-SYN01, 101010_M03249_0302_000000000-SYN02, and 101010_M03249_0302_000000000-SYN03.







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

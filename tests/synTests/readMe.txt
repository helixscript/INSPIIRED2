
source inspiired2/bin/activate

./buildSynData.py --seed 1 --outputDir U5_1000sites_seed1
./buildSynData.py --seed 2 --outputDir U5_1000sites_seed2
./buildSynData.py --seed 3 --outputDir U5_1000sites_seed3

./evalSynDataResult.R -s U5_1000sites_seed1/out/buildSites.rds -m U5_1000sites_seed1/out/buildStdFragments_multHitClusters.rds -t U5_1000sites_seed1/truth.tsv -o U5_1000sites_seed1
./evalSynDataResult.R -s U5_1000sites_seed2/out/buildSites.rds -m U5_1000sites_seed2/out/buildStdFragments_multHitClusters.rds -t U5_1000sites_seed2/truth.tsv -o U5_1000sites_seed2
./evalSynDataResult.R -s U5_1000sites_seed3/out/buildSites.rds -m U5_1000sites_seed3/out/buildStdFragments_multHitClusters.rds -t U5_1000sites_seed3/truth.tsv -o U5_1000sites_seed3


source inspiired2/bin/activate

./buildSynData.py --seed 1 --nSites 5000 --nFrags 3 --nReadsPerFrag 3 --singleSample --nReplicates 4 --outputDir U5_sites5000_seed1_error0.00
./buildSynData.py --seed 1 --nSites 5000 --nFrags 3 --nReadsPerFrag 3 --singleSample --nReplicates 4 --outputDir U5_sites5000_seed1_error0.01 --percentGenomicError 0.01
./buildSynData.py --seed 1 --nSites 5000 --nFrags 3 --nReadsPerFrag 3 --singleSample --nReplicates 4 --outputDir U5_sites5000_seed1_error0.02 --percentGenomicError 0.02
./buildSynData.py --seed 1 --nSites 5000 --nFrags 3 --nReadsPerFrag 3 --singleSample --nReplicates 4 --outputDir U5_sites5000_seed1_error0.03 --percentGenomicError 0.03

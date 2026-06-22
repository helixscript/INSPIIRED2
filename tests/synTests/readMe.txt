
source inspiired2/bin/activate

./buildSynData.py --seed 1 --nSites 100000 --nFrags 4 --singleSample --outputDir U5_sites100000_seed1_error0.00
./buildSynData.py --seed 1 --nSites 100000 --nFrags 4 --singleSample --percentGenomicError 0.01 --outputDir U5_sites100000_seed1_error0.01
./buildSynData.py --seed 1 --nSites 100000 --nFrags 4 --singleSample --percentGenomicError 0.02 --outputDir U5_sites100000_seed1_error0.02
./buildSynData.py --seed 1 --nSites 100000 --nFrags 4 --singleSample --percentGenomicError 0.03 --outputDir U5_sites100000_seed1_error0.03


./buildSynData.py --seed 1 --nSites 50000 --nFrags 4 --singleSample --outputDir U5_sites50000_seed1_error0.00

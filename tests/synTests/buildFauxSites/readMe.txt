
source inspiired2/bin/activate

./buildSynData.py --seed 1 --nSites 25000 --nFrags 1 --nReadsPerFrag 1 --singleSample --nReplicates 1 --positionChatterSD 0 --percentGenomicError 0    --refGenomePath ../../../data/referenceGenomes/hg38.2bit --refGenomeID hg38 --outputDir hg38_U5_25Ksites_singleFrag_singleRead_0error    --illuminaReadIDs

./buildSynData.py --seed 1 --nSites 25000 --nFrags 1 --nReadsPerFrag 1 --singleSample --nReplicates 1 --positionChatterSD 0 --percentGenomicError 0.01 --refGenomePath ../../../data/referenceGenomes/hg38.2bit --refGenomeID hg38 --outputDir hg38_U5_25Ksites_singleFrag_singleRead_0.01error --illuminaReadIDs

./buildSynData.py --seed 1 --nSites 25000 --nFrags 1 --nReadsPerFrag 1 --singleSample --nReplicates 1 --positionChatterSD 0 --percentGenomicError 0.03 --refGenomePath ../../../data/referenceGenomes/hg38.2bit --refGenomeID hg38 --outputDir hg38_U5_25Ksites_singleFrag_singleRead_0.03error --illuminaReadIDs

./buildSynData.py --seed 1 --nSites 25000 --nFrags 1 --nReadsPerFrag 1 --singleSample --nReplicates 1 --positionChatterSD 0 --percentGenomicError 0.05 --refGenomePath ../../../data/referenceGenomes/hg38.2bit --refGenomeID hg38 --outputDir hg38_U5_25Ksites_singleFrag_singleRead_0.05error --illuminaReadIDs

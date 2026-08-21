# System requirements  
INSPIIRED2 is designed to run on mid- to large- scale computational servers with a minimum of 30 cores and 100GB RAM. Smaller servers can be used as well though the number of computational cores used by each module should be reduced from the default of 50 with the ```--threads``` flag. Memory requirement scales with the number of computational cores employed. In order to ensure consistent results, it is recommended that INSPIIRED is executed from the provided Docker image,  eg.

```docker run --rm --shm-size=10g -v /home/everett:/workspace -w /workspace inspiired2 bash run.sh``` 

<br>

# Quick Start  
#### Download the INSPIIRED Docker image.
```wget https://bushmanlab.org/export/inspiired2_latest.tar.gz```

#### Load the image
```docker load -i inspiired2_latest.tar.gz```

#### Process a small 50 site synthetic data set from within the Docker container.
```
docker run --rm -it --shm-size=5g inspiired2 bash
%> cd /opt/INSPIIRED2/tests/synTests/U5_50sites_seed1
%> ./run.sh
```

# Preparing your data for analysis 
Paired-end Illumina formatted sequencing data, a sample data file describing samples, and a pipeline script need to be shared with the INSPIIRED2 Docker container.  The [sample data file](sampleData.tsv) provides INSPIIRED2 information about sequencing bar code and linker sequences. Sequenced amplicons are expected to have the structure defined in the [2016 INSPIIRED paper](https://pubmed.ncbi.nlm.nih.gov/28344990). The [pipeline processing script](run.sh) contains calls to pipelines modules. Each module outputs a single final result RDS formatted file which serves as input for downstream modules such that modules can be daisy-chained together: 
 
```
#!/bin/bash
set -euo pipefail

inspiired2 demultiplex --outputDir out             \
                       --sampleData sampleData.tsv \
                       --I1 I1.fastq.gz            \
                       --R1 R1.fastq.gz            \
                       --R2 R2.fastq.gz
inspiired2 prepReads         --outputDir out  --inputData out/demultiplex.rds
inspiired2 alignReads        --outputDir out  --inputData out/prepReads.rds
inspiired2 buildFragments    --outputDir out  --inputData out/alignReads.rds
inspiired2 buildStdFragments --outputDir out  --inputData out/buildFragments.rds
inspiired2 buildSites        --outputDir out  --inputData out/buildStdFragments.rds
inspiired2 nearestGenes      --outputDir out  --inputData out/buildSites.rds
inspiired2 annotateRepeats   --outputDir out  --inputData out/nearestGenes.rds
```

### Setting up the sample data file 

Reads originating from within LTR sequences and transverse genomic junctures are referred to as *anchor reads* because they anchor sequencing reads to integration positions. Reads originating from within ligated linkers at the opposite ends of fragments are referred to as *adrift reads* because their alignment positions drift due to the genome being sheared during library preparation. For each sample replicate, the sample data file will need the sequence of the adrift linker (eg. GTTAAAGGTGTTCCCTGCCGNNNNNNNNNNNNCTCCGCTTAAGGGACT) and I1 barcode (eg. ACCTAAGTCCGT).

<p align="center">
  <img src="figures/fragmentStructure.png" />
</p>

In addition to this sequence information, the sample configuration file needs information about the reference genome against which to align your data (refGenome), information about your vector (vectorFastaFile), the name of them HMM file about needed to recognize the ends of LTR sequences (leaderSeqHMM), and processing details (mode). Importantly, leaderSeqHMM is the name of a file stored within the software's  ```data/hmms``` folder. The software is provided with HMMs that recognize the ends of HIV-1's 3' and 5' LTR sequences. Custom vector sequences and HMM profiles can be imported at run time using the ```--vectorDir``` and ```--hmmDir``` flags. The design of HMMs that match your vectors is discussed in the 'Working with HMMs' section below. 


### Start processing in Docker

The pipeline writes intermediate results directly to memory. The maximum amount of allowed memory for this scratch space is defined within the Docker call ```--shm-size=10g```. Larger data sets may require the allocation of more memory. In the following Docker call, the -v flag is binding where your data is located on your server, e.g. ```/home/everett```, to the Docker container directory ```/workspace```. The -w flag tells Docker that file paths will be relative to /workspace within the container.  ```run.sh ``` is the name of the pipeline processing script.

```docker run --rm --shm-size=10g -v /home/everett:/workspace -w /workspace inspiired2 bash run.sh```

The software is provided with the following reference genomes:

 - hg38 (human)
 - hs1 (human T2T draft [2022])
 - mm9 (mouse)
 - canFam4 (dog)
 - macFas5 (crab-eating macaque) 
 - chlSab2 (green monkey)
 - sacCer3 (Brewer's yeast)
 
 <br>

# Working with HMMs
Anchor reads containing the ends of vector LTR sequences are recognized using vector specific HMMs. HMMs are used because them are particularly adept at recognizing mismatches and minor indels that can occur due to natural variation and sequencing error.  Vector HMMs are created with the HMMER software package for each vector used in your analysis. 

Developing HMMs for analyses that expect a fixed length, unchanging LTR sequence before genomic junctures is straight forward. HMMs can be created with a single expected sequence and the range of acceptable HMM scores simply depends on how much sequencing or alignment error should be accepted.  Developing HMMs designed to recognize multiple LTR variants with variable lengths (e.g. capturing wild HIV from patients) is possible as well. 
<br>
#### Creating HMMs to capture a single LTR sequence
Creating a HMM for a single  expected LTR sequence, typically the DNA sequence between where your sequencing primer binds the LTR to the terminal CA preceding genomic junctures, is straight forward. First create a FASTA file for the expected LTR sequence you expect to observe in your anchor reads sequences. 

File: mySeq.fasta
```
>mySeq
GAAAATCTCTAGCA
```

Next, use [HMMER](http://hmmer.org/) to create a HMM using this FASTA file.
```
%> hmmbuild mySeq.hmm mySeq.fasta
```
Now that we created an HMM, we need to determine how to score it. Next, create a second FASTA file containing minor variations in your sequence to see how it affects its HMM score. For example, here we create a file name *mySeqTests.fasta* and make minor changes which we would still consider valid hits.

File: mySeqTests.fasta
``` 
>mySeq
GAAAATCTCTAGCA
>mySeq_1SNP
GAAGATCTCTAGCA
>mySeq_2SNPs
GAAGATCTCAAGCA
>mySeq_1del
GAAAATTCTAGCA
>mySeq_1del_1ins
GAAATCTCTGAGCA
```
After creating  a couple of minor variations in our target sequence, we evaluate the variations with the HMM.  
```
nhmmer --F1 1 --F2 1 --F3 1 -T -5 --incT -5 --nobias --popen 0.15 --pextend 0.05 --tblout mySeqTests.tbl mySeq.hmm mySeqTests.fasta
```

Next, review the output (mySeqTests.tbl) to determine a minimum acceptable score.
  
```
# target name        accession  query name           accession  hmmfrom hmm to alifrom  ali to envfrom  env to  sq len strand   E-value  score  bias  description of target
#------------------- ---------- -------------------- ---------- ------- ------- ------- ------- ------- ------- ------- ------ --------- ------ ----- ---------------------
mySeq                -          mySeq                -                1      14       1      14       1      14      14    +      0.0063    3.7   1.1  -
mySeq_1SNP           -          mySeq                -                1      14       1      14       1      14      14    +       0.016    2.8   0.3  -
mySeq_2SNPs          -          mySeq                -                1      13       1      13       1      14      14    +        0.15    0.5   0.9  -
mySeq_1del_1ins      -          mySeq                -                3      10       2       9       1      14      14    +        0.29   -0.1   0.2  -
mySeq_1del           -          mySeq                -                4      13       3      12       1      13      13    +        0.38   -0.4   1.2  -
```
  
Examine the HMM scores in column 14 (score) and decide the range of scores associated with acceptable mismatches. In this example, a researcher may consider a score between 0.5 and 3.7 acceptable. For convenience, HMM scoring parameters can be saved in a separate file along with HMM files. 

This file needs to have the same name as the HMM file except we replace ".hmm" with ".cfg".  Here is an example:
  
```
HMMminStartPos  1
HMMmaxStartPos  5
HMMminFullBitScore      10
HMMmaxFullBitScore      30
HMMmatchEnd     TRUE
HMMmatchTerminalSeq     CA
HMMmatchEndRadius       2
```

Alternatively, for each HMM defined in your sampelData.tsv file, you can provide these parameters as a comma delimited string where each HMM is separated by a pipe character. HMM parameters provided on the command line will override parameters found in the default  .cfg files.

```
inspiired2 prepReads --outputDir out --inputData out/demultiplex.rds --HMMparam 'HIV1_1-100_U5.hmm,1,5,10,30,TRUE,CA,2|HIV1_1-100_U3_RC.hmm,1,5,30,60,TRUE,CA,2'
```
<br>

#### Creating HMMs to capture multiple LTR sequences

In order to create an HMM that can capture multiple LTR sequences, rather than starting with a single sequence as shown above, start with a FASTA formatted multiple sequence alignment (MSA) of varied LTR sequences. If the sequences contain required elements such as the terminal CA sequence found int retroviral LTR sequences, these elements should be fixed in MSAs.

The same approach for determine a range of acceptable HMM scores used with a single expected sequence can be used here. Alternatively, the testHMMs module can be used to test how HMMs score real anchor reads:

```
inspiired2 testHMMs --outputDir out  --sampleData sampleData.tsv --anchorReads  Undetermined_S0_R2_001.fastq.gz --HMMmatchEnd --HMMmatchTerminalSeq CA --HMMmatchEndRadius 2
```
This module tests sequencing data using the HMM profiles found in the sample data file and provided a graphical output useful for tuning HMM parameters:

<p align="center">
  <img src="figures/HMM_scoring.png" />
</p>

<br>

# Processing with databasing

To include databasing features in your processing, additional parameters need to be passed to Docker. Databasing is only available for the buildFragments module and is implemented by including the ```--dbConfigFile``` ```--dbConfigID``` flags. When databasing is used, fragment records are written both to the the output directory and to a local data warehouse. The data warehouse includes two parts:
  1. A SQL database that holds sample data and a pointers to a local parquet file.
  2. A collection of parquet files containing fragment details.

In order for the files to be written to the data warehouse, this flag must be included in your Docker call:

```--user $(id -u):$(id -g)```

The path to your data warehouse must also be included:

```-v /media/md0/data/inspiired:/data ```

The final Docker call for a run including databasing would be:
  
```docker run --rm --shm-size=10g --user $(id -u):$(id -g) -v /media/md0/data/inspiired:/data -v /home/everett:/workspace -w /workspace inspiired2 bash run.sh```

#### Pulling data from the data warehouse

The pullDBrecords module retrieves fragment records based on trial, subject, and sample indentifiers. Rather than starting the pipeline with the typical calls to the demultiplex, prepReads, alignReads, and buildFragments modules, the pipeline can start with a the pullDBrecords module:

```
inspiired2 pullDBrecords     --outputDir out  --outputFile fragments.rds --dbConfigFile my.cnf --dbConfigID inspiired2 --dataPath  /media/md0/data/inspiired --trials Penn_CART --subjects "p432,p663,p215"
inspiired2 buildStdFragments --outputDir out  --inputData out/buildFragments.rds
inspiired2 buildSites        --outputDir out  --inputData out/buildStdFragments.rds
```
Specific fragment records can be pulled using the --trials,  --subjects, and --samples flags which accept comma delimited lists of identifiers. One or more trial identifier must be provided. The module will treat all other flags as 'all' unless identifiers are provided.

  




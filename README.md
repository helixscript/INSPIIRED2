# INSPIIRED2

INSPIIRED2 identifies vector integration sites in reference genomes from paired-end Illumina short-read data. It is designed for linker-mediated libraries created using the original [INSPIIRED protocol](https://pubmed.ncbi.nlm.nih.gov/28344990). The structure of generated genomic fragments is shown below. The software labels reads originating from the LTR/ITR side of fragments as **anchor reads** because they are anchored to sites of vector integration.  Reads orginating from the linker side of fragments are labeled **adrift reads** because they their genomic align positions drift due random sonic sheering along the genome. 

<p align="center">
  <img src="figures/fragmentStructure.png" alt="INSPIIRED fragment and read structure" />
</p>

- **anchor reads** crosses the vector-genome junction. These reads begin with recognizable vector sequences and then transition into genomic sequences. Reads from short fragments that continue into linker sequences at the other end of fragments are automatically trimmed.

- **adrift reads** begins at the ligated linker and read into genomic sequences from the sonic-shearing boundary. Variation in this boundary provides the primary estimate of clonal abundance.

INSPIIRED2 demultiplexes and trims reads, recognizes vector-terminal sequence with a profile HMM, aligns both mates to a reference genome, constructs genomic fragments, standardizes fragment boundaries, filters likely PCR rearrangements, assembles fragments into sample-level integration sites, and adds gene and repeat annotations.

### System requirements

The supplied Docker image is the recommended execution environment because INSPIIRED2 depends on multiple R/Bioconductor packages and command-line tools including HMMER, BLAT, BLAST+, CD-HIT-EST, and UCSC sequence utilities. The versions of these packages and tools can affect results and using the provided Docker image ensures that results are reproducible. 

INSPIIRED2 uses Linux multicore processing. CPU, memory, disk, and shared-memory requirements depend on library size, read complexity, and the reference genomes used. A server with approximately 30 cores and 100 GB RAM is a reasonable starting point for substantial datasets, but these are not hard minimums. On smaller systems, reduce `--threads` values; memory use generally increases with the number of concurrent workers.

Temporary files are written below `--ramDiskPath`, which defaults to the server's RAM disk `/dev/shm`. If that location is not writable, INSPIIRED2 uses the output directory. The Docker `--shm-size` setting controls the space available in `/dev/shm`.

### Installation

Using the provided Docker image [[download here (7 GB)](https://bushmanlab.org/export/inspiired2_latest.tar.gz)] is the reccomended installation method. The documentation below will refer to this method. Alternatively, the software can be installed directly on your system by cloning this repository and calling inspiired2.R. The software requires a number of R libraries to be preinstalled as well as third party software packages (blastn, nhmmer, cd-hit-est). 


### Required inputs

Inspiired2 requires three inputs. 
  - Paired end sequencing (R1, R2, I1) with readIDs in the same order between files.
  - A [sample data](sampleData.tsv) file matching sample replicates to I1 barcode sequecnes and linker sequences.
  - A [bash script](https://github.com/helixscript/INSPIIRED2/blob/main/run.sh) calling one or more INSPIIRED2 modules.  


### Quick start and validation
Download and load the distributed Docker image:

```
wget https://bushmanlab.org/export/inspiired2_latest.tar.gz
```
```bash
docker load -i inspiired2_latest.tar.gz
```
Run the bundled synthetic test:

```
docker run --rm -it --shm-size=5g inspiired2 bash
```

```
cd /opt/INSPIIRED2/tests/synTests/U5_50sites_seed1
```
```bash
./run.sh
```

The test takes about 5 minutes to complete and report PASS at the end will if the MD5sum of the output matches the expected value. 

### Starting analyses using the Docker image

The basic INSPIIRED2 invocation command has this structure:
```
docker run --rm     \
  --shm-size=20g    \
  -v ./:/workspace  \
  -w /workspace     \
  inspiired2 bash run.sh
```
The`--shm-size` flag defines the max. amount of memory allowed to be used as scratch space during analysis. 20GB is a reasonable value for most moderate size Illumina paied-end data sets. This value should be increased for large data sets and should not reach an appreciable percentage of your total RAM. 

` -v ./:/workspace` mounts your analysis directory to `/workspace` inside of the Docker container. Here we are mounting the current directory `./` to `/workspace` within the Docker container. The analysis directory is expected to contain your sequencing data, sample data file, and processing script (described next).

`-w /workspace` instructs Docker to make all paths relative to `/workspace` within the Docker image.

`inspiired2 bash run.sh` instructs docker to run the processing script [run.sh](https://github.com/helixscript/INSPIIRED2/blob/main/run.sh), located in your analysis directory, in a Docker container created with the `inspiired2` Docker image.

By default, all output files will be owned by root. To change ownership to the user initiating the analysis, add this argument:  `--user "$(id -u):$(id -g)"`

### Resource data files

INSPIIRED2 is provided with a number of reference genomes (hg38, hs1, sacCer3, mm10, canFam4, and macFas5) as well as U3 and U5 LTR HMMs created with data from Los Alamos National laboratories. The `showResources` command can be used to list available resources provided with the Docker image. All genomes and genome annotations were created with the included `tools/buildRefGenomeObjects.R` script. This script accepts UCSC genome IDs and pulls data from UCSC data portals to build required data objects. A local install of RepeatMasker is required to create *.repeatTable.gz files required by the `annotateRepeats` module.  

```
%>docker run --rm inspiired2 bash -c 'inspiired2 showResources' 

+-- data 
+-- genomeAnnotations 
|   +-- hg38.TUs.rds 
|   +-- hg38.exons.rds 
|   +-- hg38.repeatTable.gz 
    ...
|   +-- sacCer3.TUs.rds 
|   +-- sacCer3.exons.rds 
|   \-- sacCer3.repeatTable.gz 
+-- hmms 
|   +-- HIV1_LTR_U3_v1.0.cfg 
|   +-- HIV1_LTR_U3_v1.0.hmm 
|   +-- HIV1_LTR_U5_v1.0.cfg 
|   +-- HIV1_LTR_U5_v1.0.hmm 
    ...
|   +-- validation.cfg 
|   \-- validation.hmm 
+-- referenceGenomes 
|   +-- canFam4.2bit 
|   +-- hg38.2bit 
|   +-- hs1.2bit 
|   +-- macFas5.2bit 
|   +-- mm10.2bit 
|   \-- sacCer3.2bit 
\-- vectors 
 +-- CART19.fasta 
 +-- HXB2.fasta 
 \-- synDataTest.fasta
 ```

In this data tree, reference genomes, stored in the referenceGenomes directory, are stored using the 2bit data format and are named with an identifier followed by '.2bit'. Genome annotations are stored in the genomeAnnotations directory. For each genome identifier, *.TUs.rds files stores gene transcription unit coordinates and *.exons.rds files store gene exon coordinates. Coordinates are stored as GenomicRange objects. Repeat annotations, created by RepeatMasker, are stored in *.repeatTable.gz files. These files contain the compressed tabular output created by RepeatMasker. Vector FASTA file are stored in the vectors directory and HMMs are stored in the hmms directory. Vector FASTA files are used to filter out anchor reads which read into the vector bodies rather than into flanking genomic DNA. Each HMM has a corresponding configuration file where the .hmm suffix of the HMM has been replaced with .cfg (discussed in the HMM section below).


Custom reference genomes, gene annotations, vector sequences, and HMMs can be shared with the Docker image at run time by using an additional mount flag: `-v  ~/data:/resources:ro` to overlay custom files onto INSPIIRED2's data file tree. In this example, we are mounting our local file tree  `~/data` to `/resources` within the Docker container. Custom data must be organized in the same way that data is organized within INSPIIRED's data folder. 
```
%> tree ~/data 
  ~/data  
    └── hmms  
        └── myCustomProfile.hmm
 ```

When data matching INSPIIRED2's data tree is mounted to `/resources` in the Docker container, the data is superimposed onto the INSPIIRED2's data file tree and overwrites existing entries if the same names are used, eg.

```
docker run --rm -v  ~/data:/resources:ro inspiired2 bash -c 'inspiired2 showResources'
```
```
(...)
|   +-- sacCer3.exons.rds 
|   \-- sacCer3.repeatTable.gz 
+-- hmms 
|   +-- HIV1_LTR_U3_v1.0.cfg 
|   +-- HIV1_LTR_U3_v1.0.hmm 
|   +-- HIV1_LTR_U5_v1.0.cfg 
|   +-- HIV1_LTR_U5_v1.0.hmm 
|   +-- generic_CART19_v1.0.cfg 
|   +-- generic_CART19_v1.0.hmm 
|   +-- myCustomProfile.hmm   <===== overlayed data file 
|   +-- validation.cfg 
|   \-- validation.hmm 
+-- referenceGenomes 
|   +-- canFam4.2bit 
|   +-- hg38.2bit
(...)
```


### Sample data file

A tab delimited file defining sample replicate barcode and linker sequences is required to demultiplex sequencing runs. This file is a required parameter for the demultiplex module. An example file is provided with the software [`sampleData.tsv`](sampleData.tsv).

| Column | Description |
|---|---|
| `trial` | Study, experiment, or analysis-group identifier. |
| `subject` | Biological subject identifier. Integration positions are standardized across samples and replicates within each `trial`/`subject`. |
| `sample` | Sample identifier. Final site records are assembled at this level. |
| `replicate` | Technical-replicate identifier. Values must be integers. |
| `index1Seq` | Expected Index 1 barcode sequence. |
| `adriftReadLinkerSeq` | Complete linker sequence at the beginning of the adrift read, with the UMI represented by `N` characters. |
| `refGenome` | Reference identifier matching an installed `<refGenome>.2bit` file. |
| `leaderSeqHMM` | HMM filename, including `.hmm`, used to recognize the vector-terminal sequence in anchor reads. |
| `vectorFastaFile` | Vector FASTA filename used by the internal-vector read filter. |
| `mode` | Vector-end detection mode(`U3` or `U5`). |

Example:

```text
trial	subject	sample	replicate	index1Seq	adriftReadLinkerSeq	refGenome	leaderSeqHMM	vectorFastaFile	mode
trial1	subject01	day0	1	CAGTGGGTCTAA	GAACGAGCACTAGTAAGCCCNNNNNNNNNNNNCTCCGCTTAAGGGACT	hg38	HIV1_1-100_U5.hmm	HXB2.fasta	U5
```

Every `adriftReadLinkerSeq` must contain one contiguous UMI region (Ns) and match, case-insensitively:

Additional requirements:

- Resource names are case-sensitive and must match installed HMM, vector, and reference files.
- Barcodes and linkers should remain distinguishable after the configured mismatch allowances are applied.

<br>

## Standard workflow

The pipeline is a daisy chain: the primary RDS file output of one module becomes the input to the next. Modules are chained together in a shell script which is passed to the Docker image.

```bash
#!/usr/bin/env bash
set -euo pipefail

inspiired2 demultiplex --outputDir out --threads 30 \
  --sampleData sampleData.tsv \
  --indexReads I1.fastq.gz \
  --adriftReads R1.fastq.gz \
  --anchorReads R2.fastq.gz

inspiired2 prepReads         --outputDir out --threads 30 --inputData out/demultiplex.rds
inspiired2 alignReads        --outputDir out --threads 30 --inputData out/prepReads.rds
inspiired2 buildFragments    --outputDir out --threads 30 --inputData out/alignReads.rds
inspiired2 buildStdFragments --outputDir out --threads 30 --inputData out/buildFragments.rds
inspiired2 buildSites        --outputDir out --threads 30 --inputData out/buildStdFragments.rds
inspiired2 nearestGenes      --outputDir out --threads 30 --inputData out/buildSites.rds
inspiired2 annotateRepeats   --outputDir out --threads 30 --inputData out/nearestGenes.rds
```

Setting `set -euo pipefail` at the top of processing script instructs the script to stop at the first failed module.

| Module | Main role | Primary output |
|---|---|---|
| `demultiplex` | Quality-trim and assign reads using Index 1 and linker sequences | `demultiplex.rds` |
| `prepReads` | Recognize vector-terminal sequence and prepare genomic read segments | `prepReads.rds` |
| `alignReads` | Align anchor and adrift genomic sequences to reference genomes | `alignReads.rds` |
| `buildFragments` | Pair compatible mate alignments into candidate physical fragments | `buildFragments.rds` |
| `buildStdFragments` | Standardize boundaries, handle multi-hits, filter PCR artifacts, and collapse fragments | `buildStdFragments.rds` |
| `buildSites` | Assemble sample-level integration sites and calculate abundance | `buildSites.rds` |
| `nearestGenes` | Add gene, exon, and nearest-gene annotations | `nearestGenes.rds` |
| `annotateRepeats` | Add overlapping repeat names and classes | `annotateRepeats.rds` |

`showResources`, `testHMMs`, `buildSeqDataMap`, `testDBconn`, and `pullDBrecords` are supporting commands rather than required stages of the standard chain.

<br>

## General command behavior

Core modules accept the following options:

| Flag | Default | Meaning |
|---|---:|---|
| `--outputDir` | Required | Output directory. It is normally created during module setup. |
| `--inputData` | Required | RDS output from the preceding module. `demultiplex` uses raw input flags instead. |
| `--threads` | `50` | Maximum worker or library thread count.  |
| `--fileTag` | Module name | Base name for output files. Allows modules to be run more than once by changing their output file base names.  |
| `--ramDiskPath` | `/dev/shm` | Scratch filesystem; falls back to `outputDir` when not writable. |






INSPIIRED2 is provided with an SQL database and the ability to create a data warehouse to store data from multiple experiments. Databasing and warehousing is enabled by providing database credentials to the buildFragments module arguments: `--dbConfigFile --dbConfigID`

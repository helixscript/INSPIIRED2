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

### Working with HMMs
Anchor reads containing the ends of vector LTR sequences are recognized using vector specific HMMs. HMMs are used because them are particularly adept at recognizing mismatches and minor indels that can occur due to natural variation and sequencing error.  Vector HMMs are created with the HMMER software package for each vector used in your analysis. To create a vector HMM, first create a FASTA file for the expected vector sequence you expect to observe in your R2 read sequences. This will be the expected sequence observed before transitioning into genomic DNA, eg.

```
docker run -it --rm  inspiired2 bash
```

```
echo -e ">seq\nGAAAATCTCTAGCA" > test.ff
```

Next, use HMMER to create a HMM with this FASTA file.
```
hmmbuild test.hmm test.fasta
```

Now that we created an HMM, we need to determine how to score it. Next create a FASTA file containing minor variations in your sequence to see how it affects the HMM score. For example, here we create a file name *mySeqTests.fasta* and make minor changes which we would still consider valid hits.

``` 
echo ">seq
GAAAATCTCTAGCA
>seq_1SNP
GAAGATCTCTAGCA
>seq_2SNPs
GAAGATCTCAAGCA
>seq_1del
GAAAATTCTAGCA
>seq_1del_1ins
GAAATCTCTGAGCA" > test2.ff
```
Once we create a couple of minor variations in our target sequence, we evaluate the variations with our HMM.  
First run this command to evaluate the test sequences:

```
nhmmer --F1 1 --F2 1 --F3 1 -T -5 --incT -5 --nobias --popen 0.15 --pextend 0.05 --tblout out.tbl test.hmm test2.ff
```

Next, review the output (out.tbl) to determine a minimum acceptable score:
  
```
# target name        accession  query name           accession  hmmfrom hmm to alifrom  ali to envfrom  env to  sq len strand   E-value  score  bias  description of target
#------------------- ---------- -------------------- ---------- ------- ------- ------- ------- ------- ------- ------- ------ --------- ------ ----- ---------------------
seq                  -          test                 -                1      14       1      14       1      14      14    +      0.0063    3.7   1.1  -
seq_1SNP             -          test                 -                1      14       1      14       1      14      14    +       0.016    2.8   0.3  -
seq_2SNPs            -          test                 -                1      13       1      13       1      14      14    +        0.15    0.5   0.9  -
seq_1del_1ins        -          test                 -                3      10       2       9       1      14      14    +        0.29   -0.1   0.2  -
seq_1del             -          test                 -                4      13       3      12       1      13      13    +        0.38   -0.4   1.2  -
```
  
Examine the HMM scores in column 14 (score) and make a decision about the lowest score that provides an acceptable match. In this example, we will go with 0.5. Next we will create an settings file for the new HMM. This file needs to have the same name as the HMM file except we replace ".hmm" with ".cfg". The settings file provides default scoring parameters for the HMM. Here is an example:
  
```
HMMminStartPos  1
HMMmaxStartPos  5
HMMminFullBitScore      10
HMMmaxFullBitScore      30
HMMmatchEnd     TRUE
HMMmatchTerminalSeq     CA
HMMmatchEndRadius       2
```

Once an HMM is created, it should be tested on real data. The `testHMMs` module reads in the output of the demultiplex module and runs demultiplexed reads through their associated HMMs. HMM scores and HMM alignment start positions are plotted on a grid. HMM hits that would be included in an analysis are within the blue box drawn atop of the grid. The position of the blue box is determined by the HMM processing parameters which are passed to the module using the `--HMMparams` flag. This flag accepts a comma delimited string of the processing parameters shown above. Adjust parameters until the blue box is gating scores in an manner appropriate for your work. These settings should be recorded in an HMM cfg file to be used in the future. If you are working with wild infections or anchor reads have the option of landing on multiple locations within LTRs, this configuration should be done for each sequencing run and an experiment HMM parameters should be passed to the `prepReads` module with the same  `--HMMparams` flag.



```
inspiired2 testHMMs --outputDir out --outputDir INSPIIRED2  \
--inputData INSPIIRED2/demultiplex.rds                       \
--HMMparams 'HIV1_LTR_U5_v1.0.hmm,1,6,6,26,TRUE,CA,2'        \
--scoreBinWidth 2 --startPosBinWidth 2 --minScoreBinPct 1
```

<p align="center">
  <img src="figures/testHMMs.png" alt="INSPIIRED HMM test" />
</p>

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

<br>


## Core analysis modules

INSPIIRED2 processes sequencing data through eight core modules. Each module saves an RDS object that becomes the input to the next stage, allowing an analysis to resume from an intermediate result. The first six modules identify integration sites and calculate their abundance; the final two add genomic annotations.

The examples below use `out` as the output directory and the default output names. Set `--fileTag` to change a module's output prefix, and update the next module's `--inputData` accordingly. Each core module also writes a `.log` file, a `.yml` parameter record, and a `.done` marker after successful completion. Run `inspiired2 <module> --help` for its full argument list.

### demultiplex: assign read pairs to sample replicates

`demultiplex` separates a sequencing run into the sample replicates defined in the sample data file. Assignment uses the Index 1 barcode together with the linker sequences on the adrift read. Combining these identifiers helps distinguish libraries that share an index or linker.

**Input:** synchronized Index 1, adrift-read, and anchor-read FASTQ files, plus a tab-delimited sample data file.

**Output:** `demultiplex.rds`, containing assigned read pairs, sample metadata, read sequences, and read counts; `demultiplex.tbl`, a tab-delimited copy of the sample table with `demultiplexedReads` counts.

```bash
inspiired2 demultiplex --outputDir out --threads 30 \
  --sampleData sampleData.tsv \
  --indexReads I1.fastq.gz \
  --adriftReads R1.fastq.gz \
  --anchorReads R2.fastq.gz
```

The module checks that the three FASTQ files have matching read IDs in the same order. It tests whether the Index 1 sequences should be reverse-complemented, trims low-quality read tails, and matches the barcode and the linker segments on either side of the UMI. It then removes the linker from the adrift read and trims matching poly-G tails from both mates. Reads assigned more than once are removed to avoid ambiguous sample assignments.

By default, identical combinations of UMI, anchor sequence, and adrift sequence are collapsed within each sample replicate. The number of original read pairs is retained in `nReads`.

| Option | Default | Effect |
|---|---|---|
| `--index1ReadMaxMismatch` | `1` | Allowed mismatches when matching Index 1. |
| `--adriftReadLinkerMaxMismatch` | `1` | Allowed mismatches in the linker segment before the UMI. |
| `--postUmiLinkerMaxMismatch` | `1` | Allowed mismatches in the linker segment after the UMI. |
| `--qualTrimScore` | `10` | Phred score threshold used for tail trimming. |
| `--qualTrimHalfWidth` | `3` | Half-width, in bases, of the quality-trimming window. |
| `--qualTrimEvents` | `2` | Number of low-quality events within the window required to trigger trimming. |
| `--captureUMIs` | Off | Preserve recovered UMI sequences for downstream processing. |
| `--correctGolayIndexReads` | Off | Apply correction for 12-nucleotide Golay barcodes before assignment. |
| `--disableSequenceCollapse` | Off | Retain individual assigned read pairs instead of collapsing identical records. |

**UMI behavior:** UMI sequences are extracted during demultiplexing, but the default output replaces them with a common placeholder. Enable `--captureUMIs` if UMI information is required. With the default settings, downstream abundance estimates use distinct fragment lengths, and the final UMI count fields are reported as `NA`.

Source: [modules/demultiplex.R](modules/demultiplex.R).

### prepReads: identify the vector boundary and recover genomic sequence

`prepReads` identifies the vector-terminal sequence at the beginning of each anchor read and removes it before genome alignment. This sequence, called the *leader sequence*, is retained separately for later reporting. Recognizing the leader helps establish where the read crosses from vector DNA into flanking genomic DNA.

**Input:** `demultiplex.rds`.

**Output:** `prepReads.rds`, containing the prepared genomic read pairs and recovered `leaderSeq`; `prepReads_vectorHitReads.tsv.gz`, containing reads rejected by the vector filter when that filter is enabled.

```bash
inspiired2 prepReads --outputDir out --threads 30 \
  --inputData out/demultiplex.rds
```

The module runs `nhmmer` with the HMM assigned to each library. It selects the highest-scoring forward-strand hit per anchor read and applies the configured start-position range, minimum and maximum bit scores, and optional HMM-end and terminal-sequence requirements. Accepted leaders are removed from the anchor reads.

For short fragments, a read may continue through the genomic insert and into the sequence at the opposite end. The default over-read trimming step detects this using reverse-complemented linker and leader sequences. It then requires both genomic read segments to meet the minimum length. Finally, a BLAST search tests the tail of each anchor read against its assigned vector sequence; matching pairs are removed as likely internal-vector reads.

| Option | Default | Effect |
|---|---|---|
| `--minReadLength` | `30` | Minimum length of each mate after the enabled over-read trimming step. |
| `--ORtrimPatternWidth` | `8` | Number of bases used to recognize over-reading. |
| `--ORseqMaxMismatch` | `0.10` | Mismatch fraction used to calculate the over-read pattern allowance. |
| `--vectorTestWidth` | `25` | Number of bases at the anchor-read tail tested against the vector. |
| `--vectorTestMinPercentID` | `90` | Minimum sequence identity percentage for a vector hit. |
| `--vectorTestMinCoverage` | `90` | Minimum coverage percentage of the tested sequence for a vector hit. |
| `--disableOverReadTrimming` | Off | Skip over-read trimming and its associated minimum-length filter. |
| `--disableVectorFilter` | Off | Skip the internal-vector read filter. |
| `--HMMparams` | `none` | Supply HMM-specific processing parameters instead of reading the matching `.cfg` files. |

HMM thresholds normally come from `data/hmms/<HMM-name>.cfg`. A `--HMMparams` entry contains the HMM filename followed by seven comma-separated values: minimum start, maximum start, minimum bit score, maximum bit score, match-end flag, terminal sequence, and end radius. Separate multiple entries with `|` and quote the entire argument. When this override is used, provide an entry for every HMM present in the input. The supporting `testHMMs` command can help inspect these settings.

Source: [modules/prepReads.R](modules/prepReads.R).

### alignReads: align both mates to their reference genomes

`alignReads` uses BLAT to align the prepared anchor and adrift sequences to the reference genome assigned to each library. It retains all alignments that pass the configured filters, allowing later stages to evaluate ambiguous mappings using paired-read evidence.

**Input:** `prepReads.rds` and the required `data/referenceGenomes/<refGenome>.2bit` files.

**Output:** `alignReads.rds`, an R list containing separate `anchorReads` and `adriftReads` alignment tables.

```bash
inspiired2 alignReads --outputDir out --threads 30 \
  --inputData out/prepReads.rds
```

Identical sequences are aligned once per reference genome, and their alignments are joined back to the corresponding read records. Anchor reads are aligned first; adrift reads are processed only when their anchor mate has an accepted alignment. Filters require sufficient sequence identity and query coverage and limit insertions in the query and reference. Stored alignment coordinates are converted from BLAT's format to 1-based, inclusive coordinates.

| Option | Default | Effect |
|---|---|---|
| `--minPercentID` | `95` | Minimum alignment identity percentage. |
| `--minAlignmentCoverage` | `95` | Minimum percentage of the query covered by the alignment span. |
| `--blatTileSize` | `11` | Tile size, in bases, used by BLAT. |
| `--blatStepSize` | `5` | Spacing, in bases, between BLAT seed tiles. |
| `--blatRepMatch` | `3000` | Repeat-match threshold used by BLAT. |
| `--blatMaxtNumInsert` | `1` | Maximum number of insertion events in the target. |
| `--blatMaxqNumInsert` | `1` | Maximum number of insertion events in the query. |
| `--blatMaxtBaseInsert` | `1` | Maximum total inserted bases in the target. |
| `--blatMaxqBaseInsert` | `1` | Maximum total inserted bases in the query. |
| `--dataRowChunkSize` | `2500` | Number of unique sequences sent to each alignment worker. |

Source: [modules/alignReads.R](modules/alignReads.R).

### buildFragments: reconstruct candidate genomic fragments

`buildFragments` combines anchor and adrift alignments from the same read pair into candidate physical fragments. Each fragment connects a vector-genome junction to a genomic shearing boundary.

**Input:** `alignReads.rds`.

**Output:** `buildFragments.rds`, a table of candidate fragments with genomic coordinates, strand, sequences, read counts, and sample metadata.

```bash
inspiired2 buildFragments --outputDir out --threads 30 \
  --inputData out/alignReads.rds
```

Mate alignments must lie on the same chromosome, have opposite strands, and define a fragment within the permitted length range. The anchor alignment determines the fragment strand and integration-side boundary. A read pair can produce several candidate fragments when its alignments are ambiguous; these candidates are carried forward for standardization and resolution.

| Option | Default | Effect |
|---|---|---|
| `--minFrgamentLength` | `40` | Minimum accepted fragment length in bases. |
| `--maxFrgamentLength` | `100000` | Maximum accepted fragment length in bases. |
| `--dataRowChunkSize` | `5000` | Number of read IDs processed in each fragment-building batch. |
| `--dbConfigFile` | `none` | Path to the database credential file. |
| `--dbConfigID` | `none` | Name of the credential group within the database configuration file. |
| `--overwriteDBrecords` | Off | Allow replacement of existing database entries for the same sample-replicate keys. |

The spellings `minFrgamentLength` and `maxFrgamentLength` match the current command-line interface.

Supplying both database options enables optional database and Parquet storage. The module then also writes candidate fragment data to checksum-named Parquet files under `/data`. The SQL `fragments` table stores the trial, subject, sample, replicate, reference genome, detection mode, fragment count, and corresponding filename. This preserves the fragment-level evidence for later selection and reanalysis with different downstream settings. The database path requires a writable `/data` location containing the `.inspiired` marker file. The supporting `pullDBrecords` command retrieves selected Parquet records into an RDS file that can be passed to `buildStdFragments`.

Source: [modules/buildFragments.R](modules/buildFragments.R).

### buildStdFragments: standardize boundaries and filter ambiguous fragments

`buildStdFragments` reduces small differences in fragment coordinates, resolves supported multi-mapping reads, and filters fragment patterns consistent with PCR rearrangements. It then collapses the retained evidence into standardized fragment records for site assembly.

**Input:** `buildFragments.rds`, or a compatible fragment table retrieved with `pullDBrecords`.

**Output:** `buildStdFragments.rds`, containing standardized fragment coordinates, supporting read counts and IDs, UMI lists, and representative leader sequences.

```bash
inspiired2 buildStdFragments --outputDir out --threads 30 \
  --inputData out/buildFragments.rds
```

Processing proceeds through the following steps:

1. **Standardize integration positions.** Nearby integration-side coordinates are mapped to locally supported positions using read-count and distance weighting. This is performed across samples and replicates within a trial and subject, while keeping reference genomes, detection modes, chromosomes, and strands separate.
2. **Standardize shearing boundaries.** The opposite fragment boundary is standardized within each sample replicate and integration position. This reduces inflation of fragment counts caused by small coordinate differences.
3. **Resolve multiple candidate fragments.** When all candidates for a read support the same integration position, the shortest candidate is retained. A read mapping to several integration positions can be rescued when exactly one candidate position has uniquely mapped support within the same trial, subject, reference genome, and detection mode.
4. **Summarize unresolved multi-hits.** Remaining ambiguous reads and candidate positions are grouped into connected networks. Within each network, CD-HIT-EST clusters linker-adjacent adrift sequences to estimate shearing diversity. These results are saved separately and are excluded from the main standardized-fragment output.
5. **Filter competing anchor-sequence clusters.** Similar sequences at the beginning of anchor reads are clustered. When a cluster supports multiple integration positions, the module retains a clearly dominant position based on fragment support or read-record support; otherwise, it removes the competing records. This filter is applied within samples by default.
6. **Process UMIs and collapse fragments.** When real UMIs have been retained, the module consolidates UMI variants using their support. It then combines records with the same standardized fragment identity, sums their read counts, chooses a representative leader, and applies the minimum-read requirement.

| Option | Default | Effect |
|---|---|---|
| `--intSite_sp_window` | `8` | Search radius, in bases, for integration-position standardization. |
| `--breakPoint_sp_window` | `5` | Search radius, in bases, for shearing-boundary standardization. |
| `--disableIntSitePosStd` | Off | Disable integration-position standardization. |
| `--disableBreakPointPosStd` | Off | Disable shearing-boundary standardization. |
| `--anchorReadClusterLen` | `30` | Number of genomic anchor-read bases used for sequence clustering. |
| `--anchorReadClusterGrouping` | `sample` | Apply the anchor filter within `sample`, `subject`, or `trial` groups. |
| `--anchorReadClusterMinAbundDiff` | `5` | Minimum fragment-count advantage over the next candidate for selection. |
| `--anchorReadClusterMinReadMult` | `10` | Alternative minimum ratio of supporting read records for selection. |
| `--disableAnchorReadClusteringFilter` | Off | Skip the competing anchor-sequence filter. |
| `--multiHitclusteringNTlen` | `30` | Number of linker-adjacent adrift bases used to cluster unresolved multi-hits. |
| `--minReadsPerFrag` | `1` | Minimum summed read support for a retained fragment. |

The anchor-cluster read-support comparison counts input read records; it does not sum their `nReads` values. Final fragment read totals do sum `nReads`.

Additional outputs make the filtering decisions available for inspection:

| File | Contents |
|---|---|
| `buildStdFragments_multiHitFrags.rds` | Candidate fragment records for unresolved multi-mapping reads. |
| `buildStdFragments_multiHitClusters.rds` | Multi-hit networks, candidate positions, and sequence-based shearing-diversity summaries. |
| `buildStdFragments_anchorReadClusters.rds` | Decisions for anchor-sequence clusters with competing positions, when that filter is enabled. |
| `buildStdFragments_multiHitClusterAssignments.tsv.gz` | Per-read CD-HIT assignments, written only with `--saveMultiHitClusteringDetails`. |

The current implementation requires at least one uniquely mapped integration position before attempting multi-hit rescue. The `inspiired2 buildStdFragments` command defaults to `30` threads; the other core commands default to `50`.

Sources: [modules/buildStdFragments.R](modules/buildStdFragments.R) and [lib/buildStdFragments.R](lib/buildStdFragments.R).

### buildSites: assemble integration sites and calculate abundance

`buildSites` combines standardized fragments into sample-level integration-site records. It summarizes evidence across technical replicates and uses distinct fragment lengths as a measure of independent shearing events supporting each site.

**Input:** `buildStdFragments.rds`.

**Output:** `buildSites.rds`, containing site coordinates, detection modes, abundance measures, replicate-level summaries, and representative leader sequences.

```bash
inspiired2 buildSites --outputDir out --threads 30 \
  --inputData out/buildStdFragments.rds
```

When both U3 and U5 evidence are present, the module searches for nearby, opposite-strand detections within the same trial, subject, sample, and reference genome. It combines a pair only when each site has exactly one candidate partner. Ambiguous pairings remain separate. Accepted pairs are labeled `dual detect`, and both representative leader sequences are retained.

The default processing also corrects reported positions for the genomic duplication associated with integration and reports vector orientation in `posid`. U3 strand signs are reversed relative to the original anchor alignment; U5 strand signs are retained. The coordinate correction is controlled by `--integraseCorrectionDist` and should match the integration system being analyzed.

| Option | Default | Effect |
|---|---|---|
| `--dualDetectWidth` | `6` | Maximum distance, in bases, when searching for U3/U5 partners. |
| `--integraseCorrectionDist` | `2` | Coordinate shift applied according to fragment strand during integration-site correction. |
| `--sumSonicBreaksWithin` | `replicates` | Count distinct fragment lengths within each replicate and sum them; use `samples` to count distinct lengths across the entire sample. |
| `--disableDualDetect` | Off | Keep U3 and U5 detections separate. |
| `--disableOrientationCorrection` | Off | Skip the orientation and coordinate correction step for unmerged detections. Accepted dual detections still receive their own correction. |

Key output fields are:

| Field | Meaning |
|---|---|
| `posid` | Chromosome, strand, and integration coordinate combined into one identifier, such as `chr1+123456`. Interpret it together with `refGenome`. |
| `sonicLengths` | Distinct standardized fragment lengths counted using the selected replicate- or sample-level rule. |
| `reads` | Total supporting sequencing read pairs, including counts carried forward from collapsed duplicates. |
| `UMIs` | Number of distinct retained UMIs; `NA` when UMIs were not captured. |
| `nRepsObs` | Number of replicates supporting the site; reported as `NA` for dual detections. |
| `percentSampleRelAbund` | The site's percentage of total `sonicLengths` within its trial, subject, sample, and reference genome, rounded to two decimal places. |
| `repLeaderSeq` | Representative recovered leader sequence, or paired U3/U5 leaders for a dual detection. |
| `repLeaderSeqClusters` | Number of sequence clusters among the fragment-level representative leaders. |
| `rep<N>-...` | Replicate-specific UMI counts, fragment-length counts, read counts, and representative leaders. |

For example, a site supported by two distinct lengths in replicate 1 and three in replicate 2 has `sonicLengths = 5` with the default setting, even if some lengths occur in both replicates. With `--sumSonicBreaksWithin samples`, lengths shared across replicates are counted once. These values describe recovered integration-site evidence; `percentSampleRelAbund` is not a direct measurement of the percentage of cells carrying the integration.

Source: [modules/buildSites.R](modules/buildSites.R).

### nearestGenes: annotate gene and exon context

`nearestGenes` determines whether each integration site overlaps a gene transcription unit or exon and identifies the nearest annotated gene. These annotations help interpret where integrations fall relative to known genes.

**Input:** `buildSites.rds`, plus the matching `<refGenome>.TUs.rds` and `<refGenome>.exons.rds` files in `data/genomeAnnotations`.

**Output:** `nearestGenes.rds`, preserving the site records and adding gene-context columns.

```bash
inspiired2 nearestGenes --outputDir out --threads 30 \
  --inputData out/buildSites.rds
```

Sites are annotated separately for each reference genome. Overlap and nearest-gene searches ignore strand, and tied nearest genes are retained as comma-separated values.

| Added field | Meaning |
|---|---|
| `inGene` | Whether the integration coordinate overlaps an annotated transcription unit. |
| `inExon` | Whether it overlaps an annotated exon. |
| `nearestGene` | Name or names of the nearest annotated genes. |
| `nearestGeneDist` | Distance in bases to the nearest transcription-unit interval; `0` for a site inside a gene. |
| `nearestGeneStrand` | Strand or strands of the reported nearest genes. |
| `beforeNearestGene` | Whether the site has a lower genomic coordinate than the start of the nearest gene interval; for ties, the smallest interval start is used. |

`nearestGeneDist` measures distance to the gene interval, not specifically to its transcription start site. `beforeNearestGene` describes genomic coordinate order, so it should not be interpreted as transcriptional upstream/downstream without considering gene strand. Sites without a nearest annotation on their reference sequence retain missing nearest-gene values. This module has no additional analysis-specific command-line options beyond the shared module options.

Source: [modules/nearestGenes.R](modules/nearestGenes.R).

### annotateRepeats: annotate overlapping repetitive elements

`annotateRepeats` adds repeat annotations to integration sites using the precomputed RepeatMasker-derived table for each reference genome. It reports repeats overlapping the integration coordinate while retaining the existing site and gene annotations.

**Input:** `nearestGenes.rds` in the standard workflow, plus `data/genomeAnnotations/<refGenome>.repeatTable.gz`. A site table from `buildSites` can also be used when gene annotations are not needed.

**Output:** `annotateRepeats.rds`, the final annotated site table in the standard workflow.

```bash
inspiired2 annotateRepeats --outputDir out --threads 30 \
  --inputData out/nearestGenes.rds
```

The module tests overlap at the integration coordinate, independently of strand, and adds `repeat_name` and `repeat_class`. Multiple overlapping repeat annotations are combined into comma-separated values; sites without an overlapping repeat receive `NA`. The input row count is preserved. This stage uses existing repeat annotations and does not run RepeatMasker during the analysis.

There are no additional analysis-specific command-line options beyond the shared module options.

Source: [modules/annotateRepeats.R](modules/annotateRepeats.R).

### Supporting commands

The following commands support setup, inspection, or reuse of results:

| Command | Purpose |
|---|---|
| [`showResources`](modules/showResources.R) | List available reference genomes, annotations, HMMs, and vector sequences, including resource overlays. |
| [`testHMMs`](modules/testHMMs.R) | Inspect HMM scores and read-position distributions using demultiplexed reads; produces an HMM diagnostic PDF. |
| [`buildSeqDataMap`](modules/buildSeqDataMap.R) | Create a PNG showing sequence composition across sorted and binned FASTQ reads. Set `--fileTag` explicitly because its current default is `testHMMs`. |
| [`testDBconn`](modules/testDBconn.R) | Check database connectivity using the supplied credential file and configuration group. |
| [`pullDBrecords`](modules/pullDBrecords.R) | Select stored fragments by trial and optional subject, sample, reference genome, and mode filters, read their Parquet files, and save an RDS input for downstream reanalysis. |

<br>




INSPIIRED2 is provided with an SQL database and the ability to create a data warehouse to store data from multiple experiments. Databasing and warehousing is enabled by providing database credentials to the buildFragments module arguments: `--dbConfigFile --dbConfigID`

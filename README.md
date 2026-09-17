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
<br> 

### Quick start and validation
Download and load the distributed Docker image:

```bash
wget https://bushmanlab.org/export/inspiired2_latest.tar.gz
```
```bash
docker load -i inspiired2_latest.tar.gz
```
Run the bundled synthetic test:

```bash
docker run --rm -it --shm-size=5g inspiired2 bash
```

```bash
cd /opt/INSPIIRED2/tests/synTests/U5_50sites_seed1
```
```bash
./run.sh
```

The test takes about 5 minutes to complete and report PASS at the end will if the MD5sum of the output matches the expected value. 

### Starting analyses

The basic INSPIIRED2 invocation command has this structure:
```bash
docker run --rm     \
  --shm-size=20g    \
  -v ./:/workspace  \
  -w /workspace     \
  inspiired2 bash run.sh
```
The`--shm-size` flag defines the max. amount of memory allowed to be used as scratch space during analysis. 20GB is a reasonable value for most moderate size Illumina paied-end data sets. This value should be increased for large data sets and should not reach an appreciable percentage of your total RAM. 

` -v ./:/workspace` mounts your analysis directory to `/workspace` inside of the Docker container. Here we are mounting the current directory `./`. The analysis directory is expected to contain your sequencing data, sample data file, and processing script (described next).

`-w /workspace` instructs Docker to make all paths relative to `/workspace` within the Docker image.

`inspiired2 bash run.sh` instructs docker to run the processing script `run.sh`, located in your analysis directory, in a Docker container created with the `inspiired2` Docker image. 

By default, all output files will be owned by root. To change ownership to the user initiating the analysis, add this argument:  `--user "$(id -u):$(id -g)"`

INSPIIRED2 is provided with a number of reference genomes (hg38, hs1, sacCer3, mm10, canFam4, and macFas5) as well as U3 and U5 LTR HMMs created with data from Los Alamos National laboratories. Custom reference genomes, gene annotations, vector sequences, and HMMs can be shared with the Docker image at run time by using an additional mount flag: `-v  /usr/local/customData:/resources:ro`
Custom data must be organized in the same way that data is organized withing INSPIIRED's data folder (below). Importantly, Docker's `--user` flag can not be used when overlaying local resources onto the Docker container. Doing so will cause an error. 

```
customData/
├── hmms/
│   ├── myVectorEnd.hmm
│   └── myVectorEnd.cfg
├── vectors/
│   └── myVector.fasta
├── referenceGenomes/
│   └── myAssembly.2bit
└── genomeAnnotations/
    ├── myAssembly.TUs.rds
    ├── myAssembly.exons.rds
    └── myAssembly.repeatTable.gz
```



INSPIIRED2 is provided with an SQL database and the ability to create a data warehouse to store data from multiple experiments. Databasing and warehousing is enabled by providing database credentials to the buildFragments module arguments: `--dbConfigFile --dbConfigID`

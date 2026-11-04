# Processing provenance

INSPIIRED2 carries processing history through the data-producing module chain. RDS
outputs store it in the `inspiired2.provenance` attribute. Fragment and site Parquets
store a portable YAML bundle under the same name in file metadata. Multihit RDS
archives also carry the attribute. No SQL schema change or additional data columns
are required.

Zero-cluster database records retain the existing `NULL` archive filename. Their
history is available in the local `_multiHitClusters.rds` output, but there is no
database archive file in which to embed it.

Each module inherits the history from its actual input and appends its own execution
record. `buildStdFragments --pullSubjectFragments` also reads and combines the
histories embedded in the archived fragments it introduces. Source records retain
the fragment database keys and archive file identity, so older U5 and newer U3
histories remain distinguishable. The complete inherited history is embedded in the
new output; it does not depend on an older archive remaining available.

The bundle contains:

| Field | Meaning |
| --- | --- |
| `version` | Provenance format version, currently 1. |
| `heads` | IDs of the records that produced or introduced this artifact. |
| `records` | Execution and source records, indexed by immutable IDs. Each record links to its `parents`. |
| `artifact` | Description of this output, including row count and available sample keys. |

Execution records contain `module`, `recorded_at`, `software`, `arguments_yaml`,
`details`, `resources`, and `input_files`. Software identity includes the reported
INSPIIRED2 version, R version, and checksums of the available relevant code files.
The argument snapshot includes parsed defaults; credential fields, connection
objects, and the raw command string are omitted. It is a sanitized snapshot, not a
byte-for-byte copy of the module's original `.yml`. Credential-file contents are
never embedded.

The module details supplement the arguments with decisions made during processing.
For example, demultiplex records canonical sample definitions and the selected I1
orientation. prepReads records the effective HMM parameter table, including values
resolved from configuration files. Used HMM/configuration/vector resources and
alignment reference files are identified by checksums. FASTQ and input RDS files
are described by path, size, and modification time; those descriptors are **not**
content checksums. Resource checksumming adds file-reading work, particularly for
large reference genomes.

Source records distinguish `history: embedded` from `history: unavailable`. Missing
history in an older input remains explicit; the current run's parameters are never
assigned to that older data. Malformed or unsupported embedded history is an error.
The separate `inspiired2.siteRebuild` manifest and `pullSubjectFragments` database
authorization checks retain their existing meaning. Processing provenance does not
make an otherwise ineligible legacy object eligible for a site database update.

Only actual data ancestors are included. `validateSampleData` and `testHMMs` are
independent checks, not implicit ancestors of prepReads. The software does not scan
an output directory and assume that every YAML found there belongs to the input.

## Exporting a readable history

Run the helper from the installation directory, or provide its full path:

```bash
Rscript tools/exportProvenance.R \
  --inputData INSPIIRED2_U5/annotateRepeats.rds \
  --outputFile INSPIIRED2_U5/annotateRepeats_provenance.yml

Rscript tools/exportProvenance.R \
  --inputData /data/ACTUAL_FILENAME.parquet \
  --outputFile archived_provenance.yml
```

Replace `ACTUAL_FILENAME.parquet` with the `data_file_name` from the corresponding
SQL record. The utility infers the installation root from its own location;
`--softwareRoot` can override it. RDS inspection loads the object. Parquet inspection
reads its footer without loading the data rows. Exporting does not alter the input
or append a new processing record. An old file without embedded history produces
an explicit unavailable-history record.

## Scope and retention

This feature records processing history; it does not enforce matching parameters
between runs, create multiple analysis versions, or change the database replacement
keys. Two database-enabled buildSites executions for the same scope still replace
one another, even if their `--fileTag` values differ. Run exploratory alternatives
without database credentials when the primary database result should be retained.

Metadata is embedded before the final archive checksum and filename are calculated.
Consequently, changed provenance changes the archive checksum even when data values
are identical. Unique execution IDs and timestamps also change RDS checksums between
runs. The synthetic regression check removes only processing provenance from a
temporary copy before comparing the data checksum; the original output retains its
history.

Replacing an archive retains the new result's input history, not a recoverable copy
of every previous result. Existing obsolete-file cleanup remains in effect. Full
module outputs should continue directly through the pipeline; arbitrary external
table transformations can discard custom attributes.

## U5/U3 limit-of-detection acceptance run

Use the single-clone dilution experiment to check both the biological result and
the recorded history. Define the expected integration locus or loci and the allowed
coordinate tolerance before evaluating recovery.

1. Use a distinct sample ID for each dilution. U5 and U3 measurements of the same
   biological sample must share the normalized trial, subject, sample, and reference
   genome identifiers. With the current `--sumSonicBreaksWithin replicates` behavior,
   use distinct assay replicate numbers, for example U5 1–4 and U3 5–8, to avoid
   pooling same-numbered libraries during dual detection.
2. Process U5 through the normal chain, using database credentials and
   `--pullSubjectFragments` at buildStdFragments, including the initial run with no
   prior history. Complete site building and annotations. Export the fragment and
   final site provenance from both their RDS outputs and actual archived files.
   Check the recorded HMM settings, sample keys, and software identity.
3. Process U3 with the same biological identifiers and again use
   `--pullSubjectFragments`. Complete the chain before another run for that
   trial/subject. Export the final site history. It should include both original U5
   and new U3 fragment histories, their respective HMM settings, the archive source
   identities for pulled fragments, and the new downstream processing records.
4. Check clone recovery and abundance across dilutions. Sites with suitable U5 and
   U3 support should become `dual detect`; unpaired evidence can remain U5-only or
   U3-only. Low-input stochastic dropout is not itself a provenance failure. Check
   that unrelated sample records remain intact after each replacement.
5. Run an optional comparison buildSites with `--disableDualDetect`,
   `--disableOrientationCorrection`, and a distinct `--fileTag`, **without database
   credentials**. Its history should show those flags while the primary database
   result remains unchanged.

Globally empty stages retain their existing error exits: a dilution run with no
surviving fragments or sites may stop instead of producing a successful zero-result
report. Review its logs to distinguish this behavior from a provenance error.

Automated provenance tests exercise metadata and lineage independently of the full
sequencing pipeline. A real Arrow round trip is required for the archive checks;
these tests do not substitute for this experiment or a live database run.

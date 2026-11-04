# Processing provenance tests

From the INSPIIRED2 installation root:

```bash
Rscript tests/unitTests/provenance/tests.R .
Rscript tests/unitTests/provenance/buildFragments.R .
Rscript tests/unitTests/provenance/downstream.R .
Rscript tests/unitTests/buildFragmentsDB/tests.R modules/buildFragments.R
```

The provenance tests use real RDS and Arrow Parquet files. They check data and
attribute preservation, complete U5/U3 ancestry, legacy inputs, malformed
metadata, scoped archive imports, parameter capture, resource checksums, export,
and reproducibility after removing only provenance. The main suite executes the
real `pullDBfragments` function with a database query double and runs the actual
`buildSites` body for default and alternative settings. The other scripts run
fragment construction and annotation fixtures. External services and runtime
startup are substituted where indicated in each script.

For an exact comparison against an unmodified checkout, pass its directory as
the second argument to the provenance scripts:

```bash
Rscript tests/unitTests/provenance/tests.R . /path/to/unmodified/INSPIIRED2
Rscript tests/unitTests/provenance/buildFragments.R . /path/to/unmodified/INSPIIRED2
Rscript tests/unitTests/provenance/downstream.R . /path/to/unmodified/INSPIIRED2
```

The database suite injects SQL and serialization failures while exercising real
file operations. Its Parquet/provenance writer is a double; the provenance suite
separately tests actual Arrow metadata. These tests do not replace a live
MariaDB end-to-end U5/U3 sequencing run.

The synthetic sequencing test still uses its existing `expected_md5sum`. Its
runner validates provenance, removes only the new provenance attribute in a
temporary copy, and calculates the scientific-object checksum. It does not
modify the archived output. The unit suite checks that this operation reproduces
the original compressed RDS bytes for a representative data.table with factors
and site-rebuild attributes. The complete FASTQ synthetic test still needs the
normal INSPIIRED2 environment and its supplied input files.

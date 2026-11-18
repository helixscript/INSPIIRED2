#!/usr/bin/env python3
"""Package fixed hg38 synthetic datasets as INSPIIRED1 run folders.

Run from the directory containing hg38_U5_5Ksites_seed*_*error:
    python3 package_syn_illumina_runs.py --dataset-set 2pct --dry-run
    python3 package_syn_illumina_runs.py --dataset-set 2pct

Dataset sets:
    original (default): the original 12 runs, 0001-0012 / GTSP9050-GTSP9061.
    2pct: only the new 2% runs, 0013-0015 / GTSP9062-GTSP9064 (seeds 1-3).
    all: all 15 runs; every destination must be new.
The original run mapping is unchanged, including the 5% conditions.
The new 2% inputs must contain generation.json recording the expected seed,
error rate, and zero integrase correction. Original inputs are also checked
for seed/error agreement when generation.json is available.

Alternatively specify --source-dir and, optionally, --output-dir.
Only Python's standard library is required (Python 3.8+).
Source datasets are copied, never moved or modified. Existing run folders
are never overwritten. FASTQ bytes and read headers are preserved exactly.
"""

import argparse
import csv
import json
import math
import re
import shutil
import sys
import tempfile
from pathlib import Path


SAMPLE_SHEET_HEADER = '''[Header],,,,,,,,
Experiment Name,"SynDataSet",,,,,,,
Date,10/1/26,,,,,,,
Module,GenerateFASTQ - 2.0.1,,,,,,,
Workflow,GenerateFASTQ,,,,,,,
Library Prep Kit,Custom,,,,,,,
Chemistry,Default,,,,,,,
[Reads],,,,,,,,
171,,,,,,,,
135,,,,,,,,
[Data],,,,,,,,
Sample_ID,Description,I7_Index_ID,index,Sample_Project,,,,
GeneTherapy-10012026,,TTTTTTTTTTTT,TTTTTTTTTTTT,,,,,
,,,,,,,,
[metaData],,,,,,,,
alias,linkerSequence,bcSeq,gender,primer,ltrBit,largeLTRFrag,vectorSeq,refGenome
'''

FASTQ_NAMES = {
    read + '.fastq.gz': 'Undetermined_S0_' + read + '_001.fastq.gz'
    for read in ('I1', 'R1', 'R2')
}
BASECALLS = Path('Data') / 'Intensities' / 'BaseCalls'

# Keep the original ordering fixed. New conditions are appended, never inserted.
ORIGINAL_RUNS = tuple(
    (run_number, seed, error)
    for run_number, (seed, error) in enumerate(
        ((seed, error) for seed in (1, 2, 3)
         for error in ('0.00', '0.01', '0.03', '0.05')), start=1)
)
TWO_PERCENT_RUNS = ((13, 1, '0.02'), (14, 2, '0.02'), (15, 3, '0.02'))
RUN_SETS = {
    'original': ORIGINAL_RUNS,
    '2pct': TWO_PERCENT_RUNS,
    'all': ORIGINAL_RUNS + TWO_PERCENT_RUNS,
}


def read_metadata(path):
    """Require the single U5/hg38 replicate described in the experiment plan."""
    required = {'replicate', 'index1Seq', 'adriftReadLinkerSeq', 'refGenome', 'mode'}
    with path.open('r', encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames or []
        missing = required - set(fields)
        if missing:
            raise ValueError('{}: missing columns: {}'.format(path, ', '.join(sorted(missing))))
        if len(fields) != len(set(fields)):
            raise ValueError('{}: duplicate column names'.format(path))
        rows = [row for row in reader if any(value not in ('', None) for value in row.values())]
    if len(rows) != 1:
        raise ValueError('{}: expected one sample row, found {}'.format(path, len(rows)))
    row = rows[0]
    if None in row or any(value is None for value in row.values()):
        raise ValueError('{}: malformed TSV row'.format(path))
    row = {key: value.strip() for key, value in row.items()}
    for field, expected in [('replicate', '1'), ('mode', 'U5'), ('refGenome', 'hg38')]:
        if row[field] != expected:
            raise ValueError('{}: {} must be {!r}, found {!r}'.format(path, field, expected, row[field]))
    if not re.fullmatch(r'[ACGTacgt]{12}', row['index1Seq']):
        raise ValueError('{}: index1Seq must contain 12 A/C/G/T bases'.format(path))
    if not re.fullmatch(r'[ACGTNacgtn]+', row['adriftReadLinkerSeq']):
        raise ValueError('{}: adriftReadLinkerSeq must contain A/C/G/T/N bases'.format(path))
    return row


def validate_generation(source, seed, error):
    """Reject a mislabeled dataset before copying; allow older original inputs."""
    path = source / 'generation.json'
    if not path.is_file():
        if error == '0.02':
            raise ValueError('Missing provenance for new 2% dataset: {}'.format(path))
        return
    with path.open('r', encoding='utf-8') as handle:
        provenance = json.load(handle)
    parameters = provenance.get('parameters') if isinstance(provenance, dict) else None
    if not isinstance(parameters, dict):
        raise ValueError('{}: expected a parameters object'.format(path))
    actual_seed = parameters.get('seed')
    if isinstance(actual_seed, bool) or actual_seed != seed:
        raise ValueError('{}: expected seed {}, found {!r}; check dataset directory names'.format(
            path, seed, actual_seed))
    actual_error = parameters.get('percentGenomicError')
    if (isinstance(actual_error, bool) or not isinstance(actual_error, (int, float))
            or not math.isclose(actual_error, float(error), rel_tol=0, abs_tol=1e-12)):
        raise ValueError('{}: expected percentGenomicError {}, found {!r}'.format(
            path, error, actual_error))
    if error == '0.02':
        for label, actual in (
                ('integraseCorrectionDist', provenance.get('integraseCorrectionDist')),
                ('parameters.integraseCorrectionDist', parameters.get('integraseCorrectionDist'))):
            if isinstance(actual, bool) or actual != 0:
                raise ValueError('{}: {} must be 0 for the new benchmark datasets, found {!r}'.format(
                    path, label, actual))


def plan_runs(source_root, output_root, dataset_set='original'):
    """Validate all inputs and destinations before copying any dataset."""
    jobs = []
    for run_number, seed, error in RUN_SETS[dataset_set]:
        dataset = 'hg38_U5_5Ksites_seed{}_{}error'.format(seed, error)
        run_id = '200821_MN01490_{:04d}_A000HCJNY3'.format(run_number)
        gtsp = 'GTSP{}'.format(9049 + run_number)
        source = source_root / dataset
        destination = output_root / run_id
        if not source.is_dir():
            raise ValueError('Missing source directory: {}'.format(source))
        resolved_source = source.resolve()
        if output_root == resolved_source or resolved_source in output_root.parents:
            raise ValueError('Output directory must not be inside a source dataset: {}'.format(source))
        if destination.exists() or destination.is_symlink():
            raise ValueError('Run folder already exists; nothing will be overwritten: {}'.format(destination))
        for name in list(FASTQ_NAMES) + ['sampleData.tsv', 'truth.tsv']:
            if not (source / name).is_file():
                raise ValueError('Missing source file: {}'.format(source / name))
        for renamed in FASTQ_NAMES.values():
            target = source / BASECALLS / renamed
            if target.exists() or target.is_symlink():
                raise ValueError('Source already contains a packaged FASTQ: {}'.format(target))
        metadata = read_metadata(source / 'sampleData.tsv')
        validate_generation(source, seed, error)
        jobs.append((source, destination, gtsp, metadata))
    return jobs


def package_runs(jobs, output_root):
    output_root.mkdir(parents=True, exist_ok=True)
    # Finish all copies in a temporary directory before publishing run folders.
    with tempfile.TemporaryDirectory(prefix='.faux-illumina-', dir=str(output_root)) as temporary:
        staging_root = Path(temporary)
        for source, destination, gtsp, metadata in jobs:
            staged = staging_root / destination.name
            shutil.copytree(str(source), str(staged))
            basecalls = staged / BASECALLS
            basecalls.mkdir(parents=True, exist_ok=True)
            for original, renamed in FASTQ_NAMES.items():
                (staged / original).rename(basecalls / renamed)
            with (staged / 'SampleSheet.csv').open('w', encoding='utf-8', newline='') as handle:
                handle.write(SAMPLE_SHEET_HEADER)
                writer = csv.writer(handle, lineterminator='\n')
                writer.writerow([
                    gtsp + '-1', metadata['adriftReadLinkerSeq'], metadata['index1Seq'],
                    'm', 'GAAAATC', 'TCTAGCA', 'TGCTAGAGATTTTCCACACTGACTAAAAGGGTCT',
                    'NONE.fa', metadata['refGenome'],
                ])
        for source, destination, gtsp, metadata in jobs:
            if destination.exists() or destination.is_symlink():
                raise ValueError('Destination appeared during copying; refusing to replace: {}'.format(destination))
            (staging_root / destination.name).rename(destination)
            print('Created {}  [{}-1; {}]'.format(destination, gtsp, source.name))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--source-dir', type=Path, default=Path('.'),
                        help='Directory containing the selected hg38_* datasets (default: current directory).')
    parser.add_argument('--output-dir', type=Path,
                        help='Parent for the run folders (default: same as --source-dir).')
    parser.add_argument('--dry-run', action='store_true',
                        help='Validate inputs and display the mapping without creating files.')
    parser.add_argument('--dataset-set', choices=tuple(RUN_SETS), default='original',
                        help='original: runs 1-12 (default); 2pct: only runs 13-15; all: all 15 runs.')
    args = parser.parse_args(argv)
    source_root = args.source_dir.expanduser().resolve()
    output_root = args.output_dir.expanduser().resolve() if args.output_dir else source_root
    try:
        jobs = plan_runs(source_root, output_root, args.dataset_set)
        if args.dry_run:
            for source, destination, gtsp, metadata in jobs:
                print('{} -> {}  [{}-1; barcode {}]'.format(
                    source.name, destination, gtsp, metadata['index1Seq']))
            print('Validated {} datasets; no files created.'.format(len(jobs)))
        else:
            package_runs(jobs, output_root)
            print('Done: created {} faux Illumina runs. Source datasets are unchanged.'.format(len(jobs)))
    except (OSError, ValueError, csv.Error) as error:
        print('ERROR: {}'.format(error), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())

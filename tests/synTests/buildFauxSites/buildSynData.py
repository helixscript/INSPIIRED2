#!/usr/bin/env python3
"""Generate reproducible integration-site reads, replicate metadata, and truth."""

import argparse
from array import array
from bisect import bisect_left, bisect_right
from contextlib import ExitStack
import csv
from dataclasses import dataclass
import gzip
import hashlib
import io
import math
import os
from pathlib import Path
import random
import re
import sys
import tempfile

REMNANT0 = 'GAAAATCTCTAGCA'
REMNANT = 'TCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTG'
LINKER = 'GAACGAGCACTAGTAAGCCCNNNNNNNNNNNNCTCCGCTTAAGGGACT'
DNA = 'ATCG'
COMPLEMENT = str.maketrans('ACGT', 'TGCA')
BLOCK_LENGTH = 1000
START_POSITION = 10
BASE_FRAGMENT_END = 500
CHROMOSOME_MARGIN = 10000
# Fixed context preserves both Gaussian tails without making site selection
# depend on the requested chatter SD (at most 5 nt).
CONTEXT_PADDING = 64
MAX_SITE_ATTEMPTS = 10000
OUTPUT_NAMES = ('R1.fastq.gz', 'R2.fastq.gz', 'I1.fastq.gz',
                'sampleData.tsv', 'truth.tsv')
REPEAT_TRUTH_FIELDS = (
    'repeatGroup', 'repeatClass', 'repeatName', 'repeatStart', 'repeatEnd',
    'repeatLength', 'repeatStrand', 'repeatTableLine', 'repeatDistanceToStart',
    'repeatDistanceToEnd', 'repeatBoundaryDistance', 'nR1GenomicWithinRepeat',
    'nR2GenomicWithinRepeat', 'nReadPairsFullyWithinRepeat')


class ExplicitSiteCount(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        namespace.nSites = values
        namespace.nSitesExplicit = True


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    options = [
        ('-o', '--outputDir', str, './output', 'Output directory.'),
        ('-m', '--mode', str, 'integrase', 'Mode: integrase or AAV.'),
        ('-n', '--nSites', int, 1000, 'Total number of integration sites.'),
        ('-f', '--nFrags', int, 5, 'Fragments per site across all replicates (1-100).'),
        ('-p', '--nReadsPerFrag', int, 25, 'Read pairs per fragment.'),
        ('-d', '--distBetweenFragEnds', int, 25, 'Distance between true fragment ends (nt).'),
        ('-s', '--seed', int, 1, 'Random seed.'),
        ('-r', '--refGenomePath', str, '../../../data/referenceGenomes/hg38.2bit', 'Reference 2bit file.'),
        ('-g', '--refGenomeID', str, 'hg38', 'Reference ID for sample metadata.'),
        ('-i', '--integraseHMM', str, 'validation.hmm', 'Leader HMM (integrase only).'),
        ('-t', '--anchorReadStartSeq', str, 'TCTGCGCGCT', 'Anchor filter metadata (AAV only).'),
        ('-x', '--R1_length', int, 150, 'Total R1 length, including linker and UMI.'),
        ('-y', '--R2_length', int, 150, 'Total R2 length, including leader.'),
        ('-z', '--I1_length', int, 12, 'Index barcode and I1 read length.'),
        ('-e', '--percentGenomicError', float, 0, 'Substitution fraction for R1 genomic bases and all R2 bases, including the remnant (0-1).'),
        ('-c', '--positionChatterSD', float, 0.50, 'Gaussian endpoint chatter SD (0-5 nt).'),
        ('-b', '--minSiteDistance', int, 500, 'Minimum distance between reported sites on a chromosome.'),
        ('-v', '--nReplicates', int, 4, 'Technical replicates per sample.'),
    ]
    for short, long, kind, default, help_text in options:
        extra = {'action': ExplicitSiteCount} if long == '--nSites' else {}
        parser.add_argument(short, long, type=kind, default=default,
                            help=f'{help_text} Default: {default}.', metavar='', **extra)
    parser.set_defaults(nSitesExplicit=False)
    parser.add_argument('-k', '--singleSample', action='store_true',
                        help='Consolidate all sites into one sample.')
    parser.add_argument('--illuminaReadIDs', action='store_true',
                        help='Use Illumina-formatted FASTQ headers with unique seeded '
                             'cluster coordinates and the I1 barcode. Default: descriptive IDs.')
    sources = parser.add_mutually_exclusive_group()
    sources.add_argument('--posidFile', help='Text file with one exact chr+position or chr-position per line; its row count defines nSites.')
    sources.add_argument('--repeatTable', help='RepeatMasker TSV, plain or gzip: balance sites across repeat groups and sample their center thirds.')
    parser.add_argument('--repeatClasses', nargs='+', help='Selected groups, separated by spaces or commas; e.g. SINE/Alu Satellite Simple_repeat. Alu is an alias for SINE/Alu.')
    parser.add_argument('--repeatGrouping', choices=('class', 'family'), default='class',
                        help='Without repeatClasses, balance by broad class or exact class/family label. Default: class.')
    parser.add_argument('--includeRepeatAlternates', action='store_true',
                        help='Include RepeatMasker rows marked * (overlapping, lower-scoring matches).')
    return parser


def make_rng(seed, domain, *indices):
    """Stable independent streams; never use Python's process-randomized hash()."""
    key = repr((seed, domain, *indices)).encode('utf-8')
    return random.Random(int.from_bytes(hashlib.sha256(key).digest(), 'big'))


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def chatter_radius(sd):
    # Truncate symmetrically beyond eight SDs, never clip draws onto a boundary.
    return math.ceil(8 * sd)


def gaussian_offset(rng, sd):
    if sd == 0:
        return 0
    radius = chatter_radius(sd)
    while True:
        offset = round(rng.gauss(0, sd))
        if abs(offset) <= radius:
            return offset


def unique_dna(rng, length, count):
    """Sample uniformly without replacement, including nearly full spaces.

    A sparse Fisher-Yates shuffle avoids collision retries and does not need to
    materialize 4**length sequences. Prefixes do not depend on the final count.
    """
    remaining = 4 ** length
    if count > remaining:
        raise ValueError(f'Cannot make {count} unique {length}-base sequences.')
    swaps = {}
    for _ in range(count):
        index = rng.randrange(remaining)
        value = swaps.get(index, index)
        remaining -= 1
        swaps[index] = swaps.get(remaining, remaining)
        swaps.pop(remaining, None)
        bases = []
        for _ in range(length):
            value, base = divmod(value, 4)
            bases.append(DNA[base])
        yield ''.join(reversed(bases))


def build_remnant_library(seed, mode):
    if mode == 'integrase':
        return [REMNANT0]
    rng = make_rng(seed, 'remnant-library')
    pieces = []
    for _ in range(500):
        start = rng.randint(12, len(REMNANT) - 12)
        width = rng.randint(12, min(36, len(REMNANT) - start))
        pieces.append(REMNANT[start:start + width])
    remnants = []
    for _ in range(500):
        prefix = REMNANT[:rng.randint(12, 24)]
        parts = [piece if rng.randrange(2) == 0 else reverse_complement(piece)
                 for piece in rng.sample(pieces, rng.randint(0, 3))]
        remnants.append(prefix + ''.join(parts))
    return remnants


def validate_args(args):
    if args.posidFile and args.repeatTable:
        raise ValueError('posidFile and repeatTable are mutually exclusive.')
    if not args.repeatTable and (args.repeatClasses or args.repeatGrouping != 'class'
                                 or args.includeRepeatAlternates):
        raise ValueError('Repeat selection options require --repeatTable.')
    if args.mode not in ('integrase', 'AAV'):
        raise ValueError('mode must be integrase or AAV.')
    for name in ('nSites', 'nFrags', 'nReadsPerFrag', 'distBetweenFragEnds',
                 'nReplicates', 'R1_length', 'R2_length', 'I1_length'):
        if getattr(args, name) < 1:
            raise ValueError(f'{name} must be at least 1.')
    if args.nFrags > 100:
        raise ValueError('nFrags must be no more than 100.')
    if args.nFrags * args.distBetweenFragEnds > 500:
        raise ValueError('nFrags * distBetweenFragEnds must be no more than 500.')
    if args.minSiteDistance < 0:
        raise ValueError('minSiteDistance must be non-negative.')
    for name, maximum in (('percentGenomicError', 1), ('positionChatterSD', 5)):
        value = getattr(args, name)
        if not math.isfinite(value) or not 0 <= value <= maximum:
            raise ValueError(f'{name} must be finite and between 0 and {maximum}.')
    rows = args.nReplicates * (1 if args.singleSample else 9)
    if (rows - 1).bit_length() > 2 * args.I1_length:
        raise ValueError('I1_length is too short to give every replicate a unique barcode.')


def validate_read_lengths(args, remnants):
    if args.R1_length <= len(LINKER):
        raise ValueError(f'R1_length must exceed the {len(LINKER)}-nt linker to include genomic sequence.')
    if args.R2_length <= max(map(len, remnants)):
        raise ValueError('R2_length must exceed every leader length to include genomic sequence; '
                         f'this seed/mode requires at least {max(map(len, remnants)) + 1} nt.')
    shortest_fragment = (BASE_FRAGMENT_END + args.distBetweenFragEnds - START_POSITION
                         - 2 * chatter_radius(args.positionChatterSD))
    longest_payload = max(args.R1_length - len(LINKER),
                          args.R2_length - min(map(len, remnants)))
    if longest_payload > shortest_fragment:
        raise ValueError('Requested read lengths exceed the shortest fragment allowed by '
                         'distBetweenFragEnds and positionChatterSD.')


@dataclass(frozen=True)
class RepeatAnnotation:
    group: str
    repeat_class: str
    name: str
    start: int
    end: int
    strand: str
    line: int

    @property
    def length(self):
        return self.end - self.start + 1


@dataclass(frozen=True)
class Site:
    chromosome: str
    position: int  # Reported, one-based integration position, including mode correction.
    strand: str
    sequence: str  # Oriented 1000-nt block with CONTEXT_PADDING bases at each end.
    repeat: object = None

    @property
    def posid(self):
        return f'{self.chromosome}{self.strand}{self.position}'


def block_start(position, strand, mode):
    """Invert the original posid offsets (+12/-11 integrase; +10/-9 AAV)."""
    if strand == '+':
        return position - (12 if mode == 'integrase' else 10)
    return position + (11 if mode == 'integrase' else 9) - BLOCK_LENGTH


def reference_sequence(record):
    return record.seq if hasattr(record, 'seq') else record


def reference_bounds(position, strand):
    starts = [block_start(position, strand, mode) for mode in ('integrase', 'AAV')]
    return min(starts) - CONTEXT_PADDING, max(starts) + BLOCK_LENGTH + CONTEXT_PADDING


def site_at_position(genome, chromosome, position, strand, mode, repeat=None):
    posid = f'{chromosome}{strand}{position}'
    if chromosome not in genome.keys():
        raise ValueError(f'{posid}: chromosome is absent from the reference.')
    reference = reference_sequence(genome[chromosome])
    union_start, union_end = reference_bounds(position, strand)
    if position < 1 or union_start < 0 or union_end > len(reference):
        raise ValueError(f'{posid}: insufficient reference context for the requested site.')
    window = str(reference[union_start:union_end]).upper()
    if len(window) != union_end - union_start or set(window) - set('ACGT'):
        raise ValueError(f'{posid}: reference context contains non-ACGT bases.')
    relative_start = block_start(position, strand, mode) - CONTEXT_PADDING - union_start
    sequence = window[relative_start:relative_start + BLOCK_LENGTH + 2 * CONTEXT_PADDING]
    if strand == '-':
        sequence = reverse_complement(sequence)
    return Site(chromosome, position, strand, sequence, repeat)


def open_text(path):
    path = Path(path).expanduser()
    if path.suffix.lower() in ('.gz', '.bgz'):
        return gzip.open(path, 'rt', encoding='utf-8-sig', newline='')
    return path.open('r', encoding='utf-8-sig', newline='')


def read_posids(args):
    """Preserve input order and exact sites; spacing constraints do not move them."""
    specifications, seen = [], set()
    pattern = re.compile(r'([^\s+\-]+)([+-])([1-9]\d*)')
    with open_text(args.posidFile) as handle:
        for line_number, line in enumerate(handle, 1):
            value = line.strip()
            if not value or value.startswith('#'):
                continue
            match = pattern.fullmatch(value)
            if not match:
                raise ValueError(f'posidFile line {line_number}: invalid position ID {value!r}.')
            if value in seen:
                raise ValueError(f'posidFile line {line_number}: duplicate position ID {value}.')
            seen.add(value)
            chromosome, strand, position = match.groups()
            specifications.append((chromosome, int(position), strand))
    if not specifications:
        raise ValueError('posidFile contains no position IDs.')
    if args.nSitesExplicit and args.nSites != len(specifications):
        raise ValueError(f'--nSites={args.nSites} disagrees with the {len(specifications)} sites in posidFile.')
    args.nSites = len(specifications)
    return specifications


def supplied_sites(genome, args, specifications):
    for chromosome, position, strand in specifications:
        yield site_at_position(genome, chromosome, position, strand, args.mode)


def spacing_index(positions, position, minimum_gap):
    insertion = bisect_left(positions, position)
    if ((insertion and position - positions[insertion - 1] < minimum_gap)
            or (insertion < len(positions) and positions[insertion] - position < minimum_gap)):
        return None
    return insertion


def select_sites(genome, args, max_attempts=MAX_SITE_ATTEMPTS):
    """Sample reported coordinates uniformly over eligible genomic positions.

    Both mode-specific blocks must be usable so switching modes cannot change
    placement. Reference slices stay lazy when backed by Biopython's 2bit reader.
    """
    names, cumulative = [], []
    total = 0
    for name in [f'chr{i}' for i in range(1, 100)] + ['chrX']:
        if name not in genome.keys():
            continue
        available = len(reference_sequence(genome[name])) - 2 * CHROMOSOME_MARGIN
        if available > 0:
            names.append(name)
            total += available
            cumulative.append(total)
    if not names:
        raise ValueError('No eligible chromosomes: require chr1-chr99 or chrX with '
                         f'more than {2 * CHROMOSOME_MARGIN} bases.')
    if args.nSites > total:
        raise ValueError('nSites exceeds the number of eligible genomic positions.')
    chosen = {name: [] for name in names}
    minimum_gap = max(1, args.minSiteDistance)
    for site_index in range(args.nSites):
        rng = make_rng(args.seed, 'site', site_index)
        for _ in range(max_attempts):
            draw = rng.randrange(total)
            chromosome_index = bisect_right(cumulative, draw)
            chromosome = names[chromosome_index]
            previous_total = cumulative[chromosome_index - 1] if chromosome_index else 0
            position = CHROMOSOME_MARGIN + 1 + draw - previous_total
            positions = chosen[chromosome]
            insertion = spacing_index(positions, position, minimum_gap)
            if insertion is None:
                continue
            strand = rng.choice(('+', '-'))
            try:
                site = site_at_position(genome, chromosome, position, strand, args.mode)
            except ValueError:
                continue
            positions.insert(insertion, position)
            yield site
            break
        else:
            raise ValueError(f'Could not place site {site_index + 1} after {max_attempts} attempts '
                             f'with minSiteDistance={args.minSiteDistance}. The reference may lack '
                             'usable sequence or the requested sites may be too densely packed.')


class RepeatPool:
    """Compact columns instead of millions of per-row Python objects."""
    def __init__(self):
        self.starts, self.ends, self.lines = array('Q'), array('Q'), array('Q')
        self.chromosomes, self.names, self.classes = array('I'), array('I'), array('I')
        self.strands = array('B')

    def append(self, chromosome, start, end, name, repeat_class, strand, line):
        self.chromosomes.append(chromosome)
        self.starts.append(start)
        self.ends.append(end)
        self.names.append(name)
        self.classes.append(repeat_class)
        self.strands.append(strand == 'C' or strand == '-')
        self.lines.append(line)

    def __len__(self):
        return len(self.starts)


def center_third_ranges(start, end, chromosome_length):
    trim = (end - start + 1) // 3
    lower, upper = start + trim, end - trim
    ranges = []
    for strand in ('+', '-'):
        context_start, context_end = reference_bounds(0, strand)
        left = max(lower, 1, -context_start)
        right = min(upper, chromosome_length - context_end)
        if left <= right:
            ranges.append((strand, left, right))
    return ranges


def intern(value, lookup, values):
    if value not in lookup:
        lookup[value] = len(values)
        values.append(value)
    return lookup[value]


def repeat_selectors(args):
    selectors = []
    for item in args.repeatClasses or []:
        for value in item.split(','):
            value = value.strip()
            if not value:
                raise ValueError('repeatClasses contains an empty group name.')
            key = 'sine/alu' if value.casefold() == 'alu' else value.casefold()
            if key in selectors:
                raise ValueError(f'Duplicate repeat group: {value}.')
            selectors.append(key)
    return selectors


def load_repeat_pools(genome, args):
    selectors = repeat_selectors(args)
    matched = set()
    pools = {}
    chromosomes, names, classes = [], [], []
    chromosome_ids, name_ids, class_ids = {}, {}, {}
    lengths = {name: len(reference_sequence(genome[name])) for name in genome.keys()}
    with open_text(args.repeatTable) as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader, None)
        required = ('query_seq', 'query_start', 'query_end', 'repeat_name', 'repeat_class', 'strand')
        if not header or len(header) != len(set(header)) or any(name not in header for name in required):
            raise ValueError('repeatTable requires a unique TSV header with query_seq, query_start, '
                             'query_end, repeat_name, repeat_class, and strand.')
        column = {name: header.index(name) for name in required}
        for optional in ('alt', 'strand'):
            if optional in header:
                column[optional] = header.index(optional)
        for line_number, row in enumerate(reader, 2):
            if not row or not any(row):
                continue
            if len(row) < max(column[name] for name in required) + 1:
                raise ValueError(f'repeatTable line {line_number}: missing required fields.')
            value = lambda name, default='': row[column[name]].strip() if name in column and column[name] < len(row) else default
            repeat_class = value('repeat_class')
            if not repeat_class:
                raise ValueError(f'repeatTable line {line_number}: empty repeat_class.')
            if selectors:
                keys = [key for key in selectors if repeat_class.casefold() == key
                        or repeat_class.casefold().startswith(key + '/')]
                if len(keys) > 1:
                    raise ValueError(f'Overlapping repeatClasses selections for {repeat_class}: {keys}.')
                if not keys:
                    continue
                key = keys[0]
                group = repeat_class[:len(key)]
                matched.add(key)
            else:
                group = repeat_class.split('/')[0] if args.repeatGrouping == 'class' else repeat_class
            if value('alt') == '*' and not args.includeRepeatAlternates:
                continue
            chromosome = value('query_seq')
            if chromosome not in lengths:
                continue
            try:
                start, end = int(value('query_start')), int(value('query_end'))
            except ValueError:
                raise ValueError(f'repeatTable line {line_number}: invalid query coordinates.') from None
            if not 1 <= start <= end <= lengths[chromosome]:
                raise ValueError(f'repeatTable line {line_number}: interval outside {chromosome} or start > end.')
            if not center_third_ranges(start, end, lengths[chromosome]):
                continue
            if value('strand') not in ('+', '-', 'C'):
                raise ValueError(f'repeatTable line {line_number}: strand must be +, -, or C.')
            if not value('repeat_name'):
                raise ValueError(f'repeatTable line {line_number}: empty repeat_name.')
            if group not in pools:
                pools[group] = RepeatPool()
            pools[group].append(intern(chromosome, chromosome_ids, chromosomes), start, end,
                intern(value('repeat_name'), name_ids, names), intern(repeat_class, class_ids, classes),
                value('strand'), line_number)
    if selectors and (len(matched) != len(selectors) or len(pools) != len(selectors)):
        available = {group.casefold() for group in pools}
        missing = [key for key in selectors if key not in matched or key not in available]
        raise ValueError(f'Requested repeat groups have no eligible intervals: {", ".join(missing)}.')
    if not pools:
        raise ValueError('No repeat intervals have a center-third position with sufficient reference context.')
    return pools, chromosomes, names, classes, lengths


def shuffled_indices(rng, count):
    remaining, swaps = count, {}
    while remaining:
        index = rng.randrange(remaining)
        value = swaps.get(index, index)
        remaining -= 1
        swaps[index] = swaps.get(remaining, remaining)
        swaps.pop(remaining, None)
        yield value


def repeat_group_order(pools, seed):
    groups = sorted(pools)
    make_rng(seed, 'repeat-group-order').shuffle(groups)
    return groups


def select_repeat_sites(genome, args, max_attempts=MAX_SITE_ATTEMPTS, prepared_pools=None):
    pools, chromosomes, names, classes, lengths = (
        load_repeat_pools(genome, args) if prepared_pools is None else prepared_pools)
    groups = repeat_group_order(pools, args.seed)
    if args.nSites < len(groups):
        raise ValueError(f'nSites must be at least the number of repeat groups ({len(groups)}), '
                         'or select fewer groups with --repeatClasses.')
    quotient, remainder = divmod(args.nSites, len(groups))
    for index, group in enumerate(groups):
        requested = quotient + (index < remainder)
        if requested > len(pools[group]):
            raise ValueError(f'Repeat group {group} needs {requested} sites but has only '
                             f'{len(pools[group])} eligible annotation rows. Reduce nSites or select broader groups.')
    candidates = {group: shuffled_indices(make_rng(args.seed, 'repeat-elements', group), len(pools[group]))
                  for group in groups}
    chosen = {chromosome: [] for chromosome in chromosomes}
    for site_index in range(args.nSites):
        group = groups[site_index % len(groups)]
        pool = pools[group]
        rng = make_rng(args.seed, 'repeat-site', group, site_index // len(groups))
        found = False
        for _ in range(max_attempts):
            index = next(candidates[group], None)
            if index is None:
                break
            chromosome = chromosomes[pool.chromosomes[index]]
            ranges = center_third_ranges(pool.starts[index], pool.ends[index], lengths[chromosome])
            total = sum(right - left + 1 for _, left, right in ranges)
            # Retry a few positions in long annotations before considering the next row.
            for _ in range(min(32, total)):
                draw = rng.randrange(total)
                for strand, left, right in ranges:
                    width = right - left + 1
                    if draw < width:
                        position = left + draw
                        break
                    draw -= width
                positions = chosen[chromosome]
                insertion = spacing_index(positions, position, max(1, args.minSiteDistance))
                if insertion is None:
                    continue
                annotation = RepeatAnnotation(group, classes[pool.classes[index]], names[pool.names[index]],
                    pool.starts[index], pool.ends[index], '-' if pool.strands[index] else '+', pool.lines[index])
                try:
                    site = site_at_position(genome, chromosome, position, strand, args.mode, annotation)
                except ValueError:
                    continue
                positions.insert(insertion, position)
                yield site
                found = True
                break
            if found:
                break
        if not found:
            raise ValueError(f'Could not meet the site quota for repeat group {group}. '
                             'Insufficient distinct intervals satisfy spacing and reference-sequence requirements; '
                             'reduce nSites/minSiteDistance or select broader groups.')


def interval_within_repeat(site, first, stop, mode):
    """Test a genomic slice of the oriented reference against the selected row."""
    origin = block_start(site.position, site.strand, mode) - CONTEXT_PADDING
    if site.strand == '+':
        left, right = origin + first + 1, origin + stop
    else:
        left = origin + len(site.sequence) - stop + 1
        right = origin + len(site.sequence) - first
    return site.repeat.start <= left <= right <= site.repeat.end


def repeat_truth(site):
    annotation = site.repeat
    left, right = site.position - annotation.start, annotation.end - site.position
    return dict(repeatGroup=annotation.group, repeatClass=annotation.repeat_class,
        repeatName=annotation.name, repeatStart=annotation.start, repeatEnd=annotation.end,
        repeatLength=annotation.length, repeatStrand=annotation.strand, repeatTableLine=annotation.line,
        repeatDistanceToStart=left, repeatDistanceToEnd=right, repeatBoundaryDistance=min(left, right),
        nR1GenomicWithinRepeat=0, nR2GenomicWithinRepeat=0, nReadPairsFullyWithinRepeat=0)


def simulate_error(sequence, fraction, rng):
    """Substitute round(fraction * len(sequence)) distinct, randomly chosen bases."""
    count = round(fraction * len(sequence))
    if not count:
        return sequence
    bases = list(sequence)
    for position in rng.sample(range(len(bases)), count):
        bases[position] = rng.choice([base for base in DNA if base != bases[position]])
    return ''.join(bases)


def build_sample_rows(args):
    groups = [('test', 'subjectA', 'sample1')] if args.singleSample else [
        ('test', subject, sample) for subject in ('subjectA', 'subjectB', 'subjectC')
        for sample in ('sample1', 'sample2', 'sample3')]
    barcodes = unique_dna(make_rng(args.seed, 'barcodes'), args.I1_length,
                          len(groups) * args.nReplicates)
    rows = []
    for trial, subject, sample in groups:
        for replicate in range(1, args.nReplicates + 1):
            row = dict(replicate=replicate, subject=subject, sample=sample, trial=trial,
                       index1Seq=next(barcodes), adriftReadLinkerSeq=LINKER, refGenome=args.refGenomeID)
            if args.mode == 'integrase':
                row.update(leaderSeqHMM=args.integraseHMM, mode='U5', vectorFastaFile='synDataTest.fasta')
            else:
                row.update(anchorReadStartSeq=args.anchorReadStartSeq, mode='AAV',
                           vectorFastaFile='validationVector.fasta')
            rows.append(row)
    return groups, rows


class ReplicateAllocator:
    """Balance whole fragments within sites, repeat groups, and a sample.

    strata contains (repeat_group, site_count) pairs for one sample, or a
    single (None, site_count) pair when repeat targeting is not used.
    Every fragment number visits every replicate once per nReplicates sites
    within its group. This also works when fragment and replicate counts
    share a divisor; simply advancing by nFrags would lock fragment lengths
    to a subset of replicates in that case.

    Full cycles of R/gcd(F,R) sites balance fragment totals. Shifting the
    phase between those cycles balances fragment numbers as well. Starting
    phases place each group's final remainder consecutively across the
    sample, keeping final sample totals within one fragment. Each group's
    totals and each fragment number's counts also differ by at most one.
    Assignment is deterministic and consumes no sequence-generation RNG.
    """
    def __init__(self, n_fragments, n_replicates, strata):
        self.n_fragments = n_fragments
        self.n_replicates = n_replicates
        self.divisor = math.gcd(n_fragments, n_replicates)
        self.cycle_sites = n_replicates // self.divisor
        self.offsets, self.counts, self.quotas = {}, {}, {}
        fragment_cursor = site_cursor = 0
        for group, count in strata:
            extra = (count * n_fragments) % n_replicates
            if extra:
                tail_phase = (count // self.cycle_sites) % self.divisor
                offset = (fragment_cursor - tail_phase) % n_replicates
            else:
                # Totals are already equal: stagger fragment lengths using
                # the number of sites allocated to preceding groups.
                offset = self.phase(site_cursor) % n_replicates
            self.offsets[group] = offset
            self.counts[group] = 0
            self.quotas[group] = count
            fragment_cursor = (fragment_cursor + extra) % n_replicates
            site_cursor += count

    def phase(self, site_number):
        return ((site_number % self.cycle_sites) * self.n_fragments
                + (site_number // self.cycle_sites) % self.divisor)

    def assign(self, group):
        if group not in self.counts or self.counts[group] >= self.quotas[group]:
            raise ValueError('Site order disagrees with the replicate allocation plan.')
        first = self.offsets[group] + self.phase(self.counts[group])
        self.counts[group] += 1
        return tuple((first + fragment) % self.n_replicates
                     for fragment in range(self.n_fragments))


def sample_strata(repeat_groups, first_site, site_count):
    """Count the interleaved repeat groups in one sample's contiguous sites."""
    if not repeat_groups:
        return [(None, site_count)]
    quotient, remainder = divmod(site_count, len(repeat_groups))
    return [(repeat_groups[(first_site + index) % len(repeat_groups)],
             quotient + (index < remainder))
            for index in range(len(repeat_groups))
            if quotient or index < remainder]


class IlluminaReadIDs:
    """Stream unique, randomized synthetic cluster IDs in constant memory.

    A seeded permutation of 40-bit counters addresses four lanes, 64 tiles
    per lane, and 65,536 x 65,536 integer coordinates per tile. Addition,
    multiplication by odd numbers modulo 2**40, and right-xor shifts are
    all invertible: distinct read-pair counters cannot collide. This avoids
    a growing set of used IDs or random-draw collision retries.

    Geometry is synthetic; the FASTQ field layout follows the supplied
    Illumina examples. IDs are reproducible for the same seed and record
    order. Header randomness never consumes sequence-generation streams.
    """
    CAPACITY = 1 << 40

    def __init__(self, seed, read_pairs):
        if not 0 <= read_pairs <= self.CAPACITY:
            raise ValueError('--illuminaReadIDs supports at most '
                             f'{self.CAPACITY:,} read pairs per dataset.')
        rng = make_rng(seed, 'illumina-read-ids')
        self.prefix = f'SYNTH01:1:S{rng.getrandbits(32):08X}'
        self.mask = self.CAPACITY - 1
        self.offset = rng.getrandbits(40)
        self.multiplier1 = rng.getrandbits(40) | (1 << 39) | 1
        self.multiplier2 = rng.getrandbits(40) | (1 << 39) | 1
        self.counter = 0

    def next_cluster(self):
        if self.counter >= self.CAPACITY:
            raise ValueError('Illumina cluster-coordinate space is exhausted.')
        value = (self.counter + self.offset) & self.mask
        self.counter += 1
        value ^= value >> 20
        value = (value * self.multiplier1) & self.mask
        value ^= value >> 17
        value = (value * self.multiplier2) & self.mask
        value ^= value >> 20
        y = 1000 + (value & 0xffff)
        x = 1000 + ((value >> 16) & 0xffff)
        tile_index = (value >> 32) & 63
        lane = (value >> 38) + 1
        surface, remainder = divmod(tile_index, 32)
        swath, tile = divmod(remainder, 16)
        tile_name = f'{surface + 1}{swath + 1}1{tile + 1:02d}'
        return f'{self.prefix}:{lane}:{tile_name}:{x}:{y}'


def write_fastq(handle, read_id, sequence):
    handle.write(f'@{read_id}\n{sequence}\n+\n{"?" * len(sequence)}\n')


def gzip_writer(stack, path):
    raw = stack.enter_context(path.open('wb'))
    compressed = stack.enter_context(gzip.GzipFile(filename='', mode='wb', fileobj=raw,
                                                   compresslevel=6, mtime=0))
    return stack.enter_context(io.TextIOWrapper(compressed, encoding='ascii', newline='\n'))


def generate(args, genome):
    """Stream a complete dataset into staging files before replacing prior outputs."""
    specifications = read_posids(args) if args.posidFile else None
    validate_args(args)
    remnants = build_remnant_library(args.seed, args.mode)
    validate_read_lengths(args, remnants)
    illumina_ids = (IlluminaReadIDs(args.seed, args.nSites * args.nFrags * args.nReadsPerFrag)
                   if args.illuminaReadIDs else None)
    groups, rows = build_sample_rows(args)
    output = Path(args.outputDir).expanduser()
    output.mkdir(parents=True, exist_ok=True)
    repeat_groups = None
    if specifications is not None:
        sites = supplied_sites(genome, args, specifications)
    elif args.repeatTable:
        # Load the compact index once. Its group order also defines each
        # sample's allocation quotas; selected reference windows stay lazy.
        prepared_pools = load_repeat_pools(genome, args)
        repeat_groups = repeat_group_order(prepared_pools[0], args.seed)
        sites = select_repeat_sites(genome, args, prepared_pools=prepared_pools)
    else:
        sites = select_sites(genome, args)
    sites_per_group, remainder = divmod(args.nSites, len(groups))
    with tempfile.TemporaryDirectory(prefix='.buildSynData-', dir=output) as temporary:
        staging = Path(temporary)
        with ExitStack() as stack:
            r1 = gzip_writer(stack, staging / 'R1.fastq.gz')
            r2 = gzip_writer(stack, staging / 'R2.fastq.gz')
            i1 = gzip_writer(stack, staging / 'I1.fastq.gz')
            sample_handle = stack.enter_context((staging / 'sampleData.tsv').open('w', newline='', encoding='utf-8'))
            sample_writer = csv.DictWriter(sample_handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
            sample_writer.writeheader()
            sample_writer.writerows(rows)
            truth_handle = stack.enter_context((staging / 'truth.tsv').open('w', newline='', encoding='utf-8'))
            truth_fields = ('trial', 'subject', 'sample', 'posid', 'nReads', 'nFrags', 'nUMIs', 'leaderSeq')
            if args.repeatTable:
                truth_fields += REPEAT_TRUTH_FIELDS
            truth_writer = csv.DictWriter(truth_handle, fieldnames=truth_fields,
                delimiter='\t', lineterminator='\n')
            truth_writer.writeheader()
            site_index = 0
            for group_index, (trial, subject, sample) in enumerate(groups):
                group_rows = rows[group_index * args.nReplicates:(group_index + 1) * args.nReplicates]
                sample_site_count = sites_per_group + (group_index < remainder)
                allocation = ReplicateAllocator(args.nFrags, args.nReplicates,
                    sample_strata(repeat_groups, site_index, sample_site_count))
                for _ in range(sample_site_count):
                    site = next(sites)
                    replicate_indices = allocation.assign(site.repeat.group if site.repeat else None)
                    leader = make_rng(args.seed, 'leader', site_index).choice(remnants)
                    umis = unique_dna(make_rng(args.seed, 'umis', site_index), 12, args.nFrags)
                    truth = dict(trial=trial, subject=subject, sample=sample,
                        posid=site.posid, nReads=args.nFrags * args.nReadsPerFrag,
                        nFrags=args.nFrags, nUMIs=args.nFrags, leaderSeq=leader)
                    if site.repeat is not None:
                        truth.update(repeat_truth(site))
                    read_number = 1
                    for fragment_index, umi in enumerate(umis, 1):
                        replicate_index = replicate_indices[fragment_index - 1]
                        barcode = group_rows[replicate_index]['index1Seq']
                        linker = LINKER.replace('N' * 12, umi)
                        start_rng = make_rng(args.seed, 'chatter-start', site_index, fragment_index)
                        end_rng = make_rng(args.seed, 'chatter-end', site_index, fragment_index)
                        r1_errors = make_rng(args.seed, 'errors-R1', site_index, fragment_index)
                        r2_errors = make_rng(args.seed, 'errors-R2', site_index, fragment_index)
                        for _ in range(args.nReadsPerFrag):
                            start = CONTEXT_PADDING + START_POSITION - 1 + gaussian_offset(start_rng, args.positionChatterSD)
                            end = (CONTEXT_PADDING + BASE_FRAGMENT_END + fragment_index * args.distBetweenFragEnds - 1
                                   + gaussian_offset(end_rng, args.positionChatterSD))
                            genomic_r1 = reverse_complement(site.sequence[end - (args.R1_length - len(linker)):end])
                            genomic_r2 = site.sequence[start:start + args.R2_length - len(leader)]
                            # Length validation ensures both reads fit inside this fragment.
                            if not (0 <= start < end <= len(site.sequence)
                                    and max(len(genomic_r1), len(genomic_r2)) <= end - start):
                                raise ValueError('Fragment endpoint or read length falls outside reference context.')
                            if illumina_ids is None:
                                id_r1 = id_r2 = f'{site.posid}_read{read_number}_frag{fragment_index}'
                            else:
                                cluster = illumina_ids.next_cluster()
                                id_r1 = f'{cluster} 1:N:0:{barcode}'
                                id_r2 = f'{cluster} 2:N:0:{barcode}'
                            read_number += 1
                            sequence_r1 = linker + simulate_error(genomic_r1, args.percentGenomicError, r1_errors)
                            # Mutate the assembled R2 so vector remnants can receive errors too.
                            sequence_r2 = simulate_error(leader + genomic_r2, args.percentGenomicError, r2_errors)
                            if len(sequence_r1) != args.R1_length or len(sequence_r2) != args.R2_length:
                                raise ValueError('Generated read length does not match the requested length.')
                            if site.repeat is not None:
                                truth['nR1GenomicWithinRepeat'] += interval_within_repeat(
                                    site, end - len(genomic_r1), end, args.mode)
                                truth['nR2GenomicWithinRepeat'] += interval_within_repeat(
                                    site, start, start + len(genomic_r2), args.mode)
                                truth['nReadPairsFullyWithinRepeat'] += interval_within_repeat(
                                    site, start, end, args.mode)
                            write_fastq(r1, id_r1, sequence_r1)
                            write_fastq(r2, id_r2, sequence_r2)
                            write_fastq(i1, id_r1, barcode)
                    truth_writer.writerow(truth)
                    site_index += 1
        for name in OUTPUT_NAMES:
            os.replace(staging / name, output / name)
        # Discard only obsolete, generator-owned uncompressed streams, after success.
        for name in ('R1.fastq', 'R2.fastq', 'I1.fastq'):
            (output / name).unlink(missing_ok=True)
    return [output / name for name in OUTPUT_NAMES]


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        validate_args(args)
        from Bio.SeqIO import TwoBitIO
        # Keep the file open for lazy reference access throughout generation.
        with Path(args.refGenomePath).expanduser().open('rb') as handle:
            genome = TwoBitIO.TwoBitIterator(handle)
            if args.posidFile:
                print(f'Generating {args.mode} sites from {args.posidFile}.')
            else:
                print(f'Generating {args.nSites} {args.mode} sites.')
            generate(args, genome)
    except ImportError as error:
        parser.exit(1, f'Error: Biopython is required to read the 2bit reference: {error}\n')
    except (OSError, ValueError) as error:
        parser.exit(1, f'Error: {error}\n')
    print(f'Wrote synthetic reads, sampleData.tsv, and truth.tsv to {args.outputDir}.')
    return 0


if __name__ == '__main__':
    sys.exit(main())

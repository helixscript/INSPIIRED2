#!/usr/bin/env python3
import os
import argparse
import random
import pandas
import numpy
import re
import sys
import yaml
import subprocess
import itertools
from Bio.SeqIO import TwoBitIO

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--outputDir', type=str, default = './output', help='Output directory path. Default (./output).', metavar='')
parser.add_argument('-m', '--mode', type=str, default = 'integrase', help='Mode (integrase or AAV). Default (integrase).', metavar='')
parser.add_argument('-n', '--nSites', type=int, default = 1000, help='Number of sites to generate. Default (1000).', metavar='')
parser.add_argument('-f', '--nFrags', type=int, default = 5, help='Number of fragments for each site (limit 100). Default (5).', metavar='')
parser.add_argument('-p', '--nReadsPerFrag', type=int, default = 25, help='Number of reads per fragment. Default (25).', metavar='')
parser.add_argument('-d', '--distBetweenFragEnds', type=int, default = 25, help='Distance (NT) between fragment break points. Default (25).', metavar='')
parser.add_argument('-s', '--seed', type=int, default = 1, help='Random seed. Default (1).', metavar='')
parser.add_argument('-r', '--refGenomePath', type=str, default = '../../../data/referenceGenomes/hg38.2bit', help='Path to 2bit reference genome. Default (../../../data/referenceGenomes/hg38.2bit).', metavar='')
parser.add_argument('-g', '--refGenomeID', type=str, default = 'hg38', help='Reference genome ID for sample data table. Default (hg38).', metavar='')
parser.add_argument('-i', '--integraseHMM', type=str, default = 'validation.hmm', help='Name of AAVengeR HMM to use (integrase mode only).  Default (validation.hmm).', metavar='')
parser.add_argument('-t', '--anchorReadStartSeq', type=str, default = 'TCTGCGCGCT', help='Anchor read start sequence filter (AAV mode only). Default (TCTGCGCGCT).', metavar='')
parser.add_argument('-x', '--R1_length', type=int, default = 150, help='Total length of R1 reads. Default (150).', metavar='')
parser.add_argument('-y', '--R2_length', type=int, default = 150, help='Total length of R2 reads. Default (150).', metavar='')
parser.add_argument('-z', '--I1_length', type=int, default = 12, help='Total length of I1 reads. Default (12).', metavar='')
parser.add_argument('-e', '--percentGenomicError', type=float, default = 0, help='Percent gDNA error (0.0 - 1.0) to simulate in R1 and R2 reads. Default (0).', metavar='')
parser.add_argument('-c', '--positionChatterSD', type=float, default = 0.50, help='StdDev of Gaussian centered on expected fragment ends used to simulate position chatter. (Default 0.50).', metavar='')
parser.add_argument('-k', '--singleSample', action='store_true', default=False, help='Consolidate all synthetic sites into a single sample.')
parser.add_argument('-b', '--minSiteDistance', type=int, default = 500, help='Minimum distance (NT) between simulated integration sites. Default (500).', metavar='')
parser.add_argument('-v', '--nReplicates', type=int, default = 4, help='Number of technical replicates per sample. Default (4).', metavar='')

args = parser.parse_args()
args.refGenomePath = os.path.expanduser(args.refGenomePath)

# Test inputs 
if not os.path.exists(args.refGenomePath):
  print('Error - refGenomePath does not')
  sys.exit(1)

if args.mode not in ['integrase', 'AAV']:
  print('Error - mode must be set to "integrase" or "AAV".')
  sys.exit(1)

if args.nFrags > 100:
  print('Error - nFrags must be set to a value no more than 100.')
  sys.exit(1)

if args.nReplicates < 1:
  print('Error - nReplicates must be at least 1.')
  sys.exit(1)

if args.percentGenomicError < 0 or args.percentGenomicError > 1:
  print('Error - percentGenomicError must be set to a value between 0 and 1.')
  sys.exit(1)

if args.positionChatterSD < 0 or args.positionChatterSD > 5:
  print('Error - positionChatterSD must be set to a value between 0 and 5')
  sys.exit(1)

if(args.nFrags * args.distBetweenFragEnds > 500):
  print('Error - the number of requested fragments x the distance between fragment ends can not exceed 500.')
  sys.exit(1)

# Set random seed.
random.seed(args.seed)

# Create output directory.
if not os.path.isdir(args.outputDir):
  os.mkdir(args.outputDir)

if not os.path.exists(args.outputDir):
  print('Error - could not create output directory.')
  sys.exit(1)

# Empty out output dir if it contains a previous result.
files = [os.path.join(args.outputDir, 'R1.fastq.gz'),
         os.path.join(args.outputDir, 'R2.fastq.gz'),
         os.path.join(args.outputDir, 'I1.fastq.gz'),
         os.path.join(args.outputDir, 'testVector.fasta'),
         os.path.join(args.outputDir, 'sampleData.tsv'),
         os.path.join(args.outputDir, 'config.yml'),
         os.path.join(args.outputDir, 'truth.tsv')]

for f in files:
  if os.path.exists(f):
    os.remove(f)

# Linkers to be used for faux sites.
remnant0 = 'TCTGCGCGCTCGCTCGCTCA' # To be used for integrase mode.
remnant  = 'TCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTG' 
linker   = 'GAACGAGCACTAGTAAGCCCNNNNNNNNNNNNCTCCGCTTAAGGGACT' 

# Helper functions.
def flatten(xss):
    return [x for xs in xss for x in xs]

def randomBarCode(length = 12):
  NTs = ['A', 'T', 'C', 'G']
  return(''.join(random.choices(NTs, k = length)))

def compSeq(seq):
  compBases = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
  return(''.join([compBases[x] for x in seq]))

def revCompSeq(seq):
  return(compSeq(seq)[::-1])

def buildGaussDistSample(pos, n, popSize = 1000, sd = 1):
  nums = []
  for x in range(popSize):
    nums.append(round(random.gauss(pos, sd)))
  return(random.sample(nums, n))

def mutateAtPos(string, index):
    NTs = ['A', 'T', 'C', 'G']
    exclude = list(string[index])
    NTs =[x for x in NTs if x not in exclude]
    return string[:index] + random.sample(NTs, 1)[0] + string[index + 1:]

def simulateError(seq):
  nPos = round(args.percentGenomicError * len(seq))
  orgSeq = seq
  locs = random.sample(list(range(1, len(seq))), nPos)
  for x in locs:
    seq = mutateAtPos(seq, x)
  return(seq)

# Build a collection of scrambled remnants to draw from. 
# Pieces will be between 12 and 36 NTs with a 50% of being flipped to RC.
start = []
end = []

print('Building remnant library.')

for i in range(500):
  while True:
    a = random.randint(12, len(remnant)-12)
    b = random.randint(a+12, len(remnant))
    
    if (b - a + 1) >= 12 & (b - a + 1) <= 36:
      start.append(a)
      end.append(b)
      break

pieces = []
for i in range(len(start)):
  pieces.append(remnant[start[i-1]:end[i-1]])

remnants = []
for i in range(500):
   n = random.randint(0, 3)
   if n == 0:
     remnants.append(remnant[0:random.randint(12, 24)])
   else:
      # Start remnant with a variable length fragment starting from the beginning.
      o = [remnant[0:random.randint(12, 24)]]
      # Randomly draw 1 - 3 remnant pieces.
      k = random.sample(range(1, len(pieces)), n)
      p = [pieces[x] for x in k]
      
      # Randomly flip pieces to RC.
      p2 = []
      for r in p:
         if random.randint(0, 1) == 0:
            p2.append(r)
         else:
            p2.append(revCompSeq(r))
      o.append(p2)
      o = flatten(o)
      remnants.append(''.join(o))

if args.mode == 'integrase':
   remnants = [remnant0]

# Read in genomic data pointers.
print('Reading reference genome data.')
_handle = open(args.refGenomePath, 'rb')          # keep file open for lazy-loading
g = TwoBitIO.TwoBitIterator(_handle)

# Create a list of allowed chromosomes from which to sample.
o = list(range(1, 100))
o = [str(x) for x in o]
o.extend(['X'])  # Exclude Y, often full of repeats and poorly characterized.
allowedChromosomes = ['chr' + x for x in o]
chromosomes = []
for c in g.keys():
    if c in allowedChromosomes:
        chromosomes.append(c)

# Read select chromosomes into memory, slow.
d = {}
for c in g.keys():
    if c in chromosomes:
        print('  Reading chromosome', c, 'into memory.')
        d[c] = str(g[c].seq)   # g[c] is a SeqRecord; use .seq

del g
_handle.close()

# Select chromosomes for sites.
chromosomes = random.choices(chromosomes, k = args.nSites)

# Compile data to build sites.
class site:
    def __init__(self, chr, pos, seq, strand):
        self.chr = chr
        self.pos = pos
        self.seq = seq
        self.strand = strand

sites = []
print("Compiling site building data.")

# Keep track of chosen positions per chromosome to enforce minimum distance
chosen_positions = {c: [] for c in set(chromosomes)}

for c in chromosomes:
  attempts = 0
  while True:
    # Prevent infinite loops if a chromosome is packed or very small
    attempts += 1
    if attempts > 1000:
       print(f"Warning: Struggling to find appropriately spaced sites on {c}. Proceeding with closest fit.")
       
    # Select a position within the chromosome excluding the ends.
    pos = random.randint(1+10000, len(d[c])-10000)
    
    # 1. Check distance against previously chosen sites on this chromosome
    too_close = False
    for existing_pos in chosen_positions[c]:
      if abs(pos - existing_pos) < args.minSiteDistance:
        too_close = True
        break
        
    if too_close and attempts <= 1000:
       continue # Try a new random position

    # 2. Check for Ns in the 1000 NT block
    seq = d[c][pos:pos+1000]
    if 'N' not in seq.upper():
      chosen_positions[c].append(pos)
      break # Found a valid site!

  strand = random.sample(['+', '-'], 1)[0]

  if strand == '-':
    pos = pos + 1000
    seq = revCompSeq(seq)

  sites.append(site(c, pos, seq, strand))

# Loop through site objects and build read sequences.
print("Building site reads.")
for s in sites:
  s.readIDs = []
  s.UMIs = []
  s.fragmentSeqs = []
  s.R1 = []
  s.R2 = []

  nReads = args.nFrags * args.nReadsPerFrag

  # Start within excised gDNA chunks to allow chatter around start position.
  # Create start positions around position 10.
  startPositions = buildGaussDistSample(10, nReads, sd = args.positionChatterSD)

  # Create end positions simulating sonic breaks.
  baseFragWidth = 500
  readNum = 1

  # Create alt start postion for read ids to correct for zero-based 2bit system and starting 10 NT within seq blocks.
  # Integrase vs AAV mode invoke different gDNA duplication corrections in buildSites.
  if s.strand == '+':
    if(args.mode == 'integrase'):
       s.pos2 = s.pos + 12
    else:
       s.pos2 = s.pos + 10
  else:
    if(args.mode == 'integrase'):
      s.pos2 = s.pos - 11
    else:
      s.pos2 = s.pos - 9

  endPositions = []

  for x in range(1, args.nFrags + 1):
    # Define a frag specific UMI sequence.
    UMI = ''.join(random.choices(['A', 'T', 'C', 'G'], k = 12))

    for i in range(1, args.nReadsPerFrag + 1):
      s.readIDs.append(s.chr + s.strand + str(s.pos2) + '_read' + str(readNum) + '_frag' + str(x))
      s.UMIs.append(UMI)
      readNum += 1

    offSet = buildGaussDistSample(baseFragWidth + (x * args.distBetweenFragEnds), args.nReadsPerFrag, sd = args.positionChatterSD)
    endPositions.append(offSet)

  endPositions = flatten(endPositions)

  # Extract substrings for fragment sequences.
  for x in range(len(s.readIDs)):
     s.fragmentSeqs.append(s.seq[startPositions[x]-1:endPositions[x]-1])

  # Randomly select a remnant. The same remnant will always be selected when in 'integrase' mode.
  s.remnant = remnants[random.randint(0, len(remnants)-1)]

  # Build R2 and R1 read sequences.
  for x in range(len(s.fragmentSeqs)):
    gDNA = s.fragmentSeqs[x][0:args.R2_length - len(s.remnant)]

    if(args.percentGenomicError > 0):
      saved_state = random.getstate()
      gDNA = simulateError(gDNA)
      random.setstate(saved_state)

    s.R2.append(s.remnant + gDNA)

  for x in range(len(s.fragmentSeqs)):
    gDNA = revCompSeq(s.fragmentSeqs[x]) [0:args.R1_length - len(linker)]

    if(args.percentGenomicError > 0):
      saved_state = random.getstate()
      gDNA = simulateError(gDNA)
      random.setstate(saved_state)

    s.R1.append(linker + gDNA)

# Build sample data file.
print("Building sample table.")
reps_seq = list(range(1, args.nReplicates + 1))

if args.singleSample:
  sampleData = pandas.DataFrame({
    'replicate': reps_seq,
    'subject': ['subjectA'] * args.nReplicates,
    'sample': ['sample1'] * args.nReplicates,
    'trial': 'test',
    'index1Seq': [randomBarCode() for i in range(args.nReplicates)],
    'adriftReadLinkerSeq': linker,
    'refGenome': args.refGenomeID
  })
else:
  # 3 subjects * nReplicates = total rows per subject
  sub_len = 3 * args.nReplicates
  sampleData = pandas.DataFrame({
    'replicate': reps_seq * 9,
    'subject': flatten([['subjectA'] * sub_len, ['subjectB'] * sub_len, ['subjectC'] * sub_len]),
    'sample': flatten([['sample1'] * args.nReplicates, ['sample2'] * args.nReplicates, ['sample3'] * args.nReplicates]*3),
    'trial': 'test',
    'index1Seq': [randomBarCode() for i in range(9 * args.nReplicates)],
    'adriftReadLinkerSeq': linker,
    'refGenome': args.refGenomeID
  })

# Add mode specific columns.
if(args.mode == 'integrase'):
  sampleData['leaderSeqHMM'] = args.integraseHMM
  sampleData['mode'] = 'U5'
  sampleData['vectorFastaFile'] = 'synDataTest.fasta'
else:
  sampleData['anchorReadStartSeq'] = args.anchorReadStartSeq
  sampleData['mode'] = 'AAV'
  sampleData['vectorFastaFile'] = 'validationVector.fasta'

# Split sample table and sites into sample groups.
sampleGroups = sampleData.groupby(['trial', 'subject', 'sample'])
siteGroups = numpy.array_split(sites, sampleGroups.ngroups)

z = 0
truths = []

# Create a continuous round-robin dealer for replicates
rep_cycler = itertools.cycle(range(1, args.nReplicates + 1))

print("Building I1 reads and writing all reads to output.")
for (key, group) in sampleGroups:
    # Create list of replicate bar codes associated with this sample.
    barCodes = list(group['index1Seq'])

    # Loop through sites in c
    for x in siteGroups[z]:

      # Keep fragments within replicates since they are standardized within replicates.
      fragIDs = [re.search(r'frag\d+', m).group(0) for m in x.readIDs]
      
      # Extract unique fragments in order, and deal them to the next available replicate
      unique_frags = list(dict.fromkeys(fragIDs))
      site_fragToRep = {f: next(rep_cycler) for f in unique_frags}
      
      reps = [site_fragToRep[n] for n in fragIDs]
      barCodesToUse = [barCodes[x-1] for x in reps]

      # Add entry to truths.
      truths.append({'trial': key[0], 'subject': key[1], 'sample': key[2], 'posid': x.readIDs[0].split('_')[0], 
                     'nReads': args.nFrags * args.nReadsPerFrag, 'nFrags': args.nFrags, 'nUMIs': args.nFrags, 
                     'leaderSeq': x.remnant})

      for i in range(len(x.readIDs)):
        # Replace NNN in linker with fragment UMIs, build R1 FASTQ, write.
        x.R1[i] = x.R1[i].replace('NNNNNNNNNNNN', x.UMIs[i])
        o = ['@' + x.readIDs[i], x.R1[i], '+', ''.join('?'*len(x.R1[i]))]

        with open(os.path.join(args.outputDir, 'R1.fastq'), 'a') as f:
           for line in o:
             f.write("%s\n" % line)

        # Build and write R2.
        o = ['@' + x.readIDs[i], x.R2[i], '+', ''.join('?'*len(x.R2[i]))]

        with open(os.path.join(args.outputDir, 'R2.fastq'), 'a') as f:
         for line in o:
           f.write("%s\n" % line)

        # Build and write I1.
        # Select a barcode from the ones associated with this sample.
        b = barCodesToUse[i]
        o = ['@' + x.readIDs[i], b, '+', ''.join('?'*len(b))]

        with open(os.path.join(args.outputDir, 'I1.fastq'), 'a') as f:
         for line in o:
           f.write("%s\n" % line)
    z += 1

print("Writing tabular outputs.")
# Write out sample table.
sampleData.to_csv(os.path.join(args.outputDir, 'sampleData.tsv'), sep='\t', index=False)

# Write out truth table.
t = pandas.DataFrame.from_records([truth for truth in truths])
t.to_csv(os.path.join(args.outputDir, 'truth.tsv'), sep='\t', index=False)

# Compress fastq files.
os.system('gzip ' + args.outputDir + '/*.fastq')

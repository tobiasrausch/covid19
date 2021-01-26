#! /usr/bin/env python

import csv
import argparse
import sys
import collections
import cyvcf2
import os
import json

# Parse command line
parser = argparse.ArgumentParser(description='Aggregate QC statistics')
parser.add_argument('-p', '--prefix', metavar='prefix', required=True, dest='prefix', help='file prefix (required)')
args = parser.parse_args()

qc = dict()

# Trim galore
for rn in ["1", "2"]:
    filep = args.prefix + "." + rn + ".fq.gz_trimming_report.txt"
    if (os.path.exists(filep)) and (os.path.isfile(filep)):
        with open(filep) as f:
            for line in f:
                if line.startswith("Reads with adapters"):
                    qc['AdaptersRead' + rn] = line[(line.find('(')+1):line.find(')')]

# Host reads (GRCh38 mapping)
filep = args.prefix + ".host.reads"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with open(filep) as f:
        allreads = None
        grch38reads = None
        for line in f:
            if line.strip().endswith(".all.reads"):
                allreads = int(line.strip().split(' ')[0])
            if line.strip().endswith(".remove.reads"):
                grch38reads = int(line.strip().split(' ')[0])
        qc['#ReadPairs'] = allreads
        qc['FractionGRCh38'] = float(grch38reads) / float(allreads)

# Host reads (kraken2 human DB)
filep = args.prefix + ".kraken2.report.txt"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with open(filep) as f:
        for line in f:
            if line.strip().endswith("Homo sapiens"):
                fields = ' '.join(line.split()).split(' ')
                qc['Kraken2AddHumanReads'] = fields[2]

# SARS-CoV-2 reads
filep = args.prefix + ".srt.bam.flagstat"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with open(filep) as f:
        for line in f:
            if line.strip().endswith("read1"):
                fields = ' '.join(line.split()).split(' ')
                sarsreads = int(fields[0])
                qc['FractionSars'] = float(sarsreads) / float(allreads)                

# Coverage
filep = args.prefix + ".depth"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    f_reader = csv.reader(open(filep), delimiter="\t")
    npos = 0
    zerocov = 0
    below10 = 0
    totcov = 0
    for fields in f_reader:
        cov = int(fields[2])
        if cov == 0:
            zerocov += 1
        if cov < 10:
            below10 += 1
        totcov += cov
        npos += 1
    qc['#CoverageEqual0'] = zerocov
    qc['#CoverageBelow10'] = below10
    qc['ReferenceLength'] = npos
    qc['MeanCoverage'] = float(totcov) / float(npos)

# Mutation file
filep = args.prefix + ".mutation.csv"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    f_reader = csv.reader(open(filep), delimiter=",")
    nline = 0
    for fields in f_reader:
        if nline == 1:
            ct = collections.Counter(fields)
            qc['#MissingKeyMutations'] = ct['X']
            qc['MutationString'] = ','.join(fields)
        nline += 1

# Variants
filep = args.prefix + ".bcf"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    vcf = cyvcf2.VCF(filep)
    ncount = 0
    indels = False
    mtct = collections.Counter()
    for record in vcf:
        mt = record.REF + ">" + record.ALT[0]
        if len(mt) != 3:
            indels = True
        mtct[mt] += 1
        ncount += 1
    qc['MutationTypes'] = json.dumps(mtct).replace(' ','')
    qc['#CalledVariants'] = ncount
    qc['InDelsPresent'] = indels

# Consensus composition
filep = args.prefix + ".cons.comp"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with open(filep) as f:
        ncount = 0
        ambcount = 0
        for line in f:
            fields = ' '.join(line.split()).split(' ')
            if fields[1] == 'N':
                ncount += int(fields[0])
            if (fields[1] != 'A') and (fields[1] != 'C') and (fields[1] != 'G') and (fields[1] != 'T') and (fields[1] != 'N'):
                ambcount += int(fields[0])
        qc['#ConsensusNs'] = ncount
        qc['#ConsensusAmbiguous'] = ambcount

# Freebayes & iVar diff
filep = args.prefix + ".cons.diff"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with open(filep) as f:
        diffct = 0
        for line in f:
            if line.startswith('<'):
                fields = ' '.join(line.split()).split(' ')
                # Any diff for non-N and non-ambiguous
                if fields[1] not in ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]:
                    diffct += 1
        qc['iVarFreeBayesDiff'] = diffct

# Determine success/fail for this sample
qc['outcome'] = "fail"
if qc["#CalledVariants"] < 30:
    if qc['#ConsensusNs'] < 1000:
        if qc['#MissingKeyMutations'] == 0:
            if qc['FractionGRCh38'] < 0.5:
                if qc['iVarFreeBayesDiff'] == 0:
                    qc['outcome'] = "pass"
        
# Output QC dictionary
for key in sorted(qc.keys()):
    print(key, qc[key])

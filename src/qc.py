#! /usr/bin/env python

import csv
import argparse
import collections
import os
import json
import cyvcf2

# Parse command line
parser = argparse.ArgumentParser(description='Aggregate QC statistics')
parser.add_argument('-p', '--prefix', required=True, dest='prefix', help='file prefix (required)')
args = parser.parse_args()

qc = dict()
qc['Sample'] = args.prefix

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
        qc['MappingFractionGRCh38'] = str(round(float(grch38reads) / float(allreads), 4))

# Host reads (kraken2 human DB)
kraken2addreads = 0
filep = args.prefix + ".kraken2.report.txt"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with open(filep) as f:
        for line in f:
            if line.strip().endswith("Homo sapiens"):
                fields = ' '.join(line.split()).split(' ')
                kraken2addreads = fields[2]
qc['Kraken2AddHumanReads'] = kraken2addreads

# SARS-CoV-2 reads
filep = args.prefix + ".srt.bam.flagstat"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with open(filep) as f:
        for line in f:
            if line.strip().endswith("read1"):
                fields = ' '.join(line.split()).split(' ')
                sarsreads = int(fields[0])
                qc['MappingFractionSars'] = str(round(float(sarsreads) / float(allreads), 4))

# Parse alignment statistics
filep = args.prefix + ".alfred.tsv"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with open(filep) as f:
        columns = None
        for line in f:
            if line.startswith("ME"):
                if columns is None:
                    columns = line.strip().split('\t')
                else:
                    records = line.strip().split('\t')
                    qc['FractionUnmapped'] = records[columns.index('UnmappedFraction')]
                    qc['FractionSupplementaryAlignments'] = records[columns.index('SupplementaryAlignmentFraction')]
                    qc['FractionMismatch'] = records[columns.index('MismatchRate')]
                    qc['FractionInsertion'] = records[columns.index('InsertionRate')]
                    qc['FractionDeletion'] = records[columns.index('DeletionRate')]
                    qc['FractionSoftClip'] = records[columns.index('SoftClipRate')]
                    qc['FractionSequencingErrors'] = records[columns.index('ErrorRate')]
                    qc['MedianCoverage'] = records[columns.index('MedianCoverage')]
                    qc['SDCoverage'] = records[columns.index('SDCoverage')]
                    qc['MedianInsertSize'] = records[columns.index('MedianInsertSize')]

# Coverage
filep = args.prefix + ".depth"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    f_reader = csv.reader(open(filep), delimiter="\t")
    npos = 0
    zerocov = 0
    for fields in f_reader:
        cov = int(fields[2])
        if cov == 0:
            zerocov += 1
        npos += 1
    qc['#CoverageEqual0'] = zerocov
    qc['ReferenceLength'] = npos

# Mutation file
filep = args.prefix + ".mutation.csv"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    f_reader = csv.reader(open(filep), delimiter=",")
    nline = 0
    for fields in f_reader:
        if nline == 1:
            ct = collections.Counter(fields)
            qc['#CovDropoutsKeyMutations'] = ct['X']
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
    qc['MutationTypes'] = json.dumps(mtct).replace(' ', '')
    qc['#CalledVariants'] = ncount
    qc['InDelMutationsPresent'] = indels

# Parse variants tsv file
filep = args.prefix + ".variants.tsv"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    f_reader = csv.DictReader(open(filep), delimiter="\t")
    varct = collections.Counter()
    for fields in f_reader:
        varct[fields['Consequence']] += 1
    qc['Consequence'] = json.dumps(varct).replace(' ', '')

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

# Primer trimming
filep = args.prefix + ".iVar.trim"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with open(filep) as f:
        for line in f:
            if line.startswith('Trimmed primers from'):
                fields = line.strip().split(' ')
                qc['PrimerTrimmed'] = fields[3]
            elif line.strip().endswith('and were not written to file.'):
                fields = line.strip().split(' ')
                qc['PrimerTrimmedTooShort'] = fields[0]
            elif line.strip().endswith('insert size smaller than their read length'):
                fields = line.strip().split(' ')
                qc['PrimerTrimmedISizeIssue'] = fields[0]

# Freebayes & iVar diff
filep = args.prefix + ".cons.diff"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with open(filep) as f:
        diffct = 0
        ivarct = 0
        freect = 0
        for line in f:
            if line.startswith('<'):
                ivarct += 1
            elif line.startswith('>'):
                freect += 1
            if line.startswith('<'):
                fields = ' '.join(line.split()).split(' ')
                # Any diff for non-N and non-ambiguous
                if fields[1] not in ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]:
                    diffct += 1
        absdiff = abs(ivarct - freect)
        if diffct > absdiff:
            qc['iVarFreeBayesDiff'] = diffct
        else:
            qc['iVarFreeBayesDiff'] = absdiff

# Parse pangolin lineage
filep = args.prefix + ".lineage.csv"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    f_reader = csv.DictReader(open(filep), delimiter=",")
    for fields in f_reader:
        qc['Lineage'] = fields['lineage']
        qc['PangolinStatus'] = fields['status']
        qc['LineageProb'] = fields['probability']

# Nextclade
qc['NextcladeStatus'] = None
filep = args.prefix + ".json"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with open(filep) as json_file:
        data = json.load(json_file)
        for p in data:
            if 'clade' in p.keys():
                qc['Clade'] = p['clade']
            else:
                qc['Clade'] = None
            if 'qc' in p.keys():
                qc['NextcladeStatus'] = p['qc']['overallStatus']

# Vadr
qc['VadrStatus'] = None
filep = args.prefix + ".out.vadr.pass.tbl"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    qc['VadrStatus'] = 'fail'
    with open(filep) as f:
        for line in f:
            if line.startswith('>'):
                qc['VadrStatus'] = 'pass'
                break

# Determine success/borderline/fail for this sample
qc['outcome'] = "fail"
if qc["#CalledVariants"] < 50:
    if qc['#ConsensusNs'] < 5000:
        ncstatus = 'good'
        if qc['NextcladeStatus'] is not None:
            ncstatus = qc['NextcladeStatus']
        if (qc['PangolinStatus'] == 'passed_qc') and (ncstatus == 'good'):
            qc['outcome'] = "pass"
        elif (qc['PangolinStatus'] == 'passed_qc') or (ncstatus == 'good'):
            qc['outcome'] = "borderline"
        else:
            qc['outcome'] = "fail"

# Output QC dictionary
for key in sorted(qc.keys()):
    print(key, qc[key])

#! /usr/bin/env python

from readfq import readfq
import csv
import argparse
import collections
import os
import json
import gzip

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
allreads = None
filep = args.prefix + ".host.reads"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with open(filep) as f:
        grch38reads = None
        for line in f:
            if line.strip().endswith(".all.reads"):
                allreads = int(line.strip().split(' ')[0])
            if line.strip().endswith(".remove.reads"):
                grch38reads = int(line.strip().split(' ')[0])
        qc['PercHuman'] = str(round(float(grch38reads) / float(allreads) * 100, 2)) + "%"

# SARS-CoV-2 reads
filep = args.prefix + ".srt.bam.flagstat"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with open(filep) as f:
        for line in f:
            if line.strip().endswith("read1"):
                fields = ' '.join(line.split()).split(' ')
                sarsreads = int(fields[0])
                qc['PercSars'] = str(round(float(sarsreads) / float(allreads) * 100, 2)) + "%"

# Parse alignment statistics
filep = args.prefix + ".alfred.tsv.gz"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with gzip.open(filep, "rb") as f:
        columns = None
        for line in f.read().decode("ascii").splitlines():
            if line.startswith("ME"):
                if columns is None:
                    columns = line.strip().split('\t')
                else:
                    records = line.strip().split('\t')
                    qc['PercUnmapped'] = str(round(float(records[columns.index('UnmappedFraction')]) * 100, 2)) + "%"
                    qc['PercSeqError'] = str(round(float(records[columns.index('ErrorRate')]) * 100, 2)) + "%"
                    qc['MedianCoverage'] = records[columns.index('MedianCoverage')]
                    qc['SDCoverage'] = records[columns.index('SDCoverage')]
                    qc['MedianInsertSize'] = records[columns.index('MedianInsertSize')]

# Percent identity to SARS-CoV-2 reference
filep = args.prefix + ".align.fa.gz"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with gzip.open(filep, "rb") as f:
        reflen = 0
        seqlen = 0
        nummatch = 0
        for line in f.read().decode("ascii").splitlines():
            columns = line.strip().split('\t')
            if columns[2] != '-':
                reflen += 1
            if columns[3] != '-':
                seqlen += 1
                if columns[2] == columns[3]:
                    nummatch += 1
        qc['PercIdentity'] = str(round(float(nummatch) / float(reflen) * 100, 2)) + "%"

# Mutation file
filep = args.prefix + ".mutation.csv"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    f_reader = csv.reader(open(filep), delimiter=",")
    header = None
    for fields in f_reader:
        if header is None:
            header = fields
            continue
        strout = []
        for (k, v) in zip(header[1:], fields[1:]):
            strout.append(k + ":" + v)
        qc['S_Typing'] = ','.join(strout)

# Parse variants tsv file
filep = args.prefix + ".variants.tsv"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    f_reader = csv.DictReader(open(filep), delimiter="\t")
    varsetS = []
    for fields in f_reader:
        if ('IMPACT' in fields.keys()) and ('SYMBOL' in fields.keys()):
            if (fields['IMPACT'] == 'MODERATE') or (fields['IMPACT'] == 'HIGH'):
                aa = fields['Amino_acids'].split('/')
                if len(aa) == 2:
                    if fields['SYMBOL'] == 'S':
                        aachange = aa[0] + str(fields['Protein_position']) + aa[1]
                        if aachange not in varsetS:
                            varsetS.append(aachange)
    qc['S_Variants'] = ','.join(varsetS)

# Consensus composition
filep = args.prefix + ".cons.fa"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with open(filep) as f:
        ncount = 0
        ambcount = 0
        dnacount = 0
        for seqName, seqNuc, seqQuals in readfq(f):
            for bp in seqNuc:
                if bp == 'N':
                    ncount += 1
                elif (bp == 'A') or (bp == 'C') or (bp == 'G') or (bp == 'T'):
                    dnacount += 1
                else:
                    ambcount += 1
        nuclen = ncount + dnacount + ambcount
        qc['#ConsensusNs'] = ncount
        qc['#ConsensusAmbiguous'] = ambcount
        qc['PercN'] = str(round(float(ncount) / float(nuclen) * 100, 2)) + "%"
        qc['PercACGT'] = str(round(float(dnacount) / float(nuclen) * 100, 2)) + "%"

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
        lcount = 0
        for line in f:
            lcount += 1
        if lcount > 0:
            qc['iVarFreeBayesDiff'] = True
        else:
            qc['iVarFreeBayesDiff'] = False

# Parse pangolin lineage
qc['PangolinStatus'] = None
filep = args.prefix + ".lineage.csv"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    f_reader = csv.DictReader(open(filep), delimiter=",")
    for fields in f_reader:
        qc['Lineage'] = fields['lineage']
        qc['PangolinStatus'] = fields['status']

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
qc['outcome'] = "pass"
if (qc['NextcladeStatus'] is not None) and (qc['PangolinStatus'] is not None):
    if (qc['PangolinStatus'] == 'passed_qc') and (qc['NextcladeStatus'] == 'good'):
        qc['outcome'] = "pass"
    elif (qc['PangolinStatus'] == 'passed_qc') and (qc['#ConsensusNs'] <= 5000):
        qc['outcome'] = "borderline"
    else:
        qc['outcome'] = "fail"
else:
    if qc['#ConsensusNs'] > 5000:
        qc['outcome'] = "fail"

# Check percent identity, median coverage and percent ACGT
qc['rki'] = 'fail'
if float(qc['PercIdentity'][:-1]) >= 90:
    if float(qc['PercN'][:-1]) <= 5:
        if float(qc['MedianCoverage']) >= 100:
            if float(qc['PercACGT'][:-1]) >= 90:
                qc['rki'] = "pass"

# Output QC dictionary
for key in ["Sample", "outcome", "rki", "Lineage", "Clade", "MedianCoverage", "#ConsensusAmbiguous", "#ConsensusNs", "PercACGT", "PercHuman", "PercIdentity", "PercN", "PercSars", "PercSeqError", "PercUnmapped", "S_Variants", "AdaptersRead1", "AdaptersRead2", "MedianInsertSize", "PrimerTrimmed", "PrimerTrimmedISizeIssue", "PrimerTrimmedTooShort", "SDCoverage", "VadrStatus", "NextcladeStatus", "PangolinStatus", "iVarFreeBayesDiff", "S_Typing"]:
    if key in qc.keys():
        print(key, qc[key])

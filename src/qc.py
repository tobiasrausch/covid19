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
        if allreads != 0:
            qc['PercHuman'] = str(round(float(grch38reads) / float(allreads) * 100, 2)) + "%"
        else:
            qc['PercHuman'] = "0%"

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
        qc['Typing_S'] = ','.join(strout)

# Parse variants tsv file
filep = args.prefix + ".variants.tsv"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    f_reader = csv.DictReader(open(filep), delimiter="\t")
    varsetS = []
    vartot = []
    for fields in f_reader:
        if ('IMPACT' in fields.keys()) and ('SYMBOL' in fields.keys()):
            if (fields['IMPACT'] == 'MODERATE') or (fields['IMPACT'] == 'HIGH'):
                aa = fields['Amino_acids'].split('/')
                if len(aa) == 2:
                    aachange = aa[0] + str(fields['Protein_position']) + aa[1]
                    if fields['SYMBOL'] == 'S':
                        if aachange not in varsetS:
                            varsetS.append(aachange)
                    aachange = fields['SYMBOL'] + ":" + aachange
                    if aachange not in vartot:
                        vartot.append(aachange)
    qc['Mutations_S'] = ','.join(varsetS)
    qc['Mutations'] = ','.join(vartot)

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
        if nuclen != 0:
            qc['PercN'] = str(round(float(ncount) / float(nuclen) * 100, 2)) + "%"
            qc['PercACGT'] = str(round(float(dnacount) / float(nuclen) * 100, 2)) + "%"
        else:
            qc['PercN'] = "100%"
            qc['PercACGT'] = "0%"

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
filep = args.prefix + ".cons.diff.gz"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with gzip.open(filep, "rb") as f:
        ivardiff = 0
        for line in f.read().decode("ascii").splitlines():
            columns = line.strip().split('\t')
            if columns[2] == columns[3]:
                continue
            if (columns[2] == 'N') and ((columns[3] == 'A') or (columns[3] == 'C') or (columns[3] == 'G') or (columns[3] == 'T') or (columns[3] == '-')):
                continue
            if (columns[2] == 'R') and ((columns[3] == 'A') or (columns[3] == 'G')):
                continue
            if (columns[2] == 'Y') and ((columns[3] == 'C') or (columns[3] == 'T')):
                continue
            if (columns[2] == 'S') and ((columns[3] == 'C') or (columns[3] == 'G')):
                continue
            if (columns[2] == 'W') and ((columns[3] == 'A') or (columns[3] == 'T')):
                continue
            if (columns[2] == 'K') and ((columns[3] == 'G') or (columns[3] == 'T')):
                continue
            if (columns[2] == 'M') and ((columns[3] == 'A') or (columns[3] == 'C')):
                continue
            if (columns[2] == 'V') and ((columns[3] == 'A') or (columns[3] == 'C') or (columns[3] == 'G')):
                continue
            if (columns[2] == 'H') and ((columns[3] == 'A') or (columns[3] == 'C') or (columns[3] == 'T')):
                continue
            if (columns[2] == 'D') and ((columns[3] == 'A') or (columns[3] == 'G') or (columns[3] == 'T')):
                continue
            if (columns[2] == 'B') and ((columns[3] == 'C') or (columns[3] == 'G') or (columns[3] == 'T')):
                continue
            ivardiff += 1
        qc['iVarFreeBayesDiff'] = ivardiff

# Parse pangolin lineage
qc['PangolinStatus'] = None
qc['scorpio_call'] = None
filep = args.prefix + ".lineage.csv"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    f_reader = csv.DictReader(open(filep), delimiter=",")
    for fields in f_reader:
        qc['Lineage'] = fields['lineage']
        if 'status' in fields.keys():
            qc['PangolinStatus'] = fields['status']
        else:
            qc['PangolinStatus'] = fields['qc_status']
        qc['scorpio_call'] = fields['scorpio_call']

# Nextclade
qc['NextcladeStatus'] = None
qc['Clade'] = None
filep = args.prefix + ".json"
if (os.path.exists(filep)) and (os.path.isfile(filep)):
    with open(filep) as json_file:
        data = json.load(json_file)
        data = data['results']
        for p in data:
            if 'clade' in p.keys():
                qc['Clade'] = p['clade']
            if 'qc' in p.keys():
                qc['NextcladeStatus'] = p['qc']['overallStatus']
            if qc['Lineage'] == "Unassigned":
                if 'customNodeAttributes' in p.keys():
                    if 'Nextclade_pango' in p['customNodeAttributes']:
                        qc['Lineage'] = p['customNodeAttributes']['Nextclade_pango']

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
qc['QC'] = "pass"
if (qc['NextcladeStatus'] is not None) and (qc['PangolinStatus'] is not None):
    if (qc['PangolinStatus'] == 'pass') and (qc['NextcladeStatus'] == 'good'):
        qc['QC'] = "pass"
    elif (qc['PangolinStatus'] == 'pass') and (qc['#ConsensusNs'] <= 5000):
        qc['QC'] = "borderline"
    else:
        qc['QC'] = "fail"
else:
    if qc['#ConsensusNs'] > 5000:
        qc['QC'] = "fail"

# Check percent identity, median coverage and percent ACGT
qc['RKI'] = 'fail'
if float(qc['PercIdentity'][:-1]) >= 90:
    if float(qc['PercN'][:-1]) <= 5:
        if float(qc['MedianCoverage']) >= 100:
            if float(qc['PercACGT'][:-1]) >= 90:
                if qc['#ConsensusAmbiguous'] <= 25:
                    qc['RKI'] = "pass"

# Assign simplified type
qc['Type'] = None
if qc ['RKI'] == "pass":
    if qc['Lineage'] == "B.1.1.7":
        qc['Type'] = "B.1.1.7"
    elif qc['Lineage'].startswith("Q."):
        qc['Type'] = "B.1.1.7"
    elif qc['Lineage'] == "B.1.351":
        qc['Type'] = "B.1.351"
    elif qc['Lineage'] == "B.1.1.529":
        qc['Type'] = "B.1.1.529"
    elif qc['Lineage'].startswith("BA."):
        qc['Type'] = "B.1.1.529"
    elif qc['Lineage'].startswith("BE."):
        qc['Type'] = "B.1.1.529"
    elif qc['Lineage'].startswith("BF."):
        qc['Type'] = "B.1.1.529"
    elif qc['Lineage'].startswith("BG."):
        qc['Type'] = "B.1.1.529"
    elif qc['Lineage'].startswith("BH."):
        qc['Type'] = "B.1.1.529"
    elif qc['Lineage'].startswith("BK."):
        qc['Type'] = "B.1.1.529"
    elif qc['Lineage'].startswith("BL."):
        qc['Type'] = "B.1.1.529"
    elif qc['Lineage'].startswith("BN."):
        qc['Type'] = "B.1.1.529"
    elif qc['Lineage'].startswith("BQ."):
        qc['Type'] = "B.1.1.529"
    elif qc['Lineage'].startswith("CE."):
        qc['Type'] = "B.1.1.529"
    elif qc['Lineage'].startswith("CC."):
        qc['Type'] = "B.1.1.529"
    elif qc['Lineage'] == "A.27":
        qc['Type'] = "A.27"
    elif qc['Lineage'] == "P.1":
        qc['Type'] = "P.1"
    elif qc['Lineage'].startswith("P.1."):
        qc['Type'] = "P.1"
    elif qc['Lineage'] == "P.2":
        qc['Type'] = "P.2"
    elif qc['Lineage'] == "P.3":
        qc['Type'] = "P.3"
    elif qc['Lineage'] == "B.1.617":
        qc['Type'] = "B.1.617"
    elif qc['Lineage'] == "B.1.617.1":
        qc['Type'] = "B.1.617.1"
    elif qc['Lineage'] == "B.1.617.2":
        qc['Type'] = "B.1.617.2"
    elif qc['Lineage'].startswith("AY."):
        qc['Type'] = "B.1.617.2"
    elif qc['Lineage'] == "B.1.617.3":
        qc['Type'] = "B.1.617.3"
    elif (qc['Lineage'].startswith("X")) and (qc['Clade'] == "recombinant"):
        qc['Type'] = "recombinant"
    elif (qc['Lineage'] == "None") or (qc['Lineage'] == "Unassigned"):
        if qc['scorpio_call'].startswith("Omicron"):
            qc['Type'] = "B.1.1.529"
        else:
            qc['Type'] = "None"
    else:
        qc['Type'] = "WT"


# Output QC dictionary
for key in ["Sample", "QC", "RKI", "Lineage", "Clade", "Type", "Mutations_S", "Mutations", "MedianCoverage", "#ConsensusAmbiguous", "#ConsensusNs", "PercACGT", "PercHuman", "PercIdentity", "PercN", "PercSars", "PercSeqError", "PercUnmapped", "AdaptersRead1", "AdaptersRead2", "MedianInsertSize", "PrimerTrimmed", "PrimerTrimmedISizeIssue", "PrimerTrimmedTooShort", "SDCoverage", "VadrStatus", "NextcladeStatus", "PangolinStatus", "iVarFreeBayesDiff", "Typing_S"]:
    if key in qc.keys():
        print(key, qc[key])

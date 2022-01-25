#! /usr/bin/env python

import csv
import argparse
import sys
import collections
import json
import cyvcf2

# Parse command line
parser = argparse.ArgumentParser(description='Create variant table')
parser.add_argument('-v', '--variants', metavar='sample.bcf', required=True, dest='variants', help='VCF file (required)')
parser.add_argument('-d', '--depth', metavar='sample.depth', required=True, dest='depth', help='depth file (required)')
parser.add_argument('-c', '--coverage', metavar='10', default=10, type=int, required=False, dest='thres', help='coverage threshold (optional)')
parser.add_argument('-l', '--lowest', metavar='LOW', default='LOW', required=False, dest='lowest', help='lowest VEP impact to report (optional)')
parser.add_argument('-s', '--sample', metavar='s1', required=True, dest='sample', help='sample name (required)')
parser.add_argument('-o', '--output', metavar='out.csv', required=True, dest='outfile', help='output file (required)')
args = parser.parse_args()

# Typing JSON
typing = '[{"key": [3266, 3268], "value":"t1001i"},{"key": [21614, 21616], "value":"l18f"},{"key": [22226, 22228], "value":"a222v"},{"key": [22466, 22468], "value":"t302t"},{"key": [22877, 22879], "value":"n439k"},{"key": [22916, 22918], "value":"l452r"},{"key": [22919, 22921], "value":"y453f"},{"key": [23063, 23065], "value":"n501y"},{"key": [23402, 23404], "value":"d614g"},{"key": [23519, 23521], "value":"a653v"},{"key": [23525, 23527], "value":"h655y"},{"key": [23603, 23605], "value":"p681h"},{"key": [23948, 23950], "value":"d796y"},{"key": [25217, 25219], "value":"g1219v"},{"key": [27972, 27974], "value":"q27*"},{"key": [21764, 21764], "value":"del_21765_6"}]'
typedict = dict()
for i in json.loads(typing):
    typedict[(i['key'][0], i['key'][1])] = i['value']

# Get lowest impact to report
impact = dict({'MODIFIER': 0, 'LOW': 1, 'MODERATE': 2, 'HIGH': 3})
if args.lowest not in impact:
    print("-l parameter needs to be [MODIFIER|LOW|MODERATE|HIGH]", file=sys.stderr)
    sys.exit(1)
lowestimpact = impact[args.lowest]

# Estimate threshold
invalidpos = set()
if args.depth:
    f_reader = csv.reader(open(args.depth), delimiter="\t")
    for fields in f_reader:
        cov = int(fields[2])
        pos = int(fields[1])
        if cov < args.thres:
            invalidpos.add(pos)

# VEP columns
vepcols = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID", "CLIN_SIG", "SOMATIC", "PHENO"]
descols = ["Consequence", "IMPACT", "SYMBOL", "Feature", "cDNA_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation"]
addr = dict()
for idx, val in enumerate(vepcols):
    addr[val] = idx

# VCF/BCF parsing
varstore = collections.defaultdict(collections.defaultdict)
vcf = cyvcf2.VCF(args.variants)
print("sample", "reference", "position", "REF", "ALT", '\t'.join(descols), sep='\t')
for record in vcf:
    csq = record.INFO.get('CSQ')
    if csq is None:
        print(args.sample, record.CHROM, record.POS, record.REF, ','.join(record.ALT), "unknown_variant", sep='\t')
        continue
    transcripts = csq.split(',')
    for tr in transcripts:
        fields = tr.split('|')
        if len(fields) != len(vepcols):
            print("Error: VEP annotation is corrupted!", file=sys.stderr)
            sys.exit(1)
        if impact[fields[addr['IMPACT']]] >= impact['MODERATE']:
            size = len(record.ALT[0]) - len(record.REF)   # deletions negative, insertions positive
            varstore[record.POS][fields[addr['Feature']]] = (size, fields[addr['Protein_position']], fields[addr["Amino_acids"]])
        if impact[fields[addr['IMPACT']]] >= lowestimpact:
            print(args.sample, record.CHROM, record.POS, record.REF, ','.join(record.ALT), '\t'.join([fields[addr[cname]] for cname in descols]), sep='\t')


# Create observed amino acid string (X: unobserved)
with open(args.outfile, 'w') as fout:
    obs = []
    for (start, end) in sorted(typedict.keys()):
        mt = typedict[(start, end)]
        aaobs = None
        for pos in range(start, end+1):
            if pos in invalidpos:
                aaobs = 'X'
        if aaobs is None:
            if end - start == 2:
                # Set default
                aaobs = mt[0].upper()
                aapos = mt[1:(len(mt)-1)]
                # Any called variants?
                for pos in range(start, end+1):
                    if pos in varstore.keys():
                        for key in varstore[pos]:
                            (size, dbaapos, aachange) = varstore[pos][key]
                            if (size == 0) and (dbaapos == aapos):  # Only SNVs
                                # Get target amino acid
                                aaobs = aachange[len(aachange)-1:]
            else:
                aaobs = 'ref'
                if start in varstore.keys():
                    for key in varstore[start]:
                        (size, dbaapos, aachange) = varstore[pos][key]
                        if size == -6:
                            aaobs = 'del'
        obs.append(aaobs)
    print('sample', ','.join([typedict[(start, end)] for (start, end) in sorted(typedict.keys())]), sep=',', file=fout)
    print(args.sample, ','.join(obs), sep=',', file=fout)
    fout.close()

#! /usr/bin/env python

import csv
import argparse
import sys
import collections
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
vepcols = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID", "MANE", "TSL", "APPRIS", "CLIN_SIG", "SOMATIC", "PHENO"]
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


# Lineage matching (COG-UK scheme)
mutorder = ['d614g', 'n439k', 'p4715l', 'a222v', 'y453f', 'n501y', 't1001i', 'p681h', 'q27*', 'del_21765_6']
mut = dict()
mut['d614g'] = (23402, 23404)
mut['n439k'] = (22877, 22879)
mut['p4715l'] = (14407, 14409)
mut['a222v'] = (22226, 22228)
mut['y453f'] = (22919, 22921)
mut['n501y'] = (23063, 23065)
mut['t1001i'] = (3266, 3268)
mut['p681h'] = (23603, 23605)
mut['q27*'] = (27972, 27974)
mut['del_21765_6'] = (21764, 21764)  # last aligned base for coverage
if sorted(mut.keys()) != sorted(mutorder):
    print("Mutation disctionaries are inconsisten!", file=sys.stderr)

# Create observed mutation string
with open(args.outfile, 'w') as fout:
    obs = []
    for mt in mutorder:
        (start, end) = mut[mt]
        # Check coverage
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
    print('sample', ','.join(mutorder), sep=',', file=fout)
    print(args.sample, ','.join(obs), sep=',', file=fout)
    fout.close()

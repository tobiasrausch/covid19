#! /usr/bin/env python

from __future__ import print_function
import argparse
import collections
import sys
import numpy
from readfq import readfq


# Parse command line
parser = argparse.ArgumentParser(description='FASTA/MSA filter tool')
parser.add_argument('-m', '--msa', metavar='mafft.msa', required=True, dest='msa', help='MSA file (required)')
parser.add_argument('-n', '--nperc', metavar='0.05', required=False, dest='nperc', help='max. fraction N (optional)')
parser.add_argument('-g', '--gperc', metavar='0.5', required=False, dest='gperc', help='max. fraction gaps & Ns per column (optional)')
parser.add_argument('-l', '--lperc', metavar='0.9', required=False, dest='lperc', help='min. fraction length (optional)')
parser.add_argument('-a', '--amb', metavar='15', required=False, dest='amb', help='max. ambiguous nucleotides (optional)')
parser.add_argument('-d', '--dedup', dest='dedup', action='store_true', default=False, help='remove duplicates')
args = parser.parse_args()

# Parameters
gperc = 0.5
if args.gperc:
    gperc = float(args.gperc)
nperc = 0.05
if args.nperc:
    nperc = float(args.nperc)
lperc = 0.9
if args.lperc:
    lperc = float(args.lperc)
amb = 15
if args.amb:
    amb = int(args.amb)

# Compute stats
slmed = None
fracnmed = None
ambmed = None
sizecut = None
anygaps = False
flagpos = set()
if args.msa:
    for stats in [True, False]:
        seqlen = list()
        fracn = list()
        ambiguous = list()
        gaps = collections.Counter()
        dups = set()

        # Computate stats
        f_in = open(args.msa)
        num = 0
        passcount = 0
        dupcount = 0
        for seqName, seqNuc, seqQuals in readfq(f_in):
            counter = collections.Counter(seqNuc.upper())
            sl = sum(counter.values()) - counter['-']
            acgtn = counter['A'] + counter['C'] + counter['G'] + counter['T'] + counter['N']
            fn = float(counter['N']) / float(sl)
            ab = sl - acgtn
            if stats:
                seqlen.append(sl)
                fracn.append(fn)
                ambiguous.append(ab)
                for idx, c in enumerate(seqNuc.upper()):
                    if (c == 'N') or (c == '-'):
                        if c == '-':
                            anygaps = True # MSA
                        gaps[idx] += 1
            else:
                if fn <= nperc:
                    if ab <= amb:
                        if sl >= sizecut:
                            if anygaps:
                                # Alignment
                                outnuc = ''
                                for idx, c in enumerate(seqNuc):
                                    if idx not in flagpos:
                                        outnuc += c
                                ostr = True
                                if args.dedup:
                                    cleanstr = outnuc.upper().replace('N', '').replace('-', '')
                                    if cleanstr in dups:
                                        ostr = False
                                        dupcount += 1
                                    else:
                                        dups.add(cleanstr)
                                if ostr:
                                    print(">" + seqName)
                                    print(outnuc)
                                    passcount += 1
                            else:
                                # FASTA
                                ostr = True
                                if args.dedup:
                                    cleanstr = seqNuc.upper().replace('N', '').replace('-', '')
                                    if cleanstr in dups:
                                        ostr = False
                                        dupcount += 1
                                    else:
                                        dups.add(cleanstr)
                                if ostr:
                                    print(">" + seqName)
                                    print(seqNuc)
                                    passcount += 1
            num += 1
        if stats:
            slmed = numpy.median(seqlen)
            maxlen = numpy.max(seqlen)
            sizecut = int(lperc * float(maxlen))
            fracnmed = numpy.median(fracn)
            ambmed = numpy.median(ambiguous)
            if anygaps:
                # MSA
                threshold = int(gperc * float(num))
                for k in gaps.keys():
                    if gaps[k] > threshold:
                        flagpos.add(k)
            print("Median sequence length:", slmed, file=sys.stderr)
            print("Maximum sequence length:", maxlen, file=sys.stderr)
            print("Minimum required length:", sizecut, file=sys.stderr)
            print("Median N fraction:", fracnmed, file=sys.stderr)
            print("Median number of ambigous nucleotides:", ambmed, file=sys.stderr)
            print("Gapped FASTA (alignment):", anygaps, file=sys.stderr)
            print("Number of removed alignment columns:", len(flagpos), file=sys.stderr)
        else:
            print("Total sequences:", num, file=sys.stderr)
            print("Duplicate sequences:", dupcount, file=sys.stderr)
            print("Passed sequences:", passcount, file=sys.stderr)
            print("Fraction passed:", float(passcount) / float(num), file=sys.stderr)

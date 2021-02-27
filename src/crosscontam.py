#! /usr/bin/env python

from __future__ import print_function
import cyvcf2
import argparse
import numpy
import collections

# Parse command line
parser = argparse.ArgumentParser(description='Cross-contamination assessment')
parser.add_argument('-b', '--bcf', required=True, dest='bcf', help='multi-sample BCF file (required)')
args = parser.parse_args()

if args.bcf:
    vcf = cyvcf2.VCF(args.bcf)
    samples = numpy.array(vcf.samples)
    vaf = collections.defaultdict(list)
    for record in vcf:
        dp = record.format('DP')
        ad = record.format('AD')
        for idx, s in enumerate(samples):
            vaf[idx].append(float(max(ad[idx])) / float(max(dp[idx])))
    for idx, s in enumerate(samples):
        print(s, str(round(1 - sorted(vaf[idx])[10], 2)), sep="\t")

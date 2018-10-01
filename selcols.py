#!/usr/bin/env python

import os
import sys

# Usage
if len(sys.argv) is not 2:
    exit("Usage:\n\t" + sys.argv[0] + "colname_1,colname_2,colname_3 <input.tsv> output.tsv\n\nInput was:\n\tA tab-separated values (TSV) file")

# argv
c = sys.argv[1] # column name, e.g. 'participant_id,Mutect2_passed_combined_vcf,Mutect2_passed_combined_vcf_index'

c = c.split(',')

# main
l = sys.stdin.readline()
l = l.strip('\n').split('\t')
h = dict((map(reversed, enumerate(l)))) # column name:indices
ci = [h[cx] for cx in c] # column name indices

lc = [l[cix] for cix in ci]
sys.stdout.write('\t'.join(lc) + '\n') # column name indices

for l in sys.stdin:
    l = l.strip('\n').split('\t')
    lc = [l[cix] for cix in ci]
    sys.stdout.write('\t'.join(lc) + '\n') # column name indices

### The End ###

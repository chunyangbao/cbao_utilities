#!/usr/bin/env python
## row2rows.py Ensembl2uniprot.txt 1 ','

import os
import sys

def r2rcs(l,c,s): # split one Row to Rows based on C and S 
    # ['ENSRNOT00000016175,ENSRNOT00000046192', 'P32215', 'Adcyap1r1'] -> ['ENSRNOT00000016175\tP32215\tAdcyap1r1', 'ENSRNOT00000046192\tP32215\tAdcyap1r1']
    ll = l.strip('\n').split('\t') # list of line
    l1c = ll[c].split(s) # list of col
    l0c = ['\t'.join(ll[:c] + ll[(c + 1):])] # list of others
    lc = ['\t'.join([x,y]) for x,y in [(x,y) for x in l1c for y in l0c]] # list of combined col
    return '\n'.join(lc) + '\n'

usage = """
    %prog $filename $col
"""

# argv
fi = sys.argv[1] # filename, fi = './snvs_sum/IN2CA/fmtas_snps1CA1loh.txt'
ci = int(sys.argv[2])-1 # col, ci = int('16')-1
si = sys.argv[3] # sep, si = ','
# main
if fi=="-":
    for li in sys.stdin:
        sys.stdout.write(r2rcs(li,ci,si))
else:
    f = open(fi,'r')
    for li in f:
        sys.stdout.write(r2rcs(li,ci,si))

if fi!="-":
    f.close()

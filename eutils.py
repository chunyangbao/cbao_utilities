#!/usr/bin/env python

import os
import sys
import re
import argparse
import urllib2

su1 = 'https://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=' # Suffix_Url
su2 = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=' # Suffix_Url
idl = [] # ID_list

def parse_args(arg_lst):
    parser = argparse.ArgumentParser()
    parser.add_argument('I', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('o', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument("-d", "--database", metavar='STR', type=str, default='nuccore+nucest', help="NCBI database ID")
    parser.add_argument("-t", "--rettype", metavar='STR', type=str, default='fasta', help="Return type")

    return parser.parse_args(arg_lst)

def main(argv):
    args = parse_args(argv)
    I = args.I
    dats = args.database.split('+')
    t = args.rettype
    for i, d in [(i,d) for i in I for d in dats]:
        try:
            i = i.strip()
            u1 = su1 + d + '&term=' + i
            id = re.findall('<Id>(.+?)</Id>', urllib2.urlopen(u1).read() ,re.DOTALL) # List_Pathway
            id = list(set(id) - set(idl))
            idl.extend(id)
            if len(id) == 0:
                continue
            
            for x in id:
                u2 = su2 + d + '&id=' + x + '&rettype=' + t
                o = urllib2.urlopen(u2).read().strip() + '\n'
                sys.stdout.write(o)
                
        except KeyboardInterrupt:
            pass

if __name__ == "__main__":
    main(sys.argv[1:])


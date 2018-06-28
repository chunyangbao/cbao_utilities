#!/usr/bin/env python
## ant2x.py entrezgene.txt 10116

import mygene
from Bio import Entrez
from bioservices.kegg import KEGG
import os
import sys
import re

def ld2dl(ld): # List_Dict to Dict_List
    #[{'k1': 'v11', 'k2': 'v21'}, {'k1': 'v12', 'k2': 'v22'}] -> {'k1': ['v11', 'v12'], 'k2': ['v21', 'v22']}
    ld = [ld] if type(ld)==dict else ld # Dict to List_Dict
    lk = [] # List_Key
    for dx in ld:
        lk.extend(dx.keys())
    
    lk = sorted(set(lk),key=lk.index) # unique List_Key
    dl = {}.fromkeys(lk,[])
    for kx in lk:
        lx = []
        for dy in ld:
            lx.append(dy[kx]) if dy.has_key(kx) else 'NA' # Dict(k) -> List
        
        dl[kx] = lx

    return dl

usage = """
    %prog $Filename $Taxonomy
"""
# argv
fi = sys.argv[1] # Filename, fi = 'entrezgene.txt'
ti = sys.argv[2] # Taxonomy, ti = '10116'
# load
Entrez.email = "cybao.cpu@163.com"
handle = Entrez.efetch(db="Taxonomy", id=ti, retmode="xml")
records = Entrez.read(handle)
si =records[0]["ScientificName"] # Specics_Id, si = 'Rattus norvegicus'
mg = mygene.MyGeneInfo()
kg = KEGG()
kt = str(kg.lookfor_organism(si)).replace("u'", "").split(" ")[1] # KEGG_Taxonomy
kge = kg.list(kt)
# init
lq = [] # List_Query
if fi=="-":
    for ln in sys.stdin:
        lq.append(ln.strip('\n'))
else:
    f = open(fi,'r')
    for ln in f:
        lq.append(ln.strip('\n'))

lf = ['entrezgene', 'ensembl.gene', 'symbol', 'name', 'alias', 'summary', 'refseq', 'unigene', 'ensembl.transcript', 'ensembl.protein', 'uniprot', 'interpro', 'go'] # List_Field
lr = ['query','ensembl_gene', 'ensembl_transcript', 'ensembl_protein', 'entrezgene', 'symbol', 'name', 'alias', 'unigene', 'refseq_genomic', 'refseq_rna', 'refseq_protein', 'uniprot', 'interpro', 'go.cc', 'go.mf', 'go.bp', 'kegg', 'summary'] # List_Return
mq = mg.querymany(lq, scopes=['entrezgene', 'ensembl.gene', 'refseq', 'unigene'], fields=lf, species=ti) # Mygene_Query
sys.stderr.write("\t".join([str(h) for h in lr])+"\n") # header
# main
for dq in mq: # Dict_Query dq = mq[0]
    dr = {} # Dict_Return
    # id and sum
    dr['query'] = str(dq['query']).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "")
    dr['entrezgene'] = str(dq['entrezgene']).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq.has_key('entrezgene') else 'NA'
    dr['ensembl_gene'] = str(dq['ensembl.gene']).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq.has_key('ensembl.gene') else 'NA'
    dr['ensembl_transcript'] = str(dq['ensembl.transcript']).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq.has_key('ensembl.transcript') else 'NA'
    dr['ensembl_protein'] = str(dq['ensembl.protein']).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq.has_key('ensembl.protein') else 'NA'
    dr['symbol'] = str(dq['symbol']).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq.has_key('symbol') else 'NA'
    dr['name'] = str(dq['name']).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq.has_key('name') else 'NA'
    dr['alias'] = str(dq['alias']).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq.has_key('alias') else 'NA'
    dr['summary'] =  str(dq['summary']).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq.has_key('summary') else 'NA'
    dr['unigene'] = str(dq['unigene']).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq.has_key('unigene') else 'NA'
    dr['uniprot'] = str(dq['uniprot']).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq.has_key('uniprot') else 'NA'
    dr['interpro'] = str(dq['interpro']).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq.has_key('interpro') else 'NA'
    # refseq
    if dq.has_key('refseq'):
        dr['refseq_genomic'] =  str(dq['refseq']['genomic']).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq['refseq'].has_key('genomic') else 'NA'
        dr['refseq_rna'] =  str(dq['refseq']['rna']).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq['refseq'].has_key('rna') else 'NA'
        dr['refseq_protein'] =  str(dq['refseq']['protein']).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq['refseq'].has_key('protein') else 'NA'
    else:
        dr['refseq_genomic'] =  'NA'
        dr['refseq_rna'] =  'NA'
        dr['refseq_protein'] =  'NA'
    # go
    if dq.has_key('go'):
        dr['go.cc'] = str(ld2dl(dq['go']['CC'])).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq['go'].has_key('CC') else 'NA'
        dr['go.mf'] = str(ld2dl(dq['go']['MF'])).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq['go'].has_key('MF') else 'NA'
        dr['go.bp'] = str(ld2dl(dq['go']['BP'])).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "") if dq['go'].has_key('BP') else 'NA'
    else:
        dr['go.cc'] = 'NA'
        dr['go.mf'] = 'NA'
        dr['go.bp'] = 'NA'
    # kegg
    if len(re.findall(kt+':'+dr['entrezgene']+"\t", kge))!=0:
        dr['kegg'] = str(kg.get_pathway_by_gene(dr['entrezgene'], kt)).replace('u"', "").replace('"', "").replace("u'", "").replace("'", "")
    else:
        dr['kegg'] = 'NA'
    
    sys.stderr.write("\t".join([dr[kr] for kr in lr])+"\n") # write

if fi!="-":
    f.close()

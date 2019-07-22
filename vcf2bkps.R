#!/usr/bin/env Rscript

library(optparse)
library(data.table)

V="Version: 1.0"
D="Depends: R (>= 3.5.0), optparse, data.table"

a = commandArgs(trailingOnly = TRUE)

option_list = list(
    make_option(c("-p", "--prune_colnames"), type="logical", default="FALSE", action = "store_true",
     help="Should colnames be pruned [%default].", metavar="logical")
)

opt_parser = OptionParser(usage = "usage: %prog [options] <sv.vcf>\n\t
    This required input file should be in VCF format (separator: '\\t', stdin: '-')", 
    option_list=option_list, description = paste(V, D, sep = "\n"))
opt = parse_args(opt_parser, args = a, positional_arguments = TRUE)

p = opt$options$prune_colnames

n = opt$args[1]

# main
## load packages
library(VariantAnnotation)

## read
### fread
sv_vcf = readVcf(n)
### formating
ad_mat = geno(sv_vcf)$AD
if (p) colnames(ad_mat) = gsub('^.*/', '', colnames(ad_mat))
if (p) colnames(ad_mat) = gsub('\\.bam$', '', colnames(ad_mat))
adt_mat = ad_mat[, -1, drop = FALSE]
adn_mat = ad_mat[, 1, drop = FALSE]

sr_mat = geno(sv_vcf)$SR
if (p) colnames(sr_mat) = gsub('^.*/', '', colnames(sr_mat))
if (p) colnames(sr_mat) = gsub('\\.bam$', '', colnames(sr_mat))
srt_mat = sr_mat[, -1, drop = FALSE]

## process
### build DT
sv_dt = data.table(
    ID = names(sv_vcf), 
    SVidx = gsub(':.*$', '', names(sv_vcf)), 
    Sample = colnames(adt_mat)[apply(adt_mat, 1, which.max)], 
    Group = colnames(adt_mat)[apply(adt_mat, 1, which.max)], 
    chr1 = as.character(seqnames(sv_vcf)), 
    pos1 = start(sv_vcf), 
    str1 = ifelse(grepl("^\\[", as.character(alt(sv_vcf))) | grepl("^\\]", as.character(alt(sv_vcf))), '-1', '1'), 
    TotalCount = rowSums(adt_mat), 
    SplitCount = rowSums(srt_mat), 
    maq = info(sv_vcf)$MAPQ, 
    suppSamples = rowSums(adt_mat > 0), 
    TCount = rowSums(adt_mat), 
    NCount = rowSums(adn_mat), 
    adt_mat, 
    adn_mat, 
    key = 'SVidx'
)
### split the DT to ':1' and ':2'
sv1_dt = sv_dt[ID %like% ":1$"]
sv2_dt = sv_dt[ID %like% ":2$"]
### merge 2 DTs
sv2_dt[, `:=`(chr2 = chr1, pos2 = pos1, str2 = str1), ]
sv2_dt = sv2_dt[, .(SVidx, chr2, pos2, str2)]
o_dt = merge(sv1_dt, sv2_dt)
### reorder columns
o_dt[, ID := NULL]
setcolorder(o_dt, c('SVidx', 'Sample', 'Group', 'chr1', 'pos1', 'str1', 'chr2', 'pos2', 'str2'))
### reorder rows
o_dt = o_dt[order(as.character(chr1), as.integer(pos1))]

## output
### fwrite
fwrite(o_dt, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)



# test at checkpoint
#n='/cga/bass/Chunyang/task/Matthew/WGS_EAC/SvABA_ms/filtered_somatic_sv1g_vcf/EAC-9_11.svaba.somatic.sv1g.vcf'

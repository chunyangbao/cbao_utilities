#!/usr/bin/env Rscript

library("optparse")

V="Version: 1.3"
D="Depends: R (>= 3.1.0), optparse, data.table, maftools"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-o", "--output_prefix"), type="character", default="oncoplot",
     help="Output prefix [%default].", metavar="character"),
    make_option(c("-b", "--feature_table"), type="character", default="NA",
     help="Table of features. The header of first column should be \n\t\t'Tumor_Sample_Barcode'. Could be a text file or a data.frame [%default].", metavar="character"),
    make_option(c("-l", "--feature_col"), type="character", default="NA",
     help="Column(s) in features table separated by ','. [%default].", metavar="character"),
    make_option(c("--sort_mutation"), type="logical", default="FALSE", action = "store_true",
     help="Force sort matrix according mutations. Helpful in case of MAF was read \n\t\talong with copy number data. [%default]", metavar="logical"),
    make_option(c("--sort_annotation"), type="logical", default="FALSE", action = "store_true",
     help="Logical sort oncomatrix (samples) by provided 'clinicalFeatures'. \n\t\tSorts based on first 'clinicalFeatures'. [%default]", metavar="logical")
)

opt_parser = OptionParser(usage="usage: %prog [options] <*.maf>\n\tBy default <*.maf> is '-' (stdin)", option_list=option_list, description = paste(V, D, sep="\n"))
opt = parse_args(opt_parser, args=a, positional_arguments=TRUE)

o = opt$options$output_prefix
b = opt$options$feature_table
l = unlist(strsplit(opt$options$feature_col, ",", fixed=TRUE))
sm = opt$options$sort_mutation
sa = opt$options$sort_annotation

m = opt$args[1] # a MAF file
if (is.na(m)|m=='-') m = 'file:///dev/stdin'

if (b=='NA') b=NULL
if (l[1]=='NA') l=NULL

# motif matrix
.libPaths(c("~/lib", .libPaths()))
library('maftools')

all_maf = read.maf(maf = m, clinicalData = b)

# oncoplot
cairo_pdf(filename = paste0(o, ".oncoplot.pdf"), width = 7, height = 10)
    oncoplot(maf = all_maf, clinicalFeatures = l, top =20, fontSize = 7, showTumorSampleBarcodes = TRUE, sortByMutation = sm, sortByAnnotation = sa)
dev.off()

cairo_pdf(filename = paste0(o, ".oncoplot_all.pdf"), width = 7, height = 10)
    oncoplot(maf = all_maf, clinicalFeatures = l, top =nrow(getGeneSummary(all_maf)), fontSize = 7, showTumorSampleBarcodes = TRUE, sortByMutation = sm, sortByAnnotation = sa, writeMatrix = TRUE)
    om_df <- read.delim('onco_matrix.txt', check.names=FALSE)
    write.table(t(c('gene_symbol', colnames(om_df))), paste0(o, ".oncomatrix_all.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(om_df, paste0(o, ".oncomatrix_all.tsv"), sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE, append = TRUE)
    file.remove('onco_matrix.txt')
dev.off()

### THE END ###

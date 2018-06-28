#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.0"
D="Depends: R (>= 3.4.0), optparse, data.table"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("--output_prefix"), type="character", default="annovarToMaf",
     help="If provided writes resulting MAF file to an output file [%default].", metavar="character"),
    make_option(c("--center"), type="character", default="NA",
     help="Center field in MAF file will be filled with this value [%default].", metavar="character"),
    make_option(c("--refBuild"), type="character", default="hg19",
     help="NCBI_Build field in MAF file will be filled with this value [%default].", metavar="character"),
    make_option(c("--tsbCol"), type="character", default="Tumor_Sample_Barcode",
     help="column name containing Tumor_Sample_Barcode or sample names in input file [%default].", metavar="character"),
    make_option(c("--ref_table"), type="character", default="refGene",
     help="reference table used for gene-based annotations. Can be 'ensGene' or 'refGene' [%default].", metavar="character")
)

opt_parser <- OptionParser(usage="usage: %prog [options] <input.tsv> \n\t<input.tsv> is a TSV file (separator: '\\t', stdin: '-').", option_list=option_list, description = paste(V, D, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

output_prefix <- opt$options$output_prefix
center <- opt$options$center
refBuild <- opt$options$refBuild
tsbCol <- opt$options$tsbCol
ref_table <- opt$options$ref_table

d <- ifelse(opt$args[1]=='-', 'file:///dev/stdin', opt$args[1]) # Data_path

# main
.libPaths(c("~/lib", .libPaths()))
library("maftools")

ann_maf = annovarToMaf(annovar = d, Center = center, refBuild = refBuild, 
                               tsbCol = tsbCol, table = ref_table, basename = output_prefix)

### THE END ###

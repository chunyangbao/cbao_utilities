#!/usr/bin/env Rscript

library("optparse")

V="Version: 1.1"
D="Depends: R (>= 3.1.0), optparse, data.table, maftools"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-o", "--dir_out"), type="character", default=".",
     help="Output directory [%default] \n\t\tBy default <o> is '.' (current directory).", metavar="character")
)

opt_parser = OptionParser(usage="usage: %prog [options] <*.maf>\n\tBy default <*.maf> is '-' (stdin)", option_list=option_list, description = paste(V, D, sep="\n"))
opt = parse_args(opt_parser, args=a, positional_arguments=TRUE)

o = opt$options$dir_out
o = paste0(sub("/$", "", o), "/")
if (!dir.exists(o)) dir.create(o)

m = opt$args[1] # a MAF file
if (is.na(m)|m=='-') m = 'file:///dev/stdin'

# motif matrix
library('maftools')

snv_maf = read.maf(maf = m)

## pair
for (s_id in levels(getSampleSummary(snv_maf)[[1]])) {
    # dir.create
    s_o <- paste0(o, s_id, "/")
    if (!dir.exists(s_o)) dir.create(s_o)
    setwd(s_o)
    # plot
    pdf(paste0(s_o, s_id, '.maftools_rainfallPlot.pdf'), paper = 'a4')
        rainfallPlot(maf = snv_maf, tsb = s_id, detectChangePoints = TRUE, fontSize = 12, pointSize = 0.6)
    dev.off()
    }

### THE END ###

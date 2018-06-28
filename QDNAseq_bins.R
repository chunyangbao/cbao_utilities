#!/usr/bin/env Rscript

library("optparse")

V="Version: 2.0"
L="Libraries: R (>= 3.1.0), optparse, QDNAseq, future, QDNAseq.mm10 (optional)"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-p", "--processors"), type="integer", default="1",
     help="Number of processors for parallell computing [%default].", metavar="integer"),
    make_option(c("-o", "--output"), type="character", default=".",
     help="Output directory [%default].", metavar="character"),
    make_option(c("-b", "--bin"), type="integer", default="15",
     help="Width of the bins in units of kbp (1000 base pairs) [%default].", metavar="integer"),
    make_option(c("-G", "--genome"), type="character", default="hg19",
     help="Genome version, e.g. 'hg19' or 'mm10' [%default].", metavar="character"),
    make_option(c("-y", "--seqType"), type="character", default="SR50",
     help="Experiment type, e.g. 'SR50' or 'PE100' [%default].", metavar="character"),
    make_option(c("-l", "--blacklist"), type="character", default="NA",
     help="Blacklisted Regions in a BED file [%default].", metavar="character")
)

opt_parser <- OptionParser(usage="usage: %prog [options]", option_list=option_list, description = paste(V, L, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

p <- opt$options$processors
o <- opt$options$output
b <- opt$options$bin
g <- opt$options$genome
y <- opt$options$seqType
l <- opt$options$blacklist

co <- ifelse(g %in% c("hg19"), "human", "other") # callBins_Organism
# Main
## Library
library(QDNAseq)
library(CGHcall)
library(future)
if (g == "mm10") library(QDNAseq.mm10)
## Parallel processing on the current machine
future::plan("multicore")
options(mc.cores=p)
## Bin annotations
bins <- getBinAnnotations(binSize=b, genome=g, type=y) # Bin annotations
if (l != "NA") {
    bins_ta <- data.frame(chromosome=bins$chromosome, start=bins$start, end=bins$end)
    bins$blacklist <- calculateBlacklist(bins_ta, bedFiles=c(l))
    }
# Output
if (!dir.exists(o)) dir.create(o)
setwd(o)
save(bins, file=paste0(paste("QDNAseq","bins",b,g,y,sep="_"),".RData"))
### THE END ###

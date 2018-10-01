#!/usr/bin/env Rscript

library("optparse")

V="Version: 2.0"
L="Libraries: R (>= 3.1.0), optparse, QDNAseq, future"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-p", "--processors"), type="integer", default="1",
     help="Number of processors for parallell computing [%default].", metavar="integer"),
    make_option(c("-o", "--output"), type="character", default="bin.RData",
     help="Output file [%default].", metavar="character"),
    make_option(c("-b", "--bin"), type="integer", default="10",
     help="Width of the bins in units of kbp (1000 base pairs) [%default].", metavar="integer"),
    make_option(c("-s", "--bsg"), type="character", default="BSgenome.Hsapiens.UCSC.hg19",
     help="Package name of 'BSgenome' [%default]", metavar="character"),
    make_option(c("-w", "--bigwig"), type="character", default="NA",
     help="Mappability file in the bigWig format' [%default]", metavar="character"),
    make_option(c("-e", "--bigWigAverageOverBed"), type="character", default="NA",
     help="bigWigAverageOverBed binary [%default]", metavar="character"),
    make_option(c("-l", "--blacklist"), type="character", default="NA",
     help="Blacklisted Regions in a BED file [%default].", metavar="character")
)

opt_parser <- OptionParser(usage="usage: %prog [options]", option_list=option_list, description = paste(V, L, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

p <- opt$options$processors
o <- opt$options$output
b <- opt$options$bin
s <- opt$options$bsg
w <- opt$options$bigwig
e <- opt$options$bigWigAverageOverBed
l <- opt$options$blacklist

# Main
## Library
library(QDNAseq)
library(Biobase)
library(s, character.only=T)
## Parallel processing on the current machine
future::plan("multicore")
options(mc.cores=p)
## createBins
args_createBins <- paste0("createBins(bsgenome=", s, ", binSize=", b, ")") # Load the sequencing data from bam file list.
bins <- eval(parse(text=args_createBins)) # binReadCounts
if(w != "NA" & e != "NA") bins$mappability <- calculateMappability(bins, bigWigFile=w, bigWigAverageOverBed=e)
if(l != "NA") bins$blacklist <- calculateBlacklist(bins, bedFiles=c(l))
## output
save(bins, file=o)

### THE END ###

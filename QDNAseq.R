#!/usr/bin/env Rscript

library("optparse")

V="Version: 4.1"
L="Libraries: R (>= 3.4.0), optparse, QDNAseq, future, QDNAseq.mm10 (optional)"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-p", "--processors"), type="integer", default="1",
     help="Number of processors for parallell computing [%default].", metavar="integer"),
    make_option(c("-o", "--output"), type="character", default=".",
     help="Output directory [%default].", metavar="character"),
    make_option(c("-e", "--extra_args_binReadCounts"), type="character", default="isUnmappedQuery=FALSE",
     help="Extra arguments for binReadCounts [%default].", metavar="character"),
    make_option(c("-C", "--CHR"), type="character", default="X,Y",
     help="Chromosomes to be filtered out, e.g. 'X' or 'Y' or 'NA' [%default].", metavar="character"),
    make_option(c("-f", "--transformFun"), type="character", default="sqrt",
     help="A function to transform the data with, e.g. 'sqrt' or 'log2' or 'none' [%default].", metavar="character"),
    make_option(c("-s", "--sampleType"), type="character", default="T",
     help="Sample type, e.g. 'T' or 'N' [%default].", metavar="character"),
    make_option(c("-P", "--plots"), type="logical", default="FALSE", action = "store_true",
     help="Should the plots be generated [%default].", metavar="logical")
)

opt_parser <- OptionParser(usage="usage: %prog [options] <bam> <bin>\n\t<bam> should be either a sample file or a path containing BAM files. The sample file is a TSV file for bam file of each sample (row: sample, separator: '\\t', stdin: '-') containing 'sample_id' and 'read_file' in two separate columns. The possible value of 'read_file' column: '/path/to/your.bam'. <bin> should be the directory of bin file", option_list=option_list, description = paste(V, L, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

p <- opt$options$processors
o <- opt$options$output
e <- opt$options$extra_args_binReadCounts
CHR <- opt$options$CHR
f <- opt$options$transformFun
s <- opt$options$sampleType
P <- opt$options$plots

bam <- ifelse(opt$args[1]=='-', 'file:///dev/stdin', opt$args[1]) # bam
bin <- opt$args[2] # bin

CHR <- unlist(strsplit(CHR,","))
if (s != "T" & s != "N") {cat("'s' must be one of 'T' or 'N' \n"); q('no')}
# Main
## Library
library(QDNAseq)
library(CGHcall)
library(future)
## Parallel processing on the current machine
future::plan("multicore")
options(mc.cores=p)
## Bin annotations
load(bin)
g <- attr(bins,'QDNAseq')$build
co <- ifelse(g %in% c("hg19"), "human", "other") # callBins_Organism
## Processing bam file
if (file_test("-f", bam)) {
    bam <- read.delim(bam, check.names=FALSE, stringsAsFactors=FALSE)
    args_binReadCounts <- paste0("binReadCounts(bins, bamfiles=as.character(bam$read_file), ", e, ")") # Load the sequencing data from bam file list.
} else if (file_test("-d", bam)) {
    args_binReadCounts <- paste0("binReadCounts(bins, path=bam, ", e, ")") # Load the sequencing data from bam file list.
} else {
    write("Input should be either a BAM file or a path containing BAM files.", stderr())
    q(save='no')
}

readCounts <- eval(parse(text=args_binReadCounts)) # binReadCounts

## readCounts
readCounts_af <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
readCounts_ec <- estimateCorrection(readCounts_af) # Estimate the correction for GC content and mappability
readCounts_2 <- readCounts_ec
if (sum(CHR==c("X","Y")) < 2) readCounts_2 <- applyFilters(readCounts_ec, chromosomes=C) # Sex chromosomes
## bin
copyNumbers <- correctBins(readCounts_2) # Apply the correction for GC content and mappability.
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
## Downstream analyses
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun=f) # Segmentation
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented) # Tune segmentation parameters
if (s == "T") copyNumbersCalled <- callBins(copyNumbersSegmented, organism=co, ncpus=p)
# if (s == "N") copyNumbersCalled <- callBins(copyNumbersSegmented, organism=co, method="cutoff", ncpus=p)
cgh <- makeCgh(copyNumbersCalled) # for cghCall
# Output
if (!dir.exists(o)) dir.create(o)
setwd(o)
save(cgh, file="cgh.RData")
## Files
exportBins(copyNumbersSmooth, file="copyNumbersSmooth.tsv", format="tsv")
exportBins(copyNumbersSmooth, file="copyNumbersSmooth.igv", format="igv")
exportBins(copyNumbersSmooth, file="copyNumbersSmooth.bed", format="bed")
exportBins(cgh, file="cgh.segments.tsv", type="segments")
exportBins(cgh, file="cgh.copynumber.tsv", type="copynumber")
exportBins(cgh, file="cgh.calls.tsv", type="calls")
exportBins(cgh, format="vcf")
exportBins(cgh, format="seg")
## Plots
if (P) {
    pdf('plots.pdf')
    plot(readCounts, logTransform=FALSE); highlightFilters(readCounts, logTransform=FALSE, residual=TRUE, blacklist=TRUE)
    isobarPlot(readCounts_af)
    noisePlot(readCounts_ec)
    plot(copyNumbersSmooth)
    plot(copyNumbersSegmented)
    plot(copyNumbersCalled,delcol="darkblue",losscol="blue",gaincol="red",ampcol="darkred")
    dev.off()
}
### THE END ###


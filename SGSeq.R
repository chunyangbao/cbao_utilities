#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.0"
D="Depends: R (>= 3.1.0), optparse, data.table, SGSeq, txdb"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-o", "--dir_out"), type="character", default=".",
     help="Output directory [%default] \n\t\tBy default <o> is '.' (current directory).", metavar="character"),
    make_option(c("-p", "--processors"), type="integer", default="1",
     help="Number of processors for parallell computing [%default]", metavar="integer"),
    # SGSeq
    make_option(c("-d", "--txdb"), type="character", default="TxDb.Hsapiens.UCSC.hg38.knownGene",
     help="Package name of 'TxDb' [%default]", metavar="character"),
    make_option(c("-s", "--txdb_sq"), type="character", default="/home/ref/gtf/gencode_22/gencode_22.sqlite",
     help="Directory of 'sqlite' [%default]", metavar="character"),
    make_option(c("-f", "--txf"), type="logical", default="FALSE", action = "store_true",
     help="In 'predictTxFeatures', should the features be generated? [%default]", metavar="logical"),
    make_option(c("-r", "--ref"), type="logical", default="FALSE", action = "store_true",
     help="In 'analyzeFeatures', should the analysis based on reference? [%default]", metavar="logical"),
    # fread
    make_option(c("--col_sep"), type="character", default="\t",
     help="Separator between columns [\\t] \n\t\t'auto' represent for '[,\\t |;:]'. See also 'data.table::fread'.", metavar="character"),
    make_option(c("--row_num"), type="integer", default="-1",
     help="Number of rows to read [%default] \n\t\t'-1' means all. '0' is a special case that just returns the column names.", metavar="integer"),
     make_option(c("--col_name"), type="character", default="auto",
     help="The first data line contain column names, 'auto'|TRUE|FALSE [%default] \n\t\tDefaults according to whether every non-empty field on the first data\n\t\tline is type character.", metavar="character"),
    make_option(c("--na_str"), type="character", default="NA",
     help="String(s) which are to be interpreted as NA values [%default]", metavar="character"),
    make_option(c("--row_skip"), type="character", default="0",
     help="Row can be taken as the first data row [%default] \n\t\tskip='string' searches for 'string' in the file and starts on that row.", metavar="character"),
    make_option(c("--row_name"), type="character", default="1",
     help="Column number (or name) can be taken as the row names [%default]", metavar="character")
)

opt_parser = OptionParser(usage="usage: %prog [options] <samples.tsv>\n\tBy default <samples.tsv> is '-' (stdin). This matrix should contain \n\t'sample_name' and 'file_bam' in two separate columns.", option_list=option_list, description = paste(V, D, sep="\n"))
opt = parse_args(opt_parser, args=a, positional_arguments=TRUE)

o = opt$options$dir_out
p = opt$options$processors

d = opt$options$txdb
s = opt$options$txdb_sq
f = opt$options$txf
r = opt$options$ref

cs = opt$options$col_sep
rn = opt$options$row_num
cn = as.logical(opt$options$col_name)
na = opt$options$na_str
rs = opt$options$row_skip
ra = opt$options$row_name

x = opt$args[1]
if(is.na(x)|x=='-') x = 'file:///dev/stdin'

if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)
if(grepl("^[[:digit:]]+$", ra)) ra=as.numeric(ra)
# read
xm = fread(x, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE, check.names=FALSE)
xm <- xm[, c("sample_name", "file_bam")]

# main
library(SGSeq)
library(d, character.only=T)

xm <- getBamInfo(xm, cores=p)

txdb <- loadDb(s)
txdbf <- convertToTxFeatures(txdb)

if (f) txf <- predictTxFeatures(xm, min_overhang=NULL, cores=p)

if (r) {
    sgfc_pred <- analyzeFeatures(xm, features = txdbf, cores=p)
} else {
    sgfc_pred <- analyzeFeatures(xm, cores=p)
    sgfc_pred <- annotate(sgfc_pred, txdbf)
}

sgvc_pred <- analyzeVariants(sgfc_pred, cores=p)

# Annotation Table
at <- as.data.frame(mcols(sgvc_pred))
for(i in colnames(at)) if(!is.vector(at[,i])) at[,i] <- sapply(at[,i], simplify=TRUE, toString)
# Usage Table
ut <- variantFreq(sgvc_pred)
colnames(ut) <- paste0(colnames(ut), ".usage")
# Count Table
sgv <- rowRanges(sgvc_pred)
sgvc <- getSGVariantCounts(sgv, sample_info = xm, cores=p)
ct <- counts(sgvc)
colnames(ct) <- paste0(colnames(ct), ".count")
# rEads Table
et <- t(t(ct)/(xm$lib_size/1000000)) # Reads Per Billion Mapped
colnames(et) <- gsub(".count", ".RPBM", colnames(et))
# Variant Table
vt <- at[, "variantID", drop = FALSE]
at <- cbind(vt, at[, -12])
ut <- cbind(vt, ut)
ct <- cbind(vt, ct)
et <- cbind(vt, et)
# output dir
if (!dir.exists(o)) dir.create(o)
setwd(o)
# output file
if (f) save(txf, file="txf.RData")
write.table(at, 'annotation.tsv', sep='\t', row.names=F, col.names=T, quote=F)
write.table(ut, 'usage.tsv', sep='\t', row.names=F, col.names=T, quote=F)
write.table(ct, 'count.tsv', sep='\t', row.names=F, col.names=T, quote=F)
write.table(et, 'RPBM.tsv', sep='\t', row.names=F, col.names=T, quote=F)

### THE END ###

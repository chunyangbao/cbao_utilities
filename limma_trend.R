#!/usr/bin/env Rscript

library("optparse")
library("data.table")
library("methods")

V="Version: 1.5"
D="Depends: R (>= 3.1.0), optparse, data.table, limma, statmod, Biobase"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-o", "--output"), type="character", default=".",
     help="Prefix for all output files [%default].", metavar="character"),
    make_option(c("-f", "--feature_table"), type="character", default="NA",
     help="Feature table. (row: feature; column: annotation, separator: '\\t') [%default].", metavar="character"),
    make_option(c("--filter"), type="character", default="NA",
     help="Filter for values. '1,3' represents 'keep features that \n\t\thave more than 1 in at least 3 samples.' [%default].", metavar="character"),
    make_option(c("--log2"), type="logical", default="FALSE", action = "store_true",
     help="Should the input matrix be log2-transformed? [%default]", metavar="logical"),
    make_option(c("-m", "--fit_model"), type="character", default="~s$g1",
     help="Design formula. 's' represents <sample.tsv> [%default].", metavar="character"),
    make_option(c("-r", "--fit_robust"), type="logical", default="FALSE", action = "store_true",
     help="Should the estimation of 'df.prior' and 'var.prior' be robustified \n\t\tagainst outlier sample variances? [%default]", metavar="logical"),
    make_option(c("-d", "--fit_duplicate"), type="character", default="NA",
     help="Column name corresponding to block in duplicateCorrelation [%default].", metavar="character"),
    make_option(c("-e", "--fit_coefficient"), type="character", default="NULL",
     help="Column name corresponding to coefficient in topTable and topSplice [%default].", metavar="character"),
    make_option(c("--splice_id"), type="character", default="NA",
     help="Column name corresponding to geneid and exonid in topSplice (geneid,exonid) [%default].", metavar="character"),
    make_option(c("--splice_test"), type="character", default="simes",
     help="Test in topSplice [%default].", metavar="character"),
    make_option(c("-p", "--plots"), type="logical", default="FALSE", action = "store_true",
     help="Should the plots be generated [%default].", metavar="logical"),
    # fread
    make_option(c("--col_sep"), type="character", default="\t",
     help="Separator between columns [\\t] \n\t\t'auto' represents '[,\\t |;:]'. See also 'data.table::fread'.", metavar="character"),
    make_option(c("--row_num"), type="integer", default="-1",
     help="Number of rows to read [%default] \n\t\t'-1' means all. '0' is a special case that just returns the column names.", metavar="integer"),
     make_option(c("--col_name"), type="character", default="TRUE",
     help="The first data line contain column names, 'auto'|TRUE|FALSE [%default] \n\t\tDefaults according to whether every non-empty field on the first data\n\t\tline is type character.", metavar="character"),
    make_option(c("--na_str"), type="character", default="NA",
     help="String(s) which are to be interpreted as NA values [%default]", metavar="character"),
    make_option(c("--row_skip"), type="character", default="0",
     help="Row can be taken as the first data row [%default] \n\t\tskip='string' searches for 'string' in the file and starts on that row.", metavar="character"),
    make_option(c("--row_name"), type="character", default="1",
     help="Column number (or name) can be taken as the row names [%default]", metavar="character")
)

opt_parser <- OptionParser(usage="usage: %prog [options] <exp.tsv> <sample.tsv>\n\t<exp.tsv> is a TSV file for expression (row: feature, column: sample, separator: '\\t', stdin: '-').\n\t<sample.tsv> is a TSV file for sample annotation (row: sample, column: annotation, separator: '\\t').", option_list=option_list, description = paste(V, D, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

o <- opt$options$output
f <- opt$options$feature_table
l <- opt$options$filter
lt2 <- opt$options$log2
m <- opt$options$fit_model
r <- opt$options$fit_robust
d <- opt$options$fit_duplicate
e <- opt$options$fit_coefficient
ds <- opt$options$splice_id
ss <- opt$options$splice_test
p <- opt$options$plots

cs <- opt$options$col_sep
rn <- opt$options$row_num
cn <- as.logical(opt$options$col_name)
na <- opt$options$na_str
rs <- opt$options$row_skip
ra <- opt$options$row_name

n <- ifelse(opt$args[1]=='-', 'file:///dev/stdin', opt$args[1])
s <- opt$args[2]
if (l != "NA") l <- as.numeric(unlist(strsplit(l,",")))
ds <- unlist(strsplit(ds,","))

if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)
if(grepl("^[[:digit:]]+$", ra)) ra=as.numeric(ra)
# Input_Matrix
nm <- fread(n, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE)
rownames(nm) <- nm[, ra]
nm <- nm[, -ra]
nm <- as.matrix(nm)
if(lt2) nm <- log2(nm + 1)
# Sample_Table
st <- read.delim(s, check.names=FALSE)
rownames(st) = st[, 1]
st <- st[, -1, drop=FALSE]
s <- st
# Feature_Table
if (f != "NA") {
    ft <- read.delim(f, check.names=FALSE, stringsAsFactors=FALSE)
    rownames(ft) <- ft[, 1]
    colnames(ft) <- gsub("symbol","SYMBOL",colnames(ft), ignore.case=TRUE)
}

# main
library(limma)
library(statmod)
library(Biobase)
## filter
if (f != "NA") nm <- nm[rownames(ft), ] # keeping only the features in Feature_Table
nm <- nm[, rownames(st)] # keeping only the samples in Sample_Table
if (length(l) > 1) {
    nm <- nm[rowSums(nm > l[1]) >= l[2], ]
    ft <- ft[rownames(nm),]
}
## ExpressionSet
if (f != "NA") {
    eset <- ExpressionSet(assayData=nm, phenoData=new("AnnotatedDataFrame", data=st), featureData=new("AnnotatedDataFrame", data=ft))
} else {
    eset <- ExpressionSet(assayData=nm, phenoData=new("AnnotatedDataFrame", data=st))
}
## design
design <- model.matrix(eval(parse(text=m))) # design
## fit
if (d != "NA") {
    cl <- duplicateCorrelation(eset, design, block=st[, d]) # Cor_List
    fit <- lmFit(eset, design, block=st[,d], correlation=cl$consensus) # fit
} else {
    fit <- lmFit(eset, design)
}
fit <- eBayes(fit, robust=r, trend=TRUE)

## Differently expressed genes
fit_de <- topTable(fit, coef=e, sort="p", n=Inf, adjust="BH")

# output
## fit
save(fit, file=paste0(o, ".fit.RData"))
write.table(fit_de, paste0(o, '.statistics.tsv'), sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)

## Splicing
if (ds != "NA") {
    fit_as <- diffSplice(fit, geneid=ds[1], exonid=ds[2], robust=r)
    fit_as <- topSplice(fit_as, coef=e, test=ss, number=Inf)
    write.table(fit_as, paste0(o, '.as.tsv'), sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)
}

## Plot
if (p) {
    pdf(paste0(o, '.plots.pdf'))
    par(mfrow=c(2,2))
    plotMDS(eset)
    plotSA(fit, main="Final model: Mean-variance trend", xlab="logExp")
    volcanoplot(fit, coef=e, highlight=10, names=fit$genes$SYMBOL, main="Volcano of Fit", xlab="logFC", ylab="logOdds")
    plotMD(fit, xlab="Average logExp", ylab="logFC"); abline(h=c(-1, 1), col="blue"); abline(h=0, col="darkgrey")
    dev.off()
}
### THE END ###

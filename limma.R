#!/usr/bin/env Rscript

library("optparse")
library("data.table")
library("methods")

V="Version: 1.5"
D="Depends: R (>= 3.1.0), optparse, data.table, edgeR, limma, statmod"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-o", "--output"), type="character", default=".",
     help="Prefix for all output files [%default].", metavar="character"),
    make_option(c("-f", "--feature_table"), type="character", default="NA",
     help="Feature table. (row: feature; column: annotation, separator: '\\t') [%default].", metavar="character"),
    make_option(c("-l", "--exp_limit"), type="character", default="1,3",
     help="Expression limit. '1,3' represents 'keep features that \n\t\thave more than 1 cpm in at least 3 samples.' [%default].", metavar="character"),
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

opt_parser <- OptionParser(usage="usage: %prog [options] <count.tsv> <sample.tsv>\n\t<count.tsv> is a TSV file for feature count (row: feature, column: sample, separator: '\\t', stdin: '-').\n\t<sample.tsv> is a TSV file for sample annotation (row: sample, column: annotation, separator: '\\t').", option_list=option_list, description = paste(V, D, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

o <- opt$options$output
f <- opt$options$feature_table
l <- opt$options$exp_limit
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

n <- ifelse(opt$args[1]=='-', 'file:///dev/stdin', opt$args[1]) # couNts
s <- opt$args[2] # Sample_Table
l <- as.numeric(unlist(strsplit(l,",")))
ds <- unlist(strsplit(ds,","))

if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)
if(grepl("^[[:digit:]]+$", ra)) ra=as.numeric(ra)
# read
nm <- fread(n, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE) # couNt_Matrix
rownames(nm) <- nm[, ra]
nm <- nm[, -ra]
nm <- as.matrix(nm)

st <- read.delim(s, check.names=FALSE) # Sample_Table
rownames(st) = st[, 1]
st <- st[, -1, drop=FALSE]
s <- st

if (f != "NA") {
    ft <- read.delim(f, check.names=FALSE, stringsAsFactors=FALSE) # Feature_Table
    rownames(ft) <- ft[, 1]
}
colnames(ft) <- gsub("symbol","SYMBOL",colnames(ft), ignore.case=TRUE)

# main
library(edgeR)
library(limma)
library(statmod)

## trim
nm <- nm[rownames(ft), ] # keeping only the features in Feature_Table
nm <- nm[, rownames(st)] # keeping only the samples in Sample_Table

ne <- DGEList(counts=nm)
if ("lib_size" %in% colnames(st)) ne$samples$lib.size=st[, "lib_size"]
if (f != "NA") ne$genes=ft

ne_0 <- ne # keep raw data
el <- rowSums(cpm(ne) > l[1]) >= l[2] # Genes must be expressed in at least three samples (Expr_Log)
if ("lib_size" %in% colnames(st)) {ne <- ne[el, , keep.lib.sizes=TRUE]} else ne <- ne[el, , keep.lib.sizes=FALSE] # keep expressed genes
ne_1 <- ne # keep filtered data
ne <- calcNormFactors(ne) # normalizationation
design <- model.matrix(eval(parse(text=m))) # design
ve <- voom(ne, design) # couNt_Voom
## fit
if (d != "NA") {
    cl <- duplicateCorrelation(ve, design, block=st[, d]) # Cor_List
    fit <- lmFit(ve, design, block=st[,d], correlation=cl$consensus) # fit
} else {
    fit <- lmFit(ve, design)
}
fit <- eBayes(fit, robust=r)

## Differently expressed genes
fit_de <- topTable(fit, coef=e, sort="p", n=Inf, adjust="BH")

# output
## DEGs
cM <- cpm(ne_0, log=TRUE, prior.count=0.5)
ct <- base::merge(ft, data.frame(cm_id=rownames(cM), cM), by=1, all.y=TRUE)
write.table(ct, paste0(o, '.logCPM.tsv'), sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)

vt <- base::merge(ft, data.frame(ve_id=rownames(ve$E), ve$E), by=1, all.y=TRUE)
write.table(vt, paste0(o, '.logCPM.voom.tsv'), sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)

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
    ### check the filter
    cp <- rainbow(ncol(nm)) # Color_Pattern
    par(mfrow=c(2,2))
    
    lcpm <- cpm(ne_0, log=TRUE, prior.count=0.5)
    plot(density(lcpm[,1]), col=cp[1], ylim=c(0,1), main="Raw Data", xlab="log2(CPM + 0.5)", ylab="Density")
    abline(v=0, lty=3)
    for (i in 2:ncol(nm)){den <- density(lcpm[,i]); lines(den$x, den$y, col=cp[i])}
    legend("topright", colnames(ne), text.col=cp, bty="n")
    
    lcpm <- cpm(ne_1, log=TRUE, prior.count=0.5)
    plot(density(lcpm[,1]), col=cp[1], ylim=c(0,1), main="Filtered Data", xlab="log2(CPM + 0.5)", ylab="Density")
    abline(v=0, lty=3)
    for (i in 2:ncol(nm)){den <- density(lcpm[,i]); lines(den$x, den$y, col=cp[i])}
    legend("topright", colnames(ne), text.col=cp, bty="n")
    
    voom(ne, design, plot=TRUE)
    plotSA(fit, main="Final model: Mean-variance trend", xlab="log2(CPM + 0.5)")
    ### check the normalization
    boxplot(ve$E, ylab="log2(CPM + 0.5)", las=2)
    ### check the replicates
    plotMDS(ve)
    ### check the DEGs
    volcanoplot(fit, coef=e, highlight=10, names=fit$genes$SYMBOL, main="Volcano of Fit", xlab="logFC", ylab="logOdds")
    plotMD(fit, xlab="Average logExp", ylab="logFC"); abline(h=c(-1, 1), col="blue"); abline(h=0,col="darkgrey")
    dev.off()
}
### THE END ###

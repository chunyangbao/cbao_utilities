#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.0"
D="Depends: R (>= 3.4.0), optparse, data.table"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-o", "--output"), type="character", default=".",
     help="Output directory [%default].", metavar="character"),
    make_option(c("--bins"), type="integer", default="10",
     help="Bins in the histogram [%default].", metavar="integer"),
    make_option(c("--xcut"), type="character", default="NA",
     help="Cutoff values [%default]. \n\t\te.g. '-2,-1,0,1,2'", metavar="character"),
    make_option(c("--xlim"), type="character", default="NA",
     help="Cutoff values [%default]. \n\t\te.g. '-2,2'", metavar="character"),
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

opt_parser <- OptionParser(usage="usage: %prog [options] <input.tsv> \n\t<input.tsv> is a TSV file (1st column: sample, 2nd column: feature, separator: '\\t', stdin: '-').", option_list=option_list, description = paste(V, D, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

o <- opt$options$output
bins <- opt$options$bins
xcut <- opt$options$xcut
xlim <- opt$options$xlim

cs <- opt$options$col_sep
rn <- opt$options$row_num
cn <- as.logical(opt$options$col_name)
na <- opt$options$na_str
rs <- opt$options$row_skip
ra <- opt$options$row_name

if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)
if(grepl("^[[:digit:]]+$", ra)) ra=as.numeric(ra)

if(xcut!="NA") xcut <- as.numeric(strsplit(xcut, ",")[[1]])

d <- ifelse(opt$args[1]=='-', 'file:///dev/stdin', opt$args[1]) # Data_path

# read
dta <- fread(d, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE, stringsAsFactors=TRUE) # D_Table

if(xlim=="NA") {
    xlim[1] <- floor(min(abs(dta[, 2])))
    xlim[2] <- ceiling(max(abs(dta[, 2])))
    xlim <- as.numeric(xlim)
} else {
    xlim <- as.numeric(strsplit(xlim, ",")[[1]])
    dta <- dta[dta[,2]>=xlim[1]&dta[,2]<=xlim[2],]
}

# main
pdf(o)
for(i in 1:length(levels(dta[, 1]))){
    si <- levels(dta[, 1])[i] # Sample_I
    di <- dta[dta[, 1] == si, 2] # Data_I
    hist(di, breaks=seq(xlim[1], xlim[2], (xlim[2]-xlim[1])/bins), prob=TRUE,
        main=si, xlim=c(xlim[1], xlim[2]), xlab=colnames(dta)[2])
    lines(density(di), col='blue', lwd=2)
    if(is.numeric(xcut)) abline(v = xcut, col = "red", lwd=2, lty=2)
}
dev.off()

### THE END ###

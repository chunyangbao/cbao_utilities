#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.0"
D="Depends: R (>= 3.4.0), optparse, data.table, edgeR, limma, statmod"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
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

opt_parser <- OptionParser(usage="usage: %prog [options] <count.tsv> <sample.tsv>\n\t<seg.tsv> is a TSV file for segmentations (row: seg; column: feature, chromosome, start, end and samples, separator: '\\t', stdin: '-').", option_list=option_list, description = paste(V, D, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

cs <- opt$options$col_sep
rn <- opt$options$row_num
cn <- as.logical(opt$options$col_name)
na <- opt$options$na_str
rs <- opt$options$row_skip
ra <- opt$options$row_name

s <- ifelse(opt$args[1]=='-', 'file:///dev/stdin', opt$args[1]) # Seg

if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)
if(grepl("^[[:digit:]]+$", ra)) ra=as.numeric(ra)
# read
st <- fread(s, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE) # couNt_Matrix
# main
o <- data.frame(Sample=c(), Chromosome=c(), Start=c(), End=c(), Num_Probes=c(), Segment_Mean=c())

for (i in 5:ncol(st)) {
	d <- cbind(st[,2:4], st[,i])

	rle(paste(d[,1], d[,4], sep=":")) -> rleD

	endI <- cumsum(rleD$lengths)
	posI <- c(1, endI[-length(endI)] + 1)

	chr <- d[posI,1]
	pos <- d[posI,2]
	end <- d[endI,3]
	segVal <- d[endI,4]
	bins <- rleD$lengths

	options(scipen=100) #  force R to use regular numbers

	out <- data.frame(Sample=colnames(st)[i], Chromosome=chr, Start=pos, End=end, Num_Probes=bins, Segment_Mean=segVal)

	o <- rbind(o, out)

	}
# output
write.table(o, quote=F, sep="\t", append=FALSE, col.names=TRUE, row.names=FALSE)
### THE END ###

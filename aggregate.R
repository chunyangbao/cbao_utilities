#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.0"
D="Depends: R (>= 3.1.0), optparse, data.table"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-I", "--input"), type="character", default="-",
     help="Input file [%default] \n\t\tUse '-' if passing 'I' with a UNIX pipe.", metavar="character"),
    make_option(c("-f", "--formula"), type="character", default="y~x", 
     help="a formula, such as y ~ x or cbind(y1, y2) ~ x1 + x2. [%default]", metavar="character"),
    make_option(c("-F", "--fun"), type="character", default="sum", 
     help="a function to compute the summary statistics. [%default]", metavar="character"),
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
     help="Row can be taken as the first data row [%default] \n\t\tskip='string' searches for 'string' in the file and starts on that row.", metavar="character")
)

opt_parser = OptionParser(usage="usage: %prog [options] -I <file> -f <formula> -F <function>\n\tBy default <file> is -' (stdin)", option_list=option_list, description = paste(V, D, sep="\n"))
opt = parse_args(opt_parser, args=a, positional_arguments=TRUE)

I = opt$options$input
f = opt$options$formula
F = opt$options$fun

cs = opt$options$col_sep
rn = opt$options$row_num
cn = as.logical(opt$options$col_name)
na = opt$options$na_str
rs = opt$options$row_skip

if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)

if(I=="-") I = "file:///dev/stdin"

# main
it = fread(I, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE, check.names=FALSE) # Input_Table
out_ta=aggregate(formula(f), data = it, F)
write.table(out_ta,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

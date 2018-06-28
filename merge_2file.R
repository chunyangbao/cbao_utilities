#!/usr/bin/env Rscript

library("optparse")
library("data.table")

v="Version: 1.0"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-a", "--file_a"), type="character", default="-",
     help="File 'a' [%default] \n\t\tUse '-' if passing 'a' with a UNIX pipe.", metavar="character"),
    make_option(c("-b", "--file_b"), type="character", default="-",
     help="File 'b' [%default] \n\t\tUse '-' if passing 'b' with a UNIX pipe.", metavar="character"),
    make_option(c("-x", "--key_a"), type="character", default="1",
     help="Key of File 'a' [%default] \n\t\tComma separated string is acceptable (i.e. '1,3,2' or 'A,B,C').", metavar="character"),
    make_option(c("-y", "--key_b"), type="character", default="1",
     help="Key of File 'b' [%default] \n\t\tComma separated string is acceptable (i.e. '1,3,2' or 'A,B,C').", metavar="character"),
    make_option(c("-J", "--row_intersect"), type="logical", default="FALSE", action = "store_true",
     help="Only output the rows with the shared key (intersection) [%default]", metavar="logical"),
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

opt_parser = OptionParser(usage="usage: %prog [options] -a <file> -b <file>\n\tBy default <file> is -' (stdin)", option_list=option_list, description = v)
opt = parse_args(opt_parser, args=a, positional_arguments=TRUE)

a = opt$options$file_a
b = opt$options$file_b
x = unlist(strsplit(opt$options$key_a, ",", fixed=TRUE))
y = unlist(strsplit(opt$options$key_b, ",", fixed=TRUE))
J = opt$options$row_intersect
cs = opt$options$col_sep
rn = opt$options$row_num
cn = as.logical(opt$options$col_name)
na = opt$options$na_str
rs = opt$options$row_skip

if(sum(grepl("^[[:digit:]]+$", x))==length(x)) x=as.numeric(x)
if(sum(grepl("^[[:digit:]]+$", y))==length(y)) y=as.numeric(y)
if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)

if(a=="-") a = "file:///dev/stdin"
if(b=="-") b = "file:///dev/stdin"

fa = fread(a, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE)
fb = fread(b, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE)
o = merge(fa, fb, by.x=x, by.y=y, all=!J, sort=FALSE)
write.table(o, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

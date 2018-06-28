#!/usr/bin/env Rscript

library("optparse")
library("data.table")

v="Version: 2.2"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-p", "--file_pattern"), type="character", default=".",
     help="Pattern of filename [%default] \n\t\tA regular expression in R.", metavar="character"),
    make_option(c("-k", "--col_key"), type="character", default="1",
     help="Key column(s) [%default] \n\t\tComma separated string is acceptable (i.e. '1,3,2' or 'A,B,C').", metavar="character"),
    make_option(c("-J", "--row_intersect"), type="logical", default="FALSE", action = "store_true",
     help="Only output the rows with the shared key (intersection) [%default]", metavar="logical"),
    make_option(c("-R", "--dir_recurse"), type="logical", default="FALSE", action = "store_true",
     help="Should the listing recurse into directories? [%default]", metavar="logical"),

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
    make_option(c("--col_select"), type="character", default="NULL",
     help="Vector of column names or numbers to keep, drop the rest. [%default]", metavar="character")
)

opt_parser = OptionParser(usage="usage: %prog [options] <dir>\n\tBy default <dir> is '.' (current directory)", option_list=option_list, description = v)
opt = parse_args(opt_parser, args=a, positional_arguments=TRUE)

p = opt$options$file_pattern
k = unlist(strsplit(opt$options$col_key, ",", fixed=TRUE))
J = opt$options$row_intersect
R = opt$options$dir_recurse

cs = opt$options$col_sep
rn = opt$options$row_num
cn = ifelse(opt$options$col_name=='auto', opt$options$col_name,as.logical(opt$options$col_name))
na = opt$options$na_str
rs = opt$options$row_skip
ce = opt$options$col_select

if(sum(grepl("^[[:digit:]]+$", k))==length(k)) k=as.numeric(k)
if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)
if(ce != "NULL") ce = unlist(strsplit(ce, ",", fixed=TRUE))
if(sum(grepl("^[[:digit:]]+$", ce))==length(ce)) ce=as.numeric(ce)

d = opt$args[1]
if(is.na(d)) d = "."

setwd(d)
fs = dir(recursive=R, pattern=p)

if(!"NULL" %in% ce) o = fread(fs[1], sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, select=ce, data.table=FALSE)
if("NULL" %in% ce) o = fread(fs[1], sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE)
if(is.numeric(k)) {
    nk = k # Number of Key column(s)
    nk_o = k
    hk_o = colnames(o)[k]
} else {
    hk_o = k # Header of Key column(s)
    nc_o = 1:length(colnames(o)) # Numeber of Column(s)
    names(nc_o) = colnames(o)
    nk_o = nc_o[k]
}
colnames(o) = paste(colnames(o), sub("\\.[^\\.]*$", "", basename(fs[1])), sep=".")
colnames(o)[nk_o] = hk_o

for (f in fs[-1]) {
    if(!"NULL" %in% ce) ti = fread(f, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, select=ce, data.table=FALSE)
    if("NULL" %in% ce) ti = fread(f, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE)
    if(!is.numeric(k)){
        hk = k # Header of Key column(s)
        nc = 1:length(colnames(ti))
        names(nc) = colnames(ti)
        nk = nc[k]
    }
    colnames(ti) = paste(colnames(ti), sub("\\.[^\\.]*$", "", basename(f)), sep=".")
    colnames(ti)[nk] = hk_o
    o = merge(o, ti, by.x=nk_o, by.y=nk, all=!J, sort=FALSE)
}

write.table(o, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

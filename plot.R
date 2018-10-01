#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V <- "Version: 1.0"
D <- "Depends: R (>= 3.1.0), optparse, data.table"

a <- commandArgs(trailingOnly=TRUE)

option_list <- list(
    make_option(c("-o", "--output"), type="character", default=".",
     help="Output directory [%default] \n\t\tBy default <o> is '.' (current directory).", metavar="character"),
    make_option(c("-l", "--lib_path"), type="character", default="NA",
     help="Location of the default library [%default]", metavar="character"),
     # main
    make_option(c("-b", "--lib"), type="character", default="NA",
     help="Only set if another package is required [%default] ", metavar="character"),
    make_option(c("-f", "--fun"), type="character", default="NA",
     help="function (required). Possible value is 'function1(...) + function12(...)' [%default].", metavar="character"),
    make_option(c("-s", "--plot_size"), type="character", default="7,7",
     help="Size of the plot in inches separated by ',' [%default]", metavar="character"),
     # fread
    make_option(c("-d", "--data_table"), type="logical", default="FALSE", action = "store_true",
     help="Convert input file to data.table [%default].", metavar="logical"),
    make_option(c("--out_quote"), type="character", default="FALSE",
     help="Should output be surrounded by double quotes, 'auto'|TRUE|FALSE [%default]", metavar="character"),
    make_option(c("--str2fac"), type="logical", default="FALSE", action = "store_true",
     help="Should strings be converted to factors [%default].", metavar="logical"),
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

opt_parser <- OptionParser(usage="usage: %prog [options] <input.tsv>\n\tBy default <input.tsv> is '-' (stdin).", option_list=option_list, description = paste(V, D, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

# options
o <- opt$options$output
l <- opt$options$lib_path
b <- opt$options$lib
f <- opt$options$fun
ps = as.numeric(unlist(strsplit(opt$options$plot_size, ",", fixed=TRUE)))

da <- opt$options$data_table
oq <- opt$options$out_quote
sf <- opt$options$str2fac
cs <- opt$options$col_sep
rn <- opt$options$row_num
cn <- opt$options$col_name
na <- opt$options$na_str
rs <- opt$options$row_skip

d <- ifelse(opt$args[1]=='-', 'file:///dev/stdin', opt$args[1])

if(oq != 'auto') oq <- as.logical(oq)
if(cn != 'auto') cn <- as.logical(cn)
if(grepl("^[[:digit:]]+$", rs)) rs <- as.numeric(rs)

# read
d <- fread(d, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, stringsAsFactors=sf, data.table=da)

# main
if(l != 'NA') .libPaths(c(l, .libPaths()))

## parameters
ps = c(w = ps[1], h = ps[2])
## library
if(b != 'NA') library(b, character.only=T)
## call
pdf(o, width = ps['w'], height = ps['h'])
    eval(parse(text=f))
dev.off()

### THE END ###

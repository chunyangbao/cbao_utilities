#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.1"
D="Depends: R (>= 3.4.0), optparse, data.table"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-o", "--output"), type="character", default=".",
     help="Output file [%default].", metavar="character"),
    make_option(c("-m", "--mapping"), type="character", default="NA",
     help="Mapping (required). Possible value is 'aes(...) + geom_bar(...)' [%default].", metavar="character"),
    make_option(c("-b", "--background"), type="character", default="NA",
     help="Background. Possible values are 'NA, blank, ...' [%default].", metavar="character"),
    make_option(c("-s", "--plot_size"), type="character", default="7,7",
     help="Size of the plot in inches separated by ',' [%default]", metavar="character"),
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

opt_parser <- OptionParser(usage="usage: %prog [options] <input.tsv> \n\t<input.tsv> is a TSV file (separator: '\\t', stdin: '-').", option_list=option_list, description = paste(V, D, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

o <- opt$options$output
m <- opt$options$mapping
b <- opt$options$background
ps = as.numeric(unlist(strsplit(opt$options$plot_size, ",", fixed=TRUE)))

cs <- opt$options$col_sep
rn <- opt$options$row_num
cn <- as.logical(opt$options$col_name)
na <- opt$options$na_str
rs <- opt$options$row_skip
ra <- opt$options$row_name

ps = c(w = ps[1], h = ps[2])

if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)
if(grepl("^[[:digit:]]+$", ra)) ra=as.numeric(ra)

d <- ifelse(opt$args[1]=='-', 'file:///dev/stdin', opt$args[1]) # Data_path

# read
d <- fread(d, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE, stringsAsFactors=TRUE) # D_Table

# main
library(ggplot2)

b <- switch (b,
        "NA" = NULL,
        "blank" = " + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))"
)

g <- paste0(m, b)

pdf(o, width = ps['w'], height = ps['h'])
    eval(parse(text=g))
dev.off()

### THE END ###

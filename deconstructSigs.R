#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V <- "Version: 2.0"
D <- "Depends: R (>= 3.1.0), optparse, data.table, deconstructSigs"

a <- commandArgs(trailingOnly=TRUE)

option_list <- list(
    make_option(c("-o", "--output"), type="character", default=".",
     help="Output directory [%default] \n\t\tBy default <o> is '.' (current directory).", metavar="character"),
    make_option(c("-l", "--lib_path"), type="character", default="NA",
     help="Location of the default library [%default]", metavar="character"),

     # deconstructSigs
    make_option(c("-b", "--bsgenome"), type="character", default="BSgenome.Hsapiens.UCSC.hg19",
     help="Only set if another genome build is required. Must be a name of BSgenome \n\t\tobject. [%default] ", metavar="character"),
    make_option(c("-r", "--signatures_ref"), type="character", default="signatures.cosmic",
     help="Either a data frame or location of signature text file, where rows are \n\t\tsignatures, columns are trinucleotide contexts [%default] ", metavar="character"),
    make_option(c("-s", "--sample_id"), type="character", default="NA",
     help="Name of sample. Only one item is acceptable. If not given, only the \n\t\tfirst sample is used [%default] ", metavar="character"),
    make_option(c("-m", "--tri_counts_method"), type="character", default="default",
     help="The method of normalization. \n\t\tPossible values: 'exome', 'genome', 'exome2genome', 'genome2exome' or location of scaling factor data.frame [%default] ", metavar="character"),
    make_option(c("-a", "--associated"), type="character", default="NA",
     help="Vector of associated signatures. If given, will narrow the signatures \n\t\ttested to only the ones listed. [%default] ", metavar="character"),
    make_option(c("-i", "--signatures_limit"), type="character", default="NA",
     help="Number of signatures to limit the search to [%default] ", metavar="character"),
    make_option(c("-u", "--signature_cutoff"), type="double", default="0.06",
     help="Discard any signature contributions with a weight less than this amount [%default] ", metavar="numeric"),
    make_option(c("-p", "--plots"), type="logical", default="FALSE", action = "store_true",
     help="Should the plots be generated [%default].", metavar="logical"),
     # fread
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

opt_parser <- OptionParser(usage="usage: %prog [options] <input.maf>\n\tBy default <input.maf> is '-' (stdin). A one-sample MAF file is recommended.", option_list=option_list, description = paste(V, D, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

# options
o <- opt$options$output
lp <- opt$options$lib_path
bg <- opt$options$bsgenome
sr <- opt$options$signatures_ref
si <- opt$options$sample_id
tm <- opt$options$tri_counts_method
as <- opt$options$associated
sl <- opt$options$signatures_limit
sc <- opt$options$signature_cutoff
p <- opt$options$plots

oq <- opt$options$out_quote
sf <- opt$options$str2fac
cs <- opt$options$col_sep
rn <- opt$options$row_num
cn <- opt$options$col_name
na <- opt$options$na_str
rs <- opt$options$row_skip

if(oq != 'auto') oq <- as.logical(oq)
if(cn != 'auto') cn <- as.logical(cn)
if(grepl("^[[:digit:]]+$", rs)) rs <- as.numeric(rs)

# input
## R version
rv0 <- as.numeric(R.version$major)
rv1 <- as.numeric(R.version$minor)
## stdin
if ((rv0<3)|((rv0==3)&(rv1<5))) {
    d <- ifelse(opt$args[1]=='-', 'file:///dev/stdin', opt$args[1])
} else {
    d <- ifelse(opt$args[1]=='-', 'cat /dev/stdin', opt$args[1])
}

# read
input_df <- fread(d, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, stringsAsFactors=sf, data.table=FALSE)

# main
if(lp != 'NA') .libPaths(c(lp, .libPaths()))
library(deconstructSigs)

## parameters
si <- ifelse(si == 'NA', input_df[1,'Tumor_Sample_Barcode'], si)
if(bg == 'BSgenome.Hsapiens.UCSC.hg19') {bg <- NULL} else {library(bg,character.only=T); bg <- get(bg)}
if(as == 'NA') {as <- NULL} else {as <- unlist(strsplit(as, ",", fixed=TRUE))}
if(sl == 'NA') {sl <- NA} else {sl <- as.numeric(sl)}
if(!sr %in% c('signatures.nature2013', 'signatures.cosmic')) {sr <- utils::read.table(sr, sep="\t", header=TRUE, as.is=TRUE, check.names=FALSE, row.names=1)} else {sr <- get(sr)}
if(!tm %in% c('default', 'exome', 'genome', 'exome2genome', 'genome2exome')) tm <- utils::read.table(tm, sep="\t", header=TRUE, as.is=TRUE, check.names=FALSE, row.names=1)

## Convert to deconstructSigs input
context_df <- mut.to.sigs.input(mut.ref = input_df,
                                sample.id = "Tumor_Sample_Barcode",
                                chr = "Chromosome",
                                pos = "Start_position",
                                ref = "Reference_Allele",
                                alt = "Tumor_Seq_Allele2",
                                bsg = bg)
## call sigs
sig_li = whichSignatures(tumor.ref = context_df,
                         signatures.ref = sr,
                         sample.id = si,
                         contexts.needed = TRUE,
                         tri.counts.method = tm,
                         associated = as,
                         signatures.limit = sl,
                         signature.cutoff = sc)

sig_unknown_df <- data.frame(unknown=sig_li$unknown, row.names=si)

# output
if (!dir.exists(o)) dir.create(o)
setwd(o)

fwrite(context_df, file=paste0(si, '.context_df.tsv'), quote=oq, sep=cs, na=na, row.names=TRUE, col.names=TRUE)
fwrite(as.data.frame(sig_li$weights), file=paste0(si, '.sig_weights.tsv'), quote=oq, sep=cs, na=na, row.names=TRUE, col.names=TRUE)
fwrite(as.data.frame(sig_li$tumor), file=paste0(si, '.sig_tumor.tsv'), quote=oq, sep=cs, na=na, row.names=TRUE, col.names=TRUE)
fwrite(as.data.frame(sig_li$product), file=paste0(si, '.sig_product.tsv'), quote=oq, sep=cs, na=na, row.names=TRUE, col.names=TRUE)
fwrite(as.data.frame(sig_li$diff), file=paste0(si, '.sig_diff.tsv'), quote=oq, sep=cs, na=na, row.names=TRUE, col.names=TRUE)
fwrite(sig_unknown_df, file=paste0(si, '.sig_unknown.tsv'), quote=oq, sep=cs, na=na, row.names=TRUE, col.names=TRUE)

if (p) {
    pdf(paste0(si, '.plot_bar.pdf'))
    plotSignatures(sig_li)
    dev.off()
    pdf(paste0(si, '.plot_pie.pdf'))
    makePie(sig_li)
    dev.off()
}

### THE END ###

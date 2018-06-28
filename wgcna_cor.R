#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.0"
D="Depends: R (>= 3.1.0), optparse, data.table, WGCNA, pheatmap"

A = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-a", "--file_a"), type="character", default="NA",
     help="File 'a' [%default] \n\t\tUse '-' if passing 'a' with a UNIX pipe.", metavar="character"),
    make_option(c("-b", "--file_b"), type="character", default="NA",
     help="File 'b' [%default] \n\t\tUse '-' if passing 'b' with a UNIX pipe.", metavar="character"),
    make_option(c("-o", "--dir_out"), type="character", default=".",
     help="Output directory [%default] \n\t\tBy default <o> is '.' (current directory).", metavar="character"),
    make_option(c("-l", "--log_val"), type="character", default="NA",
     help="Log transform the values in file 'a' and/or file 'b', 'a|b|both' [%default]", metavar="character"),
    make_option(c("-T", "--transpose"), type="character", default="NA",
     help="Transpose file 'a' and/or file 'b', 'a|b|both' [%default]", metavar="character"),
    make_option(c("--use"), type="character", default="all.obs",
     help="The handling of missing data [%default] \n\t\tFor other options, see R’s 'cor' or 'WGCNA::cor'.", metavar="character"),
    make_option(c("--method"), type="character", default="pearson",
     help="Method to be used [%default] \n\t\t'kendall|spearman'", metavar="character"),
    make_option(c("--quick"), type="double", default="0",
     help="Precision of handling of missing data, '0-1' [%default] \n\t\tSee also 'WGCNA::cor'.", metavar="character"),
    make_option(c("--cosine"), type="logical", default="FALSE", action = "store_true",
     help="Use the cosine calculation for file 'a' and 'b' [%default]", metavar="logical"),
     make_option(c("-t", "--nThreads"), type="integer", default="0",
     help="Number of parallel threads [%default] \n\t\t'0' means 'dynamic'.", metavar="integer"),
    make_option(c("--plot_wh"), type="character", default="NA",
     help="Width and Height of the plot in inches separated by ',', such as'7,7' [%default] \n\t\t'NA' means no plot", metavar="character"),
    make_option(c("--scale"), type="character", default="none",
     help="Values should be centered and scaled, 'row|column' [%default]", metavar="character"),
    make_option(c("--cluster_rc"), type="character", default="both",
     help="Rows and columns should be clustered, 'row|column|both' [%default]", metavar="character"),
    make_option(c("--color_lh"), type="character", default="yellow,red",
     help="Colors for the lowest or highest value [%default]", metavar="character"),
     # fread
    make_option(c("--col_sep"), type="character", default="\t",
     help="Separator between columns [\\t] \n\t\t'auto' represent for '[,\\t |;:]'. See also 'data.table::fread'.", metavar="character"),
    make_option(c("--row_num"), type="integer", default="-1",
     help="Number of rows to read [%default] \n\t\t'-1' means all.  just returns the column names.", metavar="integer"),
     make_option(c("--col_name"), type="character", default="auto",
     help="The first data line contain column names, 'auto'|TRUE|FALSE [%default] \n\t\tDefaults according to whether every non-empty field on the first data\n\t\tline is type character.", metavar="character"),
    make_option(c("--na_str"), type="character", default="NA",
     help="String(s) which are to be interpreted as NA values [%default]", metavar="character"),
    make_option(c("--row_skip"), type="character", default="0",
     help="Row can be taken as the first data row [%default] \n\t\tskip='string' searches for 'string' in the file and starts on that row.", metavar="character"),
    make_option(c("--row_name"), type="character", default="1",
     help="Column number (or name) can be taken as the row names [%default]", metavar="character")
)

opt_parser = OptionParser(usage="usage: %prog [options]\n\tBy default, the input file (a) is stdin.", option_list=option_list, description = paste(V, D, sep="\n"))
opt = parse_args(opt_parser, args=A, positional_arguments=TRUE)

a = opt$options$file_a
b = opt$options$file_b
o = opt$options$dir_out
l = opt$options$log_val
tp = opt$options$transpose
us = opt$options$use
me = opt$options$method
qu = opt$options$quick
co = opt$options$cosine
nt = opt$options$nThreads
pl = opt$options$plot_wh
sc = opt$options$scale
if(pl!='NA') {
pwh = as.numeric(unlist(strsplit(opt$options$plot_wh, ",", fixed=TRUE)))
pwh = c(w = pwh[1], h = pwh[2])
cl = opt$options$cluster_rc
lcl = switch(cl, row = c(r=TRUE, c=FALSE), column = c(r=FALSE, c=TRUE), both = c(r=TRUE, c=TRUE))
if(is.null(lcl)) lcl = c(r=FALSE, c=FALSE)
cr = as.character(unlist(strsplit(opt$options$color_lh, ",", fixed=TRUE)))
cr = c(l = cr[1], h = cr[2])
}
cs = opt$options$col_sep
rn = opt$options$row_num
cn = as.logical(opt$options$col_name)
na = opt$options$na_str
rs = opt$options$row_skip
ra = opt$options$row_name

if(a=="NA") a = opt$args[1]

ll = switch(l, a = c(a=TRUE, b=FALSE), b = c(a=FALSE, b=TRUE), both = c(a=TRUE, b=TRUE))
if(is.null(ll)) ll = c(a=FALSE, b=FALSE)
lt = switch(tp, a = c(a=TRUE, b=FALSE), b = c(a=FALSE, b=TRUE), both = c(a=TRUE, b=TRUE))
if(is.null(lt)) lt = c(a=FALSE, b=FALSE)

if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)
if(grepl("^[[:digit:]]+$", ra)) ra=as.numeric(ra)

if((a=="-")|(a=="NA")|is.na(a)) a = "file:///dev/stdin"
if(b=="-") b = "file:///dev/stdin"

fa = fread(a, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE, check.names=FALSE)
ma = as.matrix(fa[,-ra])
rownames(ma) = fa[,ra]
if(lt['a']) ma = t(ma)
if(ll['a']) ma = log2(ma + 1)

if(b!="NA") {
    fb = fread(b, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE, check.names=FALSE)
    mb = as.matrix(fb[,-ra])
    rownames(mb) = fb[,ra]
    if(lt['b']) mb = t(mb)
    if(ll['b']) mb = log2(mb + 1)
    mr = WGCNA::cor(ma, mb, use = us, method = me, quick = qu, cosine = co, nThreads = nt)
} else {
    mr = WGCNA::cor(ma, use = us, method = me, quick = qu, cosine = co, nThreads = nt)
}

if(is.null(colnames(mr))) colnames(mr) = colnames(fb)[2]
cname_a = rownames(mr)[which(abs(mr) >= 0, arr.ind = TRUE)[,1]]
cname_b = colnames(mr)[which(abs(mr) >= 0, arr.ind = TRUE)[,2]]
vr = mr[which(abs(mr) >= 0, arr.ind = TRUE)]
vp = WGCNA::corPvalueFisher(vr, nSamples = dim(ma)[1])
vq = p.adjust(vp, "BH")
tc = data.frame(colnames_a = cname_a, colnames_b = cname_b, R = vr, p_value = vp, FDR = vq)
tc = tc[order(tc[,'R'],decreasing=T),]

if (!dir.exists(o)) dir.create(o)
setwd(o)

write.table(t(c('coorelations',colnames(mr))),'coorelations.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(mr,'coorelations.tsv',sep='\t',row.names=T,col.names=F,quote=F,append=T)
write.table(tc,'statistic.tsv',sep='\t',row.names=F,col.names=T,quote=F)

if(pl!='NA') {
pdf('coorelations.pdf', width=pwh['w'], height=pwh['h'])
pheatmap::pheatmap(mr,scale=sc,cluster_rows=lcl['r'],cluster_cols=lcl['c'],color=colorRampPalette(c(cr['l'],cr['h']))(100),border_color=NA,cellwidth=10,cellheight=10)
dev.off()
}
### THE END ###

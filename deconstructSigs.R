#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.0"
D="Depends: R (>= 3.1.0), optparse, data.table, foreach, doParallel, deconstructSigs"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-o", "--dir_out"), type="character", default=".",
     help="Output directory [%default] \n\t\tBy default <o> is '.' (current directory).", metavar="character"),
    make_option(c("-p", "--processors"), type="integer", default="1",
     help="Number of processors for parallell computing [%default]", metavar="integer"),
    # deconstructSigs
    make_option(c("-s", "--sig_ref"), type="character", default="signatures.cosmic",
     help="Either a data frame or location of signature text file, where rows are \n\t\tsignatures, columns are trinucleotide contexts [%default] ", metavar="character"),
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
    make_option(c("--row_name"), type="character", default="1",
     help="Column number (or name) can be taken as the row names [%default]", metavar="character")
)

opt_parser = OptionParser(usage="usage: %prog [options] <motifs2samples.tsv>\n\tBy default <motifs2samples.tsv> is '-' (stdin). This matrix should contain \n\tmutation counts along 96 tri-nucleotide mutation contexts (rows) \n\tacross samples (columns). Rownames of the lego matrix should be \n\t4-letters. ex: T[C>G]A (C to G mutation at 5'-TCA-3'contexts)", option_list=option_list, description = paste(V, D, sep="\n"))
opt = parse_args(opt_parser, args=a, positional_arguments=TRUE)

o = opt$options$dir_out
p = opt$options$processors
s = opt$options$sig_ref

cs = opt$options$col_sep
rn = opt$options$row_num
cn = as.logical(opt$options$col_name)
na = opt$options$na_str
rs = opt$options$row_skip
ra = opt$options$row_name

x = opt$args[1]

if(is.na(x)|x=='-') x = 'file:///dev/stdin'

if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)
if(grepl("^[[:digit:]]+$", ra)) ra=as.numeric(ra)
# read
xm = fread(x, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE, check.names=FALSE)
rownames(xm) = xm[,ra]
xm = xm[,-ra]
# main
library(deconstructSigs)
library(foreach)
library(doParallel)

if(file.exists(s)) {
    s = read.delim(s, row.names=1, check.names=FALSE)
} else {
    s = get(s)
}

cl = makeCluster(p)
registerDoParallel(cl)

comb <- function(x, ...) {
  lapply(seq_along(x), function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

l = foreach (n=rownames(xm), .combine='comb', .multicombine=TRUE, 
              .packages='deconstructSigs', .init=list(list(), list(), list(), list(), list())) %dopar% {
              
    ln = whichSignatures(tumor.ref = xm, 
                        signatures.ref = s, 
                        sample.id = n, 
                        contexts.needed = TRUE,
                        tri.counts.method = 'default')

    ln$weights = cbind(ln$weights, data.frame(unknown=ln$unknown))
    ln$weights = as.matrix(ln$weights)
    lnc = round(ln$weights*sum(xm[n,])) # counts
    
    list(
    cbind(data.frame(signatures=colnames(ln$weights)), t(ln$weights)),
    cbind(data.frame(motifs=colnames(ln$tumor)), t(ln$tumor)),
    cbind(data.frame(motifs=colnames(ln$product)), t(ln$product)),
    cbind(data.frame(motifs=colnames(ln$diff)), t(ln$diff)),
    cbind(data.frame(signatures=colnames(ln$weights)), t(lnc))
    )
}

stopCluster(cl)

lo = list() # List_Output
lo[[1]] = data.frame(signatures=c(rownames(s), 'unknown'))
lo[[2]] = data.frame(motifs=colnames(s))
lo[[3]] = data.frame(motifs=colnames(s))
lo[[4]] = data.frame(motifs=colnames(s))
lo[[5]] = data.frame(signatures=c(rownames(s), 'unknown'))
for (i in 1:length(l)) {
    for (j in 1:length(l[[i]])) {
        lo[[i]] = merge(lo[[i]], l[[i]][[j]], all=TRUE)
    }
}

if (!dir.exists(o)) dir.create(o)
setwd(o)

write.table(t(colnames(lo[[1]])),'weights.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(lo[[1]],'weights.tsv',sep='\t',row.names=F,col.names=F,quote=F,append=T)
write.table(t(colnames(lo[[2]])),'tumor.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(lo[[2]],'tumor.tsv',sep='\t',row.names=F,col.names=F,quote=F,append=T)
write.table(t(colnames(lo[[3]])),'product.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(lo[[3]],'product.tsv',sep='\t',row.names=F,col.names=F,quote=F,append=T)
write.table(t(colnames(lo[[4]])),'diff.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(lo[[4]],'diff.tsv',sep='\t',row.names=F,col.names=F,quote=F,append=T)
write.table(t(colnames(lo[[5]])),'counts.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(lo[[5]],'counts.tsv',sep='\t',row.names=F,col.names=F,quote=F,append=T)

### THE END ###

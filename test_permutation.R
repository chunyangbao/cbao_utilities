#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.1"
D="Depends: R (>= 3.1.0), optparse, data.table"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-x", "--file_x"), type="character", default="-",
     help="Matrix of variable (row:feature; col: sample) [%default] \n\t\tUse '-' if passing 'x' with a UNIX pipe.", metavar="character"),
    make_option(c("-y", "--file_y"), type="character", default="-",
     help="Matrix of group (row:feature; col: sample) [%default] \n\t\tUse '-' if passing 'y' with a UNIX pipe.", metavar="character"),
    make_option(c("-p", "--processors"), type="integer", default="1",
     help="Number of processors for parallell computing [%default]", metavar="integer"),
    make_option(c("--model"), type="character", default="wilcox.test(v~g,vg)",
     help="A string deciding which statistical model will be used [%default] \n\t\tOptions include all functions in 'stats'.", metavar="character"),
    make_option(c("--returns"), type="character", default="c(unlist(r)['p.value'], effectSize(vg[vg[,'g']==1, 'v'],vg[vg[,'g']==0, 'v']))",
     help="A string deciding which statistical value will be returned [%default]", metavar="character"),
    make_option(c("--minSample"), type="integer", default="10",
     help="Min number of sample [%default]", metavar="integer"),
    make_option(c("--nPermute"), type="integer", default="0",
     help="Number of permutation [%default]", metavar="integer"),

    # fread
    make_option(c("--col_sep"), type="character", default="\t",
     help="Separator between columns [\\t] \n\t\t'auto' represent for '[,\\t |;:]'. See also 'data.table::fread'.", metavar="character"),
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

opt_parser = OptionParser(usage="usage: %prog [options] -x <file> -y <file> -z <file> -o <dir> -p 2 --cv --iter_y \n\tBy default <file> is -' (stdin)", option_list=option_list, description = paste(V, D, sep="\n"))
opt = parse_args(opt_parser, args=a, positional_arguments=TRUE)

x = opt$options$file_x
y = opt$options$file_y
p = opt$options$processors
m = opt$options$model
re = opt$options$returns
np = opt$options$nPermute
minS = opt$options$minSample

cs = opt$options$col_sep
rn = opt$options$row_num
cn = as.logical(opt$options$col_name)
na = opt$options$na_str
rs = opt$options$row_skip
ra = opt$options$row_name

if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)
if(grepl('^[[:digit:]]+$', ra)) ra=as.numeric(ra)

if(x=="-") x = "file:///dev/stdin"
if(y=="-") y = "file:///dev/stdin"

# read
vt = fread(x, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE, check.names=FALSE)
rownames(vt) = vt[,ra]
vt = as.matrix(vt[,-ra]) # Variable_Table

gt = fread(y, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE, check.names=FALSE)
rownames(gt) = gt[,ra]
gt = as.matrix(gt[,-ra]) # Group_Table

# main
library(foreach)
library(doParallel)

effectSize <- function(x, y) {
  (mean(x)-mean(y))/sqrt((((length(x)-1)*sd(x)^2)+((length(y)-1)*sd(y)^2))/(length(x)+length(y)-2))
}

stat.test <- function(v, g) {
    v=v[!is.na(v)]
    g=g[!is.na(g)]
    vg=intersect(names(v),names(g))
    v=v[vg]
    g=g[vg]
    vg=data.frame(g=g,v=v)
    vg=vg[order(rownames(vg))&order(vg[,'g']),]
    if((min(table(vg[,'g'])) < minS)|(length(table(vg[,'g']))<=1)) {
        c(NA,NA)
    } else {
        r=eval(parse(text=m))
        eval(parse(text=re))
    }
}

cl = makeCluster(p)
registerDoParallel(cl)

comb <- function(x, ...) {
  lapply(seq_along(x), function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

l = foreach (n=rownames(vt), .combine='comb', .multicombine=TRUE, .export='effectSize',
              .init=list(list(), list())) %dopar% {
              
    vtn = vt[n,]
    names(vtn) = colnames(vt)
    lr = apply(gt,1,function(gtn) stat.test(vtn,gtn))
    
    list(
    cbind(data.frame(features=colnames(lr)), matrix(lr[1,],dimnames=list(c(),n))),
    cbind(data.frame(features=colnames(lr)), matrix(lr[2,],dimnames=list(c(),n)))
    )
}
stopCluster(cl)

lo = list() # List_Output
lo[[1]] = data.frame(features=rownames(gt)) # p.value
lo[[2]] = data.frame(features=rownames(gt)) # foldchange
for (i in 1:length(l)) {
    for (j in 1:length(l[[i]])) {
        lo[[i]] = merge(lo[[i]], l[[i]][[j]], all=TRUE)
    }
    rownames(lo[[i]]) = lo[[i]][,1]
    lo[[i]] = lo[[i]][,-1]
    lo[[i]] = as.matrix(lo[[i]])
}

i = which(lo[[1]] >= 0, arr.ind = TRUE)
rvo = rownames(lo[[1]])[i[,1]]
cvo = colnames(lo[[1]])[i[,2]]
pvo = lo[[1]][i]
qvo = p.adjust(pvo, "BH")
fvo = lo[[2]][i]
to = data.frame(x = rvo, y = cvo, foldchange = fvo, p_value = pvo, FDR = qvo)
to = to[order(to[,'p_value'],decreasing=F),]

# y1 <- simulate(nullmodel(y, "curveball"), nsim=2, seed=1, burnin=10, thin=5)

write.table(to,sep='\t',row.names=F,col.names=T,quote=F)

### THE END ###

#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.0"
D="Depends: R (>= 3.1.0), optparse, data.table, foreach, doParallel"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-p", "--processors"), type="integer", default="1",
     help="Number of processors for parallell computing [%default]", metavar="integer"),
    make_option(c("-o", "--dir_out"), type="character", default=".",
     help="Output directory [%default] \n\t\tBy default <o> is '.' (current directory).", metavar="character"),
    # waterfall
    make_option(c("-t", "--type"), type="character", default="lnIC50",
     help="Type of drug response, 'lnIC50'|'AUC'|'Amax' [%default] ", metavar="character"),
    make_option(c("--nGroup"), type="integer", default="3",
     help="Number of groups for samples based on drug response [%default] ", metavar="integer"),
    make_option(c("--interFold"), type="character", default="default",
     help="Intermediate fold [%default] \n\t\tFor 'default', 'IC50':4, 'lnIC50':4, 'AUC':1.2, 'Amax':1.2", metavar="character"),
    make_option(c("--minR"), type="double", default="0.95",
     help="Minimum of Pearson correlation coefficient to the linear fit [%default] ", metavar="double"),
    make_option(c("--plot"), type="logical", default="TRUE",
     help="Plot or not [%default]", metavar="logical"),
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

opt_parser = OptionParser(usage="usage: %prog [options] <samples2drugs.tsv>\n\tBy default <motifs2samples.tsv> is '-' (stdin). This matrix should contain \n\tdrug response (columns) across samples (rows).", option_list=option_list, description = paste(V, D, sep="\n"))
opt = parse_args(opt_parser, args=a, positional_arguments=TRUE)

p = opt$options$processors
o = opt$options$dir_out
ty = opt$options$type
g = opt$options$nGroup
nf = opt$options$interFold
rl = opt$options$minR
po = opt$options$plot

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
if(ty=='lnIC50') {
    xm = exp(xm)
    ty = 'IC50'
}
if(nf=='default') {
    switch (ty,"IC50" = {nf <- 4},"AUC" = {nf <- 1.2},"Amax" = {nf <- 1.2})
} else if (grepl("^[[:digit:]]+\\.*[[:digit:]]*$", nf)) {
    nf=as.numeric(nf)
} else {stop("'--interFold' should be either 'default' or a positive numeric value")}

# main
library(foreach)
library(doParallel)

##========================================================
##  Credits:
##  Theory by Paul Bourke http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
##  Based in part on C code by Damian Coventry Tuesday, 16 July 2002
##  Based on VBA code by Brandon Crosby 9-6-05 (2 dimensions)
##  With grateful thanks for answering our needs!
##  This is an R (http://www.r-project.org) implementation by Gregoire Thomas 7/11/08
##========================================================
distancePointLine <- function(x, y, slope, intercept) {
 ## x, y is the point to test.
 ## slope, intercept is the line to check distance.
 ##
 ## Returns distance from the line.
 ##
 ## Returns 9999 on 0 denominator conditions.
 x1 <- x-10
 x2 <- x+10
 y1 <- x1*slope+intercept
 y2 <- x2*slope+intercept
 dd <- distancePointSegment(x,y, x1,y1, x2,y2)
 return(dd)
}

##========================================================
##  Credits:
##  Theory by Paul Bourke http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
##  Based in part on C code by Damian Coventry Tuesday, 16 July 2002
##  Based on VBA code by Brandon Crosby 9-6-05 (2 dimensions)
##  With grateful thanks for answering our needs!
##  This is an R (http://www.r-project.org) implementation by Gregoire Thomas 7/11/08
##========================================================
distancePointSegment <- function(px, py, x1, y1, x2, y2) {
 ## px,py is the point to test.
 ## x1,y1,x2,y2 is the line to check distance.
 ##
 ## Returns distance from the line, or if the intersecting point on the line nearest
 ## the point tested is outside the endpoints of the line, the distance to the
 ## nearest endpoint.
 ##
 ## Returns 9999 on 0 denominator conditions.
 lineMagnitude <- function(x1, y1, x2, y2) sqrt((x2-x1)^2+(y2-y1)^2)
 ans <- NULL
 ix <- iy <- 0   # intersecting point
 lineMag <- lineMagnitude(x1, y1, x2, y2)
 if( lineMag < 0.00000001) {
   stop("Short segment")
 }
 u <- (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
 u <- u / (lineMag * lineMag)
 if((u < 0.00001) || (u > 1)) {
   ## closest point does not fall within the line segment, take the shorter distance
   ## to an endpoint
   ix <- lineMagnitude(px, py, x1, y1)
   iy <- lineMagnitude(px, py, x2, y2)
   if(ix > iy)  ans <- iy
   else ans <- ix
 } else {
   ## Intersecting point is on the line, use the formula
   ix <- x1 + u * (x2 - x1)
   iy <- y1 + u * (y2 - y1)
   ans <- lineMagnitude(px, py, ix, iy)
 }
 return(ans)
}

## Drug sensitivity calling using waterfall plots
## Method:
## 1. Sensitivity calls were made using one of IC50, ActArea or Amax
## 2. Sort log IC50s (or ActArea or Amax) of the cell lines to generate a “waterfall distribution”
## 3. Identify cutoff:
##  3.1 If the waterfall distribution is non-linear (pearson cc to the linear fit <=0.95), estimate the major inflection point of the log IC50 curve as the point on the curve with the maximal distance to a line drawn between the start and end points of the distribution.
##  3.2 If the waterfall distribution appears linear (pearson cc to the linear fit > 0.95), then use the median IC50 instead.
## 4. Cell lines within a 4-fold IC50 (or within a 1.2-fold ActArea or 20% Amax difference) difference centered around this inflection point are classified as being “intermediate”,  cell lines with lower IC50s (or ActArea/Amax values) than this range are defined as sensitive, and those with IC50s (or ActArea/Amax) higher than this range are called “insensitive”.
## 5. Require at least x sensitive and x insensitive cell lines after applying these criteria (x=5 in our case).

## Input:
##  type: lnIC50, AUC, Amax
##   ic50: IC50 values in micro molar (positive values)
##   AUC: Activity Area, that is area under the drug activity curve (positive values)
##   Amax: Activity at max concentration (positive values)
##  intermediate.fold: vector of fold changes used to define the intermediate sensitivities for ic50 (4), actarea (1.2) and amax (1.2) respectively
sensitivity.calling.waterfall <- function(x, group=3, type="IC50", intermediate.fold=4, cor.min.linear=0.95, name="Drug", plot=FALSE) {
  
  # type <- match.arg(type)
  
  if (is.null(names(x))) { names(x) <- paste("X", 1:length(x), sep=".") }
  
  xx <- x[complete.cases(x)]
  switch (type,
    "IC50" = {
      xx <- -log10(xx)
      ylabel <- "-log10(IC50)"
      ## 4 fold difference around IC50 cutoff
      interfold <- log10(intermediate.fold)
    },
    "AUC" = {
      ylabel <- "Activity area"
      ## 1.2 fold difference around Activity Area cutoff
      interfold <- intermediate.fold
    },
    "Amax" = {
      ylabel <- "Amax"
      ## 1.2 fold difference around Amax
      interfold <- intermediate.fold
    }
  )
  oo <- order(xx, decreasing=TRUE)
  ## test linearity with Perason correlation
  cc <- cor.test(-xx[oo], 1:length(oo), method="pearson")
  ## line between the two extreme sensitivity values
  dd <- cbind("y"=xx[oo][c(1, length(oo))], "x"=c(1, length(oo)))
  rr <- lm(y ~ x, data=data.frame(dd))
  ## compute distance from sensitivity values and the line between the two extreme sensitivity values
  ddi <- apply(cbind(1:length(oo), xx[oo]), 1, function(x, slope, intercept) {
    return(distancePointLine(x=x[1], y=x[2], slope=slope, intercept=intercept))
  }, slope=rr$coefficients[2], intercept=rr$coefficients[1])
  if(cc$estimate > cor.min.linear){
    ## approximately linear waterfall
    cutoff <- which.min(abs(xx[oo] - median(xx[oo])))
    cutoffn <- names(cutoff)[1]
  } else {
    ## non linear waterfall
    ## identify cutoff as the maximum distance
    cutoff <- which.max(abs(ddi))
    cutoffn <- names(ddi)[cutoff]
  }
  ## identify intermediate sensitivities
  switch (type,
    "IC50" = {
      rang <- c(xx[oo][cutoff] - interfold, xx[oo][cutoff] + interfold)
    },
    "AUC" = {
     rang <- c(xx[oo][cutoff] / interfold, xx[oo][cutoff] * interfold)
    },
    "Amax" = {
      rang <- c(xx[oo][cutoff] / interfold, xx[oo][cutoff] * interfold)
    }
  )
  if(group==2) rang <- c(xx[oo][cutoff], xx[oo][cutoff])
  ## compute calls
  calls <- rep(NA, length(xx))
  names(calls) <- names(xx)
  calls[xx < rang[1]] <- "resistant"
  calls[xx > rang[2]] <- "sensitive"
  calls[xx >= rang[1] & xx <= rang[2]] <- "intermediate"
  if(group==2) rang <- calls[xx >= rang[1] & xx <= rang[2]] <- "resistant"

  if (plot) {
    pdf(file.path(o, 'waterfall_plots', paste0(gsub("\\W","_",name),'.pdf')), width=5, height=10)
    par(mfrow=c(2, 1))
    ccols <- rainbow(4)
    mycol <- rep("grey", length(xx))
    names(mycol) <- names(xx)
    mycol[calls == "sensitive"] <- ccols[2]
    mycol[calls == "intermediate"] <- ccols[3]
    mycol[calls == "resistant"] <- ccols[4]
    mycol[cutoffn] <- ccols[1]
    mypch <- rep(16, length(xx))
    names(mypch) <- names(xx)
    mypch[cutoffn] <- 19
        
    plot(xx[oo], col=mycol[oo], pch=mypch[oo], ylab=ylabel, main=sprintf("%s\nWaterfall", name))
    points(x=cutoff, y=xx[cutoffn], pch=mypch[cutoffn], col=mycol[cutoffn])
    abline(a=rr$coefficients[1], b=rr$coefficients[2], lwd=2, col="darkgrey")
    lines(x=c(cutoff, cutoff), y=c(par("usr")[3], xx[cutoffn]), col="red")
    lines(x=c(par("usr")[1], cutoff), y=c(xx[cutoffn], xx[cutoffn]), col="red")
    legend("topright", legend=c(sprintf("resistant (n=%i)", sum(!is.na(calls) & calls == "resistant")), sprintf("intermediate (n=%i)", sum(!is.na(calls) & calls == "intermediate")), sprintf("sensitive (n=%i)", sum(!is.na(calls) & calls == "sensitive")), "cutoff", sprintf("R=%.3g", cc$estimate)), col=c(rev(ccols), NA), pch=c(16, 16, 16, 19, NA), bty="n")
        
    plot(ddi, pch=mypch[oo], col=mycol[oo], ylab="Distance", main=sprintf("%s\n%s", name, "Distance from min--max line"))
    points(x=cutoff, y=ddi[cutoffn], pch=mypch[cutoffn], col=mycol[cutoffn])
    legend("topright", legend=c("resistant", "intermediate", "sensitive", "cutoff"), col=rev(ccols), pch=c(16, 16, 16, 19), bty="n")
    dev.off()
  } 
  tt <- rep(NA, length(x))
  names(tt) <- names(x)
  tt[names(calls)] <- calls
  return(tt)  
}

### Parallel ###
cl = makeCluster(p)
registerDoParallel(cl)

comb <- function(x, ...) {
  lapply(seq_along(x), function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

### Main ###
if(!dir.exists(o)) dir.create(o)
if(!dir.exists(file.path(o, 'waterfall_plots'))) dir.create(file.path(o, 'waterfall_plots'))
l = foreach (n=colnames(xm), .combine='comb', .multicombine=TRUE, 
              .packages='deconstructSigs', .init=list(list())) %dopar% {
              
    xn = xm[,n]
    names(xn) = rownames(xm)
    if(length(xn[complete.cases(xn)])>10) {
        xno = sensitivity.calling.waterfall(xn, group=g, type=ty, intermediate.fold=nf, cor.min.linear=rl, plot=po, name=n)
        list(cbind(data.frame(samples=names(xno)), matrix(xno,dimnames=list(c(),n))))
    } else {
        xno = rep(NA, length(xn))
        names(xno) <- names(xn)
        list(cbind(data.frame(samples=names(xno)), matrix(xno,dimnames=list(c(),n))))
    }

}
stopCluster(cl)

lo = list() # List_Output
lo[[1]] = data.frame(samples=rownames(xm))
for (i in 1:length(l)) {
    for (j in 1:length(l[[i]])) {
        if(is.null(l[[i]][[j]])) next
        lo[[i]] = merge(lo[[i]], l[[i]][[j]], all=TRUE)
    }
}

write.table(t(colnames(lo[[1]])),file.path(o, 'waterfall.tsv'),sep='\t',row.names=F,col.names=F,quote=F)
write.table(lo[[1]],file.path(o, 'waterfall.tsv'),sep='\t',row.names=F,col.names=F,quote=F,append=T)

### THE END ###

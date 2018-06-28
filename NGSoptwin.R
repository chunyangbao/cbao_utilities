#!/usr/bin/env Rscript

library(optparse)
library(NGSoptwin, lib.loc="~/lib")

V="Version: 1.3"
D="Depends: R (>= 3.1.0), optparse, NGSoptwin"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-o", "--output"), type="character", default=".",
     help="Output directory [%default].", metavar="character")
)

opt_parser <- OptionParser(usage="usage: %prog [options] <test.pos> <control.pos>\n\t<*.pos> is tsv file containing the first two columns of the bed files.", option_list=option_list, description = paste(V, D, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

o <- opt$options$output

tp <- opt$args[1]
if(is.na(opt$args[2])) {
    test.pos <- read.pos(tp)
    optwin <- opt.win.onesample(test.pos, win.size=c(seq(10,100,10), seq(150,950,50), seq(1000,2000,by=100))*1000)
} else {
    cp <- opt$args[2]
    control.pos <- read.pos(cp)
    optwin <- opt.win(test.pos, control.pos, win.size=c(seq(10,100,10), seq(150,950,50), seq(1000,2000,by=100))*1000)
}
# output
if (!dir.exists(o)) dir.create(o)
setwd(o)
if(is.na(opt$args[2])) {
    ot=data.frame(test_AIC=attr(optwin, "winsize")[which.min(optwin$aic)],
                  test_CV=attr(optwin, "winsize")[which.max(optwin$cv)])
    write.table(ot, 'winsize.tsv', sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)
} else {
    ot=data.frame(test_AIC=attr(optwin, "winsize")[which.min(optwin$test.aic)],
                  test_CV=attr(optwin, "winsize")[which.max(optwin$test.cv)],
                  control_AIC=attr(optwin, "winsize")[which.min(optwin$control.aic)],
                  control_CV=attr(optwin, "winsize")[which.max(optwin$control.cv)])
    write.table(ot, 'winsize.tsv', sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)
}
## Plot
pdf('plots.pdf')
plot(optwin, min=T)# AIC
plot(optwin, "cv", min=T)# CV log-likelihood
dev.off()

### THE END ###

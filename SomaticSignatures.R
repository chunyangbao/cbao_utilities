#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.1"
D="Depends: R (>= 3.1.0), optparse, data.table, SomaticSignatures, \n\tBSgenome.Hsapiens.1000genomes.hs37d5 (or other BSgenome package), \n\tggplot2, ggdendro"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-p", "--file_phe"), type="character", default="NA",
     help="A file of phenotype data, tab-separated [%default] \n\t\tUse '-' if passing the file with a UNIX pipe.", metavar="character"),
    make_option(c("-o", "--dir_out"), type="character", default=".",
     help="Output directory [%default] \n\t\tBy default <o> is '.' (current directory).", metavar="character"),
    # SomaticSignatures
    make_option(c("-b", "--pkg_bsg"), type="character", default="BSgenome.Hsapiens.1000genomes.hs37d5",
     help="Package name of 'BSgenome' [%default]", metavar="character"),
    make_option(c("-n", "--n_sigs"), type="integer", default="5",
     help="Number of signatures (n_sigs) [%default] \n\t\t'n_sigs' should not be greater than the number of samples", metavar="integer"),
    make_option(c("--mc_Unify"), type="logical", default="TRUE", action = "store_false",
     help="In 'mutationContext', should NOT the alterations be converted to have a \n\t\tC/T base pair as a reference alleles? [%default]", metavar="logical"),
    make_option(c("--mm_norm"), type="logical", default="FALSE", action = "store_true",
     help="In 'motifMatrix', should the counts be normalized to frequency? [%default]", metavar="logical"),
    make_option(c("-A", "--ans"), type="logical", default="FALSE", action = "store_true",
     help="In 'identifySignatures', should the 'n_sigs' be infered \n\t\tby 'assessNumberSignatures'? [%default]", metavar="logical"),
    make_option(c("--ans_nrep"), type="integer", default="5",
     help="In 'assessNumberSignatures', how many runs should be used for \n\t\tassessing a value of 'n_sigs'? [%default]", metavar="integer"),
    make_option(c("--dcm"), type="character", default="nmf",
     help="Function to apply for the matrix decomposition, 'nmf|pca' [%default]", metavar="character"),
    make_option(c("-C", "--cbt"), type="logical", default="FALSE", action = "store_true",
     help="Correction for Batch Effects by ComBat? [%default] \n\t\t'-p' is required", metavar="logical"),
    make_option(c("--cbt_bat"), type="character", default="1",
     help="Column number (or name) can be taken as the batch variable [%default] \n\t\t'-p' is required", metavar="logical"),
    make_option(c("-S", "--sig"), type="logical", default="FALSE", action = "store_true",
     help="Identify signatures? [%default]", metavar="logical"),
    make_option(c("-P", "--plot"), type="logical", default="FALSE", action = "store_true",
     help="Should the plot be generated? [%default]", metavar="logical"),
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
    make_option(c("--col_select"), type="character", default="5,6,7,11,13,16,16",
     help="Names or numbers of column to keep, drop the rest [%default]\n\t\tComma separated string is acceptable (i.e. '1,3,2' or 'A,B,C'). \n\t\tThe order ('chr', 'start', 'end', 'ref', 'alt', 'sample', 'group') \n\t\tand name should be considered!", metavar="character"),
    make_option(c("--row_name"), type="character", default="1",
     help="Column number (or name) can be taken as the row names [%default]", metavar="character")
)

opt_parser = OptionParser(usage="usage: %prog [options] <*.maf>\n\tBy default <*.maf> is '-' (stdin)", option_list=option_list, description = paste(V, D, sep="\n"))
opt = parse_args(opt_parser, args=a, positional_arguments=TRUE)

p = opt$options$file_phe
o = opt$options$dir_out
b = opt$options$pkg_bsg
n = opt$options$n_sigs
mU = opt$options$mc_Unify
mn = opt$options$mm_norm
A = opt$options$ans
an = opt$options$ans_nrep
d = tolower(opt$options$dcm)
C = opt$options$cbt
cbb = opt$options$cbt_bat
S = opt$options$sig
pl = opt$options$plot

cs = opt$options$col_sep
rn = opt$options$row_num
cn = as.logical(opt$options$col_name)
na = opt$options$na_str
rs = opt$options$row_skip
ce = unlist(strsplit(opt$options$col_select, ",", fixed=TRUE))
ra = opt$options$row_name

m = opt$args[1] # a MAF file

if ((is.na(m)|m=='-')&p=='-') {
    write("SomaticSignatures.R: '-m' and '-p' can not be passed with a UNIX pipe simultaneously", stderr())
    q(save='no')
}

if (d!='nmf'&d!='pca') {
    write("SomaticSignatures.R: '--dcm' must be either 'nmf' or 'pca'", stderr())
    q(save='no')
}

if (is.na(m)|m=='-') m = 'file:///dev/stdin'
if (p=='-') p = 'file:///dev/stdin'
if(grepl("^[[:digit:]]+$", cbb)) cbb=as.numeric(cbb)

if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)
if(ce[1]=="NULL") ce=as.null(ce)
if(sum(grepl("^[[:digit:]]+$", ce))==length(ce)) ce=as.numeric(ce)
if(grepl("^[[:digit:]]+$", ra)) ra=as.numeric(ra)
# read
fm = fread(m, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, select=unique(ce), data.table=FALSE, check.names=FALSE)

if(sum(grepl("^[[:digit:]]+$", ce))==length(ce)) ce=c(1:(length(ce)-1),grep(ce[duplicated(ce)],ce)[1])

if (file.exists(p)) {
    fp = fread(p, sep="\t", header=TRUE, data.table=FALSE, check.names=FALSE)
    rownames(fp) = fp[,ra]
}
# motif matrix
library(SomaticSignatures)
library(b,character.only=T)

sca_vr=VRanges(
    seqnames=Rle(fm[,ce[1]]),
    ranges=IRanges(star=fm[,ce[2]], end=fm[,ce[3]]),
    ref=fm[,ce[4]],
    alt=fm[,ce[5]],
    sampleNames=fm[,ce[6]],
    group=fm[,ce[7]])

sca_mc=mutationContext(sca_vr, get(b), unify=mU)
sca_mm=motifMatrix(sca_mc, group='group', normalize=mn)
rownames(sca_mm)=gsub(' ', '_', rownames(sca_mm))
# correction for batch effects
if (C&file.exists(p)) {
    library(sva)
    model_null = model.matrix(~1, fp)
    sca_mm = ComBat(sca_mm, batch=fp[,cbb], mod=model_null)
} 
if (C&!file.exists(p)) {
    write('SomaticSignatures.R: A file of phenotype data is required', stderr())
    q(save='no')
}
if (!dir.exists(o)) dir.create(o)
setwd(o)
# write
write.table(t(c('motifs',colnames(sca_mm))),'motifMatrix.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(sca_mm,'motifMatrix.tsv',sep='\t',row.names=T,col.names=F,quote=F,append=T)
sca_mm = sca_mm[!rowSums(sca_mm)==0,]
# assess number signatures
if(A&(d=='nmf')) {
    n_sigs = c(2:n)
    gof = assessNumberSignatures(sca_mm, n_sigs, nReplicates = an)
} else if(A&(d=='pca')) {
    n_sigs = c(2:n)
    gof = assessNumberSignatures(sca_mm, n_sigs, get(paste0(d, 'Decomposition')))
}

if(A&(d=='nmf'|d=='pca')) {
    pdf(paste0('assessNumberSignatures.pdf'))
    print(plotNumberSignatures(gof))
    dev.off()
    q(save='no')
}
# identify signatures
if(S) {
    sigs = identifySignatures(sca_mm, n, get(paste0(d, "Decomposition")))
    # write
    write.table(t(c('samples',colnames(samples(sigs)))),'samples2sigs.tsv',sep='\t',row.names=F,col.names=F,quote=F)
    write.table(samples(sigs),'samples2sigs.tsv',sep='\t',row.names=T,col.names=F,quote=F,append=T)
    write.table(t(c('motifs',colnames(signatures(sigs)))),'motifs2sigs.tsv',sep='\t',row.names=F,col.names=F,quote=F)
    write.table(signatures(sigs),'motifs2sigs.tsv',sep='\t',row.names=T,col.names=F,quote=F,append=T)
    }
# plot
if(!pl) q(save='no')

library(ggplot2)
library(ggdendro)

pdf('plot1.pdf', paper="a4")
plotMutationSpectrum(sca_mc, 'group') + ggtitle('Mutation Spectrum')
ggdendrogram(clusterSpectrum(sca_mm, 'motif'), rotate = TRUE) + ggtitle('Cluster Spectrum')
if(S) {
    plotSignatures(sigs) + ggtitle('Somatic Signatures')
    plotObservedSpectrum(sigs) + ggtitle('Observed Spectrum')
    plotFittedSpectrum(sigs) + ggtitle('Fitted Spectrum')
    plotSignatureMap(sigs) + ggtitle('Signature Map')
    plotSampleMap(sigs) + ggtitle('Sample Map')
    plotSamples(sigs) + ggtitle('Samples')
    }
dev.off()

### THE END ###

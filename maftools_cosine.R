#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.0"
D="Depends: R (>= 3.4.0), optparse, data.table"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-o", "--output_prefix"), type="character", default="cosine",
     help="If provided writes resulting MAF file to an output file [%default].", metavar="character")
)

opt_parser <- OptionParser(usage="usage: %prog [options] <input.tsv> \n\t<input.tsv> is a TSV file (separator: '\\t', stdin: '-'). This matrix \n\tshould contain mutation counts along 96 tri-nucleotide mutation \n\tcontexts (rows) across samples (columns). Rownames of the lego matrix \n\tshould be 4-letters. ex: A[C>A]A (C to A mutation at 5'-ACA-3'contexts)", option_list=option_list, description = paste(V, D, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

o <- opt$options$output_prefix

d <- ifelse(opt$args[1]=='-', 'file:///dev/stdin', opt$args[1]) # Data_path

# main
.libPaths(c("~/lib", .libPaths()))
library("maftools")

# w
w_dt <- data.table::fread(input = d, stringsAsFactors = FALSE, data.table = FALSE)
w <- w_dt[,-1]
rownames(w) <- w_dt[, 1]
# cosmic
sigs = data.table::fread(input = system.file('extdata', 'signatures.txt', package = 'maftools'), stringsAsFactors = FALSE, data.table = FALSE)
colnames(sigs) = gsub(pattern = ' ', replacement = '_', x = colnames(sigs))
rownames(sigs) = sigs$Somatic_Mutation_Type
sigs = sigs[,-c(1:3)]
# merge
sigs = sigs[rownames(w),]
coSineMat = c()
for(i in 1:ncol(w)){
  sig = w[,i]
  coSineMat = rbind(coSineMat, apply(sigs, 2, function(x){
    crossprod(sig, x)/sqrt(crossprod(x) * crossprod(sig)) #Estimate cosine similarity against all 21 signatures
  }))
}

rownames(coSineMat) = colnames(w)

cairo_pdf(filename = paste0(o, ".cosine_signatures.heatmap.pdf"), width = 7, height = 5)
    pheatmap::pheatmap(mat = coSineMat, cluster_rows = FALSE, main = "cosine similarity")
dev.off()

write.table(t(c('W', colnames(coSineMat))), paste0(o, ".cosine_signatures.heatmap.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(coSineMat, paste0(o, ".cosine_signatures.heatmap.tsv"), sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

### THE END ###

#!/usr/bin/env Rscript

# functions
## ggplot
### generate top panel of the ouput plot (optional)
top_plot = function(sv_df, contig_df, cn_max) {
    sv1_df = sv_df[sv_df[['GPOS_X']] != sv_df[['GPOS2_X']], ]
    sv2_df = data.frame(GPOS_X = c(sv_df[['GPOS_X']], sv_df[['GPOS2_X']]), 
                        SV_type = c(as.character(sv_df[['SV_type']]), as.character(sv_df[['SV_type']])), 
                        SV_size = c(as.numeric(sv_df[['SV_size']]), as.numeric(sv_df[['SV_size']])))
    p = ggplot() + scale_color_identity() + 
      geom_segment(aes(x = GPOS_X, y = cn_max, xend = GPOS_X, yend = cn_max + 1, color = SV_type), 
        data = sv2_df, size = sv2_df[['SV_size']], alpha = 0.5, show.legend = FALSE) + 
      scale_x_continuous(limits = c(0, max(contig_df[['ends']])), labels = NULL) + 
      coord_cartesian(xlim = c(0, max(contig_df[['ends']])), ylim = c(cn_max, cn_max + 4), expand=FALSE) + 
      scale_y_continuous(limits = c(cn_max, cn_max + 4), breaks = cn_max) + 
      theme(plot.margin = unit(c(2, 2, 0, 2), 'lines'), plot.background = element_blank(), 
        axis.line = element_blank(), axis.ticks = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 12, face = 'bold'), axis.text.y = element_text(size = 10), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+ 
      labs(x = '', y = '')
    if (nrow(sv1_df) > 0) {
        p = p + geom_curve(aes(x = GPOS_X, y = cn_max + 1, xend = GPOS2_X, yend = cn_max + 1, color = SV_type), 
          data = sv1_df, size = sv1_df[['SV_size']], alpha = 0.5, curvature = -0.2, show.legend = FALSE)
    }
    return(p)
}

### generate middle panel of the ouput plot
mid_plot = function(cn_df, sv_df = NULL, contig_df, cn_max) {
    p = ggplot() + scale_color_identity() + 
      geom_point(aes(x = GPOS_X, y = CN, color = CN_type), 
        data = cn_df[cn_df[['CN']] >= 2, ], alpha = 0.5, size = 0.2, show.legend = FALSE) + 
      geom_segment(aes(x = ends, y = 2, xend = ends, yend = cn_max), 
        data = contig_df, color = 'gray80', show.legend = FALSE) + 
      scale_x_continuous(limits = c(0, max(contig_df[['ends']])), labels = NULL) + 
      coord_cartesian(xlim = c(0, max(contig_df[['ends']])), ylim = c(2, cn_max), expand=FALSE) + 
      scale_y_continuous(trans = 'log2', limits = c(2, cn_max), breaks = 2^seq(2, log2(cn_max))) + 
      theme(plot.margin = unit(c(-0.19, 2, 0, 2), 'lines'), plot.background = element_blank(), 
        axis.line = element_line(colour = 'black'), axis.line.x = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title.y = element_text(size = 12, face = 'bold'), axis.text.y = element_text(size = 10), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+ 
      labs(x = '', y = 'Allelic CN')
    if (!is.null(sv_df)) {
        sv2_df = data.frame(GPOS_X = c(sv_df[['GPOS_X']], sv_df[['GPOS2_X']]), 
                            SV_type = c(as.character(sv_df[['SV_type']]), as.character(sv_df[['SV_type']])), 
                            SV_size = c(as.numeric(sv_df[['SV_size']]), as.numeric(sv_df[['SV_size']])))
        p = p + geom_segment(aes(x = GPOS_X, y = 2, xend = GPOS_X, yend = cn_max, color = SV_type), 
          data = sv2_df, size = sv2_df[['SV_size']], alpha = 0.5, show.legend = FALSE)
    }
    return(p)
}

### generate bottom panel of the ouput plot
bot_plot = function(cn_df, sv_df = NULL, contig_df, cn_max) {
    p = ggplot() + scale_color_identity() + 
      geom_point(aes(x = GPOS_X, y = CN, color = CN_type), 
        data = cn_df[cn_df[['CN']] <= 2, ], alpha = 0.5, size = 0.2, show.legend = FALSE) + 
      geom_segment(aes(x = ends, y = 0, xend = ends, yend = 2), 
        data = contig_df, color = 'gray80', show.legend = FALSE) + 
      scale_x_continuous(limits = c(0, max(contig_df[['ends']])), 
        breaks = contig_df[['centers']], labels = contig_df[['names']]) + 
      coord_cartesian(xlim = c(0, max(contig_df[['ends']])), ylim = c(0, 2), expand = FALSE) +
      scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2)) + 
      theme(plot.margin = unit(c(-0.19, 2, 0, 2), 'lines'), plot.background = element_blank(), 
        axis.line = element_line(colour = 'black'), 
        axis.title.x = element_text(size = 12, face = 'bold'), axis.text.x = element_text(size = 10), axis.ticks.x = element_blank(), 
        axis.title.y = element_text(size = 12, , face = 'bold'), axis.text.y = element_text(size = 10), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
      labs(x = 'Chromosome', y = '')
    if (!is.null(sv_df)) {
        sv2_df = data.frame(GPOS_X = c(sv_df[['GPOS_X']], sv_df[['GPOS2_X']]), 
                            SV_type = c(as.character(sv_df[['SV_type']]), as.character(sv_df[['SV_type']])), 
                            SV_size = c(as.numeric(sv_df[['SV_size']]), as.numeric(sv_df[['SV_size']])))
        p = p + geom_segment(aes(x = GPOS_X, y = 0, xend = GPOS_X, yend = 2, color = SV_type), 
          data = sv2_df, size = sv2_df[['SV_size']], alpha = 0.5, show.legend = FALSE)
    }
    return(p)
}



library(optparse)
library(data.table)

V="Version: 1.0"
D="Depends: R (>= 3.4.0), optparse, data.table, ggplot2, grid, gridExtra, gtable"

a = commandArgs(trailingOnly = TRUE)

option_list = list(
    make_option(c("-s", "--sample_name"), type = "character", default = "tumor",
     help = "Sample name for output [%default]", metavar = "character"),
    make_option(c("-d", "--ref_dict"), type = "character", default = "",
     help = "Reference file in DICT format [%default]", metavar = "character"),
    make_option(c("-n", "--sv_tsv"), type = "character", default = "",
     help = "Structural variation data in TSV format. The first four columns must be \n\t\tCONTIG, POS, POS2, SV_type and SV_size (optinal) [%default]", metavar = "character"),
    make_option(c("-u", "--purity"), type = "double", default = "1",
     help = "Purity [%default]", metavar = "double"),
    make_option(c("-l", "--ploidy"), type = "double", default = "2",
     help = "Ploidy [%default]", metavar = "double"),
    make_option(c("-g", "--region"), type = "character", default = "",
     help = "The output region [%default] \n\t\te.g. g='chr17:36000000-42000000'", metavar = "character"),
    make_option(c("-o", "--output_prefix"), type = "character", default = "plotACN",
     help = "Output prefix [%default]", metavar = "character"))

opt_parser = OptionParser(usage = "usage: %prog [options] <copy_ratio.tsv> <het_count.tsv>\n\t
    These required input files should be in TSV format (separator: '\\t', stdin: '-') \n\t
    The first four columns of <copy_ratio.tsv> must be CONTIG, START, END and LOG2_COPY_RATIO. \n\t
    The first four columns of <het_count.tsv> must be CONTIG, POSITION, REF_COUNT and ALT_COUNT.", 
    option_list=option_list, description = paste(V, D, sep = "\n"))
opt = parse_args(opt_parser, args = a, positional_arguments = TRUE)

s = opt$options$sample_name
d = opt$options$ref_dict
n = opt$options$sv_tsv
u = opt$options$purity
l = opt$options$ploidy
g = opt$options$region
o = opt$options$output_prefix

r = ifelse(opt$args[1] == '-', 'file:///dev/stdin', opt$args[1]) # Data_path
e = ifelse(opt$args[2] == '-', 'file:///dev/stdin', opt$args[2]) # Data_path
n = ifelse(n == '-', 'file:///dev/stdin', n) # Data_path

# main
## read
### fread
r_dt = fread(r, skip = 'CONTIG', sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
e_dt = fread(e, skip = 'CONTIG', sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)
if (file.exists(n)) sv_dt = fread(n, skip = 'CONTIG', sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, data.table = TRUE, showProgress = FALSE, verbose = FALSE)

## format
### remove unnecessary columns in r_dt
r_dt = r_dt[, .(CONTIG, START, END, LOG2_COPY_RATIO)]
### generate required features in e_dt
e_dt[, BAF := min(REF_COUNT, ALT_COUNT)/sum(REF_COUNT, ALT_COUNT), by = .(CONTIG, POSITION)]
e_dt[, POSITION2 := POSITION]
### remove 'chr' in CONTIG column r_dt and e_dt
r_dt[['CONTIG']] = as.character(gsub('^chr', '', r_dt[['CONTIG']]))
e_dt[['CONTIG']] = as.character(gsub('^chr', '', e_dt[['CONTIG']]))
if (file.exists(n)) sv_dt[['CONTIG']] = as.character(gsub('^chr', '', sv_dt[['CONTIG']]))
### setkey for foverlaps
setkey(r_dt, CONTIG, START, END)
setkey(e_dt, CONTIG, POSITION, POSITION2)
if (file.exists(n)) setkey(sv_dt, CONTIG, POS, POS2)
### build region data.table
if (grepl('[0-9]+', g)) {
    g = gsub('^chr', '', g)
    g = unlist(strsplit(g, '[:-]'))
    if ((length(g) == 3) & (as.integer(g[2]) <= as.integer(g[3]))) {
        g_dt = data.table(CONTIG = g[1], START_g = as.integer(g[2]), END_g = as.integer(g[3]))
        setkey(g_dt, CONTIG, START_g, END_g)
        r_dt = foverlaps(r_dt, g_dt, type = 'within', nomatch=0L, by.x=c('CONTIG', 'START', 'END'))
    } else if (length(g) == 1) {
        g_dt = data.table(CONTIG = g[1])
        setkey(g_dt, CONTIG)
        r_dt = r_dt[CONTIG == g_dt[['CONTIG']]]
    } else {
        cat("ERROR: The value of '-g' must be either 'chr1' or 'chr1:1-100'\n")
        quit(save="no", status=1, runLast=FALSE)
    }
} else {
    g_dt = data.table()
}


## allelic CN
### average ploidy
D = u*l + 2*(1-u)
### map the het sites to each bin
er_dt = foverlaps(e_dt, r_dt, type = 'within', nomatch=0L, by.x=c('CONTIG', 'POSITION', 'POSITION2'))
### call mean BAF
mer_dt = er_dt[, .(mean(BAF)), by = .(CONTIG, START, END)]
colnames(mer_dt)[colnames(mer_dt) == 'V1'] = 'MEAN_BAF'
rem_dt = merge(r_dt, mer_dt, by = c('CONTIG', 'START', 'END'))
### call MAJOR_COPY_NUMBER and MINOR_COPY_NUMBER
rem_dt[, MAJOR_COPY_NUMBER := ((2^LOG2_COPY_RATIO) * D * (1 - MEAN_BAF) - 1) / u + 1]
rem_dt[, MINOR_COPY_NUMBER := ((2^LOG2_COPY_RATIO) * D * MEAN_BAF - 1) / u + 1]
rem_dt[, MAJOR_COPY_NUMBER := ifelse(MAJOR_COPY_NUMBER < 0, 0, MAJOR_COPY_NUMBER)]
rem_dt[, MINOR_COPY_NUMBER := ifelse(MINOR_COPY_NUMBER < 0, 0, MINOR_COPY_NUMBER)]
### convert to data.frame (*checkpoint)
rem_df = as.data.frame(rem_dt[, .(CONTIG, START, END, MAJOR_COPY_NUMBER, MINOR_COPY_NUMBER, MEAN_BAF)])
### determine copy-ratio midpoints
rem_df[['MIDDLE']] = round((rem_df[['START']] + rem_df[['END']]) / 2)


## ggplot
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)

### build contig table
if (file.exists(d)) {
    contig_dict_df = read.delim(d, header=FALSE)
    contig_names = gsub('SN:', '', contig_dict_df[contig_dict_df[,1] == '@SQ', 2])
    contig_lengths = setNames(as.integer(gsub('LN:', '', contig_dict_df[contig_dict_df[,1] == '@SQ', 3])), contig_names)
    if (ncol(g_dt) == 3) {
        contig_gl = as.integer(contig_lengths[[g_dt[['CONTIG']]]])
        if (g_dt[['END_g']] > contig_gl) g_dt[['END_g']] = contig_gl
        contig_names = g_dt[['CONTIG']]
        contig_ends = g_dt[['END_g']] - g_dt[['START_g']] + 1
        contig_starts = 0
    } else if (ncol(g_dt) == 1) {
        contig_gl = as.integer(contig_lengths[[g_dt[['CONTIG']]]])
        contig_names = g_dt[['CONTIG']]
        contig_ends = contig_gl
        contig_starts = 0
    } else {
        contig_ends = cumsum(contig_lengths)
        contig_starts = c(0, head(contig_ends, -1))
    }
} else {
    contig_names = unique(rem_dt[['CONTIG']])
    contig_lengths = rem_dt[, max(END), by = CONTIG][, structure(V1, names = contig_names)]
    if (ncol(g_dt) == 3) {
        contig_names = g_dt[['CONTIG']]
        contig_ends = g_dt[['END_g']] - g_dt[['START_g']] + 1
        contig_starts = 0
    } else if (ncol(g_dt) == 1) {
        contig_names = g_dt[['CONTIG']]
        contig_ends = contig_lengths[g_dt[['CONTIG']]]
        contig_starts = 0
    } else {
        contig_ends = cumsum(contig_lengths)
        contig_starts = c(0, head(contig_ends, -1))
    }
}
contig_centers = (contig_starts + contig_ends) / 2
contig_df = data.frame(names = contig_names, starts = contig_starts, ends = contig_ends, centers = contig_centers)
### transform features to colors
rem_df[['GPOS_X']] = contig_starts[match(rem_df[['CONTIG']], contig_names)] + rem_df[['MIDDLE']]
if (ncol(g_dt) == 3) rem_df[['GPOS_X']] = rem_df[['GPOS_X']] - g_dt[['START_g']] + 1
map = setNames(c('darkorange', 'darkcyan', 'red', 'blue', 'gold', 'black'),
               c('MAJOR_COPY_NUMBER', 'MINOR_COPY_NUMBER', 'AMP', 'DEL', 'INV', 'ITC'))
### determine genomic coordinates for CN
rem_gdf = melt(rem_df, id='GPOS_X', measure=c('MAJOR_COPY_NUMBER', 'MINOR_COPY_NUMBER'), 
               variable.name = 'CN_type', value.name = 'CN')
rem_gdf[['CN_type']] = map[unlist(as.character(rem_gdf[['CN_type']]))]
### determine genomic coordinates for SV
if (exists('sv_dt')) {
    if (!'SV_size' %in% colnames(sv_dt)) sv_dt[['SV_size']] = 0.5
    if (ncol(g_dt) == 3) sv_dt = foverlaps(sv_dt, g_dt, type = 'within', nomatch=0L, by.x=c('CONTIG', 'POS', 'POS2'))
    if (ncol(g_dt) == 1) sv_dt = sv_dt[CONTIG == g_dt[['CONTIG']]]
    sv_df = as.data.frame(sv_dt[, .(CONTIG, POS, POS2, SV_type, SV_size)])
    sv_df[['GPOS_X']] = contig_starts[match(sv_df[['CONTIG']], contig_names)] + sv_df[['POS']]
    sv_df[['GPOS2_X']] = contig_starts[match(sv_df[['CONTIG']], contig_names)] + sv_df[['POS2']]
    if (ncol(g_dt) == 3) sv_df[['GPOS_X']] = sv_df[['GPOS_X']] - as.numeric(g_dt[['START_g']]) + 1
    if (ncol(g_dt) == 3) sv_df[['GPOS2_X']] = sv_df[['GPOS2_X']] - as.numeric(g_dt[['START_g']]) + 1
    sv_df[['SV_type']] = map[unlist(as.character(sv_df[['SV_type']]))]
} else {
    sv_df = data.table()
}
### ggplotGrob
cn_max = 2^ceiling(log2(max(max(rem_df[['MAJOR_COPY_NUMBER']]), max(rem_df[['MINOR_COPY_NUMBER']]))))
pdf(paste0(o, '.plotACN.pdf'), 12, 7)
    if (nrow(sv_df) > 0) {
        top_p = top_plot(sv_df, contig_df, cn_max)
        mid_p = mid_plot(cn_df = rem_gdf, sv_df = sv_df, contig_df = contig_df, cn_max = cn_max)
        bot_p = bot_plot(cn_df = rem_gdf, sv_df = sv_df, contig_df = contig_df, cn_max = cn_max)
        g0 <- ggplotGrob(top_p); g1 <- ggplotGrob(mid_p); g2 <- ggplotGrob(bot_p)
        g0$heights[unique(g0$layout[g0$layout$name == 'panel', 't'])] = unit(4, 'cm')
        g1$heights[unique(g1$layout[g1$layout$name == 'panel', 't'])] = unit((ceiling(log2(cn_max))-1)/2*2, 'cm')
        g2$heights[unique(g2$layout[g2$layout$name == 'panel', 't'])] = unit(2, 'cm')
        gg <- rbind(g0, g1, g2, size = 'first')
    } else {
        mid_p = mid_plot(cn_df = rem_gdf, contig_df = contig_df, cn_max = cn_max)
        bot_p = bot_plot(cn_df = rem_gdf, contig_df = contig_df, cn_max = cn_max)
        g1 <- ggplotGrob(mid_p); g2 <- ggplotGrob(bot_p)
        g1$heights[unique(g1$layout[g1$layout$name == 'panel', 't'])] = unit((ceiling(log2(cn_max))-1)/2*2, 'cm')
        g2$heights[unique(g2$layout[g2$layout$name == 'panel', 't'])] = unit(2, 'cm')
        gg <- rbind(g1, g2, size = 'first')
    }
    grid.newpage()
    grid.draw(gg)
dev.off()

write.table(rem_df, paste0(o, '.plotACN.tsv'), sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)


# test at checkpoint
#rem_dt = data.table(CONTIG=c(rep(1,5),2,2,2,3,3),
#    START=c(1,2,6,10,12,18,3,11,2,7),END=c(3,5,8,11,14,22,8,17,6,12),
#    MAJOR_COPY_NUMBER=c(0,0.8,1,1.1,1.8,2,2.2,4,22,54), MINOR_COPY_NUMBER=c(0,0.4,0.8,1,0.86,0.9,0.85,0.76,11,20))
#sv_dt = data.table(CONTIG=c(rep(1,4),2,2,2,2,3,3),
#    POS=c(2,3,5,10,5,7,8,11,5,7),POS2=c(2,6,5,11,12,9,15,17,12,9),
#    SV_type=c('DEL','AMP','AMP','DEL','DEL','INV','DEL','AMP','AMP','INV'))
#d=''

#s = 'EAC-9_11-EAC1'
#r = '/cga/bass/Chunyang/task/Matthew/WGS_EAC/CNVSomaticPairWorkflow_v4p0p8p0_50k/denoised_copy_ratios_tumor/9_11_EAC1.denoisedCR.tsv'
#e = '/cga/bass/Chunyang/task/Matthew/WGS_EAC/CNVSomaticPairWorkflow_v4p0p8p0_50k/het_allelic_counts_tumor/9_11_EAC1.hets.tsv'
#n = '/cga/bass/Chunyang/task/Matthew/WGS_EAC/SV_evolution/CZ_bkps__plotACV_tsv/EAC-9_11-EAC1.plotACN_sv.tsv'
#d = ''
#u = 0.4
#l = 3.85
#o = '/cga/bass/Chunyang/task/Matthew/WGS_EAC/SV_evolution/CZ_bkps__plotACV/plotACN'
#g = 'chr17:35000000-45000000'
#d = '/cga/bass/Chunyang/ref/hg19/Homo_sapiens_assembly19.dict'

### filter by absolute copy number
#D = u*l + 2*(1-u)
#r_dt[, COPY_NUMBER := ((2^LOG2_COPY_RATIO)*D-2*(1-u))/(u)]
#r_dt = r_dt[COPY_NUMBER > 0]

#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.0"
D="Depends: R (>= 3.4.0), optparse, data.table"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
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

opt_parser <- OptionParser(usage="usage: %prog [options] <D.tsv>\n\t<D.tsv> is a TSV file for GISTIC D data (row: event; column: feature, separator: '\\t', stdin: '-').", option_list=option_list, description = paste(V, D, sep="\n"))
opt <- parse_args(opt_parser, args=a, positional_arguments=TRUE)

cs <- opt$options$col_sep
rn <- opt$options$row_num
cn <- as.logical(opt$options$col_name)
na <- opt$options$na_str
rs <- opt$options$row_skip
ra <- opt$options$row_name

if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)
if(grepl("^[[:digit:]]+$", ra)) ra=as.numeric(ra)

d <- ifelse(opt$args[1]=='-', 'file:///dev/stdin', opt$args[1]) # D

# read
dta <- fread(d, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE) # D_Table
# main
add_SCNA_stats_from_D <- function(t_D, t_D_sample_colname='sample')
{
  dt_D = data.table(t_D)
  dt_D$amplitude_score = 0
  dt_D$amplitude_score[dt_D$amplitude > 0.3] = 1
  dt_D$amplitude_score[dt_D$amplitude > 1] = 2
  SCNA_stats = data.frame(sample = levels(as.factor(as.data.frame(dt_D)[, t_D_sample_colname])), N_chr = NA)
  
  # total number of chr-level events
  temp = dt_D[(dt_D$amplitude_score > 0 & arm_length >= 1.5), .N, by = c(t_D_sample_colname)]
  SCNA_stats$N_chr = as.data.frame(temp)[match(SCNA_stats[, 1], as.data.frame(temp)[, 1]), 2]
  SCNA_stats$N_chr[is.na(SCNA_stats$N_chr)] = 0
  # total number of arm-level events
  temp = dt_D[(dt_D$amplitude_score > 0 & arm_length < 1.5 & arm_length >= 0.5), .N, by = c(t_D_sample_colname)]
  SCNA_stats$N_arm = as.data.frame(temp)[match(SCNA_stats[, 1], as.data.frame(temp)[, 1]), 2]
  SCNA_stats$N_arm[is.na(SCNA_stats$N_arm)] = 0
  
  SCNA_stats$N_chr_and_arm = SCNA_stats$N_chr + SCNA_stats$N_arm
  
  # total number of arm-level amp events
  temp = dt_D[(dt_D$amplitude_score > 0 & arm_length < 0.5 & (event_type == 'amp' | event_type == 'aod')), .N, by = c(t_D_sample_colname)]
  SCNA_stats$N_focal_amp = as.data.frame(temp)[match(SCNA_stats[, 1], as.data.frame(temp)[, 1]), 2]
  SCNA_stats$N_focal_amp[is.na(SCNA_stats$N_focal_amp)] = 0
  # total number of arm-level del events
  temp = dt_D[(dt_D$amplitude_score > 0 & arm_length < 0.5 & (event_type == 'del' | event_type == 'doa')), .N, by = c(t_D_sample_colname)]
  SCNA_stats$N_focal_del = as.data.frame(temp)[match(SCNA_stats[, 1], as.data.frame(temp)[, 1]), 2]
  SCNA_stats$N_focal_del[is.na(SCNA_stats$N_focal_del)] = 0
  
  SCNA_stats$N_focal = SCNA_stats$N_focal_amp + SCNA_stats$N_focal_del
  
  # total number of focal-level amp events excluding short
  temp = dt_D[(dt_D$amplitude_score > 0 & arm_length < 0.5 & (base_end-base_start) > 3e+6), .N, by = c(t_D_sample_colname)]
  SCNA_stats$N_focal_excluding_short = as.data.frame(temp)[match(SCNA_stats[, 1], as.data.frame(temp)[, 1]), 2]
  SCNA_stats$N_focal_excluding_short[is.na(SCNA_stats$N_focal_excluding_short)] = 0
  # total number of focal-level high-amp events
  temp = dt_D[(dt_D$amplitude_score == 2 & arm_length < 0.5 & (event_type == 'amp' | event_type == 'aod')), .N, by = c(t_D_sample_colname)]
  SCNA_stats$N_focal_amp_high = as.data.frame(temp)[match(SCNA_stats[, 1], as.data.frame(temp)[, 1]), 2]
  SCNA_stats$N_focal_amp_high[is.na(SCNA_stats$N_focal_amp_high)] = 0
  # total number of focal-level deep-del events
  temp = dt_D[(dt_D$amplitude_score == 2 & arm_length < 0.5 & (event_type == 'del' | event_type == 'doa')), .N, by = c(t_D_sample_colname)]
  SCNA_stats$N_focal_del_deep = as.data.frame(temp)[match(SCNA_stats[, 1], as.data.frame(temp)[, 1]), 2]
  SCNA_stats$N_focal_del_deep[is.na(SCNA_stats$N_focal_del_deep)] = 0
  
  # total score of chr-level events
  temp = dt_D[(dt_D$amplitude_score > 0 & arm_length >= 1.5), sum(amplitude_score), by = c(t_D_sample_colname)]
  SCNA_stats$L_chr = as.data.frame(temp)[match(SCNA_stats[, 1], as.data.frame(temp)[, 1]), 2]
  SCNA_stats$L_chr[is.na(SCNA_stats$L_chr)] = 0
  # total score of arm-level events
  temp = dt_D[(dt_D$amplitude_score > 0 & arm_length < 1.5 & arm_length >= 0.5), sum(amplitude_score), by = c(t_D_sample_colname)]
  SCNA_stats$L_arm = as.data.frame(temp)[match(SCNA_stats[, 1], as.data.frame(temp)[, 1]), 2]
  SCNA_stats$L_arm[is.na(SCNA_stats$L_arm)] = 0
  
  SCNA_stats$L_chr_and_arm = SCNA_stats$L_chr + SCNA_stats$L_arm
  
  # total score of arm-level amp events
  temp = dt_D[(dt_D$amplitude_score > 0 & arm_length < 0.5 & (event_type == 'amp' | event_type == 'aod')), sum(amplitude_score), by = c(t_D_sample_colname)]
  SCNA_stats$L_focal_amp = as.data.frame(temp)[match(SCNA_stats[, 1], as.data.frame(temp)[, 1]), 2]
  SCNA_stats$L_focal_amp[is.na(SCNA_stats$L_focal_amp)] = 0
  # total score of arm-level del events
  temp = dt_D[(dt_D$amplitude_score > 0 & arm_length < 0.5 & (event_type == 'del' | event_type == 'doa')), sum(amplitude_score), by = c(t_D_sample_colname)]
  SCNA_stats$L_focal_del = as.data.frame(temp)[match(SCNA_stats[, 1], as.data.frame(temp)[, 1]), 2]
  SCNA_stats$L_focal_del[is.na(SCNA_stats$L_focal_del)] = 0
  
  SCNA_stats$L_focal = SCNA_stats$L_focal_amp + SCNA_stats$L_focal_del

  # total score of arm-level extreme-amp events
  temp = dt_D[(dt_D$amplitude > 2 & arm_length < 0.5 & (base_end-base_start) < 3e+6 & (event_type == 'amp' | event_type == 'aod')), sum((base_end-base_start)/1e+6*amplitude), by = c(t_D_sample_colname)]
  SCNA_stats$Score_focal_amp_extreme = as.data.frame(temp)[match(SCNA_stats[, 1], as.data.frame(temp)[, 1]), 2]
  SCNA_stats$Score_focal_amp_extreme[is.na(SCNA_stats$Score_focal_amp_extreme)] = 0
  
  SCNA_stats$Score_focal_amp_extreme_class = ifelse(SCNA_stats$Score_focal_amp_extreme > 2, 1, 0)
  
  return(SCNA_stats)
}

st = add_SCNA_stats_from_D(dta) # Score_Table

# output
write.table(st, quote=F, sep="\t", append=FALSE, col.names=TRUE, row.names=FALSE)

### THE END ###

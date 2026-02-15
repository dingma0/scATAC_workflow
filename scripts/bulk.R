library(Signac)
library(Seurat)
library(GenomicRanges)
library(GenomeInfoDb)
library(optparse)
library(glue)
library(stringr)

option_list = list(make_option("--objs", default = NULL),
                   make_option("--out_dir", default = NULL))
args = parse_args(OptionParser(option_list = option_list))

vec <- str_split(args$objs, " ")[[1]]
obj_list <- lapply(vec, readRDS)

all.peaks <- lapply(obj_list, granges)
combined.peaks <- reduce(do.call(c, all.peaks))
combined.peaks <- keepStandardChromosomes(combined.peaks)
w <- width(combined.peaks)
combined.peaks <- combined.peaks[w < 10000 & w > 20]

saveRDS(combined.peaks, glue("{args$out_dir}/combined_peaks.rds"))
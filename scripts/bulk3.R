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

combined <- merge(
  obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = unname(sapply(obj_list, function(obj) {obj$ss_id[1]})),
  project = "SyS_Batch_2" # todo
)

saveRDS(combined, glue("{args$out_dir}/combined_atac.rds"))
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

combined <- readRDS(vec[1])
combined <- RenameCells(combined, add.cell.id = combined$ss_id[1])

if (length(vec) > 1) {
  for (i in 2:length(vec)) {
    obj <- readRDS(vec[i])
    obj <- RenameCells(obj, add.cell.id = obj$ss_id[1])

    combined <- merge(
      x = combined,
      y = obj,
      project = "SyS_Batch_2"
    )

    rm(obj)
    gc()
  }
}

saveRDS(combined, glue("{args$out_dir}/combined_atac.rds"))
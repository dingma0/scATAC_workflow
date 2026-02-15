library(Signac)
library(Seurat)
library(GenomicRanges)
library(GenomeInfoDb)
library(optparse)
library(glue)
library(stringr)

option_list = list(make_option("--sample_id", default = NULL),
                   make_option("--obj", default = NULL),
                   make_option("--peaks_obj", default = NULL),
                   make_option("--out_dir", default = NULL))
args = parse_args(OptionParser(option_list = option_list))

obj <- readRDS(args$obj)
combined.peaks <- readRDS(args$peaks_obj)
frags <- Fragments(obj)

counts <- FeatureMatrix(
  fragments = frags,
  features = combined.peaks,
  cells = colnames(obj)
)

obj[["peaks_consensus"]] <- CreateChromatinAssay(
  counts = counts,
  fragments = frags
)

DefaultAssay(obj) <- "peaks_consensus"

saveRDS(obj, glue("{args$out_dir}/{args$sample_id}_atac_consensus.rds"))
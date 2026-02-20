library(Signac)
library(Seurat)
library(GenomicRanges)
library(GenomeInfoDb)
library(glue)
library(stringr)
library(ggplot2)
library(optparse)

option_list = list(make_option("--comb", default = NULL),
                   make_option("--out_dir", default = NULL))
args = parse_args(OptionParser(option_list = option_list))

combined <- readRDS(args$comb)

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')

umap_ss <- DimPlot(combined, group.by = 'ss_id')
ggsave(glue("{args$out_dir}/combined_atac_umap.png"), umap_ss)

combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3, resolution = 0.4)
ggsave(glue("{args$out_dir}/combined_atac_umap_clustered.png"), DimPlot(combined))

agg_peaks <- AggregateExpression(
  object = combined,
  assays = "peaks",
  group.by = "ss_id",
  slot = "counts"
)

saveRDS(combined, glue("{args$out_dir}/combined_atac_proc.rds"))
saveRDS(agg_peaks, glue("{args$out_dir}/combined_atac_agg.rds"))

library(Signac)
library(Seurat)
library(optparse)
library(glue)
library(ggplot2)

option_list = list(make_option("--sample_id", default = NULL),
                   make_option("--obj", default = NULL),
                   make_option("--out_dir", default = NULL))
args = parse_args(OptionParser(option_list = option_list))

obj <- readRDS(args$obj)

obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(obj)

obj <- RunUMAP(object = obj, reduction = 'lsi', dims = 2:30)
obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:30)
obj <- FindClusters(object = obj, verbose = FALSE, algorithm = 3, resolution = 0.4)

vln <- VlnPlot(object = obj, features = c('nCount_peaks', 'TSS.enrichment',
                                    'blacklist_ratio', 'nucleosome_signal',
                                    'pct_reads_in_peaks'),
        pt.size = 0.1, ncol = 5, group.by = "seurat_clusters")

ggsave(glue("{args$out_dir}/{args$sample_id}_atac_post_vln.png"))
ggsave(glue("{args$out_dir}/{args$sample_id}_atac_umap.png"), DimPlot(object = obj))

# gene.activities <- GeneActivity(obj)

# obj[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
# obj <- NormalizeData(
#   object = obj,
#   assay = 'ACTIVITY',
#   normalization.method = 'LogNormalize',
#   scale.factor = median(obj$nCount_ACTIVITY)
# )

saveRDS(obj, glue("{args$out_dir}/{args$sample_id}_atac.rds"))

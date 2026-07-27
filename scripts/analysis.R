library(Signac)
library(Seurat)
library(optparse)
library(glue)
library(ggplot2)

option_list = list(make_option("--sample_id", default = NULL),
                   make_option("--obj", default = NULL),
                   make_option("--out_dir", default = NULL),
                   make_option("--params", default = NULL))
args = parse_args(OptionParser(option_list = option_list))

obj <- readRDS(args$obj)
params <- read.csv(args$params)

obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(obj)

elbow <- ElbowPlot(obj, ndims = 30, reduction = "lsi")
ggsave(glue("{args$out_dir}/{args$sample_id}_atac_elbow.png"), elbow,
        units = "px", height = 400, width = 600, scale = 5)

obj <- RunUMAP(object = obj, reduction = 'lsi', dims = 2:params$n_dims)
obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:params$n_dims)
obj <- FindClusters(object = obj, verbose = FALSE, algorithm = 3, resolution = params$clust_res)

vln <- VlnPlot(object = obj, features = c('nCount_peaks', 'TSS.enrichment',
                                    'blacklist_ratio', 'nucleosome_signal',
                                    'pct_reads_in_peaks'),
        pt.size = 0.1, ncol = 5, group.by = "seurat_clusters")

ggsave(glue("{args$out_dir}/{args$sample_id}_atac_post_vln.png"), vln, 
        units = "px", height = 400, width = 800, scale = 5)

umap <- DimPlot(object = obj, label = TRUE) + NoLegend()
ggsave(glue("{args$out_dir}/{args$sample_id}_atac_umap.png"), umap,
       units = "px", height = 400, width = 400, scale = 5)

# gene.activities <- GeneActivity(obj)

# obj[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
# obj <- NormalizeData(
#   object = obj,
#   assay = 'ACTIVITY',
#   normalization.method = 'LogNormalize',
#   scale.factor = median(obj$nCount_ACTIVITY)
# )

saveRDS(obj, glue("{args$out_dir}/{args$sample_id}_atac.rds"))

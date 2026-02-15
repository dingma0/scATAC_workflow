library(Signac)
library(Seurat)
library(AnnotationHub)
library(GenomicRanges)
library(GenomeInfoDb)
library(optparse)
library(glue)
library(ggplot2)

option_list = list(make_option("--sample_id", default = NULL),
                   make_option("--bivalency", default = NULL),
                   make_option("--h5", default = NULL),
                   make_option("--meta", default = NULL),
                   make_option("--frags", default = NULL),
                   make_option("--annots", default = NULL),
                   make_option("--out_dir", default = NULL))
args = parse_args(OptionParser(option_list = option_list))

counts <- Read10X_h5(filename = args$h5)
metadata <- read.csv(
    file = args$meta,
    header = TRUE,
    row.names = 1
)

chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    fragments = args$frags,
    min.cells = 10,
    min.features = 200
)

obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
)

peaks.keep <- seqnames(granges(obj)) %in% standardChromosomes(granges(obj))
obj <- obj[as.vector(peaks.keep), ]
obj <- NucleosomeSignal(object = obj)

# ah <- AnnotationHub()
# ensdb_v98 <- ah[["AH75011"]]
# annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
# seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
# genome(annotations) <- "hg38"
annotations <- readRDS(args$annots)
Annotation(obj) <- annotations

obj <- NucleosomeSignal(object = obj)
obj <- TSSEnrichment(object = obj)

obj$pct_reads_in_peaks <- obj$peak_region_fragments / obj$passed_filters * 100
obj$blacklist_ratio <- FractionCountsInRegion(
    object = obj, 
    assay = 'peaks',
    regions = blacklist_hg38_unified
)

qc_vln <- VlnPlot(
            object = obj,
            features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
            pt.size = 0.1,
            ncol = 5)
ggsave(glue("{args$out_dir}/{args$sample_id}_atac_pre_vln.png"), qc_vln)

is_outlier <- function(adata, metric, nmads) {
    M <- adata@meta.data[[metric]]
    outlier <- (M < median(M) - nmads * mad(M)) | (median(M) + nmads * mad(M) < M)
    return(outlier)
}

obj@meta.data$outlier <- is_outlier(obj, "nCount_peaks", 5)
obj <- subset(x = obj, subset = blacklist_ratio < 0.01 &
                nucleosome_signal < 4 &
                outlier == FALSE)

obj$ss_id <- args$sample_id
obj$bivalency <- args$bivalency

saveRDS(obj, glue("{args$out_dir}/{args$sample_id}_atac_pre.rds"))


process SIGNAC_PREPO {
    publishDir "$workflow.outputDir", pattern:"signac/$sample_id/*.png", mode:'copy'

    input:
    tuple val(sample_id), val(bivalency), path(h5), path(meta), path(frags), path(frags_idx)
    each path(script)
    each path(annots)

    output:
    tuple val(sample_id), path("signac/$sample_id/${sample_id}_atac_pre.rds"), path(frags), path(frags_idx), emit: data
    path("signac/$sample_id/${sample_id}_atac_pre_vln.png"), emit: img

    script:
    """
    mkdir -p signac/$sample_id
    export HOME=\$(pwd)

    Rscript $script \
        --sample_id $sample_id \
        --bivalency $bivalency \
        --h5 $h5 \
        --meta $meta \
        --frags $frags \
        --annots $annots \
        --out_dir signac/$sample_id \
    """
}

process SIGNAC_ANALYSIS {
    publishDir "$workflow.outputDir", mode:'copy'

    input:
    tuple val(sample_id), path(prepo), path(frags), path(frags_idx)
    each path(script)

    output:
    tuple val(sample_id), path("signac/$sample_id/${sample_id}_atac.rds"), emit: data
    tuple path("signac/$sample_id/${sample_id}_atac_post_vln.png"), path("signac/$sample_id/${sample_id}_atac_umap.png"), emit: img

    script:
    """
    mkdir -p signac/$sample_id

    Rscript $script \
        --sample_id $sample_id \
        --obj $prepo \
        --out_dir signac/$sample_id
    """
}

process SIGNAC_BULK {
    input:
    path(objs)
    each path(script)

    output:
    path("signac/combined_peaks.rds")

    script:
    """
    mkdir -p signac

    Rscript $script \
        --objs "$objs" \
        --out_dir signac
    """
}

process SIGNAC_BULK2 {
    input:
    tuple val(sample_id), path(obj), path(frags), path(frags_idx)
    each path(peaks)
    each path(script)

    output:
    tuple val(sample_id), path("signac/$sample_id/${sample_id}_atac_consensus.rds"), path(frags), path(frags_idx)

    script:
    """
    mkdir -p signac/$sample_id

    Rscript $script \
        --sample_id $sample_id \
        --obj $obj \
        --peaks_obj $peaks \
        --out_dir signac/$sample_id
    """
}

process SIGNAC_BULK3 {
    publishDir "$workflow.outputDir", mode:'copy'

    input:
    path(objs)
    path(frags)
    path(frags_idx)
    each path(script)

    output:
    path("signac/combined_atac.rds")

    script:
    """
    mkdir -p signac

    Rscript $script \
        --objs "$objs" \
        --out_dir signac
    """
}

workflow {
    main:
    samples = channel.fromPath(params.sample_sheet).splitCsv(header:true)

    prepo_script = channel.fromPath("scripts/prepo.R")
    annots = channel.fromPath("data/annotations.rds")
    SIGNAC_PREPO(samples, prepo_script, annots)

    analysis_script = channel.fromPath("scripts/analysis.R")
    SIGNAC_ANALYSIS(SIGNAC_PREPO.out.data, analysis_script)

    bulk_script = channel.fromPath("scripts/bulk.R")
    SIGNAC_BULK(SIGNAC_PREPO.out.data.collect { tuple -> tuple[1] }, bulk_script)

    bulk_script2 = channel.fromPath("scripts/bulk2.R")
    SIGNAC_BULK2(SIGNAC_PREPO.out.data, SIGNAC_BULK.out, bulk_script2)

    bulk_script3 = channel.fromPath("scripts/bulk3.R")
    SIGNAC_BULK3(SIGNAC_BULK2.out.collect { tuple -> tuple[1] },
                 SIGNAC_BULK2.out.collect { tuple -> tuple[2] },
                 SIGNAC_BULK2.out.collect { tuple -> tuple[3] }, bulk_script3)
}

process SIGNAC_PREPO {
    publishDir "$workflow.outputDir", pattern:"signac/$sample_id/*.png", mode:'copy'

    input:
    tuple val(sample_id), path(h5), path(meta), path(frags), path(frags_idx)
    path script
    path annots
    path sc_params

    output:
    tuple val(sample_id), path("signac/$sample_id/${sample_id}_atac_pre.rds"), path(frags), path(frags_idx), emit: data
    path("signac/$sample_id/${sample_id}_atac_*.png"), emit: img

    script:
    """
    mkdir -p signac/$sample_id
    export HOME=\$(pwd)

    Rscript $script \
        --sample_id $sample_id \
        --h5 $h5 \
        --meta $meta \
        --frags $frags \
        --annots $annots \
        --out_dir signac/$sample_id \
        --do_qc $params.qc_threshold \
        --params $sc_params
    """
}

process SIGNAC_ANALYSIS {
    publishDir "$workflow.outputDir", mode:'copy'

    input:
    tuple val(sample_id), path(prepo), path(frags), path(frags_idx)
    path script
    path sc_params

    output:
    tuple val(sample_id), path("signac/$sample_id/${sample_id}_atac.rds"), emit: data
    path("signac/$sample_id/${sample_id}_atac_*.png"), emit: img

    script:
    """
    mkdir -p signac/$sample_id

    Rscript $script \
        --sample_id $sample_id \
        --obj $prepo \
        --out_dir signac/$sample_id \
        --params $sc_params
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

process SIGNAC_AGGRE {
    publishDir "$workflow.outputDir", mode:'copy'

    input:
    path(comb)
    path(script)

    output:
    tuple path("signac/combined_atac_agg.rds"), path("signac/combined_atac_proc.rds"), emit: data
    tuple path("signac/combined_atac_umap.png"), path("signac/combined_atac_umap_clustered.png"), emit: img

    script:
    """
    mkdir -p signac

    Rscript $script \
        --comb $comb \
        --out_dir signac
    """
}

workflow {
    main:
    samples = channel.fromPath(params.sample_sheet).splitCsv(header:true)

    prepo_script = channel.value(file("scripts/prepo.R"))
    annots = channel.value(file("data/annotations.rds"))
    sc_params = channel.value(file(params.sc_params))
    SIGNAC_PREPO(samples, prepo_script, annots, sc_params)

    if (params.single_cell) {
        analysis_script = channel.value(file("scripts/analysis.R"))
        SIGNAC_ANALYSIS(SIGNAC_PREPO.out.data, analysis_script, sc_params)
    }

    if (params.pseudobulk) {
        bulk_script = channel.fromPath("scripts/bulk.R")
        SIGNAC_BULK(SIGNAC_PREPO.out.data.collect { tuple -> tuple[1] }, bulk_script)

        bulk_script2 = channel.fromPath("scripts/bulk2.R")
        SIGNAC_BULK2(SIGNAC_PREPO.out.data, SIGNAC_BULK.out, bulk_script2)

        bulk_script3 = channel.fromPath("scripts/bulk3.R")
        SIGNAC_BULK3(SIGNAC_BULK2.out.collect { tuple -> tuple[1] },
                    SIGNAC_BULK2.out.collect { tuple -> tuple[2] },
                    SIGNAC_BULK2.out.collect { tuple -> tuple[3] }, bulk_script3)

        agg_script = channel.fromPath("scripts/aggregate.R")
        SIGNAC_AGGRE(SIGNAC_BULK3.out, agg_script)
    }
}

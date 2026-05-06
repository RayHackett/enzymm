#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/////// Expected parameters
// params.outdir
// params.input
// params.batch_size

// profile options: local or cluster
// USEAGE: nextflow run template_matcher.nf -profile <profile>

process matching {
    tag { meta.id }

    input:
    tuple(
        val(meta),
        path(files, stageAs: 'inputs/*')
    )

    output:
    tuple(
        val(meta),
        path("${ meta.id }.parquet"),
        path("${ meta.id }.residues.parquet"), 
        path("${ meta.id }.tsv"),
        path("${ meta.id }.residues.tsv"),
        path("${ meta.id }.transformations.npz")
    )

    script:
    def args = []

    // controlling boolean flags
    if (params.verbose)                args << "-v"
    if (params.warn)                   args << "-w"
    if (params.match_small_templates)  args << "--match-small-templates"
    if (params.skip_smaller_hits)      args << "--skip-smaller-hits"
    if (params.skip_annotation)        args << "--skip-annotation"
    if (params.unfiltered)             args << "--unfiltered"
    if (params.write_parquet)          args << "--write-parquet"
    if (params.save_transformations)   args << "--save-transformations"
    if (params.simple_results)         args << "--simple-results"
    if (params.per_residue_results)    args << "--per-residue-results"

    // controlling arguments
    if (params.conservation_cutoff)    args << "-c ${params.conservation_cutoff}"
    if (params.template_folder)        args << "-c ${params.template_folder}"
    if (params.n_cpus)                 args << "-j ${params.n_cpus}"
    if (params.rmsd && params.distance && params.max_dynamic_distance){
        args << " -p ${params.rmsd} ${params.distance} ${params.max_dynamic_distance}"
    }

    """
    ls inputs/ | xargs -I{} echo inputs/{} > input_list.txt
    enzymm -l input_list.txt -o ${meta.id} ${args.join(' ')}
    touch ${ meta.id }.residues.parquet
    touch ${ meta.id }.residues.tsv
    touch ${ meta.id }.parquet
    touch ${ meta.id }.tsv
    touch ${ meta.id }.transformations.npz
    """
}

process merge_parquet {

    tag { meta.id }

    publishDir params.outdir, mode: 'copy'

    input:
    tuple(
        val(meta),
        path(tables)
    )

    output:
    tuple(
        val(meta),
        path("merged.${ meta.type }.parquet")
    )

    script:
    """
    python - <<'EOF'
    import polars
    from pathlib import Path

    files = [str(p) for p in Path('.').glob("*.parquet")]

    if files:
        polars.scan_parquet(files).sink_parquet("merged.${ meta.type }.parquet")
    EOF
    """
}

process merge_npz {

    tag { meta.id }

    publishDir params.outdir, mode: 'copy'

    input:
    tuple(
        val(meta),
        path(tfms)
    )

    output:
    tuple(
        val(meta),
        path("merged.${ meta.type }.npz")
    )

    script:
    """
    python - <<'EOF'
    import numpy as np
    from pathlib import Path

    files = [str(p) for p in Path('.').glob("*${ meta.type }.npz")]

    if files:
        merged = {}
        for f in files:
            data = np.load(f)
            for key in data.files:
                if key in merged:
                    raise ValueError(f"Duplicate key found: {key}")
                merged[key] = data[key]

        np.savez("merged.${ meta.type }.npz", **merged)
    EOF
    """
}

process publish_tsv {
    tag { tsv.name }
    publishDir params.outdir, mode: 'copy'
    input:
        path(tsv)
    output:
        path(tsv)
    script:
        "true"
}

workflow {
    /*
    * run with:
    * module load bioinformatics/tools/Nextflow/25.04.6 # LUMC cluster
    * module load container/apptainer/1.4.1/gcc-8.5.0 # LUMC cluster
    * NFDIR=/exports/archive/lucid-grpzeller-primary/hackett/template_matching/nextflow
    * WORKDIR=/exports/lucid-grpzeller-work/hackett/swissprot_work
    * nextflow run $NFDIR -profile cluster_lumc -params-file local_params.yml -with-tower -resume -work-dir $WORKDIR
    */
    if (!params.input) error "Missing required parameter: input"
    if (!params.outdir) error "Missing required parameter: outdir"
    if (!params.batch_size) error "Missing required parameter: batch_size"
    /*
    * Create a channel emitting chunks of file paths from params.input,
    * each line is one file path
    * items in the channel 'ch_files' contains a tuple( meta.id, list of file objects).
    */
    Channel
        .fromPath(params.input)
        .splitText()
        .map { file(it.trim()) }
        .collate(params.batch_size)
        .map { files ->
            tuple([id: "batch_${files.hashCode()}"], files)
        }
        .set { ch_files }

    ch_files = matching(ch_files)

    merged_ch = ch_files
    .map { meta, parquet, residue_parquet, tsv, residue_tsv ->
        def key = [ id: 'merged' ]
        tuple(key, parquet, residue_parquet, tsv, residue_tsv)
    }
    .groupTuple(by: 0)

    parquet_ch = merged_ch.map { key, parquet, residue_parquet, tsv, residue_tsv, tfm ->
        tuple([id: key.id, type: 'results'], parquet)
    }

    res_parquet_ch = merged_ch.map { key, parquet, residue_parquet, tsv, residue_tsv, tfm ->
        tuple([id: key.id, type: 'residues'], residue_parquet)
    }

    tsv_ch = merged_ch.map { key, parquet, residue_parquet, tsv, residue_tsv, tfm ->
        tsv
    }

    res_tsv_ch = merged_ch.map { key, parquet, residue_parquet, tsv, residue_tsv, tfm ->
        residue_tsv
    }

    tfm_ch = merged_ch.map { key, parquet, residue_parquet, tsv, residue_tsv, tfm ->
        tuple([id: key.id, type: 'transformations'], tfm)
    }

    // merge parquet outputs
    if (params.write_parquet) {
        def merged_parquet_ch = parquet_ch
        // if per_residue_results, then mixin the residue tsv channel
        if (params.per_residue_results) {
            merged_parquet_ch = merged_parquet_ch.mix(res_parquet_ch)
        }
        merge_parquet(merged_parquet_ch)
    } else {
        // merge tsv outputs
        merged_tsv = tsv_ch.collectFile(name: 'merged.tsv', keepHeader: true, skip: 2)
        def publish_ch = merged_tsv
        // if per_residue_results, then mixin the residue tsv channel
        if (params.per_residue_results) {
            merged_residue_tsv = res_tsv_ch.collectFile(name: 'merged.residues.tsv', keepHeader: true, skip: 2)
            publish_ch = publish_ch.mix(merged_residue_tsv)
        }
        publish_tsv(publish_ch)
    }

    if (params.save_transformations) {
        merge_npz(tfm_ch)
    }
}
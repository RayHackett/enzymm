#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/////// Expected parameters
// params.output
// params.input
// params.batch_size
if (!params.input) error "Missing required parameter: --input"
if (!params.output) error "Missing required parameter: --output"
if (!params.batch_size) error "Missing required parameter: --batch_size"

//////// Optional parameters
def matcher_params = ""
matcher_params += params.template_folder ? " -t ${params.template_folder}" : ""

///////// Jess control parameters
def jess_params = ""
if (params.rmsd && params.distance && params.max_dynamic_distance){
    jess_params = " -j ${params.rmsd} ${params.distance} ${params.max_dynamic_distance}"
}


///////// optional control parameters
matcher_params += params.conservation_cutoff ? " -c ${params.conservation_cutoff}" : ""
matcher_params += params.warn ? " -w" : ""
matcher_params += params.verbose ? " -v" : ""
matcher_params += params.match_small_templates ? " --match-small-templates" : ""
matcher_params += params.skip_smaller_hits ? " --skip-smaller-hits" : ""
matcher_params += params.skip_annotation ? " --skip-annotation": ""
matcher_params += params.unfiltered ? " --unfiltered" : ""
matcher_params += params.n_cpus ? " -n ${params.n_cpus}" : ""

// profile options: local or cluster
// USEAGE: nextflow run template_matcher.nf -profile <profile>

workflow {
    /*
    * Create a channel emitting chunks of file paths from params.input,
    * each line is one file path
    * where each chunk contains batch_size lines written to a temp file
    * the channel 'ch_files' contains these temp files.
    */
    Channel
        .fromPath(params.input)
        .splitText( by: params.batch_size , file:true )
        .set { ch_files }

    /*
     * Execute matching process for each chunk emitted by the 'ch_fasta' channel
     * emit all the tsv files generated
     */
    
    ch_tsv_files = matching(ch_files)

    // merge all tsv files into one
    merged_tsv = ch_tsv_files.collectFile(name: 'merged.tsv', keepHeader: true, skip: 2)

    // move the merged file to output
    merged_tsv.view { f ->
        println "Moving merged TSV to ${params.output}"
        f.moveTo(params.output)
    }
}

process matching {

    input:
    path input_file

    output:
    path "output.tsv"

    script:
    """
    echo $SLURM_JOB_ID
    export PYTHONPATH='/exports/archive/lucid-grpzeller-primary/hackett/template_matching/'
    python -m enzymm -l ${input_file} -o output.tsv ${jess_params}${matcher_params}
    """
}
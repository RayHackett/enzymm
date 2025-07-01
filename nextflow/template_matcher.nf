#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/////// Expected parameters
// params.output
// params.input
// params.batch_size

//////// Optional parameters
// params.template_folder
params.template_folder = params.template_folder ? " -t ${params.template_folder}" : ""

///////// Jess control parameters
if (params.rmsd && params.distance && params.max_dynamic_distance){
    jess_params = " -j ${params.rmsd} ${params.distance} ${params.max_dynamic_distance}"
} else {
    jess_params = ""
}

// params.conservation_cutoff
conservation_cutoff = params.conservation_cutoff ? "-c ${params.conservation_cutoff}" : ""

///////// optional control parameters
// params.warn
warn = params.warn ? " -w" : ""
// params.verbose
verbose = params.verbose ? " -v" : ""
// params.match_smaller_templates
match_small_templates = params.match_small_templates ? " --match-small-templates" : ""
// params.skip_smaller_hits
skip_smaller_hits = params.skip_smaller_hits ? " -skip-smaller-hits" : ""
// params.skip_annotation
skip_annotation = params.skip_annotation ? " --skip-annotation": ""
// params.unfiltered
unfiltered = params.unfiltered ? " --unfiltered" : ""

optional_param_flags = warn + verbose + match_small_templates + skip_smaller_hits + skip_annotation


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

    /*
     * Collect all the tsv files into a single file via the concatendate_tsvs process
     * and print the resulting file contents when complete.
    */
    concatenate_tsvs(ch_tsvs)
}

process matching {

    input:
    path input_file

    output:
    // path("./results/*.pdb")
    path "output.tsv"

    script:
    """
    export PYTHONPATH='/exports/archive/lucid-grpzeller-primary/hackett/template_matching/'
    python -m matcher.jess_run -l ${input_file} -o output.tsv ${params.template_folder}${jess_params}${conservation_cutoff}${optional_param_flags}
    """
}

process concatenate_tsvs {

    input:
    // tsv_files will be a list of paths
    path tsv_files

    output:
    path(params.output)

    // head -n 1 ${tsv_files[0]} > ${params.output} takes the header from the first tsv
    // for f in ${tsv_files[@]} loops over tsv files. @ expands the channel to an array
    // while skipping subsequent header lines with tail -n +2

    script:
    """
    # Use the header from the first file, then skip headers in the rest
    head -n 1 ${tsv_files[0]} > ${params.output}
    for f in ${tsv_files[@]}; do
        tail -n +2 "$f" >> ${params.output}
    done
    echo "Collected template matching results to ${params.output}"
    """
}

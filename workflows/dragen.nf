/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

nextflow.enable.types = true

include { DRAGEN                } from '../modules/local/dragen'
include { DRAGEN_BUILDHASHTABLE } from '../modules/local/dragen_buildhashtable'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

workflow DRAGEN_WORKFLOW {
    take:
    ch_input: Channel<DragenSample>
    index: Path?
    fasta: Path?

    main:
    ch_index_built = channel.empty()
    if (!index) {
        // Built from the samples so that no samples means no build
        ch_index_built = DRAGEN_BUILDHASHTABLE(ch_input.map { _sample -> fasta }.unique())
        // One hash table; reduce makes it a value that every sample reads (null when no samples, so nothing reads it)
        ch_index = ch_index_built.reduce(null) { _acc, index_dir -> index_dir }
    }
    else {
        ch_index = channel.value(index)
    }

    //
    // MODULE: Run DRAGEN on each sample
    //
    ch_results = DRAGEN(ch_input, ch_index)

    emit:
    index       = ch_index
    index_built = ch_index_built // only a hash table built here, not one the user supplied
    results     = ch_results
}

// One sample: every lane's FASTQs, in lane order (fastq_2 is empty for single-end)
record DragenSample {
    id: String
    seq_type: String
    lanes: List<String>
    fastq_1: List<Path>
    fastq_2: List<Path>
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

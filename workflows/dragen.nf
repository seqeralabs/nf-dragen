/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

include { DRAGEN                } from '../modules/local/dragen.nf'
include { DRAGEN_BUILDHASHTABLE } from '../modules/local/dragen_buildhashtable.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

workflow DRAGEN_WORKFLOW {
    take:
    ch_input
    index
    fasta

    main:
    ch_index_built = channel.empty()
    if (!index) {
        ch_fasta = channel.fromPath(fasta, checkIfExists: true, type: 'file')
        // Use combine to not run the process if ch_input is empty
        DRAGEN_BUILDHASHTABLE(
            ch_input.combine(ch_fasta).map { meta, fastq_1, fastq_2, _fasta -> _fasta }.first()
        )
        ch_index = DRAGEN_BUILDHASHTABLE.out.index
        ch_index_built = ch_index
    }
    else {
        ch_index = channel.fromPath(index, checkIfExists: true, type: 'dir')
    }

    //
    // MODULE: Run DRAGEN on DNA samples to generate BAM from FastQ
    //
    DRAGEN(
        ch_input,
        ch_index.collect(),
    )

    emit:
    index        = ch_index
    index_built  = ch_index_built          // only a hash table built here, not one the user supplied
    bam          = DRAGEN.out.bam
    fastq        = DRAGEN.out.fastq
    vcf          = DRAGEN.out.vcf
    tbi          = DRAGEN.out.tbi
    vcf_filtered = DRAGEN.out.vcf_filtered
    tbi_filtered = DRAGEN.out.tbi_filtered
    results      = DRAGEN.out.bam.mix(
        DRAGEN.out.fastq,
        DRAGEN.out.vcf,
        DRAGEN.out.tbi,
        DRAGEN.out.vcf_filtered,
        DRAGEN.out.tbi_filtered,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

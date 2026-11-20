/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DRAGEN                } from '../modules/local/dragen.nf'
include { DRAGEN_BUILDHASHTABLE } from '../modules/local/dragen_buildhashtable.nf'

include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'

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
        ch_versions = Channel.empty()
        if (!index && !fasta) {
            // No reference for this seq_type. Only an error if a sample of this type
            // is actually present, so DNA-only runs need no RNA reference (and vice versa).
            ch_index = ch_input.map { meta, _fastq_1, _fastq_2 ->
                error("No DRAGEN reference for ${meta.seq_type} sample '${meta.id}': set --dragen_index_${meta.seq_type} or --fasta")
            }
        } else if (!index) {
            ch_fasta = Channel.fromPath(fasta, checkIfExists: true, type: 'file')
            // Use combine to not run the process if ch_input is empty
            DRAGEN_BUILDHASHTABLE (
                ch_input.combine(ch_fasta).map { meta, fastq_1, fastq_2, _fasta -> _fasta }.first()
            )
            ch_index    = DRAGEN_BUILDHASHTABLE.out.index
            ch_versions = ch_versions.mix(DRAGEN_BUILDHASHTABLE.out.versions)
        } else {
            ch_index = Channel.fromPath(index, checkIfExists: true, type: 'dir')
        }

        //
        // MODULE: Run DRAGEN on DNA samples to generate BAM from FastQ
        //
        DRAGEN (
            ch_input,
            ch_index.collect()
        )
        ch_versions = ch_versions.mix(DRAGEN.out.versions.first())


    emit:
        index        = ch_index
        bam          = DRAGEN.out.bam
        fastq        = DRAGEN.out.fastq
        vcf          = DRAGEN.out.vcf
        tbi          = DRAGEN.out.tbi
        vcf_filtered  = DRAGEN.out.vcf_filtered
        tbi_filtered  = DRAGEN.out.tbi_filtered
        versions     = DRAGEN.out.versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

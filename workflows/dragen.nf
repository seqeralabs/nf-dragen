/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                } from '../modules/nf-core/fastqc/main'
include { MULTIQC               } from '../modules/nf-core/multiqc/main'
include { DRAGEN                } from '../modules/local/dragen.nf'
include { DRAGEN_BUILDHASHTABLE } from '../modules/local/dragen_buildhashtable.nf'

include { paramsSummaryMap          } from 'plugin/nf-validation'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_dragen_pipeline'

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
        ch_versions     = Channel.empty()
        ch_multiqc_files = Channel.empty()

        if (!index) {
            ch_fasta = Channel.fromPath(fasta, checkIfExists: true, type: 'file')
            DRAGEN_BUILDHASHTABLE (
                ch_input.combine(ch_fasta).map { meta, fastq_1, fastq_2, fasta -> fasta }.first()
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
        index       = ch_index
        bam         = DRAGEN.out.bam
        fastq       = DRAGEN.out.fastq
        vcf         = DRAGEN.out.vcf
        tbi         = DRAGEN.out.tbi
        vcf_filtered = DRAGEN.out.vcf_filtered
        tbi_filtered = DRAGEN.out.tbi_filtered
        versions    = DRAGEN.out.versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

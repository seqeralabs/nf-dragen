/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DRAGEN_BUILDHASHTABLE as DRAGEN_BUILDHASHTABLE_DNA } from '../modules/local/dragen_buildhashtable.nf'
include { DRAGEN_BUILDHASHTABLE as DRAGEN_BUILDHASHTABLE_RNA } from '../modules/local/dragen_buildhashtable.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DRAGEN_INDEX {
    take:
        ch_fasta
        ch_dna_index
        ch_rna_index

    main:
        ch_versions = Channel.empty()

        if (!ch_dna_index) {
            DRAGEN_BUILDHASHTABLE_DNA (
                ch_fasta
            )
            ch_dna_index = DRAGEN_BUILDHASHTABLE_DNA.out.index
            ch_versions = ch_versions.mix(DRAGEN_BUILDHASHTABLE_DNA.out.versions)
        }

        if (!ch_rna_index) {
            DRAGEN_BUILDHASHTABLE_RNA (
                ch_fasta
            )
            ch_rna_index = DRAGEN_BUILDHASHTABLE_RNA.out.index
            ch_versions = ch_versions.mix(DRAGEN_BUILDHASHTABLE_RNA.out.versions)
        }
    
    emit:
        dna     = ch_dna_index
        rna     = ch_rna_index
        version = ch_versions
}

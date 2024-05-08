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
        fasta
        dna_index
        rna_index

    main:
        ch_versions = Channel.empty()

        if (!dna_index) {
            DRAGEN_BUILDHASHTABLE_DNA (
                fasta
            )
            ch_dna_index = DRAGEN_BUILDHASHTABLE_DNA.out.index
            ch_versions = ch_versions.mix(DRAGEN_BUILDHASHTABLE_DNA.out.versions)
        } else {
            ch_dna_index = Channel.fromPath(dna_index, checkIfExists: true, type: 'dir')
        }

        if (!rna_index) {
            DRAGEN_BUILDHASHTABLE_RNA (
                fasta
            )
            ch_rna_index = DRAGEN_BUILDHASHTABLE_RNA.out.index
            ch_versions = ch_versions.mix(DRAGEN_BUILDHASHTABLE_RNA.out.versions)
        } else {
            ch_rna_index = Channel.fromPath(rna_index, checkIfExists: true, type: 'dir')
        }
    
    emit:
        dna     = ch_dna_index
        rna     = ch_rna_index
        version = ch_versions
}

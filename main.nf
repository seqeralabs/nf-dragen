#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    seqeralabs/nf-dragen
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/seqeralabs/nf-dragen
----------------------------------------------------------------------------------------
*/

include { samplesheetToList             } from 'plugin/nf-schema'
include { validateParameters            } from 'plugin/nf-schema'
include { DRAGEN_WORKFLOW as DRAGEN_DNA } from './workflows/dragen'
include { DRAGEN_WORKFLOW as DRAGEN_RNA } from './workflows/dragen'

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow SEQERALABS_DRAGEN {

    take:
        input     // tuple: [ meta, [fastq_1, ...], [fastq_2, ...] ]
        fasta     // str: path to fasta file
        dna_index // str: path to dragen index directory for dna
        rna_index // str: path to dragen index directory for rna

    main:

        DRAGEN_DNA (
            input.filter { meta, fastq_1, fastq_2 -> meta.seq_type == 'dna' && !params.skip_dragen },
            dna_index,
            fasta
        )

        DRAGEN_RNA (
            input.filter { meta, fastq_1, fastq_2 -> meta.seq_type == 'rna' && !params.skip_dragen },
            rna_index,
            fasta
        )

        ch_versions = DRAGEN_DNA.out.versions.mix(DRAGEN_RNA.out.versions)

    emit:
        index          = DRAGEN_DNA.out.index
        bam            = DRAGEN_DNA.out.bam
        fastq          = DRAGEN_DNA.out.fastq
        vcf            = DRAGEN_DNA.out.vcf
        tbi            = DRAGEN_DNA.out.tbi
        vcf_filtered   = DRAGEN_DNA.out.vcf_filtered
        tbi_filtered   = DRAGEN_DNA.out.tbi_filtered
        versions       = ch_versions
}

workflow {

    main:
    validateParameters()

    ch_samplesheet = channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { meta, fastq_1, fastq_2 ->
            return [
                meta + [ seq_type: meta.seq_type ?: params.seq_type],
                fastq_1,
                fastq_2
            ]
        }
        .groupTuple()
        .map { meta, fastq_1, fastq_2 ->
            return [
                meta,
                fastq_1.flatten(),
                fastq_2.flatten()
            ]
        }

    SEQERALABS_DRAGEN (
        ch_samplesheet,
        params.fasta,
        params.dragen_index_dna,
        params.dragen_index_rna
    )

    SEQERALABS_DRAGEN.out.versions
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_versions.yml', sort: true, newLine: true)
}

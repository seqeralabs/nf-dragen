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

params {
    // Input options
    input: Path
    seq_type: String = 'dna'

    // References
    fasta: Path?
    dragen_index_dna: Path?
    dragen_index_rna: Path?

    // Tools
    skip_dragen: Boolean = false

    // Boilerplate options
    outdir: String
    publish_dir_mode: String = 'copy'
    monochrome_logs: Boolean = false
}

workflow {

    main:
    validateParameters()

    ch_samplesheet = channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .groupTuple()
        .map { meta, fastq_1, fastq_2 ->
            return [meta, fastq_1.flatten(), fastq_2.flatten()]
        }

    SEQERALABS_DRAGEN(
        ch_samplesheet,
        params.fasta,
        params.dragen_index_dna,
        params.dragen_index_rna,
        params.seq_type,
        params.skip_dragen,
    )

    // One entry per distinct process/tool version, from every process
    ch_versions = channel.topic('versions')
        .map { versions -> versions.text }
        .unique()
        .collectFile(name: 'software_versions.yml', sort: true, newLine: true)

    publish:
    dragen_index_dna = SEQERALABS_DRAGEN.out.dna_index
    dragen_index_rna = SEQERALABS_DRAGEN.out.rna_index
    dragen_dna       = SEQERALABS_DRAGEN.out.dna
    dragen_rna       = SEQERALABS_DRAGEN.out.rna
    versions         = ch_versions
}

output {
    dragen_index_dna {
        path 'genome/index'
    }
    dragen_index_rna {
        path 'genome/index'
    }
    dragen_dna {
        path 'dragen/dna_fastq_to_vcf'
    }
    dragen_rna {
        path 'dragen/rna_fastq_to_bam'
    }
    versions {
        path 'pipeline_info'
    }
}

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow SEQERALABS_DRAGEN {
    take:
    input       // tuple: [ meta, [fastq_1, ...], [fastq_2, ...] ]
    fasta       // path: fasta file
    dna_index   // path: dragen index directory for dna
    rna_index   // path: dragen index directory for rna
    seq_type    // str: seq_type for samplesheet rows without one
    skip_dragen // bool: skip all DRAGEN steps

    main:

    ch_input = input
        .map { meta, fastq_1, fastq_2 -> [meta + [seq_type: meta.seq_type ?: seq_type], fastq_1, fastq_2] }
        .filter { _meta, _fastq_1, _fastq_2 -> !skip_dragen }

    DRAGEN_DNA(
        ch_input.filter { meta, _fastq_1, _fastq_2 -> meta.seq_type == 'dna' },
        dna_index,
        fasta,
    )

    DRAGEN_RNA(
        ch_input.filter { meta, _fastq_1, _fastq_2 -> meta.seq_type == 'rna' },
        rna_index,
        fasta,
    )

    emit:
    dna_index = DRAGEN_DNA.out.index_built
    rna_index = DRAGEN_RNA.out.index_built
    dna       = DRAGEN_DNA.out.results
    rna       = DRAGEN_RNA.out.results
}

#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seqeralabs/nf-dragen
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Github : https://github.com/seqeralabs/nf-dragen
----------------------------------------------------------------------------------------
 */

include { validateParameters            } from 'plugin/nf-schema'
include { PARSE_SAMPLESHEET             } from './subworkflows/local/parse_samplesheet'
include { SELECT_DRAGEN_MODE            } from './subworkflows/local/select_dragen_mode'
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
    target_bed: Path?

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

    PARSE_SAMPLESHEET(params.input, params.seq_type)

    SEQERALABS_DRAGEN(
        PARSE_SAMPLESHEET.out.samples,
        params.fasta,
        params.dragen_index_dna,
        params.dragen_index_rna,
        params.target_bed,
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
    dragen_somatic   = SEQERALABS_DRAGEN.out.somatic
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
    dragen_somatic {
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
    input       // channel: Sample records from PARSE_SAMPLESHEET
    fasta       // path: fasta file
    dna_index   // path: dragen index directory for dna
    rna_index   // path: dragen index directory for rna
    target_bed  // path: exome targets for germline and somatic runs, or null
    skip_dragen // bool: skip all DRAGEN steps

    main:

    ch_input = input.filter { _sample -> !skip_dragen }

    SELECT_DRAGEN_MODE(ch_input)

    // Germline and somatic runs share one DNA hash table
    DRAGEN_DNA(
        SELECT_DRAGEN_MODE.out.germline,
        SELECT_DRAGEN_MODE.out.somatic,
        dna_index,
        fasta,
        target_bed,
    )

    DRAGEN_RNA(
        SELECT_DRAGEN_MODE.out.rna,
        channel.of(), // not channel.empty(): a typed workflow cannot read that twice (Nextflow 26.10.6)
        rna_index,
        fasta,
        null,
    )

    emit:
    dna_index = DRAGEN_DNA.out.index_built
    rna_index = DRAGEN_RNA.out.index_built
    dna       = DRAGEN_DNA.out.results
    somatic   = DRAGEN_DNA.out.somatic
    rna       = DRAGEN_RNA.out.results
}

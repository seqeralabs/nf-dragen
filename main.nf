#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seqeralabs/nf-dragen
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Github : https://github.com/seqeralabs/nf-dragen
----------------------------------------------------------------------------------------
 */

include { validateParameters                                 } from 'plugin/nf-schema'
include { PARSE_SAMPLESHEET                                  } from './subworkflows/local/parse_samplesheet'
include { SELECT_DRAGEN_MODE                                 } from './subworkflows/local/select_dragen_mode'
include { DRAGEN as DRAGEN_DNA                               } from './modules/local/dragen'
include { DRAGEN as DRAGEN_SOMATIC                           } from './modules/local/dragen'
include { DRAGEN as DRAGEN_RNA                               } from './modules/local/dragen'
include { DRAGEN_BUILDHASHTABLE as DRAGEN_BUILDHASHTABLE_DNA } from './modules/local/dragen_buildhashtable'
include { DRAGEN_BUILDHASHTABLE as DRAGEN_BUILDHASHTABLE_RNA } from './modules/local/dragen_buildhashtable'
include { MULTIQC_REPORT                                     } from './subworkflows/local/multiqc_report'

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
    skip_multiqc: Boolean = false

    // MultiQC options
    multiqc_config: Path?
    multiqc_title: String?
    multiqc_logo: Path?

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

    // One entry per distinct process/tool version, from every process; the Workflow entry keeps the report non-empty when no task runs
    ch_versions_text = channel.topic('versions')
        .map { versions -> versions.text }
        .unique()
        .mix(channel.of("Workflow:\n    Nextflow: ${workflow.nextflow.version}\n    ${workflow.manifest.name}: ${workflow.manifest.version}"))

    ch_multiqc = channel.empty()
    if (!params.skip_multiqc) {
        MULTIQC_REPORT(
            SEQERALABS_DRAGEN.out.dna.mix(SEQERALABS_DRAGEN.out.somatic, SEQERALABS_DRAGEN.out.rna).map { result -> result.metrics.toList() }, // a Set, which collect() would not flatten
            ch_versions_text,
            params.multiqc_config,
            params.multiqc_logo,
        )
        ch_multiqc = MULTIQC_REPORT.out.report
        // MULTIQC cannot report to the versions topic, as its input waits on that topic
        ch_versions_text = ch_versions_text.mix(MULTIQC_REPORT.out.versions)
    }

    ch_versions = ch_versions_text.collectFile(name: 'software_versions.yml', sort: true, newLine: true)

    publish:
    dragen_index_dna = SEQERALABS_DRAGEN.out.dna_index
    dragen_index_rna = SEQERALABS_DRAGEN.out.rna_index
    dragen_dna       = SEQERALABS_DRAGEN.out.dna
    dragen_somatic   = SEQERALABS_DRAGEN.out.somatic
    dragen_rna       = SEQERALABS_DRAGEN.out.rna
    multiqc          = ch_multiqc
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
    multiqc {
        path 'multiqc'
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

    // Hash tables are built only when no index is given, from the samples so that no samples means no build.
    // reduce makes the one table a value that every sample reads (null when no samples, so nothing reads it).
    ch_dna_index_built = channel.empty()
    if (!dna_index) {
        // Germline and somatic runs share one DNA hash table
        ch_dna_index_built = DRAGEN_BUILDHASHTABLE_DNA(
            SELECT_DRAGEN_MODE.out.germline.mix(SELECT_DRAGEN_MODE.out.somatic).map { sample -> referenceFasta(sample, fasta, 'dna') }.unique()
        )
        ch_dna_index = ch_dna_index_built.reduce(null) { _acc, index_dir -> index_dir }
    }
    else {
        ch_dna_index = channel.value(dna_index)
    }

    ch_rna_index_built = channel.empty()
    if (!rna_index) {
        ch_rna_index_built = DRAGEN_BUILDHASHTABLE_RNA(SELECT_DRAGEN_MODE.out.rna.map { sample -> referenceFasta(sample, fasta, 'rna') }.unique())
        ch_rna_index = ch_rna_index_built.reduce(null) { _acc, index_dir -> index_dir }
    }
    else {
        ch_rna_index = channel.value(rna_index)
    }

    ch_dna = DRAGEN_DNA(SELECT_DRAGEN_MODE.out.germline, ch_dna_index, target_bed)
    ch_somatic = DRAGEN_SOMATIC(SELECT_DRAGEN_MODE.out.somatic, ch_dna_index, target_bed)
    // RNA runs get no target_bed
    ch_rna = DRAGEN_RNA(SELECT_DRAGEN_MODE.out.rna, ch_rna_index, null)

    emit:
    dna_index = ch_dna_index_built // only hash tables built here, not ones the user supplied
    rna_index = ch_rna_index_built
    dna       = ch_dna
    somatic   = ch_somatic
    rna       = ch_rna
}

// A sample whose sequencing type has neither --fasta nor a DRAGEN index cannot get a hash table; stop before any task
def referenceFasta(sample, fasta, seq_type) {
    if (!fasta) {
        error("Sample '${sample.id}' (seq_type '${seq_type}') needs a reference: give --fasta to build the DRAGEN hash table, or --dragen_index_${seq_type} to use an existing one")
    }
    return fasta
}

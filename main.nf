#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seqeralabs/nf-dragen
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Github : https://github.com/seqeralabs/nf-dragen
----------------------------------------------------------------------------------------
 */

nextflow.enable.types = true

include { validateParameters                                 } from 'plugin/nf-schema'
include { PARSE_SAMPLESHEET ; Sample } from './subworkflows/local/parse_samplesheet'
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
}

workflow {

    main:
    validateParameters()

    ch_samples = PARSE_SAMPLESHEET(params.input, params.seq_type, file("${projectDir}/assets/schema_input.json"))

    dragen = SEQERALABS_DRAGEN(
        ch_samples,
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
        .mix(channel.of("Workflow:\n    Nextflow: ${nextflow.version}\n    ${workflow.manifest.name}: ${workflow.manifest.version}"))

    ch_multiqc = channel.empty()
    if (!params.skip_multiqc) {
        multiqc_report = MULTIQC_REPORT(
            dragen.dna.mix(dragen.somatic).mix(dragen.rna).flatMap { result -> result.metrics },
            ch_versions_text,
            params.multiqc_config,
            params.multiqc_logo,
            file("${projectDir}/assets/multiqc_config.yml"),
        )
        ch_multiqc = multiqc_report.report
        // MULTIQC cannot report to the versions topic, as its input waits on that topic
        ch_versions_text = ch_versions_text.mix(multiqc_report.versions)
    }

    ch_versions = ch_versions_text.collectFile(name: 'software_versions.yml', sort: true, newLine: true)

    publish:
    dragen_index_dna = dragen.dna_index
    dragen_index_rna = dragen.rna_index
    dragen_dna       = dragen.dna
    dragen_somatic   = dragen.somatic
    dragen_rna       = dragen.rna
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
    input: Channel<Sample>
    fasta: Path?
    dna_index: Path?       // dragen index directory for dna
    rna_index: Path?       // dragen index directory for rna
    target_bed: Path?      // exome targets for germline and somatic runs
    skip_dragen: Boolean   // skip all DRAGEN steps

    main:

    ch_input = input.filter { _sample -> !skip_dragen }

    modes = SELECT_DRAGEN_MODE(ch_input)

    // Hash tables are built only when no index is given, from the samples so that no samples means no build.
    // reduce makes the one table a value that every sample reads (null when no samples, so nothing reads it).
    ch_dna_index_built = channel.empty()
    if (!dna_index) {
        // Germline and somatic runs share one DNA hash table
        ch_dna_fasta = modes.germline
            .map { sample -> sample.id }
            .mix(modes.somatic.map { run -> run.id })
            .map { id ->
                if (!fasta) {
                    error(missingReference(id, 'dna'))
                }
                fasta
            }
        ch_dna_index_built = DRAGEN_BUILDHASHTABLE_DNA(ch_dna_fasta.unique())
        ch_dna_index = ch_dna_index_built.reduce(null) { _acc, index_dir -> index_dir }
    }
    else {
        ch_dna_index = channel.value(dna_index)
    }

    ch_rna_index_built = channel.empty()
    if (!rna_index) {
        ch_rna_fasta = modes.rna.map { sample ->
            if (!fasta) {
                error(missingReference(sample.id, 'rna'))
            }
            fasta
        }
        ch_rna_index_built = DRAGEN_BUILDHASHTABLE_RNA(ch_rna_fasta.unique())
        ch_rna_index = ch_rna_index_built.reduce(null) { _acc, index_dir -> index_dir }
    }
    else {
        ch_rna_index = channel.value(rna_index)
    }

    ch_dna = DRAGEN_DNA(modes.germline, ch_dna_index, target_bed)
    ch_somatic = DRAGEN_SOMATIC(modes.somatic, ch_dna_index, target_bed)
    // RNA runs get no target_bed
    ch_rna = DRAGEN_RNA(modes.rna, ch_rna_index, null)

    emit:
    dna_index: Channel<Path> = ch_dna_index_built // only hash tables built here, not ones the user supplied
    rna_index: Channel<Path> = ch_rna_index_built
    dna: Channel<Record> = ch_dna
    somatic: Channel<Record> = ch_somatic
    rna: Channel<Record> = ch_rna
}

// A sample whose sequencing type has neither --fasta nor a DRAGEN index cannot get a hash table; stop before any task.
// error() is called in the operator closure: under static typing an error() thrown inside a function called from an operator shows only 'InvocationTargetException'
def missingReference(id: String, seq_type: String) -> String {
    return "Sample '${id}' (seq_type '${seq_type}') needs a reference: give --fasta to build the DRAGEN hash table, or --dragen_index_${seq_type} to use an existing one"
}

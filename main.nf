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
    dragen_germline_args: String = '--enable-duplicate-marking true'
    dragen_somatic_args: String = '--enable-duplicate-marking true --vc-target-vaf 0.03 --enable-sv true'
    dragen_rna_args: String = '--enable-duplicate-marking true'
    dragen_hashtable_dna_args: String = '--ht-build-cnv-hashtable true'
    dragen_hashtable_rna_args: String = '--ht-build-rna-hashtable true'
    dragen_license: String = 'marketplace'

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
        params.dragen_germline_args,
        params.dragen_somatic_args,
        params.dragen_rna_args,
        params.dragen_hashtable_dna_args,
        params.dragen_hashtable_rna_args,
        params.dragen_license,
        params.skip_dragen,
    )

    // One entry per distinct process/tool version, from every process; the Workflow entry keeps the report non-empty when no task runs
    ch_versions_text = channel.topic('versions')
        .map { v -> "\"${v.process}\":\n    ${v.tool}: ${v.version}" }
        .unique()
        .mix(channel.of("Workflow:\n    Nextflow: ${nextflow.version}\n    ${workflow.manifest.name}: ${workflow.manifest.version}"))

    ch_multiqc = channel.empty()
    if (!params.skip_multiqc) {
        // MULTIQC reports no version: it reads the versions topic, and its report footer shows its own version
        ch_multiqc = MULTIQC_REPORT(
            dragen.dna.mix(dragen.somatic).mix(dragen.rna).flatMap { result -> result.metrics },
            ch_versions_text,
            params.multiqc_config,
            params.multiqc_logo,
            file("${projectDir}/assets/multiqc_config.yml"),
            params.multiqc_title,
        )
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
    // DRAGEN records: every file in a record is published; the index lists each record's sample fields and published paths
    dragen_dna: Channel<Record> {
        path 'dragen/germline/'
        index {
            path 'dragen/germline/germline_index.csv'
            header true
        }
    }
    dragen_somatic: Channel<Record> {
        path 'dragen/somatic/'
        index {
            path 'dragen/somatic/somatic_index.csv'
            header true
        }
    }
    dragen_rna: Channel<Record> {
        path 'dragen/rna/'
        index {
            path 'dragen/rna/rna_index.csv'
            header true
        }
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
    dna_index: Path?                  // dragen index directory for dna
    rna_index: Path?                  // dragen index directory for rna
    target_bed: Path?                 // exome targets for germline and somatic runs
    dragen_germline_args: String      // tunable DRAGEN options for germline runs
    dragen_somatic_args: String       // tunable DRAGEN options for somatic runs
    dragen_rna_args: String           // tunable DRAGEN options for RNA runs
    dragen_hashtable_dna_args: String // DRAGEN options for the DNA hash table build
    dragen_hashtable_rna_args: String // DRAGEN options for the RNA hash table build
    dragen_license: String            // marketplace or byol
    skip_dragen: Boolean              // skip all DRAGEN steps

    main:

    ch_input = input.filter { _sample -> !skip_dragen }

    modes = SELECT_DRAGEN_MODE(ch_input)

    // Single quotes: the shell, not Nextflow, expands the secrets
    def lic_server = dragen_license == 'byol' ? '--lic-server=$DRAGEN_USERNAME:$DRAGEN_PASSWORD@license.edicogenome.com' : null

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
        ch_dna_index_built = DRAGEN_BUILDHASHTABLE_DNA(
            ch_dna_fasta.unique().map { ref -> record(id: 'dna', fasta: ref, dragen_args: dragen_hashtable_dna_args, dragen_lic_server: lic_server, output_prefix: 'dragen_index_dna') }
        ).map { table -> table.index }
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
        ch_rna_index_built = DRAGEN_BUILDHASHTABLE_RNA(
            ch_rna_fasta.unique().map { ref -> record(id: 'rna', fasta: ref, dragen_args: dragen_hashtable_rna_args, dragen_lic_server: lic_server, output_prefix: 'dragen_index_rna') }
        ).map { table -> table.index }
        ch_rna_index = ch_rna_index_built.reduce(null) { _acc, index_dir -> index_dir }
    }
    else {
        ch_rna_index = channel.value(rna_index)
    }

    def map_align = '--enable-map-align true --enable-map-align-output true --enable-bam-indexing true --enable-sort=true'
    def dna_args = "${map_align} --enable-variant-caller true --output-format=BAM ${dragen_germline_args}"
    def somatic_args = "${map_align} --enable-variant-caller true --output-format=BAM ${dragen_somatic_args}"
    def rna_args = "--enable-rna true ${map_align} --output-format=BAM ${dragen_rna_args}"

    ch_dna = DRAGEN_DNA(
        modes.germline.combine(index: ch_dna_index).map { run -> run + record(target_bed: target_bed, dragen_args: dna_args, dragen_lic_server: lic_server) }
    )
    ch_somatic = DRAGEN_SOMATIC(
        modes.somatic.combine(index: ch_dna_index).map { run -> run + record(target_bed: target_bed, dragen_args: somatic_args, dragen_lic_server: lic_server) }
    )
    // RNA runs get no target_bed
    ch_rna = DRAGEN_RNA(
        modes.rna.combine(index: ch_rna_index).map { run -> run + record(target_bed: null, dragen_args: rna_args, dragen_lic_server: lic_server) }
    )

    // DRAGEN results with the samples they came from, for the published index files
    ch_dna_published = modes.germline
        .map { sample -> record(id: sample.id, patient: sample.patient, sample: sample.id) }
        .join(ch_dna, by: 'id')
    ch_somatic_published = modes.somatic
        .map { run -> record(id: "${run.tumour_id}_vs_${run.id}", patient: run.patient, tumour: run.tumour_id, normal: run.id) }
        .join(ch_somatic, by: 'id')
    ch_rna_published = modes.rna
        .map { sample -> record(id: sample.id, patient: sample.patient, sample: sample.id) }
        .join(ch_rna, by: 'id')

    emit:
    dna_index: Channel<Path> = ch_dna_index_built   // only hash tables built here, not ones the user supplied
    rna_index: Channel<Path> = ch_rna_index_built
    dna: Channel<Record> = ch_dna_published
    somatic: Channel<Record> = ch_somatic_published
    rna: Channel<Record> = ch_rna_published
}

// A sample whose sequencing type has neither --fasta nor a DRAGEN index cannot get a hash table; stop before any task.
// error() is called in the operator closure: under static typing an error() thrown inside a function called from an operator shows only 'InvocationTargetException'
def missingReference(id: String, seq_type: String) -> String {
    return "Sample '${id}' (seq_type '${seq_type}') needs a reference: give --fasta to build the DRAGEN hash table, or --dragen_index_${seq_type} to use an existing one"
}

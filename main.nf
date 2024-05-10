#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    seqeralabs/nf-dragen
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/seqeralabs/nf-dragen
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap              } from 'plugin/nf-validation'
include { getGenomeAttribute            } from './subworkflows/local/utils_nfcore_dragen_pipeline'
include { softwareVersionsToYAML        } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { paramsSummaryMultiqc          } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText        } from './subworkflows/local/utils_nfcore_dragen_pipeline'

include { PIPELINE_INITIALISATION       } from './subworkflows/local/utils_nfcore_dragen_pipeline'
include { FASTQC                        } from './modules/nf-core/fastqc/main'
include { DRAGEN_WORKFLOW as DRAGEN_DNA } from './workflows/dragen'
include { DRAGEN_WORKFLOW as DRAGEN_RNA } from './workflows/dragen'
include { MULTIQC                       } from './modules/nf-core/multiqc/main'
include { PIPELINE_COMPLETION           } from './subworkflows/local/utils_nfcore_dragen_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta = getGenomeAttribute('fasta')


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main nf-core/dragen analysis pipeline
//
workflow SEQERALABS_DRAGEN {
    take:
        input
        fasta
        dna_index
        rna_index

    main:

        ch_versions     = Channel.empty()
        ch_multiqc_files = Channel.empty()

        //
        // FASTQC
        //
        input
            .filter { params.fastqc }
            // Coerce input to format for FASTQC input
            .map { meta, fastq_1, fastq_2 -> [meta, fastq_1 + fastq_2]}
            | FASTQC
        ch_versions     = ch_versions.mix(FASTQC.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { it[1] })

        // Branch DNA and RNA samples
        input
            .branch { meta, fastq_1, fastq_2 ->
                dna: meta.seq_type == 'dna'
                rna: meta.seq_type == 'rna'
            }
            .set { input_branched}

        //
        // WORKFLOW: Run DNA pipeline
        //
        DRAGEN_DNA (
            input_branched.dna,
            dna_index,
            fasta
        )
        ch_versions = ch_versions.mix(DRAGEN_DNA.out.versions)

        //
        // WORKFLOW: Run RNA pipeline
        //
        DRAGEN_RNA (
            input_branched.rna,
            rna_index,
            fasta
        )
        ch_versions = ch_versions.mix(DRAGEN_RNA.out.versions)

        if ( params.multiqc ) {
            //
            // MODULE: MultiQC
            //
            softwareVersionsToYAML(ch_versions)
                .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
                .set { ch_collated_versions }

            // Build custom multiqc content for the repot
            ch_multiqc_config                      = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
            ch_multiqc_custom_config               = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
            ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
            summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
            ch_workflow_summary                    = Channel.value(paramsSummaryMultiqc(summary_params))
            ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
            ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
            ch_multiqc_files                       = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
            ch_multiqc_files                       = ch_multiqc_files.mix(ch_collated_versions)
            ch_multiqc_files                       = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

            MULTIQC (
                ch_multiqc_files.collect(),
                ch_multiqc_config.toList(),
                ch_multiqc_custom_config.toList(),
                ch_multiqc_logo.toList()
            )
            multiqc_report = MULTIQC.out.report.toList()
            ch_versions    = ch_versions.mix(MULTIQC.out.versions)
        } else {
            multiqc_report = Channel.empty()
        }

    emit:
        index          = DRAGEN_DNA.out.index
        bam            = DRAGEN_DNA.out.bam
        fastq          = DRAGEN_DNA.out.fastq
        vcf            = DRAGEN_DNA.out.vcf
        tbi            = DRAGEN_DNA.out.tbi
        vcf_filtered    = DRAGEN_DNA.out.vcf_filtered
        tbi_filtered    = DRAGEN_DNA.out.tbi_filtered
        multiqc_report = multiqc_report
        versions       = DRAGEN_DNA.out.versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    SEQERALABS_DRAGEN (
        PIPELINE_INITIALISATION.out.samplesheet,
        params.fasta,
        params.dragen_index_dna,
        params.dragen_index_rna
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        SEQERALABS_DRAGEN.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

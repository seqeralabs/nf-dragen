include { MULTIQC } from '../../../modules/local/multiqc/main'

workflow MULTIQC_REPORT {
    take:
    metrics        // channel: DRAGEN metrics CSVs
    versions       // channel: versions.yml text, one entry per process
    multiqc_config // path: extra MultiQC config, or null
    multiqc_logo   // path: report logo, or null

    main:
    // MultiQC reads *_mqc_versions.yml files into its software versions section
    ch_versions_mqc = versions.collectFile(name: 'nf_dragen_software_mqc_versions.yml', sort: true, newLine: true)
    ch_config = [file("${projectDir}/assets/multiqc_config.yml"), multiqc_config].findAll { config -> config }

    MULTIQC(
        metrics.mix(ch_versions_mqc).collect().map { files -> record(id: 'multiqc', files: files.toSet(), configs: ch_config.toSet(), logo: multiqc_logo) }
    )

    emit:
    report   = MULTIQC.out.flatMap { run -> [run.report, run.data, run.plots].findAll { f -> f } } // HTML, multiqc_data/ and plots
    versions = MULTIQC.out.map { run -> "\"${run.process}\":\n    multiqc: ${run.version}" }
}

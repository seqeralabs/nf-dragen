include { MULTIQC } from '../../../modules/nf-core/multiqc/main'

workflow MULTIQC_REPORT {
    take:
    metrics        // channel: DRAGEN metrics CSVs
    versions       // channel: versions.yml text, one entry per process
    multiqc_config // path: extra MultiQC config, or null
    multiqc_logo   // path: report logo, or null

    main:
    // MultiQC reads *_mqc_versions.yml files into its software versions section
    ch_versions_mqc = versions.collectFile(name: 'nf_dragen_software_mqc_versions.yml', sort: true, newLine: true)
    ch_config = [file("${moduleDir}/../../../assets/multiqc_config.yml"), multiqc_config].findAll { config -> config }

    MULTIQC(
        metrics.mix(ch_versions_mqc).collect().map { files -> [[id: 'multiqc'], files, ch_config, multiqc_logo ?: [], [], []] }
    )

    emit:
    report   = MULTIQC.out.report.mix(MULTIQC.out.data, MULTIQC.out.plots).map { _meta, report_file -> report_file } // HTML, multiqc_data/ and plots
    versions = MULTIQC.out.versions.map { process, tool, version -> "\"${process}\":\n    ${tool}: ${version}" }
}

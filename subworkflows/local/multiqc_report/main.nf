nextflow.enable.types = true

include { MULTIQC } from '../../../modules/local/multiqc/main'

workflow MULTIQC_REPORT {
    take:
    metrics: Channel<Path>    // DRAGEN metrics CSVs
    versions: Channel<String> // software versions YAML text, one entry per process
    multiqc_config: Path?     // extra MultiQC config
    multiqc_logo: Path?       // report logo
    bundled_config: Path      // the pipeline's MultiQC config, overridden by multiqc_config
    multiqc_title: String?    // report title

    main:
    // MultiQC reads *_mqc_versions.yml files into its software versions section
    ch_versions_mqc = versions.collectFile(name: 'nf_dragen_software_mqc_versions.yml', sort: true, newLine: true)
    ch_config = multiqc_config ? [bundled_config, multiqc_config] : [bundled_config]
    multiqc_args = multiqc_title ? "--title \"${multiqc_title}\"" : null

    ch_multiqc = MULTIQC(
        metrics.mix(ch_versions_mqc).collect().map { files -> record(id: 'multiqc', files: files.toSet(), configs: ch_config, logo: multiqc_logo, multiqc_args: multiqc_args) }
    )

    emit:
    report: Channel<Path> = ch_multiqc.flatMap { run -> [run.report, run.data, run.plots].findAll { f -> f } } // HTML, multiqc_data/ and plots
}

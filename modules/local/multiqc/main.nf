nextflow.enable.types = true

process MULTIQC {
    tag "${id}"
    label 'process_single'

    container 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'

    input:
    record(
        id: String,
        files: Set<Path>,
        configs: List<Path>,
        logo: Path?,
        multiqc_args: String?
    )

    stage:
    stageAs files, '?/*'
    stageAs configs, '?/*'

    output:
    record(
        id: id,
        report: file("*.html"),
        data: file("*_data"),
        plots: file("*_plots", optional: true),
    )

    script:
    def args = multiqc_args ?: ''
    def config = configs.collect { config_file -> "--config ${config_file}" }.join(' ')
    def logo_arg = logo ? "--cl-config 'custom_logo: \"${logo}\"'" : ''
    """
    multiqc \\
        --force \\
        ${args} \\
        ${config} \\
        ${logo_arg} \\
        .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_data/.stub
    mkdir multiqc_plots
    touch multiqc_plots/.stub
    touch multiqc_report.html
    """
}

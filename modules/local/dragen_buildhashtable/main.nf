nextflow.enable.types = true

process DRAGEN_BUILDHASHTABLE {
    tag "${fasta.name}"
    label 'dragen'

    input:
    record(
        id: String,
        fasta: Path,
        dragen_args: String?,
        dragen_lic_server: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        index: file(prefix),
    )

    topic:
    record(process: "${task.process}", tool: 'dragen', version: eval(version_cmd)) >> 'versions'

    script:
    def args = dragen_args ?: ''
    version_cmd = 'echo $(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g"'
    prefix = output_prefix ?: 'dragen'
    """
    mkdir -p ${prefix}

    /opt/edico/bin/dragen \\
        --build-hash-table true \\
        --output-directory ${prefix} \\
        --ht-reference ${fasta} \\
        ${dragen_lic_server ?: ''} \\
        ${args}
    """

    stub:
    def args = dragen_args ?: ''
    // the stub never calls the binary, which exists only on DRAGEN hosts
    version_cmd = 'echo stub'
    prefix = output_prefix ?: 'dragen'
    // Quoted, so the echo prints the variable names, not the credentials
    def lic_server = dragen_lic_server ? "'${dragen_lic_server}'" : ''
    """
    mkdir -p ${prefix}

    echo /opt/edico/bin/dragen \\
        --build-hash-table true \\
        --output-directory ${prefix} \\
        --ht-reference ${fasta} \\
        ${lic_server} \\
        ${args}

    touch ${prefix}/dragen.ht
    """
}

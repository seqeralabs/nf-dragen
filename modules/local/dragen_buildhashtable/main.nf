nextflow.enable.types = true

process DRAGEN_BUILDHASHTABLE {
    tag "${fasta.name}"
    label 'dragen'

    input:
    record(
        id: String,
        fasta: Path,
        dragen_args: String?,
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
    // -profile byol provides the DRAGEN_USERNAME/DRAGEN_PASSWORD secrets as env vars; the script holds only their names
    """
    mkdir -p ${prefix}

    /opt/edico/bin/dragen \\
        --build-hash-table true \\
        --output-directory ${prefix} \\
        --ht-reference ${fasta} \\
        \${DRAGEN_USERNAME:+--lic-server=\$DRAGEN_USERNAME:\$DRAGEN_PASSWORD@license.edicogenome.com} \\
        ${args}
    """

    stub:
    def args = dragen_args ?: ''
    // the stub never calls the binary, which exists only on DRAGEN hosts
    version_cmd = 'echo stub'
    prefix = output_prefix ?: 'dragen'
    // single-quoted inside the expansion, so the echo prints the secret names, not the credentials
    """
    mkdir -p ${prefix}

    echo /opt/edico/bin/dragen \\
        --build-hash-table true \\
        --output-directory ${prefix} \\
        --ht-reference ${fasta} \\
        \${DRAGEN_USERNAME:+'--lic-server=\$DRAGEN_USERNAME:\$DRAGEN_PASSWORD@license.edicogenome.com'} \\
        ${args}

    touch ${prefix}/dragen.ht
    """
}

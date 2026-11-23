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
    file('versions.yml') >> 'versions'

    script:
    def args = dragen_args ?: ''
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: '\$(echo \$(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g")'
    END_VERSIONS
    """

    stub:
    def args = dragen_args ?: ''
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: stub
    END_VERSIONS
    """
}

nextflow.enable.types = true

process DRAGEN_BUILDHASHTABLE {
    tag "${fasta.name}"
    label 'dragen'

    secret 'DRAGEN_USERNAME'
    secret 'DRAGEN_PASSWORD'

    input:
    fasta: Path

    output:
    file(prefix)

    topic:
    file('versions.yml') >> 'versions'

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'dragen'
    """
    mkdir -p ${prefix}

    /opt/edico/bin/dragen \\
        --build-hash-table true \\
        --output-directory ${prefix} \\
        --ht-reference ${fasta} \\
        --lic-server=\$DRAGEN_USERNAME:\$DRAGEN_PASSWORD@license.edicogenome.com \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: '\$(echo \$(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g")'
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'dragen'
    """
    mkdir -p ${prefix}

    echo /opt/edico/bin/dragen \\
        --build-hash-table true \\
        --output-directory ${prefix} \\
        --ht-reference ${fasta} \\
        --lic-server=\$DRAGEN_USERNAME:\$DRAGEN_PASSWORD@license.edicogenome.com \\
        ${args}

    touch ${prefix}/dragen.ht

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: stub
    END_VERSIONS
    """
}

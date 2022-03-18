process DRAGEN_BUILDHASHTABLE {
    tag "$fasta"
    label 'process_medium'
    label 'dragen'

    secret 'DRAGEN_USERNAME'
    secret 'DRAGEN_PASSWORD'

    input:
    path fasta

    output:
    path "$prefix"     , emit: index
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'dragen'
    """
    mkdir -p $prefix

    /opt/edico/bin/dragen_reset

    /opt/edico/bin/dragen \\
        --build-hash-table true \\
        --output-directory $prefix \\
        --ht-reference $fasta \\
        --lic-server=$DRAGEN_USERNAME:$DRAGEN_PASSWORD@license.edicogenome.com \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo \$(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g")
    END_VERSIONS
    """
}

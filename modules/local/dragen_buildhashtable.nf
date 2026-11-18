process DRAGEN_BUILDHASHTABLE {
    tag "$fasta"
    label 'dragen'

    input:
    path(fasta)

    output:
    path "$prefix"     , emit: index
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    // BYOL needs the license server credentials (supplied as pipeline secrets, see conf/byol.config).
    // Marketplace AMIs are metered by AWS and take no license arguments.
    def lic = params.dragen_license == 'byol' ? "--lic-server=\$DRAGEN_USERNAME:\$DRAGEN_PASSWORD@license.edicogenome.com" : ''
    prefix = task.ext.prefix ?: 'dragen'
    """
    mkdir -p $prefix

    /opt/edico/bin/dragen \\
        --build-hash-table true \\
        --output-directory $prefix \\
        --ht-reference $fasta \\
        $lic \\
        $args \\
        $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: '\$(echo \$(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g")'
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    // BYOL needs the license server credentials (supplied as pipeline secrets, see conf/byol.config).
    // Marketplace AMIs are metered by AWS and take no license arguments.
    def lic = params.dragen_license == 'byol' ? "--lic-server=\$DRAGEN_USERNAME:\$DRAGEN_PASSWORD@license.edicogenome.com" : ''
    prefix = task.ext.prefix ?: 'dragen'
    """
    mkdir -p $prefix

    echo /opt/edico/bin/dragen \\
        --build-hash-table true \\
        --output-directory $prefix \\
        --ht-reference $fasta \\
        $lic \\
        $args \\
        $args2

    touch ${prefix}/dragen.ht
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: '\$(echo \$(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g")'
    END_VERSIONS
    """
}

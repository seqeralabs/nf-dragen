process DRAGEN_BUILDHASHTABLE {
    tag "$fasta"
    label 'dragen'

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

    /opt/edico/bin/dragen \\
        --build-hash-table true \\
        --output-directory $prefix \\
        --ht-reference $fasta \\
        --lic-server=\$DRAGEN_USERNAME:\$DRAGEN_PASSWORD@license.edicogenome.com \\
        $args

    /opt/edico/bin/dragen --version > version.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(cat version.txt | sed -e "s/dragen Version //g")
    END_VERSIONS
    """
}

include { executeOnJarvice } from './executeOnJarvice.nf'

process DRAGEN_BUILDHASHTABLE {
    tag "$fasta"
    label 'dragen'

    cpus 1
    memory 1.GB
    container 'us-docker.pkg.dev/jarvice/images/jarvice-dragen-service:1.2-seqera-gcsfuse'

    input:
    path(fasta)

    output:
    path "$prefix"     , emit: index
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'dragen'
    def dragen_cmd = """
    mkdir -p $prefix

    /opt/edico/bin/dragen \\
        --build-hash-table true \\
        --output-directory $prefix \\
        --ht-reference $fasta \\
        --lic-server=\$DRAGEN_USERNAME:\$DRAGEN_PASSWORD@license.edicogenome.com \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: '\$(echo \$(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g")'
    END_VERSIONS
    """
    executeOnJarvice(dragen_cmd)
}

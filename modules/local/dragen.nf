process DRAGEN {
    tag "$meta.id"
    label 'dragen'

    secret 'DRAGEN_USERNAME'
    secret 'DRAGEN_PASSWORD'

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)
    path index

    output:
    tuple val(meta), path('*.bam')                             , emit: bam        , optional:true
    tuple val(meta), path('*fastq.gz')                         , emit: fastq      , optional:true
    tuple val(meta), path("${prefix}.vcf.gz")                   , emit: vcf        , optional:true
    tuple val(meta), path("${prefix}.vcf.gz.tbi")               , emit: tbi        , optional:true
    tuple val(meta), path("${prefix}.hard-filtered.vcf.gz")      , emit: vcf_filtered, optional:true
    tuple val(meta), path("${prefix}.hard-filtered.vcf.gz.tbi")  , emit: tbi_filtered, optional:true
    path  "versions.yml"                                       , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def ref = index ? "-r $index" : ''

    // Generate appropriate parameter for input files
    def input = "-1 ${fastq_1} -2 ${fastq_2}"
    def rgid  = meta.rgid ? "--RGID ${meta.rgid}" : "--RGID ${meta.id}"
    def rgsm  = meta.rgsm ? "--RGSM ${meta.rgsm}" : "--RGSM ${meta.id}"
    """
    /opt/edico/bin/dragen \\
        $ref \\
        --output-directory ./ \\
        --output-file-prefix $prefix \\
        --lic-server=\$DRAGEN_USERNAME:\$DRAGEN_PASSWORD@license.edicogenome.com \\
        $input \\
        $rgid \\
        $rgsm \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo \$(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def ref = index ? "-r $index" : ''

    // Generate appropriate parameter for input files
    def input = "-1 ${fastq_1} -2 ${fastq_2}"
    def rgid  = meta.rgid ? "--RGID ${meta.rgid}" : "--RGID ${meta.id}"
    def rgsm  = meta.rgsm ? "--RGSM ${meta.rgsm}" : "--RGSM ${meta.id}"
    """
    echo /opt/edico/bin/dragen \\
        $ref \\
        --output-directory ./ \\
        --output-file-prefix $prefix \\
        --lic-server=\$DRAGEN_USERNAME:\$DRAGEN_PASSWORD@license.edicogenome.com \\
        $input \\
        $rgid \\
        $rgsm \\
        $args
    touch ${prefix}.bam
    touch ${prefix}.1.fastq.gz
    touch ${prefix}.2.fastq.gz
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.hard-filtered.vcf.gz
    touch ${prefix}.hard-filtered.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: v1
    END_VERSIONS
    """
}

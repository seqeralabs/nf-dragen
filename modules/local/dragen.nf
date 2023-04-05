process DRAGEN {
    tag "$meta.id"
    label 'dragen'

    input:
    tuple val(meta), path(files_in)
    path index

    output:
    tuple val(meta), path('*.bam')                             , emit: bam         , optional:true
    tuple val(meta), path('*fastq.gz')                         , emit: fastq       , optional:true
    tuple val(meta), path("${prefix}.vcf.gz")                  , emit: vcf         , optional:true
    tuple val(meta), path("${prefix}.vcf.gz.tbi")              , emit: tbi         , optional:true
    tuple val(meta), path("${prefix}.hard-filtered.vcf.gz")    , emit: vcf_filtered, optional:true
    tuple val(meta), path("${prefix}.hard-filtered.vcf.gz.tbi"), emit: tbi_filtered, optional:true
    path  "versions.yml"                                       , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def ref = index ? "-r $index" : ''

    // Generate appropriate parameter for input files
    def input = ''
    def rgid = ''
    def rgdm = ''
    def file_list = files_in.collect { it.toString() }
    if (file_list[0].endsWith('.bam')) {
        input = "-b ${files_in}"
    } else {
        input = meta.single_end ? "-1 ${files_in}" : "-1 ${files_in[0]} -2 ${files_in[1]}"
        rgid = meta.rgid ? "--RGID ${meta.rgid}" : "--RGID ${meta.id}"
        rgsm = meta.rgsm ? "--RGSM ${meta.rgsm}" : "--RGSM ${meta.id}"
    }
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
}

process DRAGEN {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(files_in)
    path index

    output:
    tuple val(meta), path('*.bam')    , emit: bam  , optional:true
    tuple val(meta), path('*.vcf.gz') , emit: vcf  , optional:true
    tuple val(meta), path('*.tbi')    , emit: tbi  , optional:true
    tuple val(meta), path('*fastq.gz'), emit: fastq, optional:true
    path  "versions.yml"              , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Generate appropriate parameter for input files
    def file_list = files_in.collect { it.toString() }
    def input = ''
    if (file_list[0].endsWith('.bam')) {
        input = "-b $files_in"
    } else {
        input = meta.single_end ? "-1 $files_in" : "-1 ${files_in[0]} -2 ${files_in[1]}"
    }
    def ref = index ? "-r $index" : ''
    def rgid = meta.rgid ? "--RGID ${meta.id}" : ''
    def rgsm = meta.rgsm ? "--RGSM ${meta.id}" : ''
    """
    dragen \\
        $ref \\
        --output-directory ./ \\
        --output-file-prefix $prefix \\
        $input \\
        $rgid \\
        $rgsm \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo \$(dragen --version))
    END_VERSIONS
    """
}

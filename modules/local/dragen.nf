process DRAGEN {
    tag "$meta.id"
    label 'dragen'

    secret 'DRAGEN_USERNAME'
    secret 'DRAGEN_PASSWORD'

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
    tuple val(meta), path("${prefix}.*.csv")                   , emit: csv , optional:true
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

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Generate stub files
    // Include one possible csv file, and faked data to prompt multiqc
    """
    touch ${prefix}.bam
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.hard-filtered.vcf.gz
    touch ${prefix}.hard-filtered.vcf.gz.tbi

    # Define the header row of the CSV file
    echo "RUN TIME,,Time loading reference,00:00.1,0.09" > ${prefix}.time_metrics.csv
    echo "RUN TIME,,Time aligning reads,00:00.5,0.53" >> ${prefix}.time_metrics.csv
    echo "RUN TIME,,Time sorting,00:00.4,0.38" >> ${prefix}.time_metrics.csv
    echo "RUN TIME,,Time DRAGStr calibration,00:00.0,0" >> ${prefix}.time_metrics.csv
    echo "RUN TIME,,Time saving map/align output,00:05.6,5.57" >> ${prefix}.time_metrics.csv
    echo "RUN TIME,,Time variant calling,00:05.5,5.45" >> ${prefix}.time_metrics.csv
    echo "RUN TIME,,Time partitioning,00:00.1,0.12" >> ${prefix}.time_metrics.csv
    echo "RUN TIME,,Time structural variant calling,00:02.8,2.83" >> ${prefix}.time_metrics.csv
    echo "RUN TIME,,Time accessing license server,00:13.2,13.21" >> ${prefix}.time_metrics.csv
    echo "RUN TIME,,Total runtime,00:27.2,27.17" >> ${prefix}.time_metrics.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo "vSTUB")
    END_VERSIONS
    """
}

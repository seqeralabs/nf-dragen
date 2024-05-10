process DRAGEN {
    tag "$meta.id"
    label 'dragen'

    secret 'DRAGEN_USERNAME'
    secret 'DRAGEN_PASSWORD'

    input:
    tuple val(meta), path(fastq_1, stageAs: "input_S1_L001_R1_00*.fastq.gz"), path(fastq_2, stageAs: "input_S1_L001_R2_00*.fastq.gz")
    path index

    output:
    tuple val(meta), path('*.{bam,sam,cram}')                   , emit: bam        , optional:true
    tuple val(meta), path('*fastq.gz')                          , emit: fastq      , optional:true
    tuple val(meta), path("${prefix}.vcf.gz")                   , emit: vcf        , optional:true
    tuple val(meta), path("${prefix}.vcf.gz.tbi")               , emit: tbi        , optional:true
    tuple val(meta), path("${prefix}.hard-filtered.vcf.gz")      , emit: vcf_filtered, optional:true
    tuple val(meta), path("${prefix}.hard-filtered.vcf.gz.tbi")  , emit: tbi_filtered, optional:true
    path  "versions.yml"                                       , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def ref = index ? "-r $index" : ''


    // Check FASTQ numbers match
    def num_r1 = fastq_1 instanceof List ? fastq_1.size() : 1
    def num_r2 = fastq_2 instanceof List ? fastq_2.size() : 1

    if ( fastq_2 && num_r1 != num_r2 ) {
        error "Number of R1 and R2 FASTQ files do not match"
    }

    // Generate appropriate parameter for input files
    r1_in = fastq_1 ? "-1 input_S1_L001_R1_001.fastq.gz" : "" 
    r2_in = fastq_2 ? "-2 input_S1_L001_R2_001.fastq.gz" : ""
    def rgid  = meta.rgid ? "--RGID ${meta.rgid}" : "--RGID ${meta.id}"
    def rgsm  = meta.rgsm ? "--RGSM ${meta.rgsm}" : "--RGSM ${meta.id}"
    """
    /opt/edico/bin/dragen \\
        $ref \\
        --output-directory ./ \\
        --output-file-prefix $prefix \\
        --lic-server=\$DRAGEN_USERNAME:\$DRAGEN_PASSWORD@license.edicogenome.com \\
        $r1_in \\
        $r2_in \\
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


    // Check FASTQ numbers match
    def num_r1 = fastq_1 instanceof List ? fastq_1.size() : 1
    def num_r2 = fastq_2 instanceof List ? fastq_2.size() : 1

    if ( fastq_2 && num_r1 != num_r2 ) {
        error "Number of R1 and R2 FASTQ files do not match for sample ${prefix}"
    }

    // Generate appropriate parameter for input files
    r1_in = fastq_1 ? "-1 input_S1_L001_R1_001.fastq.gz" : "" 
    r2_in = fastq_2 ? "-2 input_S1_L001_R2_001.fastq.gz" : ""
    def rgid  = meta.rgid ? "--RGID ${meta.rgid}" : "--RGID ${meta.id}"
    def rgsm  = meta.rgsm ? "--RGSM ${meta.rgsm}" : "--RGSM ${meta.id}"
    """
    echo /opt/edico/bin/dragen \\
        $ref \\
        --output-directory ./ \\
        --output-file-prefix $prefix \\
        --lic-server=\$DRAGEN_USERNAME:\$DRAGEN_PASSWORD@license.edicogenome.com \\
        $r1_in \\
        $r2_in \\
        $rgid \\
        $rgsm \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo \$(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g")
    END_VERSIONS
    """
}

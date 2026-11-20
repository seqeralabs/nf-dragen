process DRAGEN {
    tag "$meta.id"
    label 'dragen'

    input:
    tuple val(meta), path(fastq_1, stageAs: "input_S1_L001_R1_???.fastq.gz"), path(fastq_2, stageAs: "input_S1_L001_R2_???.fastq.gz")
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
    // BYOL needs the license server credentials (supplied as pipeline secrets, see conf/byol.config).
    // Marketplace AMIs are metered by AWS and take no license arguments.
    def lic = params.dragen_license == 'byol' ? "--lic-server=\$DRAGEN_USERNAME:\$DRAGEN_PASSWORD@license.edicogenome.com" : ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def ref = index ? "-r $index" : ''


    // Check FASTQ numbers match
    def num_r1 = fastq_1 instanceof List ? fastq_1.size() : 1
    def num_r2 = fastq_2 instanceof List ? fastq_2.size() : 1

    if ( fastq_2 && num_r1 != num_r2 ) {
        error "Number of R1 and R2 FASTQ files do not match"
    }

    // DRAGEN only reads the FASTQ files named on its command line, so every staged pair goes into a
    // fastq-list CSV (one row per pair). With more than one pair each row gets its own read group ID.
    def rgid = meta.rgid ?: meta.id
    def rgsm = meta.rgsm ?: meta.id
    def rglb = meta.rglb ?: 'UnknownLibrary'
    def row_rgid = num_r1 > 1 ? "${rgid}.\${lane}" : "${rgid}"
    def inputs = "--fastq-list fastq_list.csv --fastq-list-sample-id ${rgsm}"
    """
    echo "RGID,RGSM,RGLB,Lane,Read1File,Read2File" > fastq_list.csv
    lane=0
    for r1 in input_S1_L001_R1_*.fastq.gz; do
        lane=\$((lane + 1))
        r2=\${r1/_R1_/_R2_}
        [ -e "\$r2" ] && r2="\$PWD/\$r2" || r2=""
        echo "${row_rgid},${rgsm},${rglb},\${lane},\$PWD/\${r1},\${r2}" >> fastq_list.csv
    done
    cat fastq_list.csv
    /opt/edico/bin/dragen \\
        $ref \\
        --output-directory ./ \\
        --output-file-prefix $prefix \\
        $lic \\
        $inputs \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: '\$(echo \$(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g")'
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    // BYOL needs the license server credentials (supplied as pipeline secrets, see conf/byol.config).
    // Marketplace AMIs are metered by AWS and take no license arguments.
    def lic = params.dragen_license == 'byol' ? "--lic-server=\$DRAGEN_USERNAME:\$DRAGEN_PASSWORD@license.edicogenome.com" : ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def ref = index ? "-r $index" : ''


    // Check FASTQ numbers match
    def num_r1 = fastq_1 instanceof List ? fastq_1.size() : 1
    def num_r2 = fastq_2 instanceof List ? fastq_2.size() : 1

    if ( fastq_2 && num_r1 != num_r2 ) {
        error "Number of R1 and R2 FASTQ files do not match for sample ${prefix}"
    }

    // DRAGEN only reads the FASTQ files named on its command line, so every staged pair goes into a
    // fastq-list CSV (one row per pair). With more than one pair each row gets its own read group ID.
    def rgid = meta.rgid ?: meta.id
    def rgsm = meta.rgsm ?: meta.id
    def rglb = meta.rglb ?: 'UnknownLibrary'
    def row_rgid = num_r1 > 1 ? "${rgid}.\${lane}" : "${rgid}"
    def inputs = "--fastq-list fastq_list.csv --fastq-list-sample-id ${rgsm}"
    """
    echo "RGID,RGSM,RGLB,Lane,Read1File,Read2File" > fastq_list.csv
    lane=0
    for r1 in input_S1_L001_R1_*.fastq.gz; do
        lane=\$((lane + 1))
        r2=\${r1/_R1_/_R2_}
        [ -e "\$r2" ] && r2="\$PWD/\$r2" || r2=""
        echo "${row_rgid},${rgsm},${rglb},\${lane},\$PWD/\${r1},\${r2}" >> fastq_list.csv
    done
    cat fastq_list.csv
    echo /opt/edico/bin/dragen \\
        $ref \\
        --output-directory ./ \\
        --output-file-prefix $prefix \\
        $lic \\
        $inputs \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: '\$(echo \$(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g")'
    END_VERSIONS
    """
}

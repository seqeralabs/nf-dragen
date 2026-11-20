nextflow.enable.types = true

process DRAGEN {
    tag "${id}"
    label 'dragen'

    input:
    record(
        id: String,
        lanes: List<String>,
        fastq_1: List<Path>,
        fastq_2: List<Path>
    )
    index: Path
    target_bed: Path?

    output:
    record(
        id: id,
        bam: file("${prefix}.bam", optional: true),
        bai: file("${prefix}.bam.bai", optional: true),
        vcf: file("${prefix}.vcf.gz", optional: true),
        tbi: file("${prefix}.vcf.gz.tbi", optional: true),
        vcf_filtered: file("${prefix}.hard-filtered.vcf.gz", optional: true),
        tbi_filtered: file("${prefix}.hard-filtered.vcf.gz.tbi", optional: true),
        metrics: files("${prefix}.*metrics.csv", optional: true),
    )

    topic:
    file('versions.yml') >> 'versions'

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: id
    if (fastq_2 && fastq_1.size() != fastq_2.size()) {
        error("Number of R1 and R2 FASTQ files do not match for sample ${id}")
    }
    def targets = target_bed ? "--vc-target-bed ${target_bed}" : ''
    // One row, and so one read group, per lane; DRAGEN needs a number in Lane, so the lane name goes into RGID
    def fastq_list = fastq_1.withIndex().collect { r1, i -> "'${id}.${lanes[i]},${id},UnknownLibrary,${i + 1},${r1.name},${fastq_2 ? fastq_2[i].name : ''}'" }.join(' ')
    """
    printf '%s\\n' 'RGID,RGSM,RGLB,Lane,Read1File,Read2File' ${fastq_list} > fastq_list.csv

    /opt/edico/bin/dragen \\
        -r ${index} \\
        --output-directory ./ \\
        --output-file-prefix ${prefix} \\
        --fastq-list fastq_list.csv \\
        ${targets} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: '\$(echo \$(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g")'
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: id
    if (fastq_2 && fastq_1.size() != fastq_2.size()) {
        error("Number of R1 and R2 FASTQ files do not match for sample ${id}")
    }
    def targets = target_bed ? "--vc-target-bed ${target_bed}" : ''
    def fastq_list = fastq_1.withIndex().collect { r1, i -> "'${id}.${lanes[i]},${id},UnknownLibrary,${i + 1},${r1.name},${fastq_2 ? fastq_2[i].name : ''}'" }.join(' ')
    """
    printf '%s\\n' 'RGID,RGSM,RGLB,Lane,Read1File,Read2File' ${fastq_list} > fastq_list.csv

    echo /opt/edico/bin/dragen \\
        -r ${index} \\
        --output-directory ./ \\
        --output-file-prefix ${prefix} \\
        --fastq-list fastq_list.csv \\
        ${targets} \\
        ${args}

    touch ${prefix}.bam ${prefix}.bam.bai
    touch ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
    touch ${prefix}.hard-filtered.vcf.gz ${prefix}.hard-filtered.vcf.gz.tbi
    touch ${prefix}.mapping_metrics.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: stub
    END_VERSIONS
    """
}

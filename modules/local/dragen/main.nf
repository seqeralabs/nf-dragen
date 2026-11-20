nextflow.enable.types = true

process DRAGEN {
    tag "${runId(id, tumour_id)}"
    label 'dragen'

    input:
    record(
        id: String,
        lanes: List<String>,
        fastq_1: List<Path>,
        fastq_2: List<Path>,
        tumour_id: String?,
        tumour_lanes: List<String>?,
        tumour_fastq_1: List<Path>?,
        tumour_fastq_2: List<Path>?
    )
    index: Path
    target_bed: Path?

    output:
    record(
        id: runId(id, tumour_id),
        bam: files("${prefix}*.bam", optional: true),
        bai: files("${prefix}*.bam.bai", optional: true),
        vcf: file("${prefix}.vcf.gz", optional: true),
        tbi: file("${prefix}.vcf.gz.tbi", optional: true),
        vcf_filtered: file("${prefix}.hard-filtered.vcf.gz", optional: true),
        tbi_filtered: file("${prefix}.hard-filtered.vcf.gz.tbi", optional: true),
        sv_vcf: file("${prefix}.sv.vcf.gz", optional: true),
        sv_tbi: file("${prefix}.sv.vcf.gz.tbi", optional: true),
        cnv_vcf: file("${prefix}.cnv.vcf.gz", optional: true),
        cnv_tbi: file("${prefix}.cnv.vcf.gz.tbi", optional: true),
        metrics: files("${prefix}.*metrics.csv", optional: true),
    )

    topic:
    file('versions.yml') >> 'versions'

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: runId(id, tumour_id)
    def fastq_lists = fastqLists(id, lanes, fastq_1, fastq_2, tumour_id, tumour_lanes, tumour_fastq_1, tumour_fastq_2)
    def targets = targetArgs(target_bed, tumour_id)
    """
    ${fastq_lists.tumour}
    ${fastq_lists.normal}

    /opt/edico/bin/dragen \\
        -r ${index} \\
        --output-directory ./ \\
        --output-file-prefix ${prefix} \\
        ${fastq_lists.args} \\
        ${targets} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: '\$(echo \$(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g")'
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: runId(id, tumour_id)
    def fastq_lists = fastqLists(id, lanes, fastq_1, fastq_2, tumour_id, tumour_lanes, tumour_fastq_1, tumour_fastq_2)
    def targets = targetArgs(target_bed, tumour_id)
    """
    ${fastq_lists.tumour}
    ${fastq_lists.normal}

    echo /opt/edico/bin/dragen \\
        -r ${index} \\
        --output-directory ./ \\
        --output-file-prefix ${prefix} \\
        ${fastq_lists.args} \\
        ${targets} \\
        ${args}

    touch ${prefix}.bam ${prefix}.bam.bai
    touch ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
    touch ${prefix}.hard-filtered.vcf.gz ${prefix}.hard-filtered.vcf.gz.tbi
    ${tumour_id ? "touch ${prefix}_tumor.bam ${prefix}_tumor.bam.bai ${prefix}.sv.vcf.gz ${prefix}.sv.vcf.gz.tbi" : ''}
    ${tumour_id && !target_bed ? "touch ${prefix}.cnv.vcf.gz ${prefix}.cnv.vcf.gz.tbi" : ''}
    touch ${prefix}.mapping_metrics.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: stub
    END_VERSIONS
    """
}

// A tumour/normal run is named after both samples; `id` is then the normal
def runId(id, tumour_id) {
    return tumour_id ? "${tumour_id}_vs_${id}" : id
}

// Commands that write the FASTQ list(s) and the DRAGEN arguments that read them
def fastqLists(id, lanes, fastq_1, fastq_2, tumour_id, tumour_lanes, tumour_fastq_1, tumour_fastq_2) {
    def header = "printf '%s\\n' 'RGID,RGSM,RGLB,Lane,Read1File,Read2File'"
    if (!tumour_id) {
        return [tumour: '', normal: "${header} ${fastqListRows(id, lanes, fastq_1, fastq_2)} > fastq_list.csv", args: '--fastq-list fastq_list.csv']
    }
    return [
        tumour: "${header} ${fastqListRows(tumour_id, tumour_lanes, tumour_fastq_1, tumour_fastq_2)} > tumour_fastq_list.csv",
        normal: "${header} ${fastqListRows(id, lanes, fastq_1, fastq_2)} > normal_fastq_list.csv",
        args: "--tumor-fastq-list tumour_fastq_list.csv --tumor-fastq-list-sample-id ${tumour_id} --fastq-list normal_fastq_list.csv --fastq-list-sample-id ${id}",
    ]
}

// fastq_list.csv rows of one sample: one row, and so one read group, per lane; DRAGEN needs a number in Lane, so the lane name goes into RGID
def fastqListRows(id, lanes, fastq_1, fastq_2) {
    if (fastq_2 && fastq_1.size() != fastq_2.size()) {
        error("Number of R1 and R2 FASTQ files do not match for sample ${id}")
    }
    return fastq_1.withIndex().collect { r1, i -> "'${id}.${lanes[i]},${id},UnknownLibrary,${i + 1},${r1.name},${fastq_2 ? fastq_2[i].name : ''}'" }.join(' ')
}

def targetArgs(target_bed, tumour_id) {
    if (!tumour_id) {
        return target_bed ? "--vc-target-bed ${target_bed}" : ''
    }
    // Exome CNV needs a panel of normals and self-normalization is WGS-only, so CNV runs only without a target BED
    return target_bed ? "--vc-target-bed ${target_bed} --sv-exome true --sv-call-regions-bed ${target_bed}" : '--enable-cnv true --cnv-use-somatic-vc-baf true --cnv-enable-self-normalization true'
}

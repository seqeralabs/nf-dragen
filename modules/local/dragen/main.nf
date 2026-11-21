nextflow.enable.types = true

process DRAGEN {
    tag "${runId(id, tumour_id)}"
    label 'dragen'

    input:
    record(
        id: String,
        fastq_1: Path,
        fastq_2: Path?,
        tumour_id: String?,
        tumour_fastq_1: Path?,
        tumour_fastq_2: Path?
    )
    index: Path
    target_bed: Path?

    output:
    record(
        id: runId(id, tumour_id),
        bam: files("${prefix}*.bam", optional: true),
        bai: files("${prefix}*.bam.bai", optional: true),
        // DRAGEN 4.4.4 germline and RNA runs write no unfiltered VCF; kept for somatic runs (unverified)
        vcf: file("${prefix}.vcf.gz", optional: true),
        tbi: file("${prefix}.vcf.gz.tbi", optional: true),
        vcf_filtered: file("${prefix}.hard-filtered.vcf.gz", optional: true),
        tbi_filtered: file("${prefix}.hard-filtered.vcf.gz.tbi", optional: true),
        sv_vcf: file("${prefix}.sv.vcf.gz", optional: true),
        sv_tbi: file("${prefix}.sv.vcf.gz.tbi", optional: true),
        cnv_vcf: file("${prefix}.cnv.vcf.gz", optional: true),
        cnv_tbi: file("${prefix}.cnv.vcf.gz.tbi", optional: true),
        metrics: files("${prefix}.{*metrics*,*_hist*,*_mean_cov*}.csv", optional: true),
        // Germline: named in the DRAGEN 4.4.4 run log. BAF: the log names `.baf`, the v4.4 docs `.baf.bw` (somatic too, unverified)
        ploidy_vcf: file("${prefix}.ploidy.vcf.gz", optional: true),
        roh_bed: file("${prefix}.roh.bed", optional: true),
        baf: files("${prefix}{,.hard-filtered}.baf{,.bw}", optional: true),
        // RNA splice junctions (v4.4 docs, unverified names): filtered and unfiltered SJ.out.tab, chimeric junctions
        junctions: files("${prefix}.{SJ.out.tab*,Chimeric.out.junction}", optional: true),
        dragen_files: file("${prefix}.dragen_files.txt"),
    )

    topic:
    file('versions.yml') >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def lic = task.ext.lic_server ?: ''
    prefix = task.ext.prefix ?: runId(id, tumour_id)
    def fastq_lists = fastqLists(id, fastq_1, fastq_2, tumour_id, tumour_fastq_1, tumour_fastq_2)
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
        ${lic} \\
        ${args}

    ls -lR > ${prefix}.dragen_files.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: '\$(echo \$(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g")'
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    // quoted below, so the echo prints the secret names, not the credentials
    def lic = task.ext.lic_server ?: ''
    prefix = task.ext.prefix ?: runId(id, tumour_id)
    def fastq_lists = fastqLists(id, fastq_1, fastq_2, tumour_id, tumour_fastq_1, tumour_fastq_2)
    def targets = targetArgs(target_bed, tumour_id)
    // Germline and RNA names are those DRAGEN 4.4.4 wrote in a real run; somatic names come from the DRAGEN docs (unverified).
    // Coverage reports (RNA too): whole genome, plus target_bed with --vc-target-bed; tumour/normal runs write one per sample
    // The coverage histogram and mean-coverage names touched below are unverified until a real run
    def suffixes = tumour_id ? ['_tumor', '_normal'] : ['']
    def coverage = (target_bed ? ['wgs', 'target_bed'] : ['wgs']).collectMany { region -> suffixes.collect { suffix -> "${prefix}.${region}_coverage_metrics${suffix}.csv" } }.join(' ')
    """
    ${fastq_lists.tumour}
    ${fastq_lists.normal}

    echo /opt/edico/bin/dragen \\
        -r ${index} \\
        --output-directory ./ \\
        --output-file-prefix ${prefix} \\
        ${fastq_lists.args} \\
        ${targets} \\
        '${lic}' \\
        ${args}

    touch ${prefix}.bam ${prefix}.bam.bai ${coverage} ${prefix}.mapping_metrics.csv ${prefix}.time_metrics.csv
    ${args.contains('--enable-variant-caller') ? "touch ${prefix}.hard-filtered.vcf.gz ${prefix}.hard-filtered.vcf.gz.tbi ${prefix}.vc_metrics.csv" : ''}
    ${tumour_id ? "touch ${prefix}_tumor.bam ${prefix}_tumor.bam.bai ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi ${prefix}.sv.vcf.gz ${prefix}.sv.vcf.gz.tbi" : "touch ${prefix}.fastqc_metrics.csv ${prefix}.trimmer_metrics.csv ${prefix}.wgs_overall_mean_cov.csv ${prefix}.wgs_contig_mean_cov.csv ${prefix}.wgs_fine_hist.csv ${prefix}.wgs_hist.csv ${prefix}.fragment_length_hist.csv"}
    ${tumour_id && !target_bed ? "touch ${prefix}.cnv.vcf.gz ${prefix}.cnv.vcf.gz.tbi ${prefix}.cnv_metrics.csv" : ''}
    ${!tumour_id && args.contains('--enable-variant-caller') ? "touch ${prefix}.vc_hethom_ratio_metrics.csv ${prefix}.ploidy_estimation_metrics.csv ${prefix}.roh_metrics.csv ${prefix}.allele_transition_noise_metrics.csv ${prefix}.ploidy.vcf.gz ${prefix}.roh.bed ${prefix}.baf ${prefix}.hard-filtered.baf" : ''}
    ${tumour_id && args.contains('--enable-variant-caller') ? "touch ${prefix}.baf.bw ${prefix}.hard-filtered.baf.bw" : ''}
    ${args.contains('--enable-rna') ? "touch ${prefix}.SJ.out.tab ${prefix}.Chimeric.out.junction" : ''}

    ls -lR > ${prefix}.dragen_files.txt

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
def fastqLists(id, fastq_1, fastq_2, tumour_id, tumour_fastq_1, tumour_fastq_2) {
    def header = "printf '%s\\n' 'RGID,RGSM,RGLB,Lane,Read1File,Read2File'"
    if (!tumour_id) {
        return [tumour: '', normal: "${header} ${fastqListRow(id, fastq_1, fastq_2)} > fastq_list.csv", args: '--fastq-list fastq_list.csv']
    }
    return [
        tumour: "${header} ${fastqListRow(tumour_id, tumour_fastq_1, tumour_fastq_2)} > tumour_fastq_list.csv",
        normal: "${header} ${fastqListRow(id, fastq_1, fastq_2)} > normal_fastq_list.csv",
        args: "--tumor-fastq-list tumour_fastq_list.csv --tumor-fastq-list-sample-id ${tumour_id} --fastq-list normal_fastq_list.csv --fastq-list-sample-id ${id}",
    ]
}

// The one fastq_list.csv row of a sample: one read group named after the sample
def fastqListRow(id, fastq_1, fastq_2) {
    return "'${id},${id},UnknownLibrary,1,${fastq_1.name},${fastq_2 ? fastq_2.name : ''}'"
}

def targetArgs(target_bed, tumour_id) {
    if (!tumour_id) {
        return target_bed ? "--vc-target-bed ${target_bed}" : ''
    }
    // Exome CNV needs a panel of normals and self-normalization is WGS-only, so CNV runs only without a target BED
    return target_bed ? "--vc-target-bed ${target_bed} --sv-exome true --sv-call-regions-bed ${target_bed}" : '--enable-cnv true --cnv-use-somatic-vc-baf true --cnv-enable-self-normalization true'
}

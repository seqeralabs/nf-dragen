nextflow.enable.types = true

process DRAGEN {
    tag "${runId(id, tumour_id)}"
    label 'dragen'

    input:
    record(
        id: String,
        fastq_1: List<Path>,
        fastq_2: List<Path>,
        tumour_id: String?,
        tumour_fastq_1: List<Path>?,
        tumour_fastq_2: List<Path>?,
        index: Path,
        target_bed: Path?,
        dragen_args: String?,
        dragen_lic_server: String?,
        output_prefix: String?
    )

    stage:

    // FASTQ pair n goes to r1/<n>/ and r2/<n>/ (tumour: tumour/...): lanes, mates and a tumour and its normal often share file names
    stageAs fastq_1, 'r1/?/*'
    stageAs fastq_2, 'r2/?/*'
    stageAs tumour_fastq_1, 'tumour/r1/?/*'
    stageAs tumour_fastq_2, 'tumour/r2/?/*'

    output:
    record(
        id: runId(id, tumour_id),
        bam: files("${prefix}*.bam"),
        bai: files("${prefix}*.bam.bai"),
        vcf_filtered: file("${prefix}.hard-filtered.vcf.gz", optional: true),
        tbi_filtered: file("${prefix}.hard-filtered.vcf.gz.tbi", optional: true),
        sv_vcf: file("${prefix}.sv.vcf.gz", optional: true),
        sv_tbi: file("${prefix}.sv.vcf.gz.tbi", optional: true),
        cnv_vcf: file("${prefix}.cnv.vcf.gz", optional: true),
        cnv_tbi: file("${prefix}.cnv.vcf.gz.tbi", optional: true),
        mapping_metrics: file("${prefix}.mapping_metrics.csv"), // also matched by metrics below; required, so a missing file fails the task
        metrics: files("${prefix}.{*metrics*,*_hist*,*_mean_cov*}.csv", optional: true),
        ploidy_vcf: file("${prefix}.ploidy.vcf.gz", optional: true),
        ploidy_tbi: file("${prefix}.ploidy.vcf.gz.tbi", optional: true),
        roh_bed: file("${prefix}.roh.bed", optional: true),
        baf: files("${prefix}{,.hard-filtered}.baf.bw", optional: true),
        // Germline targeted callers (e.g. CYP2D6, SMN)
        targeted_vcf: file("${prefix}.targeted.vcf.gz", optional: true),
        targeted_tbi: file("${prefix}.targeted.vcf.gz.tbi", optional: true),
        targeted_json: file("${prefix}.targeted.json", optional: true),
        callable_bed: file("${prefix}.wgs_somatic_callable_regions.bed", optional: true),
        junctions: files("${prefix}.{SJ.out.tab,unfiltered.SJ.out.tab,Chimeric.out.junction}", optional: true),
    )

    topic:
    record(process: "${task.process}", tool: 'dragen', version: eval(version_cmd)) >> 'versions'

    script:
    def args = dragen_args ?: ''
    version_cmd = 'echo $(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g"'
    prefix = output_prefix ?: runId(id, tumour_id)
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
        ${dragen_lic_server ?: ''} \\
        ${args}
    """

    stub:
    def args = dragen_args ?: ''
    // the stub never calls the binary, which exists only on DRAGEN hosts
    version_cmd = 'echo stub'
    prefix = output_prefix ?: runId(id, tumour_id)
    def fastq_lists = fastqLists(id, fastq_1, fastq_2, tumour_id, tumour_fastq_1, tumour_fastq_2)
    def targets = targetArgs(target_bed, tumour_id)
    // Quoted, so the echo prints the variable names, not the credentials
    def lic_server = dragen_lic_server ? "'${dragen_lic_server}'" : ''
    // File names are those DRAGEN 4.4.4 wrote in real germline WGS, somatic WGS and RNA runs; target_bed names come from the v4.4 docs
    def rna = args.contains('--enable-rna')
    def vc = args.contains('--enable-variant-caller')
    def suffixes = tumour_id ? ['_tumor', '_normal'] : ['']
    def coverage = (target_bed ? ['wgs', 'target_bed'] : ['wgs']).collectMany { region -> ['coverage_metrics', 'overall_mean_cov', 'contig_mean_cov', 'fine_hist', 'hist'].collectMany { kind -> suffixes.collect { suffix -> "${prefix}.${region}_${kind}${suffix}.csv" } } }.join(' ')
    """
    ${fastq_lists.tumour}
    ${fastq_lists.normal}

    echo /opt/edico/bin/dragen \\
        -r ${index} \\
        --output-directory ./ \\
        --output-file-prefix ${prefix} \\
        ${fastq_lists.args} \\
        ${targets} \\
        ${lic_server} \\
        ${args}

    touch ${prefix}.bam ${prefix}.bam.bai ${coverage} ${prefix}.mapping_metrics.csv ${prefix}.time_metrics.csv ${prefix}.fastqc_metrics.csv ${prefix}.trimmer_metrics.csv
    ${rna ? "touch ${prefix}.SJ.out.tab ${prefix}.unfiltered.SJ.out.tab ${prefix}.Chimeric.out.junction" : "touch ${prefix}.fragment_length_hist.csv"}
    ${vc ? "touch ${prefix}.hard-filtered.vcf.gz ${prefix}.hard-filtered.vcf.gz.tbi ${prefix}.vc_metrics.csv ${prefix}.ploidy.vcf.gz ${prefix}.ploidy.vcf.gz.tbi ${prefix}.ploidy_estimation_metrics.csv ${prefix}.allele_transition_noise_metrics.csv ${prefix}.baf.bw ${prefix}.hard-filtered.baf.bw" : ''}
    ${vc && !tumour_id ? "touch ${prefix}.vc_hethom_ratio_metrics.csv ${prefix}.roh_metrics.csv ${prefix}.roh.bed ${prefix}.targeted.vcf.gz ${prefix}.targeted.vcf.gz.tbi ${prefix}.targeted.json" : ''}
    ${tumour_id ? "touch ${prefix}_tumor.bam ${prefix}_tumor.bam.bai ${prefix}.sv.vcf.gz ${prefix}.sv.vcf.gz.tbi ${prefix}.sv_metrics.csv" : ''}
    ${tumour_id && !target_bed ? "touch ${prefix}.cnv.vcf.gz ${prefix}.cnv.vcf.gz.tbi ${prefix}.cnv_metrics.csv ${prefix}.wgs_somatic_callable_regions.bed" : ''}
    """
}

// A tumour/normal run is named after both samples; `id` is then the normal
def runId(id: String, tumour_id: String?) -> String {
    return tumour_id ? "${tumour_id}_vs_${id}" : id
}

// Commands that write the FASTQ list(s) and the DRAGEN arguments that read them
def fastqLists(id: String, fastq_1: List<Path>, fastq_2: List<Path>, tumour_id: String?, tumour_fastq_1: List<Path>?, tumour_fastq_2: List<Path>?) -> Map<String, String> {
    def header = "printf '%s\\n' 'RGID,RGSM,RGLB,Lane,Read1File,Read2File'"
    if (!tumour_id) {
        return [tumour: '', normal: "${header} ${fastqListRows(id, '', fastq_1, fastq_2)} > fastq_list.csv", args: '--fastq-list fastq_list.csv']
    }
    return [
        tumour: "${header} ${fastqListRows(tumour_id, 'tumour/', tumour_fastq_1, tumour_fastq_2)} > tumour_fastq_list.csv",
        normal: "${header} ${fastqListRows(id, '', fastq_1, fastq_2)} > normal_fastq_list.csv",
        args: "--tumor-fastq-list tumour_fastq_list.csv --tumor-fastq-list-sample-id ${tumour_id} --fastq-list normal_fastq_list.csv --fastq-list-sample-id ${id}",
    ]
}

// fastq_list.csv rows of one sample: one read group per FASTQ pair; a single pair keeps RGID = sample.
// Paths follow the stage: patterns, not the staged variables, which map a file given twice to one name (and whose .name is the staged path)
def fastqListRows(id: String, dir: String, fastq_1: List<Path>, fastq_2: List<Path>) -> String {
    return fastq_1.withIndex().collect { r1, i -> "'${fastq_1.size() == 1 ? id : "${id}.${i + 1}"},${id},UnknownLibrary,${i + 1},${dir}r1/${i + 1}/${r1.name.tokenize("/").last()},${fastq_2 ? "${dir}r2/${i + 1}/${fastq_2[i].name.tokenize("/").last()}" : ''}'" }.join(' ')
}

def targetArgs(target_bed: Path?, tumour_id: String?) -> String {
    if (!tumour_id) {
        return target_bed ? "--vc-target-bed ${target_bed}" : ''
    }
    // Exome CNV needs a panel of normals and self-normalization is WGS-only, so CNV runs only without a target BED
    return target_bed ? "--vc-target-bed ${target_bed} --sv-exome true --sv-call-regions-bed ${target_bed}" : '--enable-cnv true --cnv-use-somatic-vc-baf true --cnv-enable-self-normalization true'
}

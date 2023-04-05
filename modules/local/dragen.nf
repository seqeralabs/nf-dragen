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
    yum -y update
    yum -y install kernel-devel kernel-headers

    yum -y install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm

    curl -L https://www.xilinx.com/bin/public/openDownload?filename=xrt_202020.2.8.832_7.4.1708-x86_64-xrt.rpm -o xrt_202020.2.8.832_7.4.1708-x86_64-xrt.rpm
    curl -L https://www.xilinx.com/bin/public/openDownload?filename=xrt_202020.2.8.832_7.4.1708-x86_64-azure.rpm -o xrt_202020.2.8.832_7.4.1708-x86_64-azure.rpm

    yum -y install xrt*.rpm

    curl -L https://www.xilinx.com/bin/public/openDownload?filename=xilinx-u250-gen3x16-xdma-platform-2.1-3.noarch.rpm.tar.gz  -o xilinx-u250-gen3x16-xdma-platform-2.1-3.noarch.rpm.tar.gz
    curl -L https://www.xilinx.com/bin/public/openDownload?filename=xilinx-u250-gen3x16-xdma-validate-2.1-3005608.1.noarch.rpm -o xilinx-u250-gen3x16-xdma-validate-2.1-3005608.1.noarch.rpm
    tar -zxvf xilinx-u250-gen3x16-xdma-platform-2.1-3.noarch.rpm.tar.gz

    yum -y install xilinx*.rpm

    systemctl enable mpd
    systemctl start mpd

    /opt/xilinx/xrt/bin/xbutil host_mem -d 0 --enable --size 1G

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

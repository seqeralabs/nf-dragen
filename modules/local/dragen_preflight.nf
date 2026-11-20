process DRAGEN_PREFLIGHT {
    tag "dragen"
    label 'dragen'

    output:
    path "dragen_preflight.txt", emit: report
    path "versions.yml"        , emit: versions

    script:
    // Fail fast, before any FASTQ or reference staging, if the DRAGEN host
    // mounts or the FPGA device are not visible inside the task container.
    """
    fail() {
        cat dragen_preflight.txt >&2
        echo "ERROR: DRAGEN pre-flight failed: \$1. Check aws.batch.volumes in the Launchpad Nextflow config (see docs/launchpad_aws_f2.config) and that the task runs on the DRAGEN queue." >&2
        exit 1
    }

    echo "== DRAGEN pre-flight on \$(hostname) ==" > dragen_preflight.txt
    ls -l /opt/edico/bin/dragen >> dragen_preflight.txt 2>&1 || fail "/opt/edico/bin/dragen not mounted"
    ls -l /lib64/libdragen.so.* >> dragen_preflight.txt 2>&1 || fail "/lib64/libdragen.so.* not mounted"
    ls -l /dev/dragen*          >> dragen_preflight.txt 2>&1 || fail "no /dev/dragen FPGA device visible"
    /opt/edico/bin/dragen --version >> dragen_preflight.txt 2>&1 || fail "dragen --version failed"
    cat dragen_preflight.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: '\$(echo \$(/opt/edico/bin/dragen --version 2>&1) | sed -e "s/dragen Version //g")'
    END_VERSIONS
    """

    stub:
    """
    echo "DRAGEN pre-flight (stub)" > dragen_preflight.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: 'stub'
    END_VERSIONS
    """
}

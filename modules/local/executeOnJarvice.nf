String executeOnJarvice(cmd) {
    String taskset_cmd = "/sys/bus/pci/devices/\$(/opt/xilinx/xrt/bin/xbutil scan | tail -5 | grep xilinx | cut -d ' ' -f3)/local_cpulist \\\n"
    String full_cmd = cmd
    """
    # Replace the shell script
    cp .command.sh .backup.command.sh
    echo "Working dir: \$NXF_CHDIR"
    echo '$full_cmd' > .command.sh

    echo "about to run Jarvice CLI"
    /usr/local/bin/entrypoint \\
        --api-host https://illumina-beta.nimbix.net/api \\
        --machine n1-2 \\
        --dragen-app illumina-dragen_4_3_6 \\
        --google-sa default \\
        --job-priority normal \\
        --s3-access-key null \\
        --s3-secret-key null \\
        -- \\
        "\$NXF_CHDIR/.command.run"
    echo "finished running command"
    """.stripIndent()
}

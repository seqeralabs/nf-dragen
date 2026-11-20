# seqeralabs/nf-dragen: Usage

> Every pipeline parameter is documented in [parameters.md](parameters.md).

## Introduction

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with a header row, as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

Columns can be in any order. Samplesheets written for nf-core/sarek work unchanged, and columns the pipeline does not know (for example nf-core/rnaseq's `strandedness`) are ignored.

| Column     | Required | Description                                                                                                                                                 |
| ---------- | -------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`   | yes      | Sample name, without spaces. Rows with the same `sample` are lanes (or re-sequencing runs) of that sample and are analysed together.                        |
| `fastq_1`  | yes      | Path to the read 1 FASTQ file. Must exist, be gzipped and end in `.fastq.gz` or `.fq.gz`.                                                                   |
| `fastq_2`  | no       | Path to the read 2 FASTQ file, for paired-end data. Give it for every row of a sample or for none.                                                          |
| `patient`  | no       | Patient ID, without spaces. Defaults to the sample name.                                                                                                    |
| `status`   | no       | `0` for normal, `1` for tumour. Defaults to `0`.                                                                                                            |
| `lane`     | no       | Lane name, for example `1` or `L001`. Each lane is its own read group, `<sample>.<lane>`. Defaults to the row's position within the sample (`1`, `2`, ...). |
| `sex`      | no       | Sex of the patient, for example `XX` or `XY`.                                                                                                               |
| `seq_type` | no       | `dna` or `rna`. Defaults to `--seq_type` (default `dna`). One samplesheet can hold both DNA and RNA samples.                                                |

The pipeline stops before running anything, with a message naming the problem, when:

- a required column is missing, or a FASTQ file does not exist;
- the rows of one sample disagree on `patient`, `status`, `sex` or `seq_type`;
- some rows of one sample have `fastq_2` and others do not.

### Minimal samplesheet

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
SAMPLE1,SAMPLE1_L001_R1.fastq.gz,SAMPLE1_L001_R2.fastq.gz
SAMPLE1,SAMPLE1_L002_R1.fastq.gz,SAMPLE1_L002_R2.fastq.gz
SAMPLE2,SAMPLE2_L001_R1.fastq.gz,
```

`SAMPLE1` is paired-end and sequenced on two lanes; `SAMPLE2` is single-end.

### Full samplesheet

```csv title="samplesheet.csv"
patient,sex,status,sample,lane,fastq_1,fastq_2,seq_type
PATIENT1,XX,0,PATIENT1_NORMAL,L001,P1N_L001_R1.fastq.gz,P1N_L001_R2.fastq.gz,dna
PATIENT1,XX,1,PATIENT1_TUMOUR,L001,P1T_L001_R1.fastq.gz,P1T_L001_R2.fastq.gz,dna
PATIENT1,XX,1,PATIENT1_TUMOUR,L002,P1T_L002_R1.fastq.gz,P1T_L002_R2.fastq.gz,dna
PATIENT2,XY,0,PATIENT2_RNA,L001,P2R_L001_R1.fastq.gz,P2R_L001_R2.fastq.gz,rna
```

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Germline, somatic and RNA runs

The pipeline picks the DRAGEN run for each DNA sample from `patient` and `status`:

- A patient with a normal and a tumour gets one **somatic** (tumour/normal) run. The normal gets no germline run.
- A patient with several tumours and one normal gets one somatic run per tumour, each with that normal.
- A normal without a tumour gets a **germline** run.
- A tumour without a normal, or a patient with two or more normals, stops the pipeline before any DNA task runs, with a message naming the patient. Tumour-only calling is not supported.

RNA samples ignore `patient` and `status` and each get an **RNA** run.

Somatic runs call small variants, structural variants (SV) and copy number variants (CNV). HRD, TMB, HLA and MSI are off; switch them on by adding the DRAGEN options (for example `--enable-hrd true`, `--enable-tmb true`, `--enable-hla true`, `--msi-command tumor-normal`) to `ext.args` of `.*:DRAGEN_DNA:DRAGEN_SOMATIC` in a custom config (see [Custom Tool Arguments](#custom-tool-arguments)).

For exome or panel data, give the target regions with `--target_bed targets.bed`. Germline runs get it as `--vc-target-bed`; somatic runs also as `--sv-call-regions-bed` and `--cnv-target-bed`. RNA runs never get it.

## MultiQC report

MultiQC runs after DRAGEN on the DRAGEN metrics and the software versions (see [Output](output.md)). It runs in a Docker container, so the pipeline needs Docker on the machine that runs MultiQC; the DRAGEN processes have no container and run on the host. Skip it with `--skip_multiqc`.

## Configuring the reference

We can configure the reference in two ways:

1. Supply a `--fasta`, from which Dragen will build a genome index for alignment.
2. If you have already performed a pipeline run using `1.`, the reference index will be written to a results directory at `${params.outdir}/genome/index/dragen_index_dna` or `${params.outdir}/genome/index/dragen_index_rna`. You can reuse these indexes by supplying the parameter `--dragen_index_dna` and/or `--dragen_index_rna` which will skip 1 for DNA and RNA respectively.

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run seqeralabs/nf-dragen --input samplesheet.csv --fasta genome.fa --outdir results
```

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run seqeralabs/nf-dragen -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
fasta: './genome.fa'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull seqeralabs/nf-dragen
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [seqeralabs/nf-dragen releases page](https://github.com/seqeralabs/nf-dragen/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Note that multiple profiles can be loaded; the order of arguments is important, as they are loaded in sequence.

- `test`
  - A profile with a complete configuration for automated testing
  - Uses the small synthetic data in `tests/data/`, so needs no other parameters except `--outdir`. Only works with `-stub`, because DRAGEN is not available outside its licensed environment.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

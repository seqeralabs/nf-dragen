# seqeralabs/nf-dragen: Usage

> Every pipeline parameter is documented in [parameters.md](parameters.md).

## Introduction

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with a header row, as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

Columns can be in any order, and columns the pipeline does not know (for example nf-core/sarek's `lane` or nf-core/rnaseq's `strandedness`) are ignored. Each sample is one row, with one FASTQ pair or one single-end FASTQ. A nf-core/sarek samplesheet with several lanes per sample must first be merged into one row per sample, for example by concatenating each sample's lane FASTQs.

| Column     | Required | Description                                                                                                  |
| ---------- | -------- | ------------------------------------------------------------------------------------------------------------ |
| `sample`   | yes      | Sample name, without spaces, unique within the samplesheet.                                                  |
| `fastq_1`  | yes      | Path to the read 1 FASTQ file. Must exist, be gzipped and end in `.fastq.gz` or `.fq.gz`.                    |
| `fastq_2`  | no       | Path to the read 2 FASTQ file, for paired-end data.                                                          |
| `patient`  | no       | Patient ID, without spaces. Defaults to the sample name.                                                     |
| `status`   | no       | `0` for normal, `1` for tumour. Defaults to `0`.                                                             |
| `sex`      | no       | Sex of the patient, for example `XX` or `XY`.                                                                |
| `seq_type` | no       | `dna` or `rna`. Defaults to `--seq_type` (default `dna`). One samplesheet can hold both DNA and RNA samples. |

The pipeline stops before running anything, with a message naming the problem, when:

- a required column is missing, or a FASTQ file does not exist;
- a sample name is on more than one row.

### Minimal samplesheet

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
SAMPLE1,SAMPLE1_R1.fastq.gz,SAMPLE1_R2.fastq.gz
SAMPLE2,SAMPLE2_R1.fastq.gz,
```

`SAMPLE1` is paired-end; `SAMPLE2` is single-end.

### Full samplesheet

```csv title="samplesheet.csv"
patient,sex,status,sample,fastq_1,fastq_2,seq_type
PATIENT1,XX,0,PATIENT1_NORMAL,P1N_R1.fastq.gz,P1N_R2.fastq.gz,dna
PATIENT1,XX,1,PATIENT1_TUMOUR,P1T_R1.fastq.gz,P1T_R2.fastq.gz,dna
PATIENT2,XY,0,PATIENT2_RNA,P2R_R1.fastq.gz,P2R_R2.fastq.gz,rna
```

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Germline, somatic and RNA runs

The pipeline picks the DRAGEN run for each DNA sample from `patient` and `status`:

- A patient with a normal and a tumour gets one **somatic** (tumour/normal) run. The normal gets no germline run.
- A patient with several tumours and one normal gets one somatic run per tumour, each with that normal.
- A normal without a tumour gets a **germline** run.
- A tumour without a normal, or a patient with two or more normals, stops the pipeline before any DNA task runs, with a message naming the patient. Tumour-only calling is not supported.

RNA samples ignore `patient` and `status` and each get an **RNA** run.

Somatic runs call small variants, structural variants (SV) and, for whole-genome data (no `--target_bed`), copy number variants (CNV) with self-normalization. HRD, TMB, HLA and MSI are off; switch them on by adding the DRAGEN options (for example `--enable-hrd true`, `--enable-tmb true`, `--enable-hla true`, `--msi-command tumor-normal`) to `ext.args` of `.*:DRAGEN_SOMATIC` in a custom config (see [Custom Tool Arguments](#custom-tool-arguments)).

For exome or panel data, give the target regions with `--target_bed targets.bed`. Germline runs get it as `--vc-target-bed`; somatic runs also as `--sv-call-regions-bed`. Somatic exome runs skip CNV calling, because DRAGEN needs a panel of normals for exome CNV, which the pipeline does not support. RNA runs never get it.

## MultiQC report

MultiQC runs after DRAGEN on the DRAGEN metrics and the software versions (see [Output](output.md)). It runs in a Docker container, so the pipeline needs Docker on the machine that runs MultiQC; the DRAGEN processes have no container and run on the host. Skip it with `--skip_multiqc`.

## Configuring the reference

We can configure the reference in two ways:

1. Supply a `--fasta`, from which Dragen will build a genome index for alignment. The FASTA must not contain ALT contigs, for example use the GRCh38 no-ALT analysis set: DRAGEN refuses to build from ALT contigs without a liftover file, and the pipeline does not pass one. The DNA hash table is built with CNV k-mers, so CNV calling can use it; a DNA index built by an earlier version of the pipeline lacks them and must be rebuilt for CNV.
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

Each process has default CPUs, memory and time requests in `conf/base.config`. Processes labelled `dragen` have fixed resources (8 CPUs, 96 GB of memory and 24 hours) that do not change between attempts. Most other labels multiply their memory or time by the attempt number, so a retry asks for twice the resources of the first attempt.

A task that exits with code 130 to 145 or 104, which usually means it was killed, for example for running out of memory, is retried once. Any other exit code is not retried. The pipeline then uses `finish`: tasks that are already running complete, no new tasks start, and the run exits with an error naming the sample.

To change the resources of a process, pass a config file with `-c`. Select processes by label with `withLabel`, or by name with `withName`:

```groovy
process {
    withLabel: dragen {
        memory = 120.GB
    }
}
```

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

[DRAGEN options](dragen_options.md) lists every DRAGEN option the pipeline passes, checked against the Illumina DRAGEN v4.4 docs, and the recommended options it leaves out.

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

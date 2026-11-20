[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A526.10.6-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)

> This is a proof-of-concept pipeline under active development. Syntax, layout and parameters may change without notice.

## Introduction

**seqeralabs/nf-dragen** runs the [Illumina DRAGEN](https://emea.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html) licensed suite of tools on FASTQ samples. It produces germline calls for DNA samples, tumour/normal somatic calls, and RNA alignments, and summarises the DRAGEN metrics in a MultiQC report.

The pipeline is built with [Nextflow](https://www.nextflow.io). DRAGEN runs from the AWS Marketplace pay-as-you-go DRAGEN AMI on AWS F2 instances, so a real run needs a compute environment that uses that AMI (see [Requirements](#requirements)).

## Pipeline steps

1. **Samplesheet parsing.** Reads the samplesheet, which uses the same columns as nf-core/sarek, and checks it before any task runs.
2. **Run selection.** Groups DNA samples by patient: a normal without a tumour gets a germline run, and each tumour gets a somatic run paired with the patient's normal. RNA samples get an RNA run.
3. **Hash tables.** Builds the DRAGEN DNA and RNA hash tables from `--fasta`, unless you supply them with `--dragen_index_dna` and `--dragen_index_rna`.
4. **DRAGEN.** Runs germline and somatic DNA calling (somatic runs also call structural and copy number variants) and RNA alignment.
5. **MultiQC.** Summarises the DRAGEN metrics and the software versions in one report. Skip it with `--skip_multiqc`.

```mermaid
flowchart LR
    SS[Samplesheet] --> PS[Parse and validate]
    PS --> SEL[Select run per patient]
    FA[fasta] -.->|built when no index is given| HT[Hash tables]
    SEL -->|normal only| GER[DRAGEN germline]
    SEL -->|tumour and normal| SOM[DRAGEN somatic]
    SEL -->|RNA samples| RNA[DRAGEN RNA]
    HT --> GER
    HT --> SOM
    HT --> RNA
    GER --> MQ[MultiQC]
    SOM --> MQ
    RNA --> MQ
```

## Quickstart

Run the test profile in stub mode from the repository root. Stub mode does not run DRAGEN, so it needs neither the AMI nor a licence:

```bash
nextflow run . -profile test -stub --outdir results
```

A real run needs the DRAGEN AMI on F2 instances (see [Requirements](#requirements)):

```bash
nextflow run . --input samplesheet.csv --fasta genome.fa --outdir results
```

Read [docs/usage.md](docs/usage.md) for the samplesheet format and run options.

## Requirements

- Nextflow 26.10.6 or later.
- A DRAGEN-capable compute environment: the AWS Marketplace pay-as-you-go DRAGEN AMI on AWS F2 instances, in a region where F2 instances are available. AWS meters the DRAGEN licence, so no credentials are needed. The instances need outbound internet access for DRAGEN's licence check. The DRAGEN processes carry the `dragen` label, which sends them to those instances.
- Docker on the machine that runs MultiQC. The DRAGEN processes run on the host and need no container.

### Seqera Platform setup

You can use Seqera Batch Forge to create an AWS Batch compute environment with F2 instances. In the Seqera UI, go to `Compute Environments` -> `Add Compute Environment`, fill in your AWS Batch settings and paste the AWS Marketplace pay-as-you-go DRAGEN AMI id:

![Seqera enable DRAGEN](docs/images/tower_ce_enable_dragen.png)

## Documentation

- [Usage](docs/usage.md): samplesheet format, run modes, configuration and running the pipeline.
- [Output](docs/output.md): the files the pipeline writes to `--outdir`.
- [Parameters](docs/parameters.md): every pipeline parameter.

## Credits

nf-dragen was originally written by [Harshil Patel](https://github.com/drpatelh), [Graham Wright](https://github.com/gwright99) and [Paolo Di Tommaso](https://github.com/pditommaso), [Seqera Labs](https://seqera.io/).

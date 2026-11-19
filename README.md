[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)

> THIS IS A PROOF-OF-CONCEPT REPOSITORY THAT IS UNDER ACTIVE DEVELOPMENT. SYNTAX, ORGANISATION AND LAYOUT MAY CHANGE WITHOUT NOTICE!

## Introduction

**nf-dragen** is a simple, proof-of-concept pipeline to run the [Illumina DRAGEN](https://emea.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html) licensed suite of tools.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker containers making installation trivial and results highly reproducible. This pipeline has only been tested on AWS Batch.

## Integration with Seqera Platform

Seqera Batch Forge can create an AWS Batch compute environment with a separate queue of DRAGEN F2 instances. Any process with `label 'dragen'` runs on that queue.

### 1. Choose a licensing model

| `--dragen_license`      | AMI                                                     | Credentials                                                                                       |
| ----------------------- | ------------------------------------------------------- | ------------------------------------------------------------------------------------------------- |
| `marketplace`           | Illumina DRAGEN from AWS Marketplace (subscribe first)  | None. AWS meters the DRAGEN software per instance-hour                                            |
| `byol` (default)        | Private AMI supplied by Illumina                        | Pipeline secrets `DRAGEN_USERNAME` and `DRAGEN_PASSWORD`, passed to `--lic-server` by the pipeline |

For BYOL, add both secrets in Seqera (`Secrets` -> `Add Pipeline Secret`) and select them when adding the pipeline:

![Seqera Secrets DRAGEN username](docs/images/tower_secrets_dragen_username.png)

![Seqera Secrets DRAGEN password](docs/images/tower_secrets_dragen_password.png)

BYOL is the default, as in earlier releases. The run fails at launch with `Required secrets are missing` if the secrets are not selected. For a Marketplace AMI, set `--dragen_license marketplace`: no secrets are needed.

### 2. Create the compute environment

In Seqera, go to `Compute Environments` -> `Add Compute Environment` -> AWS Batch (Batch Forge) and set:

- **Region** with F2 instances, e.g. `us-east-1`
- **Enable DRAGEN**: the DRAGEN AMI ID and an F2 instance type (`f2.6xlarge` is the smallest and cheapest)
- **Wave**, **Fusion v2** and **fast instance storage (NVMe)** enabled
- **Provisioning model**: on-demand is recommended. The DRAGEN software fee is the same on Spot, and a Spot reclaim means paying it again for the retry
- **Nextflow config**: leave empty (see step 3)

![Seqera enable DRAGEN](docs/images/tower_ce_enable_dragen.png)

### 3. Add and launch the pipeline

In Seqera, go to `Launchpad` -> `Add Pipeline`, select the compute environment and, for BYOL only, the two secrets. In **Nextflow config**, paste [`docs/launchpad_aws_f2.config`](docs/launchpad_aws_f2.config), with the `libdragen.so` version matching the AMI's DRAGEN version.

This config mounts the DRAGEN install from the F2 host into the task container and keeps DRAGEN's scratch on the local NVMe disk. These settings depend on the AMI, so they live in the Launchpad rather than in the pipeline.

> [!WARNING]
> Do not put `aws.batch.volumes` in the compute environment's Nextflow config. For DRAGEN-enabled compute environments, Seqera Platform appends its own `aws.batch.volumes` after that field, and because the setting is a single value, Platform's line replaces yours without libdragen or `/ephemeral`. The Launchpad config is applied last, so set it there.

Before any data is staged, the pipeline runs `DRAGEN_PREFLIGHT` on the DRAGEN queue. It checks that `/opt/edico/bin/dragen`, `/lib64/libdragen.so.*` and `/dev/dragen` are visible and that `dragen --version` works, and writes `pipeline_info/dragen_preflight.txt`. Skip it with `--skip_dragen_preflight`.

Useful options:

- `--dragen_index_dna` / `--dragen_index_rna`: reuse a pre-built hash table instead of building one (the hash table must match the DRAGEN version)

### Smoke test

To check a new compute environment as cheaply as possible:

```bash
nextflow run seqeralabs/nf-dragen -profile test_smoke --outdir <OUTDIR>
```

`test_smoke` runs one small sample against a yeast reference and skips FastQC and MultiQC. Cap the compute environment's maximum CPUs at one F2 instance (24 for `f2.6xlarge`) so the test can only start one DRAGEN instance.

## Credits

nf-dragen was originally written by [Harshil Patel](https://github.com/drpatelh) and [Graham Wright](https://github.com/gwright99) and [Paolo Di Tommasso](https://github.com/pditommaso), [Seqera Labs](https://seqera.io/).

## Citations

The nf-core pipeline template was used to create the skeleton of this pipeline but there are no plans to contribute it to nf-core at this point.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

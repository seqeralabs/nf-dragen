# seqeralabs/nf-dragen: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of seqeralabs/nf-dragen, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- (10)[https://github.com/seqeralabs/nf-dragen/pull/10] - Use nf-core template v3.1.1
- `--dragen_license marketplace|byol` (default `marketplace`). Only BYOL passes `--lic-server` and requires the `DRAGEN_USERNAME`/`DRAGEN_PASSWORD` secrets; AWS Marketplace AMIs need no credentials
- `DRAGEN_PREFLIGHT` check (host mounts, `/dev/dragen`, `dragen --version`) that gates all DRAGEN tasks; skip with `--skip_dragen_preflight`
- `test_smoke` profile: one sample, no QC, DRAGEN tasks serial with no retries, for minimal-cost compute environment checks
- `docs/launchpad_aws_f2.config`: Launchpad Nextflow config for DRAGEN on AWS Batch F2 (host mounts, `/ephemeral` scratch, Fusion settings)
- `--skip_pangenome_validation` and default `--ht-num-threads <task.cpus> --shared-thread-odirect-output true` for DRAGEN mapping

### `Fixed`

- Fixed error `Access to 'FASTQC.out' is undefined since the workflow 'FASTQC' has not been invoked before accessing the output attribute` when `-skip_fastqc` enabled by adjusting channel generation

### `Dependencies`

### `Deprecated`

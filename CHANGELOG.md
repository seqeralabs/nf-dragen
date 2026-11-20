# seqeralabs/nf-dragen: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of seqeralabs/nf-dragen, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- (10)[https://github.com/seqeralabs/nf-dragen/pull/10] - Use nf-core template v3.1.1
- `--dragen_license marketplace|byol` (default `byol`, as before). Only BYOL passes `--lic-server` and requires the `DRAGEN_USERNAME`/`DRAGEN_PASSWORD` secrets; AWS Marketplace AMIs need no credentials
- `DRAGEN_PREFLIGHT` check (host mounts, `/dev/dragen`, `dragen --version`) that gates all DRAGEN tasks; skip with `--skip_dragen_preflight`
- `test_smoke` profile: one sample, no QC, for minimal-cost compute environment checks
- `docs/launchpad_aws_f2.config`: Launchpad Nextflow config for DRAGEN on AWS Batch F2 (host mounts, `/ephemeral` scratch, Fusion settings, DRAGEN retries for host start-up)

### `Fixed`

- DRAGEN now reads every FASTQ pair of a sample: all pairs go into a `--fastq-list` CSV. Before, a sample with several samplesheet rows only had its first pair analysed (DRAGEN does not pick up the other staged files by name)
- Fixed error `Access to 'FASTQC.out' is undefined since the workflow 'FASTQC' has not been invoked before accessing the output attribute` when `-skip_fastqc` enabled by adjusting channel generation

### `Dependencies`

### `Deprecated`

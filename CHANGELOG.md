# seqeralabs/nf-dragen: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of seqeralabs/nf-dragen, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- (10)[https://github.com/seqeralabs/nf-dragen/pull/10] - Use nf-core template v3.1.1

### `Removed`

- FastQC and MultiQC steps, and the params `skip_fastqc`, `skip_multiqc`, `multiqc_config`, `multiqc_title`, `multiqc_logo`, `max_multiqc_email_size` and `multiqc_methods_description`. The MultiQC report and `samplesheet.csv` are no longer published.
- iGenomes support and the `genome`, `igenomes_base` and `igenomes_ignore` params. Use `--fasta` or `--dragen_index_dna`/`--dragen_index_rna` instead.
- The `help`, `help_full`, `show_hidden`, `version`, `validate_params`, `pipelines_testdata_base_path` and `trace_report_suffix` params. Parameters are always validated.
- The container, conda and `debug` profiles, and the `test_full` profile. No process uses a container.
- The timeline, report, trace and DAG files in `pipeline_info/`; Seqera Platform records these itself.

### `Changed`

- The pipeline requires Nextflow 26.10.6 or later.
- The software versions report is written to `pipeline_info/software_versions.yml` and lists the versions of every process, including the RNA steps.
- The `test` profile uses small synthetic data in `tests/data/` and only works with `-stub`.

### `Fixed`

- Fixed error `Access to 'FASTQC.out' is undefined since the workflow 'FASTQC' has not been invoked before accessing the output attribute` when `-skip_fastqc` enabled by adjusting channel generation

### `Dependencies`

- nf-schema 3.1.0

### `Deprecated`

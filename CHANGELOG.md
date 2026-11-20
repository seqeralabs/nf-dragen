# seqeralabs/nf-dragen: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of seqeralabs/nf-dragen, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- (10)[https://github.com/seqeralabs/nf-dragen/pull/10] - Use nf-core template v3.1.1
- Somatic tumour/normal calling: for DNA samples, a patient with a normal and one or more tumours gets one somatic run per tumour (small variants, SV and CNV; HRD, TMB, HLA and MSI can be switched on through `ext.args`), published to `dragen/dna_fastq_to_vcf/<tumour>_vs_<normal>.*`. A normal without a tumour gets a germline run. A tumour without a normal, or two normals for one patient, stop the pipeline with an error naming the patient.
- `--target_bed`: exome or panel target regions for germline and somatic runs (not RNA runs).

### `Removed`

- FastQC and MultiQC steps, and the params `skip_fastqc`, `skip_multiqc`, `multiqc_config`, `multiqc_title`, `multiqc_logo`, `max_multiqc_email_size` and `multiqc_methods_description`. The MultiQC report and `samplesheet.csv` are no longer published.
- iGenomes support and the `genome`, `igenomes_base` and `igenomes_ignore` params. Use `--fasta` or `--dragen_index_dna`/`--dragen_index_rna` instead.
- The `help`, `help_full`, `show_hidden`, `version`, `validate_params`, `pipelines_testdata_base_path` and `trace_report_suffix` params. Parameters are always validated.
- The container, conda and `debug` profiles, and the `test_full` profile. No process uses a container.
- The timeline, report, trace and DAG files in `pipeline_info/`; Seqera Platform records these itself.
- The `config_profile_name` and `config_profile_description` params. They had no effect since the nf-core utils were removed. Launch configs that still pass them now fail parameter validation.

### `Changed`

- `--input` accepts nf-core/sarek samplesheets unchanged, in any column order: the optional columns `patient`, `status`, `lane` and `sex` are new, and unknown columns are ignored. Existing `sample,fastq_1,fastq_2` samplesheets still work. When `lane` is given, each lane's read group is `<sample>.<lane>`.
- A samplesheet is rejected before anything runs when the rows of one sample disagree on `patient`, `status`, `sex` or `seq_type`, or when some rows of a sample have `fastq_2` and others do not.

- The pipeline requires Nextflow 26.10.6 or later.
- The software versions report is written to `pipeline_info/software_versions.yml` and lists the versions of every process, including the RNA steps.
- The `test` profile uses small synthetic data in `tests/data/` and only works with `-stub`.
- DRAGEN reads every lane of a sample, each as its own read group (`<sample>.<lane>`). The BAM index and DRAGEN's `*metrics.csv` files are now published; FASTQ files are not.

### `Fixed`

- Fixed error `Access to 'FASTQC.out' is undefined since the workflow 'FASTQC' has not been invoked before accessing the output attribute` when `-skip_fastqc` enabled by adjusting channel generation

### `Dependencies`

- nf-schema 3.1.0

### `Deprecated`

# seqeralabs/nf-dragen: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of seqeralabs/nf-dragen, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [#10](https://github.com/seqeralabs/nf-dragen/pull/10) - Use nf-core template v3.1.1
- Samplesheet columns, in any column order: the optional columns `patient`, `tissue` (`germline` or `tumor`, default `germline`) and `sex` (`male`, `female` or `unknown`, default `unknown`) are new, and unknown columns (for example `lane`) are ignored. A samplesheet with a `status` column is rejected; use `tissue`. `sample,fastq_1,fastq_2` samplesheets still work.
- Several rows with the same `sample` are that sample's FASTQ pairs, for germline, tumour/normal and RNA runs. Each row is one DRAGEN read group: `<sample>.1`, `<sample>.2`, ... in samplesheet order; a sample with one row has the read group ID `<sample>`. The rows of a sample must agree on `patient`, `tissue`, `sex` and `seq_type` and be all paired-end or all single-end, otherwise the pipeline stops before any task with an error naming the sample and the column.
- Tumour/normal calling for DNA samples: a patient with one normal and one tumour gets one somatic run (small variants and SV, plus CNV for whole-genome data, not for runs with `--target_bed`), published as `dragen/somatic/<tumour>_vs_<normal>.*`. A normal without a tumour gets a germline run. The pipeline stops before any task, with an error naming the patient or samples, for a tumour without a normal, two normals or two tumours for one patient.
- `--dragen_germline_args`, `--dragen_somatic_args` and `--dragen_rna_args`: the tunable DRAGEN options of each run type (for example to switch on HRD, TMB, HLA or MSI for somatic runs). Setting one replaces its default, so repeat the default options you want to keep.
- `--dragen_hashtable_dna_args` (default `--ht-build-cnv-hashtable true`) and `--dragen_hashtable_rna_args` (default `--ht-build-rna-hashtable true`): the DRAGEN options of the DNA and RNA hash table builds, used when no `--dragen_index_dna`/`--dragen_index_rna` is given. Setting one replaces its default.
- `--target_bed`: exome or panel target regions for germline and somatic runs (not RNA runs). Somatic exome runs use DRAGEN's exome SV mode.
- The DNA hash table is built with CNV k-mers. `--fasta` must not contain ALT contigs, for example use the GRCh38 no-ALT analysis set.
- `--dragen_license` (`marketplace` by default, or `byol`) for Illumina's private DRAGEN AMI. Setting it to `byol`, or using `-profile byol`, declares the Nextflow secrets `DRAGEN_USERNAME` and `DRAGEN_PASSWORD`; every DRAGEN task then passes them to `--lic-server`. It is untested on real DRAGEN.
- A run with DNA or RNA samples but neither `--fasta` nor the matching `--dragen_index_dna`/`--dragen_index_rna` stops before any DRAGEN task with an error naming a sample and both parameters.
- Published DRAGEN files per run: BAM and index, `*metrics*.csv`, coverage histogram and mean-coverage CSVs, VCFs and SV/CNV calls where the mode produces them, the ploidy VCF, targeted-caller calls, runs of homozygosity and B-allele frequencies (germline), the callable-regions BED (somatic), and the splice junctions (RNA). FASTQ files are not published.
- Index files `dragen/germline/germline_index.csv`, `dragen/somatic/somatic_index.csv` and `dragen/rna/rna_index.csv`: one row per DRAGEN run with its patient and sample (somatic: patient, tumour and normal) and the published path of each file.
- MultiQC report on the DRAGEN metrics and the software versions, published to `multiqc/`, with the params `skip_multiqc`, `multiqc_config`, `multiqc_title` and `multiqc_logo`. MultiQC runs in a container, so the pipeline enables Docker; DRAGEN processes run on the host.
- The software versions report `pipeline_info/software_versions.yml` lists the DRAGEN version of every DRAGEN and hash-table process, and the Nextflow and pipeline versions. It does not list MultiQC.
- Seqera Platform's Reports tab shows the MultiQC report and the DRAGEN QC metrics, listed in `tower.yml`.
- README: the Seqera Platform section has the Batch Forge settings and the Nextflow config to paste into the launch form for AWS F2 DRAGEN runs (volumes, `/ephemeral` scratch, Fusion settings, retry rule that excludes exit 134). Add it to your launch config.

### `Changed`

- The pipeline requires Nextflow 26.10.6 or later and nf-schema 3.1.0.
- Results go to the output directory given by `-output-dir <dir>` (one dash), or `outputDir` in a config file or a Seqera Platform launch. `--outdir` fails the launch.
- DRAGEN runs with the AWS Marketplace pay-as-you-go licence, so it needs no Nextflow secrets and no `--lic-server`. Use `--dragen_license byol` or `-profile byol` for a private AMI.
- DRAGEN options are not read from `ext.args` or `ext.prefix` in custom configs; use `--dragen_germline_args`, `--dragen_somatic_args` or `--dragen_rna_args`.
- RNA runs do not call variants; `dragen/rna/` holds no VCFs.
- DRAGEN tasks that abort with exit 134 (SIGABRT) are not retried; exit 143 and the other resource-class codes are.
- DNA hash tables built by an earlier version lack CNV k-mers and must be rebuilt for CNV calling.
- The `test` profile uses small public data from nf-core/test-datasets, read from raw.githubusercontent.com (needs network access), and only works with `-stub`.

### `Removed`

- The FastQC step, and the params `skip_fastqc`, `max_multiqc_email_size` and `multiqc_methods_description`. `samplesheet.csv` is not published.
- iGenomes support and the `genome`, `igenomes_base` and `igenomes_ignore` params. Use `--fasta` or `--dragen_index_dna`/`--dragen_index_rna` instead.
- The `help`, `help_full`, `show_hidden`, `version`, `validate_params`, `pipelines_testdata_base_path` and `trace_report_suffix` params. Parameters are always validated.
- The container, conda and `debug` profiles, and the `test_full` profile. Docker is always on; no other container runtime is supported.
- The timeline, report, trace and DAG files in `pipeline_info/`; Seqera Platform records these itself.
- The params `config_profile_name` and `config_profile_description`; launch configs that still pass them fail parameter validation.
- The `outdir` param (see `Changed`).
- The `publish_dir_mode` param: results are always copied. Set `workflow.output.mode` in your own config to change that.
- The `monochrome_logs` param: set `validation.monochromeLogs = true` in your own config.
- The `lane` column is not used.
- The `DRAGEN_USERNAME` and `DRAGEN_PASSWORD` secrets and `--lic-server`, unless `--dragen_license byol` or `-profile byol` is set.

### `Dependencies`

- nf-schema 3.1.0

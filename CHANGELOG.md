# seqeralabs/nf-dragen: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of seqeralabs/nf-dragen, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- (10)[https://github.com/seqeralabs/nf-dragen/pull/10] - Use nf-core template v3.1.1
- DRAGEN coverage histogram and mean-coverage CSVs (`*_hist.csv`, `*_mean_cov.csv`) are published next to the `*metrics*.csv` files and passed to MultiQC when DRAGEN writes them.
- More DRAGEN outputs are published next to each sample's DRAGEN outputs when DRAGEN writes them: for germline runs the ploidy VCF (`<prefix>.ploidy.vcf.gz`), runs of homozygosity (`<prefix>.roh.bed`) and B-allele frequencies (`<prefix>.baf.bw`, `<prefix>.hard-filtered.baf.bw`, also for somatic runs); for RNA runs the splice junctions (`<prefix>.SJ.out.tab`, `<prefix>.Chimeric.out.junction`).
- Published DRAGEN outputs now match what DRAGEN 4.4.4 writes: the RNA unfiltered splice-junction table (`<prefix>.unfiltered.SJ.out.tab`), the ploidy VCF index (`<prefix>.ploidy.vcf.gz.tbi`, germline and somatic), the germline targeted-caller calls (`<prefix>.targeted.vcf.gz(.tbi)`, `<prefix>.targeted.json`) and the whole-genome somatic callable regions (`<prefix>.wgs_somatic_callable_regions.bed`) are now published, and the somatic SV metrics show in the Seqera Platform Reports tab. Germline B-allele frequencies are published as `<prefix>.baf.bw` and `<prefix>.hard-filtered.baf.bw` (BigWig), as for somatic runs.
- Somatic tumour/normal calling: for DNA samples, a patient with one normal and one tumour gets one somatic run (small variants and SV, plus CNV for whole-genome data, not for runs with `--target_bed`; HRD, TMB, HLA and MSI can be switched on through `ext.args`), published to `dragen/dna_fastq_to_vcf/<tumour>_vs_<normal>.*`. A normal without a tumour gets a germline run. A tumour without a normal, two normals for one patient, or two tumours for one patient, stop the pipeline with an error naming the patient.
- The DNA hash table is built with CNV k-mers (`--ht-build-cnv-hashtable true`), so CNV calling can use it. A DNA index built by an earlier version lacks them and must be rebuilt for CNV. `--fasta` must not contain ALT contigs, for example use the GRCh38 no-ALT analysis set.
- `--target_bed`: exome or panel target regions for germline and somatic runs (not RNA runs).
- MultiQC report on the DRAGEN metrics and the software versions, published to `multiqc/`, with the params `skip_multiqc`, `multiqc_config`, `multiqc_title` and `multiqc_logo`. MultiQC runs in a container, so the pipeline now enables Docker (`docker.enabled = true`); DRAGEN processes still run on the host. The software versions report also lists the Nextflow and pipeline versions.
- Seqera Platform's Reports tab shows the MultiQC report and the DRAGEN QC metrics (mapping, variant calling, whole-genome and target region coverage, CNV and run time), listed in `tower.yml`. DRAGEN's tumour/normal coverage metrics (`*_coverage_metrics_tumor.csv`, `*_coverage_metrics_normal.csv`) are now published and passed to MultiQC.

- `-profile byol` (without it, DRAGEN uses the AWS Marketplace licence). It is for Illumina's private DRAGEN AMI: every DRAGEN task then needs the Nextflow secrets `DRAGEN_USERNAME` and `DRAGEN_PASSWORD` and passes them to `--lic-server`. `byol` is untested on real DRAGEN.
- A run with DNA or RNA samples but neither `--fasta` nor the matching `--dragen_index_dna`/`--dragen_index_rna` now stops before any DRAGEN task with an error naming a sample and both parameters.

- README: the Seqera Platform section now has the Batch Forge settings and the Nextflow config to paste into the launch form for AWS F2 DRAGEN runs (volumes, `/ephemeral` scratch, Fusion settings, extra retry). Add it to your launch config.

### `Removed`

- The FastQC step, and the params `skip_fastqc`, `max_multiqc_email_size` and `multiqc_methods_description`. `samplesheet.csv` is no longer published.
- iGenomes support and the `genome`, `igenomes_base` and `igenomes_ignore` params. Use `--fasta` or `--dragen_index_dna`/`--dragen_index_rna` instead.
- The `help`, `help_full`, `show_hidden`, `version`, `validate_params`, `pipelines_testdata_base_path` and `trace_report_suffix` params. Parameters are always validated.
- The container, conda and `debug` profiles, and the `test_full` profile. Docker is always on (see `Added`); no other container runtime is supported.
- The timeline, report, trace and DAG files in `pipeline_info/`; Seqera Platform records these itself.
- The `config_profile_name` and `config_profile_description` params. They had no effect since the nf-core utils were removed. Launch configs that still pass them now fail parameter validation.
- The `publish_dir_mode` param: results are always copied to the output directory. To change that, set `workflow.output.mode` in your own config.
- The `outdir` param: `--outdir` now fails the launch. Pass `-output-dir <dir>` (one dash) on every run instead, or set `outputDir` in a config file or a Seqera Platform launch (the **Output directory** field, or `tw launch --output-dir`).
- The `monochrome_logs` param: set `validation.monochromeLogs = true` in your own config instead.
- The `dragen_license` param: use `-profile byol` instead of `--dragen_license byol`, and drop `--dragen_license marketplace` (the default).
- The `lane` column: a samplesheet takes one row per sample, with one FASTQ pair or one single-end FASTQ. A sample name on more than one row is now an error, so merge the lanes of a sample (for example by concatenating its FASTQs) into one row first. A `lane` column is ignored like any unknown column.
- The `DRAGEN_USERNAME` and `DRAGEN_PASSWORD` Nextflow secrets and the `--lic-server` argument are no longer used by default. DRAGEN runs from the AWS Marketplace pay-as-you-go AMI, which meters the licence through AWS. To keep using a private AMI with your own licence, run with `-profile byol` (see `Added`).

### `Changed`

- RNA runs no longer call variants; `dragen/rna_fastq_to_bam/` holds no VCFs.
- DRAGEN process names: custom configs that set `ext.args` or `ext.prefix` must now target `.*:DRAGEN_DNA` (germline), `.*:DRAGEN_SOMATIC` (tumour/normal), `.*:DRAGEN_RNA` (RNA), `.*:DRAGEN_BUILDHASHTABLE_DNA` and `.*:DRAGEN_BUILDHASHTABLE_RNA`, instead of `.*:DRAGEN_DNA:DRAGEN`, `.*:DRAGEN_DNA:DRAGEN_SOMATIC`, `.*:DRAGEN_RNA:DRAGEN` and `.*:DRAGEN_{DNA,RNA}:DRAGEN_BUILDHASHTABLE`.

- Somatic runs no longer pass `--Aligner.clip-pe-overhang 2`, so paired-end overhangs are not clipped (it likely caused a DRAGEN 4.4.4 SV caller abort on a 1-base alignment, and Illumina's v4.4 tumour/normal recipes do not use it).
- DRAGEN tasks that abort with exit 134 (SIGABRT) are no longer retried; exit 143 and the other resource-class codes still are. If you copied the README's F2 launch config, change its `errorStrategy` to `task.exitStatus in ((130..145) - 134) + 104 + Integer.MAX_VALUE`.
- Somatic exome runs (`--target_bed`) now pass `--sv-exome true`, as the DRAGEN v4.4 WES recipe does, so SV calling uses exome mode.
- `--input` accepts nf-core/sarek samplesheets unchanged, in any column order: the optional columns `patient`, `status` and `sex` are new, and unknown columns are ignored. Existing `sample,fastq_1,fastq_2` samplesheets with one row per sample still work.
- A samplesheet is rejected before anything runs when a sample name is on more than one row; the message names the sample.

- The pipeline requires Nextflow 26.10.6 or later.
- The software versions report is written to `pipeline_info/software_versions.yml` and lists the versions of every process, including the RNA steps.
- The `test` profile uses small public data from nf-core/test-datasets, read from raw.githubusercontent.com (needs network access), and only works with `-stub`.
- DRAGEN reads each sample as one read group named after the sample. The BAM index and DRAGEN's `*metrics.csv` files are now published; FASTQ files are not.

### `Dependencies`

- nf-schema 3.1.0

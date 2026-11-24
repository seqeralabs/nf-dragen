# seqeralabs/nf-dragen: Output

The pipeline writes its results to the output directory given with `-output-dir` (one dash) or `outputDir` in a config file. It is required on every run.

```
<output-dir>/
├── genome/index/      # hash tables built from --fasta
├── dragen/
│   ├── germline/      # germline DNA runs
│   ├── somatic/       # tumour/normal runs
│   └── rna/           # RNA runs
├── multiqc/           # MultiQC report, unless --skip_multiqc
└── pipeline_info/     # software versions
```

Files are copied to the output directory. To change that, set `workflow.output.mode` in your own config, for example `workflow.output.mode = 'symlink'`.

On Seqera Platform, the MultiQC report and the DRAGEN QC files appear in the run's Reports tab (listed in `tower.yml`).

## Hash tables

| Path                             | Contents                                  |
| -------------------------------- | ----------------------------------------- |
| `genome/index/dragen_index_dna/` | DNA hash table, when built from `--fasta` |
| `genome/index/dragen_index_rna/` | RNA hash table, when built from `--fasta` |

`genome/index/` holds a hash table only when the pipeline built it. Pass it back with `--dragen_index_dna` or `--dragen_index_rna` to skip the build in a later run. An index you supply is not copied.

## DRAGEN

`dragen/` holds the outputs of each run in `germline/`, `somatic/` and `rna/`. In the tables below, `<prefix>` is the sample name, or `<tumour>_vs_<normal>` for a somatic run.

Each samplesheet row of a sample is one read group in the BAM. The read group ID is the sample name for a sample with one row, and `<sample>.1`, `<sample>.2`, ... for a sample with several rows.

Each file is published when DRAGEN writes it.

### Germline

| File                                                       | Contents                                                                                     |
| ---------------------------------------------------------- | -------------------------------------------------------------------------------------------- |
| `germline_index.csv`                                       | One row per germline sample (see [Index files](#index-files))                                |
| `<prefix>.bam`, `<prefix>.bam.bai`                         | Alignments                                                                                   |
| `<prefix>.hard-filtered.vcf.gz(.tbi)`                      | Small variant calls                                                                          |
| `<prefix>.ploidy.vcf.gz(.tbi)`                             | Ploidy estimator calls                                                                       |
| `<prefix>.baf.bw`, `<prefix>.hard-filtered.baf.bw`         | B-allele frequencies of the small variant calls as BigWig (the second holds PASS calls only) |
| `<prefix>.roh.bed`                                         | Runs of homozygosity                                                                         |
| `<prefix>.targeted.vcf.gz(.tbi)`, `<prefix>.targeted.json` | Calls of DRAGEN's targeted callers (genes such as CYP2D6 and SMN1/2)                         |
| `<prefix>.*metrics*.csv`                                   | QC metrics (see [QC metrics](#qc-metrics))                                                   |

### Somatic

A somatic run writes the outputs of a tumour/normal pair, named `<tumour>_vs_<normal>`, with a BAM for each of the two samples. It writes the same files as a germline run except `.roh.bed` and `.targeted.*`, plus:

| File                                        | Contents                                                                                      |
| ------------------------------------------- | --------------------------------------------------------------------------------------------- |
| `somatic_index.csv`                         | One row per tumour/normal run (see [Index files](#index-files))                               |
| `<prefix>.bam(.bai)`                        | Normal BAM                                                                                    |
| `<prefix>_tumor.bam(.bai)`                  | Tumour BAM                                                                                    |
| `<prefix>.sv.vcf.gz(.tbi)`                  | Structural variants                                                                           |
| `<prefix>.cnv.vcf.gz(.tbi)`                 | Copy number variants, whole-genome runs only                                                  |
| `<prefix>.wgs_somatic_callable_regions.bed` | Regions with enough coverage in both samples to call somatic variants, whole-genome runs only |

### RNA

RNA runs write no VCFs.

| File                               | Contents                                                 |
| ---------------------------------- | -------------------------------------------------------- |
| `rna_index.csv`                    | One row per RNA sample (see [Index files](#index-files)) |
| `<prefix>.bam`, `<prefix>.bam.bai` | Alignments                                               |
| `<prefix>.SJ.out.tab`              | Filtered splice junctions                                |
| `<prefix>.unfiltered.SJ.out.tab`   | Unfiltered splice junctions                              |
| `<prefix>.Chimeric.out.junction`   | Chimeric alignments                                      |
| `<prefix>.*metrics*.csv`           | QC metrics (see [QC metrics](#qc-metrics))               |

### QC metrics

| Runs                                      | Files                                                                                                                                                                                                                                                                                                                                |
| ----------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| All                                       | `<prefix>.mapping_metrics.csv` (for RNA runs it includes the RNA-specific mapping metrics), `.time_metrics.csv`, `.fastqc_metrics.csv`, `.trimmer_metrics.csv`, and the whole-genome coverage reports `.wgs_coverage_metrics.csv`, `.wgs_overall_mean_cov.csv`, `.wgs_contig_mean_cov.csv`, `.wgs_fine_hist.csv` and `.wgs_hist.csv` |
| Germline and somatic                      | `<prefix>.fragment_length_hist.csv`, `.vc_metrics.csv`, `.ploidy_estimation_metrics.csv` and `.allele_transition_noise_metrics.csv`                                                                                                                                                                                                  |
| Germline and somatic, with `--target_bed` | The same five coverage reports for the target regions: `<prefix>.target_bed_coverage_metrics.csv`, `_overall_mean_cov`, `_contig_mean_cov`, `_fine_hist` and `_hist`                                                                                                                                                                 |
| Germline                                  | `<prefix>.vc_hethom_ratio_metrics.csv` and `.roh_metrics.csv`                                                                                                                                                                                                                                                                        |
| Somatic                                   | `<prefix>.sv_metrics.csv`, and the coverage reports once per sample, ending in `_tumor` and `_normal` (for example `<prefix>.wgs_coverage_metrics_tumor.csv`, `<prefix>.wgs_hist_normal.csv`)                                                                                                                                        |
| Whole-genome somatic                      | `<prefix>.cnv_metrics.csv`                                                                                                                                                                                                                                                                                                           |

### Index files

`germline_index.csv`, `somatic_index.csv` and `rna_index.csv` list each DRAGEN run with the samples it came from and the published path of each of its files, as CSV with a header.

| Columns                                                                                                                                                                                                                                          | Contents                                                     |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------ |
| `id`, `patient`, `sample`                                                                                                                                                                                                                        | Germline and RNA: the file prefix, patient and sample        |
| `id`, `patient`, `tumour`, `normal`                                                                                                                                                                                                              | Somatic: `<tumour>_vs_<normal>`, patient and the two samples |
| `bam`, `bai`, `vcf_filtered`, `tbi_filtered`, `sv_vcf`, `sv_tbi`, `cnv_vcf`, `cnv_tbi`, `mapping_metrics`, `metrics`, `ploidy_vcf`, `ploidy_tbi`, `roh_bed`, `baf`, `targeted_vcf`, `targeted_tbi`, `targeted_json`, `callable_bed`, `junctions` | Then one column per DRAGEN output                            |

A column is empty (or `[]`) when the run did not write that file. Columns that can hold several files (`bam`, `bai`, `metrics`, `baf`, `junctions`) hold a list such as `[<path>, <path>]`. A somatic run's `bam` and `bai` hold the normal and the tumour file.

### Not published

- DRAGEN's logs, replay and usage files, and MD5 sums.
- `<prefix>.metrics.json` and `<prefix>.insert-stats.tab`, which hold the same metrics as the CSVs.
- The CNV intermediates: `.tn.tsv.gz`, `tumor.*`, `.cnv.gff3` and `.cnv.igv_session.xml`.

## MultiQC

`multiqc/` is written unless you pass `--skip_multiqc`. It summarises every DRAGEN `*metrics*.csv` and the software versions. With `--skip_dragen` it holds only the software versions.

| Path                          | Contents                                |
| ----------------------------- | --------------------------------------- |
| `multiqc/multiqc_report.html` | MultiQC report                          |
| `multiqc/multiqc_data/`       | The report's data, as text files        |
| `multiqc/multiqc_plots/`      | The report's plots, as PNG, SVG and PDF |

Set `--multiqc_title` to name the report (the HTML file name then starts with the title), `--multiqc_config` to add a MultiQC config and `--multiqc_logo` to add a logo.

## Pipeline info

| Path                                  | Contents                                                                    |
| ------------------------------------- | --------------------------------------------------------------------------- |
| `pipeline_info/software_versions.yml` | DRAGEN version of every DRAGEN process, plus Nextflow and pipeline versions |

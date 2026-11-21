# seqeralabs/nf-dragen: Output

The pipeline writes these files to the output directory given with `-output-dir` (one dash) or `outputDir` in a config file. It is required on every run:

```
<output-dir>/
├── genome/
│   └── index/
│       ├── dragen_index_dna/          # DNA hash table, when built from --fasta
│       └── dragen_index_rna/          # RNA hash table, when built from --fasta
├── dragen/
│   ├── dna_fastq_to_vcf/              # DNA samples
│   │   ├── <sample>.bam
│   │   ├── <sample>.bam.bai
│   │   ├── <sample>.*metrics*.csv     # DRAGEN QC metrics, see below
│   │   ├── <sample>.hard-filtered.vcf.gz
│   │   ├── <sample>.hard-filtered.vcf.gz.tbi
│   │   ├── <tumour>_vs_<normal>.*        # somatic runs: the same files as above, plus
│   │   ├── <tumour>_vs_<normal>.vcf.gz(.tbi)      # unfiltered variants
│   │   ├── <tumour>_vs_<normal>.sv.vcf.gz(.tbi)   # structural variants
│   │   └── <tumour>_vs_<normal>.cnv.vcf.gz(.tbi)  # copy number variants, whole-genome somatic runs only
│   └── rna_fastq_to_bam/              # RNA samples: <sample>.bam, .bam.bai and *metrics*.csv, no VCFs
├── multiqc/                           # unless --skip_multiqc
│   ├── multiqc_report.html            # MultiQC report: DRAGEN metrics and software versions
│   ├── multiqc_data/                  # the report's data, as text files
│   └── multiqc_plots/                 # the report's plots, as PNG, SVG and PDF
└── pipeline_info/
    └── software_versions.yml          # versions of the software every process ran
```

- `genome/index/` holds a hash table only when the pipeline built it. Pass it back with `--dragen_index_dna` or `--dragen_index_rna` to skip the build in a later run. An index you supply is not copied.
- `dragen/` holds the DRAGEN outputs of each sample, named by the sample name. A somatic run writes the outputs of the tumour/normal pair under `dragen/dna_fastq_to_vcf/`, named `<tumour>_vs_<normal>`, with a BAM for each of the two samples. Each sample is one read group in the BAM, with the sample name as its ID.
- DRAGEN QC metrics (`<prefix>` is the sample name, or `<tumour>_vs_<normal>` for a somatic run):
  - every run: `<prefix>.mapping_metrics.csv` (for RNA runs it includes the RNA-specific mapping metrics), `<prefix>.time_metrics.csv` and `<prefix>.wgs_coverage_metrics.csv`;
  - germline and RNA runs: `<prefix>.fastqc_metrics.csv` and `<prefix>.trimmer_metrics.csv`;
  - germline and somatic runs: `<prefix>.vc_metrics.csv`, plus `<prefix>.target_bed_coverage_metrics.csv` with `--target_bed`;
  - germline runs: `<prefix>.vc_hethom_ratio_metrics.csv`, `<prefix>.ploidy_estimation_metrics.csv`, `<prefix>.roh_metrics.csv` and `<prefix>.allele_transition_noise_metrics.csv`;
  - somatic runs write the coverage metrics once per sample, ending in `_tumor` and `_normal` (e.g. `<prefix>.wgs_coverage_metrics_tumor.csv`);
  - whole-genome somatic runs: `<prefix>.cnv_metrics.csv`.
  - when DRAGEN writes them, the coverage histograms and mean coverage: `<prefix>.<region>_overall_mean_cov.csv`, `<prefix>.<region>_contig_mean_cov.csv`, `<prefix>.<region>_fine_hist.csv`, `<prefix>.<region>_hist.csv` and `<prefix>.fragment_length_hist.csv`.
- `<prefix>.dragen_files.txt` (in both `dragen/` folders): a recursive `ls -l` listing of the DRAGEN task directory after DRAGEN finished, so it also shows the files the pipeline does not publish.
- On Seqera Platform, the MultiQC report and these DRAGEN QC files appear in the run's Reports tab (listed in `tower.yml`).
- `multiqc/` summarises every DRAGEN `*metrics*.csv` and the software versions. With `--skip_dragen` it holds only the software versions. Set `--multiqc_title` to name the report (the HTML file name then starts with the title), `--multiqc_config` to add a MultiQC config and `--multiqc_logo` to add a logo.
- Files are copied to the output directory. To change that, set `workflow.output.mode` in your own config (for example `workflow.output.mode = 'symlink'`).

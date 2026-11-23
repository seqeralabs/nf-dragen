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
│   │   ├── <sample>.ploidy.vcf.gz(.tbi)  # ploidy estimator calls
│   │   ├── <sample>.baf.bw, .hard-filtered.baf.bw  # B-allele frequencies
│   │   ├── <sample>.roh.bed, .targeted.vcf.gz(.tbi), .targeted.json  # germline runs, see below
│   │   ├── <tumour>_vs_<normal>.*        # somatic runs: the same files as above (not .roh.bed or .targeted.*), plus
│   │   ├── <tumour>_vs_<normal>_tumor.bam(.bai)   # tumour BAM; <tumour>_vs_<normal>.bam is the normal BAM
│   │   ├── <tumour>_vs_<normal>.sv.vcf.gz(.tbi)   # structural variants
│   │   ├── <tumour>_vs_<normal>.cnv.vcf.gz(.tbi)  # copy number variants, whole-genome somatic runs only
│   │   └── <tumour>_vs_<normal>.wgs_somatic_callable_regions.bed  # whole-genome somatic runs only
│   └── rna_fastq_to_bam/              # RNA samples: <sample>.bam, .bam.bai, metrics and splice junctions, no VCFs
├── multiqc/                           # unless --skip_multiqc
│   ├── multiqc_report.html            # MultiQC report: DRAGEN metrics and software versions
│   ├── multiqc_data/                  # the report's data, as text files
│   └── multiqc_plots/                 # the report's plots, as PNG, SVG and PDF
└── pipeline_info/
    └── software_versions.yml          # DRAGEN version of every DRAGEN process, plus Nextflow and pipeline versions
```

- `genome/index/` holds a hash table only when the pipeline built it. Pass it back with `--dragen_index_dna` or `--dragen_index_rna` to skip the build in a later run. An index you supply is not copied.
- `dragen/` holds the DRAGEN outputs of each sample, named by the sample name. A somatic run writes the outputs of the tumour/normal pair under `dragen/dna_fastq_to_vcf/`, named `<tumour>_vs_<normal>`, with a BAM for each of the two samples. Each samplesheet row of a sample is one read group in the BAM: the ID is the sample name for a sample with one row, and `<sample>.1`, `<sample>.2`, ... for a sample with several rows.
- DRAGEN QC metrics (`<prefix>` is the sample name, or `<tumour>_vs_<normal>` for a somatic run):
  - every run: `<prefix>.mapping_metrics.csv` (for RNA runs it includes the RNA-specific mapping metrics), `<prefix>.time_metrics.csv`, `<prefix>.fastqc_metrics.csv`, `<prefix>.trimmer_metrics.csv`, and the whole-genome coverage reports `<prefix>.wgs_coverage_metrics.csv`, `<prefix>.wgs_overall_mean_cov.csv`, `<prefix>.wgs_contig_mean_cov.csv`, `<prefix>.wgs_fine_hist.csv` and `<prefix>.wgs_hist.csv`;
  - germline and somatic runs: `<prefix>.fragment_length_hist.csv`, `<prefix>.vc_metrics.csv`, `<prefix>.ploidy_estimation_metrics.csv` and `<prefix>.allele_transition_noise_metrics.csv`, plus the same five coverage reports for the target regions (`<prefix>.target_bed_coverage_metrics.csv`, `_overall_mean_cov`, `_contig_mean_cov`, `_fine_hist`, `_hist`) with `--target_bed`;
  - germline runs: `<prefix>.vc_hethom_ratio_metrics.csv` and `<prefix>.roh_metrics.csv`;
  - somatic runs: `<prefix>.sv_metrics.csv`, and the coverage reports once per sample, ending in `_tumor` and `_normal` (e.g. `<prefix>.wgs_coverage_metrics_tumor.csv`, `<prefix>.wgs_hist_normal.csv`);
  - whole-genome somatic runs: `<prefix>.cnv_metrics.csv`.
- Germline and somatic runs also write `<prefix>.ploidy.vcf.gz` and its index (ploidy estimator calls) and the B-allele frequencies of the small variant calls as BigWig, `<prefix>.baf.bw` and `<prefix>.hard-filtered.baf.bw` (PASS calls only).
- Germline runs also write `<prefix>.roh.bed` (runs of homozygosity) and the calls of DRAGEN's targeted callers (genes such as CYP2D6 and SMN1/2), `<prefix>.targeted.vcf.gz(.tbi)` and `<prefix>.targeted.json`.
- Whole-genome somatic runs also write `<prefix>.wgs_somatic_callable_regions.bed`, the regions with enough coverage in both samples to call somatic variants.
- RNA runs also write the splice junctions: `<prefix>.SJ.out.tab` (filtered), `<prefix>.unfiltered.SJ.out.tab` and `<prefix>.Chimeric.out.junction` (chimeric alignments).
- Each of these files is published when DRAGEN writes it. DRAGEN's logs, replay and usage files, MD5 sums, `<prefix>.metrics.json` and `<prefix>.insert-stats.tab` (the same metrics as the CSVs), and the CNV intermediates (`.tn.tsv.gz`, `tumor.*`, `.cnv.gff3`, `.cnv.igv_session.xml`) are not published.
- On Seqera Platform, the MultiQC report and these DRAGEN QC files appear in the run's Reports tab (listed in `tower.yml`).
- `multiqc/` summarises every DRAGEN `*metrics*.csv` and the software versions. With `--skip_dragen` it holds only the software versions. Set `--multiqc_title` to name the report (the HTML file name then starts with the title), `--multiqc_config` to add a MultiQC config and `--multiqc_logo` to add a logo.
- Files are copied to the output directory. To change that, set `workflow.output.mode` in your own config (for example `workflow.output.mode = 'symlink'`).

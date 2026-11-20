# seqeralabs/nf-dragen: Output

The pipeline writes these files to the directory given with `--outdir`:

```
<outdir>/
├── genome/
│   └── index/
│       ├── dragen_index_dna/          # DNA hash table, when built from --fasta
│       └── dragen_index_rna/          # RNA hash table, when built from --fasta
├── dragen/
│   ├── dna_fastq_to_vcf/              # DNA samples
│   │   ├── <sample>.bam
│   │   ├── <sample>.vcf.gz
│   │   ├── <sample>.vcf.gz.tbi
│   │   ├── <sample>.hard-filtered.vcf.gz
│   │   └── <sample>.hard-filtered.vcf.gz.tbi
│   └── rna_fastq_to_bam/              # RNA samples, same files as for DNA
└── pipeline_info/
    └── software_versions.yml          # versions of the software every process ran
```

- `genome/index/` holds a hash table only when the pipeline built it. Pass it back with `--dragen_index_dna` or `--dragen_index_rna` to skip the build in a later run. An index you supply is not copied.
- `dragen/` holds the DRAGEN outputs of each sample, named by the sample name. FASTQ files that DRAGEN writes are published here too.
- `--publish_dir_mode` sets how files are written to `--outdir` (default `copy`).

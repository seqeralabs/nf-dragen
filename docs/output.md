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
│   │   ├── <sample>.bam.bai
│   │   ├── <sample>.*metrics.csv      # DRAGEN metrics, e.g. <sample>.mapping_metrics.csv
│   │   ├── <sample>.vcf.gz
│   │   ├── <sample>.vcf.gz.tbi
│   │   ├── <sample>.hard-filtered.vcf.gz
│   │   ├── <sample>.hard-filtered.vcf.gz.tbi
│   │   ├── <tumour>_vs_<normal>.*        # somatic runs: the same files as above, plus
│   │   ├── <tumour>_vs_<normal>.sv.vcf.gz(.tbi)   # structural variants
│   │   └── <tumour>_vs_<normal>.cnv.vcf.gz(.tbi)  # copy number variants
│   └── rna_fastq_to_bam/              # RNA samples, same files as for DNA
└── pipeline_info/
    └── software_versions.yml          # versions of the software every process ran
```

- `genome/index/` holds a hash table only when the pipeline built it. Pass it back with `--dragen_index_dna` or `--dragen_index_rna` to skip the build in a later run. An index you supply is not copied.
- `dragen/` holds the DRAGEN outputs of each sample, named by the sample name. A somatic run writes the outputs of the tumour/normal pair under `dragen/dna_fastq_to_vcf/`, named `<tumour>_vs_<normal>`, with a BAM for each of the two samples. Each lane of a sample is its own read group in the BAM, with ID `<sample>.<lane>`, where `<lane>` counts the sample's lanes from 1 in samplesheet order.
- `--publish_dir_mode` sets how files are written to `--outdir` (default `copy`).

# seqeralabs/nf-dragen: Parameters

The output directory is not a parameter: set it with `-output-dir <dir>` (one dash) or `outputDir` in a config file. It is required on every run.

See [usage.md](usage.md) for how the DRAGEN args, `--target_bed` and the hash-table options behave.

## Input/output options

| Name         | Description                                                                                                  | Type     | Default | Required |
| ------------ | ------------------------------------------------------------------------------------------------------------ | -------- | ------- | :------: |
| `--input`    | Path to comma-separated file containing information about the samples in the experiment.                     | `string` | n/a     |   yes    |
| `--seq_type` | Default type of sequencing (dna, rna etc). Can be overridden in the samplesheet using the column `seq_type`. | `string` | `dna`   |    no    |

## Tools

| Name                          | Description                                                                   | Type      | Default                                                                 | Required |
| ----------------------------- | ----------------------------------------------------------------------------- | --------- | ----------------------------------------------------------------------- | :------: |
| `--dragen_germline_args`      | DRAGEN options for germline runs. Replaces the default.                       | `string`  | `--enable-duplicate-marking true`                                       |    no    |
| `--dragen_somatic_args`       | DRAGEN options for somatic (tumour/normal) runs. Replaces the default.        | `string`  | `--enable-duplicate-marking true --vc-target-vaf 0.03 --enable-sv true` |    no    |
| `--dragen_license`            | How DRAGEN is licensed.                                                       | `string`  | `marketplace`                                                           |    no    |
| `--dragen_rna_args`           | DRAGEN options for RNA runs. Replaces the default.                            | `string`  | `--enable-duplicate-marking true`                                       |    no    |
| `--dragen_hashtable_dna_args` | DRAGEN options for the DNA hash table build. Replaces the default.            | `string`  | `--ht-build-cnv-hashtable true`                                         |    no    |
| `--dragen_hashtable_rna_args` | DRAGEN options for the RNA hash table build. Replaces the default.            | `string`  | `--ht-build-rna-hashtable true`                                         |    no    |
| `--skip_dragen`               | Skip all DRAGEN steps, including building the hash tables. Only MultiQC runs. | `boolean` | `False`                                                                 |    no    |
| `--skip_multiqc`              | Skip the MultiQC report                                                       | `boolean` | `False`                                                                 |    no    |

## Reference genome options

| Name                 | Description                                                                                                                                 | Type     | Default | Required |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- | -------- | ------- | :------: |
| `--fasta`            | Reference genome FASTA file. The DRAGEN hash tables are built from it when no index is given with --dragen_index_dna or --dragen_index_rna. | `string` | n/a     |    no    |
| `--dragen_index_dna` | Path to DNA index created by Dragen.                                                                                                        | `string` | n/a     |    no    |
| `--dragen_index_rna` | Path to RNA index created by Dragen.                                                                                                        | `string` | n/a     |    no    |
| `--target_bed`       | BED file of exome or panel target regions.                                                                                                  | `string` | n/a     |    no    |

## MultiQC options

| Name               | Description                                                                                 | Type     | Default | Required |
| ------------------ | ------------------------------------------------------------------------------------------- | -------- | ------- | :------: |
| `--multiqc_config` | Custom config file to supply to MultiQC.                                                    | `string` | n/a     |    no    |
| `--multiqc_title`  | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` | n/a     |    no    |
| `--multiqc_logo`   | Custom logo file to supply to MultiQC                                                       | `string` | n/a     |    no    |

---

_This pipeline was built with [Nextflow](https://nextflow.io).
Documentation generated by [nf-docs](https://github.com/ewels/nf-docs) v0.4.0 on 2026-11-24 11:50:10 UTC._

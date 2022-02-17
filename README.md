
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)

> THIS IS A PROOF-OF-CONCEPT REPOSITORY THAT IS UNDER ACTIVE DEVELOPMENT. SYNTAX, ORGANISATION AND LAYOUT MAY CHANGE WITHOUT NOTICE!

## Introduction

**nf-core/dragen** is a simple, proof-of-concept pipeline to run the [Illumina DRAGEN](https://emea.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html) suite of tools.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker  containers making installation trivial and results highly reproducible.

DRAGEN is licensed software that you can purchase in the [AWS Marketplace](https://aws.amazon.com/marketplace/pp/prodview-ypz2tpzy6f5xq). This pipeline has only been tested using AWS Batch and doesn't come with any bundled containers.

## Credits

nf-core/dragen was originally written by Harshil Patel, Graham Wright.

## Citations

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

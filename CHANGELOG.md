# nf-core/dragen: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of nf-core/dragen, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Publish csv outputs from `DRAGEN` for `multiQC` module
- Stub runs for `DRAGEN`, including creation of a file recognised by the `multiQC` module

### `Fixed`

- Fixed error `Access to 'FASTQC.out' is undefined since the workflow 'FASTQC' has not been invoked before accessing the output attribute` when `-skip_fastqc` enabled by adjusting channel generation
- Added `arm` profile for docker runs on Mac OS

### `Dependencies`

### `Deprecated`

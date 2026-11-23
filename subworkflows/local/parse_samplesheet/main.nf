nextflow.enable.types = true

include { samplesheetToList } from 'plugin/nf-schema'

workflow PARSE_SAMPLESHEET {
    take:
    samplesheet: Path
    default_seq_type: String // for rows without a seq_type
    schema: Path             // samplesheet JSON schema

    main:
    // Validated and grouped here, so a bad samplesheet fails before any task runs
    def rows = samplesheetToList(samplesheet, schema)
    def sample_records = rows.groupBy { row -> row[0].id }.collect { id, sample_rows -> toSample(id, sample_rows, default_seq_type) }

    // One tumour per DNA patient, checked here so it fails before any task (SELECT_DRAGEN_MODE runs in an operator)
    sample_records
        .findAll { sample -> sample.seq_type == 'dna' && sample.tissue == 'tumor' }
        .groupBy { sample -> sample.patient }
        .each { patient, tumours ->
            if (tumours.size() > 1) {
                error("Patient '${patient}' has ${tumours.size()} tumour DNA samples (${tumours*.id.join(', ')}); give one tumour per patient")
            }
        }

    // DRAGEN allows each FASTQ file once in its fastq-list, in fastq_1 and fastq_2 together
    rows
        .collectMany { row -> [row[1], row[2]].findAll { fastq -> fastq }.collect { fastq -> [path: fastq.toUriString(), sample: row[0].id] } }
        .groupBy { entry -> entry.path }
        .each { path, entries ->
            if (entries.size() > 1) {
                error("Samplesheet: FASTQ '${path}' is used ${entries.size()} times (samples ${entries*.sample.unique().join(', ')}); each FASTQ file may be listed only once")
            }
        }

    emit:
    samples: Channel<Sample> = channel.fromList(sample_records)
}

// All rows (FASTQ pairs) of one sample as one record; rows must agree on the sample-level columns
def toSample(id, rows, default_seq_type) {
    def metas = rows.collect { row ->
        def meta = row[0]
        [
            patient: meta.patient ?: id,
            tissue: meta.tissue,
            sex: meta.sex,
            seq_type: meta.seq_type ?: default_seq_type,
        ]
    }
    ['patient', 'tissue', 'sex', 'seq_type'].each { column ->
        def values = metas.collect { meta -> meta[column] }.unique()
        if (values.size() > 1) {
            error("Samplesheet: rows of sample '${id}' disagree on column '${column}' (${values.join(', ')})")
        }
    }
    def fastq_2 = rows.collect { row -> row[2] }.findAll { fastq -> fastq }
    if (fastq_2 && fastq_2.size() != rows.size()) {
        error("Samplesheet: sample '${id}' mixes paired-end and single-end rows; give fastq_2 for every row or for none")
    }
    return record(
        id: id,
        patient: metas[0].patient,
        tissue: metas[0].tissue,
        sex: metas[0].sex,
        seq_type: metas[0].seq_type,
        fastq_1: rows.collect { row -> row[1] },
        fastq_2: fastq_2,
    )
}

// One sample with all its FASTQ pairs, in samplesheet order (fastq_2 is empty for single-end)
record Sample {
    id: String
    patient: String
    tissue: String
    sex: String
    seq_type: String
    fastq_1: List<Path>
    fastq_2: List<Path>
}

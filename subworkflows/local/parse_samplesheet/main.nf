nextflow.enable.types = true

include { samplesheetToList } from 'plugin/nf-schema'

workflow PARSE_SAMPLESHEET {
    take:
    samplesheet: Path
    seq_type: String  // for rows without a seq_type

    main:
    // Validated and grouped here, so a bad samplesheet fails before any task runs
    def rows = samplesheetToList(samplesheet, "${moduleDir}/../../../assets/schema_input.json")
    def sample_records = rows.groupBy { row -> row[0].id }.collect { id, sample_rows -> toSample(id, sample_rows, seq_type) }

    emit:
    samples: Channel<Sample> = channel.fromList(sample_records)
}

// All rows (lanes) of one sample as one record; rows must agree on the sample-level columns
def toSample(id, rows, default_seq_type) {
    def metas = rows.collect { row ->
        def meta = row[0]
        [
            patient: meta.patient ?: id,
            status: meta.status ?: 0,
            sex: meta.sex ?: null,
            seq_type: meta.seq_type ?: default_seq_type,
        ]
    }
    ['patient', 'status', 'sex', 'seq_type'].each { column ->
        def values = metas.collect { meta -> meta[column] }.unique()
        if (values.size() > 1) {
            error("Samplesheet: rows of sample '${id}' disagree on column '${column}' (${values.join(', ')})")
        }
    }
    def fastq_2 = rows.collect { row -> row[2] }.findAll { fastq -> fastq }
    if (fastq_2 && fastq_2.size() != rows.size()) {
        error("Samplesheet: sample '${id}' mixes paired-end and single-end rows; give fastq_2 for every row or for none")
    }
    // Without a lane column, a lane is named by its position in the sample
    def lanes = rows.withIndex().collect { row, i -> (row[0].lane ?: (i + 1)).toString() }
    return record(
        id: id,
        patient: metas[0].patient,
        status: metas[0].status,
        sex: metas[0].sex,
        seq_type: metas[0].seq_type,
        lanes: lanes,
        fastq_1: rows.collect { row -> row[1] },
        fastq_2: fastq_2,
    )
}

// One sample with all its lanes, in samplesheet order (fastq_2 is empty for single-end)
record Sample {
    id: String
    patient: String
    status: Integer
    sex: String?
    seq_type: String
    lanes: List<String>
    fastq_1: List<Path>
    fastq_2: List<Path>
}

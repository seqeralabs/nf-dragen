nextflow.enable.types = true

include { samplesheetToList } from 'plugin/nf-schema'

workflow PARSE_SAMPLESHEET {
    take:
    samplesheet: Path
    seq_type: String  // for rows without a seq_type

    main:
    // Validated here (nf-schema also rejects a sample name on two rows), so a bad samplesheet fails before any task runs
    def rows = samplesheetToList(samplesheet, "${moduleDir}/../../../assets/schema_input.json")
    def sample_records = rows.collect { meta, fastq_1, fastq_2 ->
        record(
            id: meta.id,
            patient: meta.patient ?: meta.id,
            status: meta.status ?: 0,
            sex: meta.sex ?: null,
            seq_type: meta.seq_type ?: seq_type,
            fastq_1: fastq_1,
            fastq_2: fastq_2 ?: null,
        )
    }

    emit:
    samples: Channel<Sample> = channel.fromList(sample_records)
}

// One samplesheet row (fastq_2 is null for single-end)
record Sample {
    id: String
    patient: String
    status: Integer
    sex: String?
    seq_type: String
    fastq_1: Path
    fastq_2: Path?
}

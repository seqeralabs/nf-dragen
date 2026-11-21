nextflow.enable.types = true

workflow SELECT_DRAGEN_MODE {
    take:
    samples: Channel<ModeSample>

    main:
    // Pairing needs every DNA sample of a patient, so it waits for all samples (they come from the samplesheet at once)
    ch_dna_runs = samples
        .filter { sample -> sample.seq_type == 'dna' }
        .collect()
        .map { dna -> dnaRuns(dna) }

    emit:
    germline = ch_dna_runs.flatMap { runs -> runs.germline }
    somatic  = ch_dna_runs.flatMap { runs -> runs.somatic }
    rna      = samples.filter { sample -> sample.seq_type == 'rna' }
}

// Per patient: a normal alone gets a germline run; a normal with a tumour gets one somatic run (PARSE_SAMPLESHEET rejects two or more tumours)
def dnaRuns(samples) {
    def germline = []
    def somatic = []
    samples
        .groupBy { sample -> sample.patient }
        .each { patient, patient_samples ->
            def normals = patient_samples.findAll { sample -> sample.status == 0 }
            def tumours = patient_samples.findAll { sample -> sample.status == 1 }
            if (normals.size() > 1) {
                error("Patient '${patient}' has ${normals.size()} normal DNA samples (${normals*.id.join(', ')}); give one normal per patient")
            }
            if (!normals) {
                error("Patient '${patient}' has tumour sample(s) ${tumours*.id.join(', ')} but no normal; tumour-only calling is not supported")
            }
            def normal = normals[0]
            if (!tumours) {
                germline << normal
            }
            else {
                def tumour = tumours[0]
                somatic << record(
                    id: normal.id,
                    seq_type: normal.seq_type,
                    fastq_1: normal.fastq_1,
                    fastq_2: normal.fastq_2,
                    tumour_id: tumour.id,
                    tumour_fastq_1: tumour.fastq_1,
                    tumour_fastq_2: tumour.fastq_2,
                )
            }
        }
    return [germline: germline, somatic: somatic]
}

// The Sample fields that pick the mode; records keep their other fields
record ModeSample {
    id: String
    patient: String
    status: Integer
    seq_type: String
}

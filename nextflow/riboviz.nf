#!/usr/bin/env nextflow

rrna_fasta_file = Channel.fromPath(params.rrna_fasta_file,
                                   checkIfExists: true)

process buildIndicesRrna {
    input:
        file fasta_file from rrna_fasta_file
    output:
        file "${params.rrna_index_prefix}.*.ht2*" into rrna_index_files
    when:
        params.build_indices
    shell:
        """
        hisat2-build --version
        hisat2-build ${fasta_file} ${params.rrna_index_prefix}
        """
}

// rrna_index_files.subscribe { println "RRNA Index files: ${it}" }

orf_fasta_file = Channel.fromPath(params.orf_fasta_file,
                                  checkIfExists: true)

process buildIndicesOrf {
    input:
        file fasta_file from orf_fasta_file
    output:
        file "${params.orf_index_prefix}.*.ht2*" into orf_index_files
    when:
        params.build_indices
    shell:
        """
        hisat2-build --version
        hisat2-build ${fasta_file} ${params.orf_index_prefix}
        """
}

// orf_index_files.subscribe { println "ORF index files: ${it}" }

// Alternative one task approach.
alt_rrna_fasta_file = Channel.fromPath(params.rrna_fasta_file,
                                       checkIfExists: true)
alt_rrna_index_prefix = Channel.of(params.rrna_index_prefix)
alt_orf_fasta_file = Channel.fromPath(params.orf_fasta_file,
                                      checkIfExists: true)
alt_orf_index_prefix = Channel.of(params.orf_index_prefix)
index_fasta_files = alt_rrna_fasta_file.concat(alt_orf_fasta_file)
index_prefixes = alt_rrna_index_prefix.concat(alt_orf_index_prefix)

process buildIndices {
    input:
        file index_fasta_file from index_fasta_files
        val index_prefix from index_prefixes
    output:
        file "${index_prefix}.*.ht2*" into index_files
    when:
        params.build_indices
    shell:
        """
        hisat2-build --version
        hisat2-build ${index_fasta_file} ${index_prefix}
        """
}

// index_files.subscribe { println "Index files: ${it}" }

/*
 * cutadapt using concurrent sample ID and filename channels.
 */

sample_ids = Channel.fromList(params.fq_files.keySet())
sample_files = Channel.fromPath(
    params.fq_files.values().collect({"${params.dir_in}/${it}"}))

process cutAdapters {
    input:
        val adapters from params.adapters
        val sample_id from sample_ids
        file sample_file from sample_files
    output:
        file "${sample_id}_trim.fq" into trimmed_sample_fastq
    shell:
        // TODO configure -j 0 in a more Nextflow-esque way.
        """
        cutadapt --trim-n -O 1 -m 5 -a ${adapters} -o ${sample_id}_trim.fq ${sample_file} -j 0
        """
}

/*
 * Alternative cutadapt using channel of (sample ID, filename) tuples.
 */

sample_list = params.fq_files.collect(
    {key, value -> [key, file("${params.dir_in}/${value}")]})

process cutAdaptersTuple {
    input:
        val adapters from params.adapters
        tuple val(sample_id), file(sample_file) from sample_list
    output:
        file "${sample_id}_trim.fq" into trimmed_sample_fastq_tuple
    shell:
        // TODO configure -j 0 in a more Nextflow-esque way.
        """
        cutadapt --trim-n -O 1 -m 5 -a ${adapters} -o ${sample_id}_trim.fq ${sample_file} -j 0
        """
}

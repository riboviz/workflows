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
fq_files: # fastq files to be processed, relative to dir_in
  WTnone: SRR1042855_s1mi.fastq.gz
  WT3AT: SRR1042864_s1mi.fastq.gz
  NotHere: example_missing_file.fastq.gz # Test case for missing file
*/

/*
fq_file = 'SRR1042855_s1mi.fastq.gz'
sample = Channel.fromPath("${params.dir_in}/${fq_file}")
process cutAdapters {
    input:
        file sample from sample
        val adapters from params.adapters
    output:
        file 'trim.fq' into trimmed_fastq
    shell:
        """
        cutadapt --trim-n -O 1 -m 5 -a ${adapters} -o trim.fq ${sample} -j 0
        """
}
*/

/**
params.fq_files.each({key, value -> println "Key-value: $key $value"})
samples = Channel.of(params.fq_files.each({key, value -> [key, value]}))
samples.subscribe { println "Sample: ${it}" }
sample_ids = params.fq_files.each({key, value -> key}).collect()
println "${sample_ids}\n"
sample_files = params.fq_files.each({key, value -> value}).collect()
println "${sample_files}\n"
xxx = params.fq_files.entrySet()
yyy = xxx.each { entry -> entry.key }
println "YYY ${yyy}\n"
*/


sample_ids = Channel.fromList(params.fq_files.keySet())
sample_files = Channel.fromPath(params.fq_files.values().collect({"${params.dir_in}/${it}"}))

sample_list = params.fq_files.collect({key, value -> [key, file(value)]})

process cutAdapters {
    input:
        val sample_id from sample_ids
        file sample_file from sample_files

        // 1
        // cutadapt --trim-n -O 1 -m 5 -a CTGTAGGCACC -o WT3AT_trim.fq SRR1042864_s1mi.fastq.gz -j 0
        // OSError: pigz: skipping: SRR1042864_s1mi.fastq.gz does not exist as file has not been staged.
        // tuple val(sample_id), sample_file from sample_list
	// 2
        // cutadapt --trim-n -O 1 -m 5 -a CTGTAGGCACC -o WTnone_trim.fq input.1 -j 0
        // tuple val(sample_id), file(sample_file) from sample_list

        val adapters from params.adapters
    output:
        file "${sample_id}_trim.fq" into trimmed_sample_fastq
    shell:
        // TODO configure -j 0 in a more Nextflow-esque way.
        """
        echo ${sample_id}
        echo ${sample_file}
        cutadapt --trim-n -O 1 -m 5 -a ${adapters} -o ${sample_id}_trim.fq ${sample_file} -j 0
        """
}

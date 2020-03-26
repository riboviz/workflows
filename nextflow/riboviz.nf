#!/usr/bin/env nextflow

rrna_fasta_file = Channel.fromPath(params.rrna_fasta_file,
                                   checkIfExists: true)

process buildIndicesrRNA {
    input:
        file fasta_file from rrna_fasta_file
    output:
        file "${params.rrna_index_prefix}.*.ht2" into rrna_index_files
    when:
        params.build_indices
    shell:
        """
        hisat2-build --version
        hisat2-build ${fasta_file} ${params.rrna_index_prefix}
        """
}

orf_fasta_file = Channel.fromPath(params.orf_fasta_file,
                                  checkIfExists: true)

process buildIndicesORF {
    input:
        file fasta_file from orf_fasta_file
    output:
        file "${params.orf_index_prefix}.*.ht2" into orf_index_files
    when:
        params.build_indices
    shell:
        """
        hisat2-build --version
        hisat2-build ${fasta_file} ${params.orf_index_prefix}
        """
}

// Filter samples down to those whose files exist.
sample_list = []
for (entry in params.fq_files) {
    sample_file = file("${params.dir_in}/${entry.value}")
    if (sample_file.exists()) {
        sample_list.add([entry.key, sample_file])
    } else {
        println "Missing file ($entry.key): $entry.value"
    }
}

process cutAdapters {
    input:
        tuple val(sample_id), file(sample_file) from sample_list
        val adapters from params.adapters
    output:
        tuple val("${sample_id}"), file("${sample_id}_trim.fq") into trim_sample_fastq_tuple
    shell:
        // TODO configure -j 0 in a more Nextflow-esque way.
        """
        cutadapt --trim-n -O 1 -m 5 -a ${adapters} -o ${sample_id}_trim.fq ${sample_file} -j 0
        """
}

process hisat2rRNA {
    input:
        tuple val(sample_id), file(trim_fq) from trim_sample_fastq_tuple
        each file("${params.rrna_index_prefix}.*.ht2") from rrna_index_files
    output:
        tuple val("${sample_id}"), file("${sample_id}_nonrRNA.fq") into non_rrna_fq
        tuple val("${sample_id}"), file("${sample_id}_rRNA_map.sam") into rrna_map_sam
    shell:
        """
        hisat2 --version
        hisat2 -p ${params.num_processes} -N 1 -k 1 --un ${sample_id}_nonrRNA.fq -x ${params.rrna_index_prefix} -S ${sample_id}_rRNA_map.sam -U ${trim_fq}
        """
}

process hisat2ORF {
    input:
        tuple val(sample_id), file(non_rrna_fq) from non_rrna_fq
        each file("${params.orf_index_prefix}.*.ht2") from orf_index_files
    output:
        tuple val("${sample_id}"), file("${sample_id}_unaligned.fq") into unaligned_fq
        tuple val("${sample_id}"), file("${sample_id}_orf_map.sam") into orf_map_sam
    shell:
        """
        hisat2 --version
        hisat2 -p ${params.num_processes} -k 2 --no-spliced-alignment --rna-strandness F --no-unal --un ${sample_id}_unaligned.fq -x ${params.orf_index_prefix} -S ${sample_id}_orf_map.sam -U ${non_rrna_fq}
        """
}

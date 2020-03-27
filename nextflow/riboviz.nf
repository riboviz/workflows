#!/usr/bin/env nextflow

rrna_fasta = Channel.fromPath(params.rrna_fasta_file,
                              checkIfExists: true)

process buildIndicesrRNA {
    tag "${params.rrna_index_prefix}"
    publishDir "${params.dir_index}"
    input:
        file fasta from rrna_fasta
    output:
        file "${params.rrna_index_prefix}.*.ht2" into rrna_indices
    when:
        params.build_indices
    shell:
        """
        hisat2-build --version
        hisat2-build ${fasta} ${params.rrna_index_prefix}
        """
}

orf_fasta = Channel.fromPath(params.orf_fasta_file,
                             checkIfExists: true)

process buildIndicesORF {
    tag "${params.orf_index_prefix}"
    publishDir "${params.dir_index}"
    input:
        file fasta from orf_fasta
    output:
        file "${params.orf_index_prefix}.*.ht2" into orf_indices
    when:
        params.build_indices
    shell:
        """
        hisat2-build --version
        hisat2-build ${fasta} ${params.orf_index_prefix}
        """
}

// Filter samples down to those whose files exist.
samples = []
for (entry in params.fq_files) {
    sample_file = file("${params.dir_in}/${entry.value}")
    if (sample_file.exists()) {
        samples.add([entry.key, sample_file])
    } else {
        println "Missing file ($entry.key): $entry.value"
    }
}

process cutAdapters {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}"
    input:
        tuple val(sample_id), file(sample_file) from samples
        val adapters from params.adapters
    output:
        tuple val(sample_id), file("${sample_id}_trim.fq") into cut_adapters
    shell:
        // TODO configure -j 0 in a more Nextflow-esque way.
        """
        echo Trimming ${sample_id}
        cutadapt --trim-n -O 1 -m 5 -a ${adapters} -o ${sample_id}_trim.fq ${sample_file} -j 0
        """
}

process hisat2rRNA {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}"
    input:
        tuple val(sample_id), file(fastq) from cut_adapters
        each file(indices) from rrna_indices
    output:
        tuple val(sample_id), file("${sample_id}_nonrRNA.fq") into non_rrnas
        tuple val(sample_id), file("${sample_id}_rRNA_map.sam") into rrna_maps
    shell:
        """
        hisat2 --version
        hisat2 -p ${params.num_processes} -N 1 -k 1 --un ${sample_id}_nonrRNA.fq -x ${params.rrna_index_prefix} -S ${sample_id}_rRNA_map.sam -U ${fastq}
        """
}

process hisat2ORF {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}"
    input:
        tuple val(sample_id), file(fastq) from non_rrnas
        each file(indices) from orf_indices
    output:
        tuple val(sample_id), file("${sample_id}_unaligned.fq") into unaligneds
        tuple val(sample_id), file("${sample_id}_orf_map.sam") into orf_maps
    shell:
        """
        hisat2 --version
        hisat2 -p ${params.num_processes} -k 2 --no-spliced-alignment --rna-strandness F --no-unal --un ${sample_id}_unaligned.fq -x ${params.orf_index_prefix} -S ${sample_id}_orf_map.sam -U ${fastq}
        """
}

process trim5pMismatches {
    tag "${sample_id}"
    publishDir "${params.dir_tmp}/${sample_id}"
    input:
        tuple val(sample_id), file(sam) from orf_maps
    output:
        tuple val(sample_id), file("${sample_id}_orf_map_clean.sam") into clean_orf_maps
        tuple val(sample_id), file("${sample_id}_trim_5p_mismatch.tsv") into trim_summaries
    shell:
        """
        python -m riboviz.tools.trim_5p_mismatch -m 2 -i ${sam} -o ${sample_id}_orf_map_clean.sam -s ${sample_id}_trim_5p_mismatch.tsv
        """
}

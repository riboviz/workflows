#!/usr/bin/env nextflow

process buildIndicesRrna {
    input:
        file fasta_file from file(params.rrna_fasta_file)
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

process buildIndicesOrf {
    input:
        file fasta_file from file(params.orf_fasta_file)
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

// Alternative one task approach.
alt_rrna_fasta_file = Channel.fromPath(params.rrna_fasta_file)
alt_rrna_index_prefix = Channel.of(params.rrna_index_prefix)
alt_orf_fasta_file = Channel.fromPath(params.orf_fasta_file)
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

rrna_index_files.subscribe { println "RRNA Index files: ${it}" }
orf_index_files.subscribe { println "ORF index files: ${it}" }
index_files.subscribe { println "Index files: ${it}" }

input_dir = Channel.fromPath(params.dir_in)
adapters = params.adapters

/*
lambda wildcards: os.path.join(config['dir_in'], SAMPLES[wildcards.sample])
os.path.join(config['dir_tmp'], "{sample}", "trim.fq")
log: os.path.join(config['dir_logs'], TIMESTAMP, "{sample}", "cutadapt.log")
cutadapt --trim-n -O 1 -m 5 -a	{config[adapters]} -o {output} {input} -j 0 &> {log}
*/
/*
process cutAdapters {
    input:
    file sequence from TODO
    val adapters from adapters
    output:
    shell:
    // TODO configure -j 0 in a more Nextflow-esque way.
    """
    cutadapt --trim-n -O 1 -m 5 -a ${adapters} -o ${trimmed_sequence} ${sequence} -j 0
    """
}
*/


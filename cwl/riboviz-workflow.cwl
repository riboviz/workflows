#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {} # Allow use of "valueFrom" in "step".
inputs:
  rrna_fasta_file: File
  rrna_index_prefix: string
  orf_fasta_file: File
  orf_index_prefix: string
  sample: string
  sample_file: File
  adapter: string
outputs:
  indexed_rrna_files:
    type:
      type: array
      items: File
    outputSource: hisat2_build_rrna/index_files
  indexed_orf_files:
    type:
      type: array
      items: File
    outputSource: hisat2_build_orf/index_files
  trim_file:
    type: File
    outputSource: cutadapt/output_file
steps:
  hisat2_build_rrna:
    run: hisat2_index.cwl
    in:
      reference_fasta: rrna_fasta_file
      index_basename: rrna_index_prefix
    out: [index_files]
  hisat2_build_orf:
    run: hisat2_index.cwl
    in:
      reference_fasta: orf_fasta_file
      index_basename: orf_index_prefix
    out: [index_files]
  cutadapt:
    run: cutadapt.cwl
    in:
      trim:
        default: true
      overlap:
        default: 1
      min_length:
        default: 5
      adapter: adapter
      cores:
        default: 0
      input: sample_file
      output:
        source: sample # Get value of "sample" for use below.
        valueFrom: $(self + '.trim.fq')
    out: [output_file]

cwlVersion: v1.0
class: CommandLineTool
label: "cutaddapt"
doc: "...cutadapt...TODO...add URL..."
requirements:
  InlineJavascriptRequirement: {}
baseCommand: [cutadapt]
inputs:
  trim:
    type: boolean
    inputBinding:
      position: 1
      prefix: --trim-n
  overlap:
    type: int
    inputBinding:
      position: 2
      prefix: -O
  min_length:
    type: int
    inputBinding:
      position: 3
      prefix: -m
  adapter:
    type: string
    inputBinding:
      position: 4
      prefix: -a
  output:
     type: string
     inputBinding:
       position: 5
       prefix: -o
  cores:
    type: int
    inputBinding:
      position: 6
      prefix: -j
  input:
    type: File
    inputBinding:
      position: 7
outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output)

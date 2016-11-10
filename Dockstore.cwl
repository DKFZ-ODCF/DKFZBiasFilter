#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "DKFZBiasFilter"
label: "DKFZ Bias Filter"
cwlVersion: v1.0 
description: |
    A Docker container for the DKFZ Bias Filter.
    ```
    Usage:
    # fetch CWL
    $> dockstore tool cwl --entry quay.io/jwerner_dkfz/dkfzbiasfilter:1.0 > Dockstore.cwl
    # make a runtime JSON template and edit it (or use the content of sample_configs.json in this git repo)
    $> dockstore tool convert cwl2json --cwl Dockstore.cwl > Dockstore.json
    # run it locally with the Dockstore CLI
    $> dockstore tool launch --entry quay.io/jwerner_dkfz/dkfzbiasfilter:1.0 \
        --json Dockstore.json
    ```

dct:creator:
  foaf:name: Ivo Buchhalter
  foaf:mbox: "mailto:i.buchhalter@dkfz-heidelberg.de"

requirements:
  - class: DockerRequirement
    # dockerPull: "quay.io/jwerner_dkfz/dkfzbiasfilter:1.0"
    dockerPull: "wernerjo/bias"
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.input_vcf)
      - $(inputs.input_bam)
      - $(inputs.input_bam_index)
      - $(inputs.reference_sequence)
      - $(inputs.reference_sequence_index)

inputs:
  write_qc:
    type: boolean
    default: true
    doc: "Write quality control? If true, then a folder is created within the same folder as the filtered vcf file storing bias plots and qc files"
    inputBinding:
      position: 1
      prefix: -q
  mapq:
    type: int
    default: 1
    doc: "Minimal mapping quality of a read to be considered for error count calculation"
    inputBinding:
      position: 2
      prefix: --mapq=
      separate: false
  baseq:
    type: int
    default: 1
    doc: "Minimal base quality to be considered for error count calculation"
    inputBinding:
      position: 3
      prefix: --baseq=
      separate: false
  temp_folder:
    type: string
    default: "/var/spool/cwl/"
    doc: "Path to the folder where temporary files are stored"
    inputBinding:
      position: 4
      prefix: --tempFolder=
      separate: false
  input_vcf:
    type: File
    default: "/home/pcawg/input.vcf"
    doc: "Absolute filename of input vcf file"
    inputBinding:
      position: 5
      valueFrom: $(self.basename)
  input_bam:
    type: File
    default: "/home/pcawg/tumor.bam"
    doc: "Absolute filename of tumor bam file"
    inputBinding:
      position: 6
      valueFrom: $(self.basename)
  input_bam_index:
    type: File
    default: "/home/pcawg/tumor.bam.bai"
    doc: "Absolute filename of tumor bam file index"
  reference_sequence:
    type: File
    default: "/home/pcawg/hs37d5.fa"
    doc: "Absolute filename of reference sequence file"
    inputBinding:
      position: 7
      valueFrom: $(self.basename)
  reference_sequence_index:
    type: File
    default: "/home/pcawg/hs37d5.fa.fai"
    doc: "Absolute filename of reference sequence file index"
  output_vcf_filename:
    type: string
    default: "/var/spool/cwl/filtered.vcf"
    doc: "Absolute filename of filtered vcf file"
    inputBinding:
      position: 8

outputs:
  output_vcf_file:
    type: File
    outputBinding:
      glob: filtered.vcf
    doc: "The filtered vcf file"

baseCommand: ["python", "/usr/local/bin/biasFilter.py"]

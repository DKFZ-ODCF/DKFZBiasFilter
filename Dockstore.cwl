#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "DKFZBiasFilter"
label: "DKFZ Bias Filter"
cwlVersion: v1.0

doc: |
    A Docker container for the DKFZ Bias Filter.
    
    Usage:
    
    Clone the git repository and build the docker.
    ```
    git clone https://github.com/eilslabs/DKFZBiasFilter.git
    cd DKFZBiasFilter/
    docker build -t "DKFZBiasFilter" .
    ```
    
    Run the docker:
    ```
    docker run \
        -v /path/to/tumor.bam:/home/pcawg/tumor.bam \
        -v /path/to/tumor.bam.bai:/home/pcawg/tumor.bam.bai \
        -v /path/to/reference.fa:/home/pcawg/hs37d5.fa \
        -v /path/to/somatic.vcf:/home/pcawg/input.vcf \
        -v /path/to/results_directory/:/home/pcawg/results \
        DKFZBiasFilter
    ```

description: |
    A Docker container for the DKFZ Bias Filter.

dct:creator:
  foaf:name: Ivo Buchhalter
  foaf:mbox: "mailto:i.buchhalter@dkfz-heidelberg.de"

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/jwerner_dkfz/dkfzbiasfilter:1.2.2"
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
  input_vcf:
    type: File
    default: "/home/pcawg/input.vcf"
    doc: "Absolute filename of input vcf file"
    inputBinding:
      position: 2
      valueFrom: $(self.basename)
  input_bam:
    type: File
    default: "/home/pcawg/tumor.bam"
    doc: "Absolute filename of tumor bam file"
    inputBinding:
      position: 3
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
      position: 4
      valueFrom: $(self.basename)
  reference_sequence_index:
    type: File
    default: "/home/pcawg/hs37d5.fa.fai"
    doc: "Absolute filename of reference sequence file index"

outputs:
  output_vcf_file:
    type: File
    outputBinding:
      glob: filtered.vcf
    doc: "The filtered vcf file"
  output_qc_folder:
    type: Directory
    outputBinding:
      glob: filtered_qcSummary
    doc: "The qc folder"

baseCommand: ["/usr/local/bin/run_biasfilter.sh"]

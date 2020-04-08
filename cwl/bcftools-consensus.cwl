#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.1
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.vcf)
      - $(inputs.vcf_index)
baseCommand: bcftools
arguments:
  - consensus
  - -i'QUAL>1 && (GT="AA" || GT="Aa")'
  - -Hla
  - -f
  - $(inputs.ref_fasta)
  - $(inputs.vcf)
inputs:
  - id: ref_fasta
    type: File
  - id: vcf
    type: File
  - id: vcf_index
    type: File
outputs:
  - id: out_fasta
    type: stdout
stdout: out.fasta

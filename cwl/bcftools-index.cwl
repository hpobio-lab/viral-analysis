#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.1
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
baseCommand: bcftools
arguments:
  - index
  - $(inputs.bcf)
inputs:
  - id: bcf
    type: File
outputs:
  - id: index
    type: File
    outputBinding:
      glob: "$(inputs.bcf).csi"

# Reference:
#  https://github.com/VGP/vgp-assembly/blob/33cd6236a68a1aee5f282e365dfe6b97e0b4ebb7/pipeline/freebayes-polish/freebayes.sh
#  https://github.com/VGP/vgp-assembly/blob/33cd6236a68a1aee5f282e365dfe6b97e0b4ebb7/pipeline/freebayes-polish/consensus.sh
class: Workflow
cwlVersion: v1.1
id: bam2fasta
label: bam2fasta
requirements: []

inputs:
  bam:
    type: File
  fasta:
    type: File
  threads:
    type: int
    default: 4
  freebayes_max_coverage:
    type: int
    default: 200

outputs:
  out_fasta:
    type: File
    outputSource: bcftools_consensus/out_fasta

steps:
  freebayes:
    in:
      bam: bam
      ref_fasta: fasta
      max_coverage: freebayes_max_coverage
    out: [vcf]
    run: freebayes.cwl
  bcftools_view_exclude_ref:
    in:
      vcf: freebayes/vcf
      threads: threads
    out: [bcf]
    run: bcftools-view-exclude-ref.cwl
  bcftools_norm:
    in:
      ref_fasta: fasta
      bcf: bcftools_view_exclude_ref/bcf
      threads: threads
    out: [normalized_bcf]
    run: bcftools-norm.cwl
  bcftools_index_after_normalization:
    in:
      bcf: bcftools_norm/normalized_bcf
    out: [index]
    run: bcftools-index.cwl
  bcftools_view_qc:
    in:
      bcf: bcftools_norm/normalized_bcf
      bcf_index: bcftools_index_after_normalization/index
      threads: threads
    out: [vcf]
    run: bcftools-view-qc.cwl
  bcftools_index_after_qc:
    in:
      bcf: bcftools_view_qc/vcf
    out: [index]
    run: bcftools-index.cwl
  bcftools_consensus:
    in:
      ref_fasta: fasta
      vcf: bcftools_view_qc/vcf
      vcf_index: bcftools_index_after_qc/index
    out: [out_fasta]
    run: bcftools-consensus.cwl

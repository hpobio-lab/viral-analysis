cwlVersion: v1.1
class: Workflow
requirements:
  SubworkflowFeatureRequirement: {}

inputs:
  ref_fasta:
    type: File
  fastq_forward:
    type: File
  fastq_reverse:
    type: File?
  threads:
    type: int
    default: 4
  freebayes_max_coverage:
    type: int
    default: 200

outputs:
  out_fasta:
    type: File
    outputSource: bam2fasta/out_fasta

steps:
  bwa-index:
    in:
      input_fasta: ref_fasta
    out: [amb, ann, bwt, pac, sa]
    run: bwa-index.cwl
  bwa-mem:
    in:
      threads: threads
      fastq_forward: fastq_forward
      fastq_reverse: fastq_reverse
      index_base: ref_fasta
      amb: bwa-index/amb
      ann: bwa-index/ann
      bwt: bwa-index/bwt
      pac: bwa-index/pac
      sa: bwa-index/sa
    out: [output]
    run: bwa-mem.cwl
  samtools-view:
    in:
      threads: threads
      input_file: bwa-mem/output
    out: [bam]
    run: samtools-view.cwl
  samtools-sort:
    in:
      input_bamfile: samtools-view/bam
      threads: threads
    out: [sorted_bam]
    run: samtools-sort.cwl
  bam2fasta:
    in:
      bam: samtools-sort/sorted_bam
      fasta: ref_fasta
      threads: threads
      freebayes_max_coverage: freebayes_max_coverage
    out: [out_fasta]
    run: bam2fasta.cwl

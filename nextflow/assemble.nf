#!/usr/bin/env nextflow

/*

minia -kmer-size 51 -abundance-min 20 -in sample1.fastq.gz -out sample1.k51.a20
cat sample1.k51.a20.contigs.fa ../NC_045512.fa >seqs.fa
minimap2 -c -x asm20 ../NC_045512.fa sample1.k51.a20.contigs.fa >seqs.paf
seqwish -s seqs.fa -p seqs.paf -g sample1_vs_ref.gfa
vg convert -x -g sample1_vs_ref.gfa >sample1_vs_ref.xg
vg deconstruct -p NC_045512 sample1_vs_ref.xg >sample1_vs_ref.vcf
*/

params.dir = "${baseDir}/../data"
params.reference = "${baseDir}/../ref/SARS-CoV2-NC_045512.2.fasta"

reference = file("${params.reference}")

fastqFiles = "${params.dir}/**.fastq.gz"
fastqs = Channel.fromPath(fastqFiles).map { path -> tuple(path.simpleName, path) }

process assembly {
  tag { sample }
  publishDir "$sample", mode: 'copy'
  container "quay.io/biocontainers/minia:3.2.1--he513fc3_0"

  input:
    set sample, file(fastq) from fastqs
  output:
    set sample, file("${sample}.fasta") into assemblies

  """
  minia -kmer-size 51 -abundance-min 20 -in $fastq -out ${sample}.k51.a20
  cat ${sample}.k51.a20.contigs.fa > ${sample}.fasta
  """
}

process overlapReads {
  tag { sample }
  publishDir "$sample", mode: 'copy'
  container "quay.io/biocontainers/minimap2:2.17--h8b12597_1"

  input:
    set sample, file(fasta) from assemblies
    file reference from reference
  output:
    set sample, file(fasta), file("${sample}.paf") into alignments

  """
  minimap2 -c -x asm20 $reference $fasta > ${sample}.paf
  """
}

process induceGraph {
  tag { sample }
  publishDir "$sample", mode: 'copy'
  container "quay.io/biocontainers/seqwish:0.4.1--h8b12597_0"

  input:
    set sample, file(fasta), file(alignment) from alignments
    file reference from reference
  output:
    set sample, file("${sample}.gfa") into graphs

  """
  cat $fasta $reference > ${sample}.withRef.fasta
  seqwish -s ${sample}.withRef.fasta -p $alignment -g ${sample}.gfa
  """
}

process callVariants {
  tag { sample }
  publishDir "$sample", mode: 'copy'
  container "heuermh/vg:latest"

  input:
    set sample, file(graph) from graphs
  output:
    set sample, file("${sample}.xg"), file("${sample}.vcf") into variants

  """
  vg convert -x -g $graph > ${sample}.xg
  vg deconstruct -p NC_045512.2 ${sample}.xg > ${sample}.vcf
  """
}

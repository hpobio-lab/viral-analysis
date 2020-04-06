#!/usr/bin/env nextflow

params.dir = "${baseDir}/data"

fastaFiles = "${params.dir}/**.fasta"
sequences = Channel.fromPath(fastaFiles).map { path -> tuple(path.simpleName, path) }

process minimap2 {
  tag { sample }
  container "quay.io/biocontainers/minimap2:2.17--h8b12597_1"

  input:
    set sample, file(fasta) from sequences
  output:
    set sample, file(fasta), file("${sample}.paf") into alignments

  """
  minimap2 -cx asm20 -w 1 -t ${task.cpus} $fasta $fasta > ${sample}.paf
  """
}

process seqwish {
  tag { sample }
  container "quay.io/biocontainers/seqwish:0.2.1--h8b12597_0"

  input:
    set sample, file(fasta), file(alignment) from alignments
  output:
    set sample, file("${sample}.gfa") into graphs

  """
  seqwish -t ${task.cpus} -k 16 -s $fasta -p $alignment -g ${sample}.gfa
  """
}

process odgiBuild {
  tag { sample }
  container "quay.io/biocontainers/odgi:v0.3--py37h8b12597_0"

  input:
    set sample, file(graph) from graphs
  output:
    set sample, file("${sample}.odgi") into sortedGraphs

  """
  odgi build -g $graph -s -o - |\
     odgi sort -i - -p s -o ${sample}.odgi
  """
}

process odgiViz {
  tag { sample }
  container "quay.io/biocontainers/odgi:v0.3--py37h8b12597_0"

  input:
    set sample, file(sortedGraph) from sortedGraphs
  output:
    set sample, file("${sample}.png") into visualizations

  """
  odgi viz -i $sortedGraph -o ${sample}.png -x 50000 -y 500 -R -P 4 -R
  """
}


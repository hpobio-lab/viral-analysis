version 1.0

task minimap2{
  input {
    File readsFA
    Int diskGB
    String outbase = basename(basename(basename(readsFA, ".gz"), ".fasta"), ".fa")
  }
  command{
    minimap2 -cx asm20 -X ${readsFA} ${readsFA} > ${outbase}.paf
  }
  runtime{
    docker : "hpobiolab/minimap2"
    memory : "12GB"
    cpu : 4
    disks : "local-disk " + diskGB + " HDD"
    preemptible : 4
  }
  output{
    File readsPAF = "${outbase}.paf"
  }
}

task seqwish{
  input {
    File readsFA
    File readsPAF
    Int diskGB
    String outbase = basename(readsPAF, ".paf")
  }
  command {
    seqiwsh -s ${readsFA} -p ${readsPAF} -g ${outbase}.gfa
  }
  runtime {
    docker : "hpobiolab/seqwish"
    memory : "12GB"
    cpu : 4
    preemptible : 4
  }
  output{
    File seqwishGFA = "${outbase}.gfa"
  }
  
}

task odgiBuild{
  input {
    File inputGFA
    Int diskGB
    String outbase = basename(inputGFA, ".gfa")
  }
  command {
    odgi build -g ${inputGFA} -s -o ${outbase}.odgi
  }
  runtime {
    docker : "hpobiolab/odgi"
    memory : "32GB"
    cpu : 4
    preemptible : 3
  }
  output {
    File odgiGraph = "${outbase}.odgi"
  }
}

task odgiViz{
  input{
    File inputODGI
    Int diskGB
    String outbase = basename(inputODGI, ".odgi")
  }
  command {
    odgi viz -i ${inputODGI} -o ${outbase}.png -x 4000 -y 500 -R -P 5
  }
  runtime {
    docker : "hpobiolab/odgi"
    memory : "10GB"
    cpu : 2
    preemptible : 2
  }
  output{
    File odgiPNG = "${outbase}.png"
  }

}

workflow PangenomeGenerate{
  input {
    File inputReads
  }

  Int mapGB = ceil(size(inputReads, "GB") + 20)
  call minimap2 as overlapReads{
    input:
      readsFA=inputReads,
      diskGB=mapGB
  }

  Int induceGB = ceil(size(inputReads, "GB") + size(overlapReads.readsPAF, "GB") * 2)
  call seqwish as induceGraph{
    input:
      readsFA=inputReads,
      readsPAF=overlapReads.readsPAF
  }

  Int buildGB = ceil(size(induceGraph.seqwishGFA, "GB") * 2)
  call odgiBuild as buildGraph{
    input:
      inputGFA=induceGraph.seqwishGFA,
      diskGB=buildGB

  }

  Int vizGB = ceil(size(buildGraph.odgiGraph, "GB") + 20)
  call odgiViz as vizGraph{
    input:
      inputODGI=buildGraph.odgiGraph,
      diskGB=vizGB
  }
}


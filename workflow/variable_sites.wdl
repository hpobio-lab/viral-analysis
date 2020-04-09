version 1.0

task rkmhFilter{
  input{
    File refFA
    File readsFA
    Int diskGB
    Int? kmerSize = 16
    Int? sketchSize = 10000
    Int? readBatch = 4000
    Int? refBatch = 5
    Int? minLength = 50
    Int? threshold = 25
    #File? readPairFA
    Int? threads = 4
   
    String outbase = basename(basename(readsFA, ".fasta"), ".fa")
  }
  command{
    rkmh2 filter -t ${threads} \
      -r ~{refFA} \
      -f ~{readsFA} \
      -t ~{threads} \
      ~{"-R " + refBatch} \
      ~{"-F " + readBatch} \
      ~{"-k " + kmerSize} \
      -t ~{threads} \
      ~{"-l " + minLength} \
      ~{"-s " + sketchSize} \
      ~{"-m " + threshold} > ${outbase}.filtered.fa
  }
  runtime{
    docker : "hpobiolab/rkmh2"
    cpu : threads
    memory : "12GB"
  }
  output{
    File filteredReads = "${outbase}.filtered.fa"
  }
}


task bwaMem{
    input{
        File reads
        File refFA
        File refFAI
        File refSA
        File refBWT
        File refPAC
        File refAMB
        File refANN
        #Boolean paired
        Int? bwaThreads = 20
        Int? samtoolsThreads = 10
        String outbase = basename(reads, ".bam")
        Int totalThreads = bwaThreads + samtoolsThreads
    }
    command{
      bwa mem -t ${bwaThreads} ${refFA} ${reads} | \
        samtools sort -T tmp -O bam -m 2G -@ ${samtoolsThreads} > ${outbase}.sorted.bam
    }
    runtime{
      docker : "erictdawson/bwa"
      memory : "38GB"
      cpu : totalThreads
      preemptible : 2
    }
    output{
      File mappedBam = "${outbase}.sorted.bam"
    }
}

task freebayesCall{
  input{
    File refFA
    File readsBAM
    String outbase = basename(readsBAM, ".bam")
  }
  command {
    freebayes --pooled-continuous --ploidy 1 --bam ${readsBAM} -f ${refFA} > ${outbase}.vcf
  }
  runtime {
    docker : "quay.io/biocontainers/freebayes:1.3.2--py27hc088bd4_0"
    cpu : 2
    memory : "13GB"
  }
  output{
    File freebayesVCF = "${outbase}.vcf"
  }

}

task samtoolsMpileup{
  input{
    File refFA
    File readsBAM
    Int? preemptAttempts = 2

    String outbase = basename(readsBAM, ".bam")
  }
  command{
    samtools mpileup --output-QNAME ${refFA} ${readsBAM} > ${outbase}.mpileup
  }
  runtime{
    docker : "erictdawson/samtools"
    cpu : 2
    memory : "7GB"
    preemptible : preemptAttempts
  }
  output{
    File readsPileup = "${outbase}.mpileup"
  }
}

task translate{
  input{
    File readsVARS
    File refGTF
  }
  command{
    echo
  }
  runtime{

  }
  output{

  }
}

task plotVariableSites{
  input{

  }
  command{
    echo 
  }
  runtime{
    docker : "rocker/tidyverse"
    cpu : 1
    memory : "7GB"
    preemptible : 0
  }
  output{

  }
}

task gffToTidy{
  input{

  }
  command{

  }
  runtime{

  }
  output{

  }
}

workflow findVariableSites{
  input{
    File refFA
    File refFAI
    File refGFF
    File readsFA

    File bwaSA
    File bwaPAC
    File bwaANN
    File bwaAMB

    Int? bwaThreads
    Int? samtoolsThreads
    Int? rkmhThreads
    Int? rkmhThreshold
    Int? rkmhKmerSize
    Int? rkmhSketchSize

    Int rkmhDiskGB = ceil(size(refFA, "gb") + (2 * size(readsFA, "gb")))
  }

  call rkmhFilter{
    input:
      refFA=refFA,
      readsFA=readsFA,
      threads=rkmhThreads,
      diskGB=rkmhDiskGB
  }

  #call gff_to_tidy{
    
  #}

  call bwaMem{
    input:
      reads=rkmhFilter.filteredReads,
      refFA=refFA,
      refFAI=refFAI,
      refSA=bwaSA,
      refPAC=bwaPAC,
      refANN=bwaANN,
      refAMB=bwaAMB

  }
  call freebayesCall{

  }
  call translate{

  }
  call plotVariableSites{

  }
}

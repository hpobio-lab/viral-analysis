version 1.0

task filterViralReads{
  input{
    File refFA
    File readsFA
    File? readPairFA
  }
  command{

  }
  runtime{

  }
  output{

  }
}


task bwaMem{
    input{
        File reads
        File ref
        File refSA
        File refBWT
        File refPAC
        File refAMB
        File refANN
        File refFAI
        #Boolean paired
        Int? bwaThreads = 20
        Int? samtoolsThreads = 10
        String outbase = basename(reads, ".bam")
        Int totalThreads = bwaThreads + samtoolsThreads
    }
    command{
      bwa mem -t ${bwaThreads} ${ref} ${reads} | \
        samtools -T tmp -O bam -m 2G -@ ${samtoolsThreads} > ${outbase}.sorted.bam
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

task samtoolsCall{
input{

  }
  command{

  }
  runtime{

  }
  output{

  }
}

task translate{
input{
    File readsVARS
    File refGTF
  }
  command{

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

workflow findVariableSites{

}

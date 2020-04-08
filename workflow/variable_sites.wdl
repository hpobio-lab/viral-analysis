version 1.0

task filterViralReads{
  input{
    File refFA
    File readsFA
    Int diskGB
    Int? kmerSize = 16
    Int? sketchSize = 10000
    Int? threshold = 50
    #File? readPairFA
    Int? threads = 4
   
    String outbase = basename(basename(readsFA, ".fasta"), ".fa")
  }
  command{
    rkmh stream -t ${threads} -r ${refFA} -f ${readsFA} -k 16 -s 5000 | \
    awk '{ if($3 >= ${threshold}) { print }}' | cut -f 2 > reads.lst && \
    seqtk ${readsFA} reads.lst > ${outbase}.filtered.fa
  }
  runtime{
    docker : "hpobiolab/rkmh"
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
    String outbase = basename(readsBam, ".bam")
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
    File bwaSAI
    File bwaPAC
    File bwaANN
    File bwaAMB

    Int? bwaThreads
    Int? samtoolsThreads
    Int? rkmhThreads
    Int? rkmhSimilarityThreshold
    Int? rkmhKmerSize
    Int? rkmhSketchSize
  }
}

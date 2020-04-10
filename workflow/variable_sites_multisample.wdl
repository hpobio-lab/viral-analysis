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
docker : "hpobiolab/rkmh2-ivybridge"
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
      samtools sort -T tmp -O bam -m 1G -@ ${samtoolsThreads} > ${outbase}.sorted.bam
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
      File refFAI
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
    File readsVars
    File refTidyGFF
    File refFA
    File refFAI

  }

  String outbase = basename(readsVars, ".vcf")
  command{
    python translate.py -g ~{refTidyGFF} -r ~{refFA} -i ~{readsVars} > ~{outbase}.tidytranslate.tsv
  }
  runtime {
    docker : "hpobiolab/variable-sites-translate"
    cpu : 1
    memory : "1GB"
    preemptible : 2
  }
  output{
    File tidyTranslatedVars = "~{outbase}.tidytranslate.tsv"
  }
}

task plotVariableSites{
  input{
    File tidyVars
  }
  String outbase = basename(tidyVars, ".tsv")
  command{
    Rscript --vanilla plot_variable_sites.R -i ~{tidyVars} -o ~{outbase}
  }
  runtime{
      docker : "rocker/tidyverse"
      cpu : 1
      memory : "7GB"
      preemptible : 2
  }
  output{
     File perFeatureCountPlot = "~{outbase}.per_feature_count_plot.pdf"
     File featureVarsPlot = "~{outbase}.feature_variants_plot.pdf"
     File genomeVarsPlot = "~{outbase}.genome_variants_plot.pdf"
     #File report = "~{outbase}.report.pdf"
  }
}

task gffToTidy{
  input{
    File refGFF
  }
  command{

  }
  runtime{

  }
  output{

  }
}

task collapseCalls{
  input{
    File fileList
    String outbase = "multisample_calls"
  }
  command{
    cat_file_of_files ~{fileList} > ~{outbase}.tidyVars.tsv
  }
  runtime{
    docker : "erictdawson/base"
    cpu : 1
    memory : "2GB"
    preemptible : 1
  }
  output {
    File collapsedVars = "~{outbase}.tidyVars.tsv"
  }
}



workflow findVariableSites{
  input{
    File refFA
      File refFAI
      File refTidyGFF

      File readsFAManifest

      File bwaSA
      File bwaPAC
      File bwaANN
      File bwaAMB
      File bwaBWT

      Int? bwaThreads
      Int? samtoolsThreads
      Int? rkmhThreads
      Int? rkmhThreshold
      Int? rkmhKmerSize
      Int? rkmhSketchSize

  }
  Int rkmhDiskGB = ceil(size(refFA, "GB") + (2 * size(readsFA, "GB")))

  Array[File] readsFA = read_lines(readsFAManifest)

  scatter (r in readsFA){

    call rkmhFilter{
      input:
        refFA=refFA,
        readsFA=r,
        threads=rkmhThreads,
        diskGB=rkmhDiskGB
    }

    call bwaMem{
      input:
        reads=rkmhFilter.filteredReads,
        refFA=refFA,
        refFAI=refFAI,
        refSA=bwaSA,
        refPAC=bwaPAC,
        refANN=bwaANN,
        refAMB=bwaAMB,
        refBWT=bwaBWT

    }
    call freebayesCall{
      input:
        refFA=refFA,
        refFAI=refFAI,
        readsBAM=bwaMem.mappedBam
    }
    call translate{
      input:
        refFA=refFA,
        refFAI=refFAI,
        refTidyGFF=refTidyGFF,
        readsVars = freebayesCall.freebayesVCF
    }
  }

  File tidyVarFiles = write_lines(translate.tidyTranslatedVars)
  call collapseCalls{
    input:
      fileList=tidyVarFiles      
  }

  call plotVariableSites{
    input:
      tidyVars=collapseCalls.collapsedVars
  }
}

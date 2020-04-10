FROM rocker/tidyverse

RUN apt-get update -q -yy && \
    apt-get install -q -yy \
    libgsl-dev \
    tabix && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages("data.table")'
RUN Rscript -e 'install.packages("seqinr")'

RUN apt-get update -q -yy && \
    apt-get install -q -yy \
    libbz2-dev \
    samtools && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager");'

RUN apt-get update -q -yy && \
    apt-get install -q -yy \
    liblzma-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Rhtslib");'
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Rsamtools");'
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("TFBSTools");'
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings");'
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("rtracklayer");'
RUN Rscript -e 'install.packages("foreach")'
RUN Rscript -e 'install.packages("doParallel")'


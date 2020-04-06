from rocker/verse

RUN Rscript -e 'install.packages("ggplot2")'

RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("DESeq2"); BiocManager::install("EdgeR"); BiocManager::install("limma")'
RUN Rscript -e 'install.packages("devtools"); library(devtools); install_github("PF2-pasteur-fr/SARTools", build_opts="--no-resave-data")'

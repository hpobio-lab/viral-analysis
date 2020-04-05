FROM ubuntu:18.04

RUN apt-get update &&  apt-get install -y \
    gcc \
    make \
    wget && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://mafft.cbrc.jp/alignment/software/mafft-7.453-with-extensions-src.tgz && \
    tar xvzf mafft-7.453-with-extensions-src.tgz && \
    cd mafft-7.453-with-extensions/core/ && \
    make clean && make -j 4 && make install

FROM erictdawson/base
LABEL maintainer="Eric T Dawson"
WORKDIR /app

RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.4.1/bowtie2-2.4.1-linux-x86_64.zip && \
    unzip bowtie2-2.4.1-linux-x86_64.zip && \
    cd bowtie2-2.4.1-linux-x86_64 && \
    cp bowtie2 /usr/bin/ && \
    cp bowtie2-align-s /usr/bin/ && \
    cp bowtie2-build-s /usr/bin/ && \
    cp bowtie2-align-l /usr/bin/ && \
    cp bowtie2-build-l /usr/bin/ && \
    cp bowtie2-inspect /usr/bin/ && \
    cp bowtie2-inspect-s /usr/bin/ && \
    cp bowtie2-inspect-l /usr/bin/ 

    #bowtie2, bowtie2-align-s, bowtie2-align-l, bowtie2-build, bowtie2-build-s, bowtie2-build-l, bowtie2-inspect, bowtie2-inspect-s and bowtie2-inspect-l

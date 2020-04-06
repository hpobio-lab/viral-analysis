FROM erictdawson/base
LABEL maintainer="Eric T Dawson"
WORKDIR /app

RUN wget https://github.com/lh3/minimap2/releases/download/v2.16/minimap2-2.16.tar.bz2 && \
    tar xjf minimap2-2.16.tar.bz2 && \
    cd minimap2-2.16 && \
    make && \
    cp minimap2 /usr/bin

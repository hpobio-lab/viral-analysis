FROM erictdawson/base
LABEL maintainer="Eric T Dawson"
WORKDIR /app

RUN apt install -y build-essential cmake && rm -rf /var/lib/apt/lists/*

RUN git clone --recursive https://github.com/ekg/seqwish.git && \
    cd seqwish && \
    cmake -H. -Bbuild && cmake --build build -- -j3 && \
    cp bin/seqwish /usr/bin/seqwish

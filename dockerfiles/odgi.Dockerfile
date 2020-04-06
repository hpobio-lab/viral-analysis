FROM erictdawson/base
LABEL maintainer="Eric T Dawson"
WORKDIR /app

RUN apt-get update &&  apt-get install -y \
    python3-dev  && \
    rm -rf /var/lib/apt/lists/*
RUN git clone --recursive https://github.com/vgteam/odgi.git && \
    cd odgi && \
    cmake -DBUILD_STATIC=1 -H. -Bbuild && cmake --build build -- -j 3 && \
    cp bin/odgi /usr/bin/odgi

FROM ubuntu:18.04
COPY --from=0 /usr/bin/odgi /usr/bin/odgi

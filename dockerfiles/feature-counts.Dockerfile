FROM ubuntu:18.04
WORKDIR /app
COPY subread-2.0.0-Linux-x86_64.tar.gz /app
RUN tar xvzf subread-2.0.0-Linux-x86_64.tar.gz && \
    cd  subread-2.0.0-Linux-x86_64 && \
    cp -r  bin/* /usr/bin/ && \
    cp -r bin/utilities/* /usr/bin/

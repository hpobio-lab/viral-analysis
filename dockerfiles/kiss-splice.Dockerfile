FROM erictdawson/base
WORKDIR /app
COPY kissplice-2.5.0.tar.gz /app
RUN tar xvzf kissplice-2.5.0.tar.gz && \
    cd kissplice-2.5.0 && \
    cmake .

FROM erictdawson/base
LABEL maintainer="Eric T Dawson"
WORKDIR /app

RUN apt-get update && \
    apt-get install -q -yy python3-setuptools && \
    rm -rf /var/lib/apt/lists/*
RUN pip3 install pyfaidx intervaltree

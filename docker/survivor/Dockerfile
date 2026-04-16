FROM ubuntu:22.04
LABEL maintainer="Ayan Malakar" \
      version="1.0.7" \
      description="SURVIVOR SV merger for ensemble pipeline"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        git \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        bcftools \
        tabix \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/fritzsedlazeck/SURVIVOR.git /tmp/SURVIVOR

WORKDIR /tmp/SURVIVOR/Debug
RUN make && cp SURVIVOR /usr/local/bin/

WORKDIR /workspace
RUN rm -rf /tmp/SURVIVOR

CMD ["SURVIVOR", "merge"]
FROM ubuntu:22.04
LABEL maintainer="Ayan Malakar" \
      version="1.0.7" \
      description="SURVIVOR SV merger + bcftools for ensemble pipeline"

ENV DEBIAN_FRONTEND=noninteractive \
    SURVIVOR_VERSION=1.0.7

RUN apt-get update && apt-get install -y \
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

# Build SURVIVOR from source
RUN git clone https://github.com/fritzsedlazeck/SURVIVOR.git \
    && cd SURVIVOR/Debug \
    && make \
    && cp SURVIVOR /usr/local/bin/ \
    && cd ../.. && rm -rf SURVIVOR

WORKDIR /workspace
CMD ["SURVIVOR", "merge"]

FROM ubuntu:22.04
LABEL maintainer="Ayan Malakar" \
      version="1.4.5" \
      description="BreakDancer SV caller for WGS ensemble pipeline"

ENV DEBIAN_FRONTEND=noninteractive \
    BREAKDANCER_VERSION=1.4.5

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget git build-essential cmake \
        python3 python3-pip \
        libncurses5-dev zlib1g-dev \
        libbz2-dev liblzma-dev \
        libcurl4-openssl-dev \
        libgsl-dev libboost-all-dev \
        perl cpanminus \
        bcftools tabix \
    && rm -rf /var/lib/apt/lists/*

# Perl modules required by bam2cfg.pl
RUN cpanm Statistics::Descriptive

# Install samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 \
    && tar -xjf samtools-1.17.tar.bz2 \
    && cd samtools-1.17 && ./configure --prefix=/usr/local && make -j4 && make install \
    && cd .. && rm -rf samtools-1.17*

# Build BreakDancer from source
RUN git clone https://github.com/genome/breakdancer.git \
    && cd breakdancer \
    && mkdir build && cd build \
    && cmake .. && make -j4 \
    && cp breakdancer-max /usr/local/bin/ \
    && cp ../perl/bam2cfg.pl /usr/local/bin/ \
    && cd ../.. && rm -rf breakdancer

# breakdancer2vcf conversion script
COPY breakdancer2vcf.py /usr/local/bin/breakdancer2vcf.py
RUN chmod +x /usr/local/bin/breakdancer2vcf.py

WORKDIR /workspace
CMD ["breakdancer-max", "--help"]

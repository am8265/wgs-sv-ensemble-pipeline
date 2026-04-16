FROM ubuntu:22.04
LABEL maintainer="Ayan Malakar" \
      version="0.4.1" \
      description="CNVnator read-depth CNV caller for WGS ensemble pipeline"

ENV DEBIAN_FRONTEND=noninteractive \
    CNVNATOR_VERSION=0.4.1 \
    ROOT_VERSION=6.28.04

RUN apt-get update && apt-get install -y \
        wget git build-essential cmake \
        libncurses5-dev zlib1g-dev \
        libbz2-dev liblzma-dev \
        libcurl4-openssl-dev \
        libx11-dev libxpm-dev libxft-dev \
        libssl-dev libgsl-dev \
        perl bcftools tabix \
    && rm -rf /var/lib/apt/lists/*

# Install ROOT (required by CNVnator)
RUN wget -q https://root.cern/download/root_v${ROOT_VERSION}.Linux-ubuntu22-x86_64-gcc11.3.tar.gz \
    && tar -xzf root_v${ROOT_VERSION}.Linux-ubuntu22-x86_64-gcc11.3.tar.gz \
    && mv root /opt/root \
    && rm root_v${ROOT_VERSION}.Linux-ubuntu22-x86_64-gcc11.3.tar.gz

ENV ROOTSYS=/opt/root
ENV PATH="${ROOTSYS}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${ROOTSYS}/lib:${LD_LIBRARY_PATH}"

# Install samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 \
    && tar -xjf samtools-1.17.tar.bz2 \
    && cd samtools-1.17 && ./configure --prefix=/usr/local && make -j4 && make install \
    && cd .. && rm -rf samtools-1.17*

# Install CNVnator
RUN git clone https://github.com/abyzovlab/CNVnator.git \
    && cd CNVnator \
    && ln -s /usr/local samtools \
    && make HTSDIR=/usr/local \
    && cp cnvnator /usr/local/bin/ \
    && cp src/cnvnator2VCF.pl /usr/local/bin/ \
    && cd .. && rm -rf CNVnator

WORKDIR /workspace
CMD ["cnvnator", "--help"]

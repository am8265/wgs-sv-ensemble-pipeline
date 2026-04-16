FROM ubuntu:22.04
LABEL maintainer="Ayan Malakar" \
      version="1.6.0" \
      description="Manta SV caller for WGS ensemble pipeline"

ENV DEBIAN_FRONTEND=noninteractive \
    MANTA_VERSION=1.6.0

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        python3 \
        bzip2 \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        tabix \
        bcftools \
    && rm -rf /var/lib/apt/lists/*

RUN wget -q https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.centos6_x86_64.tar.bz2 \
    && tar -xjf manta-${MANTA_VERSION}.centos6_x86_64.tar.bz2 \
    && mv manta-${MANTA_VERSION}.centos6_x86_64 /opt/manta \
    && rm manta-${MANTA_VERSION}.centos6_x86_64.tar.bz2

RUN wget -q https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 \
    && tar -xjf samtools-1.17.tar.bz2 \
    && rm samtools-1.17.tar.bz2

WORKDIR /tmp/samtools-1.17
RUN ./configure --prefix=/usr/local && make -j4 && make install

WORKDIR /workspace
RUN rm -rf /tmp/samtools-1.17

ENV PATH="/opt/manta/bin:${PATH}"

CMD ["configManta.py", "--help"]
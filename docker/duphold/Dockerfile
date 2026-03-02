FROM ubuntu:22.04
LABEL maintainer="Ayan Malakar" \
      version="0.2.3" \
      description="Duphold depth fold-change QC for SV calls"

ENV DEBIAN_FRONTEND=noninteractive \
    DUPHOLD_VERSION=0.2.3

RUN apt-get update && apt-get install -y \
        wget \
        bcftools \
        tabix \
        zlib1g-dev \
        libbz2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN wget -q https://github.com/brentp/duphold/releases/download/v${DUPHOLD_VERSION}/duphold \
    -O /usr/local/bin/duphold \
    && chmod +x /usr/local/bin/duphold

WORKDIR /workspace
CMD ["duphold", "--help"]

FROM ubuntu:22.04
LABEL maintainer="Ayan Malakar" \
      version="1.1.6" \
      description="Delly SV caller for WGS ensemble pipeline"

ENV DEBIAN_FRONTEND=noninteractive \
    DELLY_VERSION=1.1.6

RUN apt-get update && apt-get install -y \
        wget \
        bcftools \
        tabix \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

RUN wget -q https://github.com/dellytools/delly/releases/download/v${DELLY_VERSION}/delly_v${DELLY_VERSION}_linux_x86_64bit \
    -O /usr/local/bin/delly \
    && chmod +x /usr/local/bin/delly

WORKDIR /workspace
CMD ["delly", "--help"]

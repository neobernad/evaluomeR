FROM bioconductor/bioconductor_docker:RELEASE_3_19

LABEL org.opencontainers.image.source="https://github.com/neobernad/evaluomeR"
LABEL org.opencontainers.image.description="evaluomeR — evaluation of bioinformatics metrics (built from source)"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake \
    gfortran \
    liblapack-dev \
    libblas-dev \
    libnlopt-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    pandoc \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /work

# Default: interactive R shell with repo mounted at /work via docker-compose
CMD ["R"]

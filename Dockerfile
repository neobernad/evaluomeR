FROM rocker/r-ver:4.4.2

LABEL org.opencontainers.image.source="https://github.com/neobernad/evaluomeR"
LABEL org.opencontainers.image.description="evaluomeR — evaluation of bioinformatics metrics (built from source)"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
  && rm -rf /var/lib/apt/lists/*

RUN R -e 'install.packages(c("remotes", "BiocManager"), repos = "https://packagemanager.posit.co/cran/__linux__/jammy/latest")'

WORKDIR /evaluomeR

COPY DESCRIPTION NAMESPACE LICENSE ./
COPY R/ R/
COPY data/ data/
COPY man/ man/
COPY inst/ inst/
COPY vignettes/ vignettes/

RUN R -e 'BiocManager::install(c("SummarizedExperiment", "MultiAssayExperiment"), ask = FALSE, update = FALSE)' \
 && R -e 'remotes::install_local(".", dependencies = TRUE, upgrade = FALSE)'

CMD ["R"]

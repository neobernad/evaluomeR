FROM rocker/r-ver:latest

LABEL org.opencontainers.image.source="https://github.com/neobernad/evaluomeR"
LABEL org.opencontainers.image.description="evaluomeR — evaluation of bioinformatics metrics (built from source)"

ENV DEBIAN_FRONTEND=noninteractive

# System libraries required by evaluomeR's compiled dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake \
    gfortran \
    liblapack-dev \
    libblas-dev \
    libnlopt-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    pandoc \
    libuv1-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg-dev \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /work

# Bootstrap BiocManager and docs toolchain.
# install.packages() exits 0 even on failure, so we verify afterwards.
RUN R -e ' \
  pkgs <- c("BiocManager", "remotes", "xfun", "knitr", "rmarkdown", "pkgdown"); \
  install.packages(pkgs, repos = "https://cloud.r-project.org", Ncpus = 2L); \
  ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE); \
  if (!all(ok)) stop("Failed to install: ", paste(names(ok)[!ok], collapse = ", ")) \
'

# Pre-install all evaluomeR Depends + Imports via BiocManager.
# (covers Bioconductor packages and compiled CRAN deps like lme4 / RcppEigen)
RUN R -e ' \
  pkgs <- c( \
    "SummarizedExperiment", "MultiAssayExperiment", \
    "dplyr", "cluster", "fpc", "randomForest", "flexmix", "RSKC", "sparcl", \
    "ggrepel", "corrplot", "reshape2", "ggplot2", "ggdendro", "plotrix", \
    "matrixStats", "Rdpack", "MASS", "prabclus", "mclust", "kableExtra", \
    "dendextend", "fossil", "psych", "FactoMineR", "factoextra", "nFactors", \
    "RcppEigen", "lme4" \
  ); \
  BiocManager::install(pkgs, update = FALSE, ask = FALSE); \
  ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE); \
  if (!all(ok)) stop("Failed to install: ", paste(names(ok)[!ok], collapse = ", ")) \
'

CMD ["R"]

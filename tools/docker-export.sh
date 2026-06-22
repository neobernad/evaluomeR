#!/usr/bin/env bash
# Install evaluomeR and dependencies inside the Docker image, then export demo JSON.
# Usage: docker-compose run --rm evaluomeR bash tools/docker-export.sh

set -euo pipefail

apt-get update -qq
apt-get install -y -qq \
  cmake gfortran liblapack-dev libblas-dev libnlopt-dev \
  libcurl4-openssl-dev libssl-dev libxml2-dev pandoc
# Expected runtime: ~20-40 min (NCI-60, k=3..8, bs=100)

R -e 'install.packages(c(
  "jsonlite", "remotes", "RcppEigen", "lme4", "xfun", "knitr", "kableExtra"
), repos = "https://cloud.r-project.org")'

R -e 'remotes::install_local("/work", dependencies = TRUE, upgrade = FALSE, build = FALSE)'

Rscript tools/export-demo-data.R

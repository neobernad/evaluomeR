# Docker usage for evaluomeR

Run **evaluomeR** in a reproducible Linux container with R and all Bioconductor/CRAN dependencies pre-installed. The image builds the package from this repository source.

## Prerequisites

- [Docker Engine](https://docs.docker.com/engine/install/) or Docker Desktop

## Build the image

From the repository root:

```bash
docker build -t evaluomer:local .
```

The first build may take 10–20 minutes while R packages are compiled. Subsequent builds reuse cached layers when only source code changes.

## Verify installation

```bash
docker run --rm evaluomer:local R -e 'packageVersion("evaluomeR")'
```

```bash
docker run --rm evaluomer:local R -e 'library(evaluomeR); data("rnaMetrics"); print(rnaMetrics)'
```

## Interactive R session

```bash
docker run -it --rm evaluomer:local
```

Inside R:

```r
library(evaluomeR)
data("rnaMetrics")
stability(rnaMetrics, k = 2, bs = 20, seed = 100, getImages = FALSE)
```

## Example one-shot commands

### Correlation analysis

```bash
docker run --rm evaluomer:local R -e '
  library(evaluomeR)
  data("rnaMetrics")
  corSE <- metricsCorrelations(rnaMetrics)
  print(assay(corSE, 1)[1:3, 1:3])
'
```

### Stability index (single k)

```bash
docker run --rm evaluomer:local R -e '
  library(evaluomeR)
  data("rnaMetrics")
  stab <- stability(rnaMetrics, k = 2, bs = 20, seed = 100, getImages = FALSE)
  print(assay(stab$stability_mean))
'
```

### Stability + quality range and optimal k

```bash
docker run --rm evaluomer:local R -e '
  library(evaluomeR)
  data("rnaMetrics")
  s <- stabilityRange(rnaMetrics, k.range = c(3, 4), bs = 20, seed = 100, getImages = FALSE)
  q <- qualityRange(rnaMetrics, k.range = c(3, 4), seed = 100, getImages = FALSE)
  print(getOptimalKValue(s, q, k.range = c(3, 4)))
'
```

### Ontology metrics (bundled dataset)

```bash
docker run --rm evaluomer:local R -e '
  library(evaluomeR)
  data("ontMetrics")
  stab <- stability(ontMetrics, k = 3, bs = 20, seed = 100, getImages = FALSE)
  print(assay(stab$stability_mean))
'
```

### ATSC pipeline (reduced parameters for a quick demo)

ATSC runs many bootstrap iterations; use small `bs` and a narrow `k.range` for demos:

```bash
docker run --rm evaluomer:local R -e '
  library(evaluomeR)
  data("rnaMetrics")
  result <- ATSC(rnaMetrics, k.range = c(2, 3), bs = 20, seed = 100)
  print(result)
'
```

## Custom data via volume mount

Mount a host directory and run your own R script:

```bash
mkdir -p work
docker run --rm -v "$PWD/work:/work" -w /work evaluomer:local Rscript my_analysis.R
```

Example `work/my_analysis.R` using a CSV (first column = sample ID, remaining columns = metrics):

```r
library(evaluomeR)
library(SummarizedExperiment)

df <- read.csv("metrics.csv", check.names = FALSE)
se <- SummarizedExperiment(assays = list(counts = data.matrix(df[, -1])),
                           rowData = data.frame(Description = df[, 1]))

stab <- stability(se, k = 2, bs = 20, seed = 100, getImages = FALSE)
print(assay(stab$stability_mean))
```

## Docker Compose

```bash
mkdir -p work
docker compose run --rm evaluomeR
```

This starts an interactive R session with `./work` mounted at `/work`.

## Development smoke test (optional)

To run integration scripts from a mounted checkout:

```bash
docker run --rm -v "$PWD:/repo" -w /repo evaluomer:local \
  Rscript -e 'source("tests/testAll.R")'
```

Note: some test scripts reference datasets not bundled in the package.

## Parallel analyses

Functions such as `stability(..., numCores = 2)` spawn R child processes. Allocate CPUs when needed:

```bash
docker run --rm --cpus 4 evaluomer:local R -e '
  library(evaluomeR)
  data("rnaMetrics")
  stability(rnaMetrics, k = 2, bs = 20, seed = 100, numCores = 2, getImages = FALSE)
'
```

## Alternative: Bioconductor release image

If you only need the stable Bioconductor release (not this repo’s source), a minimal Dockerfile is:

```dockerfile
FROM rocker/r-ver:4.4.2
RUN R -e 'install.packages("BiocManager"); BiocManager::install("evaluomeR")'
CMD ["R"]
```

## See also

- [README](../README.md) — Bioconductor and Conda install options
- [Bioconductor evaluomeR](https://bioconductor.org/packages/release/bioc/html/evaluomeR.html)
- Conda: `conda install -c bioconda bioconductor-evaluomer`

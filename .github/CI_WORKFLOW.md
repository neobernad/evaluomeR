# GitHub Actions CI Workflow

This document describes the continuous integration (CI) workflows for evaluomeR.

## Overview

CI validates the installable R package. `R CMD check` runs the testthat suite
in `tests/testthat/` (~10 s). The slow timing scripts in `tests/benchmarks/` are
left out of the check and run weekly in the Benchmarks workflow.

## Workflows

### R CMD check — `.github/workflows/test.yaml`

**Triggers:** push or pull request to `main`, `master`, or `develop`; manual
(`workflow_dispatch`).

**Matrix:**

| OS | R version |
|----|-----------|
| ubuntu-latest | release |

**Steps:** checkout → setup-r (RSPM binaries) → setup-r-dependencies →
`r-lib/actions/check-r-package`.

This runs `R CMD check`, which exercises examples, vignette code, the testthat
suite and package metadata. A failing check blocks merge. A new push to the
same branch cancels the run it supersedes.

### Benchmarks — `.github/workflows/benchmarks.yaml`

**Triggers:** weekly (Mondays 03:00 UTC); manual dispatch.

Runs `testL1_clustering.R`, `testPCA.R` and `testCH_index.R` from
`tests/benchmarks/`. They have no assertions, so this fails only on a crash.

### pkgdown — `.github/workflows/pkgdown.yaml`

**Triggers:** push to `main`, `master`, or `develop`; manual dispatch.

Builds the documentation site (and demo SPA) on every trigger. **Deploy** to
GitHub Pages at https://neobernad.github.io/evaluomeR/ runs only on push to
`main` or `master` — `develop` validates the build but does not publish
(must match the `github-pages` environment branch protection).

## Local testing

```r
devtools::test()                                        # whole suite
testthat::test_file("tests/testthat/test-stability.R")  # one area, see AGENTS.md
```

Full package check locally:

```bash
R CMD build . && R CMD check evaluomeR_*.tar.gz
```

## Badge

```markdown
[![R-CMD-check](https://github.com/neobernad/evaluomeR/actions/workflows/test.yaml/badge.svg)](https://github.com/neobernad/evaluomeR/actions/workflows/test.yaml)
```

## Related files

- `.github/workflows/test.yaml` — R CMD check workflow
- `.github/workflows/pkgdown.yaml` — documentation deploy
- `AGENTS.md` — area-specific test mapping for agents
- `.github/workflows/benchmarks.yaml` — weekly benchmark run
- `tests/testthat/` — test suite run by R CMD check
- `tests/benchmarks/` — timing scripts (not run by R CMD check)

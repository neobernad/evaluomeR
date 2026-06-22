# GitHub Actions CI Workflow

This document describes the continuous integration (CI) workflows for evaluomeR.

## Overview

CI validates the installable R package. The standalone scripts in `tests/` are
**not** run by GitHub Actions; run those locally after code changes (see
`AGENTS.md` → Agent Workflow).

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

This runs `R CMD check`, which exercises examples, vignette code, and package
metadata. A failing check blocks merge.

### pkgdown — `.github/workflows/pkgdown.yaml`

**Triggers:** push to `main`, `master`, or `develop`; manual dispatch.

Builds the documentation site (and demo SPA) on every trigger. **Deploy** to
GitHub Pages at https://neobernad.github.io/evaluomeR/ runs only on push to
`main` or `master` — `develop` validates the build but does not publish
(must match the `github-pages` environment branch protection).

## Local testing

```r
library(evaluomeR)
source("tests/testAll.R")       # smoke test
source("tests/testStability.R") # area-specific — see AGENTS.md
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
- `tests/` — standalone integration scripts (local only)

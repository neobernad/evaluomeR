# GitHub Actions CI Workflow

This document describes the continuous integration (CI) workflow for the evaluomeR package.

## Overview

The `.github/workflows/test.yaml` workflow automatically runs your R package tests whenever code is pushed or a pull request is created. This ensures code quality and prevents regressions before merging to main branches.

## Workflow Triggers

The workflow runs on the following events:

- **Push**: Any push to `main`, `master`, or `develop` branches
- **Pull Request**: Any pull request targeting `main`, `master`, or `develop` branches
- **Manual**: Manually triggered via GitHub Actions UI (workflow_dispatch)

## Test Matrix

Tests run across multiple R versions to ensure compatibility:

- **R 4.3** - Latest stable release from that series
- **Latest R** - Most current R version available

All testing occurs on **Ubuntu Latest** runner.

## Workflow Steps

### 1. Checkout
Clones your repository code.

### 2. Set up R
Installs the specified R version using the official r-lib setup action.

### 3. Install System Dependencies
Installs required system libraries:
- `libcurl4-openssl-dev` - For HTTP communication
- `libssl-dev` - For SSL/TLS support
- `libxml2-dev` - For XML parsing

### 4. Install R Dependencies
Automatically installs all package dependencies declared in `DESCRIPTION`:
- **Depends**: SummarizedExperiment, MultiAssayExperiment, dplyr, cluster, fpc, randomForest, flexmix, RSKC, sparcl
- **Imports**: All packages listed in the Imports field
- **Extra packages**: devtools (for package installation)

### 5. Install Package
Builds and installs your package locally using `devtools::install_local()`.

### 6. Run Tests
Executes all test scripts in the `tests/` directory matching the pattern `test*.R`:
- testARI.R
- testATSC.R
- testATSC_GoldStandard.R
- testATSC_parallel.R
- testAll.R
- testAlpha_parallel.R
- testAnalysis.R
- testCBI.R
- testCH_index.R
- testDataPreprocessing.R
- testK_optimo_parallel.R
- testL1_clustering.R
- testMetricsRelevancy.R
- testPCA.R
- testQuality.R
- testQuality_parallel.R
- testStability.R
- testStability_parallel.R

Each test is sourced in sequence, with detailed output showing which test is running. If any test fails, the workflow exits with status 1 and reports the error.

### 7. Upload Artifacts
Test results and scripts are uploaded as workflow artifacts for 30 days, allowing you to review test output even after the workflow completes.

## Viewing Results

### On GitHub
1. Navigate to your repository
2. Click the **Actions** tab
3. Click on the workflow run to see detailed output
4. Expand any step to see its console output
5. Download artifacts from the "Artifacts" section

### Badges
You can add a status badge to your README:

```markdown
[![R Package Tests](https://github.com/neobernad/evaluomeR/actions/workflows/test.yaml/badge.svg)](https://github.com/neobernad/evaluomeR/actions/workflows/test.yaml)
```

## Troubleshooting

### Tests Fail in CI but Pass Locally
- Ensure your local R version matches one of the CI versions (R 4.3 or latest)
- Check for platform-specific issues (paths, environment variables)
- Verify all dependencies are correctly listed in DESCRIPTION

### Missing Dependencies
- Add missing packages to DESCRIPTION under either `Depends` or `Imports`
- For system packages, add them to the "Install system dependencies" step
- Re-run the workflow after updating DESCRIPTION

### Long Test Execution Times
- Consider splitting large parallel tests into separate workflows if needed
- Optimize test code to reduce runtime
- Check if tests are downloading large data files

## Future Enhancements

Possible improvements to the CI workflow:

1. **Code Coverage**: Add code coverage reports using `covr` package
2. **R CMD Check**: Add `R CMD check` for comprehensive package validation
3. **Platform Matrix**: Extend to include macOS and Windows runners
4. **Documentation Builds**: Generate and validate vignettes/documentation
5. **Performance Benchmarks**: Track test execution times over time

## Related Files

- `.github/workflows/test.yaml` - The main workflow definition
- `DESCRIPTION` - Package metadata and dependencies
- `tests/` - Directory containing all test scripts

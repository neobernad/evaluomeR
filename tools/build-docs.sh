#!/usr/bin/env bash
# Build the full evaluomeR docs site (pkgdown) and the React demo app.
#
# Usage (from repo root):
#   bash tools/build-docs.sh
#
# What it does:
#   1. Always rebuilds the Docker image (picks up Dockerfile / base-image changes).
#   2. Verifies remotes + pkgdown are available in the image.
#   3. Inside Docker: installs evaluomeR, force-cleans docs/, builds pkgdown site.
#   4. On the host: runs `npm run build` to write the React demo to docs/app/.
#
# Prerequisites:
#   - Docker + docker-compose available
#   - Node.js + npm available (for the React app step)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

echo "──────────────────────────────────────────────────"
echo " Step 1: Build Docker image"
echo "──────────────────────────────────────────────────"
docker-compose build --no-cache

echo ""
echo "──────────────────────────────────────────────────"
echo " Step 2: Verify image has remotes + pkgdown"
echo "──────────────────────────────────────────────────"
docker-compose run --rm evaluomeR Rscript -e '
  pkgs <- c("remotes", "pkgdown")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing in image after build: ", paste(missing, collapse = ", "),
         ". Re-run: docker-compose build --no-cache")
  }
'

echo ""
echo "──────────────────────────────────────────────────"
echo " Step 3: Install evaluomeR + build pkgdown site"
echo "──────────────────────────────────────────────────"
docker-compose run --rm evaluomeR bash -c "
  set -euo pipefail
  R -e 'remotes::install_local(\"/work\", dependencies=TRUE, upgrade=FALSE, build=FALSE)'
  Rscript -e 'stopifnot(requireNamespace(\"evaluomeR\", quietly=TRUE))'
  Rscript -e 'pkgdown::clean_site(force=TRUE); pkgdown::build_site_github_pages()'
"

echo ""
echo "──────────────────────────────────────────────────"
echo " Step 4: Build React demo app → docs/app/"
echo "──────────────────────────────────────────────────"
cd "$REPO_ROOT/ui"
npm ci
npm run build

echo ""
echo "──────────────────────────────────────────────────"
echo " Done. Docs are in docs/"
echo " pkgdown site : docs/index.html"
echo " React demo   : docs/app/index.html"
echo "──────────────────────────────────────────────────"

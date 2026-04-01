#!/usr/bin/env bash
# ─── Run CI check locally in Docker ──────────────────────────────────────────
# Usage:
#   ./ci-local.sh               → full R CMD check, output to stdout
#   ./ci-local.sh logs          → same, but also saves Rcheck dir + log file
#   ./ci-local.sh shell         → interactive container for debugging
#   ./ci-local.sh build         → (re)build image only, keep cache
#   ./ci-local.sh rebuild       → force full rebuild (no cache)
#
# Error collection:
#   - Plain stdout:   ./ci-local.sh 2>&1 | tee ci-local.log
#   - Full Rcheck:    ./ci-local.sh logs   → writes ./ci-check-output/
#     The directory contains semseeker.Rcheck/ with all log files:
#       00check.log        → full R CMD check output
#       00install.out      → package installation log
#       semseeker/tests/testthat.Rout → all test output (every test_that block)
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

IMAGE="semseeker-ci"
MODE="${1:-check}"
OUTDIR="$(pwd)/ci-check-output"

_build() {
  docker build -f Dockerfile.ci -t "$IMAGE" "$@" .
}

case "$MODE" in

  build)
    echo "==> Building $IMAGE (cached) ..."
    _build
    ;;

  rebuild)
    echo "==> Force-rebuilding $IMAGE (no cache) ..."
    _build --no-cache
    ;;

  shell)
    echo "==> Building $IMAGE (if needed) ..."
    _build -q
    echo "==> Starting interactive shell in /pkg ..."
    echo "    Useful commands inside:"
    echo "      R CMD check --no-manual --as-cran /pkg"
    echo "      Rscript -e 'devtools::load_all(); testthat::test_file(\"tests/testthat/test-7-association_analysis.R\")'"
    docker run --rm -it -w /pkg "$IMAGE" bash
    ;;

  logs)
    echo "==> Building $IMAGE (if needed) ..."
    _build -q
    echo ""
    echo "==> Running R CMD check, collecting artefacts in $OUTDIR ..."
    mkdir -p "$OUTDIR"

    # Write the R script to a temp file to avoid shell-escaping issues
    RSCRIPT_TMP="$(mktemp /tmp/ci-check-XXXXXX.R)"
    cat > "$RSCRIPT_TMP" <<'REOF'
res <- rcmdcheck::rcmdcheck(
  args       = c("--no-manual", "--as-cran"),
  build_args = "--no-manual",
  error_on   = "never",
  check_dir  = "/check-output"
)
cat("\n=== SUMMARY ===\n")
print(res)
if (length(res$errors) > 0) quit(status = 1)
REOF

    docker run --rm \
      -v "$OUTDIR":/check-output \
      -v "$RSCRIPT_TMP":/tmp/ci-check.R:ro \
      "$IMAGE" \
      Rscript /tmp/ci-check.R

    EXIT_CODE=$?
    rm -f "$RSCRIPT_TMP"

    echo ""
    echo "==> Artefacts saved to: $OUTDIR"
    echo "    Key files:"
    echo "      $OUTDIR/semseeker.Rcheck/00check.log"
    echo "      $OUTDIR/semseeker.Rcheck/00install.out"
    echo "      $OUTDIR/semseeker.Rcheck/semseeker/tests/testthat.Rout"
    exit $EXIT_CODE
    ;;

  check|*)
    echo "==> Building $IMAGE (if needed) ..."
    _build
    echo ""
    echo "==> Running R CMD check ..."
    echo "    Tip: pipe through 'tee' to save: ./ci-local.sh 2>&1 | tee ci-local.log"
    echo "    Tip: use './ci-local.sh logs' to also extract the full .Rcheck directory"
    echo ""
    docker run --rm "$IMAGE"
    ;;

esac

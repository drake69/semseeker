#!/usr/bin/env bash
# ─── Run CI check locally in Docker ──────────────────────────────────────────
# Usage:
#   ./ci-local.sh               → full R CMD check, output to stdout
#   ./ci-local.sh check         → same as above
#   ./ci-local.sh coverage      → run covr::package_coverage() (mirrors test-coverage CI)
#   ./ci-local.sh logs          → R CMD check + saves .Rcheck dir to ./ci-check-output/
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
#   - Coverage:       ./ci-local.sh coverage
#     On failure prints the full .Rout.fail content so you can see which test failed.
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
    RSCRIPT_TMP="/tmp/ci-check-$$.R"
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

  coverage)
    echo "==> Building $IMAGE (if needed) ..."
    _build -q
    echo ""
    echo "==> Running covr::package_coverage() (mirrors test-coverage CI) ..."
    echo ""

    RSCRIPT_TMP="/tmp/ci-coverage-$$.R"
    cat > "$RSCRIPT_TMP" <<'REOF'
Sys.setenv(NOT_CRAN = "true", RENV_CONFIG_SANDBOX_ENABLED = "false")
if (!requireNamespace("covr", quietly = TRUE))
  install.packages("covr", repos = "https://packagemanager.posit.co/cran/__linux__/noble/latest")
install_path <- "/tmp/covr-pkg"
tryCatch({
  cov <- covr::package_coverage(
    path        = "/pkg",
    quiet       = FALSE,
    clean       = FALSE,
    install_path = install_path
  )
  cat("\n=== COVERAGE SUMMARY ===\n")
  print(cov)
}, error = function(e) {
  cat("\n=== COVERAGE FAILED ===\n", conditionMessage(e), "\n")
  # Print .Rout.fail if present
  fail_files <- list.files(install_path, pattern = "\\.Rout\\.fail$",
                           recursive = TRUE, full.names = TRUE)
  for (f in fail_files) {
    cat("\n=== FAIL FILE:", f, "===\n")
    cat(readLines(f), sep = "\n")
  }
  quit(status = 1)
})
REOF

    docker run --rm \
      -v "$RSCRIPT_TMP":/tmp/ci-coverage.R:ro \
      -e NOT_CRAN=true \
      -e RENV_CONFIG_SANDBOX_ENABLED=false \
      -e RENV_ACTIVATE_PROJECT=0 \
      "$IMAGE" \
      Rscript /tmp/ci-coverage.R

    EXIT_CODE=$?
    rm -f "$RSCRIPT_TMP"
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

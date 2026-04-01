#!/usr/bin/env bash
# ─── Run CI check locally in Docker ──────────────────────────────────────────
# Usage:
#   ./ci-local.sh          → full R CMD check (same as GitHub CI)
#   ./ci-local.sh shell    → drop into container for interactive debugging
#   ./ci-local.sh build    → (re)build image only
#   ./ci-local.sh rebuild  → force full rebuild (no cache)
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

IMAGE="semseeker-ci"
MODE="${1:-check}"

case "$MODE" in
  build)
    echo "==> Building $IMAGE ..."
    docker build -f Dockerfile.ci -t "$IMAGE" .
    ;;
  rebuild)
    echo "==> Force-rebuilding $IMAGE (no cache) ..."
    docker build --no-cache -f Dockerfile.ci -t "$IMAGE" .
    ;;
  shell)
    echo "==> Building $IMAGE (if needed) ..."
    docker build -f Dockerfile.ci -t "$IMAGE" . -q
    echo "==> Starting interactive shell ..."
    docker run --rm -it "$IMAGE" bash
    ;;
  check|*)
    echo "==> Building $IMAGE (if needed) ..."
    docker build -f Dockerfile.ci -t "$IMAGE" .
    echo ""
    echo "==> Running R CMD check ..."
    docker run --rm "$IMAGE"
    ;;
esac

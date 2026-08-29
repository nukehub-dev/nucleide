#!/usr/bin/env bash
# Build (layer-cached) and run the full cross-code validation harness in a
# container.
#
# Usage:
#   ./validation/run_container.sh
#
# The image bundles PyNE 0.7.5 (conda-forge; numerically identical to 0.7.8
# for the modules exercised here) and OpenMC v0.16.0 (built from the upstream
# release tag). See Containerfile for why these channels were chosen. The
# Nucleide abi3 wheel is rebuilt each run via the PyO3 maturin container so
# the harness always tests the current sources.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
IMAGE="${VALIDATION_IMAGE:-localhost/nucleide-validation:latest}"

# Always rebuild: the layer cache makes unchanged rebuilds near-instant, and
# Containerfile edits (e.g. new oracle packages) take effect without having to
# delete the image by hand.
echo "Building validation image $IMAGE (layer-cached)..."
podman build -t "$IMAGE" -f "$SCRIPT_DIR/Containerfile" "$REPO_ROOT"

echo "Building Nucleide release wheel (maturin container, cargo-cached)..."
mkdir -p "$REPO_ROOT/target/container-wheels"
podman run --rm -v "$REPO_ROOT:/io" ghcr.io/pyo3/maturin:latest \
    build --release -m bindings/python/Cargo.toml -o /io/target/container-wheels

echo ""
echo "Running validation harness..."
podman run --rm -v "$REPO_ROOT:/io" -w /io "$IMAGE" bash -c '
    set -e
    pip install --quiet --force-reinstall --no-deps \
        target/container-wheels/nucleide-*-manylinux*x86_64.whl
    ./validation/run_all.sh python3
'

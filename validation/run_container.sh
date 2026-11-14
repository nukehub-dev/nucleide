#!/usr/bin/env bash
# Build (layer-cached) and run the full cross-code validation harness in a
# container.
#
# Usage:
#   [CONTAINER_ENGINE=docker|podman] ./validation/run_container.sh
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
ENGINE="${CONTAINER_ENGINE:-}"
if [[ -z "$ENGINE" ]]; then
    if command -v podman >/dev/null 2>&1; then
        ENGINE=podman
    elif command -v docker >/dev/null 2>&1; then
        ENGINE=docker
    else
        echo "Install Docker or Podman, or use run_all.sh with a prepared Python environment." >&2
        exit 127
    fi
fi
case "$ENGINE" in
    docker|podman) ;;
    *) echo "CONTAINER_ENGINE must be docker or podman" >&2; exit 2 ;;
esac
if ! command -v "$ENGINE" >/dev/null 2>&1; then
    echo "Requested container engine is not installed: $ENGINE" >&2
    exit 127
fi
RUN_LIMITS=(--memory=6g --memory-swap=6g --pids-limit=128 --cpus=2)
echo "Using container engine: $ENGINE"

# Always rebuild: the layer cache makes unchanged rebuilds near-instant, and
# Containerfile edits (e.g. new oracle packages) take effect without having to
# delete the image by hand.
echo "Building validation image $IMAGE (layer-cached)..."
# The image copies no repository inputs. Restrict context to validation/
# rather than uploading target/, .venv/, node_modules/ and research caches.
"$ENGINE" build --build-arg BUILD_JOBS=2 -t "$IMAGE" \
    -f "$SCRIPT_DIR/Containerfile" "$SCRIPT_DIR"

echo "Building Nucleide release wheel (maturin container, cargo-cached)..."
mkdir -p "$REPO_ROOT/target/container-wheels"
# Stale wheels from previous runs would collide in the install glob below.
rm -f "$REPO_ROOT"/target/container-wheels/nucleide-*.whl
"$ENGINE" run --rm "${RUN_LIMITS[@]}" -v "$REPO_ROOT:/io" ghcr.io/pyo3/maturin:latest \
    build --release --jobs 2 -m bindings/python/Cargo.toml -o /io/target/container-wheels

echo ""
echo "Running validation harness..."
"$ENGINE" run --rm "${RUN_LIMITS[@]}" \
    -e OMP_NUM_THREADS=2 -e OPENBLAS_NUM_THREADS=2 -e MKL_NUM_THREADS=2 \
    -v "$REPO_ROOT:/io" -w /io "$IMAGE" bash -c '
    set -e
    pip install --quiet --force-reinstall --no-deps \
        target/container-wheels/nucleide-*-manylinux*x86_64.whl
    ./validation/run_all.sh python3
'

#!/usr/bin/env bash
# Run the full cross-code validation harness.
#
# Usage:
#   ./validation/run_all.sh [PYTHON]
#
# The optional PYTHON argument defaults to the conda env used for JOSS
# validation. Every validation/*_vs_*.py comparison script is auto-discovered
# and run in sorted order, followed by timings.py and render_results.py.
# All scripts are run with -e so the first failure stops the suite.

set -euo pipefail

PYTHON="${1:-/home/tahmid/.conda/envs/nuke-validation/bin/python}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Using Python: $PYTHON"
$PYTHON --version

for script in "$SCRIPT_DIR"/*_vs_*.py; do
    echo ""
    echo "===== $(basename "$script") ====="
    $PYTHON "$script"
done

for script in timings.py render_results.py; do
    echo ""
    echo "===== $script ====="
    $PYTHON "$SCRIPT_DIR/$script"
done

echo ""
echo "All validation scripts completed successfully."

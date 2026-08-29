#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2023-2026 Nucleide Developers
# SPDX-License-Identifier: BSD-2-Clause

# Bump the Nucleide workspace version.
#
# The git tag (vX.Y.Z) is the release source of truth; this script updates the
# checked-in version source so a release is one command:
#
#   scripts/bump-version.sh 0.2.0
#
# Updates:
#   Cargo.toml              - [workspace.package] version
#   Cargo.lock              - workspace member versions (committed lockfile)
#   website/package.json    - website package version
#   website/package-lock.json - website version fields (CI uses `npm ci`)
#   CHANGELOG.md            - stamps [Unreleased] with the new version + date
#
# The script never commits or tags; it prints the follow-up git commands.
# Re-running with the same version is a no-op.

set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

die() {
    echo "error: $*" >&2
    exit 1
}

_version="${1:-}"
[[ "$_version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] \
    || die "usage: scripts/bump-version.sh <X.Y.Z>  (e.g. scripts/bump-version.sh 0.2.0)"

_date="$(date +%F)"

# Cargo.toml: update [workspace.package] version only if it differs.
_cargo_toml="$DIR/Cargo.toml"
_current_version="$(sed -n 's/^version = "\([^"]*\)".*/\1/p' "$_cargo_toml" | head -n1)"
if [[ "$_current_version" == "$_version" ]]; then
    echo "note: Cargo.toml already at $_version; left unchanged"
else
    sed -i "s/^version = \"[^\"]*\"/version = \"$_version\"/" "$_cargo_toml"
    echo "  Cargo.toml -> $_version"
    # Keep the committed lockfile in sync (workspace member versions).
    if (cd "$DIR" && cargo metadata --no-deps --format-version 1 >/dev/null 2>&1); then
        echo "  Cargo.lock -> $_version"
    else
        echo "warning: could not refresh Cargo.lock (offline?); run any cargo command to sync it" >&2
    fi
fi

# website/package.json: update version only if it differs.
_website_pkg="$DIR/website/package.json"
_current_website_version="$(sed -n 's/^[[:space:]]*"version": "\([^"]*\)".*/\1/p' "$_website_pkg" | head -n1)"
if [[ "$_current_website_version" == "$_version" ]]; then
    echo "note: website/package.json already at $_version; left unchanged"
else
    sed -i "s/\"version\": \"[^\"]*\"/\"version\": \"$_version\"/" "$_website_pkg"
    echo "  website/package.json -> $_version"
    # Keep the website lockfile in sync without touching node_modules.
    if [[ -f "$DIR/website/package-lock.json" ]]; then
        (cd "$DIR/website" && npm install --package-lock-only --ignore-scripts >/dev/null)
        echo "  website/package-lock.json -> $_version"
    fi
fi

# CHANGELOG: stamp [Unreleased] with the new version and date (skip when the
# version heading already exists, so re-runs stay idempotent).
_changelog="$DIR/CHANGELOG.md"
if grep -q "^## \[$_version\]" "$_changelog"; then
    echo "note: CHANGELOG.md already has [$_version]; left unchanged"
elif grep -q '^## \[Unreleased\]' "$_changelog"; then
    # Insert the new release heading immediately after the [Unreleased] heading.
    sed -i "0,/^## \[Unreleased\]$/s//## [Unreleased]\n\n## [$_version] - $_date/" "$_changelog"
    echo "  CHANGELOG.md -> [$_version] - $_date"
else
    echo "warning: no [Unreleased] section in CHANGELOG.md; left unchanged" >&2
fi

echo
echo "Bumped to $_version."
echo
echo "Next steps:"
echo "  git diff Cargo.toml Cargo.lock website/package.json website/package-lock.json CHANGELOG.md"
echo "  git add Cargo.toml Cargo.lock website/package.json website/package-lock.json CHANGELOG.md"
echo "  git commit -m \"chore: bump version to $_version\""
echo "  git tag v$_version"
echo "  git push origin main --tags   # CI publishes wheels from the tag"

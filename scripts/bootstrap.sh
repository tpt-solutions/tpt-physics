#!/usr/bin/env bash
# bootstrap.sh — verify the local sibling dependency workspaces exist.
#
# tpt-physics consumes tpt-math and tpt-fem through Cargo path dependencies
# (see the [workspace.dependencies] table in the root Cargo.toml). Those
# crates are expected to live in sibling directories. This script checks they
# are present and points you at the override env vars if they are not.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

missing=0
for dep in tpt-math tpt-fem; do
    if [ ! -d "$ROOT/../$dep" ]; then
        echo "ERROR: sibling workspace '../$dep' not found." >&2
        echo "       Clone it next to this repo, or set ${dep^^}_PATH to its location." >&2
        missing=1
    else
        echo "OK: found ../$dep"
    fi
done

if [ "$missing" -ne 0 ]; then
    exit 1
fi

echo ""
echo "All sibling dependencies present. Next steps:"
echo "  just setup   # re-run this check"
echo "  just test    # build & run the test suite"
echo "  just bench   # run the criterion benchmarks"

#!/usr/bin/env bash
# Bootstrap the tpt-physics workspace: verify sibling crates are present and
# that the Rust toolchain is available, then run a full check.
#
# tpt-physics depends on sibling crates tpt-math and tpt-fem via path
# dependencies. Clone them next to this repo:
#
#   $HOME/src/tpt-math
#   $HOME/src/tpt-fem
#   $HOME/src/tpt-physics   <-- this repo
#
# Usage:  ./scripts/bootstrap.sh            # check + cargo check
#         ./scripts/bootstrap.sh --build    # also run cargo build
#         ./scripts/bootstrap.sh --test     # also run cargo test

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SIBLINGS=(tpt-math tpt-fem)

for s in "${SIBLINGS[@]}"; do
  if [ -d "$ROOT/../$s" ]; then
    echo "[ok] found sibling: $s"
  else
    echo "error: missing sibling crate '$s' at $ROOT/../$s" >&2
    echo "clone it next to this repo (see header comment)" >&2
    exit 1
  fi
done

if ! command -v cargo >/dev/null 2>&1; then
  echo "error: cargo not found on PATH. Install Rust from https://rustup.rs (>= 1.84)." >&2
  exit 1
fi
echo "[ok] $(cargo --version)"

echo "==> cargo check --workspace"
cargo check --workspace
if [ "${1:-}" = "--build" ]; then
  echo "==> cargo build --workspace"
  cargo build --workspace
fi
if [ "${1:-}" = "--test" ]; then
  echo "==> cargo test --workspace"
  cargo test --workspace
fi
echo "Bootstrap complete."

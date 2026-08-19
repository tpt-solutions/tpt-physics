#!/usr/bin/env bash
# run_gallery.sh — run every example in the tpt-physics gallery.
set -euo pipefail

run() {
    local crate="$1" example="$2"
    echo "================================================================"
    echo ">>> $crate :: $example"
    echo "================================================================"
    cargo run -q --release --example "$example" -p "$crate" 2>&1
    echo ""
}

run tpt-phys-cfd cavity
run tpt-phys-dem granular_pile
run tpt-phys-orchestrator rl_pendulum
run tpt-phys-core uq_cantilever --features uq

echo "Gallery complete."

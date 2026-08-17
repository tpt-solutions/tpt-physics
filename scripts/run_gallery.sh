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

run tpt-physics-fea beam
run tpt-physics-cfd cavity
run tpt-physics-dem granular_pile
run tpt-physics-ai rl_pendulum
run tpt-physics-fea pile_cage_spacer
run tpt-physics-fea spacer_benchmark

echo "Gallery complete."

#!/usr/bin/env bash
# Build the tpt-physics WebAssembly web playground. See build_wasm.ps1 for details.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WASM_CRATE="$ROOT/crates/tpt-physics-wasm"
OUT_DIR="$WASM_CRATE/www/pkg"

if command -v wasm-pack >/dev/null 2>&1; then
  echo "==> wasm-pack build (target web)"
  ( cd "$WASM_CRATE" && wasm-pack build --target web --out-dir www/pkg --out-name tpt_physics_wasm )
elif command -v wasm-bindgen >/dev/null 2>&1; then
  echo "==> cargo build + wasm-bindgen CLI"
  cargo build -p tpt-physics-wasm --target wasm32-unknown-unknown --release
  WASM="$ROOT/target/wasm32-unknown-unknown/release/tpt_physics_wasm.wasm"
  mkdir -p "$OUT_DIR"
  wasm-bindgen "$WASM" --out-dir "$OUT_DIR" --out-name tpt_physics_wasm --target web
else
  echo "error: neither wasm-pack nor wasm-bindgen-cli found on PATH" >&2
  exit 1
fi

echo "==> bindings written to $OUT_DIR"
echo "==> serve the playground:"
echo "   cd $WASM_CRATE/www && python -m http.server 8080   # open http://localhost:8080"

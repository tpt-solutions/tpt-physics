# Build the tpt-physics WebAssembly web playground.
#
# Produces the JS/wasm bindings consumed by crates/tpt-physics-wasm/www/playground.js
# (imports "./tpt_physics_wasm.js", generated into www/pkg/).
#
# Requires one of:
#   * wasm-pack  (https://rustwasm.github.io/wasm-pack/)
#   * cargo + wasm-bindgen-cli
#
# The wasm32-unknown-unknown target must be installed:
#   rustup target add wasm32-unknown-unknown

$ErrorActionPreference = "Stop"

$Root = Resolve-Path (Join-Path $PSScriptRoot "..")
$WasmCrate = Join-Path $Root "crates/tpt-physics-wasm"
$OutDir = Join-Path $WasmCrate "www/pkg"

if (Get-Command wasm-pack -ErrorAction SilentlyContinue) {
    Write-Host "==> wasm-pack build (target web)"
    Push-Location $WasmCrate
    try {
        wasm-pack build --target web --out-dir www/pkg --out-name tpt_physics_wasm
    } finally {
        Pop-Location
    }
} elseif (Get-Command wasm-bindgen -ErrorAction SilentlyContinue) {
    Write-Host "==> cargo build + wasm-bindgen CLI"
    cargo build -p tpt-physics-wasm --target wasm32-unknown-unknown --release
    $Wasm = Join-Path $Root "target/wasm32-unknown-unknown/release/tpt_physics_wasm.wasm"
    New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
    wasm-bindgen "$Wasm" --out-dir $OutDir --out-name tpt_physics_wasm --target web
} else {
    Write-Error "Neither wasm-pack nor wasm-bindgen-cli found on PATH. Install one (see rustwasm docs)."
    exit 1
}

Write-Host "==> bindings written to $OutDir"
Write-Host "==> serve the playground:"
Write-Host "   cd $WasmCrate/www; python -m http.server 8080   # then open http://localhost:8080"

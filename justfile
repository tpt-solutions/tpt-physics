# justfile — convenient shortcuts for the tpt-physics workspace.
#
# Run `just` to see the list, e.g. `just check`, `just test`, `just bench`.

# List available recipes.
default:
    @just --list

# Format the whole workspace (and check in CI).
fmt:
    cargo fmt --all

# Verify formatting is clean (CI gate).
fmt-check:
    cargo fmt --all -- --check

# Clippy with warnings denied (CI gate).
clint:
    cargo clippy --workspace --all-targets -- -D warnings

# Deny copyleft / unexpected licenses (CI gate).
deny:
    cargo deny check

# Build every crate.
build:
    cargo build --workspace

# Run the full test suite.
test:
    cargo test --workspace

# Run all criterion benchmarks.
bench:
    cargo bench --workspace

# Seed sibling dependencies and check the workspace.
setup:
    pwsh scripts/bootstrap.ps1

# End-to-end check used by CI.
ci: fmt-check clint deny build test

# Build the WebAssembly web playground (DEM + CFD) and its JS bindings.
# Output lands in `crates/tpt-physics-wasm/www/pkg`; serve `www/` with any
# static file server (e.g. `python -m http.server`) to use the playground.
# Requires `wasm-pack` (or `cargo` + `wasm-bindgen-cli`) on PATH.
wasm:
    pwsh scripts/build_wasm.ps1

# Serve the WebGL playground locally (after `just wasm`).
serve-wasm:
    pwsh -Command "cd crates/tpt-physics-wasm/www; python -m http.server 8080"

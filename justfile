# tpt-physics task runner.
#
# Recipes for the common developer workflows. Requires `just`
# (https://github.com/casey/just) and a Rust toolchain.

# List available recipes.
default:
    @just --list

# Verify the sibling dependency workspaces (tpt-math, tpt-fem) are present.
setup:
    ./scripts/bootstrap.sh

# Type-check the whole workspace.
check:
    cargo check --workspace --all-targets

# Apply rustfmt.
fmt:
    cargo fmt --all

# Lint with clippy, denying warnings (mirrors CI).
lint:
    cargo clippy --workspace --all-targets -- -D warnings

# Run the full test suite.
test:
    cargo test --workspace

# Run only the fast (non-ignored) tests.
test-fast:
    cargo test --workspace -- --skip large_scale

# Run the criterion benchmark harnesses.
bench:
    cargo bench --workspace

# Build the API documentation.
doc:
    cargo doc --no-deps --workspace

# Build every example.
examples:
    cargo build --workspace --examples

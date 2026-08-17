#!/usr/bin/env pwsh
# Bootstrap the tpt-physics workspace: verify sibling crates are present and
# that the Rust toolchain is available, then run a full check.
#
# tpt-physics depends on sibling crates tpt-math and tpt-fem via path
# dependencies (see the workspace [workspace.dependencies] table). Those repos
# must be cloned next to this one:
#
#   C:\Programming\tpt-math
#   C:\Programming\tpt-fem
#   C:\Programming\tpt-physics   <-- this repo
#
# Usage:  .\scripts\bootstrap.ps1            # check + cargo check
#         .\scripts\bootstrap.ps1 -Build     # also run cargo build
#         .\scripts\bootstrap.ps1 -Test      # also run cargo test

[CmdletBinding()]
param(
    [switch]$Build,
    [switch]$Test
)

$ErrorActionPreference = "Stop"

$siblings = @("tpt-math", "tpt-fem")
$root = Resolve-Path (Join-Path $PSScriptRoot "..")

foreach ($s in $siblings) {
    $p = Join-Path $root $s
    if (-not (Test-Path $p)) {
        Write-Error "Missing sibling crate '$s' at $p. Clone it next to this repo (see header)."
    } else {
        Write-Host "[ok] found sibling: $s"
    }
}

if (-not (Get-Command cargo -ErrorAction SilentlyContinue)) {
    Write-Error "cargo not found on PATH. Install Rust from https://rustup.rs (>= 1.84)."
}
Write-Host "[ok] cargo $(cargo --version)"

Write-Host "==> cargo check --workspace"
cargo check --workspace
if ($Build) {
    Write-Host "==> cargo build --workspace"
    cargo build --workspace
}
if ($Test) {
    Write-Host "==> cargo test --workspace"
    cargo test --workspace
}
Write-Host "Bootstrap complete."

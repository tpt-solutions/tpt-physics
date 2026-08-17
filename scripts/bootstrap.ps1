# bootstrap.ps1 — verify the local sibling dependency workspaces exist.
#
# tpt-physics consumes tpt-math and tpt-fem through Cargo path dependencies
# (see the [workspace.dependencies] table in the root Cargo.toml). Those
# crates are expected to live in sibling directories. This script checks they
# are present and points you at the override env vars if they are not.
$ErrorActionPreference = "Stop"

$Root = Resolve-Path (Join-Path $PSScriptRoot "..")

$missing = $false
foreach ($dep in @("tpt-math", "tpt-fem")) {
    $path = Join-Path $Root $dep
    if (-not (Test-Path $path)) {
        Write-Error "sibling workspace '$path' not found. Clone it next to this repo, or set $($dep.ToUpper())_PATH to its location."
        $missing = $true
    } else {
        Write-Host "OK: found $path"
    }
}

if ($missing) {
    exit 1
}

Write-Host ""
Write-Host "All sibling dependencies present. Next steps:"
Write-Host "  just setup   # re-run this check"
Write-Host "  just test    # build & run the test suite"
Write-Host "  just bench   # run the criterion benchmarks"

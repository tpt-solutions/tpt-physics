# run_gallery.ps1 — run every example in the tpt-physics gallery.
$ErrorActionPreference = "Continue"

function Run-Example {
    param([string]$Crate, [string]$Example, [string]$Features = "")
    Write-Host "================================================================" -ForegroundColor Cyan
    Write-Host ">>> $Crate :: $Example" -ForegroundColor Cyan
    Write-Host "================================================================"
    if ($Features -eq "") {
        cargo run -q --release --example $Example -p $Crate 2>&1
    } else {
        cargo run -q --release --example $Example -p $Crate --features $Features 2>&1
    }
    Write-Host ""
}

Run-Example tpt-phys-cfd cavity
Run-Example tpt-phys-dem granular_pile
Run-Example tpt-phys-orchestrator rl_pendulum
Run-Example tpt-phys-core uq_cantilever uq

Write-Host "Gallery complete."

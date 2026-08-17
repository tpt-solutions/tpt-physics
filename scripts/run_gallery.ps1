# run_gallery.ps1 — run every example in the tpt-physics gallery.
$ErrorActionPreference = "Continue"

function Run-Example {
    param([string]$Crate, [string]$Example)
    Write-Host "================================================================" -ForegroundColor Cyan
    Write-Host ">>> $Crate :: $Example" -ForegroundColor Cyan
    Write-Host "================================================================"
    cargo run -q --release --example $Example -p $Crate 2>&1
    Write-Host ""
}

Run-Example tpt-physics-fea beam
Run-Example tpt-physics-cfd cavity
Run-Example tpt-physics-dem granular_pile
Run-Example tpt-physics-ai rl_pendulum
Run-Example tpt-physics-fea pile_cage_spacer
Run-Example tpt-physics-fea spacer_benchmark

Write-Host "Gallery complete."

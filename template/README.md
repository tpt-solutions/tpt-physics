# {{project-name}}

A minimal [`tpt-physics`](https://github.com/tpt-solutions/tpt-physics)
application scaffolded with `cargo generate`.

It drops a self-weight granular pile under gravity using the material database
(`tpt-phys-core`) and the discrete-element `World` (`tpt-phys-dem`).

## Build

```sh
cargo run
```

Requires the sibling `tpt-math`, `tpt-fem`, and `tpt-science` workspaces to be
checked out alongside `tpt-physics` (the generated `Cargo.toml` references them
by relative path).

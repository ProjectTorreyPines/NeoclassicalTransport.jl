# NeoclassicalTransport.jl

Neoclassical transport for IMAS data (used by FUSE):

- **Hirshman-Sigmar** and **Chang-Hinton** analytic models
- **FACIT** impurity transport
- `run_neo`: calls a locally installed Fortran NEO executable

## Analytic models

```julia
geom = NeoclassicalTransport.get_equilibrium_geometry(eqt, cp1d)
prof = NeoclassicalTransport.get_plasma_profiles(eqt, cp1d)
sol  = NeoclassicalTransport.hirshmansigmar(gridpoint, eqt, cp1d, prof, geom)   # GACODE.FluxSolution
sol  = NeoclassicalTransport.changhinton(eqt, cp1d, rho, 1)                      # ion energy flux only
```

## Fortran NEO (`run_neo`)

`run_neo(ineo)` writes `input.neo`, runs the GACODE NEO executable and returns a
`GACODE.FluxSolution`; it is kept for cross-checks and is not differentiable. It needs
GACODE installed locally: https://fuse.help/install.html#Install-GACODE (see the note after
step 6 — on macOS you may need to replace `mpif90-openmpi-mp` with `mpif90-openmpi-gcc12`
in `$GACODE_ROOT/platform/build/make.inc.OSX_MONTEREY` and `mpirun-openmpi-mp` with
`mpirun-openmpi-gcc12` in `$GACODE_ROOT/platform/exec/exec.OSX_MONTEREY`). On login nodes
where the `neo -e` wrapper cannot launch, source your GACODE environment and call
`NeoclassicalTransport.use_serial_neo!()`, which builds a serial no-MPI NEO from
`utilities/serial_neo/` and points `run_neo` at it.

## Online documentation
For more details, see the [online documentation](https://projecttorreypines.github.io/NeoclassicalTransport.jl/dev).

![Docs](https://github.com/ProjectTorreyPines/NeoclassicalTransport.jl/actions/workflows/make_docs.yml/badge.svg)

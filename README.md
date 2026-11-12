# NeoclassicalTransport.jl

Neoclassical transport for IMAS data (used by FUSE):

- **Hirshman-Sigmar** and **Chang-Hinton** analytic models
- **FACIT** impurity transport
- **yancc**: drift kinetic solver with the linearized Fokker-Planck collision operator, through
  [Yancc.jl](https://github.com/ProjectTorreyPines/Yancc.jl) (see below)
- `run_neo`: calls a locally installed Fortran NEO executable (see below)

## Analytic models

```julia
geom = NeoclassicalTransport.get_equilibrium_geometry(eqt, cp1d)
prof = NeoclassicalTransport.get_plasma_profiles(eqt, cp1d)
sol  = NeoclassicalTransport.hirshmansigmar(gridpoint, eqt, cp1d, prof, geom)   # GACODE.FluxSolution
sol  = NeoclassicalTransport.changhinton(eqt, cp1d, rho, 1)                      # ion energy flux only
```

## Drift kinetic solver (`yancc`)

```julia
gridpoints = [argmin(abs.(cp1d.grid.rho_tor_norm .- rho)) for rho in 0.25:0.1:0.85]
geometry = NeoclassicalTransport.get_yancc_geometry(eqt, cp1d, gridpoints)         # Miller flux surfaces
caches = [NeoclassicalTransport.YanccCache() for _ in gridpoints]                  # optional
sol = NeoclassicalTransport.yancc(gridpoints[k], eqt, cp1d, geometry.fields[k]; geometry.a, cache=caches[k])   # GACODE.FluxSolution
```

Electron and ion energy and particle fluxes from
[Yancc.jl](https://github.com/ProjectTorreyPines/Yancc.jl), the Julia port of
[yancc](https://github.com/f0uriest/yancc) by Rory Conlin, on the Miller parametrization of
the flux surface and without rotation or radial electric field. The species are those of
`InputNEO` and `hirshmansigmar`: every ion with its total density (thermal + fast) and its
own temperature, and the electrons. A cache kept from one call to
the next at the same grid point reuses the preconditioner and the last solution, as between
the iterations of a transport solver. Keywords `nx`, `nalpha` and `ntheta` (of
`get_yancc_geometry`) set the resolution; the defaults (6, 65, 41) are within a few percent of
the converged energy fluxes at moderate collisionality and about 10 % above them in the
collisionless core of a reactor.

If you use this model, please cite R. Conlin and M. Landreman, "yancc: a GPU-accelerated,
differentiable solver for neoclassical transport in tokamaks and stellarators",
*Plasma Phys. Control. Fusion* **68**, 115001 (2026),
[doi:10.1088/1361-6587/aeb628](https://doi.org/10.1088/1361-6587/aeb628).

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

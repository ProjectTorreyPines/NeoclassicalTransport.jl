# NeoclassicalTransport.jl

Neoclassical transport for IMAS data (used by FUSE):

- **NEO**: a native Julia port of the GACODE drift-kinetic solver, no Fortran executable needed, differentiable with ForwardDiff
- **NEO-NN**: neural-network surrogates of NEO
- **Hirshman-Sigmar** and **Chang-Hinton** analytic models
- `run_neo`: calls a locally installed Fortran NEO executable

## NEO

The port covers NEO's local, axisymmetric solve: collision models 1–5 (default 4,
full linearized Fokker-Planck), s-alpha / large-aspect-ratio / Miller geometry, strong
rotation, adiabatic or kinetic electrons. Results match Fortran NEO to ~1e-8.

```julia
using NeoclassicalTransport
ineo = NeoclassicalTransport.InputNEO(eqt, cp1d, gridpoint)   # from an IMAS dd (same input as run_neo)
sol  = NeoclassicalTransport.run_neo_native(ineo)              # -> NEOSolution
sol.pflux, sol.eflux, sol.mflux        # per-species fluxes, NEO normalizations
sol.jpar, sol.uparB, sol.vpol_th0      # bootstrap current, parallel flows, poloidal velocity
GACODE.FluxSolution(sol)               # the lumped gyroBohm fluxes run_neo returns

sols = NeoclassicalTransport.run_neo_native(ineos)   # a batch, threaded over the inputs (julia -t N)
p    = NeoclassicalTransport.NEOParams("input.neo")  # or from an input.neo file
sol  = NeoclassicalTransport.solve_neo(p)
```

Options: `ineo.COLLISION_MODEL = 1..5`; `solve_neo(p; serial=true)` runs plain loops
(bit-identical, for debugging); `keep_g=true` keeps the distribution function.
Unsupported NEO features (`profile_model=2`, 3D, Spitzer, anisotropic species) raise an
explicit error.

A solve is dominated by the sparse factorization (~0.7 s for 5 species on the default
grid, ~0.2 s for 3). When the same surfaces are solved repeatedly, keep caches to reuse
the factorizations; a repeated 5-species point then costs ~0.05 s:

```julia
caches = [NeoclassicalTransport.NEOFactorCache(; refine=true) for _ in ineos]
sols = NeoclassicalTransport.run_neo_native(ineos; caches)   # call again as the profiles evolve
```

Differentiation: an `InputNEO` built from a Dual-valued `dd` (FUSE's
`jacobian_method=:forward_ad`), or `InputNEO{D}(ineo)` with Duals on its fields, gives
a `NEOSolution{Dual}`; the sparse solve uses an implicit-function rule, so derivatives
cost one triangular solve per partial on top of the Float64 factorization.

Verification: `test/runtests_neo_native.jl` compares every intermediate array, the
assembled system, the solution vector and all transport outputs against Fortran NEO
reference data in `test/neo_reference/` (regenerate with `test/neo_reference/generate.sh`
and a gacode installation). `NEO_NATIVE_FULL=1` adds gacode's regression cases reg01–15.
The one known discrepancy is NEO's own: its double-precision recursion for the
field-particle integrals is unstable for `0.1 <= (vth_a/vth_b)^2 <= 10`, where Fortran
and Julia agree only to ~1e-6; `solve_neo(p; fcoll_exact=true)` evaluates those integrals
in extended precision.

## Analytic models

```julia
geom = NeoclassicalTransport.get_equilibrium_geometry(eqt, cp1d)
prof = NeoclassicalTransport.get_plasma_profiles(eqt, cp1d)
sol  = NeoclassicalTransport.hirshmansigmar(gridpoint, eqt, cp1d, prof, geom)   # GACODE.FluxSolution
sol  = NeoclassicalTransport.changhinton(eqt, cp1d, rho, 1)                      # ion energy flux only
```

## NEO-NN

Ensemble neural-network surrogates (20 members each) trained on Fokker-Planck NEO
databases, shipped in `models/` via Git LFS:

| model | devices | outputs |
|---|---|---|
| `neonn_d3d+mastu+nstx_flux` (default) | DIII-D + MAST-U + NSTX | 9 gyroBohm fluxes (p/e/m × ion1/ion2/elec) |
| `neonn_d3d+mastu+nstx_flow` (default) | DIII-D + MAST-U + NSTX | vpol × 3 species + jpar |
| `neonn_d3d_{flux,flow}` | DIII-D only | same |
| `neonn_mastu+nstx_{flux,flow}` | MAST-U + NSTX | same |
| `neonn_d3dedge_{flux,flow}` | DIII-D, edge radii (rho 0.80-0.99) | same |
| `neonn_d3dnearedge_{flux,flow}` | DIII-D, near-edge radii (rho 0.68-0.94) | same |
| `neonn_d3dnegdedge_{flux,flow}` | DIII-D, negative triangularity, edge radii (rho 0.80-0.99) | same |
| `neonn_d3d_withnegD_{flux,flow}` | DIII-D, ± triangularity, core radii (rho 0.10-0.90); blends radially with the two nets below | same |
| `neonn_d3dnearedge_withnegD_{flux,flow}` | DIII-D, ± triangularity, near-edge radii (rho 0.68-0.94) | same |
| `neonn_d3dedge_withnegD_{flux,flow}` | DIII-D, ± triangularity, edge radii (rho 0.80-0.99) | same |
| `neonn_mastu+nstx_withnegD_{flux,flow}` | MAST-U + NSTX, ± triangularity, core radii (rho 0.10-0.90); blends radially with the two nets below | same |
| `neonn_mastunearedge+nstxnearedge_withnegD_{flux,flow}` | MAST-U + NSTX, ± triangularity, near-edge radii (rho 0.68-0.94) | same |
| `neonn_mastuedge+nstxedge_withnegD_{flux,flow}` | MAST-U + NSTX, ± triangularity, edge radii (rho 0.80-0.99) | same |

```julia
using NeoclassicalTransport
input_neos = [NeoclassicalTransport.InputNEO(eqt, cp1d, gp) for gp in gridpoints]
sols = NeoclassicalTransport.run_neonn(input_neos)               # Vector{GACODE.FluxSolution}, gyroBohm units
flows = NeoclassicalTransport.run_neonn_flow(input_neos)         # Vector{NEOFlowSolution}, bulk-ion v_norm units
sols_u = NeoclassicalTransport.run_neonn(input_neos; uncertain=true)  # ensemble mean ± std (Measurements.jl)
```

Fluxes are in the same tgyro-block gyroBohm normalization as `run_neo`, so the two are
drop-in comparable. Flow quantities are in NEO bulk-ion normalized units
(`vpol` in `sqrt(k*T_1/m_1)`, `jpar` in `e*n_1*sqrt(k*T_1/m_1)`, with `T_1`/`n_1` the
bulk-ion values). Model selection is by `model_filename`; see
`NeoclassicalTransport.available_models()`.

### Radial electric field from the flow nets

`neoclassical_Er` closes radial force balance for one species with the network's
poloidal flow, returning the outboard-midplane `E_R` [V/m] and its three terms:

```julia
gps = [argmin(abs.(cp1d.grid.rho_tor_norm .- r)) for r in (0.90, 0.95, 0.99)]
v_tor = [8.0e4, 6.0e4, 3.0e4]   # m/s, measured toroidal rotation of that species
sols = NeoclassicalTransport.neoclassical_Er(eqt, cp1d, gps, v_tor; species=:impurity)
sols[1].Er, sols[1].Er_pressure, sols[1].Er_vtor, sols[1].Er_vpol
```

    E_R = (dp_s/dR)/(Z_s e n_s) - v_φ,s B_Z + v_θ,s B_φ    (at θ = 0)

`v_tor` is required: the toroidal velocity is set by momentum transport and
torque, not by neoclassical theory — in the standard diagnostic application it
is the CER rotation of the same impurity whose `v_θ` the network supplies. The
densities, temperatures, charge and gradients come from the same `InputNEONN`
the network is evaluated on, so the 3-species lumping matches the flow
prediction exactly; `B_p` at θ = 0 is taken from the 2D ψ map (required), since
the 1D `r B_unit/q` shortcut misses `|∇r|` at θ = 0 by 1.5-2.4x on a shaped
equilibrium. `species` is `:bulk`, `:impurity` (default) or `:electron` — with a
*consistent* `v_tor` per species they must give the same `E_R`, so disagreement
between two species is a diagnostic on the supplied rotation. Caveats: the
signed `v_θ` inherits the training helicity convention (see below), the training
NEO runs carried no equilibrium-scale radial electric field so `v_θ` is the
unsqueezed neoclassical flow, and `OMEGA_ROT`/`OMEGA_ROT_DERIV` *are* network
inputs — deriving them from an assumed `E_r` means iterating to a fixed point.

Radial blending: selecting a family's core net (`neonn_d3d_*`, the joint
± triangularity `neonn_d3d_withnegD_*`, or the spherical-tokamak
`neonn_mastu+nstx_withnegD_*`; `_withnegD` = trained on the positive and
negative triangularity DBs jointly) blends radially — points with
`RMIN_OVER_A >= 0.881` are evaluated with the family's near-edge net and points
with `RMIN_OVER_A >= 0.975` with its edge net, mirroring TurbulentTransport.jl's
region switching for the `sat3_em_d3d_azf-1_withnegD` TGLF-NN model.

Caveats:
- The nets use a 3-species reduction: bulk hydrogenic ion, one lumped impurity, electrons.
  Extra hydrogenic ions (e.g. T in DT) are folded into the bulk via quasineutrality —
  outside the D+C training distribution, watch the extrapolation warnings
  (`warn_nn_train_bounds=true`, default).
- `PARTICLE_FLUX_i` has length 2 (`[bulk, lumped impurity]`), unlike `run_neo`'s one
  entry per plasma ion.
- Signed outputs (momentum flux, vpol, jpar) follow the training devices' helicity
  convention (IPCCW/BTCCW are not net inputs).
- **Differentiable end to end.** `InputNEO{T}` carries the element type, so a
  `ForwardDiff.Dual` set on a plasma quantity propagates through the species
  lumping, the electron → bulk-ion normalization, the log10 feature transform and
  the ensemble forward pass — `run_neonn` and `run_neonn_flow`, radial-family
  blending included. `InputNEO(eqt, cp1d, gp)` inherits the element type of the
  dd, so a dd carrying Duals (FUSE's AD path) needs nothing extra; otherwise
  start from `InputNEO{D}(input_neo)`:

  ```julia
  ineo0 = NeoclassicalTransport.InputNEO(eqt, cp1d, gridpoint)
  function Qi(x::AbstractVector{D}) where {D<:Real}
      ineo = NeoclassicalTransport.InputNEO{D}(ineo0)
      ineo.DLNTDR_1, ineo.TEMP_1 = x[1], x[2]
      return NeoclassicalTransport.run_neonn(ineo; warn_nn_train_bounds=false).ENERGY_FLUX_i
  end
  ForwardDiff.gradient(Qi, [ineo0.DLNTDR_1, ineo0.TEMP_1])
  ```

  `neoclassical_Er` differentiates w.r.t. the profiles (Dual `cp1d`) at frozen
  geometry — `eqt` must stay Float64, because IMAS's 2D ψ interpolant is not
  Dual-capable. Derivatives that are structurally zero are real: `DLNNDR_1` is
  one, because the bulk-ion density gradient feature is rebuilt from the electron
  and impurity gradients by quasineutrality rather than read. Two entry points are
  not differentiable and say so: `run_neo` (shells out to NEO) and
  `uncertain=true` (the ensemble spread).
- Validation against NEO on a machine with GACODE installed:
  compare `run_neonn(input_neo)` vs `run_neo(input_neo)` at a few radii.
  (Cross-checked on actual training inputs: ≤1.5% per channel at mid-radius
  against full Fokker-Planck NEO.) On login nodes where the `neo -e` wrapper
  cannot launch, see `use_serial_neo!()` below; the comparison notebook calls it
  automatically when it detects `GACODE_ROOT`.

Notebooks:
- `examples/run_NEONN.ipynb` — run NEO-NN and compare against full NEO,
  Hirshman-Sigmar and Chang-Hinton on the same plasma, with radial profiles and
  ensemble uncertainty bands.
- `utilities/convert_nn.ipynb` — export the BSON ensembles to ONNX. Each model
  has a `models/<name>/` directory (Git LFS) with one `.onnx` per ensemble member
  plus `xnames/ynames/xm/xsigma/ym/ysigma/xbounds_*` sidecar text files, for
  consumers outside Julia (the ONNX graphs are the raw networks — apply the
  log10/standardization input pipeline and output de-standardization yourself;
  batch-first `[N, features]` layout).

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

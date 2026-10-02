# NeoclassicalTransport.jl

Calls the drift-kinetic solver NEO for high-accuracy neoclassical calculations.
It also implements Chang-Hinton and Hirshman-Sigmar neoclassical calculations,
and ships NEO-NN neural-network surrogates of NEO (no NEO executable needed).

NOTE: Running NEO requires GACODE executables to be locally installed. 

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
  cannot launch (srun/mpirun dispatch), source your GACODE environment and call
  `NeoclassicalTransport.use_serial_neo!()` — this builds the serial no-MPI NEO
  from `utilities/serial_neo/` against your gacode tree (once) and points
  `run_neo` at it via `NEO_EXECUTABLE`. The comparison notebook does this
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

Link to instructions on GACODE installation: https://fuse.help/install.html#Install-GACODE

See the note following step 6 - you may need to replace `mpif90-openmpi-mp` with `mpif90-openmpi-gcc12`
in the platform-specific make file found in `$GACODE_ROOT/platform/build/make.inc.OSX_MONTEREY` and
`mpirun-openmpi-mp` with `mpirun-openmpi-gcc12` in the platform exec file found in
`$GACODE_ROOT/platform/exec/exec.OSX_MONTEREY`.

## Native NEO (Julia port)

`src/neo/` is a Julia port of the local (`profile_model=1`), axisymmetric
drift-kinetic solve of GACODE NEO with the full linearized Fokker-Planck
collision operator (`collision_model=4`): the Laguerre/Legendre energy basis,
the closed-form field-particle integrals (`neo_compute_fcoll`), s-alpha /
large-aspect-ratio / Miller (MXH) geometry, strong rotation
(`rotation_model=2`, poloidal density variation from quasi-neutrality), the
sparse kinetic matrix, the UMFPACK solve, and every transport moment NEO
writes (fluxes, gyroviscous fluxes, bootstrap current, parallel flows,
poloidal/toroidal velocities). It runs in process, so there is no file I/O,
no external binary and no MPI:

```julia
ineo = NeoclassicalTransport.InputNEO(eqt, cp1d, gridpoint)   # same input as run_neo
sol  = NeoclassicalTransport.run_neo_native(ineo)              # -> NEOSolution
sol.pflux, sol.eflux, sol.jpar, sol.vpol_th0                   # NEO normalizations, per species
GACODE.FluxSolution(sol)                                       # the lumped fluxes run_neo returns

sols = NeoclassicalTransport.run_neo_native(ineos)   # a batch: threaded over the inputs
p    = NeoclassicalTransport.NEOParams("input.neo")  # or straight from an input.neo file
sol  = NeoclassicalTransport.solve_neo(p; serial=true)  # plain loops everywhere (debugging)
```

Parallelism: the species-pair collision build and the row-wise assembly use
`Threads.@threads`; a batch of inputs is split over tasks with BLAS
single-threaded for the duration. `serial=true` gives bit-identical results
on plain loops. Start Julia with `-t N` to use it.

Cost and factorization reuse: the UMFPACK factorization is 85–95 % of a
solve (5 species on the default 6/17/17 grid: 10 710 rows, ~0.66 s for the
factorization, ~0.02 s for everything else at 8 threads). UMFPACK's own
strategy choice is the slow one for these matrices, so from 4 species on the
unsymmetric strategy is used (4 species 1.18 → 0.41 s, 5 species 1.09 →
0.66 s, and a 1e3–1e4 smaller residual). Pass caller-owned caches to reuse
factorizations across calls:

```julia
caches = [NeoclassicalTransport.NEOFactorCache(; refine=true) for _ in ineos]
sols = NeoclassicalTransport.run_neo_native(ineos; caches)   # call again as the profiles evolve
```

A solve with the same matrix values (the `ForwardDiff` passes of a Jacobian
at the point just solved) skips the factorization; with `refine=true` a solve
on a slightly changed matrix (the next flux-matcher evaluation) uses the old
factorization as a GMRES preconditioner, 4–14 iterations for 0.1–10 % changes,
and only refactorizes when that does not converge. Without `caches` each task
uses a temporary cache; `refine` is off by default so results stay
bit-reproducible.

Differentiation: `NEOParams{<:ForwardDiff.Dual}` (e.g. from an `InputNEO`
built from a Dual-valued `dd`) flows through every stage; the sparse solve
applies the implicit-function rule, `g₀ = A₀⁻¹b₀`, `ġ = A₀⁻¹(ḃ − Ȧ g₀)`, on
the Float64 factorization, so a `NEOSolution{Dual}` costs one factorization
plus one triangular solve per partial. FUSE's `ActorFluxMatcher` uses this for
`jacobian_method=:forward_ad` with `model=:neo, neo_backend=:julia`.

Supported: `collision_model` 1–5 (4 is the default full Fokker-Planck),
`equilibrium_model` 0/1/2, `rotation_model` 1/2, `laguerre_method` 1–4,
adiabatic or kinetic electrons, `sim_model` 1/2. Anything else
(`profile_model=2`, 3D, Spitzer, anisotropic species, `coll_uncoupled*`) is
rejected with an explicit error. `utilities/profile_native_neo.jl` prints the
per-stage timings and the UMFPACK variants on the reference cases.

Validation (`test/runtests_neo_native.jl`): every intermediate array
(basis, collision matrices, geometry, rotation), the assembled system
(`A·g_fortran ≈ b`), the solution vector and all transport outputs are
compared against Fortran NEO reference data generated by
`test/neo_reference/generate.sh` (`utilities/serial_neo/neo_dump.f90`, a
full-precision dump driver linked against your gacode `neo_lib.a`; the
element-wise dumps are kept for the small CI cases and reg12, the other
cases are checked on NEO's standard outputs). Two small Miller/rotation
cases run in CI; `NEO_NATIVE_FULL=1` adds the gacode
regression cases reg04/08/12/13/14/15 and a FUSE-built case in all four
`(BTCCW, IPCCW)` sign conventions. Agreement is ~1e-8 in the fluxes, except
where NEO's own double-precision recursion for the negative-index
field-particle integrals is unstable (`0.1 <= (vth_a/vth_b)^2 <= 10`, high
Legendre orders): there Fortran and Julia both carry amplified round-off
and the fluxes agree to ~1e-6 (the same spread as between two Fortran
builds). `solve_neo(p; fcoll_exact=true)` evaluates those integrals in
extended precision instead.

## Online documentation
For more details, see the [online documentation](https://projecttorreypines.github.io/NeoclassicalTransport.jl/dev).

![Docs](https://github.com/ProjectTorreyPines/NeoclassicalTransport.jl/actions/workflows/make_docs.yml/badge.svg)

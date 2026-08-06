# NeoclassicalTransport.jl

Calls the drift-kinetic solver NEO for high-accuracy neoclassical calculations.
It also implements Chang-Hinton and Hirshman-Sigmar neoclassical calculations,
and ships NEO-NN neural-network surrogates of NEO (no NEO executable needed).

NOTE: Running NEO requires GACODE executables to be locally installed. 

## NEO-NN

Ensemble neural-network surrogates (20 members each) trained on Fokker-Planck NEO
databases generated with [runNEOdb](https://github.com/ProjectTorreyPines/runNEOdb),
shipped in `models/` via Git LFS:

| model | devices | outputs |
|---|---|---|
| `neonn_tgyro_d3d+mastu+nstx_flux_v1` (default) | DIII-D + MAST-U + NSTX | 9 gyroBohm fluxes (p/e/m × ion1/ion2/elec) |
| `neonn_tgyro_d3d+mastu+nstx_flow_v1` (default) | DIII-D + MAST-U + NSTX | vpol × 3 species + jpar |
| `neonn_tgyro_d3d_{flux,flow}_v1` | DIII-D only | same |
| `neonn_tgyro_mastu+nstx_{flux,flow}_v1` | MAST-U + NSTX | same |

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

Caveats:
- The nets use a 3-species reduction: bulk hydrogenic ion, one lumped impurity, electrons.
  Extra hydrogenic ions (e.g. T in DT) are folded into the bulk via quasineutrality —
  outside the D+C training distribution, watch the extrapolation warnings
  (`warn_nn_train_bounds=true`, default).
- `PARTICLE_FLUX_i` has length 2 (`[bulk, lumped impurity]`), unlike `run_neo`'s one
  entry per plasma ion.
- Signed outputs (momentum flux, vpol, jpar) follow the training devices' helicity
  convention (IPCCW/BTCCW are not net inputs).
- Validation against NEO on a machine with GACODE installed:
  compare `run_neonn(input_neo)` vs `run_neo(input_neo)` at a few radii.
  (Cross-checked on actual training inputs: ≤1.5% per channel at mid-radius
  against full Fokker-Planck NEO.)

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

## Online documentation
For more details, see the [online documentation](https://projecttorreypines.github.io/NeoclassicalTransport.jl/dev).

![Docs](https://github.com/ProjectTorreyPines/NeoclassicalTransport.jl/actions/workflows/make_docs.yml/badge.svg)

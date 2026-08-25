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
| `neonn_d3d_withnegD_{flux,flow}` | DIII-D, ± triangularity, core radii (rho 0.10-0.90); blends radially with the near-edge/edge nets of this family once all are shipped | same |
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

## Online documentation
For more details, see the [online documentation](https://projecttorreypines.github.io/NeoclassicalTransport.jl/dev).

![Docs](https://github.com/ProjectTorreyPines/NeoclassicalTransport.jl/actions/workflows/make_docs.yml/badge.svg)

# NeoclassicalTransport.jl

Calls the drift-kinetic solver NEO for high-accuracy neoclassical calculations.
It also implements Chang-Hinton and Hirshman-Sigmar neoclassical calculations.

NOTE: Running NEO requires GACODE executables to be locally installed. 

Link to instructions on GACODE installation: https://fuse.help/install.html#Install-GACODE

See the note following step 6 - you may need to replace `mpif90-openmpi-mp` with `mpif90-openmpi-gcc12`
in the platform-specific make file found in `$GACODE_ROOT/platform/build/make.inc.OSX_MONTEREY` and
`mpirun-openmpi-mp` with `mpirun-openmpi-gcc12` in the platform exec file found in
`$GACODE_ROOT/platform/exec/exec.OSX_MONTEREY`.

## Online documentation
For more details, see the [online documentation](https://projecttorreypines.github.io/NeoclassicalTransport.jl/dev).

![Docs](https://github.com/ProjectTorreyPines/NeoclassicalTransport.jl/actions/workflows/make_docs.yml/badge.svg)

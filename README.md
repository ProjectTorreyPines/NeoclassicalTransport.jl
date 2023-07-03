# NEO.jl

Interface between FUSE and NEO, the drift-kinetic solver for high-accuracy neoclassical
calculations from the GA Code suite. 

NOTE: Requires NEO executable to be locally installed. 
Link to instructions on GA code installation: https://fuse.help/install.html#Install-GACODE

See the note following step 6 - you may need to replace mpif90-openmpi-mp with mpif90-openmpi-gcc12
in the platform-specific make file found in $GACODE_ROOT/platform/build/make.inc.OSX_MONTEREY and
mpirun-openmpi-mp with mpirun-openmpi-gcc12 in the platform exec file found in
$GACODE_ROOT/platform/exec/exec.OSX_MONTEREY. 

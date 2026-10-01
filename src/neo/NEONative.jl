"""
    NEONative

Native Julia port of the axisymmetric, local (`profile_model=1`) drift-kinetic
solve of GACODE NEO, with the full linearized Fokker-Planck collision operator
(`collision_model=4`). Everything up to the sparse solve is written to read
line-for-line against the Fortran (`\$GACODE_ROOT/neo/src`), which is what the
reference tests in `test/runtests_neo_native.jl` check against.

Entry points: [`NEOParams`](@ref), [`solve_neo`](@ref), [`run_neo_native`](@ref).
"""
module NEONative

using LinearAlgebra
using SparseArrays
using OffsetArrays
import ForwardDiff
import GACODE
import ..NeoclassicalTransport: InputNEO

include("neo_params.jl")
include("neo_grid.jl")
include("neo_fcoll.jl")
include("neo_collision.jl")
include("neo_geo_miller.jl")
include("neo_geometry.jl")
include("neo_rotation.jl")
include("neo_assemble.jl")
include("neo_solve.jl")
include("neo_transport.jl")
include("neo_driver.jl")

export NEOParams, NEOBasis, NEOCollision, NEOGeometry, NEORotation, NEOSolution, NEOFactorCache
export solve_neo, run_neo_native

end

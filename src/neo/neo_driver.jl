# Drivers: one flux surface (solve_neo), one or a batch of InputNEOs
# (run_neo_native, threaded over the batch), and the lossy conversion to
# GACODE.FluxSolution that FUSE consumes.

"""
    solve_neo(p::NEOParams; serial=false, keep_g=false, cache=nothing) -> NEOSolution

The full native NEO solve of one flux surface: basis, equilibrium, rotation,
collision matrices, assembly, sparse LU, transport moments.

- `serial=true` runs every stage on plain loops (the threaded stages are the
  species-pair collision build and the row-wise assembly; results are
  identical either way).
- `keep_g=true` stores the solution vector in the result.
- `cache::NEOFactorCache` reuses the UMFPACK symbolic analysis between solves
  on the same grid sizes (one cache per task).
- `fcoll_exact=true` builds the field-particle collision integrals in extended
  precision instead of NEO's (unstable) double-precision recursion, see
  [`compute_fcoll`](@ref). Default `false` reproduces Fortran NEO.
"""
function solve_neo(p::NEOParams{T}; serial::Bool=false, keep_g::Bool=false, cache::Union{Nothing,NEOFactorCache}=nothing,
    fcoll_exact::Bool=false) where {T<:Real}
    basis = NEOBasis(p)
    geo = equilibrium(p)
    rot = rotation_phi(p, geo)
    coll = collision_ints(p, basis; serial, fcoll_exact)
    pattern = NEOPattern(p)
    A, b, coef = assemble(p, basis, coll, geo, rot, pattern; serial)
    g = solve_system(A, b, pattern; cache)
    return transport(p, basis, geo, rot, coef, g, pattern; keep_g)
end

"""
    run_neo_native(input_neo::InputNEO; kw...) -> NEOSolution
    run_neo_native(inputs::AbstractVector{<:InputNEO}; serial=false, kw...) -> Vector{NEOSolution}
    run_neo_native(params::AbstractVector{<:NEOParams}; serial=false, kw...) -> Vector{NEOSolution}

Native-Julia replacement for [`run_neo`](@ref): the same `InputNEO`, solved in
process. A vector of inputs (e.g. one per transport grid point) is solved as
a batch, split over threads with BLAS single-threaded for the duration and one
factorization cache per task; `serial=true` solves them one after the other.
Keyword arguments are passed to [`solve_neo`](@ref).
"""
run_neo_native(input_neo::InputNEO; kw...) = solve_neo(NEOParams(input_neo); kw...)
run_neo_native(inputs::AbstractVector{<:InputNEO}; kw...) = run_neo_native([NEOParams(inp) for inp in inputs]; kw...)

function run_neo_native(params::AbstractVector{<:NEOParams}; serial::Bool=false, kw...)
    n = length(params)
    T = isempty(params) ? Float64 : promote_type(map(eltype, params)...)
    out = Vector{NEOSolution{T}}(undef, n)
    nchunks = (serial || Threads.nthreads() == 1) ? 1 : min(n, Threads.nthreads())
    if nchunks <= 1
        let cache = NEOFactorCache()
            for i in 1:n
                out[i] = solve_neo(params[i]; serial, cache, kw...)
            end
        end
        return out
    end
    blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        chunks = [i:nchunks:n for i in 1:nchunks]
        @sync for chunk in chunks
            Threads.@spawn let cache = NEOFactorCache()
                # `let`: the cache must be task-local; a plain assignment here would
                # rebind one shared variable and make every task factorize into the
                # same UmfpackLU
                for i in chunk
                    out[i] = solve_neo(params[i]; serial, cache, kw...)
                end
            end
        end
    finally
        BLAS.set_num_threads(blas_threads)
    end
    return out
end

"""
    GACODE.FluxSolution(sol::NEOSolution)

The lumped flux set FUSE uses, built exactly as [`run_neo`](@ref) builds it
from `out.neo.transport_flux`: GB-normalised tgyro totals, electrons from the
kinetic electron species, ion energy and momentum fluxes summed over ions
(ion particle fluxes kept per species, in input order). Requires kinetic
electrons (`ae_flag=0`).
"""
function GACODE.FluxSolution(sol::NEOSolution{T}) where {T<:Real}
    p = sol.params
    p.ae_flag == 0 || error("NEONative: FluxSolution needs kinetic electrons (ae_flag=0); this solve used adiabatic electrons")
    pflux, eflux, mflux = tgyro_fluxes(sol)
    e = p.is_ele
    ions = [is for is in 1:p.n_species if is != e]
    return GACODE.FluxSolution{T}(eflux[e], sum(eflux[ions]), pflux[e], pflux[ions], sum(mflux[ions]))
end

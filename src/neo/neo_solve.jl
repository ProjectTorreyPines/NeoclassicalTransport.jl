# Sparse direct solve of the kinetic system (replaces SOLVE_sparse/UMFPACK
# of NEO with SparseArrays' UMFPACK). The symbolic analysis is reused across
# solves on the same pattern through a per-task NEOFactorCache.

"""
    NEOFactorCache

Holds the UMFPACK factorization of the previous solve so that the next solve
on the same [`NEOPattern`](@ref) only redoes the numeric factorization
(`lu!`). One cache per task; never share one across threads.
"""
mutable struct NEOFactorCache
    pattern::Union{Nothing,NEOPattern}
    F::Any
end
NEOFactorCache() = NEOFactorCache(nothing, nothing)

"""
    solve_system(A, b, pattern; cache=nothing) -> g

`A \\ b` through UMFPACK. Only `Float64` systems can be factorized; a
`ForwardDiff.Dual`-valued matrix is rejected with an explanation.
"""
function solve_system(A::SparseMatrixCSC{Float64,Int}, b::AbstractVector{Float64}, pattern::NEOPattern; cache::Union{Nothing,NEOFactorCache}=nothing)
    # UMFPACK's dense kernels call BLAS; with OpenBLAS at its default thread
    # count (e.g. 128 on a Perlmutter login node) a 6000-row factorization
    # takes a minute instead of a fraction of a second
    blas_threads = BLAS.get_num_threads()
    blas_threads == 1 || BLAS.set_num_threads(1)
    try
        if cache === nothing
            F = lu(A)
        elseif cache.F !== nothing && cache.pattern === pattern
            F = lu!(cache.F, A)
        else
            F = lu(A)
            cache.F = F
            cache.pattern = pattern
        end
        return F \ b
    finally
        blas_threads == 1 || BLAS.set_num_threads(blas_threads)
    end
end

function solve_system(A::SparseMatrixCSC{T}, b::AbstractVector, pattern::NEOPattern; kw...) where {T<:Real}
    error("NEONative: the sparse LU solve is Float64-only (got element type $T); " *
          "differentiate through the native NEO solve with an implicit-function rule rather than by propagating Duals through UMFPACK")
end

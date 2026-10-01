# Sparse direct solve of the kinetic system (replaces SOLVE_sparse/UMFPACK
# of NEO with SparseArrays' UMFPACK). The factorization is reused across
# solves through a per-task NEOFactorCache: the symbolic analysis whenever the
# pattern is the same, the numeric factorization whenever the matrix values
# are unchanged as well (e.g. the Dual passes of a ForwardDiff Jacobian at the
# same primal point), and optionally a stale factorization as the
# preconditioner of an iterative refinement when the values changed a little
# (successive flux-matcher iterations). ForwardDiff.Dual systems are solved
# with the implicit-function rule on top of the Float64 factorization.

"""
    NEOFactorCache(; refine=false, refine_tol=1e-10, max_refine=20)

Holds the UMFPACK factorization of the previous solve so that the next solve
on the same [`NEOPattern`](@ref) only redoes the numeric factorization
(`lu!`), and a solve with the very same matrix values does not refactorize at
all. One cache per task; never share one across threads.

With `refine=true` a solve whose matrix differs from the factorized one first
tries the old factorization as the preconditioner of a GMRES iteration, and
only refactorizes when that does not reach, within `max_refine` iterations, a
relative residual below `max(refine_tol, 10 × the residual the fresh
factorization achieved)`. Solves on slowly changing matrices (successive
flux-matcher evaluations) then cost a handful of triangular solves instead of
a factorization: 4–14 GMRES iterations for 0.1–10 % parameter changes on the
reference cases. (Plain iterative refinement diverges here: the kinetic matrix
is ill-conditioned enough that a 0.1 % change already gives a Richardson
factor of ~10.) Multi-column right-hand sides (the partials of a
`ForwardDiff` solve) refactorize instead, which is cheaper than a GMRES per
column and makes the following Dual passes at the same point cache hits. Off
by default so results are bit-reproducible.
"""
mutable struct NEOFactorCache
    pattern::Union{Nothing,NEOPattern}
    F::Any
    nzval0::Vector{Float64}   # values of the matrix F was computed from
    refine::Bool
    refine_tol::Float64
    max_refine::Int
    resid0::Float64           # relative residual of the fresh factorization's solve
    n_factorizations::Int     # statistics
    n_refined::Int
end
NEOFactorCache(; refine::Bool=false, refine_tol::Float64=1e-10, max_refine::Int=20) =
    NEOFactorCache(nothing, nothing, Float64[], refine, refine_tol, max_refine, 0.0, 0, 0)

# UMFPACK's dense kernels call BLAS; with OpenBLAS at its default thread
# count (e.g. 128 on a Perlmutter login node) a 6000-row factorization
# takes a minute instead of a fraction of a second
function _blas1(f)
    blas_threads = BLAS.get_num_threads()
    blas_threads == 1 || BLAS.set_num_threads(1)
    try
        return f()
    finally
        blas_threads == 1 || BLAS.set_num_threads(blas_threads)
    end
end

"""
    UMFPACK_UNSYMMETRIC_MIN_SPECIES

UMFPACK's automatic strategy picks the symmetric one for the kinetic matrix
(its pattern is 99.9 % structurally symmetric), which is the slow choice for
the larger systems. Measured on the reference cases (6/17/17 grid, Perlmutter
login node, one BLAS thread), factorization time symmetric → unsymmetric:
3 species 0.18 → 0.21 s, 4 species 1.18 → 0.41 s, 5 species 1.09 → 0.66 s,
with a 1e3–1e4 smaller residual for the unsymmetric strategy. It is selected
from this species count on; `Ref` so it can be changed at run time.
"""
const UMFPACK_UNSYMMETRIC_MIN_SPECIES = Ref(4)

function _umfpack_control(pattern::NEOPattern)
    control = SparseArrays.UMFPACK.get_umfpack_control(Float64, Int)
    if pattern.n_species >= UMFPACK_UNSYMMETRIC_MIN_SPECIES[]
        control[SparseArrays.LibSuiteSparse.UMFPACK_STRATEGY+1] = SparseArrays.LibSuiteSparse.UMFPACK_STRATEGY_UNSYMMETRIC
    end
    return control
end

"""
    _factorize(A::SparseMatrixCSC{Float64,Int}, pattern, cache) -> UmfpackLU

`lu(A)`, or with a cache: the cached factorization when the values are
unchanged, `lu!` (symbolic reuse) when only the values changed, `lu` on a new
pattern.
"""
function _factorize(A::SparseMatrixCSC{Float64,Int}, pattern::NEOPattern, cache::Union{Nothing,NEOFactorCache})
    cache === nothing && return _blas1(() -> lu(A; control=_umfpack_control(pattern)))
    if cache.F !== nothing && cache.pattern === pattern
        if cache.nzval0 == A.nzval
            return cache.F
        end
        F = _blas1(() -> lu!(cache.F, A))
    else
        F = _blas1(() -> lu(A; control=_umfpack_control(pattern)))
        cache.pattern = pattern
    end
    cache.F = F
    cache.nzval0 = copy(A.nzval)
    cache.n_factorizations += 1
    return F
end

_relres(A, X, B) = (nb = norm(B); nb == 0 ? 0.0 : norm(B - A * X) / nb)

# right-preconditioned GMRES with the stale factorization as preconditioner;
# nothing when it does not reach the target within max_refine iterations
function _refine(A::SparseMatrixCSC{Float64,Int}, b::AbstractVector{Float64}, cache::NEOFactorCache)
    F = cache.F
    target = max(cache.refine_tol, 10 * cache.resid0)
    n = length(b)
    β = norm(b)
    β == 0 && return zeros(n)
    m = cache.max_refine
    V = Matrix{Float64}(undef, n, m + 1)
    H = zeros(m + 1, m)
    g = zeros(m + 1)
    cs = zeros(m)
    sn = zeros(m)
    V[:, 1] = b / β
    g[1] = β
    k = 0
    converged = false
    @views for j in 1:m
        k = j
        w = A * _blas1(() -> F \ V[:, j])
        for i in 1:j
            H[i, j] = dot(V[:, i], w)
            w .-= H[i, j] .* V[:, i]
        end
        H[j+1, j] = norm(w)
        if H[j+1, j] != 0
            V[:, j+1] = w / H[j+1, j]
        end   # else: exact breakdown, the rotation below gives g[j+1] = 0 and the loop exits converged
        for i in 1:j-1   # apply the previous Givens rotations to the new column
            t = cs[i] * H[i, j] + sn[i] * H[i+1, j]
            H[i+1, j] = -sn[i] * H[i, j] + cs[i] * H[i+1, j]
            H[i, j] = t
        end
        ρ = hypot(H[j, j], H[j+1, j])
        cs[j] = H[j, j] / ρ
        sn[j] = H[j+1, j] / ρ
        H[j, j] = ρ
        H[j+1, j] = 0.0
        g[j+1] = -sn[j] * g[j]
        g[j] = cs[j] * g[j]
        if abs(g[j+1]) / β < target
            converged = true
            break
        end
    end
    converged || return nothing
    y = UpperTriangular(view(H, 1:k, 1:k)) \ view(g, 1:k)
    x = _blas1(() -> F \ (view(V, :, 1:k) * y))
    # the recurrence estimate can drift from the true residual; check it
    return _relres(A, x, b) <= target ? x : nothing
end

# A X = B for a vector or matrix right-hand side
function _solve(A::SparseMatrixCSC{Float64,Int}, B::AbstractVecOrMat{Float64}, pattern::NEOPattern, cache::Union{Nothing,NEOFactorCache})
    if B isa AbstractVector && cache !== nothing && cache.refine && cache.F !== nothing && cache.pattern === pattern && cache.nzval0 != A.nzval
        X = _refine(A, B, cache)
        if X !== nothing
            cache.n_refined += 1
            return X
        end
    end
    fresh = cache === nothing || cache.F === nothing || cache.pattern !== pattern || cache.nzval0 != A.nzval
    F = _factorize(A, pattern, cache)
    X = _blas1(() -> F \ B)
    if cache !== nothing && cache.refine && fresh
        cache.resid0 = _relres(A, X, B)
    end
    return X
end

"""
    solve_system(A, b, pattern; cache=nothing) -> g

`A \\ b` through UMFPACK. A `ForwardDiff.Dual`-valued system (one level of
Duals over `Float64`) is solved with the implicit-function rule: the values
`g₀ = A₀⁻¹ b₀` and the partials `ġ = A₀⁻¹ (ḃ − Ȧ g₀)` share one factorization
of `A₀`; the partials are one multi-right-hand-side triangular solve.
"""
function solve_system(A::SparseMatrixCSC{Float64,Int}, b::AbstractVector{Float64}, pattern::NEOPattern; cache::Union{Nothing,NEOFactorCache}=nothing)
    return _solve(A, b, pattern, cache)
end

function solve_system(A::SparseMatrixCSC{D,Int}, b::AbstractVector{D}, pattern::NEOPattern;
    cache::Union{Nothing,NEOFactorCache}=nothing) where {Tag,N,D<:ForwardDiff.Dual{Tag,Float64,N}}
    n = size(A, 1)
    colptr, rowval, nzval = A.colptr, A.rowval, A.nzval
    A0 = SparseMatrixCSC(n, n, colptr, rowval, ForwardDiff.value.(nzval))   # shares the cached colptr/rowval
    b0 = ForwardDiff.value.(b)
    g0 = _solve(A0, b0, pattern, cache)

    # R[:, k] = ḃ_k − Ȧ_k g₀
    R = Matrix{Float64}(undef, n, N)
    @inbounds for i in 1:n
        pb = ForwardDiff.partials(b[i])
        for k in 1:N
            R[i, k] = pb[k]
        end
    end
    @inbounds for j in 1:n
        gj = g0[j]
        for idx in colptr[j]:colptr[j+1]-1
            pa = ForwardDiff.partials(nzval[idx])
            i = rowval[idx]
            for k in 1:N
                R[i, k] -= pa[k] * gj
            end
        end
    end
    G = _solve(A0, R, pattern, cache)
    return D[D(g0[i], ForwardDiff.Partials(ntuple(k -> G[i, k], Val(N)))) for i in 1:n]
end

function solve_system(A::SparseMatrixCSC{T}, b::AbstractVector, pattern::NEOPattern; kw...) where {T<:Real}
    error("NEONative: the sparse LU solve supports Float64 and one level of ForwardDiff.Dual{Tag,Float64,N} " *
          "(got element type $T); nested Duals and other number types are not supported")
end

# Per-stage timing of the native NEO solve and of the UMFPACK factorization
# variants, on the reference cases in test/neo_reference.
#
#   julia --project=/path/to/NeoclassicalTransport -t 1 utilities/profile_native_neo.jl [case ...]
#
# Default cases: fuse_miller (5 species, 10710 rows) and reg12 (3 species, 6426 rows).

using LinearAlgebra, SparseArrays, Printf, NeoclassicalTransport
const NN = NeoclassicalTransport.NEONative
const U = SparseArrays.UMFPACK
const L = SparseArrays.LibSuiteSparse   # UMFPACK control/info constants
BLAS.set_num_threads(1)
tmin(f; n=3) = minimum(@elapsed(f()) for _ in 1:n)
const REF = joinpath(dirname(@__DIR__), "test", "neo_reference")

function prof(name)
    p = NN.NEOParams(joinpath(REF, name, "input.neo"))
    println("== $name  ns=$(p.n_species) ne=$(p.n_energy) nxi=$(p.n_xi) nth=$(p.n_theta) n_row=$(NN.n_row(p))  julia threads=$(Threads.nthreads())")
    NN.solve_neo(p)   # warm-up / compile
    basis = NN.NEOBasis(p); geo = NN.equilibrium(p); rot = NN.rotation_phi(p, geo)
    coll = NN.collision_ints(p, basis); pat = NN.NEOPattern(p)
    A, b, coef = NN.assemble(p, basis, coll, geo, rot, pat); F = lu(A); g = F \ b
    NN.transport(p, basis, geo, rot, coef, g, pat)
    key = (p.n_species, p.n_energy, p.n_xi, p.n_theta, Int(NN.constraint_all_ie(p)))
    rows = Any[("NEOBasis (cached)", tmin(() -> NN.NEOBasis(p))), ("equilibrium", tmin(() -> NN.equilibrium(p))),
        ("rotation_phi", tmin(() -> NN.rotation_phi(p, geo))), ("collision_ints", tmin(() -> NN.collision_ints(p, basis))),
        ("NEOPattern (cached)", tmin(() -> NN.NEOPattern(p))), ("_build_pattern cold", tmin(() -> NN._build_pattern(key...))),
        ("assemble", tmin(() -> NN.assemble(p, basis, coll, geo, rot, pat))),
        ("lu cold", tmin(() -> lu(A))), ("lu! warm (symbolic reuse)", tmin(() -> lu!(F, A))), ("F \\ b", tmin(() -> F \ b)),
        ("transport", tmin(() -> NN.transport(p, basis, geo, rot, coef, g, pat))),
        ("solve_neo (no cache)", tmin(() -> NN.solve_neo(p)))]
    c = NN.NEOFactorCache(); NN.solve_neo(p; cache=c)
    push!(rows, ("solve_neo (cache, same values)", tmin(() -> NN.solve_neo(p; cache=c))))
    for (k, v) in rows
        @printf("  %-34s %9.4f s\n", k, v)
    end
    I, J, _ = findnz(A)
    P = SparseMatrixCSC(size(A)..., A.colptr, A.rowval, ones(nnz(A)))
    lnz, unz = U.umf_lunz(F)[1:2]
    @printf("  nnz(A)=%d  nnz/row=%.1f  maxband=%d  structural symmetry=%.3f\n", nnz(A), nnz(A) / size(A, 1),
        maximum(abs.(I .- J)), (2nnz(A) - nnz(P + P')) / nnz(A))
    @printf("  default LU: fill lnz+unz=%d (%.1fx nnz)  flops=%.3g  strategy_used=%g  ordering_used=%g\n", lnz + unz, (lnz + unz) / nnz(A),
        F.info[L.UMFPACK_FLOPS+1], F.info[L.UMFPACK_STRATEGY_USED+1], F.info[L.UMFPACK_ORDERING_USED+1])
    println("  UMFPACK variants (strategy: 0 auto, 1 unsymmetric, 3 symmetric; ordering: 0 cholmod, 1 amd, 3 metis, 4 best):")
    for (nm, s, o) in [("default", nothing, nothing), ("unsym", L.UMFPACK_STRATEGY_UNSYMMETRIC, nothing),
        ("sym", L.UMFPACK_STRATEGY_SYMMETRIC, nothing), ("AMD", nothing, L.UMFPACK_ORDERING_AMD),
        ("METIS", nothing, L.UMFPACK_ORDERING_METIS), ("BEST", nothing, L.UMFPACK_ORDERING_BEST),
        ("sym+AMD", L.UMFPACK_STRATEGY_SYMMETRIC, L.UMFPACK_ORDERING_AMD),
        ("sym+METIS", L.UMFPACK_STRATEGY_SYMMETRIC, L.UMFPACK_ORDERING_METIS)]
        try
            ctl = U.get_umfpack_control(Float64, Int)
            s === nothing || (ctl[L.UMFPACK_STRATEGY+1] = s)
            o === nothing || (ctl[U.JL_UMFPACK_ORDERING] = o)
            Fv = lu(A; control=ctl)
            tc = tmin(() -> lu(A; control=ctl)); tw = tmin(() -> lu!(Fv, A)); ts = tmin(() -> Fv \ b)
            l, u = U.umf_lunz(Fv)[1:2]
            @printf("  %-10s cold=%.4f  warm=%.4f  solve=%.4f  fill=%d  flops=%.3g  strat=%g  ord=%g  relres=%.2e\n", nm, tc, tw, ts,
                l + u, Fv.info[L.UMFPACK_FLOPS+1], Fv.info[L.UMFPACK_STRATEGY_USED+1], Fv.info[L.UMFPACK_ORDERING_USED+1],
                norm(A * (Fv \ b) - b) / norm(b))
        catch e
            println("  $nm ERROR: ", sprint(showerror, e))
        end
    end
    println("  BLAS threads for the default LU:")
    for nt in (1, 2, 4)
        BLAS.set_num_threads(nt)
        @printf("    blas=%d  lu cold=%.4f  lu! warm=%.4f  solve=%.4f\n", nt, tmin(() -> lu(A)), tmin(() -> lu!(F, A)), tmin(() -> F \ b))
    end
    BLAS.set_num_threads(1)
end

cases = isempty(ARGS) ? ["fuse_miller", "reg12"] : ARGS
foreach(prof, cases)

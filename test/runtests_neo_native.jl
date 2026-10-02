# Reference tests of the native NEO port (src/neo/) against Fortran NEO.
#
# The reference data in test/neo_reference/<case>/ is produced by
# test/neo_reference/generate.sh (utilities/serial_neo/neo_dump.f90) and holds
# every intermediate array at full precision (out.neo.dump), the solution
# vector (out.neo.f) and NEO's standard outputs. The `small*` cases always
# run; the gacode regression cases run with NEO_NATIVE_FULL=1.

using NeoclassicalTransport
using NeoclassicalTransport: NEOParams, NEOSolution, NEONative
using NeoclassicalTransport.NEONative: NEOBasis, gamma2, compute_fcoll, collision_ints, collision_ints_mono,
    equilibrium, rotation_phi, NEOPattern, assemble, solve_system, transport, solve_neo, run_neo_native, tgyro_fluxes, E_ALPHA, NEOFactorCache
# stdlibs and GACODE through the package: under Pkg.test only test/Project.toml is on the load path
using NeoclassicalTransport.NEONative.LinearAlgebra
using NeoclassicalTransport.NEONative.SparseArrays
using NeoclassicalTransport.NEONative.OffsetArrays
using Test
const GACODE = NeoclassicalTransport.GACODE
import ForwardDiff

const SF = NeoclassicalTransport.SpecialFunctions
const NEO_REFDIR = joinpath(@__DIR__, "neo_reference")
const NEO_FULL = get(ENV, "NEO_NATIVE_FULL", "") in ("1", "true", "yes")

# ---- readers --------------------------------------------------------------

"out.neo.dump: `# name dims...` headers followed by the values in column-major order"
function read_neo_dump(path)
    d = Dict{String,Any}()
    lines = readlines(path)
    i = 1
    while i <= length(lines)
        hdr = split(lines[i])
        hdr[1] == "#" || error("read_neo_dump: bad header '$(lines[i])'")
        name = String(hdr[2])
        dims = parse.(Int, hdr[3:end])
        n = prod(dims)
        vals = [parse(Float64, lines[i+k]) for k in 1:n]
        d[name] = (length(dims) == 1 && dims[1] == 1) ? vals[1] : reshape(vals, dims...)
        i += n + 1
    end
    return d
end

# 1-element arrays come back as scalars from read_neo_dump (n_species=1)
_vec(x) = x isa Number ? [x] : vec(x)

read_neo_numbers(path) = [parse(Float64, w) for line in eachline(path) if !startswith(strip(line), "#") for w in split(line)]

"out.neo.transport: r, d_phi_sqavg, jpar, vtor_0order_th0, uparB_0order, then 8 values per species"
function read_neo_transport(path, ns)
    v = read_neo_numbers(path)
    length(v) == 5 + 8ns || error("read_neo_transport: expected $(5+8ns) values, got $(length(v))")
    per = reshape(v[6:end], 8, ns)
    return (r=v[1], d_phi_sqavg=v[2], jpar=v[3], vtor_0order_th0=v[4], uparB_0order=v[5],
        pflux=per[1, :], eflux=per[2, :], mflux=per[3, :], uparB=per[4, :],
        klittle_upar=per[5, :], kbig_upar=per[6, :], vpol_th0=per[7, :], vtor_th0=per[8, :])
end

"out.neo.diagnostic_coll{test,field}: per (is,js), per ix, one line of (ne+1)^2 values (ie outer, je inner)"
function read_neo_diagnostic_coll(path, ns, ne, nxi; upper::Bool=false)
    out = OffsetArray(zeros(ns, ns, ne + 1, ne + 1, nxi + 1), 1:ns, 1:ns, 0:ne, 0:ne, 0:nxi)
    lines = filter(l -> !isempty(strip(l)), readlines(path))
    k = 1
    for is in 1:ns, js in 1:ns, ix in 0:nxi
        while !startswith(strip(lines[k]), "# ix")
            k += 1
        end
        m = match(r"ix\s*=\s*(\d+)\s*\(is,js\)=\(\s*(\d+),\s*(\d+)\)", lines[k])
        @assert parse(Int, m[1]) == ix && parse(Int, m[2]) == is && parse(Int, m[3]) == js
        vals = parse.(Float64, split(lines[k+1]))
        n = 0
        for ie in 0:ne, je in (upper ? ie : 0):ne
            n += 1
            out[is, js, ie, je, ix] = vals[n]
        end
        k += 2
    end
    return out
end

relnorm(a, b) = norm(a - b) / max(norm(b), floatmin())
# max over (is,js,ix) blocks of the relative Frobenius error of the (ie,je) matrix
function blockwise_relerr(a, b)
    ns, _, _, _, nx = size(b)
    return maximum(norm(a[is, js, :, :, ix] - b[is, js, :, :, ix]) / max(norm(b[is, js, :, :, ix]), 1e-300) for is in 1:ns, js in 1:ns, ix in 1:nx)
end

# ---- one reference case ---------------------------------------------------

"compare a solution with NEO's standard outputs (out.neo.transport*, out.neo.prec; e16.8 = 8 digits)"
function test_neo_transport_outputs(sol, dir, tol)
    p = sol.params
    ns = p.n_species
    r8 = max(tol.transport, 1e-7)
    tr = read_neo_transport(joinpath(dir, "out.neo.transport"), ns)
    scale = maximum(abs.(tr.eflux))
    @test isapprox(sol.pflux, tr.pflux; rtol=r8, atol=r8 * scale)
    @test isapprox(sol.eflux, tr.eflux; rtol=r8, atol=r8 * scale)
    @test isapprox(sol.mflux, tr.mflux; rtol=r8, atol=r8 * scale)
    @test isapprox(sol.jpar, tr.jpar; rtol=r8, atol=r8 * scale)
    @test isapprox(sol.uparB, tr.uparB; rtol=r8, atol=r8 * maximum(abs.(tr.uparB)))
    @test isapprox(sol.klittle_upar, tr.klittle_upar; rtol=r8, atol=r8 * maximum(abs.(tr.klittle_upar)))
    @test isapprox(sol.kbig_upar, tr.kbig_upar; rtol=r8, atol=r8 * maximum(abs.(tr.kbig_upar)))
    @test isapprox(sol.vpol_th0, tr.vpol_th0; rtol=r8, atol=r8 * maximum(abs.(tr.vpol_th0)))
    @test isapprox(sol.vtor_th0, tr.vtor_th0; rtol=r8, atol=r8 * maximum(abs.(tr.vtor_th0)))
    @test isapprox(sol.d_phi_sqavg, tr.d_phi_sqavg; rtol=r8, atol=r8 * abs(tr.d_phi_sqavg) + 1e-300)
    @test isapprox(sol.vtor_0order_th0, tr.vtor_0order_th0; rtol=r8, atol=1e-14)
    @test isapprox(sol.uparB_0order, tr.uparB_0order; rtol=r8, atol=1e-14)
    gvf = read_neo_numbers(joinpath(dir, "out.neo.transport_gv"))
    gv = reshape(gvf[2:end], 3, ns)
    @test isapprox(sol.pflux_gv, gv[1, :]; rtol=r8, atol=r8 * scale)
    @test isapprox(sol.eflux_gv, gv[2, :]; rtol=r8, atol=r8 * scale)
    @test isapprox(sol.mflux_gv, gv[3, :]; rtol=r8, atol=r8 * scale)
    prec = only(read_neo_numbers(joinpath(dir, "out.neo.prec")))
    @test isapprox(sol.check_sum, prec; rtol=tol.prec)
    # the tgyro block of out.neo.transport_flux (e14.5) and its lumped FluxSolution
    if p.ae_flag == 0
        v = read_neo_numbers(joinpath(dir, "out.neo.transport_flux"))
        tg = reshape(v[(8ns+1):end], 4, ns)   # Z, pflux, eflux, mflux per species
        pf, ef, mf = tgyro_fluxes(sol)
        fscale = maximum(abs.(tg[3, :]))
        @test isapprox(pf, tg[2, :]; rtol=1e-4, atol=1e-4 * fscale)
        @test isapprox(ef, tg[3, :]; rtol=1e-4, atol=1e-4 * fscale)
        @test isapprox(mf, tg[4, :]; rtol=1e-4, atol=1e-4 * fscale)
        fs = GACODE.FluxSolution(sol)
        @test fs.ENERGY_FLUX_e == ef[p.is_ele]
        @test fs.ENERGY_FLUX_i == sum(ef[is] for is in 1:ns if is != p.is_ele)
    end
    return nothing
end

# Tolerances. The port follows the Fortran operation by operation, but NEO's
# backward recursion for the negative-index fcoll entries (neo_compute_fcoll
# stage C) is unstable in double precision for 0.1 <= lambda <= 10, so the
# high-Legendre-order (ix >~ n_xi/2) field-particle blocks of such species
# pairs carry round-off amplified to O(1e-6) (lambda=1) or O(1) (reg12's C-on-D
# pair, lambda=0.17); Fortran and Julia both sit that far from a BigFloat
# evaluation of the same recursion. The low moments barely feel it: reg12's
# fluxes still agree to 3e-7, which is also how far two Fortran builds of reg12
# are apart (out.neo.prec 50.220128 here vs 50.220145 shipped with gacode).
const NEO_TOL_DEFAULT = (basis=1e-10, coll_test=5e-11, coll_field_lowix=1e-8, coll_field=1e-5, resid=1e-9, g=1e-8, transport=1e-8, prec=1e-7)
const NEO_TOL_SMALL = (basis=1e-10, coll_test=1e-11, coll_field_lowix=1e-8, coll_field=1e-8, resid=1e-9, g=1e-8, transport=1e-8, prec=1e-7)
const NEO_TOL_CASE = Dict(
    "small" => NEO_TOL_SMALL, "small_norot" => NEO_TOL_SMALL,
    "small_cm1" => NEO_TOL_SMALL, "small_cm2" => NEO_TOL_SMALL, "small_cm3" => NEO_TOL_SMALL, "small_cm5" => NEO_TOL_SMALL,
    "reg12" => (basis=1e-10, coll_test=1e-11, coll_field_lowix=1e-4, coll_field=2.0, resid=1e-2, g=1e-4, transport=1e-5, prec=1e-6))
# FUSE-built 5-species case (D, T, He, C, e at n_xi=17): several ion pairs have lambda in the
# unstable range, milder than reg12; the four (BTCCW, IPCCW) variants give identical magnitudes
const NEO_TOL_FUSE = (basis=1e-10, coll_test=1e-10, coll_field_lowix=1e-6, coll_field=1e-3, resid=1e-8, g=1e-6, transport=1e-6, prec=1e-7)
for c in ("fuse_miller", "fuse_bp1_ip1", "fuse_bm1_ip1", "fuse_bp1_im1", "fuse_bm1_im1")
    NEO_TOL_CASE[c] = NEO_TOL_FUSE
end

function test_neo_reference_case(case::String)
    tol = get(NEO_TOL_CASE, case, NEO_TOL_DEFAULT)
    dir = joinpath(NEO_REFDIR, case)
    p = NEOParams(joinpath(dir, "input.neo"))
    ns, ne, nxi, nth = p.n_species, p.n_energy, p.n_xi, p.n_theta
    if !isfile(joinpath(dir, "out.neo.dump"))
        # standard-output-only case: end-to-end comparison
        sol = solve_neo(p)
        test_neo_transport_outputs(sol, dir, get(NEO_TOL_CASE, case, NEO_TOL_DEFAULT))
        return sol
    end
    D = read_neo_dump(joinpath(dir, "out.neo.dump"))
    @test (ns, ne, nxi, nth) == Int.((D["n_species"], D["n_energy"], D["n_xi"], D["n_theta"]))

    @testset "profiles" begin
        @test p.z ≈ _vec(D["z"]) rtol = 1e-14
        @test p.mass ≈ _vec(D["mass"]) rtol = 1e-14
        @test p.dens ≈ _vec(D["dens"]) rtol = 1e-14
        @test p.temp ≈ _vec(D["temp"]) rtol = 1e-14
        @test p.nu ≈ _vec(D["nu"]) rtol = 1e-13
        @test p.vth ≈ _vec(D["vth"]) rtol = 1e-14
        @test p.q ≈ D["q"] && p.rho ≈ D["rho"] && p.rmin ≈ D["r"] && p.rmaj ≈ D["rmaj"]
        @test p.sign_q == D["sign_q"] && p.sign_bunit == D["sign_bunit"]
        @test p.dphi0dr == D["dphi0dr"] && p.omega_rot == D["omega_rot"] && p.omega_rot_deriv == D["omega_rot_deriv"]
        @test p.is_ele == Int(D["is_ele"]) && p.ae_flag == Int(D["ae_flag"])
    end

    basis = NEOBasis(p)
    @testset "energy basis" begin
        @test parent(basis.e_lag) == Int.(vec(D["e_lag"]))
        @test parent(basis.xi_beta_l) == Int.(vec(D["xi_beta_l"]))
        @test basis.mygamma2 ≈ vec(D["mygamma2"]) rtol = 1e-14
        for name in (:evec_e0, :evec_e1, :evec_e2, :evec_e05, :evec_e105, :emat_e05, :emat_en05, :emat_e05de, :emat_e0, :emat_e1)
            @test relnorm(parent(getfield(basis, name)), D[String(name)]) < tol.basis
        end
    end

    coll = collision_ints(p, basis; serial=true)
    @testset "collision matrices" begin
        # normwise per (is,js,ix) block: the alternating Laguerre sums cancel digits
        @test blockwise_relerr(parent(coll.test), D["emat_coll_test"]) < tol.coll_test
        lowix = 1:(nxi ÷ 2 + 1)
        @test blockwise_relerr(parent(coll.field)[:, :, :, :, lowix], D["emat_coll_field"][:, :, :, :, lowix]) < tol.coll_field_lowix
        @test blockwise_relerr(parent(coll.field), D["emat_coll_field"]) < tol.coll_field
        if Threads.nthreads() > 1
            collt = collision_ints(p, basis; serial=false)
            @test collt.test == coll.test && collt.field == coll.field
        end
        # density and momentum conservation of like-species collisions (full FP operators)
        for is in (p.collision_model in (4, 5) ? (1:ns) : (1:0))
            blk0 = coll.test[is, is, :, :, 0] + coll.field[is, is, :, :, 0]
            @test maximum(abs.(blk0[:, 0])) < 1e-12 * maximum(abs.(blk0))
            if nxi >= 1
                blk1 = coll.test[is, is, :, :, 1] + coll.field[is, is, :, :, 1]
                @test maximum(abs.(blk1[:, 0])) < 1e-12 * maximum(abs.(blk1))
            end
        end
        if isfile(joinpath(dir, "out.neo.diagnostic_colltest"))
            tm, fm = collision_ints_mono(p, basis)
            tref = read_neo_diagnostic_coll(joinpath(dir, "out.neo.diagnostic_colltest"), ns, ne, nxi)
            fref = read_neo_diagnostic_coll(joinpath(dir, "out.neo.diagnostic_collfield"), ns, ne, nxi)
            @test maximum(abs.(tm - tref)) < 5e-8 * maximum(abs.(tref))
            @test maximum(abs.(fm - fref)) < 5e-8 * maximum(abs.(fref))
        end
    end

    geo = equilibrium(p)
    @testset "equilibrium" begin
        @test geo.theta ≈ vec(D["theta"]) atol = 1e-15
        @test geo.d_theta ≈ D["d_theta"] rtol = 1e-15
        for name in (:k_par, :v_drift_x, :v_drift_th, :gradr, :gradpar_gradr, :w_theta, :Btor, :Bpol, :Bmag,
            :Bmag_rderiv, :gradpar_Bmag, :bigR, :bigR_rderiv, :gradpar_bigR)
            @test relnorm(getfield(geo, name), vec(D[String(name)])) < 1e-12
        end
        for name in (:bigR_th0, :bigR_th0_rderiv, :gradr_th0, :Btor_th0, :Bpol_th0, :Bmag_th0, :Bmag_th0_rderiv,
            :I_div_psip, :Bmag2_avg, :Bmag2inv_avg, :Btor2_avg, :bigRinv_avg, :gradpar_Bmag2_avg, :ftrap)
            @test isapprox(getfield(geo, name), D[String(name)]; rtol=1e-12, atol=1e-14)
        end
    end

    rot = rotation_phi(p, geo)
    @testset "rotation" begin
        @test isapprox(rot.phi_rot, vec(D["phi_rot"]); rtol=1e-12, atol=1e-14)
        @test isapprox(rot.phi_rot_deriv, vec(D["phi_rot_deriv"]); rtol=1e-11, atol=1e-13)
        @test isapprox(rot.phi_rot_rderiv, vec(D["phi_rot_rderiv"]); rtol=1e-11, atol=1e-13)
        @test isapprox(rot.phi_rot_avg, D["phi_rot_avg"]; rtol=1e-12, atol=1e-14)
        @test relnorm(rot.dens_fac, reshape(D["dens_fac"], ns, nth)) < 1e-12
    end

    pattern = NEOPattern(p)
    A, b, coef = assemble(p, basis, coll, geo, rot, pattern; serial=true)
    g_fortran = read_neo_numbers(joinpath(dir, "out.neo.f"))
    @testset "assembly" begin
        @test pattern.n_row == Int(D["n_row"]) == length(g_fortran)
        # g_fortran carries 16 digits; the residual is limited by cond(A)
        @test norm(A * g_fortran - b) < tol.resid * norm(b)
        if Threads.nthreads() > 1
            At, bt, _ = assemble(p, basis, coll, geo, rot, pattern; serial=false)
            @test At.nzval == A.nzval && bt == b
        end
    end

    g = solve_system(A, b, pattern)
    @testset "solve" begin
        @test relnorm(g, g_fortran) < tol.g
    end

    sol = transport(p, basis, geo, rot, coef, g, pattern)
    @testset "transport" begin
        dke = reshape(D["neo_dke_out"], ns, 6)   # pflux, eflux, mflux, eflux-omega*mflux, vpol_th0, vtor_th0+vtor_0order_th0
        gv = reshape(D["neo_gv_out"], ns, 4)
        d1 = vec(D["neo_dke_1d_out"])
        rt = tol.transport
        scale = maximum(abs.(dke[:, 2]))   # the energy flux sets the flux scale (pflux ~ 0 by ambipolarity)
        @test isapprox(sol.pflux, dke[:, 1]; rtol=rt, atol=rt * scale)
        @test isapprox(sol.eflux, dke[:, 2]; rtol=rt, atol=rt * scale)
        @test isapprox(sol.mflux, dke[:, 3]; rtol=rt, atol=rt * scale)
        @test isapprox(sol.vpol_th0, dke[:, 5]; rtol=rt, atol=rt * maximum(abs.(dke[:, 5])))
        @test isapprox(sol.vtor_th0 .+ sol.vtor_0order_th0, dke[:, 6]; rtol=rt, atol=rt * maximum(abs.(dke[:, 6])))
        # gv fluxes can be round-off (s-alpha: <gradpar_Bmag/B^4> = 0); compare on the dke flux scale
        @test isapprox(sol.pflux_gv, gv[:, 1]; rtol=rt, atol=rt * scale)
        @test isapprox(sol.eflux_gv, gv[:, 2]; rtol=rt, atol=rt * scale)
        @test isapprox(sol.mflux_gv, gv[:, 3]; rtol=rt, atol=rt * scale)
        @test isapprox(sol.jpar, d1[1]; rtol=rt, atol=rt * scale)
        @test isapprox(sol.jtor, d1[2]; rtol=rt, atol=rt * scale)

        test_neo_transport_outputs(sol, dir, tol)
    end

    @testset "driver" begin
        sol2 = solve_neo(p; serial=true)
        @test sol2.check_sum == sol.check_sum
        if Threads.nthreads() > 1
            sol3 = solve_neo(p; serial=false)
            @test sol3.check_sum == sol.check_sum && sol3.pflux == sol.pflux
        end
        if p.ae_flag == 0
            fs = GACODE.FluxSolution(sol)
            pf, ef, mf = tgyro_fluxes(sol)
            @test fs.ENERGY_FLUX_e == ef[p.is_ele]
            @test length(fs.PARTICLE_FLUX_i) == ns - 1
        end
    end
    return sol
end

# ---- unit tests -----------------------------------------------------------

@testset "NEONative" begin
    @testset "gamma2" begin
        for n in 1:40
            @test gamma2(n) ≈ SF.gamma(n / 2) rtol = 1e-14
        end
    end

    @testset "fcoll" begin
        fdump = joinpath(NEO_REFDIR, "small", "out.neo.dump_fcoll")
        if isfile(fdump)
            F = read_neo_dump(fdump)
            m0 = Int(F["m0"])
            lams = vec(F["lambda"])
            # read_neo_dump keeps the last block of a repeated name; re-read sequentially
            lines = readlines(fdump)
            blocks = Vector{Matrix{Float64}}()
            i = 1
            while i <= length(lines)
                hdr = split(lines[i])
                dims = parse.(Int, hdr[3:end])
                n = prod(dims)
                if hdr[2] in ("fcoll", "fcoll_bar")
                    push!(blocks, reshape([parse(Float64, lines[i+k]) for k in 1:n], dims...))
                end
                i += n + 1
            end
            @test length(blocks) == 2 * length(lams)
            for (il, lam) in enumerate(lams)
                f, fb = compute_fcoll(m0, lam)
                fref = blocks[2il-1]
                fbref = blocks[2il]
                # NEO's backward recursions for the negative-index entries are unstable in
                # double precision for 0.1 <= lambda <= 10 (Fortran and Julia both sit ~1e-5
                # from a BigFloat evaluation of the same algorithm); the well-conditioned
                # m,n >= -1 entries agree to round-off
                well = [m >= -1 && n >= -1 for m in -m0:m0, n in -m0:m0]
                @test isapprox(parent(f)[well], fref[well]; rtol=1e-12, atol=1e-14 * maximum(abs.(fref)))
                @test isapprox(parent(fb)[well], fbref[well]; rtol=1e-12, atol=1e-14 * maximum(abs.(fbref)))
                @test isapprox(parent(f), fref; rtol=1e-4, atol=1e-14 * maximum(abs.(fref)))
                @test isapprox(parent(fb), fbref; rtol=1e-4, atol=1e-14 * maximum(abs.(fbref)))
                fB, fbB = compute_fcoll(m0, big(lam))
                @test isapprox(parent(f), Float64.(parent(fB)); rtol=1e-4, atol=1e-14 * maximum(abs.(fref)))
                @test isapprox(parent(fb), Float64.(parent(fbB)); rtol=1e-4, atol=1e-14 * maximum(abs.(fbref)))
            end
        end
        # identities G + GB = B((m+1)/2,(n+1)/2) and G(m,n,λ) = GB(n,m,1/λ) for m,n >= 0;
        # NEO only forms (and the operator only reads) entries with m+n odd
        for lam in (0.03, 0.5, 1.0, 7.0, 25.0, 4000.0)
            m0 = 20
            f, fb = compute_fcoll(m0, lam)
            f2, fb2 = compute_fcoll(m0, 1 / lam)
            for m in 0:m0, n in 0:m0
                isodd(m + n) || continue
                gm = f[m, n] / (0.25 * SF.gamma((m + n + 2) / 2) / lam^((n + 1) / 2))
                gbm = fb[m, n] / (0.25 * SF.gamma((m + n + 2) / 2) / lam^((n + 1) / 2))
                beta = SF.gamma((m + 1) / 2) * SF.gamma((n + 1) / 2) / SF.gamma((m + n + 2) / 2)
                @test isapprox(gm + gbm, beta; rtol=1e-12)
                gb2 = fb2[n, m] / (0.25 * SF.gamma((m + n + 2) / 2) / (1 / lam)^((m + 1) / 2))
                @test isapprox(gm, gb2; rtol=1e-11, atol=1e-13)
            end
        end
    end

    @testset "batch driver" begin
        # several surfaces per task, so the per-task factorization cache is reused
        # across different matrices; threaded must equal serial exactly
        ps = [NEOParams(joinpath(NEO_REFDIR, c, "input.neo")) for c in ("small", "small_norot", "small_cm3", "small_cm5")]
        batch = repeat(ps, 3)
        ser = run_neo_native(batch; serial=true)
        thr = run_neo_native(batch)
        @test [s.check_sum for s in thr] == [s.check_sum for s in ser]
        @test all(thr[i].jpar == solve_neo(batch[i]).jpar for i in eachindex(batch))
        @test length(run_neo_native(ps[1:1])) == 1
    end

    @testset "ForwardDiff through solve_neo" begin
        # NEOParams{Dual}: everything up to the sparse solve propagates the Duals,
        # solve_system applies the implicit-function rule (values on A0, partials
        # from A0 g' = b' - A' g0 with the same factorization)
        d = NEONative.read_input_neo(joinpath(NEO_REFDIR, "small", "input.neo"))
        keys_ = ["DLNTDR_1", "DLNNDR_2", "TEMP_1", "OMEGA_ROT", "DENS_3"]
        x0 = [Float64(d[k]) for k in keys_]
        function neo_params(x)
            dx = copy(d)
            for (k, v) in zip(keys_, x)
                dx[k] = v
            end
            return NEOParams{eltype(x)}(dx)
        end
        function neo_outputs(x)
            sol = solve_neo(neo_params(x))
            return vcat(sol.pflux, sol.eflux, sol.mflux, sol.jpar, sol.uparB)
        end

        # 1. the rule against a dense generic LU of the Dual matrix (exact Dual propagation)
        xd = [ForwardDiff.Dual{:neo_test}(x0[k], ForwardDiff.Partials(ntuple(j -> Float64(j == k), 5))) for k in 1:5]
        D = eltype(xd)
        pd = neo_params(xd)
        basis = NEOBasis(pd)
        geo = NEONative.equilibrium(pd)
        rot = NEONative.rotation_phi(pd, geo)
        coll = collision_ints(pd, basis)
        pattern = NEONative.NEOPattern(pd)
        A, b, _ = NEONative.assemble(pd, basis, coll, geo, rot, pattern)
        @test eltype(A) === D
        g_rule = NEONative.solve_system(A, b, pattern)
        g_dense = Matrix(A) \ b
        @test ForwardDiff.value.(g_rule) ≈ ForwardDiff.value.(g_dense) rtol = 1e-10
        for k in 1:5
            @test ForwardDiff.partials.(g_rule, k) ≈ ForwardDiff.partials.(g_dense, k) rtol = 1e-10
        end

        # 2. end to end: ForwardDiff Jacobian of the transport outputs vs central differences
        y0 = neo_outputs(x0)
        J = ForwardDiff.jacobian(neo_outputs, x0)
        @test size(J) == (length(y0), 5)
        for k in 1:5
            h = 1e-6 * max(abs(x0[k]), 1.0)
            xp = copy(x0)
            xp[k] += h
            xm = copy(x0)
            xm[k] -= h
            Jfd = (neo_outputs(xp) - neo_outputs(xm)) / (2h)
            # the gradient columns are small (they only enter the RHS), so the
            # finite-difference round-off eps*|y|/h sets an absolute floor
            @test isapprox(J[:, k], Jfd; rtol=1e-4, atol=1e-7 * norm(y0))
        end

        # 3. Dual results, the Float64-only fcoll_exact path, the threaded batch
        sold = solve_neo(pd; keep_g=true)
        @test sold isa NEOSolution{D}
        @test GACODE.FluxSolution(sold) isa GACODE.FluxSolution{D}
        @test_throws ErrorException solve_neo(pd; fcoll_exact=true)
        batch = run_neo_native([pd, pd, pd])
        @test all(s.check_sum === sold.check_sum for s in batch)
        @test all(s.jpar === sold.jpar for s in batch)

        # 4. factorization reuse: a Dual solve at the primal point just solved does
        # not refactorize, and gives exactly the Float64 solution as its values
        cache = NEOFactorCache()
        s0 = solve_neo(NEOParams{Float64}(d); cache, keep_g=true)
        F0 = cache.F
        sd = solve_neo(pd; cache, keep_g=true)
        @test cache.F === F0
        @test ForwardDiff.value.(sd.g) == s0.g
        # changed values on the same pattern refactorize in place (symbolic reuse)
        s1 = solve_neo(neo_params(x0 .* 1.01); cache)
        @test cache.F === F0
        @test s1.check_sum != s0.check_sum

        # 5. caller-owned caches in the batch driver persist across calls
        ps = [NEOParams(joinpath(NEO_REFDIR, c, "input.neo")) for c in ("small", "small_norot")]
        caches = [NEOFactorCache() for _ in ps]
        r1 = run_neo_native(ps; caches)
        Fs = [c.F for c in caches]
        r2 = run_neo_native(ps; caches)
        @test all(caches[i].F === Fs[i] for i in eachindex(ps))
        @test [s.check_sum for s in r2] == [s.check_sum for s in r1]
        @test_throws ErrorException run_neo_native(ps; caches=caches[1:1])
    end

    @testset "stale-factorization GMRES refinement" begin
        d = NEONative.read_input_neo(joinpath(NEO_REFDIR, "small", "input.neo"))
        function params_scaled(f)
            dx = copy(d)
            for k in ("TEMP_1", "TEMP_2", "TEMP_3", "DENS_1", "DENS_2", "DLNTDR_1", "OMEGA_ROT")
                dx[k] = d[k] * f
            end
            return NEOParams{Float64}(dx)
        end
        cache = NEOFactorCache(; refine=true)
        s0 = solve_neo(params_scaled(1.0); cache, keep_g=true)
        @test cache.n_factorizations == 1
        @test cache.resid0 < 1e-8
        # changed matrices: GMRES on the old factorization, no refactorization
        for f in (1.01, 1.1)
            p1 = params_scaled(f)
            s1 = solve_neo(p1; cache, keep_g=true)
            @test cache.n_factorizations == 1
            s1_fresh = solve_neo(p1; keep_g=true)
            @test isapprox(s1.g, s1_fresh.g; rtol=1e-7)
            @test isapprox(s1.check_sum, s1_fresh.check_sum; rtol=1e-8)
            @test all(isapprox.(s1.eflux, s1_fresh.eflux; rtol=1e-7))
        end
        @test cache.n_refined == 2
        # a very different matrix is still solved correctly (refined or refactorized)
        p2 = params_scaled(4.0)
        s2 = solve_neo(p2; cache, keep_g=true)
        @test isapprox(s2.g, solve_neo(p2; keep_g=true).g; rtol=1e-7)
        # Dual solve at a nearby point: the values go through GMRES, the partial
        # block refactorizes, and the next Dual pass at the same point is a cache hit
        nfact = cache.n_factorizations
        xd = ForwardDiff.Dual{:neo_refine}(1.0, ForwardDiff.Partials((1.0, 0.0)))
        dx = copy(d)
        dx["TEMP_1"] = d["TEMP_1"] * 4.0 * (1 + 1e-3 * xd)
        pd = NEOParams{eltype(xd)}(dx)
        sd = solve_neo(pd; cache)
        @test cache.n_factorizations == nfact + 1
        sd2 = solve_neo(pd; cache)
        @test cache.n_factorizations == nfact + 1
        @test sd2.check_sum === sd.check_sum
        dx0 = copy(d)
        dx0["TEMP_1"] = d["TEMP_1"] * 4.0 * (1 + 1e-3)
        @test ForwardDiff.value(sd.check_sum) ≈ solve_neo(NEOParams{Float64}(dx0)).check_sum rtol = 1e-9
        # the default cache never refines
        plain = NEOFactorCache()
        solve_neo(params_scaled(1.0); cache=plain)
        solve_neo(params_scaled(1.01); cache=plain)
        @test plain.n_factorizations == 2 && plain.n_refined == 0
    end

    cases = ["small", "small_norot", "small_cm1", "small_cm2", "small_cm3", "small_cm5"]
    if NEO_FULL
        append!(cases, ["reg01", "reg02", "reg03", "reg04", "reg05", "reg06", "reg07", "reg08", "reg09", "reg10", "reg11",
            "reg12", "reg13", "reg14", "reg15",
            "fuse_miller", "fuse_bp1_ip1", "fuse_bm1_ip1", "fuse_bp1_im1", "fuse_bm1_im1"])
    end
    for case in cases
        isdir(joinpath(NEO_REFDIR, case)) || continue
        @testset "reference $case" begin
            test_neo_reference_case(case)
        end
    end
end

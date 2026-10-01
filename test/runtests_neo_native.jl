# Reference tests of the native NEO port (src/neo/) against Fortran NEO.
#
# The reference data in test/neo_reference/<case>/ is produced by
# test/neo_reference/generate.sh (utilities/serial_neo/neo_dump.f90) and holds
# every intermediate array at full precision (out.neo.dump), the solution
# vector (out.neo.f) and NEO's standard outputs. The `small*` cases always
# run; the gacode regression cases run with NEO_NATIVE_FULL=1.

using NeoclassicalTransport
using NeoclassicalTransport: NEOParams, NEONative
using NeoclassicalTransport.NEONative: NEOBasis, gamma2, compute_fcoll, collision_ints, collision_ints_mono,
    equilibrium, rotation_phi, NEOPattern, assemble, solve_system, transport, solve_neo, tgyro_fluxes, E_ALPHA
using LinearAlgebra
using SparseArrays
using OffsetArrays
using Test
import GACODE

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

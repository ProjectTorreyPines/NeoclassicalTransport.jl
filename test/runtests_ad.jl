using NeoclassicalTransport
using NeoclassicalTransport.IMAS
using ForwardDiff
using Test

dd = IMAS.json2imas(joinpath(@__DIR__, "highbetap_fpp_325_ods.json"))
eqt = dd.equilibrium.time_slice[]
cp1d = dd.core_profiles.profiles_1d[]

ir = 5  # use an interior radial index

pp = NeoclassicalTransport.get_plasma_profiles(eqt, cp1d)
eg = NeoclassicalTransport.get_equilibrium_geometry(eqt, cp1d)

# Baseline evaluation
sol0 = NeoclassicalTransport.hirshmansigmar(ir, eqt, cp1d, pp, eg)

"""
    convert_plasma_profiles(::Type{T}, pp::NeoclassicalTransport.PlasmaProfiles{S}) where {T<:Real, S}

Create a PlasmaProfiles{T} from an existing PlasmaProfiles by converting all fields to type T.
"""
function convert_plasma_profiles(::Type{T}, pp::NeoclassicalTransport.PlasmaProfiles{S}) where {T<:Real,S}
    return NeoclassicalTransport.PlasmaProfiles{T}(;
        Z=T.(pp.Z),
        mass=T.(pp.mass),
        dens=T.(pp.dens),
        temp=T.(pp.temp),
        nu=T.(pp.nu),
        dlnndr=T.(pp.dlnndr),
        dlntdr=T.(pp.dlntdr),
        vth=T.(pp.vth),
    )
end

"""
    convert_equilibrium_geometry(::Type{T}, eg::NeoclassicalTransport.EquilibriumGeometry{S}) where {T<:Real, S}

Create an EquilibriumGeometry{T} from an existing EquilibriumGeometry by converting all fields to type T.
"""
function convert_equilibrium_geometry(::Type{T}, eg::NeoclassicalTransport.EquilibriumGeometry{S}) where {T<:Real,S}
    return NeoclassicalTransport.EquilibriumGeometry{T}(;
        time=eg.time,
        rmin=T.(eg.rmin),
        rmaj=T.(eg.rmaj),
        a=T(eg.a),
        q=T.(eg.q),
        ftrap=T.(eg.ftrap),
        Bmag2_avg=T.(eg.Bmag2_avg),
        f=T.(eg.f),
    )
end

@testset "Hirshman-Sigmar AD" begin
    h = 1e-5  # finite-difference step
    fd_rtol = 1e-3  # FD truncation error tolerance

    # ---- Test derivative w.r.t. dlntdr (ion 1) ----
    @testset "d(ENERGY_FLUX_i)/d(dlntdr_ion1)" begin
        x0 = pp.dlntdr[ir, 1]

        function hs_Qi_dlntdr(x::T) where {T<:Real}
            pp_T = convert_plasma_profiles(T, pp)
            eg_T = convert_equilibrium_geometry(T, eg)
            pp_T.dlntdr[ir, 1] = x
            pflux, eflux = NeoclassicalTransport.compute_HS(ir, eqt, cp1d, pp_T, eg_T)
            return sum(eflux[1:end-1])
        end

        dQi_ad = ForwardDiff.derivative(hs_Qi_dlntdr, x0)
        dQi_fd = (hs_Qi_dlntdr(x0 + h) - hs_Qi_dlntdr(x0 - h)) / (2h)

        @test isfinite(dQi_ad)
        @test isapprox(dQi_ad, dQi_fd; rtol=fd_rtol)
        println("d(Qi)/d(dlntdr_i1): AD=$dQi_ad  FD=$dQi_fd")
    end

    # ---- Test derivative w.r.t. dlnndr (electron) ----
    @testset "d(PARTICLE_FLUX_e)/d(dlnndr_e)" begin
        x0 = pp.dlnndr[ir, end]

        function hs_Ge_dlnndr(x::T) where {T<:Real}
            pp_T = convert_plasma_profiles(T, pp)
            eg_T = convert_equilibrium_geometry(T, eg)
            pp_T.dlnndr[ir, end] = x
            pflux, eflux = NeoclassicalTransport.compute_HS(ir, eqt, cp1d, pp_T, eg_T)
            return pflux[end]
        end

        dGe_ad = ForwardDiff.derivative(hs_Ge_dlnndr, x0)
        dGe_fd = (hs_Ge_dlnndr(x0 + h) - hs_Ge_dlnndr(x0 - h)) / (2h)

        @test isfinite(dGe_ad)
        @test isapprox(dGe_ad, dGe_fd; rtol=fd_rtol)
        println("d(Γe)/d(dlnndr_e): AD=$dGe_ad  FD=$dGe_fd")
    end

    # ---- Test derivative w.r.t. collision frequency (nu, ion 1) ----
    @testset "d(ENERGY_FLUX_e)/d(nu_ion1)" begin
        x0 = pp.nu[ir, 1]

        function hs_Qe_nu(x::T) where {T<:Real}
            pp_T = convert_plasma_profiles(T, pp)
            eg_T = convert_equilibrium_geometry(T, eg)
            pp_T.nu[ir, 1] = x
            pflux, eflux = NeoclassicalTransport.compute_HS(ir, eqt, cp1d, pp_T, eg_T)
            return eflux[end]
        end

        dQe_ad = ForwardDiff.derivative(hs_Qe_nu, x0)
        dQe_fd = (hs_Qe_nu(x0 + h) - hs_Qe_nu(x0 - h)) / (2h)

        @test isfinite(dQe_ad)
        @test isapprox(dQe_ad, dQe_fd; rtol=fd_rtol)
        println("d(Qe)/d(ν_i1):    AD=$dQe_ad  FD=$dQe_fd")
    end

    # ---- Test gradient w.r.t. multiple inputs simultaneously ----
    @testset "gradient w.r.t. [dlntdr_i1, dlnndr_i1, dlntdr_e]" begin
        x0 = [pp.dlntdr[ir, 1], pp.dlnndr[ir, 1], pp.dlntdr[ir, end]]

        function hs_Qi_vec(x::AbstractVector{T}) where {T<:Real}
            pp_T = convert_plasma_profiles(T, pp)
            eg_T = convert_equilibrium_geometry(T, eg)
            pp_T.dlntdr[ir, 1] = x[1]
            pp_T.dlnndr[ir, 1] = x[2]
            pp_T.dlntdr[ir, end] = x[3]
            pflux, eflux = NeoclassicalTransport.compute_HS(ir, eqt, cp1d, pp_T, eg_T)
            return sum(eflux[1:end-1])
        end

        grad_ad = ForwardDiff.gradient(hs_Qi_vec, x0)
        grad_fd = similar(x0)
        for k in eachindex(x0)
            xp = copy(x0); xp[k] += h
            xm = copy(x0); xm[k] -= h
            grad_fd[k] = (hs_Qi_vec(xp) - hs_Qi_vec(xm)) / (2h)
        end

        @test all(isfinite, grad_ad)
        @test isapprox(grad_ad, grad_fd; rtol=fd_rtol)
        println("∇Qi (AD): $grad_ad")
        println("∇Qi (FD): $grad_fd")
    end
end

# ---- NEO-NN: ForwardDiff through the ensemble forward pass ----
# The NN layer is generic in eltype (Float64 weights x Dual inputs), so
# gradients of the surrogate fluxes w.r.t. any feature come for free.
#
# Differentiating from a *plasma quantity* rather than a raw feature (i.e.
# through InputNEO -> InputNEONN -> the nets) is covered by the "neo_nn.jl AD"
# testset in runtests.jl; what follows is the heavier full-Jacobian check that
# does not belong in CI.
@testset "NEO-NN AD" begin
    ens = NeoclassicalTransport.loadmodelonce("neonn_d3d+mastu+nstx_flux")
    ineo = NeoclassicalTransport.InputNEO(eqt, cp1d, ir)
    nn = NeoclassicalTransport.InputNEONN(ineo)
    xnames_val = NeoclassicalTransport._get_xnames_without_log10_suffix(ens)
    x0 = NeoclassicalTransport._extract_fields!(zeros(length(ens.xnames)), nn, xnames_val)

    iy = findfirst(==("OUT_eflux_tgyro_ion1"), ens.ynames)
    nn_Qi(x) = NeoclassicalTransport.flux_array(ens, x; warn_nn_train_bounds=false)[iy]

    grad_ad = ForwardDiff.gradient(nn_Qi, x0)
    @test all(isfinite, grad_ad)
    @test any(!iszero, grad_ad)

    # central-difference cross-check on the bulk-ion temperature gradient feature
    k = findfirst(==("DLNTDR_1"), ens.xnames)
    hx = 1e-4 * max(abs(x0[k]), 1.0)
    xp = copy(x0); xp[k] += hx
    xm = copy(x0); xm[k] -= hx
    fd = (nn_Qi(xp) - nn_Qi(xm)) / (2hx)
    @test isapprox(grad_ad[k], fd; rtol=1e-4)
    println("dQi/dDLNTDR_1 (AD): $(grad_ad[k])  (FD): $fd")
end

# ---- NEO-NN: full Jacobian of every output channel w.r.t. every feature ----
@testset "NEO-NN Jacobian" begin
    ens = NeoclassicalTransport.loadmodelonce("neonn_d3d+mastu+nstx_flux")
    ineo = NeoclassicalTransport.InputNEO(eqt, cp1d, ir)
    nn = NeoclassicalTransport.InputNEONN(ineo)
    xnames_val = NeoclassicalTransport._get_xnames_without_log10_suffix(ens)
    x0 = NeoclassicalTransport._extract_fields!(zeros(length(ens.xnames)), nn, xnames_val)

    f(x) = NeoclassicalTransport.flux_array(ens, x; warn_nn_train_bounds=false)
    J = ForwardDiff.jacobian(f, x0)
    @test size(J) == (length(ens.ynames), length(ens.xnames))
    @test all(isfinite, J)

    # every feature must move at least one channel: a column of exact zeros
    # would mean the net ignores an input it was trained on
    dead = [ens.xnames[k] for k in axes(J, 2) if all(iszero, J[:, k])]
    @test isempty(dead)

    # central-difference the whole Jacobian. The step has to scale with the
    # feature: several inputs are O(1e-3) (NU_1, RHO_STAR) and are fed to the
    # net through log10, so an absolute step of 1e-5 is a ~1% excursion and the
    # truncation error alone reaches 0.5%.
    Jfd = similar(J)
    for k in eachindex(x0)
        h = 1e-6 * max(abs(x0[k]), 1e-4)
        xp = copy(x0); xp[k] += h
        xm = copy(x0); xm[k] -= h
        Jfd[:, k] = (f(xp) - f(xm)) / (2h)
    end
    @test isapprox(J, Jfd; rtol=1e-4, atol=1e-6 * maximum(abs, J))

    worst = argmax(abs.(J .- Jfd) ./ max.(abs.(Jfd), 1e-12))
    println("largest AD-FD Jacobian deviation at ",
            ens.ynames[worst[1]], " / ", ens.xnames[worst[2]],
            ": AD=", J[worst], " FD=", Jfd[worst])
end

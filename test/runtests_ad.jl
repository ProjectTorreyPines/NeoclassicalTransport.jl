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


@testset "Yancc AD" begin
    NT = NeoclassicalTransport
    rho_cp = cp1d.grid.rho_tor_norm
    i = argmin(abs.(rho_cp .- 0.75))
    geometry = NT.get_yancc_geometry(eqt, cp1d, [i]; ntheta=21)
    field = geometry.fields[1]
    rho_s = NT.GACODE.rho_s(cp1d, eqt)
    rmin = NT.GACODE.r_min_core_profiles(eqt.profiles_1d, rho_cp)

    # profiles of the element type of the arguments: the ion temperatures tilted by xT
    # (value and gradient change), all densities scaled by xn
    function profiles(xT::T, xn::T) where {T<:Real}
        cp = IMAS.core_profiles__profiles_1d{T}()
        cp.grid.rho_tor_norm = T.(rho_cp)
        cp.electrons.temperature = T.(cp1d.electrons.temperature)
        cp.electrons.density_thermal = xn .* cp1d.electrons.density_thermal
        cp.electrons.density = xn .* cp1d.electrons.density
        resize!(cp.ion, length(cp1d.ion))
        for (k, ion) in enumerate(cp1d.ion)
            resize!(cp.ion[k].element, 1)
            cp.ion[k].element[1].a = T(ion.element[1].a)
            cp.ion[k].element[1].z_n = T(ion.element[1].z_n)
            cp.ion[k].temperature = ion.temperature .* (1 .+ (xT - 1) .* rho_cp)
            cp.ion[k].density_thermal = xn .* ion.density_thermal
            cp.ion[k].density = xn .* ion.density
        end
        return cp
    end
    fluxes(sol) = [sol.ENERGY_FLUX_e, sol.ENERGY_FLUX_i, sol.PARTICLE_FLUX_e, sol.PARTICLE_FLUX_i[1]]
    solve(xT, xn; kw...) = fluxes(NT.yancc(i, eqt, profiles(xT, xn), field; geometry.a, rho_s, rmin, nx=4, nalpha=25, kw...))

    # both derivatives from one call with two partials
    cache = NT.YanccCache()
    xT = ForwardDiff.Dual{Nothing}(1.0, 1.0, 0.0)
    xn = ForwardDiff.Dual{Nothing}(1.0, 0.0, 1.0)
    ad = solve(xT, xn; cache, rtol=1e-10, tangent_rtol=1e-9)
    reference = solve(1.0, 1.0; rtol=1e-11)
    @test ForwardDiff.value.(ad) ≈ reference rtol = 1e-7
    h = 1e-4
    for (k, (hT, hn)) in enumerate(((h, 0.0), (0.0, h)))
        fd = (solve(1 + hT, 1 + hn; rtol=1e-11) .- solve(1 - hT, 1 - hn; rtol=1e-11)) ./ 2h
        @test ForwardDiff.partials.(ad, k) ≈ fd rtol = 1e-4
    end
    println("d(Qe, Qi, Γe, Γi1)/d(xT): AD=$(ForwardDiff.partials.(ad, 1))")
    # the second call at the same profiles reuses the derivative solves
    solves = cache.tangent_solves
    again = solve(xT, xn; cache, rtol=1e-10, tangent_rtol=1e-9)
    @test ForwardDiff.partials.(again, 1) ≈ ForwardDiff.partials.(ad, 1) rtol = 1e-10
    @test cache.tangent_solves == solves && solves <= 2 && cache.rebuilds == 1
end

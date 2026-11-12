using NeoclassicalTransport
using NeoclassicalTransport.IMAS
using ForwardDiff
using Test

dd = IMAS.json2imas(joinpath(@__DIR__, "highbetap_fpp_325_ods.json"))
eqt = dd.equilibrium.time_slice[];
cp1d = dd.core_profiles.profiles_1d[];

@testset "NeoclassicalTransport" begin
    @testset "chang_hinton.jl" begin
        rho_fluxmatch = 0.5
        iion = 1
        sol = NeoclassicalTransport.changhinton(eqt, cp1d, rho_fluxmatch, iion)
        @test isapprox(sol.ENERGY_FLUX_i, 0.05211, rtol=0.10)
    end

    @testset "hirshman_sigmar.jl" begin
        ir = 1

        parameter_matrices = NeoclassicalTransport.get_plasma_profiles(eqt, cp1d)
        equilibrium_geometry = NeoclassicalTransport.get_equilibrium_geometry(eqt, cp1d)

        sol = NeoclassicalTransport.hirshmansigmar(ir, eqt, cp1d, parameter_matrices, equilibrium_geometry)
        @test isapprox(sol.PARTICLE_FLUX_e, 0.0005395, rtol=0.10)
        @test isapprox(sol.PARTICLE_FLUX_i[1], 0.00032, rtol=0.10)
        @test isapprox(sol.PARTICLE_FLUX_i[2], 0.000398, rtol=0.10)
        @test isapprox(sol.PARTICLE_FLUX_i[3], -0.00015, rtol=0.10)
        @test isapprox(sol.PARTICLE_FLUX_i[4], -1.232e-6, rtol=0.10)
        @test isapprox(sol.ENERGY_FLUX_e, 0.002116, rtol=0.10)
        @test isapprox(sol.ENERGY_FLUX_i, 0.0163, rtol=0.10)
    end

    @testset "input_neo.jl" begin
        gridpoint_cp = 1
        input_neo = NeoclassicalTransport.InputNEO(eqt, cp1d, gridpoint_cp)
        @test true
    end
    @testset "yancc.jl" begin
        NT = NeoclassicalTransport
        Yancc = NT.Yancc
        rho_cp = cp1d.grid.rho_tor_norm

        # Miller surface: circular limit and the safety factor carried by the field
        f = NT.miller_field(; r=0.5, a=1.0, R0=3.0, kappa=1.0, delta=0.0, q=2.0, Bunit=5.0, ntheta=32)
        R = 3.0 .+ 0.5 .* cos.(f.theta)
        F = 5.0 * 3.0 * sqrt(1 - (0.5 / 3.0)^2)   # q = r F / (R0² Bp(R0)) √(1 - ε²) for circles
        @test vec(f.B_sub_z) ≈ fill(F, 32) rtol = 1e-12
        @test vec(f.Bmag) ≈ hypot.(0.5 * 5.0 / 2.0 ./ R, F ./ R) rtol = 1e-12
        @test Yancc.flux_surface_average(f, f.B_sup_t) / Yancc.flux_surface_average(f, f.B_sup_z) ≈ -1 / 2.0 rtol = 1e-12
        shaped = (; r=0.6, a=1.1, R0=3.1, kappa=1.7, delta=0.35, zeta=-0.05, shift=-0.2, s_kappa=0.3, s_delta=0.4, s_zeta=-0.1, q=3.3, Bunit=4.0)
        f = NT.miller_field(; shaped..., ntheta=64)
        @test Yancc.flux_surface_average(f, f.B_sup_t) / Yancc.flux_surface_average(f, f.B_sup_z) ≈ -1 / 3.3 rtol = 1e-10
        # ∇·B = 0: √g B^θ is constant on the surface
        @test maximum(abs, f.sqrtg .* f.B_sup_t .- f.sqrtg[1] * f.B_sup_t[1]) < 1e-12 * abs(f.sqrtg[1] * f.B_sup_t[1])

        # the test case at rho = 0.75, posed to Yancc as InputNEO poses it to NEO
        i = argmin(abs.(rho_cp .- 0.75))
        geometry = NT.get_yancc_geometry(eqt, cp1d, [i]; ntheta=41)
        field = geometry.fields[1]
        input_neo = NT.InputNEO(eqt, cp1d, i)
        @test field.rho ≈ input_neo.RMIN_OVER_A
        @test abs(field.iota) ≈ 1 / input_neo.Q
        species, lnlambda = NT.yancc_species(input_neo, cp1d.electrons.density[i], cp1d.electrons.temperature[i], geometry.a)
        @test lnlambda ≈ IMAS.lnΛ_ei(cp1d.electrons.density[i], cp1d.electrons.temperature[i]) rtol = 1e-3
        solve(field; nx=6, nalpha=65, rtol=1e-6) = Yancc.solve_dke_tokamak(Yancc.DKEOperator(field, Yancc.UniformPitchAngleGrid(nalpha),
            Yancc.MaxwellSpeedGrid(nx), species, 0.0; p1="4z", lnlambda); rtol)
        raw = solve(field)
        sol = NT.yancc_to_GB(raw, 1:4, eqt, cp1d, i, geometry.a)

        # Fortran NEO without rotation at N_ENERGY=8, N_XI=25, N_THETA=25 (it is still
        # rising with resolution: 5.50e-3, 0.1769, 1.54e-3 at the default 6, 17, 17)
        @test isapprox(sol.ENERGY_FLUX_e, 5.7635e-3, rtol=0.05)
        @test isapprox(sol.ENERGY_FLUX_i, 0.17978, rtol=0.03)
        @test isapprox(sol.PARTICLE_FLUX_e, 1.6738e-3, rtol=0.08)
        @test isapprox(sol.PARTICLE_FLUX_i[4], 6.596e-5, rtol=0.08)
        # ambipolar up to the discretization error of the ion fluxes, which nearly cancel
        Z = [s.species.charge / Yancc.ELEMENTARY_CHARGE for s in species]
        @test abs(sum(Z[1:4] .* sol.PARTICLE_FLUX_i) - sol.PARTICLE_FLUX_e) < 0.15 * sum(abs, Z[1:4] .* sol.PARTICLE_FLUX_i)

        # fluxes do not depend on the directions of the field and of the current
        coarse = solve(NT.get_yancc_geometry(eqt, cp1d, [i]; ntheta=21).fields[1]; nx=4, nalpha=25, rtol=1e-9)
        Bunit = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, NT.GACODE.bunit(eqt.profiles_1d))(rho_cp[i])
        for (bt, ip) in ((1, 1), (1, -1), (-1, 1), (-1, -1))
            input_neo.BTCCW, input_neo.IPCCW = bt, ip
            flipped = solve(NT.miller_field(input_neo, geometry.a, Bunit; ntheta=21); nx=4, nalpha=25, rtol=1e-9)
            @test flipped.heat_flux ≈ coarse.heat_flux rtol = 1e-5
            @test flipped.particle_flux[end] ≈ coarse.particle_flux[end] rtol = 1e-5
            @test abs(flipped.JparB) ≈ abs(coarse.JparB) rtol = 1e-5
        end

        # from the profiles, with a cache; the conversion to IMAS returns the SI fluxes
        rho_s = NT.GACODE.rho_s(cp1d, eqt)
        cache = NT.YanccCache()
        sol = NT.yancc(i, eqt, cp1d, field; geometry.a, cache)
        @test sol.ENERGY_FLUX_e > 0 && sol.ENERGY_FLUX_i > 0
        @test length(sol.PARTICLE_FLUX_i) == length(cp1d.ion)
        iterations = cache.iterations
        again = NT.yancc(i, eqt, cp1d, field; geometry.a, cache)
        @test again.ENERGY_FLUX_i == sol.ENERGY_FLUX_i
        @test cache.iterations == iterations && cache.rebuilds == 1 && cache.solves == 2

        # the species from the profiles are those of InputNEO: total densities (thermal + fast),
        # every ion at its own temperature
        dd_fast = deepcopy(dd)
        cp_fast = dd_fast.core_profiles.profiles_1d[]
        cp_fast.ion[1].density_fast = 0.1 .* cp_fast.ion[1].density_thermal
        empty!(cp_fast.ion[1], :density)   # stored in the file; now thermal + fast
        cp_fast.ion[2].temperature = 1.3 .* cp_fast.ion[2].temperature
        rmin_over_a = NT.GACODE.r_min_core_profiles(eqt.profiles_1d, rho_cp) ./ (geometry.a * IMAS.cgs.m_to_cm)
        species, kept = NT.yancc_species(cp_fast, rmin_over_a, i)
        as_neo, _ = NT.yancc_species(NT.InputNEO(eqt, cp_fast, i), cp_fast.electrons.density[i], cp_fast.electrons.temperature[i], geometry.a)
        @test kept == 1:length(cp_fast.ion)
        @test species[1].density ≈ 1.1 * cp_fast.ion[1].density_thermal[i] rtol = 1e-12
        @test species[2].temperature == cp_fast.ion[2].temperature[i]
        for (sp, ref) in zip(species, as_neo)
            @test sp.temperature ≈ ref.temperature rtol = 1e-10
            @test sp.density ≈ ref.density rtol = 1e-10
            @test sp.dTdrho ≈ ref.dTdrho rtol = 1e-8
            @test sp.dndrho ≈ ref.dndrho rtol = 1e-8
            @test sp.species.mass ≈ ref.species.mass rtol = 1e-2   # atomic mass unit here, proton mass in InputNEO
        end

        model = resize!(dd.core_transport.model, :neoclassical; wipe=false)
        m1d = resize!(model.profiles_1d)
        m1d.grid_flux.rho_tor_norm = [rho_cp[i]]
        NT.GACODE.flux_gacode_to_imas((:electron_energy_flux, :ion_energy_flux, :electron_particle_flux), [sol], m1d, eqt, cp1d)
        ieq = argmin(abs.(eqt.profiles_1d.rho_tor_norm .- rho_cp[i]))
        vprime = NT.GACODE.volume_prime_miller_correction(eqt)[ieq]
        raw = cache.solution
        @test m1d.electrons.energy.flux[1] ≈ geometry.a * raw.heat_flux[end] * vprime rtol = 1e-10
        @test m1d.total_ion_energy.flux[1] ≈ geometry.a * sum(raw.heat_flux[1:end-1]) * vprime rtol = 1e-10
        @test m1d.electrons.particles.flux[1] ≈ geometry.a * raw.particle_flux[end] * vprime rtol = 1e-10
    end
end

include("runtests_ad.jl")

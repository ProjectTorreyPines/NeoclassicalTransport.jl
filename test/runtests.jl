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
        # a solve that stops above `max_residual` is an error, not fluxes
        @test_throws ErrorException NT.yancc(i, eqt, cp1d, field; geometry.a, max_residual=0.0)

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

        # the same case in the conventions of NEO's output: with the electrons normalizing the
        # input, the gyro-Bohm fluxes are those of `yancc_to_GB`
        input_neo = NT.InputNEO(eqt, cp1d, i)
        ne, Te = cp1d.electrons.density[i], cp1d.electrons.temperature[i]
        op = NT.yancc_operator(input_neo, ne, Te, geometry.a; rotation=false, nx=6)
        raw = Yancc.solve_dke_tokamak(op; rtol=1e-6)
        out = NT.yancc_neo_outputs(op, raw, input_neo, ne, Te, geometry.a)
        ref = NT.yancc_to_GB(raw, 1:4, eqt, cp1d, i, geometry.a)
        @test out.eflux_tgyro[end] ≈ ref.ENERGY_FLUX_e rtol = 2e-3   # a and ρ_s of InputNEO and of GACODE
        @test sum(out.eflux_tgyro[1:end-1]) ≈ ref.ENERGY_FLUX_i rtol = 2e-3
        @test out.pflux_tgyro[end] ≈ ref.PARTICLE_FLUX_e rtol = 2e-3
        @test out.eflux_dke == out.eflux_tgyro && all(iszero, out.pflux_gv) && all(iszero, out.mflux_dke)
        # they do not depend on the values that normalize the input
        op2 = NT.yancc_operator(input_neo, 3e19, 700.0, 0.6; rotation=false, nx=6)
        out2 = NT.yancc_neo_outputs(op2, Yancc.solve_dke_tokamak(op2; rtol=1e-6), input_neo, 3e19, 700.0, 0.6)
        for name in (:pflux_tgyro, :eflux_tgyro, :uparB, :vpol)
            @test getfield(out2, name) ≈ getfield(out, name) rtol = 1e-4
        end
        @test out2.jpar ≈ out.jpar rtol = 1e-4
    end

    # needs a Yancc.jl with toroidal rotation
    NeoclassicalTransport.yancc_has_rotation() && @testset "yancc.jl with sonic rotation" begin
        NT = NeoclassicalTransport
        Yancc = NT.Yancc
        # a spherical tokamak surface at r/a = 0.3, Mach number 0.5 at the outboard midplane,
        # normalized to the main ion as in a NEO input file
        input_neo = NT.InputNEO(; N_SPECIES=3, IPCCW=1, BTCCW=-1, ROTATION_MODEL=2,
            RMIN_OVER_A=0.3, RMAJ_OVER_A=1.7536662411465478, RHO_STAR=7.5578865445764644e-3, SHIFT=-0.15070487035015717,
            S_ZMAG=0.0, Q=2.4209627337754567, SHEAR=0.11331244035462903, KAPPA=2.5811717959967213, S_KAPPA=1.6289195249156868e-2,
            DELTA=0.12702113731164183, S_DELTA=9.5310575457983993e-2, ZETA=-1.2833695904231810e-2, S_ZETA=0.0,
            OMEGA_ROT=-0.23941167743090741, OMEGA_ROT_DERIV=8.8488439698855101e-2, NU_1=5.3152002561413416e-3,
            Z_1=1, MASS_1=1.00026045, DENS_1=1.0, TEMP_1=1.0, DLNNDR_1=2.2483299820280460, DLNTDR_1=0.43181336523337882,
            Z_2=6, MASS_2=5.9650047842589409, DENS_2=0.10288811961885454, TEMP_2=1.0010177587237705,
            DLNNDR_2=-0.48314066249864296, DLNTDR_2=0.46765248541359755,
            Z_3=-1, MASS_3=2.72308510742333e-4, DENS_3=1.6173287177131272, TEMP_3=0.83007795899279391,
            DLNNDR_3=1.2057373145816577, DLNTDR_3=1.6013811000907354)
        n_norm, T_norm, a = 1e19, 1e3, 1.0
        rotation = NT.yancc_rotation(input_neo, T_norm, a)
        @test rotation.omega ≈ -0.23941167743090741 * sqrt(1e3 * Yancc.ELEMENTARY_CHARGE / (IMAS.cgs.md / 1e3)) rtol = 1e-12
        @test rotation.domega_drho / rotation.omega ≈ 8.8488439698855101e-2 / -0.23941167743090741 rtol = 1e-12
        op = NT.yancc_operator(input_neo, n_norm, T_norm, a; nx=8)
        @test op.rotation !== nothing
        sol = Yancc.solve_dke_tokamak(op; rtol=1e-8)
        out = NT.yancc_neo_outputs(op, sol, input_neo, n_norm, T_norm, a)

        # Fortran NEO, ROTATION_MODEL=2, COLLISION_MODEL=4 at N_ENERGY=12, N_XI=35, N_THETA=35
        # (within 0.5% of its fluxes at 12, 19, 19)
        @test out.pflux_tgyro ≈ [0.63741, -0.097443, 0.052754] rtol = 0.02
        @test out.eflux_tgyro ≈ [1.6981, -0.17850, 0.14579] rtol = 0.02
        @test out.eflux_dke ≈ [1.8944, -0.35923, 0.14579] rtol = 0.02
        @test out.mflux_tgyro[1:2] ≈ [-0.74701, 0.68775] rtol = 0.02
        @test out.mflux_tgyro[3] ≈ -0.14060e-4 rtol = 0.05
        @test out.uparB ≈ [-0.079794848, -0.066067769, -0.014306856] rtol = 0.02
        @test out.vpol ≈ [-0.014156700, 0.027618043, 0.049470210] rtol = 0.03
        @test out.jpar ≈ -0.095451423 rtol = 0.02
        # the fluxes are ambipolar and nothing is gyroviscous on this up-down symmetric surface
        Z = [1, 6, -1]
        @test abs(sum(Z .* out.pflux_tgyro)) < 1e-2 * sum(abs, Z .* out.pflux_tgyro)
        @test maximum(abs, out.pflux_gv) < 1e-8 && maximum(abs, out.mflux_gv) < 1e-8

        # without the rotation the momentum flux vanishes and the energy fluxes of the two frames agree
        op0 = NT.yancc_operator(input_neo, n_norm, T_norm, a; rotation=false, nx=8)
        out0 = NT.yancc_neo_outputs(op0, Yancc.solve_dke_tokamak(op0; rtol=1e-8), input_neo, n_norm, T_norm, a)
        @test out0.eflux_dke == out0.eflux_tgyro
        @test !isapprox(out0.eflux_tgyro[1], out.eflux_tgyro[1]; rtol=0.05)

    end
end

include("runtests_ad.jl")

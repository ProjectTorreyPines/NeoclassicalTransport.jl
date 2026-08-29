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

    @testset "neo_nn.jl" begin
        # --- model loading: all 24 shipped ensembles
        shipped = [(dev, grp) for dev in ("d3d", "d3dedge", "d3dnearedge", "d3dnegdedge", "mastu+nstx", "d3d+mastu+nstx",
                                          "mastu+nstx_withnegD", "mastunearedge+nstxnearedge_withnegD", "mastuedge+nstxedge_withnegD",
                                          "d3d_withnegD", "d3dnearedge_withnegD", "d3dedge_withnegD")
                              for grp in ("flux", "flow")]
        for (dev, grp) in shipped
            name = "neonn_$(dev)_$(grp)"
            @test name in NeoclassicalTransport.available_models()
            ens = NeoclassicalTransport.loadmodelonce(name)
            @test length(ens.models) == 20
            @test length(ens.xnames) == (startswith(dev, "d3d") && dev != "d3d+mastu+nstx" ? 24 : 25)  # d3d-only nets lack TEMP_2
            @test length(ens.ynames) == (grp == "flux" ? 9 : 4)
            # every feature must be constructible from InputNEONN
            for x in ens.xnames
                @test Symbol(replace(x, "_log10" => "")) in fieldnames(NeoclassicalTransport.InputNEONN)
            end
        end

        # --- feature construction at mid-radius gridpoints
        rho = cp1d.grid.rho_tor_norm
        gps = [argmin(abs.(rho .- r)) for r in (0.3, 0.5, 0.7)]
        input_neos = [NeoclassicalTransport.InputNEO(eqt, cp1d, gp) for gp in gps]
        for (gp, ineo) in zip(gps, input_neos)
            nn = NeoclassicalTransport.InputNEONN(ineo)
            for f in fieldnames(NeoclassicalTransport.InputNEONN)
                @test isfinite(getfield(nn, f))
            end
            # bulk-ion normalization: TEMP_3 = Te/T_1
            @test isapprox(nn.TEMP_3, cp1d.electrons.temperature[gp] / cp1d.ion[1].temperature[gp]; rtol=1e-10)
            # training helicity convention: signed q, negative
            @test nn.Q < 0
        end

        # --- golden values, joint (default) models on the test ODS.
        # Regenerate after any model or conversion change with:
        #   julia --project -e 'using NeoclassicalTransport, IMAS; dd = IMAS.json2imas("test/highbetap_fpp_325_ods.json");
        #     eqt = dd.equilibrium.time_slice[]; cp1d = dd.core_profiles.profiles_1d[]; rho = cp1d.grid.rho_tor_norm;
        #     gps = [argmin(abs.(rho .- r)) for r in (0.3, 0.5, 0.7)];
        #     ineos = [NeoclassicalTransport.InputNEO(eqt, cp1d, gp) for gp in gps];
        #     display(NeoclassicalTransport.run_neonn(ineos; warn_nn_train_bounds=false));
        #     display(NeoclassicalTransport.run_neonn_flow(ineos; warn_nn_train_bounds=false))'
        sols = NeoclassicalTransport.run_neonn(input_neos; warn_nn_train_bounds=true)  # mid-radius must be in-distribution: warnings on
        golden_flux = [
            (Qe=0.0017178346833441264, Qi=0.025556752585612613, Ge=0.00012735449643806496, Gi=[-0.0014022468789535898, 0.0002596746436104749], Pi=7.451180100801163e-5),
            (Qe=0.0038533148862862797, Qi=0.06928415085842432, Ge=0.0006699483897839101, Gi=[-5.847681304467084e-5, 0.00012696371720343812], Pi=2.6813228672561452e-5),
            (Qe=0.006351501902653892, Qi=0.13854652844309343, Ge=0.0014658874150405535, Gi=[0.0015188631077603938, -5.002779301155509e-6], Pi=0.00019910699840269505),
        ]
        for (sol, g) in zip(sols, golden_flux)
            @test isapprox(sol.ENERGY_FLUX_e, g.Qe; rtol=1e-5)
            @test isapprox(sol.ENERGY_FLUX_i, g.Qi; rtol=1e-5)
            @test isapprox(sol.PARTICLE_FLUX_e, g.Ge; rtol=1e-5)
            @test length(sol.PARTICLE_FLUX_i) == 2  # [bulk ion, lumped impurity]
            @test isapprox(sol.PARTICLE_FLUX_i[1], g.Gi[1]; rtol=1e-5)
            @test isapprox(sol.PARTICLE_FLUX_i[2], g.Gi[2]; rtol=1e-5)
            @test isapprox(sol.STRESS_TOR_i, g.Pi; rtol=1e-5)
        end

        flows = NeoclassicalTransport.run_neonn_flow(input_neos; warn_nn_train_bounds=false)
        golden_flow = [
            (v1=-0.003047472114153771, v2=-0.00028947016414673646, ve=0.0035875258110280835, j=-0.009058944858801718),
            (v1=-0.0033977939772391608, v2=1.9522040705542046e-5, ve=0.004621595736978964, j=-0.011749084373043311),
            (v1=-0.003392818187325397, v2=8.847442323627779e-5, ve=0.004537927313130893, j=-0.01462742379790355),
        ]
        for (flow, g) in zip(flows, golden_flow)
            @test isapprox(flow.vpol_ion1, g.v1; rtol=1e-5)
            @test isapprox(flow.vpol_ion2, g.v2; rtol=1e-5)
            @test isapprox(flow.vpol_elec, g.ve; rtol=1e-5)
            @test isapprox(flow.jpar, g.j; rtol=1e-5)
        end

        # --- ensemble uncertainty
        su = NeoclassicalTransport.run_neonn(input_neos[2]; uncertain=true, warn_nn_train_bounds=false)
        @test su.ENERGY_FLUX_i isa NeoclassicalTransport.Measurements.Measurement
        @test NeoclassicalTransport.Measurements.uncertainty(su.ENERGY_FLUX_i) > 0
        @test isapprox(NeoclassicalTransport.Measurements.value(su.ENERGY_FLUX_i), golden_flux[2].Qi; rtol=1e-5)

        # --- joint (25-feature) vs d3d (24-feature) nets: same input, consistent physics.
        # Smoke check on feature extraction/conventions for both signatures, not accuracy.
        sols_d3d = NeoclassicalTransport.run_neonn(input_neos; model_filename="neonn_d3d_flux", warn_nn_train_bounds=false)
        for (s_joint, s_d3d) in zip(sols, sols_d3d)
            @test sign(s_joint.ENERGY_FLUX_i) == sign(s_d3d.ENERGY_FLUX_i)
            @test isapprox(s_joint.ENERGY_FLUX_i, s_d3d.ENERGY_FLUX_i; rtol=1.0)
        end

        # --- flux/flow model mix-up guard
        @test_throws ErrorException NeoclassicalTransport.run_neonn(input_neos[1]; model_filename="neonn_d3d+mastu+nstx_flow")
        @test_throws ErrorException NeoclassicalTransport.run_neonn_flow(input_neos[1]; model_filename="neonn_d3d+mastu+nstx_flux")

        # --- neoclassical Er closure (force balance on the NEO-NN poloidal flow)
        gps_er = [argmin(abs.(rho .- r)) for r in (0.5, 0.9, 0.95)]
        input_neos_er = [NeoclassicalTransport.InputNEO(eqt, cp1d, gp) for gp in gps_er]
        Rout = NeoclassicalTransport.IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.r_outboard).(rho)
        vtor = fill(5.0e4, length(gps_er))  # m/s, stand-in for a measured toroidal rotation
        ers = NeoclassicalTransport.neoclassical_Er(eqt, cp1d, gps_er, vtor; species=:impurity, warn_nn_train_bounds=false)
        @test length(ers) == length(gps_er)
        for (gp, s) in zip(gps_er, ers)
            # the three terms are the whole story
            @test isapprox(s.Er, s.Er_pressure + s.Er_vtor + s.Er_vpol; rtol=1e-12)
            @test all(isfinite, (s.Er, s.Er_pressure, s.Er_vtor, s.Er_vpol, s.vpol, s.grad_r))
            # geometry at theta=0: |grad r| > 1 on this shaped equilibrium, B_tor ~ F/R
            @test s.grad_r > 1.0
            @test isapprox(s.B_tor, NeoclassicalTransport.IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.f).(rho)[gp] / s.R_omp; rtol=1e-10)
            @test isapprox(s.R_omp, Rout[gp]; rtol=1e-10)
            @test abs(s.B_pol) < abs(s.B_tor)
            # decreasing profiles + positive charge -> diamagnetic term digs the well
            @test s.Er_pressure < 0
            # the toroidal term is exactly -v_tor*B_pol
            @test isapprox(s.Er_vtor, -5.0e4 * s.B_pol; rtol=1e-12)
        end

        # Er is linear in the supplied toroidal velocity, slope -B_pol
        ers2 = NeoclassicalTransport.neoclassical_Er(eqt, cp1d, gps_er, 2 .* vtor; species=:impurity, warn_nn_train_bounds=false)
        for (s1, s2) in zip(ers, ers2)
            @test isapprox(s2.Er - s1.Er, -5.0e4 * s1.B_pol; rtol=1e-10)
        end

        # scalar method, and the species -> network-output mapping shares one v_norm
        s_one = NeoclassicalTransport.neoclassical_Er(eqt, cp1d, gps_er[2], 5.0e4; species=:impurity, warn_nn_train_bounds=false)
        @test isapprox(s_one.Er, ers[2].Er; rtol=1e-12)
        s_bulk = NeoclassicalTransport.neoclassical_Er(eqt, cp1d, gps_er[2], 5.0e4; species=:bulk, warn_nn_train_bounds=false)
        s_elec = NeoclassicalTransport.neoclassical_Er(eqt, cp1d, gps_er[2], 5.0e4; species=:electron, warn_nn_train_bounds=false)
        flow2 = NeoclassicalTransport.run_neonn_flow(input_neos_er[2]; warn_nn_train_bounds=false)
        @test isapprox(s_one.vpol / flow2.vpol_ion2, s_bulk.vpol / flow2.vpol_ion1; rtol=1e-10)
        @test isapprox(s_bulk.vpol / flow2.vpol_ion1, s_elec.vpol / flow2.vpol_elec; rtol=1e-10)
        # electrons: Z < 0 flips the diamagnetic term
        @test s_elec.Er_pressure > 0

        # ensemble uncertainty propagates through the closure
        s_unc = NeoclassicalTransport.neoclassical_Er(eqt, cp1d, gps_er[2], 5.0e4; uncertain=true, warn_nn_train_bounds=false)
        @test s_unc.Er isa NeoclassicalTransport.Measurements.Measurement
        @test NeoclassicalTransport.Measurements.uncertainty(s_unc.Er) > 0
        @test isapprox(NeoclassicalTransport.Measurements.value(s_unc.Er), ers[2].Er; rtol=1e-5)

        # input validation
        @test_throws ErrorException NeoclassicalTransport.neoclassical_Er(eqt, cp1d, gps_er, vtor; species=:carbon)
        @test_throws ErrorException NeoclassicalTransport.neoclassical_Er(eqt, cp1d, gps_er, vtor[1:2])

        # --- radial-family blending: a family's core net switches to its
        # near-edge / edge nets at the family's switch r/a (0.881 / 0.975, the
        # TurbulentTransport withnegD boundaries). A blended point must match
        # the region net evaluated directly; a core point must be unaffected.
        cand = [NeoclassicalTransport.InputNEO(eqt, cp1d, gp) for gp in findall(rho .>= 0.80)]
        for (core, variants) in NeoclassicalTransport._RADIAL_BLEND_VARIANTS
            nearedge, edge, r_ne, r_ed = variants
            i_ne = findfirst(ineo -> r_ne <= ineo.RMIN_OVER_A < r_ed, cand)
            i_ed = findfirst(ineo -> ineo.RMIN_OVER_A >= r_ed, cand)
            @test i_ne !== nothing  # test ODS must reach into the near-edge region
            c_flux = NeoclassicalTransport.run_neonn(input_neos[2]; model_filename="$(core)_flux", warn_nn_train_bounds=false)
            for (idx, variant) in ((i_ne, nearedge), (i_ed, edge))
                idx === nothing && continue
                ineo = cand[idx]
                b_flux = NeoclassicalTransport.run_neonn([input_neos[2], ineo]; model_filename="$(core)_flux", warn_nn_train_bounds=false)
                d_flux = NeoclassicalTransport.run_neonn(ineo; model_filename="$(variant)_flux", warn_nn_train_bounds=false)
                @test isapprox(b_flux[2].ENERGY_FLUX_e, d_flux.ENERGY_FLUX_e; rtol=1e-12)
                @test isapprox(b_flux[2].ENERGY_FLUX_i, d_flux.ENERGY_FLUX_i; rtol=1e-12)
                @test isapprox(b_flux[1].ENERGY_FLUX_i, c_flux.ENERGY_FLUX_i; rtol=1e-12)  # core point untouched
                b_flow = NeoclassicalTransport.run_neonn_flow([ineo]; model_filename="$(core)_flow", warn_nn_train_bounds=false)
                d_flow = NeoclassicalTransport.run_neonn_flow(ineo; model_filename="$(variant)_flow", warn_nn_train_bounds=false)
                @test isapprox(only(b_flow).jpar, d_flow.jpar; rtol=1e-12)
                @test isapprox(only(b_flow).vpol_ion1, d_flow.vpol_ion1; rtol=1e-12)
            end
        end
    end

    @testset "neo_nn.jl AD" begin
        # The surrogates are differentiable end-to-end: InputNEO carries the
        # element type, so a Dual set on a plasma quantity propagates through the
        # species lumping, the electron -> bulk-ion normalization, the log10
        # feature transform and the ensemble forward pass.
        rho = cp1d.grid.rho_tor_norm
        gps = [argmin(abs.(rho .- r)) for r in (0.3, 0.5, 0.7)]
        ineo0 = NeoclassicalTransport.InputNEO(eqt, cp1d, gps[2])
        @test ineo0 isa NeoclassicalTransport.InputNEO{Float64}
        @test eltype(ineo0) === Float64
        @test NeoclassicalTransport.InputNEO{Float32}(ineo0) isa NeoclassicalTransport.InputNEO{Float32}
        @test NeoclassicalTransport.InputNEO{Float32}(ineo0).DELTA ≈ ineo0.DELTA rtol = 1e-6
        @test NeoclassicalTransport.InputNEO{Float32}(ineo0).N_SPECIES === ineo0.N_SPECIES  # ints stay ints

        # d(flux)/d(plasma quantity), against a central difference
        function Qi(x::AbstractVector{D}) where {D<:Real}
            ineo = NeoclassicalTransport.InputNEO{D}(ineo0)
            ineo.DLNTDR_1, ineo.DLNTDR_3, ineo.TEMP_1 = x[1], x[2], x[3]
            return NeoclassicalTransport.run_neonn(ineo; warn_nn_train_bounds=false).ENERGY_FLUX_i
        end
        x0 = [ineo0.DLNTDR_1, ineo0.DLNTDR_3, ineo0.TEMP_1]
        grad = ForwardDiff.gradient(Qi, x0)
        @test all(isfinite, grad)
        @test any(!iszero, grad)
        for k in eachindex(x0)
            h = 1e-6 * max(abs(x0[k]), 1.0)
            xp = copy(x0); xp[k] += h
            xm = copy(x0); xm[k] -= h
            @test isapprox(grad[k], (Qi(xp) - Qi(xm)) / (2h); rtol=1e-5)
        end

        # ... and through the flow nets
        function jpar(x::AbstractVector{D}) where {D<:Real}
            ineo = NeoclassicalTransport.InputNEO{D}(ineo0)
            ineo.NU_1 = x[1]
            return NeoclassicalTransport.run_neonn_flow(ineo; warn_nn_train_bounds=false).jpar
        end
        y0 = [ineo0.NU_1]
        gj = only(ForwardDiff.gradient(jpar, y0))
        h = 1e-6 * max(abs(y0[1]), 1.0)
        @test isfinite(gj) && !iszero(gj)
        @test isapprox(gj, (jpar([y0[1] + h]) - jpar([y0[1] - h])) / (2h); rtol=1e-4)

        # structurally zero: the bulk-ion density gradient is never read — the
        # feature is rebuilt from the electron and impurity gradients by
        # quasineutrality (see InputNEONN)
        dQi_dLn1 = ForwardDiff.derivative(x -> begin
                ineo = NeoclassicalTransport.InputNEO{typeof(x)}(ineo0)
                ineo.DLNNDR_1 = x
                NeoclassicalTransport.run_neonn(ineo; warn_nn_train_bounds=false).ENERGY_FLUX_i
            end, ineo0.DLNNDR_1)
        @test dQi_dLn1 == 0

        # AD survives the radial-family blending path (masked region nets)
        i_ne = findfirst(gp -> begin
                ineo = NeoclassicalTransport.InputNEO(eqt, cp1d, gp)
                NeoclassicalTransport._RADIAL_BLEND_NEAREDGE_RMIN <= ineo.RMIN_OVER_A < NeoclassicalTransport._RADIAL_BLEND_EDGE_RMIN
            end, findall(rho .>= 0.80))
        if i_ne !== nothing
            ineo_ne = NeoclassicalTransport.InputNEO(eqt, cp1d, findall(rho .>= 0.80)[i_ne])
            g_ne = ForwardDiff.derivative(x -> begin
                    ineo = NeoclassicalTransport.InputNEO{typeof(x)}(ineo_ne)
                    ineo.DLNTDR_1 = x
                    NeoclassicalTransport.run_neonn(ineo; model_filename="neonn_d3d_flux", warn_nn_train_bounds=false).ENERGY_FLUX_i
                end, ineo_ne.DLNTDR_1)
            @test isfinite(g_ne) && !iszero(g_ne)
        end

        # --- a dd carrying Duals (FUSE's AD path): InputNEO inherits the element
        # type and the derivative of a flux w.r.t. a profile value must match a
        # central difference over the Float64 dd (calc_z couples neighbouring
        # grid points; the FD sees the same coupling)
        D = ForwardDiff.Dual{Nothing,Float64,1}
        dd_ad = IMAS.dd{D}()
        IMAS.fill!(dd_ad, dd)
        eqt_ad = dd_ad.equilibrium.time_slice[]
        cp1d_ad = dd_ad.core_profiles.profiles_1d[]
        gp = gps[2]
        Ti = cp1d_ad.ion[1].temperature
        cp1d_ad.ion[1].temperature = [D(ForwardDiff.value(t), ForwardDiff.Partials((i == gp ? 1.0 : 0.0,))) for (i, t) in enumerate(Ti)]
        ineo_ad = NeoclassicalTransport.InputNEO(eqt_ad, cp1d_ad, gp)
        @test ineo_ad isa NeoclassicalTransport.InputNEO{D}
        @test ineo_ad.N_SPECIES === ineo0.N_SPECIES
        @test ineo_ad.BTCCW === ineo0.BTCCW
        Qi_ad = NeoclassicalTransport.run_neonn(ineo_ad; warn_nn_train_bounds=false).ENERGY_FLUX_i
        @test Qi_ad isa D
        @test ForwardDiff.value(Qi_ad) ≈ NeoclassicalTransport.run_neonn(ineo0; warn_nn_train_bounds=false).ENERGY_FLUX_i rtol = 1e-10
        function Qi_of_Ti(δ)
            dd_p = deepcopy(dd)
            cp = dd_p.core_profiles.profiles_1d[]
            cp.ion[1].temperature[gp] += δ
            ineo = NeoclassicalTransport.InputNEO(dd_p.equilibrium.time_slice[], cp, gp)
            return NeoclassicalTransport.run_neonn(ineo; warn_nn_train_bounds=false).ENERGY_FLUX_i
        end
        h = 1e-4 * ForwardDiff.value(Ti[gp])
        dQi_fd = (Qi_of_Ti(h) - Qi_of_Ti(-h)) / (2h)
        @test isapprox(ForwardDiff.partials(Qi_ad, 1), dQi_fd; rtol=1e-4)
        # Float64 equilibrium + Dual profiles (frozen geometry) is the setting
        # neoclassical_Er supports
        ineo_mixed = NeoclassicalTransport.InputNEO(eqt, cp1d_ad, gp)
        @test ineo_mixed isa NeoclassicalTransport.InputNEO{D}
        er_ad = NeoclassicalTransport.neoclassical_Er(eqt, cp1d_ad, gp, 5.0e4; warn_nn_train_bounds=false)
        @test er_ad.Er isa D
        @test isfinite(ForwardDiff.partials(er_ad.Er, 1)) && !iszero(ForwardDiff.partials(er_ad.Er, 1))

        # the two non-differentiable entry points refuse Duals explicitly
        @test_throws ErrorException NeoclassicalTransport.run_neonn(ineo_ad; uncertain=true, warn_nn_train_bounds=false)
        @test_throws ErrorException NeoclassicalTransport.save_inputneo(ineo_ad, tempname())
    end
end

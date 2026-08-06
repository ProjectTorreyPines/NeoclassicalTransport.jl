using NeoclassicalTransport
using NeoclassicalTransport.IMAS
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
        # --- model loading: all 6 shipped ensembles
        for dev in ("d3d", "mastu+nstx", "d3d+mastu+nstx"), grp in ("flux", "flow")
            name = "neonn_tgyro_$(dev)_$(grp)_v1"
            @test name in NeoclassicalTransport.available_models()
            ens = NeoclassicalTransport.loadmodelonce(name)
            @test length(ens.models) == 20
            @test length(ens.xnames) == (dev == "d3d" ? 24 : 25)  # d3d nets lack TEMP_2
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
        sols_d3d = NeoclassicalTransport.run_neonn(input_neos; model_filename="neonn_tgyro_d3d_flux_v1", warn_nn_train_bounds=false)
        for (s_joint, s_d3d) in zip(sols, sols_d3d)
            @test sign(s_joint.ENERGY_FLUX_i) == sign(s_d3d.ENERGY_FLUX_i)
            @test isapprox(s_joint.ENERGY_FLUX_i, s_d3d.ENERGY_FLUX_i; rtol=1.0)
        end

        # --- flux/flow model mix-up guard
        @test_throws ErrorException NeoclassicalTransport.run_neonn(input_neos[1]; model_filename="neonn_tgyro_d3d+mastu+nstx_flow_v1")
        @test_throws ErrorException NeoclassicalTransport.run_neonn_flow(input_neos[1]; model_filename="neonn_tgyro_d3d+mastu+nstx_flux_v1")
    end
end
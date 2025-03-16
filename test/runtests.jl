using NEO 
using Test 
using IMAS

dd = IMAS.json2imas(joinpath(@__DIR__,"highbetap_fpp_325_ods.json"))
eqt = dd.equilibrium.time_slice[] 
cp1d = dd.core_profiles.profiles_1d[]

@testset "NEO" begin
    @testset "chang_hinton.jl" begin
        rho_fluxmatch = 0.5
        iion = 1
        works = try
            sol = NEO.changhinton(eqt, cp1d, rho_fluxmatch, iion)
            @test isapprox(sol.ENERGY_FLUX_i, 0.05211, rtol=0.10) 
            true
        catch
            false
        end
        @test works==true
    end 

    @testset "hirshman_sigmar.jl" begin 
        ir = 1

        parameter_matrices = NEO.get_ion_electron_parameters(eqt, cp1d)
        equilibrium_geometry = NEO.get_equilibrium_parameters(eqt, cp1d)

        works = try
            sol = NEO.hirshmansigmar(ir, eqt, cp1d, parameter_matrices, equilibrium_geometry)
            @test isapprox(sol.PARTICLE_FLUX_e, 0.0005395, rtol=0.10)
            @test isapprox(sol.PARTICLE_FLUX_i[1], 0.00032, rtol=0.10)
            @test isapprox(sol.PARTICLE_FLUX_i[2], 0.000398, rtol=0.10)
            @test isapprox(sol.PARTICLE_FLUX_i[3], -0.00015, rtol=0.10)
            @test isapprox(sol.PARTICLE_FLUX_i[4], -1.232e-6, rtol=0.10)
            @test isapprox(sol.ENERGY_FLUX_e, 0.002116, rtol=0.10)
            @test isapprox(sol.ENERGY_FLUX_i, 0.0163, rtol=0.10)
            true
        catch
            false
        end
        @test works==true
    end

    @testset "input_neo.jl" begin
        gridpoint_cp = 1 
        works = try
            input_neo = NEO.InputNEO(eqt, cp1d, gridpoint_cp)
            true 
        catch
            false
        end 
        @test works==true
    end
end
using NEO 
using Test 
using IMAS

dd = IMAS.json2imas(joinpath(@__DIR__,"highbetap_fpp_325_ods.json"))
eqt = dd.equilibrium.time_slice[] 
cp1d = dd.core_profiles.profiles_1d[]

@testset "NEO" begin
    @testset "chang_hinton.jl" begin
        rho_fluxmatch = 1
        iion = 1
        works = try
            NEO.changhinton(eqt, cp1d, rho_fluxmatch, iion)
            true
        catch
            false
        end
        @test works==true
    end 

    @testset "hirshman_sigmar.jl" begin 
        ir = 1

        parameter_matrices = NEO.get_ion_electron_parameters(dd)
        equilibrium_geometry = NEO.get_equilibrium_parameters(dd)

        works = try
            NEO.hirshmansigmar(ir, dd, parameter_matrices, equilibrium_geometry)
            true
        catch
            false
        end
        @test works==true
    end

    @testset "input_neo.jl" begin
        gridpoint_cp = 1 
        works = try
            input_neo = NEO.InputNEO(dd, gridpoint_cp)
            true 
        catch
            false
        end 
        @test works==true
    end
end
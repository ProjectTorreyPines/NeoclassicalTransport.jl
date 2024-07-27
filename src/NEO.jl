module NEO

Base.@kwdef mutable struct flux_solution
    PARTICLE_FLUX_1::Union{Float64,Missing} = missing
    ENERGY_FLUX_1::Union{Float64,Missing} = missing
    MOMENTUM_FLUX_1::Union{Float64,Missing} = missing

    PARTICLE_FLUX_2::Union{Float64,Missing} = missing
    ENERGY_FLUX_2::Union{Float64,Missing} = missing
    MOMENTUM_FLUX_2::Union{Float64,Missing} = missing

    PARTICLE_FLUX_3::Union{Float64,Missing} = missing
    ENERGY_FLUX_3::Union{Float64,Missing} = missing
    MOMENTUM_FLUX_3::Union{Float64,Missing} = missing

    PARTICLE_FLUX_4::Union{Float64,Missing} = missing
    ENERGY_FLUX_4::Union{Float64,Missing} = missing
    MOMENTUM_FLUX_4::Union{Float64,Missing} = missing

    PARTICLE_FLUX_5::Union{Float64,Missing} = missing
    ENERGY_FLUX_5::Union{Float64,Missing} = missing
    MOMENTUM_FLUX_5::Union{Float64,Missing} = missing

    PARTICLE_FLUX_6::Union{Float64,Missing} = missing
    ENERGY_FLUX_6::Union{Float64,Missing} = missing
    MOMENTUM_FLUX_6::Union{Float64,Missing} = missing

    PARTICLE_FLUX_7::Union{Float64,Missing} = missing
    ENERGY_FLUX_7::Union{Float64,Missing} = missing
    MOMENTUM_FLUX_7::Union{Float64,Missing} = missing

    PARTICLE_FLUX_8::Union{Float64,Missing} = missing
    ENERGY_FLUX_8::Union{Float64,Missing} = missing
    MOMENTUM_FLUX_8::Union{Float64,Missing} = missing

    PARTICLE_FLUX_9::Union{Float64,Missing} = missing
    ENERGY_FLUX_9::Union{Float64,Missing} = missing
    MOMENTUM_FLUX_9::Union{Float64,Missing} = missing

    PARTICLE_FLUX_10::Union{Float64,Missing} = missing
    ENERGY_FLUX_10::Union{Float64,Missing} = missing
    MOMENTUM_FLUX_10::Union{Float64,Missing} = missing

    PARTICLE_FLUX_11::Union{Float64,Missing} = missing
    ENERGY_FLUX_11::Union{Float64,Missing} = missing
    MOMENTUM_FLUX_11::Union{Float64,Missing} = missing
end

Base.@kwdef mutable struct equilibrium_geometry
    rmin::Union{Vector{Float64},Missing} = missing
    rmaj::Union{Vector{Float64},Missing} = missing
    a::Union{Float64,Missing} = missing
    q::Union{Vector{Float64},Missing} = missing
    ftrap::Union{Vector{Float64},Missing} = missing
    Bmag2_avg::Union{Vector{Float64},Missing} = missing
    f::Union{Vector{Float64},Missing} = missing
end

include("input_neo.jl")

include("hirshman_sigmar.jl")

include("chang_hinton.jl")

"""
    run_neo(input_neo::InputNEO)

Saves input.neo file to a temporary directory, runs NEO on that directory and parses output 
"""

function run_neo(input_neo::InputNEO)
    folder = mktempdir()
    save_inputneo(input_neo, joinpath(folder, "input.neo"))

    open(joinpath(folder, "command.sh"), "w") do io
        return write(
            io,
            """
          neo -e &> command.log
          """)
    end

    run(Cmd(`bash command.sh`; dir=folder))

    ### parse outputs ###
    tmp_fluxes = Float64[]
    open(joinpath(folder, "out.neo.transport_flux"), "r") do io
        for line in eachline(io)
            if !startswith(line, "#")
                for word in split(line)
                    val = tryparse(Float64, word)
                    if val !== nothing
                        push!(tmp_fluxes, val)
                    end
                end
            end
        end
    end
    loc_first_tgyro = (4 * input_neo.N_SPECIES * 2) + 1
    tgyro_fluxes = tmp_fluxes[loc_first_tgyro:end]

    # figure out indexes
    e_index = [input_neo.N_SPECIES]
    i_index = collect(1:input_neo.N_SPECIES-1)
    particle_index(index) = 2 .+ ((index .- 1) .* 4)
    energy_index(index) = 3 .+ ((index .- 1) .* 4)
    momentum_index(index) = 4 .+ ((index .- 1) .* 4)

    # sort fluxes
    electrons_energy_flux = only(tgyro_fluxes[energy_index(e_index)])
    electrons_particle_flux = only(tgyro_fluxes[particle_index(e_index)])
    ion_particle_flux = tgyro_fluxes[particle_index(i_index)]
    ion_total_energy_flux = sum(tgyro_fluxes[energy_index(i_index)])
    ion_total_momentum_flux = sum(tgyro_fluxes[momentum_index(i_index)])

    # assign fluxes to flux_solution structure
    sol = IMAS.flux_solution(electrons_energy_flux, ion_total_energy_flux, electrons_particle_flux, ion_particle_flux, ion_total_momentum_flux)
    return sol
end

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end

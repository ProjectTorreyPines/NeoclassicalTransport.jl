module NeoclassicalTransport

using IMAS
import GACODE

include("input_neo.jl")

include("hirshman_sigmar.jl")

include("chang_hinton.jl")

include("FACIT.jl")

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

    # assign fluxes to FluxSolution structure
    T = Float64
    sol = GACODE.FluxSolution{T}(electrons_energy_flux, ion_total_energy_flux, electrons_particle_flux, ion_particle_flux, ion_total_momentum_flux)

    rm(folder; force=true, recursive=true)

    return sol
end

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end

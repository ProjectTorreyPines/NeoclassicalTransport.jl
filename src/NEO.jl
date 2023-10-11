module NEO

Base.@kwdef mutable struct flux_solution
	PARTICLE_FLUX_1::Union{Float64, Missing} = missing
	ENERGY_FLUX_1::Union{Float64, Missing} = missing
	MOMENTUM_FLUX_1::Union{Float64, Missing} = missing

	PARTICLE_FLUX_2::Union{Float64, Missing} = missing
	ENERGY_FLUX_2::Union{Float64, Missing} = missing
	MOMENTUM_FLUX_2::Union{Float64, Missing} = missing

	PARTICLE_FLUX_3::Union{Float64, Missing} = missing
	ENERGY_FLUX_3::Union{Float64, Missing} = missing
	MOMENTUM_FLUX_3::Union{Float64, Missing} = missing

	PARTICLE_FLUX_4::Union{Float64, Missing} = missing
	ENERGY_FLUX_4::Union{Float64, Missing} = missing
	MOMENTUM_FLUX_4::Union{Float64, Missing} = missing

	PARTICLE_FLUX_5::Union{Float64, Missing} = missing
	ENERGY_FLUX_5::Union{Float64, Missing} = missing
	MOMENTUM_FLUX_5::Union{Float64, Missing} = missing

	PARTICLE_FLUX_6::Union{Float64, Missing} = missing
	ENERGY_FLUX_6::Union{Float64, Missing} = missing
	MOMENTUM_FLUX_6::Union{Float64, Missing} = missing

	PARTICLE_FLUX_7::Union{Float64, Missing} = missing
	ENERGY_FLUX_7::Union{Float64, Missing} = missing
	MOMENTUM_FLUX_7::Union{Float64, Missing} = missing

	PARTICLE_FLUX_8::Union{Float64, Missing} = missing
	ENERGY_FLUX_8::Union{Float64, Missing} = missing
	MOMENTUM_FLUX_8::Union{Float64, Missing} = missing

	PARTICLE_FLUX_9::Union{Float64, Missing} = missing
	ENERGY_FLUX_9::Union{Float64, Missing} = missing
	MOMENTUM_FLUX_9::Union{Float64, Missing} = missing

	PARTICLE_FLUX_10::Union{Float64, Missing} = missing
	ENERGY_FLUX_10::Union{Float64, Missing} = missing
	MOMENTUM_FLUX_10::Union{Float64, Missing} = missing

	PARTICLE_FLUX_11::Union{Float64, Missing} = missing
	ENERGY_FLUX_11::Union{Float64, Missing} = missing
	MOMENTUM_FLUX_11::Union{Float64, Missing} = missing
end

Base.@kwdef mutable struct equilibrium_geometry
    rmin::Union{Vector{Float64}, Missing} = missing 
    rmaj::Union{Vector{Float64}, Missing} = missing 
    a::Union{Float64, Missing} = missing 
    q::Union{Vector{Float64}, Missing} = missing 
    ftrap::Union{Vector{Float64}, Missing} = missing
    Bmag2_avg::Union{Vector{Float64}, Missing} = missing 
    f::Union{Vector{Float64}, Missing} = missing
end

include("input_neo.jl")
include("hirshman_sigmar.jl")
include("chang_hinton.jl")

function run_neo(input_neo::InputNEO)
	folder = mktempdir()
	save_inputneo(input_neo, joinpath(folder, "input.neo"))

	open(joinpath(folder, "command.sh"), "w") do io
		write(
			io,
			"""
	neo -e &> command.log
	""")
	end

	run(Cmd(`bash command.sh`, dir = folder))

	### parse outputs ###
	fluxes = Float64[]

	tmp = open(joinpath(folder, "out.neo.transport_flux"), "r") do io
		for line in eachline(io)
			if !startswith(line, "#")
				for word in split(line)
					val = tryparse(Float64, word)
					if val !== nothing
						push!(fluxes, val)
					end
				end
			end
		end
	end

	loc_first_tgyro = (4 * input_neo.N_SPECIES * 2) + 1
	tgyro_fluxes = fluxes[loc_first_tgyro:end]

	flux_solution = NEO.flux_solution()

	for i in range(0, input_neo.N_SPECIES - 1)
		species = i + 1
		setfield!(flux_solution, Symbol("PARTICLE_FLUX_$species"), tgyro_fluxes[2+(i*4)])
		setfield!(flux_solution, Symbol("ENERGY_FLUX_$species"), tgyro_fluxes[3+(i*4)])
		setfield!(flux_solution, Symbol("MOMENTUM_FLUX_$species"), tgyro_fluxes[4+(i*4)])
	end

    energy_flux_electrons = getfield(flux_solution, Symbol("ENERGY_FLUX_$(input_neo.N_SPECIES)")) # electrons are always the last species in the list
    total_ion_energy_flux = -energy_flux_electrons # exclude electron energy flux from total ion energy flux
    particle_flux_electrons = getfield(flux_solution, Symbol("PARTICLE_FLUX_$(input_neo.N_SPECIES)"))

    for field in fieldnames(NEO.flux_solution)
        if ismissing(getfield(flux_solution, field))
            setfield!(flux_solution, field, 0.0)
        end

        if occursin("ENERGY", string(field))
            total_ion_energy_flux += getfield(flux_solution, field) 
        end
    end

    total_fluxes = IMAS.flux_solution(particle_flux_electrons, 0.0, energy_flux_electrons, total_ion_energy_flux)

    return total_fluxes
end

end

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

include("input_neo.jl")

end # module NEO

# Input handling for the native NEO solve: neo_parse.py defaults, the
# input.neo key=value parser, the InputNEO mapping, and the local-profile
# setup of neo_make_profiles.f90 (case profile_model=1) plus the checks of
# neo_check.f90 that apply to the supported subset.

const NEO_SPECIES_MAX = 11

"""
    neo_input_defaults() -> Dict{String,Any}

Defaults of every `input.neo` key the native solver reads, as set by
`\$GACODE_ROOT/neo/bin/neo_parse.py` (integer-valued keys are `Int`).
"""
function neo_input_defaults()
    d = Dict{String,Any}(
        "N_ENERGY" => 6, "N_XI" => 17, "N_THETA" => 17, "N_RADIAL" => 1,
        "RMIN_OVER_A" => 0.5, "RMAJ_OVER_A" => 3.0,
        "SIM_MODEL" => 2, "EQUILIBRIUM_MODEL" => 0, "COLLISION_MODEL" => 4,
        "PROFILE_MODEL" => 1, "IPCCW" => -1, "BTCCW" => -1,
        "ROTATION_MODEL" => 1, "OMEGA_ROT" => 0.0, "OMEGA_ROT_DERIV" => 0.0,
        "SPITZER_MODEL" => 0, "COLL_UNCOUPLEDEI_MODEL" => 0, "COLL_UNCOUPLEDANISO_MODEL" => 0,
        "AE_FLAG" => 0, "DENS_AE" => 1.0, "TEMP_AE" => 1.0, "DLNNDR_AE" => 1.0, "DLNTDR_AE" => 1.0,
        "N_SPECIES" => 1, "NU_1" => 0.1,
        "DPHI0DR" => 0.0, "EPAR0" => 0.0, "Q" => 2.0, "RHO_STAR" => 0.001,
        "SHEAR" => 1.0, "SHIFT" => 0.0, "ZMAG_OVER_A" => 0.0, "S_ZMAG" => 0.0,
        "KAPPA" => 1.0, "S_KAPPA" => 0.0, "DELTA" => 0.0, "S_DELTA" => 0.0,
        "ZETA" => 0.0, "S_ZETA" => 0.0, "BETA_STAR" => 0.0,
        "THREED_MODEL" => 0, "LAGUERRE_METHOD" => 1, "WRITE_CMOMENTS_FLAG" => 0)
    for i in 1:NEO_SPECIES_MAX
        d["Z_$i"] = 1
        d["MASS_$i"] = 1.0
        d["DENS_$i"] = i == 1 ? 1.0 : 0.0
        d["TEMP_$i"] = 1.0
        d["DLNNDR_$i"] = 1.0
        d["DLNTDR_$i"] = 1.0
        d["ANISO_MODEL_$i"] = 1
    end
    for k in 3:6
        d["SHAPE_SIN$k"] = 0.0
        d["SHAPE_S_SIN$k"] = 0.0
    end
    for k in 0:6
        d["SHAPE_COS$k"] = 0.0
        d["SHAPE_S_COS$k"] = 0.0
    end
    return d
end

# InputNEO uses a few non-NEO field names for the adiabatic-electron inputs
const INPUTNEO_KEY_MAP = Dict("NE_ADE" => "DENS_AE", "TE_ADE" => "TEMP_AE", "DLNNDRE_ADE" => "DLNNDR_AE", "DLNTDRE_ADE" => "DLNTDR_AE")

function _set_input_key!(d::Dict{String,Any}, key::AbstractString, value)
    key = get(INPUTNEO_KEY_MAP, key, key)
    haskey(d, key) || return d # keys the native solver does not use
    d[key] = d[key] isa Int ? Int(value) : value
    return d
end

"""
    read_input_neo(path) -> Dict{String,Any}

Parse an `input.neo` file (`KEY=value` lines, `#` comments) on top of
[`neo_input_defaults`](@ref). Keys the native solver does not read are ignored.
"""
function read_input_neo(path::AbstractString)
    d = neo_input_defaults()
    for line in eachline(path)
        line = strip(first(split(line, '#'; limit=2)))
        isempty(line) && continue
        kv = split(line, '='; limit=2)
        length(kv) == 2 || error("read_input_neo: cannot parse line '$line' in $path")
        key = strip(kv[1])
        value = strip(kv[2])
        if haskey(d, key)
            _set_input_key!(d, key, d[key] isa Int ? parse(Int, value) : parse(Float64, value))
        end
    end
    return d
end

"""
    input_dict(input_neo::InputNEO) -> Dict{String,Any}

The non-missing fields of an [`InputNEO`](@ref) on top of [`neo_input_defaults`](@ref).
"""
function input_dict(input_neo::InputNEO)
    d = neo_input_defaults()
    for field in fieldnames(typeof(input_neo))
        value = getfield(input_neo, field)
        ismissing(value) && continue
        _set_input_key!(d, String(field), value)
    end
    return d
end

"""
    NEOParams{T}

Everything the native NEO solve needs for one flux surface, in NEO's own
normalizations (`neo_make_profiles.f90`, `profile_model=1`): grid sizes and
model switches, the signed geometry scalars, and the per-species profile
vectors. Species index 1 is the reference species (`NU_1`, `z_1`); electrons
are wherever `z < 0` (`is_ele`, or `-1` with adiabatic electrons).

Build one with `NEOParams(::InputNEO)`, `NEOParams(path)` (an `input.neo`
file) or `NEOParams(::AbstractDict)`.
"""
struct NEOParams{T<:Real}
    # sizes and switches
    n_species::Int
    n_energy::Int
    n_xi::Int
    n_theta::Int
    sim_model::Int
    equilibrium_model::Int
    collision_model::Int
    rotation_model::Int
    ae_flag::Int
    laguerre_method::Int
    is_ele::Int
    sign_q::Float64
    sign_bunit::Float64
    # geometry (local)
    rmin::T
    rmaj::T
    q::T        # abs(Q) * sign_q
    rho::T      # abs(RHO_STAR) * sign_bunit
    shear::T
    shift::T
    zmag::T
    s_zmag::T
    kappa::T
    s_kappa::T
    delta::T
    s_delta::T
    zeta::T
    s_zeta::T
    shape_sin::OffsetVector{T,Vector{T}}   # 0:6 (only 3:6 used)
    shape_s_sin::OffsetVector{T,Vector{T}}
    shape_cos::OffsetVector{T,Vector{T}}   # 0:6
    shape_s_cos::OffsetVector{T,Vector{T}}
    beta_star::T
    # fields and rotation
    dphi0dr::T          # zeroed when rotation_model == 2
    epar0::T
    omega_rot::T        # zeroed when rotation_model == 1
    omega_rot_deriv::T
    # adiabatic electrons
    dens_ae::T
    temp_ae::T
    dlnndr_ae::T
    dlntdr_ae::T
    nu_1::T
    # species
    z::Vector{T}
    mass::Vector{T}
    dens::Vector{T}
    temp::Vector{T}
    dlnndr::Vector{T}
    dlntdr::Vector{T}
    nu::Vector{T}
    vth::Vector{T}
end

NEOParams(path::AbstractString) = NEOParams(read_input_neo(path))
NEOParams(input_neo::InputNEO{T}) where {T<:Real} = NEOParams{T}(input_dict(input_neo))
NEOParams(d::AbstractDict) = NEOParams{Float64}(d)

function NEOParams{T}(d::AbstractDict) where {T<:Real}
    n_species = d["N_SPECIES"]
    n_energy = d["N_ENERGY"]
    n_xi = d["N_XI"]
    n_theta = d["N_THETA"]

    # ---- supported subset (neo_check.f90 + what this port implements)
    1 <= n_species <= NEO_SPECIES_MAX || error("NEONative: n_species must be in 1:$NEO_SPECIES_MAX")
    isodd(n_theta) || error("NEONative: n_theta must be odd")
    d["N_RADIAL"] == 1 || error("NEONative: only n_radial=1 (local) is supported")
    d["PROFILE_MODEL"] == 1 || error("NEONative: only profile_model=1 (local input.neo profiles) is supported")
    d["THREED_MODEL"] == 0 || error("NEONative: threed_model=1 is not supported")
    d["SPITZER_MODEL"] == 0 || error("NEONative: spitzer_model=1 is not supported")
    d["COLL_UNCOUPLEDEI_MODEL"] == 0 || error("NEONative: coll_uncoupledei_model != 0 is not supported")
    d["COLL_UNCOUPLEDANISO_MODEL"] == 0 || error("NEONative: coll_uncoupledaniso_model != 0 is not supported")
    sim_model = d["SIM_MODEL"]
    sim_model in (1, 2) || error("NEONative: sim_model=$sim_model has no kinetic solve (only 1 and 2 are supported)")
    equilibrium_model = d["EQUILIBRIUM_MODEL"]
    equilibrium_model in (0, 1, 2) || error("NEONative: equilibrium_model=$equilibrium_model invalid")
    collision_model = d["COLLISION_MODEL"]
    collision_model in 1:5 || error("NEONative: collision_model=$collision_model invalid (1 Connor, 2 reduced HS, 3 full HS, 4 full FP, 5 FP + ad-hoc field)")
    rotation_model = d["ROTATION_MODEL"]
    rotation_model in (1, 2) || error("NEONative: rotation_model=$rotation_model invalid")
    laguerre_method = d["LAGUERRE_METHOD"]
    laguerre_method in 1:4 || error("NEONative: laguerre_method=$laguerre_method invalid")
    for is in 1:n_species
        d["ANISO_MODEL_$is"] == 1 || error("NEONative: anisotropic species (aniso_model=2) are not supported")
    end
    d["RHO_STAR"] >= 0 || error("NEONative: rho_unit must be positive")

    # ---- neo_make_profiles.f90, case profile_model=1
    ae_flag = d["AE_FLAG"]
    if n_species == 1
        ae_flag = 1
    end
    sign_bunit = d["BTCCW"] > 0 ? -1.0 : 1.0
    sign_q = d["IPCCW"] > 0 ? -sign_bunit : sign_bunit

    z = T[d["Z_$is"] for is in 1:n_species]
    mass = T[d["MASS_$is"] for is in 1:n_species]
    dens = T[d["DENS_$is"] for is in 1:n_species]
    temp = T[d["TEMP_$is"] for is in 1:n_species]
    dlnndr = T[d["DLNNDR_$is"] for is in 1:n_species]
    dlntdr = T[d["DLNTDR_$is"] for is in 1:n_species]
    nu_1 = T(d["NU_1"])
    nu = T[nu_1 * z[is]^4 / z[1]^4 * dens[is] / dens[1] * sqrt(mass[1] / mass[is]) * (temp[1] / temp[is])^1.5 for is in 1:n_species]
    vth = T[sqrt(temp[is] / mass[is]) for is in 1:n_species]

    for is in 1:n_species
        dens[is] > 0 || error("NEONative: density must be positive (species $is)")
        temp[is] > 0 || error("NEONative: temperature must be positive (species $is)")
        nu[is] > 0 || error("NEONative: collision frequency must be positive (species $is)")
        abs(z[is]) > eps() || error("NEONative: charge must be non-zero (species $is)")
    end

    num_ele = count(<(0), z)
    is_ele = -1
    if ae_flag == 1
        num_ele == 0 || error("NEONative: electron species specified with adiabatic electron flag")
    else
        num_ele == 0 && error("NEONative: no electron species specified")
        is_ele = findfirst(<(0), z)
    end
    num_ele <= 1 || error("NEONative: only one electron species allowed")

    dphi0dr = T(d["DPHI0DR"])
    omega_rot = T(d["OMEGA_ROT"])
    omega_rot_deriv = T(d["OMEGA_ROT_DERIV"])
    if rotation_model == 2
        # strong-rotation limit: phi_(-1) (set by omega_rot) replaces phi_(0)
        dphi0dr = zero(T)
    else
        # ROT_solve_phi zeroes these for rotation_model=1
        omega_rot = zero(T)
        omega_rot_deriv = zero(T)
    end

    shape_sin = OffsetVector(T[k >= 3 ? d["SHAPE_SIN$k"] : 0.0 for k in 0:6], 0:6)
    shape_s_sin = OffsetVector(T[k >= 3 ? d["SHAPE_S_SIN$k"] : 0.0 for k in 0:6], 0:6)
    shape_cos = OffsetVector(T[d["SHAPE_COS$k"] for k in 0:6], 0:6)
    shape_s_cos = OffsetVector(T[d["SHAPE_S_COS$k"] for k in 0:6], 0:6)

    return NEOParams{T}(
        n_species, n_energy, n_xi, n_theta,
        sim_model, equilibrium_model, collision_model, rotation_model, ae_flag, laguerre_method, is_ele,
        sign_q, sign_bunit,
        d["RMIN_OVER_A"], d["RMAJ_OVER_A"], abs(d["Q"]) * sign_q, abs(d["RHO_STAR"]) * sign_bunit,
        d["SHEAR"], d["SHIFT"], d["ZMAG_OVER_A"], d["S_ZMAG"],
        d["KAPPA"], d["S_KAPPA"], d["DELTA"], d["S_DELTA"], d["ZETA"], d["S_ZETA"],
        shape_sin, shape_s_sin, shape_cos, shape_s_cos, d["BETA_STAR"],
        dphi0dr, d["EPAR0"], omega_rot, omega_rot_deriv,
        d["DENS_AE"], d["TEMP_AE"], d["DLNNDR_AE"], d["DLNTDR_AE"], nu_1,
        z, mass, dens, temp, dlnndr, dlntdr, nu, vth)
end

Base.eltype(::Type{NEOParams{T}}) where {T} = T
Base.eltype(p::NEOParams) = eltype(typeof(p))

"""
    n_row(p::NEOParams)

Size of the kinetic system: `n_species*(n_energy+1)*(n_xi+1)*n_theta`.
"""
n_row(p::NEOParams) = p.n_species * (p.n_energy + 1) * (p.n_xi + 1) * p.n_theta

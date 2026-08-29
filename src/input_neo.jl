mutable struct InputNEO{T<:Real}
    BETA_STAR::Union{T,Missing}
    DELTA::Union{T,Missing}
    DLNNDRE_ADE::Union{T,Missing}
    DLNTDRE_ADE::Union{T,Missing}
    DPHI0DR::Union{T,Missing}
    EPAR0::Union{T,Missing}
    EPAR0_SPITZER::Union{T,Missing}
    KAPPA::Union{T,Missing}
    NE_ADE::Union{T,Missing}
    NU_1::Union{T,Missing}
    OMEGA_ROT::Union{T,Missing}
    OMEGA_ROT_DERIV::Union{T,Missing}
    Q::Union{T,Missing}
    RHO_STAR::Union{T,Missing}
    RMAJ_OVER_A::Union{T,Missing}
    RMIN_OVER_A::Union{T,Missing}
    S_DELTA::Union{T,Missing}
    S_KAPPA::Union{T,Missing}
    S_ZETA::Union{T,Missing}
    S_ZMAG::Union{T,Missing}
    SHEAR::Union{T,Missing}
    SHIFT::Union{T,Missing}
    TE_ADE::Union{T,Missing}
    THREED_EXB_DPHI0DR::Union{T,Missing}
    ZETA::Union{T,Missing}
    ZMAG_OVER_A::Union{T,Missing}


    #moment parameters 
    SHAPE_COS0::Union{T,Missing}
    SHAPE_S_COS0::Union{T,Missing}
    SHAPE_COS1::Union{T,Missing}
    SHAPE_S_COS1::Union{T,Missing}
    SHAPE_COS2::Union{T,Missing}
    SHAPE_S_COS2::Union{T,Missing}
    SHAPE_COS3::Union{T,Missing}
    SHAPE_S_COS3::Union{T,Missing}
    SHAPE_SIN3::Union{T,Missing}
    SHAPE_S_SIN3::Union{T,Missing}


    #species-specific parameters 
    ANISO_MODEL_1::Union{Int,Missing}
    DENS_1::Union{T,Missing}
    DLNNDR_1::Union{T,Missing}
    DLNTDR_1::Union{T,Missing}
    DLNTDR_PARA_1::Union{T,Missing}
    DLNTDR_PERP_1::Union{T,Missing}
    MASS_1::Union{T,Missing}
    PROFILE_DLNNDR_1_SCALE::Union{T,Missing}
    PROFILE_DLNTDR_1_SCALE::Union{T,Missing}
    TEMP_1::Union{T,Missing}
    TEMP_PARA_1::Union{T,Missing}
    TEMP_PERP_1::Union{T,Missing}
    Z_1::Union{Int64,Missing}

    ANISO_MODEL_2::Union{Int,Missing}
    DENS_2::Union{T,Missing}
    DLNNDR_2::Union{T,Missing}
    DLNTDR_2::Union{T,Missing}
    DLNTDR_PARA_2::Union{T,Missing}
    DLNTDR_PERP_2::Union{T,Missing}
    MASS_2::Union{T,Missing}
    PROFILE_DLNNDR_2_SCALE::Union{T,Missing}
    PROFILE_DLNTDR_2_SCALE::Union{T,Missing}
    TEMP_2::Union{T,Missing}
    TEMP_PARA_2::Union{T,Missing}
    TEMP_PERP_2::Union{T,Missing}
    Z_2::Union{Int64,Missing}

    ANISO_MODEL_3::Union{Int,Missing}
    DENS_3::Union{T,Missing}
    DLNNDR_3::Union{T,Missing}
    DLNTDR_3::Union{T,Missing}
    DLNTDR_PARA_3::Union{T,Missing}
    DLNTDR_PERP_3::Union{T,Missing}
    MASS_3::Union{T,Missing}
    PROFILE_DLNNDR_3_SCALE::Union{T,Missing}
    PROFILE_DLNTDR_3_SCALE::Union{T,Missing}
    TEMP_3::Union{T,Missing}
    TEMP_PARA_3::Union{T,Missing}
    TEMP_PERP_3::Union{T,Missing}
    Z_3::Union{Int64,Missing}

    ANISO_MODEL_4::Union{Int,Missing}
    DENS_4::Union{T,Missing}
    DLNNDR_4::Union{T,Missing}
    DLNTDR_4::Union{T,Missing}
    DLNTDR_PARA_4::Union{T,Missing}
    DLNTDR_PERP_4::Union{T,Missing}
    MASS_4::Union{T,Missing}
    PROFILE_DLNNDR_4_SCALE::Union{T,Missing}
    PROFILE_DLNTDR_4_SCALE::Union{T,Missing}
    TEMP_4::Union{T,Missing}
    TEMP_PARA_4::Union{T,Missing}
    TEMP_PERP_4::Union{T,Missing}
    Z_4::Union{Int64,Missing}

    ANISO_MODEL_5::Union{Int,Missing}
    DENS_5::Union{T,Missing}
    DLNNDR_5::Union{T,Missing}
    DLNTDR_5::Union{T,Missing}
    DLNTDR_PARA_5::Union{T,Missing}
    DLNTDR_PERP_5::Union{T,Missing}
    MASS_5::Union{T,Missing}
    PROFILE_DLNNDR_5_SCALE::Union{T,Missing}
    PROFILE_DLNTDR_5_SCALE::Union{T,Missing}
    TEMP_5::Union{T,Missing}
    TEMP_PARA_5::Union{T,Missing}
    TEMP_PERP_5::Union{T,Missing}
    Z_5::Union{Int64,Missing}

    ANISO_MODEL_6::Union{Int,Missing}
    DENS_6::Union{T,Missing}
    DLNNDR_6::Union{T,Missing}
    DLNTDR_6::Union{T,Missing}
    DLNTDR_PARA_6::Union{T,Missing}
    DLNTDR_PERP_6::Union{T,Missing}
    MASS_6::Union{T,Missing}
    PROFILE_DLNNDR_6_SCALE::Union{T,Missing}
    PROFILE_DLNTDR_6_SCALE::Union{T,Missing}
    TEMP_6::Union{T,Missing}
    TEMP_PARA_6::Union{T,Missing}
    TEMP_PERP_6::Union{T,Missing}
    Z_6::Union{Int64,Missing}

    ANISO_MODEL_7::Union{Int,Missing}
    DENS_7::Union{T,Missing}
    DLNNDR_7::Union{T,Missing}
    DLNTDR_7::Union{T,Missing}
    DLNTDR_PARA_7::Union{T,Missing}
    DLNTDR_PERP_7::Union{T,Missing}
    MASS_7::Union{T,Missing}
    PROFILE_DLNNDR_7_SCALE::Union{T,Missing}
    PROFILE_DLNTDR_7_SCALE::Union{T,Missing}
    TEMP_7::Union{T,Missing}
    TEMP_PARA_7::Union{T,Missing}
    TEMP_PERP_7::Union{T,Missing}
    Z_7::Union{Int64,Missing}

    ANISO_MODEL_8::Union{Int,Missing}
    DENS_8::Union{T,Missing}
    DLNNDR_8::Union{T,Missing}
    DLNTDR_8::Union{T,Missing}
    DLNTDR_PARA_8::Union{T,Missing}
    DLNTDR_PERP_8::Union{T,Missing}
    MASS_8::Union{T,Missing}
    PROFILE_DLNNDR_8_SCALE::Union{T,Missing}
    PROFILE_DLNTDR_8_SCALE::Union{T,Missing}
    TEMP_8::Union{T,Missing}
    TEMP_PARA_8::Union{T,Missing}
    TEMP_PERP_8::Union{T,Missing}
    Z_8::Union{Int64,Missing}

    ANISO_MODEL_9::Union{Int,Missing}
    DENS_9::Union{T,Missing}
    DLNNDR_9::Union{T,Missing}
    DLNTDR_9::Union{T,Missing}
    DLNTDR_PARA_9::Union{T,Missing}
    DLNTDR_PERP_9::Union{T,Missing}
    MASS_9::Union{T,Missing}
    PROFILE_DLNNDR_9_SCALE::Union{T,Missing}
    PROFILE_DLNTDR_9_SCALE::Union{T,Missing}
    TEMP_9::Union{T,Missing}
    TEMP_PARA_9::Union{T,Missing}
    TEMP_PERP_9::Union{T,Missing}
    Z_9::Union{Int64,Missing}

    ANISO_MODEL_10::Union{Int,Missing}
    DENS_10::Union{T,Missing}
    DLNNDR_10::Union{T,Missing}
    DLNTDR_10::Union{T,Missing}
    DLNTDR_PARA_10::Union{T,Missing}
    DLNTDR_PERP_10::Union{T,Missing}
    MASS_10::Union{T,Missing}
    PROFILE_DLNNDR_10_SCALE::Union{T,Missing}
    PROFILE_DLNTDR_10_SCALE::Union{T,Missing}
    TEMP_10::Union{T,Missing}
    TEMP_PARA_10::Union{T,Missing}
    TEMP_PERP_10::Union{T,Missing}
    Z_10::Union{Int64,Missing}

    ANISO_MODEL_11::Union{Int,Missing}
    DENS_11::Union{T,Missing}
    DLNNDR_11::Union{T,Missing}
    DLNTDR_11::Union{T,Missing}
    DLNTDR_PARA_11::Union{T,Missing}
    DLNTDR_PERP_11::Union{T,Missing}
    MASS_11::Union{T,Missing}
    PROFILE_DLNNDR_11_SCALE::Union{T,Missing}
    PROFILE_DLNTDR_11_SCALE::Union{T,Missing}
    TEMP_11::Union{T,Missing}
    TEMP_PARA_11::Union{T,Missing}
    TEMP_PERP_11::Union{T,Missing}
    Z_11::Union{Int64,Missing}

    #switches 
    BTCCW::Union{Int,Missing}
    COLLISION_MODEL::Union{Int,Missing}
    EQUILIBRIUM_MODEL::Union{Int,Missing}
    IPCCW::Union{Int,Missing}
    N_ENERGY::Union{Int,Missing}
    N_RADIAL::Union{Int,Missing}
    N_SPECIES::Union{Int,Missing}
    N_THETA::Union{Int,Missing}
    N_XI::Union{Int,Missing}
    PROFILE_EQUILIBRIUM_MODEL::Union{Int,Missing}
    PROFILE_ERAD0_MODEL::Union{Int,Missing}
    PROFILE_MODEL::Union{Int,Missing}
    ROTATION_MODEL::Union{Int,Missing}
    SILENT_FLAG::Union{Int,Missing}
    SIM_MODEL::Union{Int,Missing}
    SPITZER_MODEL::Union{Int,Missing}
    THREED_MODEL::Union{Int,Missing}
    THREED_EXB_MODEL::Union{Int,Missing}
end

# Keyword constructors, generated so that every field defaults to `missing`.
# (`Base.@kwdef` cannot be used on this struct: for a parametric type it also
# emits a parameter-free keyword method that routes through the positional
# UnionAll constructor, which Julia cannot resolve here because no field has
# bare type `T` — they are all `Union{T,Missing}` — and that method cannot be
# replaced without a method-overwrite error during precompilation.)
#
# The type parameter exists so the plasma quantities can carry
# `ForwardDiff.Dual`s: see [`InputNEO{T}(::InputNEO)`](@ref) and the AD-capable
# `run_neonn` path. Integer fields (species charges, switches) stay `Int`
# whatever `T` is — nothing is differentiated with respect to them.
@eval function (::Type{InputNEO{T}})(; $(map(f -> Expr(:kw, f, :missing), fieldnames(InputNEO{Float64}))...)) where {T<:Real}
    return InputNEO{T}($(fieldnames(InputNEO{Float64})...))
end

"""
    InputNEO(; kwargs...)

`InputNEO{Float64}` — the element type an unparameterized construction defaults to.
"""
InputNEO(; kwargs...) = InputNEO{Float64}(; kwargs...)

Base.eltype(::Type{InputNEO{T}}) where {T<:Real} = T
Base.eltype(input_neo::InputNEO) = eltype(typeof(input_neo))

"""
    InputNEO{T}(other::InputNEO)

Copy `other` with every real-valued field converted to `T` (`missing` stays
`missing`, integer fields stay `Int`). This is how a differentiable evaluation
starts from an already-populated input:

```julia
ineo0 = InputNEO(eqt, cp1d, gridpoint)
function Qi(x::AbstractVector{D}) where {D<:Real}
    ineo = InputNEO{D}(ineo0)
    ineo.DLNTDR_1 = x[1]
    ineo.TEMP_1 = x[2]
    return run_neonn(ineo; warn_nn_train_bounds=false).ENERGY_FLUX_i
end
ForwardDiff.gradient(Qi, [ineo0.DLNTDR_1, ineo0.TEMP_1])
```
"""
function InputNEO{T}(other::InputNEO) where {T<:Real}
    input_neo = InputNEO{T}()
    for field in fieldnames(InputNEO{T})
        value = getfield(other, field)
        ismissing(value) && continue
        setfield!(input_neo, field, value isa Integer ? value : convert(T, value))
    end
    return input_neo
end

InputNEO{T}(other::InputNEO{T}) where {T<:Real} = deepcopy(other)

"""
    InputNEO(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, gridpoint_cp)

Populates InputNEO structure with quantities from eqt and cp1d using NEO normalizations.

The element type is inherited from `eqt`/`cp1d`, so a dd carrying `ForwardDiff.Dual`s
produces an `InputNEO{<:Dual}` and the NEO-NN surrogates differentiate straight through.
"""
function InputNEO(eqt::IMAS.equilibrium__time_slice{Teq}, cp1d::IMAS.core_profiles__profiles_1d{Tcp}, gridpoint_cp) where {Teq<:Real,Tcp<:Real}
    # element type follows the dd's: under FUSE's AD path both carry Duals, and
    # every quantity below is built from them
    input_neo = InputNEO{promote_type(Teq, Tcp)}()

    eqt1d = eqt.profiles_1d
    ions = cp1d.ion

    e = IMAS.cgs.e # statcoul
    k = IMAS.cgs.k # erg/eV
    mp = IMAS.cgs.mp # g
    me = IMAS.cgs.me # g
    md = IMAS.cgs.md # g
    m_to_cm = IMAS.cgs.m_to_cm
    m³_to_cm³ = IMAS.cgs.m³_to_cm³

    rmin = GACODE.r_min_core_profiles(eqt1d, cp1d.grid.rho_tor_norm)
    a = rmin[end]

    temp_1 = ions[1].temperature
    T1 = temp_1[gridpoint_cp]
    dens_1 = ions[1].density[gridpoint_cp] ./ m³_to_cm³

    dens_e = cp1d.electrons.density ./ m³_to_cm³
    dlnnedr = -IMAS.calc_z(rmin ./ a, dens_e, :backward)
    ne = dens_e[gridpoint_cp]
    dlnnedr = dlnnedr[gridpoint_cp]

    temp_e = cp1d.electrons.temperature
    dlntedr = -IMAS.calc_z(rmin ./ a, temp_e, :backward)
    Te = temp_e[gridpoint_cp]
    dlntedr = dlntedr[gridpoint_cp]

    n_norm = ne
    t_norm = Te
    v_norm = sqrt(k .* t_norm ./ md)

    input_neo.RMIN_OVER_A = rmin[gridpoint_cp] / a

    input_neo.DELTA = IMAS.interp1d(eqt1d.rho_tor_norm, 0.5 * (eqt1d.triangularity_lower + eqt1d.triangularity_upper)).(cp1d.grid.rho_tor_norm)[gridpoint_cp]

    kappa = IMAS.interp1d(eqt1d.rho_tor_norm, eqt1d.elongation).(cp1d.grid.rho_tor_norm)
    input_neo.KAPPA = kappa[gridpoint_cp]

    loglam = IMAS.lnΛ_ei(cp1d.electrons.density[gridpoint_cp], cp1d.electrons.temperature[gridpoint_cp])
    Z1 = IMAS.avgZ(ions[1].element[1].z_n, T1)
    m1 = ions[1].element[1].a * mp
    nu1 = sqrt(2) * pi * dens_1 * Z1^4.0 * e^4.0 * loglam ./ (sqrt(m1) * (k * temp_1) .^ 1.5)

    input_neo.NU_1 = (nu1./(v_norm./a))[gridpoint_cp]

    w0 = -cp1d.rotation_frequency_tor_sonic[gridpoint_cp]
    w0p = -IMAS.gradient(rmin, cp1d.rotation_frequency_tor_sonic)[gridpoint_cp]

    input_neo.OMEGA_ROT = w0 / (v_norm / a)
    input_neo.OMEGA_ROT_DERIV = w0p * a^2 / v_norm

    ####################

    q_profile = IMAS.interp1d(eqt1d.rho_tor_norm, eqt1d.q).(cp1d.grid.rho_tor_norm)
    q = q_profile[gridpoint_cp]
    input_neo.Q = abs(q)

    input_neo.RHO_STAR = GACODE.rho_s(cp1d, eqt)[gridpoint_cp] / a

    Rmaj = IMAS.interp1d(eqt1d.rho_tor_norm, m_to_cm * 0.5 * (eqt1d.r_outboard .+ eqt1d.r_inboard)).(cp1d.grid.rho_tor_norm)
    input_neo.RMAJ_OVER_A = Rmaj[gridpoint_cp] / a

    delta = IMAS.interp1d(eqt1d.rho_tor_norm, 0.5 * (eqt1d.triangularity_lower + eqt1d.triangularity_upper)).(cp1d.grid.rho_tor_norm)
    sdelta = rmin .* IMAS.gradient(rmin, delta)
    input_neo.S_DELTA = sdelta[gridpoint_cp]

    skappa = rmin .* IMAS.gradient(rmin, kappa) ./ kappa
    input_neo.S_KAPPA = skappa[gridpoint_cp]

    zeta =
        IMAS.interp1d(
            eqt1d.rho_tor_norm,
            0.25 * (eqt1d.squareness_lower_inner .+ eqt1d.squareness_lower_outer .+ eqt1d.squareness_upper_inner .+ eqt1d.squareness_upper_outer)
        ).(cp1d.grid.rho_tor_norm)
    input_neo.ZETA = zeta[gridpoint_cp]
    szeta = rmin .* IMAS.gradient(rmin, zeta)
    input_neo.S_ZETA = szeta[gridpoint_cp]

    dqdr = IMAS.gradient(rmin, q_profile)[gridpoint_cp]
    input_neo.SHEAR = rmin[gridpoint_cp] / q * dqdr

    drmaj = IMAS.gradient(rmin, Rmaj)
    input_neo.SHIFT = drmaj[gridpoint_cp]

    Z0 = IMAS.interp1d(eqt1d.rho_tor_norm, eqt1d.geometric_axis.z * 1e2).(cp1d.grid.rho_tor_norm)
    input_neo.ZMAG_OVER_A = Z0[gridpoint_cp] / a
    sZ0 = IMAS.gradient(rmin, Z0)
    input_neo.S_ZMAG = sZ0[gridpoint_cp]

    # setproperty! (not setfield!) so every RHS is converted to the field type:
    # under AD the struct is InputNEO{Dual} and constants like me/md are Float64
    for iion in eachindex(ions)
        species = iion
        setproperty!(input_neo, Symbol("ANISO_MODEL_$species"), 1)
        setproperty!(input_neo, Symbol("MASS_$species"), ions[iion].element[1].a .* mp / md)
        setproperty!(input_neo, Symbol("Z_$species"), Int(ions[iion].element[1].z_n / ions[1].element[1].z_n))

        Ti = ions[iion].temperature ./ t_norm
        dlntidr = -IMAS.calc_z(rmin ./ a, Ti, :backward)
        Ti = Ti[gridpoint_cp]
        dlntidr = dlntidr[gridpoint_cp]

        ni = ions[iion].density ./ m³_to_cm³ / n_norm
        dlnnidr = -IMAS.calc_z(rmin ./ a, ni, :backward)
        ni = ni[gridpoint_cp]
        dlnnidr = dlnnidr[gridpoint_cp]

        setproperty!(input_neo, Symbol("TEMP_$species"), Ti)
        setproperty!(input_neo, Symbol("DENS_$species"), ni)
        setproperty!(input_neo, Symbol("DLNNDR_$species"), dlnnidr)
        setproperty!(input_neo, Symbol("DLNTDR_$species"), dlntidr)
    end

    for i in range(1, 11)
        density_val = getfield(input_neo, Symbol("DENS_$i"))
        if ismissing(density_val)
            setproperty!(input_neo, Symbol("DENS_$i"), ne / n_norm)
            setproperty!(input_neo, Symbol("TEMP_$i"), Te / t_norm)
            setproperty!(input_neo, Symbol("ANISO_MODEL_$i"), 1)
            setproperty!(input_neo, Symbol("MASS_$i"), me / md)
            setproperty!(input_neo, Symbol("Z_$i"), -1)
            setproperty!(input_neo, Symbol("DLNNDR_$i"), dlnnedr)
            setproperty!(input_neo, Symbol("DLNTDR_$i"), dlntedr)
            break
        end
    end

    # fix to PROFILE_MODEL 1, N_RADIAL must always be 1
    input_neo.PROFILE_MODEL = 1
    input_neo.N_RADIAL = 1

    input_neo.THREED_MODEL = 0
    input_neo.EQUILIBRIUM_MODEL = 2
    input_neo.ROTATION_MODEL = 2

    input_neo.N_SPECIES = length(ions) + 1 # add 1 to include electrons

    # setting sign conventions
    Bt = eqt.global_quantities.vacuum_toroidal_field.b0
    input_neo.BTCCW = sign(Bt)
    input_neo.IPCCW = sign(Bt) * sign(q)

    return input_neo
end

function save_inputneo(input_neo::InputNEO, filename::String)
    eltype(input_neo) <: ForwardDiff.Dual && error(
        "save_inputneo/run_neo are not differentiable (NEO is an external code): " *
        "evaluate them on the primal InputNEO{Float64}, not on an InputNEO{Dual}")
    open(filename, "w") do io
        for key in fieldnames(typeof(input_neo))
            if startswith(String(key), "_")
                continue
            end
            try
                value = getfield(input_neo, key)
                if ismissing(value)
                    continue
                elseif isa(value, Int)
                    println(io, "$(key)=$(convert(Int, value))")
                else
                    println(io, "$(key)=$(convert(Float64, value))")
                end
            catch
                println("Error writing $key to input.neo file")
                rethrow()
            end
        end
    end
end

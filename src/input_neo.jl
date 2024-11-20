using IMAS

Base.@kwdef mutable struct InputNEO
    BETA_STAR::Union{Float64,Missing} = missing
    DELTA::Union{Float64,Missing} = missing
    DLNNDRE_ADE::Union{Float64,Missing} = missing
    DLNTDRE_ADE::Union{Float64,Missing} = missing
    DPHI0DR::Union{Float64,Missing} = missing
    EPAR0::Union{Float64,Missing} = missing
    EPAR0_SPITZER::Union{Float64,Missing} = missing
    KAPPA::Union{Float64,Missing} = missing
    NE_ADE::Union{Float64,Missing} = missing
    NU_1::Union{Float64,Missing} = missing
    OMEGA_ROT::Union{Float64,Missing} = missing
    OMEGA_ROT_DERIV::Union{Float64,Missing} = missing
    Q::Union{Float64,Missing} = missing
    RHO_STAR::Union{Float64,Missing} = missing
    RMAJ_OVER_A::Union{Float64,Missing} = missing
    RMIN_OVER_A::Union{Float64,Missing} = missing
    S_DELTA::Union{Float64,Missing} = missing
    S_KAPPA::Union{Float64,Missing} = missing
    S_ZETA::Union{Float64,Missing} = missing
    S_ZMAG::Union{Float64,Missing} = missing
    SHEAR::Union{Float64,Missing} = missing
    SHIFT::Union{Float64,Missing} = missing
    TE_ADE::Union{Float64,Missing} = missing
    THREED_EXB_DPHI0DR::Union{Float64,Missing} = missing
    ZETA::Union{Float64,Missing} = missing
    ZMAG_OVER_A::Union{Float64,Missing} = missing


    #moment parameters 
    SHAPE_COS0::Union{Float64,Missing} = missing
    SHAPE_S_COS0::Union{Float64,Missing} = missing
    SHAPE_COS1::Union{Float64,Missing} = missing
    SHAPE_S_COS1::Union{Float64,Missing} = missing
    SHAPE_COS2::Union{Float64,Missing} = missing
    SHAPE_S_COS2::Union{Float64,Missing} = missing
    SHAPE_COS3::Union{Float64,Missing} = missing
    SHAPE_S_COS3::Union{Float64,Missing} = missing
    SHAPE_SIN3::Union{Float64,Missing} = missing
    SHAPE_S_SIN3::Union{Float64,Missing} = missing


    #species-specific parameters 
    ANISO_MODEL_1::Union{Int,Missing} = missing
    DENS_1::Union{Float64,Missing} = missing
    DLNNDR_1::Union{Float64,Missing} = missing
    DLNTDR_1::Union{Float64,Missing} = missing
    DLNTDR_PARA_1::Union{Float64,Missing} = missing
    DLNTDR_PERP_1::Union{Float64,Missing} = missing
    MASS_1::Union{Float64,Missing} = missing
    PROFILE_DLNNDR_1_SCALE::Union{Float64,Missing} = missing
    PROFILE_DLNTDR_1_SCALE::Union{Float64,Missing} = missing
    TEMP_1::Union{Float64,Missing} = missing
    TEMP_PARA_1::Union{Float64,Missing} = missing
    TEMP_PERP_1::Union{Float64,Missing} = missing
    Z_1::Union{Int64,Missing} = missing

    ANISO_MODEL_2::Union{Int,Missing} = missing
    DENS_2::Union{Float64,Missing} = missing
    DLNNDR_2::Union{Float64,Missing} = missing
    DLNTDR_2::Union{Float64,Missing} = missing
    DLNTDR_PARA_2::Union{Float64,Missing} = missing
    DLNTDR_PERP_2::Union{Float64,Missing} = missing
    MASS_2::Union{Float64,Missing} = missing
    PROFILE_DLNNDR_2_SCALE::Union{Float64,Missing} = missing
    PROFILE_DLNTDR_2_SCALE::Union{Float64,Missing} = missing
    TEMP_2::Union{Float64,Missing} = missing
    TEMP_PARA_2::Union{Float64,Missing} = missing
    TEMP_PERP_2::Union{Float64,Missing} = missing
    Z_2::Union{Int64,Missing} = missing

    ANISO_MODEL_3::Union{Int,Missing} = missing
    DENS_3::Union{Float64,Missing} = missing
    DLNNDR_3::Union{Float64,Missing} = missing
    DLNTDR_3::Union{Float64,Missing} = missing
    DLNTDR_PARA_3::Union{Float64,Missing} = missing
    DLNTDR_PERP_3::Union{Float64,Missing} = missing
    MASS_3::Union{Float64,Missing} = missing
    PROFILE_DLNNDR_3_SCALE::Union{Float64,Missing} = missing
    PROFILE_DLNTDR_3_SCALE::Union{Float64,Missing} = missing
    TEMP_3::Union{Float64,Missing} = missing
    TEMP_PARA_3::Union{Float64,Missing} = missing
    TEMP_PERP_3::Union{Float64,Missing} = missing
    Z_3::Union{Int64,Missing} = missing

    ANISO_MODEL_4::Union{Int,Missing} = missing
    DENS_4::Union{Float64,Missing} = missing
    DLNNDR_4::Union{Float64,Missing} = missing
    DLNTDR_4::Union{Float64,Missing} = missing
    DLNTDR_PARA_4::Union{Float64,Missing} = missing
    DLNTDR_PERP_4::Union{Float64,Missing} = missing
    MASS_4::Union{Float64,Missing} = missing
    PROFILE_DLNNDR_4_SCALE::Union{Float64,Missing} = missing
    PROFILE_DLNTDR_4_SCALE::Union{Float64,Missing} = missing
    TEMP_4::Union{Float64,Missing} = missing
    TEMP_PARA_4::Union{Float64,Missing} = missing
    TEMP_PERP_4::Union{Float64,Missing} = missing
    Z_4::Union{Int64,Missing} = missing

    ANISO_MODEL_5::Union{Int,Missing} = missing
    DENS_5::Union{Float64,Missing} = missing
    DLNNDR_5::Union{Float64,Missing} = missing
    DLNTDR_5::Union{Float64,Missing} = missing
    DLNTDR_PARA_5::Union{Float64,Missing} = missing
    DLNTDR_PERP_5::Union{Float64,Missing} = missing
    MASS_5::Union{Float64,Missing} = missing
    PROFILE_DLNNDR_5_SCALE::Union{Float64,Missing} = missing
    PROFILE_DLNTDR_5_SCALE::Union{Float64,Missing} = missing
    TEMP_5::Union{Float64,Missing} = missing
    TEMP_PARA_5::Union{Float64,Missing} = missing
    TEMP_PERP_5::Union{Float64,Missing} = missing
    Z_5::Union{Int64,Missing} = missing

    ANISO_MODEL_6::Union{Int,Missing} = missing
    DENS_6::Union{Float64,Missing} = missing
    DLNNDR_6::Union{Float64,Missing} = missing
    DLNTDR_6::Union{Float64,Missing} = missing
    DLNTDR_PARA_6::Union{Float64,Missing} = missing
    DLNTDR_PERP_6::Union{Float64,Missing} = missing
    MASS_6::Union{Float64,Missing} = missing
    PROFILE_DLNNDR_6_SCALE::Union{Float64,Missing} = missing
    PROFILE_DLNTDR_6_SCALE::Union{Float64,Missing} = missing
    TEMP_6::Union{Float64,Missing} = missing
    TEMP_PARA_6::Union{Float64,Missing} = missing
    TEMP_PERP_6::Union{Float64,Missing} = missing
    Z_6::Union{Int64,Missing} = missing

    ANISO_MODEL_7::Union{Int,Missing} = missing
    DENS_7::Union{Float64,Missing} = missing
    DLNNDR_7::Union{Float64,Missing} = missing
    DLNTDR_7::Union{Float64,Missing} = missing
    DLNTDR_PARA_7::Union{Float64,Missing} = missing
    DLNTDR_PERP_7::Union{Float64,Missing} = missing
    MASS_7::Union{Float64,Missing} = missing
    PROFILE_DLNNDR_7_SCALE::Union{Float64,Missing} = missing
    PROFILE_DLNTDR_7_SCALE::Union{Float64,Missing} = missing
    TEMP_7::Union{Float64,Missing} = missing
    TEMP_PARA_7::Union{Float64,Missing} = missing
    TEMP_PERP_7::Union{Float64,Missing} = missing
    Z_7::Union{Int64,Missing} = missing

    ANISO_MODEL_8::Union{Int,Missing} = missing
    DENS_8::Union{Float64,Missing} = missing
    DLNNDR_8::Union{Float64,Missing} = missing
    DLNTDR_8::Union{Float64,Missing} = missing
    DLNTDR_PARA_8::Union{Float64,Missing} = missing
    DLNTDR_PERP_8::Union{Float64,Missing} = missing
    MASS_8::Union{Float64,Missing} = missing
    PROFILE_DLNNDR_8_SCALE::Union{Float64,Missing} = missing
    PROFILE_DLNTDR_8_SCALE::Union{Float64,Missing} = missing
    TEMP_8::Union{Float64,Missing} = missing
    TEMP_PARA_8::Union{Float64,Missing} = missing
    TEMP_PERP_8::Union{Float64,Missing} = missing
    Z_8::Union{Int64,Missing} = missing

    ANISO_MODEL_9::Union{Int,Missing} = missing
    DENS_9::Union{Float64,Missing} = missing
    DLNNDR_9::Union{Float64,Missing} = missing
    DLNTDR_9::Union{Float64,Missing} = missing
    DLNTDR_PARA_9::Union{Float64,Missing} = missing
    DLNTDR_PERP_9::Union{Float64,Missing} = missing
    MASS_9::Union{Float64,Missing} = missing
    PROFILE_DLNNDR_9_SCALE::Union{Float64,Missing} = missing
    PROFILE_DLNTDR_9_SCALE::Union{Float64,Missing} = missing
    TEMP_9::Union{Float64,Missing} = missing
    TEMP_PARA_9::Union{Float64,Missing} = missing
    TEMP_PERP_9::Union{Float64,Missing} = missing
    Z_9::Union{Int64,Missing} = missing

    ANISO_MODEL_10::Union{Int,Missing} = missing
    DENS_10::Union{Float64,Missing} = missing
    DLNNDR_10::Union{Float64,Missing} = missing
    DLNTDR_10::Union{Float64,Missing} = missing
    DLNTDR_PARA_10::Union{Float64,Missing} = missing
    DLNTDR_PERP_10::Union{Float64,Missing} = missing
    MASS_10::Union{Float64,Missing} = missing
    PROFILE_DLNNDR_10_SCALE::Union{Float64,Missing} = missing
    PROFILE_DLNTDR_10_SCALE::Union{Float64,Missing} = missing
    TEMP_10::Union{Float64,Missing} = missing
    TEMP_PARA_10::Union{Float64,Missing} = missing
    TEMP_PERP_10::Union{Float64,Missing} = missing
    Z_10::Union{Int64,Missing} = missing

    ANISO_MODEL_11::Union{Int,Missing} = missing
    DENS_11::Union{Float64,Missing} = missing
    DLNNDR_11::Union{Float64,Missing} = missing
    DLNTDR_11::Union{Float64,Missing} = missing
    DLNTDR_PARA_11::Union{Float64,Missing} = missing
    DLNTDR_PERP_11::Union{Float64,Missing} = missing
    MASS_11::Union{Float64,Missing} = missing
    PROFILE_DLNNDR_11_SCALE::Union{Float64,Missing} = missing
    PROFILE_DLNTDR_11_SCALE::Union{Float64,Missing} = missing
    TEMP_11::Union{Float64,Missing} = missing
    TEMP_PARA_11::Union{Float64,Missing} = missing
    TEMP_PERP_11::Union{Float64,Missing} = missing
    Z_11::Union{Int64,Missing} = missing

    #switches 
    BTCCW::Union{Int,Missing} = missing
    COLLISION_MODEL::Union{Int,Missing} = missing
    EQUILIBRIUM_MODEL::Union{Int,Missing} = missing
    IPCCW::Union{Int,Missing} = missing
    N_ENERGY::Union{Int,Missing} = missing
    N_RADIAL::Union{Int,Missing} = missing
    N_SPECIES::Union{Int,Missing} = missing
    N_THETA::Union{Int,Missing} = missing
    N_XI::Union{Int,Missing} = missing
    PROFILE_EQUILIBRIUM_MODEL::Union{Int,Missing} = missing
    PROFILE_ERAD0_MODEL::Union{Int,Missing} = missing
    PROFILE_MODEL::Union{Int,Missing} = missing
    ROTATION_MODEL::Union{Int,Missing} = missing
    SILENT_FLAG::Union{Int,Missing} = missing
    SIM_MODEL::Union{Int,Missing} = missing
    SPITZER_MODEL::Union{Int,Missing} = missing
    THREED_MODEL::Union{Int,Missing} = missing
    THREED_EXB_MODEL::Union{Int,Missing} = missing
end
"""
    InputNEO(dd::IMAS.dd, gridpoint_cp)

Populates InputNEO structure with quantities from dd using NEO normalizations
"""
function InputNEO(dd::IMAS.dd, gridpoint_cp)
    input_neo = InputNEO()

    eq = dd.equilibrium
    eqt = eq.time_slice[]
    eq1d = eqt.profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]
    ions = cp1d.ion

    e = IMAS.cgs.e # statcoul
    k = IMAS.cgs.k # erg/eV
    mp = IMAS.cgs.mp # g
    me = IMAS.cgs.me # g
    md = IMAS.cgs.md # g
    m_to_cm = IMAS.cgs.m_to_cm
    m³_to_cm³ = IMAS.cgs.m³_to_cm³

    rmin = IMAS.r_min_core_profiles(cp1d, eqt)
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

    input_neo.DELTA = IMAS.interp1d(eq1d.rho_tor_norm, 0.5 * (eq1d.triangularity_lower + eq1d.triangularity_upper)).(cp1d.grid.rho_tor_norm)[gridpoint_cp]

    kappa = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.elongation).(cp1d.grid.rho_tor_norm)
    input_neo.KAPPA = kappa[gridpoint_cp]

    loglam = IMAS.lnΛ_ei(cp1d.electrons.density[gridpoint_cp], cp1d.electrons.temperature[gridpoint_cp])
    Z1 = IMAS.avgZ(ions[1].element[1].z_n, T1)
    m1 = ions[1].element[1].a * mp
    nu1 = sqrt(2) * pi * dens_1 * Z1^4.0 * e^4.0 * loglam ./ (sqrt(m1) * (k * temp_1).^1.5)

    input_neo.NU_1 = (nu1 ./ (v_norm ./ a))[gridpoint_cp]

    w0 = -cp1d.rotation_frequency_tor_sonic[gridpoint_cp]
    w0p = -IMAS.gradient(rmin, cp1d.rotation_frequency_tor_sonic)[gridpoint_cp]

    input_neo.OMEGA_ROT = w0 / (v_norm / a)
    input_neo.OMEGA_ROT_DERIV = w0p * a^2 / v_norm

    ####################

    q_profile = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.q).(cp1d.grid.rho_tor_norm)
    q = q_profile[gridpoint_cp]
    input_neo.Q = abs(q)

    input_neo.RHO_STAR = IMAS.rho_s(cp1d, eqt)[gridpoint_cp] / a

    Rmaj = IMAS.interp1d(eq1d.rho_tor_norm, m_to_cm * 0.5 * (eq1d.r_outboard .+ eq1d.r_inboard)).(cp1d.grid.rho_tor_norm)
    input_neo.RMAJ_OVER_A = Rmaj[gridpoint_cp] / a

    delta = IMAS.interp1d(eq1d.rho_tor_norm, 0.5 * (eq1d.triangularity_lower + eq1d.triangularity_upper)).(cp1d.grid.rho_tor_norm)
    sdelta = rmin .* IMAS.gradient(rmin, delta)
    input_neo.S_DELTA = sdelta[gridpoint_cp]

    skappa = rmin .* IMAS.gradient(rmin, kappa) ./ kappa
    input_neo.S_KAPPA = skappa[gridpoint_cp]

    zeta =
        IMAS.interp1d(
            eq1d.rho_tor_norm,
            0.25 * (eq1d.squareness_lower_inner .+ eq1d.squareness_lower_outer .+ eq1d.squareness_upper_inner .+ eq1d.squareness_upper_outer)
        ).(cp1d.grid.rho_tor_norm)
    input_neo.ZETA = zeta[gridpoint_cp]
    szeta = rmin .* IMAS.gradient(rmin, zeta)
    input_neo.S_ZETA = szeta[gridpoint_cp]

    dqdr = IMAS.gradient(rmin, q_profile)[gridpoint_cp]
    input_neo.SHEAR = rmin[gridpoint_cp] / q * dqdr

    drmaj = IMAS.gradient(rmin, Rmaj)
    input_neo.SHIFT = drmaj[gridpoint_cp]

    Z0 = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.geometric_axis.z * 1e2).(cp1d.grid.rho_tor_norm)
    input_neo.ZMAG_OVER_A = Z0[gridpoint_cp] / a
    sZ0 = IMAS.gradient(rmin, Z0)
    input_neo.S_ZMAG = sZ0[gridpoint_cp]

    for iion in eachindex(ions)
        species = iion
        setfield!(input_neo, Symbol("ANISO_MODEL_$species"), 1)
        setfield!(input_neo, Symbol("MASS_$species"), ions[iion].element[1].a .* mp / md)
        setfield!(input_neo, Symbol("Z_$species"), Int(ions[iion].element[1].z_n / ions[1].element[1].z_n))

        Ti = ions[iion].temperature ./ t_norm
        dlntidr = -IMAS.calc_z(rmin ./ a, Ti, :backward)
        Ti = Ti[gridpoint_cp]
        dlntidr = dlntidr[gridpoint_cp]

        ni = ions[iion].density ./ m³_to_cm³ / n_norm
        dlnnidr = -IMAS.calc_z(rmin ./ a, ni, :backward)
        ni = ni[gridpoint_cp]
        dlnnidr = dlnnidr[gridpoint_cp]

        setfield!(input_neo, Symbol("TEMP_$species"), Ti)
        setfield!(input_neo, Symbol("DENS_$species"), ni)
        setfield!(input_neo, Symbol("DLNNDR_$species"), dlnnidr)
        setfield!(input_neo, Symbol("DLNTDR_$species"), dlntidr)
    end

    for i in range(1, 11)
        density_val = getfield(input_neo, Symbol("DENS_$i"))
        if ismissing(density_val)
            setfield!(input_neo, Symbol("DENS_$i"), ne / n_norm)
            setfield!(input_neo, Symbol("TEMP_$i"), Te / t_norm)
            setfield!(input_neo, Symbol("ANISO_MODEL_$i"), 1)
            setfield!(input_neo, Symbol("MASS_$i"), me / md)
            setfield!(input_neo, Symbol("Z_$i"), -1)
            setfield!(input_neo, Symbol("DLNNDR_$i"), dlnnedr)
            setfield!(input_neo, Symbol("DLNTDR_$i"), dlntedr)
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

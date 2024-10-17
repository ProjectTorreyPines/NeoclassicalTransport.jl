using SpecialFunctions

Base.@kwdef mutable struct parameter_matrices
    Z::Union{Vector{Float64},Missing} = missing
    mass::Union{Vector{Float64},Missing} = missing
    dens::Union{Matrix{Float64},Missing} = missing
    temp::Union{Matrix{Float64},Missing} = missing
    nu::Union{Matrix{Float64},Missing} = missing
    dlnndr::Union{Matrix{Float64},Missing} = missing
    dlntdr::Union{Matrix{Float64},Missing} = missing
    vth::Union{Matrix{Float64},Missing} = missing
end

"""
    get_equilibrium_parameters(dd::IMAS.dd)

Populates equilibrium_geometry structure with equilibrium quantities from dd
"""
function get_equilibrium_parameters(dd::IMAS.dd)
    equilibrium_geometry = NEO.equilibrium_geometry()

    eqt = dd.equilibrium.time_slice[]
    eq1d = eqt.profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]

    m_to_cm = IMAS.gacode_units.m_to_cm

    rmin = IMAS.r_min_core_profiles(cp1d, eqt)
    a = rmin[end]
    rmaj = IMAS.interp1d(eq1d.rho_tor_norm, m_to_cm * 0.5 * (eq1d.r_outboard .+ eq1d.r_inboard)).(cp1d.grid.rho_tor_norm) ./ a
    q = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.q).(cp1d.grid.rho_tor_norm)

    ftrap = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.trapped_fraction).(cp1d.grid.rho_tor_norm)

    rmin_eqt = 0.5 * (eqt.profiles_1d.r_outboard - eqt.profiles_1d.r_inboard)
    bunit_eqt = IMAS.gradient(2 * pi * rmin_eqt, eqt.profiles_1d.phi) ./ rmin_eqt

    Bmag2_avg_eq = eqt.profiles_1d.gm5 ./ bunit_eqt .^ 2
    Bmag2_avg = IMAS.interp1d(eq1d.rho_tor_norm, Bmag2_avg_eq).(cp1d.grid.rho_tor_norm)

    f_cp = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.f).(cp1d.grid.rho_tor_norm)
    bunit_cp = IMAS.interp1d(eq1d.rho_tor_norm, IMAS.bunit(eqt)).(cp1d.grid.rho_tor_norm)
    f = f_cp .* m_to_cm ./ bunit_cp

    equilibrium_geometry.rmin = rmin
    equilibrium_geometry.rmaj = rmaj
    equilibrium_geometry.a = a
    equilibrium_geometry.q = q
    equilibrium_geometry.ftrap = ftrap
    equilibrium_geometry.Bmag2_avg = Bmag2_avg
    equilibrium_geometry.f = f

    return equilibrium_geometry

end

"""
    get_ion_electron_parameters(dd::IMAS.dd)

Populates parameter_matrices structure with profile data from dd using NEO normalizations
"""
function get_ion_electron_parameters(dd::IMAS.dd)
    parameter_matrices = NEO.parameter_matrices()

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    n = length(cp1d.grid.rho_tor_norm)

    rmin = IMAS.r_min_core_profiles(cp1d, eqt)
    a = rmin[end]

    num_ions = length(cp1d.ion)

    e = IMAS.gacode_units.e # statcoul
    k = IMAS.gacode_units.k # erg/eV
    mp = IMAS.gacode_units.mp # g
    me = IMAS.gacode_units.me # g
    md = IMAS.gacode_units.md # g
    m³_to_cm³ = IMAS.gacode_units.m³_to_cm³

    Te = cp1d.electrons.temperature # ev
    ne = cp1d.electrons.density / m³_to_cm³ # cm^-3

    loglam = 24.0 .- log.(sqrt.(ne) ./ Te)

    n_norm = ne
    m_norm = md
    t_norm = Te
    nu_norm = sqrt.(k .* Te ./ md) ./ a

    Z = Vector{Float64}(undef, num_ions + 1)
    mass = Vector{Float64}(undef, num_ions + 1)
    dens = zeros(Float64, n, num_ions + 1)
    temp = zeros(Float64, n, num_ions + 1)
    vth = zeros(Float64, n, num_ions + 1)
    nu = zeros(Float64, n, num_ions + 1)
    dlnndr = zeros(Float64, n, num_ions + 1)
    dlntdr = zeros(Float64, n, num_ions + 1)
    for i in 1:num_ions
        T1 = cp1d.ion[i].temperature[Int(ceil(end / 2))]
        Z[i] = IMAS.avgZ(cp1d.ion[i].element[1].z_n, T1)
        mass[i] = cp1d.ion[i].element[1].a * mp / m_norm

        dens[:, i] = cp1d.ion[i].density ./ m³_to_cm³ ./ n_norm
        temp[:, i] = cp1d.ion[i].temperature ./ t_norm

        nu[:, i] =
            (@. sqrt(2) * pi * (cp1d.ion[i].density ./ m³_to_cm³) * Z[i]^4.0 * e^4.0 * loglam / sqrt(cp1d.ion[i].element[1].a * mp) / (k * cp1d.ion[i].temperature)^1.5) ./ nu_norm

        dlnndr[:, i] = -IMAS.calc_z(rmin / a, cp1d.ion[i].density, :backward)
        dlntdr[:, i] = -IMAS.calc_z(rmin / a, cp1d.ion[i].temperature, :backward)

        vth[:, i] = sqrt.((cp1d.ion[i].temperature ./ t_norm) ./ (cp1d.ion[i].element[1].a * mp / m_norm))
    end

    # tacking on electron parameters at the end 
    Z[end] = -1.0
    mass[end] = me / md
    dens[:, end] = ne ./ n_norm
    temp[:, end] = Te ./ t_norm
    nu[:, end] = (@. sqrt(2) * pi * ne * (Z[end] * e)^4.0 * loglam / sqrt(me) / (k .* Te)^1.5) ./ nu_norm
    dlnndr[:, end] = -IMAS.calc_z(rmin / a, ne, :backward)
    dlntdr[:, end] = -IMAS.calc_z(rmin / a, Te, :backward)
    vth[:, end] = sqrt.((Te ./ t_norm) ./ (me ./ md))

    parameter_matrices.Z = Z
    parameter_matrices.mass = mass
    parameter_matrices.dens = dens
    parameter_matrices.temp = temp
    parameter_matrices.nu = nu
    parameter_matrices.dlnndr = dlnndr
    parameter_matrices.dlntdr = dlntdr
    parameter_matrices.vth = vth

    return parameter_matrices
end

function gauss_legendre(x1::Int, x2::Int, n::Int)
    eps = 1e-15

    xm = 0.5 * (x2 + x1)
    xl = 0.5 * (x2 - x1)

    x = zeros(Float64, n)
    w = zeros(Float64, n)

    # Exception for n=1 is required:
    if n == 1
        x[1] = xm
        w[1] = 2.0 * xl
        return x, w
    end

    # Roots are symmetric. We only need to find half of them.
    m = (n + 1) / 2

    # Initialize to fail first do test
    z1 = -1.0
    pp = 0.0

    # Loop over roots.
    for i in 1:m
        i = Int(i)
        z = cos(pi * (i - 0.25) / (n + 0.5))

        while abs(z - z1) > eps
            p1 = 1.0
            p2 = 0.0

            for j in 1:n
                p3 = p2
                p2 = p1
                p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j
            end

            # p1 is the Legendre polynomial. Now compute its derivative, pp.
            pp = n * (z * p1 - p2) / (z * z - 1.0)
            z1 = z
            z = z1 - p1 / pp
        end

        x[i] = xm - xl * z
        x[n+1-i] = xm + xl * z
        w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp)
        w[n+1-i] = w[i]
    end

    return x, w
end

function gauss_integ(
    xmin::Float64,
    xmax::Float64,
    func::Function,
    order::Int,
    n_subdiv::Int,
    parameter_matrices::NEO.parameter_matrices,
    ietype::Int,
    equilibrium_geometry::NEO.equilibrium_geometry,
    is_globalHS::Int,
    ir_global::Int)

    x0, w0 = gauss_legendre(0, 1, order)

    dx = (xmax - xmin) / n_subdiv

    n_node = n_subdiv * order

    x = zeros(Float64, n_subdiv * order)
    w = zeros(Float64, n_subdiv * order)

    for i in 1:n_subdiv
        for j in 1:order
            p = (i - 1) * order + j
            x[p] = xmin + ((i - 1) + x0[j]) * dx
            w[p] = w0[j] * dx
        end
    end

    answer = 0.0
    for p in 1:n_node
        answer = answer + w[p] * func(x[p], parameter_matrices, ietype, equilibrium_geometry, is_globalHS, ir_global)
    end

    return answer
end

function get_coll_freqs(ir_loc::Int, is_loc::Int, js_loc::Int, ene::Float64, parameter_matrices::NEO.parameter_matrices)
    Z = parameter_matrices.Z
    dens = parameter_matrices.dens
    vth = parameter_matrices.vth

    fac = @views (1.0 * Z[js_loc])^2 / (1.0 * Z[is_loc])^2 * (dens[ir_loc, js_loc] / dens[ir_loc, is_loc])

    xa = sqrt(ene)
    xb = @views xa * (vth[ir_loc, is_loc] / vth[ir_loc, js_loc])

    if xb < 1e-4
        nu_d =
            fac * (1.0 / sqrt(pi)) *
            (
                4.0 / 3.0 * (vth[ir_loc, is_loc] / vth[ir_loc, js_loc]) - 4.0 / 15.0 * (vth[ir_loc, is_loc] / vth[ir_loc, js_loc])^3 * ene +
                2.0 / 35.0 * (vth[ir_loc, is_loc] / vth[ir_loc, js_loc])^5 * ene^2 -
                2.0 / 189.0 * (vth[ir_loc, is_loc] / vth[ir_loc, js_loc])^7 * ene^3
            )
    else
        Hd_coll = exp(-xb * xb) / (xb * sqrt(pi)) + (1 - (1 / (2 * xb * xb))) * erf(xb)
        Xd_coll = 1 / xa
        nu_d = fac * Hd_coll * Xd_coll
    end

    return nu_d
end

function myHSenefunc(x::Float64, parameter_matrices::NEO.parameter_matrices, ietype::Int, equilibrium_geometry::NEO.equilibrium_geometry, is_globalHS::Int, ir_global::Int)
    rmin = equilibrium_geometry.rmin
    rmaj = equilibrium_geometry.rmaj
    a = equilibrium_geometry.a
    q = equilibrium_geometry.q
    ftrap = equilibrium_geometry.ftrap

    nu = parameter_matrices.nu
    vth = parameter_matrices.vth

    emin = 0.0
    emax = 16.0

    xa = 2.0 / (1.0 - sqrt(emin / emax))
    xb = -(1.0 + sqrt(emin / emax)) / (1.0 - sqrt(emin / emax))
    ene = emax * ((x - xb) / xa)^2
    de = 2.0 * sqrt(emax) / xa
    val = de * exp(-ene)

    eps = rmin[ir_global] / (rmaj[ir_global] .* a)

    nu_d_tot = 0.0
    for js in eachindex(parameter_matrices.Z)
        nu_d = get_coll_freqs(ir_global, is_globalHS, js, ene, parameter_matrices)
        nu_d_tot += nu_d * nu[ir_global, is_globalHS]
    end

    ft_star = (3.0 * pi / 16.0) * eps^2 * vth[ir_global, is_globalHS] * sqrt(2.0) * ene^1.5 / (rmaj[ir_global] * abs(q[ir_global]) * nu_d_tot)
    ft_fac = 1.0 / (1.0 + ftrap[ir_global] / ft_star)

    if ietype == 1
        out = val * nu_d_tot * ene * ft_fac
    elseif ietype == 2
        out = val * nu_d_tot * ene * ene * ft_fac
    elseif ietype == 3
        out = val * nu_d_tot * ene * ene * ene * ft_fac
    end

    return out
end

function compute_HS(ir::Int, dd::IMAS.dd, parameter_matrices::NEO.parameter_matrices, equilibrium_geometry::NEO.equilibrium_geometry)
    rmin = equilibrium_geometry.rmin
    a = equilibrium_geometry.a
    q = equilibrium_geometry.q
    ftrap = equilibrium_geometry.ftrap
    Bmag2_avg = equilibrium_geometry.Bmag2_avg
    f = equilibrium_geometry.f

    Z = parameter_matrices.Z
    mass = parameter_matrices.mass
    dens = parameter_matrices.dens
    temp = parameter_matrices.temp
    dlnndr = parameter_matrices.dlnndr
    dlntdr = parameter_matrices.dlntdr

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    n_species = length(cp1d.ion) + 1

    rho = IMAS.rho_s(cp1d, eqt) ./ a

    Nx = 10 # Can be lowered to speed up calculation time
    integ_order = 1
    omega_fac = 1.0 / Bmag2_avg[ir]
    HS_I_div_psip = -f[ir] * q[ir] / rmin[ir]

    nux0 = zeros(Float64, n_species)
    nux2 = zeros(Float64, n_species)
    nux4 = zeros(Float64, n_species)

    cp1d = dd.core_profiles.profiles_1d[]

    for is_global in 1:n_species
        for ietype in 1:3
            ir_global = ir
            is_globalHS = is_global

            eii_val = gauss_integ(-1.0, 1.0, NEO.myHSenefunc, integ_order, Nx, parameter_matrices, ietype, equilibrium_geometry, is_globalHS, ir_global)

            if ietype == 1
                nux0[is_global] = eii_val * 4.0 / (3.0 * sqrt(pi))
            elseif ietype == 2
                nux2[is_global] = eii_val * 4.0 / (3.0 * sqrt(pi))
            elseif ietype == 3
                nux4[is_global] = eii_val * 4.0 / (3.0 * sqrt(pi))
            end
        end
    end

    sum_nm = 0.0
    for is_global in 1:n_species
        sum_nm += mass[is_global] * dens[ir, is_global] * nux0[is_global]
    end

    pflux_multi = zeros(Float64, n_species)
    eflux_multi = zeros(Float64, n_species)
    for is_global in 1:n_species
        A1 = -dlnndr[ir, is_global] + (1.5 * dlntdr[ir, is_global])
        A2 = -dlntdr[ir, is_global]

        pflux_multi[is_global] = 0.0
        eflux_multi[is_global] = 0.0

        L_a =
            nux0[is_global] * omega_fac * HS_I_div_psip^2 * rho[ir]^2 * ftrap[ir] * dens[ir, is_global] * temp[ir, is_global] * mass[is_global] /
            (Z[is_global] * Z[is_global] * 1.0)

        for js in 1:n_species

            L_b = nux0[js] * omega_fac * HS_I_div_psip^2 * rho[ir]^2 * ftrap[ir] * dens[ir, js] * temp[ir, js] * mass[js] / (Z[js] * Z[js] * 1.0)

            if is_global == js
                L11 = -L_a * (sum_nm - mass[is_global] * dens[ir, is_global] * nux0[is_global]) / sum_nm
                L12 = L11 * (nux2[is_global] / nux0[is_global])
                L22 = (nux2[is_global] / nux0[is_global]) * (L12 + L_a * (nux2[is_global] / nux0[is_global] - nux4[is_global] / nux2[is_global]))
                L21 = L12

            else
                L11 = L_a * (Z[is_global] * temp[ir, js]) / (Z[js] * temp[ir, is_global]) * mass[js] * dens[ir, js] * nux0[js] / sum_nm
                L12 = (nux2[js] / nux0[js]) * L11
                L21 = (nux2[is_global] / nux0[is_global]) * Z[js] / (1.0 * Z[is_global]) * (mass[is_global] * dens[ir, is_global] * nux0[is_global] / sum_nm) * L_b
                L22 = (nux2[is_global] / nux0[is_global]) * L12

            end

            pflux_multi[is_global] += L11 * A1 + L12 * A2
            eflux_multi[is_global] += L21 * A1 + L22 * A2

        end
    end

    return pflux_multi, eflux_multi

end

function HS_to_GB(HS_solution::Tuple{Vector{Float64},Vector{Float64}}, dd::IMAS.dd, rho::Int)
    pflux_multi, eflux_multi = HS_solution

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    rmin = IMAS.r_min_core_profiles(cp1d, eqt)
    a = rmin[end]

    neo_rho_star = (IMAS.rho_s(cp1d, eqt)./a)[rho]

    temp_e = 1.0 # electron temperature is 1 since all NEO temps are normalized against electron temp
    dens_e = 1.0 # electron density is 1 since all NEO densities are normalized against electron density

    Gamma_neo_GB = dens_e * temp_e^1.5 * neo_rho_star .^ 2
    Q_neo_GB = dens_e * temp_e^2.5 * neo_rho_star .^ 2

    pflux_norm = pflux_multi ./ Gamma_neo_GB
    eflux_norm = eflux_multi ./ Q_neo_GB

    energy_flux_e = eflux_norm[end]
    energy_flux_i = sum(eflux_norm[1:end-1])
    particle_flux_e = pflux_norm[end]
    particle_flux_i = pflux_norm[1:end-1]

    # assign fluxes to flux_solution structure
    sol = IMAS.flux_solution(energy_flux_e, energy_flux_i, particle_flux_e, particle_flux_i, 0.0)
    return sol
end

"""
    hirshmansigmar(ir::Int, dd::IMAS.dd, parameter_matrices::NEO.parameter_matrices, equilibrium_geometry::NEO.equilibrium_geometry)

Calculates neoclassical fluxes according to Hirshman-Sigmar model
Ref: S.P. Hirshman, D.J. Sigmar, J.F. Clarke, Phys. Fluids 19, 656–666 (1976), https://doi.org/10.1063/1.861524
"""
function hirshmansigmar(ir::Int, dd::IMAS.dd, parameter_matrices::NEO.parameter_matrices, equilibrium_geometry::NEO.equilibrium_geometry)
    hs = compute_HS(ir, dd, parameter_matrices, equilibrium_geometry)
    return HS_to_GB(hs, dd, ir)
end
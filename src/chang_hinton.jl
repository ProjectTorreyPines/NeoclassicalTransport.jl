"""
    changhinton(
        eqt1d::IMAS.equilibrium__time_slice___profiles_1d, 
        cp1d::IMAS.core_profiles__profiles_1d,
        rho_fluxmatch::Real,
        iion::Integer)

Calculates the neoclassical flux using Chang-Hinton model which has been modified assuming Zi = 1, and ni=ne

Ref: C.S. Chang, F.L. Hinton, Phys. Fluids 25, 1493–1494 (1982), https://doi.org/10.1063/1.863934
"""
function changhinton(
    eqt::IMAS.equilibrium__time_slice,
    cp1d::IMAS.core_profiles__profiles_1d,
    rho_fluxmatch::Real,
    iion::Integer)

    eqt1d = eqt.profiles_1d
    ions = cp1d.ion

    e = IMAS.cgs.e # statcoul
    k = IMAS.cgs.k # erg/eV
    mp = IMAS.cgs.mp # g
    m_to_cm = IMAS.cgs.m_to_cm
    m³_to_cm³ = IMAS.cgs.m³_to_cm³

    Rmaj0 = eqt1d.geometric_axis.r[1] * m_to_cm

    rho_eq = eqt1d.rho_tor_norm
    rho_cp = cp1d.grid.rho_tor_norm
    gridpoint_cp = argmin(abs.(rho_cp .- rho_fluxmatch))

    ne = cp1d.electrons.density_thermal / m³_to_cm³
    Te = cp1d.electrons.temperature[gridpoint_cp]
    Ti = ions[iion].temperature

    Rmaj = IMAS.interp1d(eqt1d.rho_tor_norm, 0.5 * (eqt1d.r_outboard .+ eqt1d.r_inboard)).(cp1d.grid.rho_tor_norm)
    rmin = IMAS.interp1d(rho_eq, 0.5 * (eqt1d.r_outboard - eqt1d.r_inboard)).(rho_cp)
    drmaj = IMAS.gradient(rmin, Rmaj)
    shift = -drmaj[gridpoint_cp]

    eps = rmin[gridpoint_cp] * m_to_cm / Rmaj0
    q = IMAS.interp1d(rho_eq, eqt1d.q).(rho_cp)[gridpoint_cp]

    a = eqt.boundary.minor_radius * m_to_cm
    dlntdr = -IMAS.gradient(rmin, Ti)[gridpoint_cp] / Ti[gridpoint_cp]
    dlntdr = dlntdr * a / m_to_cm
    dlnndr = -IMAS.gradient(rmin, ne)[gridpoint_cp] / ne[gridpoint_cp]
    dlnndr = dlnndr * a / m_to_cm
    Ti = Ti[gridpoint_cp]
    ne = ne[gridpoint_cp]

    mi = ions[iion].element[1].a * mp

    c_s = sqrt(k * Te / mi)
    loglam = 24.0 - log(sqrt(ne / Te))

    k0 = 0.66
    a0 = 1.03
    b0 = 0.31
    c0 = 0.74

    Zi = IMAS.avgZ(ions[iion].element[1].z_n, Ti)
    alpha = Zi - 1.0
    nui = sqrt(2.0) * pi * ne * (Zi * e)^4 * loglam / sqrt(mi) / (k * Ti)^1.5
    nu = nui * a / c_s / sqrt(Ti / Ti) * (Ti / Te)^1.5
    nui_HH = nu * (4.0 / 3.0) / sqrt(2.0 * pi)
    nui_star_HH = nui_HH * (Rmaj0 / a) * abs(q / (sqrt(eps) * sqrt(eps) * sqrt(eps)))
    mu_star = (1.0 + 1.54 * alpha) * nui_star_HH

    CH_Bmag2inv_avg = ((1.0 + 1.5 * (eps * eps + eps * shift)
                        + 0.375 * eps * eps * eps * shift)
                       /
                       (1.0 + 0.5 * eps * shift))
    CH_Bmag2avg_inv = ((sqrt(1.0 - eps * eps) * (1.0 + 0.5 * eps * shift))
                       /
                       (1.0 + (shift / eps) * (sqrt(1.0 - eps * eps) - 1.0)))
    CH_I_div_psip = q / eps

    F2 = (0.5 / sqrt(eps)) * (CH_Bmag2inv_avg - CH_Bmag2avg_inv)

    K1 = (-alpha * (0.83 + 0.42 * alpha) / (0.58 + alpha)
          *
          (c0 * mu_star * sqrt(eps * eps * eps) * F2)
          /
          (1.0 + c0 * mu_star * sqrt(eps * eps * eps)))

    K2 = (
        (k0 * (1.0 + 1.54 * alpha)
         +
         (1.88 * sqrt(eps) - 1.54 * eps) * (1.0 + 3.75 * alpha)) * CH_Bmag2inv_avg
        /
        (1 + a0 * sqrt(mu_star) + b0 * mu_star)
        +
        k0 * sqrt(eps * eps * eps) * (c0 * c0 / b0) * mu_star * F2
        * (1.0 + 1.33 * alpha * (1.0 + 0.6 * alpha) / (1.0 + 1.79 * alpha))
        /
        (1.0 + c0 * sqrt(eps * eps * eps) * mu_star)
    )

    neo_rho_star_in = 0.001

    efluxi = (CH_I_div_psip^2
              *
              Ti / Te
              * (ions[iion].element[1].a / Zi^2 * neo_rho_star_in^2 * sqrt(eps) * nui_HH)
              * ((K2 + K1) * dlntdr + K1 * dlnndr))

    qneo_gb = neo_rho_star_in^2

    # assign fluxes to FluxSolution structure
    sol = GACODE.FluxSolution(0.0, efluxi / qneo_gb, 0.0, Float64[], 0.0)
    return sol
end
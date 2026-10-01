# Transport moments of the solution (TRANSP_do and compute_velocity of
# neo_transport.f90): particle, energy and momentum fluxes, parallel flows,
# bootstrap current, poloidal/toroidal velocities at theta=0, the Sugama
# gyroviscous "H" fluxes, and NEO's check_sum.

"""
    NEOSolution{T}

Everything NEO writes for one flux surface, per species and unlumped, in
NEO's normalizations (`neo_transport.f90` header):

- `pflux`, `eflux`, `mflux`: dke fluxes Γ/(n0 vt0), Q/(n0 vt0 T0), Π/(n0 a T0)
  (`eflux` includes `+omega_rot*mflux` like NEO)
- `pflux_gv`, `eflux_gv`, `mflux_gv`: gyroviscous fluxes
- `uparB` (<u_par B>/(vt0 B0)), `uparBN`, `klittle_upar`, `kbig_upar`
- `vpol_th0`, `vtor_th0` (first-order velocities at theta=0), `vtor_0order_th0`, `uparB_0order`
- `jpar` (<j_par B>/(e n0 vt0 B_unit)), `jtor`
- `d_phi(theta)`, `d_phi_sqavg` (first-order potential, rotation_model=1 only)
- `upar`, `vpol`, `vtor`: (n_species, n_theta) velocity fields
- `check_sum`: NEO's regression checksum (`out.neo.prec`)
- `g`: the solution vector in `mindx` order (`nothing` unless `keep_g=true`)
- `params`: the [`NEOParams`](@ref) solved

[`tgyro_fluxes`](@ref) gives the GB-normalised totals of `out.neo.transport_flux`,
`GACODE.FluxSolution(sol)` the lumped form FUSE uses.
"""
struct NEOSolution{T<:Real}
    params::NEOParams{T}
    pflux::Vector{T}
    eflux::Vector{T}
    mflux::Vector{T}
    pflux_gv::Vector{T}
    eflux_gv::Vector{T}
    mflux_gv::Vector{T}
    uparB::Vector{T}
    uparBN::Vector{T}
    klittle_upar::Vector{T}
    kbig_upar::Vector{T}
    vpol_th0::Vector{T}
    vtor_th0::Vector{T}
    vtor_0order_th0::T
    uparB_0order::T
    jpar::T
    jtor::T
    d_phi::Vector{T}
    d_phi_sqavg::T
    upar::Matrix{T}
    vpol::Matrix{T}
    vtor::Matrix{T}
    check_sum::T
    g::Union{Nothing,Vector{T}}
end

"""
    transport(p, basis, geo, rot, coef, g, pattern; keep_g=false) -> NEOSolution

`TRANSP_do` + `compute_velocity` on the solution vector `g`.
"""
function transport(p::NEOParams{T}, basis::NEOBasis, geo::NEOGeometry, rot::NEORotation, coef::NEOCoefficients,
    g::AbstractVector, pattern::NEOPattern; keep_g::Bool=false) where {T<:Real}
    ns, ne, nxi, nth = p.n_species, p.n_energy, p.n_xi, p.n_theta
    Z, dens, temp, vth, mass = p.z, p.dens, p.temp, p.vth, p.mass
    omega_rot, omega_rot_deriv, rho = p.omega_rot, p.omega_rot_deriv, p.rho
    w_theta, Bmag, Btor, bigR = geo.w_theta, geo.Bmag, geo.Btor, geo.bigR
    dens_fac, phi_rot, phi_rot_avg = rot.dens_fac, rot.phi_rot, rot.phi_rot_avg
    evec_e0, evec_e1, evec_e2, evec_e05, evec_e105 = basis.evec_e0, basis.evec_e1, basis.evec_e2, basis.evec_e05, basis.evec_e105
    driftx, driftxrot1, driftxrot2, driftxrot3 = coef.driftx, coef.driftxrot1, coef.driftxrot2, coef.driftxrot3
    I_div_psip = geo.I_div_psip

    pflux = zeros(T, ns)
    eflux = zeros(T, ns)
    mflux = zeros(T, ns)
    upar = zeros(T, ns, nth)
    uparB = zeros(T, ns)
    uparBN = zeros(T, ns)
    d_phi = zeros(T, nth)

    for i in 1:pattern.n_row
        is, ie, ix, it = row_indices(i, pattern)
        gi = g[i]
        nd = dens[is] * dens_fac[is, it]

        # lambda_a, the effective potential energy
        rfac = Z[is] / temp[is] * (phi_rot[it] - phi_rot_avg) - (omega_rot * bigR[it] / vth[is])^2 * 0.5

        if ix == 0
            pflux[is] = pflux[is] + w_theta[it] * 4.0 / sqrt(pi) * nd * gi *
                                    (driftx[is, it] * (4.0 / 3.0) * evec_e1[ie, ix] + driftxrot1[is, it] * evec_e0[ie, ix])
            eflux[is] = eflux[is] + w_theta[it] * temp[is] * 4.0 / sqrt(pi) * nd * gi *
                                    (driftx[is, it] * (4.0 / 3.0) * (evec_e2[ie, ix] + rfac * evec_e1[ie, ix]) +
                                     driftxrot1[is, it] * (evec_e1[ie, ix] + rfac * evec_e0[ie, ix]))
            mflux[is] = mflux[is] + w_theta[it] * temp[is] * nd * 4.0 / sqrt(pi) * gi * bigR[it] / vth[is] *
                                    (1.0 / 3.0 * driftxrot2[is, it] * sqrt(2.0) * Btor[it] / Bmag[it] * evec_e1[ie, ix] +
                                     4.0 / 3.0 * driftx[is, it] * omega_rot * bigR[it] / vth[is] * evec_e1[ie, ix] +
                                     driftxrot1[is, it] * omega_rot * bigR[it] / vth[is] * evec_e0[ie, ix])
            d_phi[it] = d_phi[it] + Z[is] * nd * gi * 4.0 / sqrt(pi) * evec_e0[ie, ix]
        elseif ix == 1
            pflux[is] = pflux[is] + w_theta[it] * nd * 4.0 / sqrt(pi) * gi * driftxrot2[is, it] * (1.0 / 3.0) * evec_e05[ie, ix]
            eflux[is] = eflux[is] + w_theta[it] * temp[is] * nd * 4.0 / sqrt(pi) * gi * driftxrot2[is, it] * (1.0 / 3.0) *
                                    (evec_e105[ie, ix] + rfac * evec_e05[ie, ix])
            mflux[is] = mflux[is] + w_theta[it] * 4.0 / sqrt(pi) * nd * gi * temp[is] / vth[is] * bigR[it] *
                                    (8.0 / 15.0 * driftx[is, it] * Btor[it] / Bmag[it] * sqrt(2.0) * evec_e105[ie, ix] +
                                     1.0 / 3.0 * driftxrot1[is, it] * sqrt(2.0) * Btor[it] / Bmag[it] * evec_e05[ie, ix] +
                                     1.0 / 3.0 * driftxrot2[is, it] * omega_rot * bigR[it] / vth[is] * evec_e05[ie, ix] +
                                     2.0 / 15.0 * driftxrot3[is, it] * evec_e105[ie, ix])
            # uparB = < B * 1/n * int vpar * (F0 g)>
            uparB[is] = uparB[is] + w_theta[it] * Bmag[it] * sqrt(2.0) * vth[is] * (1.0 / 3.0) * gi * 4.0 / sqrt(pi) * evec_e05[ie, ix]
            uparBN[is] = uparBN[is] + w_theta[it] * Bmag[it] * sqrt(2.0) * vth[is] * (1.0 / 3.0) * gi * 4.0 / sqrt(pi) * nd * evec_e05[ie, ix]
            upar[is, it] = upar[is, it] + sqrt(2.0) * vth[is] * (1.0 / 3.0) * gi * 4.0 / sqrt(pi) * evec_e05[ie, ix]
        elseif ix == 2
            pflux[is] = pflux[is] + w_theta[it] * 4.0 / sqrt(pi) * nd * gi * driftx[is, it] * (2.0 / 15.0) * evec_e1[ie, ix]
            eflux[is] = eflux[is] + w_theta[it] * temp[is] * 4.0 / sqrt(pi) * nd * gi * driftx[is, it] * (2.0 / 15.0) *
                                    (evec_e2[ie, ix] + rfac * evec_e1[ie, ix])
            mflux[is] = mflux[is] + w_theta[it] * temp[is] * nd * 4.0 / sqrt(pi) * gi * bigR[it] / vth[is] *
                                    (2.0 / 15.0 * driftxrot2[is, it] * sqrt(2.0) * Btor[it] / Bmag[it] * evec_e1[ie, ix] +
                                     2.0 / 15.0 * driftx[is, it] * omega_rot * bigR[it] / vth[is] * evec_e1[ie, ix])
        elseif ix == 3
            mflux[is] = mflux[is] + w_theta[it] * 4.0 / sqrt(pi) * nd * gi * temp[is] / vth[is] * bigR[it] * 2.0 / 35.0 * evec_e105[ie, ix] *
                                    (driftx[is, it] * Btor[it] / Bmag[it] * sqrt(2.0) - driftxrot3[is, it])
        end
    end

    for is in 1:ns
        eflux[is] = eflux[is] + omega_rot * mflux[is]
    end

    # d_phi: sum_s Z_s e int f_s = 0 (first-order potential, no rotation only)
    if p.rotation_model == 2
        fill!(d_phi, zero(T))
        d_phi_sqavg = zero(T)
    else
        poisson_F0fac = zero(T)
        for is in 1:ns
            poisson_F0fac = poisson_F0fac + dens[is] * Z[is] * Z[is] / temp[is]
        end
        if p.ae_flag == 1
            poisson_F0fac = poisson_F0fac + p.dens_ae / p.temp_ae
        end
        d_phi_sqavg = zero(T)
        for it in 1:nth
            d_phi[it] = d_phi[it] / poisson_F0fac
            d_phi_sqavg = d_phi_sqavg + d_phi[it] * d_phi[it] * w_theta[it]
        end
    end

    # bootstrap current = sum <Z*n*upar B>
    jpar = zero(T)
    for is in 1:ns
        jpar = jpar + Z[is] * uparBN[is]
    end

    # U_parallel coefficients
    kbig_upar = zeros(T, ns)
    klittle_upar = zeros(T, ns)
    for is in 1:ns
        B2_div_dens = zero(T)
        bigR2_avg = zero(T)
        for it in 1:nth
            B2_div_dens = B2_div_dens + w_theta[it] * Bmag[it]^2 / (dens[is] * dens_fac[is, it])
            bigR2_avg = bigR2_avg + w_theta[it] * bigR[it]^2
        end
        kbig_upar[is] = (1.0 / B2_div_dens) *
                        (uparB[is] - (I_div_psip * rho * temp[is] / (Z[is] * 1.0)) *
                                     (p.dlnndr[is] - (Z[is] * 1.0) / temp[is] * p.dphi0dr +
                                      p.dlntdr[is] * (1.0 + Z[is] / temp[is] * phi_rot_avg + omega_rot^2 * 0.5 / vth[is]^2 * (geo.bigR_th0^2 - bigR2_avg)) +
                                      omega_rot_deriv * omega_rot / vth[is]^2 * (geo.bigR_th0^2 - bigR2_avg) +
                                      omega_rot^2 * geo.bigR_th0 / vth[is]^2 * geo.bigR_th0_rderiv))
        if abs(_val(p.dlntdr[is])) > eps()
            klittle_upar[is] = -kbig_upar[is] * B2_div_dens / (p.dlntdr[is] * I_div_psip * rho * temp[is] / (Z[is] * 1.0))
        else
            klittle_upar[is] = zero(T)
        end
    end

    # poloidal and toroidal velocities
    vpol, vtor, vpol_th0, vtor_th0, vtor_0order_th0, uparB_0order = compute_velocity(p, geo, rot, kbig_upar, upar)

    # toroidal component of the bootstrap current = sum <Z*n*utor/R>/<1/R>
    jtor = zero(T)
    for is in 1:ns, it in 1:nth
        jtor = jtor + Z[is] * dens[is] * dens_fac[is, it] * vtor[is, it] * w_theta[it] / bigR[it]
    end
    jtor = jtor / geo.bigRinv_avg

    # Sugama gyro-viscosity "H" fluxes
    pflux_gv = zeros(T, ns)
    eflux_gv = zeros(T, ns)
    mflux_gv = zeros(T, ns)
    for is in 1:ns
        fac1 = zero(T)
        fac2 = zero(T)
        for it in 1:nth
            rfac = Z[is] / temp[is] * (phi_rot[it] - phi_rot_avg) - (omega_rot * bigR[it] / vth[is])^2 * 0.5
            gfac = w_theta[it] * dens[is] * dens_fac[is, it] / Bmag[it]^3 *
                   (2.0 * geo.gradr[it] * geo.gradpar_gradr[it] - 1.0 / Bmag[it] * geo.gradr[it]^2 * geo.gradpar_Bmag[it])
            fac1 = fac1 + gfac
            fac2 = fac2 + gfac * (1.0 + rfac)
        end
        pre = -0.5 * rho^2 * mass[is] * temp[is] / (Z[is] * 1.0)^2 * I_div_psip * p.rmin / p.q
        pflux_gv[is] = pre * fac1 * omega_rot_deriv
        mflux_gv[is] = -0.5 * temp[is] * rho^2 * mass[is] * temp[is] / (Z[is] * 1.0)^2 * I_div_psip * p.rmin / p.q *
                       (fac2 * p.dlntdr[is] + fac1 * (p.dlnndr[is] - Z[is] / temp[is] * p.dphi0dr +
                                                      omega_rot * geo.bigR_th0^2 / vth[is]^2 * omega_rot_deriv +
                                                      omega_rot^2 * geo.bigR_th0 / vth[is]^2 * geo.bigR_th0_rderiv +
                                                      p.dlntdr[is] * (1.0 + 0.5 * omega_rot^2 * geo.bigR_th0^2 / vth[is]^2) +
                                                      p.dlntdr[is] * Z[is] / temp[is] * phi_rot_avg))
        eflux_gv[is] = -0.5 * temp[is] * rho^2 * mass[is] * temp[is] / (Z[is] * 1.0)^2 * I_div_psip * p.rmin / p.q * fac2 * omega_rot_deriv +
                       2.5 * pflux_gv[is] * temp[is] + omega_rot * mflux_gv[is]
    end

    check_sum = zero(T)
    for is in 1:ns
        check_sum = check_sum + (abs(pflux[is]) + abs(eflux[is]) + abs(mflux[is])) / rho^2 + abs(uparB[is]) / rho
    end

    return NEOSolution{T}(p, pflux, eflux, mflux, pflux_gv, eflux_gv, mflux_gv, uparB, uparBN, klittle_upar, kbig_upar,
        vpol_th0, vtor_th0, vtor_0order_th0, uparB_0order, jpar, jtor, d_phi, d_phi_sqavg, upar, vpol, vtor, check_sum,
        keep_g ? collect(T, g) : nothing)
end

# compute_velocity: vpol, vtor on the grid; theta=0 values from the truncated cosine series
function compute_velocity(p::NEOParams{T}, geo::NEOGeometry, rot::NEORotation, kbig_upar::Vector{T}, upar::Matrix{T}) where {T<:Real}
    ns, nth = p.n_species, p.n_theta
    Z, dens, temp, vth = p.z, p.dens, p.temp, p.vth
    omega_rot, omega_rot_deriv, rho = p.omega_rot, p.omega_rot_deriv, p.rho
    w_theta, Bpol, Btor, bigR, theta = geo.w_theta, geo.Bpol, geo.Btor, geo.bigR, geo.theta
    dens_fac, phi_rot = rot.dens_fac, rot.phi_rot
    I_div_psip, bigR_th0, bigR_th0_rderiv = geo.I_div_psip, geo.bigR_th0, geo.bigR_th0_rderiv

    # 0th-order toroidal flow (at theta=0) and <u_par B>
    vtor_0order_th0 = omega_rot * bigR_th0
    RBt = zero(T)
    for it in 1:nth
        RBt = RBt + w_theta[it] * bigR[it] * Btor[it]
    end
    uparB_0order = omega_rot * RBt

    vpol = zeros(T, ns, nth)
    vtor = zeros(T, ns, nth)
    vpol_th0 = zeros(T, ns)
    vtor_th0 = zeros(T, ns)
    m_theta = (nth - 1) ÷ 2 - 1

    for is in 1:ns
        for it in 1:nth
            vpol[is, it] = kbig_upar[is] * Bpol[it] / (dens[is] * dens_fac[is, it])
            vtor[is, it] = kbig_upar[is] * Btor[it] / (dens[is] * dens_fac[is, it]) +
                           (I_div_psip * rho * temp[is] / (Z[is] * 1.0)) * (1.0 / Btor[it]) *
                           (p.dlnndr[is] - (1.0 * Z[is]) / temp[is] * p.dphi0dr +
                            p.dlntdr[is] * (1.0 + Z[is] / temp[is] * phi_rot[it] + omega_rot^2 * 0.5 / vth[is]^2 * (bigR_th0^2 - bigR[it]^2)) +
                            omega_rot^2 * bigR_th0 / vth[is]^2 * bigR_th0_rderiv +
                            omega_rot_deriv * omega_rot / vth[is]^2 * (bigR_th0^2 - bigR[it]^2))
        end

        for jt in 0:m_theta
            vp = zero(T)
            vt = zero(T)
            for it in 1:nth
                vp = vp + vpol[is, it] * cos(jt * theta[it])
                vt = vt + vtor[is, it] * cos(jt * theta[it])
            end
            if jt == 0
                vp = vp / (1.0 * nth)
                vt = vt / (1.0 * nth)
            else
                vp = vp / (0.5 * nth)
                vt = vt / (0.5 * nth)
            end
            # vel(theta=0) = sum of the cosine coefficients
            vpol_th0[is] = vpol_th0[is] + vp
            vtor_th0[is] = vtor_th0[is] + vt
        end
    end

    return vpol, vtor, vpol_th0, vtor_th0, vtor_0order_th0, uparB_0order
end

"""
    tgyro_fluxes(sol::NEOSolution) -> (pflux, eflux, mflux)

Per-species GB-normalised totals of `out.neo.transport_flux`'s tgyro block:
`(pflux+pflux_gv)/pgb`, `(eflux+eflux_gv-omega*(mflux+mflux_gv))/egb`,
`(mflux+mflux_gv)/mgb`, with the electron density and temperature of the
kinetic (or adiabatic) electrons setting the GB units.
"""
function tgyro_fluxes(sol::NEOSolution)
    p = sol.params
    if p.ae_flag == 0
        dens_ele = p.dens[p.is_ele]
        temp_ele = p.temp[p.is_ele]
    else
        dens_ele = p.dens_ae
        temp_ele = p.temp_ae
    end
    pgb = dens_ele * p.rho^2 * temp_ele^1.5
    egb = dens_ele * p.rho^2 * temp_ele^2.5
    mgb = dens_ele * p.rho^2 * temp_ele^2
    omega = p.omega_rot
    pflux = (sol.pflux .+ sol.pflux_gv) ./ pgb
    eflux = (sol.eflux .+ sol.eflux_gv .- omega .* sol.mflux .- omega .* sol.mflux_gv) ./ egb
    mflux = (sol.mflux .+ sol.mflux_gv) ./ mgb
    return pflux, eflux, mflux
end

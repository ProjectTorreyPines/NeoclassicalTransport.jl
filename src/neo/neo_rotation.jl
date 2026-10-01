# Poloidal variation of the equilibrium potential and densities with strong
# toroidal rotation: ROT_solve_phi of neo_rotation.f90 (isotropic species).
#
# With rotation_model=1 everything is trivial (phi_rot=0, dens_fac=1);
# the omega_rot zeroing that ROT_solve_phi does in that case is already
# applied in NEOParams.

"""
    NEORotation{T}

`phi_rot(theta)` = phi(r,theta) - phi(r,theta=0) from quasi-neutrality, its
theta and radial derivatives, its flux-surface average, and
`dens_fac[is,it]` = n_s(theta)/n_s(theta=0).
"""
struct NEORotation{T<:Real}
    phi_rot::Vector{T}
    phi_rot_deriv::Vector{T}
    phi_rot_rderiv::Vector{T}
    phi_rot_avg::T
    dens_fac::Matrix{T}   # (n_species, n_theta)
end

"""
    rotation_phi(p::NEOParams, geo::NEOGeometry) -> NEORotation

`ROT_solve_phi`: Newton solve of the quasi-neutrality relation at every theta
point (warm-started along theta from `x=0.05`, `|dx|<1e-12`, at most 200
iterations).
"""
function rotation_phi(p::NEOParams{T}, geo::NEOGeometry) where {T<:Real}
    ns, n_theta = p.n_species, p.n_theta
    if p.rotation_model == 1
        return NEORotation{T}(zeros(T, n_theta), zeros(T, n_theta), zeros(T, n_theta), zero(T), ones(T, ns, n_theta))
    end

    z, dens, temp, vth, dlnndr, dlntdr = p.z, p.dens, p.temp, p.vth, p.dlnndr, p.dlntdr
    omega_rot, omega_rot_deriv = p.omega_rot, p.omega_rot_deriv
    bigR, bigR_rderiv, w_theta = geo.bigR, geo.bigR_rderiv, geo.w_theta
    bigR_th0, bigR_th0_rderiv = geo.bigR_th0, geo.bigR_th0_rderiv
    ae = p.ae_flag == 1

    # check for equilibrium-scale QN at theta=0
    sum_zn = zero(T)
    for is in 1:ns
        sum_zn = sum_zn + z[is] * dens[is]
    end
    if ae
        sum_zn = (sum_zn - p.dens_ae) / p.dens_ae
    else
        sum_zn = sum_zn / dens[p.is_ele]
    end
    if abs(_val(sum_zn)) > 1.0e-3
        @warn "NEONative: rotation is being run without quasi-neutral densities (sum Z n = $(_val(sum_zn)) n_e)" maxlog = 1
    end

    # partial component of n/n(theta0) -- the phi_rot component is added after the QN solve
    dens_fac = zeros(T, ns, n_theta)
    for is in 1:ns, it in 1:n_theta
        dens_fac[is, it] = exp(omega_rot^2 * 0.5 / vth[is]^2 * (bigR[it]^2 - bigR_th0^2))
    end

    phi_rot = zeros(T, n_theta)
    phi_rot_avg = zero(T)
    nmax = 200
    x = T(0.05)   # initial guess for phi_rot(1)
    for it in 1:n_theta
        n = 1
        while true
            # Newton's method on the quasi-neutrality relation
            sum_zn = zero(T)
            dsum_zn = zero(T)
            for is in 1:ns
                fac = z[is] * dens[is] * dens_fac[is, it] * exp(-z[is] / temp[is] * x)
                sum_zn = sum_zn + fac
                dsum_zn = dsum_zn - z[is] / temp[is] * fac
            end
            if ae
                fac = -p.dens_ae * exp(1.0 / p.temp_ae * x)
                sum_zn = sum_zn + fac
                dsum_zn = dsum_zn + 1.0 / p.temp_ae * fac
            end
            x0 = x
            x = x0 - sum_zn / dsum_zn
            abs(_val(x - x0)) < 1.0e-12 && break
            n += 1
            n > nmax && break
        end
        n > nmax && error("NEONative: rotation density computation failed to converge")
        phi_rot[it] = x
        phi_rot_avg = phi_rot_avg + w_theta[it] * phi_rot[it]
    end

    # n(theta)/n(0)
    for is in 1:ns, it in 1:n_theta
        dens_fac[is, it] = dens_fac[is, it] * exp(-z[is] / temp[is] * phi_rot[it])
    end

    # dln n(theta)/dr
    dlnndr_fac = zeros(T, ns, n_theta)
    for is in 1:ns, it in 1:n_theta
        dlnndr_fac[is, it] = -(-dlnndr[is] + omega_rot_deriv * omega_rot / vth[is]^2 * (bigR[it]^2 - bigR_th0^2) +
                               omega_rot^2 / vth[is]^2 * (bigR[it] * bigR_rderiv[it] - bigR_th0 * bigR_th0_rderiv) -
                               phi_rot[it] * z[is] / temp[is] * dlntdr[is] +
                               0.5 * omega_rot^2 / vth[is]^2 * dlntdr[is] * (bigR[it]^2 - bigR_th0^2))
    end

    # radial derivative of QN -> dphi*/dr
    phi_rot_rderiv = zeros(T, n_theta)
    for it in 1:n_theta
        acc = zero(T)
        sum_zn = zero(T)
        for is in 1:ns
            fac = z[is] * dens[is] * dens_fac[is, it]
            sum_zn = sum_zn + z[is] / temp[is] * fac
        end
        if ae
            fac = -p.dens_ae * exp(1.0 / p.temp_ae * phi_rot[it])
            sum_zn = sum_zn - 1.0 / p.temp_ae * fac
        end
        for is in 1:ns
            fac = z[is] * dens[is] * dens_fac[is, it]
            acc = acc + fac * (-dlnndr_fac[is, it])
        end
        if ae
            fac = -p.dens_ae * exp(1.0 / p.temp_ae * phi_rot[it])
            acc = acc + fac * (-p.dlnndr_ae + phi_rot[it] / p.temp_ae * p.dlntdr_ae)
        end
        phi_rot_rderiv[it] = acc / sum_zn
    end

    # d phi*/d theta
    phi_rot_deriv = zeros(T, n_theta)
    for it in 1:n_theta
        acc = zero(T)
        for id in -2:2
            jt = thcyc(it + id, n_theta)
            acc = acc + phi_rot[jt] * cderiv(id) / (12.0 * geo.d_theta)
        end
        phi_rot_deriv[it] = acc
    end

    return NEORotation{T}(phi_rot, phi_rot_deriv, phi_rot_rderiv, phi_rot_avg, dens_fac)
end

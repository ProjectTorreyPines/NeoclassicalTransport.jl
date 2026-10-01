# Miller extended-harmonic (MXH) flux-surface geometry: a port of geo_do and
# the linear geo_interp of $GACODE_ROOT/f2py/geo/geo.f90 (geo_model_in=0),
# restricted to the functions NEO's EQUIL_do reads. NEO evaluates geo_do on a
# fixed 2001-point internal grid and interpolates linearly onto its own theta
# grid; both are reproduced here so the port agrees with the Fortran to
# round-off rather than to the ~1e-6 interpolation error.

"""
    MillerGeo{T}

`geo_do` on its internal theta grid (`-pi:pi`, `n_theta` points including both
ends): the vector functions NEO uses plus the scalars `f` (R B_t / B_unit)
and `ffprime`. Build with [`miller_geo`](@ref), sample with [`geo_interp`](@ref).
"""
struct MillerGeo{T<:Real}
    n_theta::Int
    theta::Vector{Float64}
    b::Vector{T}
    dbdt::Vector{T}
    bp::Vector{T}
    bt::Vector{T}
    gsin::Vector{T}
    gcos1::Vector{T}
    gcos2::Vector{T}
    g_theta::Vector{T}
    grad_r::Vector{T}
    l_t::Vector{T}
    nsin::Vector{T}
    bigr::Vector{T}
    bigr_r::Vector{T}
    bigr_t::Vector{T}
    theta_nc::Vector{T}
    f::T
    ffprime::T
end

const GEO_NTHETA_NEO = 2001  # neo_equilibrium.f90: GEO_ntheta_in = 2001

"""
    miller_geo(p::NEOParams; n_theta=GEO_NTHETA_NEO, beta_star=0) -> MillerGeo

`geo_do` for the Miller parameters of `p` (signed `q`, `sign_bunit` as
`geo_signb_in`). `beta_star` only enters for anisotropic species in NEO and is
therefore zero here.
"""
function miller_geo(p::NEOParams{T}; n_theta::Int=GEO_NTHETA_NEO, beta_star_in=zero(T)) where {T<:Real}
    signb = p.sign_bunit
    rmin, rmaj, drmaj, zmag, dzmag = p.rmin, p.rmaj, p.shift, p.zmag, p.s_zmag
    q, s = p.q, p.shear
    kappa, s_kappa, delta, s_delta, zeta, s_zeta = p.kappa, p.s_kappa, p.delta, p.s_delta, p.zeta, p.s_zeta
    cs, s_cs, sn, s_sn = p.shape_cos, p.shape_s_cos, p.shape_sin, p.shape_s_sin

    abs(signb) < 1e-10 && error("NEONative (geo_do): bad value for geo_signb_in")
    abs(_val(delta)) > 1.0 && error("NEONative (geo_do): |delta| > 1")

    n = n_theta
    pi_2 = 8.0 * atan(1.0)
    d_theta = pi_2 / (n - 1)
    ic(i) = mod(i - 1, n - 1) + 1   # periodic index, period n-1

    theta = zeros(n)
    bigr = zeros(T, n)
    bigr_r = zeros(T, n)
    bigr_t = zeros(T, n)
    bigz_t = zeros(T, n)
    bigz_r = zeros(T, n)
    jac_r = zeros(T, n)
    grad_r = zeros(T, n)
    l_t = zeros(T, n)
    r_c = zeros(T, n)
    bigz_l = zeros(T, n)
    nsin = zeros(T, n)
    beta_star = zeros(T, n)

    x = asin(delta)
    for i in 1:n
        th = -0.5 * pi_2 + (i - 1) * d_theta
        theta[i] = th

        # Miller extended harmonic (MHX) parameterization
        a = th + cs[0] + cs[1] * cos(th) + cs[2] * cos(2 * th) + cs[3] * cos(3 * th) + cs[4] * cos(4 * th) +
            cs[5] * cos(5 * th) + cs[6] * cos(6 * th) + x * sin(th) - zeta * sin(2 * th) + sn[3] * sin(3 * th) +
            sn[4] * sin(4 * th) + sn[5] * sin(5 * th) + sn[6] * sin(6 * th)
        a_t = 1.0 - cs[1] * sin(th) - 2 * cs[2] * sin(2 * th) - 3 * cs[3] * sin(3 * th) - 4 * cs[4] * sin(4 * th) -
              5 * cs[5] * sin(5 * th) - 6 * cs[6] * sin(6 * th) + x * cos(th) - 2 * zeta * cos(2 * th) +
              3 * sn[3] * cos(3 * th) + 4 * sn[4] * cos(4 * th) + 5 * sn[5] * cos(5 * th) + 6 * sn[6] * cos(6 * th)
        a_tt = -cs[1] * cos(th) - 4 * cs[2] * cos(2 * th) - 9 * cs[3] * cos(3 * th) - 16 * cs[4] * cos(4 * th) -
               25 * cs[5] * cos(5 * th) - 36 * cs[6] * cos(6 * th) - x * sin(th) + 4 * zeta * sin(2 * th) -
               9 * sn[3] * sin(3 * th) - 16 * sn[4] * sin(4 * th) - 25 * sn[5] * sin(5 * th) - 36 * sn[6] * sin(6 * th)

        bigr[i] = rmaj + rmin * cos(a)
        bigr_r[i] = drmaj + cos(a) - sin(a) * (s_cs[0] + s_cs[1] * cos(th) + s_cs[2] * cos(2 * th) + s_cs[3] * cos(3 * th) +
                                                s_cs[4] * cos(4 * th) + s_cs[5] * cos(5 * th) + s_cs[6] * cos(6 * th) +
                                                s_delta / cos(x) * sin(th) - s_zeta * sin(2 * th) + s_sn[3] * sin(3 * th) +
                                                s_sn[4] * sin(4 * th) + s_sn[5] * sin(5 * th) + s_sn[6] * sin(6 * th))
        bigr_t[i] = -rmin * a_t * sin(a)
        bigr_tt = -rmin * a_t^2 * cos(a) - rmin * a_tt * sin(a)

        a = th
        a_t = 1.0
        a_tt = 0.0
        # Z(theta)
        bigz_r[i] = dzmag + kappa * (1.0 + s_kappa) * sin(a)
        bigz_t[i] = kappa * rmin * cos(a) * a_t
        bigz_tt = -kappa * rmin * sin(a) * a_t^2 + kappa * rmin * cos(a) * a_tt

        g_tt = bigr_t[i]^2 + bigz_t[i]^2
        jac_r[i] = bigr[i] * (bigr_r[i] * bigz_t[i] - bigr_t[i] * bigz_r[i])
        grad_r[i] = bigr[i] * sqrt(g_tt) / jac_r[i]
        l_t[i] = sqrt(g_tt)
        # 1/(du/dl)
        r_c[i] = l_t[i]^3 / (bigr_t[i] * bigz_tt - bigz_t[i] * bigr_tt)
        # cos(u)
        bigz_l[i] = bigz_t[i] / l_t[i]
        nsin[i] = (bigr_r[i] * bigr_t[i] + bigz_r[i] * bigz_t[i]) / l_t[i]
        beta_star[i] = beta_star_in
    end

    # loop integral (1 to n-1) to compute f
    c = zero(T)
    for i in 1:(n-1)
        c = c + l_t[i] / (bigr[i] * grad_r[i])
    end
    f = rmin / (c * d_theta / pi_2)

    bt = zeros(T, n)
    bp = zeros(T, n)
    b = zeros(T, n)
    for i in 1:n
        bt[i] = f / bigr[i]
        bp[i] = (rmin / q) * grad_r[i] / bigr[i]
        b[i] = signb * sqrt(bt[i]^2 + bp[i]^2)
    end

    # 5-point stencils for the derivatives of b
    dbdt = zeros(T, n)
    gsin = zeros(T, n)
    gcos1 = zeros(T, n)
    gcos2 = zeros(T, n)
    g_theta = zeros(T, n)
    for i in 1:n
        b5 = b[ic(i + 2)]
        b4 = b[ic(i + 1)]
        b2 = b[ic(i - 1)]
        b1 = b[ic(i - 2)]
        dbdt[i] = (-b5 + 8.0 * b4 - 8.0 * b2 + b1) / (12.0 * d_theta)
        dbdl = dbdt[i] / l_t[i]
        gsin[i] = bt[i] * rmaj * dbdl / b[i]^2
        gcos1[i] = (bt[i]^2 / bigr[i] * bigz_l[i] + bp[i]^2 / r_c[i]) * rmaj / b[i]^2
        gcos2[i] = 0.5 * (rmaj / b[i]^2) * grad_r[i] * (-beta_star[i])
        g_theta[i] = bigr[i] * b[i] * l_t[i] / (rmin * rmaj * grad_r[i])
    end

    # integrands and integrals E1..E4 (only needed for ffprime here)
    ei = zeros(T, n, 4)
    for i in 1:n
        c = d_theta * l_t[i] / (bigr[i] * grad_r[i])
        ei[i, 1] = c * 2.0 * bt[i] / bp[i] * (rmin / r_c[i] - rmin * bigz_l[i] / bigr[i])
        ei[i, 2] = c * b[i]^2 / bp[i]^2
        ei[i, 3] = c * grad_r[i] * 0.5 / bp[i]^2 * (bt[i] / bp[i]) * beta_star[i]
        ei[i, 4] = -c * grad_r[i] * (bt[i] / bp[i])
    end
    e = zeros(T, n, 4)
    i0 = n ÷ 2 + 1
    for k in 1:4
        e[i0, k] = 0.0
    end
    for i in (i0+1):n, k in 1:4
        e[i, k] = e[i-1, k] + 0.5 * (ei[i-1, k] + ei[i, k])
    end
    for i in (i0-1):-1:1, k in 1:4
        e[i, k] = e[i+1, k] - 0.5 * (ei[i+1, k] + ei[i, k])
    end
    loop = [e[n, k] - e[1, k] for k in 1:4]
    f_prime = (pi_2 * q * s / rmin - loop[1] / rmin + loop[3]) / loop[2]
    ffprime = f * f_prime

    # GS2/NCLASS angle
    theta_nc = zeros(T, n)
    theta_nc[1] = theta[1]
    for i in 2:n
        theta_nc[i] = theta_nc[i-1] + 0.5 * (g_theta[i] + g_theta[i-1]) * d_theta
    end
    theta_nc .= -0.5 * pi_2 .+ pi_2 .* (0.5 * pi_2 .+ theta_nc) ./ (0.5 * pi_2 + theta_nc[n])

    return MillerGeo{T}(n, theta, b, dbdt, bp, bt, gsin, gcos1, gcos2, g_theta, grad_r, l_t, nsin,
        bigr, bigr_r, bigr_t, theta_nc, f, ffprime)
end

"""
    geo_interp(mg::MillerGeo, theta_in) -> NamedTuple

`geo_interp` (general case): linear interpolation of every `MillerGeo` vector
function onto `theta_in`, with the Fortran's end-point clamp.
"""
function geo_interp(mg::MillerGeo{T}, theta_in::AbstractVector) where {T<:Real}
    n_theta = mg.n_theta
    n = length(theta_in)
    names = (:b, :dbdt, :bp, :bt, :gsin, :gcos1, :gcos2, :g_theta, :grad_r, :l_t, :nsin, :bigr, :bigr_r, :bigr_t, :theta_nc)
    out = NamedTuple{names}(ntuple(_ -> zeros(T, n), length(names)))
    dx = mg.theta[2] - mg.theta[1]
    for itheta in 1:n
        theta_0 = theta_in[itheta]
        x0 = theta_0 - mg.theta[1]
        i1 = trunc(Int, x0 / dx) + 1
        i2 = i1 + 1
        x1 = (i1 - 1) * dx
        z = (x0 - x1) / dx
        if i2 > n_theta
            i2 = n_theta
        end
        for name in names
            v = getfield(mg, name)
            getfield(out, name)[itheta] = v[i1] + (v[i2] - v[i1]) * z
        end
    end
    return out
end

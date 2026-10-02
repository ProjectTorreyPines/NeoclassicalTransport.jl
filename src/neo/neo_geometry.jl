# Flux-surface geometry on NEO's theta grid: EQUIL_alloc/EQUIL_do and
# compute_fractrap of neo_equilibrium.f90, for the s-alpha (0),
# large-aspect-ratio (1) and Miller (2) equilibrium models.

# 4th-order centred theta derivative: coefficients cderiv(-2:2), periodic index thcyc
const CDERIV = (1, -8, 0, 8, -1)
@inline cderiv(id::Int) = CDERIV[id+3]
@inline thcyc(it::Int, n_theta::Int) = mod1(it, n_theta)

"""
    NEOGeometry{T}

Equilibrium quantities on the `n_theta` grid `theta = -pi .+ (0:n_theta-1)*2pi/n_theta`
(all normalised to `a` and `B_unit`, sign conventions of `neo_equilibrium.f90`:
`Bmag`, `k_par` carry `sign_bunit`), their values at `theta=0`, the
flux-surface averages NEO keeps, and the trapped fraction.
"""
struct NEOGeometry{T<:Real}
    n_theta::Int
    theta::Vector{Float64}
    d_theta::Float64
    k_par::Vector{T}          # bhat dot grad/a
    v_drift_x::Vector{T}      # radial curvature drift
    v_drift_th::Vector{T}
    gradr::Vector{T}          # |grad r|
    gradpar_gradr::Vector{T}
    w_theta::Vector{T}        # flux-surface-average weights (sum 1)
    Btor::Vector{T}
    Bpol::Vector{T}
    Bmag::Vector{T}
    Bmag_rderiv::Vector{T}
    gradpar_Bmag::Vector{T}
    bigR::Vector{T}
    bigR_rderiv::Vector{T}
    gradpar_bigR::Vector{T}
    theta_nc::Vector{T}
    bigR_th0::T
    bigR_th0_rderiv::T
    gradr_th0::T
    Btor_th0::T
    Bpol_th0::T
    Bmag_th0::T
    Bmag_th0_rderiv::T
    I_div_psip::T             # I(psi)/psi' = q f / r
    Bmag2_avg::T
    Bmag2inv_avg::T
    Btor2_avg::T
    bigRinv_avg::T
    gradpar_Bmag2_avg::T
    ftrap::T
end

"""
    equilibrium(p::NEOParams) -> NEOGeometry

`EQUIL_do` for the flux surface of `p`.
"""
function equilibrium(p::NEOParams{T}) where {T<:Real}
    n_theta = p.n_theta
    r, rmaj, q, rho, sign_bunit = p.rmin, p.rmaj, p.q, p.rho, p.sign_bunit
    d_theta = 2 * pi / n_theta
    theta = [-pi + (it - 1) * d_theta for it in 1:n_theta]

    k_par = zeros(T, n_theta)
    v_drift_x = zeros(T, n_theta)
    v_drift_th = zeros(T, n_theta)
    gradr = zeros(T, n_theta)
    gradpar_gradr = zeros(T, n_theta)
    w_theta = zeros(T, n_theta)
    Btor = zeros(T, n_theta)
    Bpol = zeros(T, n_theta)
    Bmag = zeros(T, n_theta)
    Bmag_rderiv = zeros(T, n_theta)
    gradpar_Bmag = zeros(T, n_theta)
    bigR = zeros(T, n_theta)
    bigR_rderiv = zeros(T, n_theta)
    gradpar_bigR = zeros(T, n_theta)
    theta_nc = zeros(T, n_theta)
    sum = zero(T)

    if p.equilibrium_model == 2
        mg = miller_geo(p)

        # geo params at theta=0
        g0 = geo_interp(mg, [0.0])
        bigR_th0 = g0.bigr[1]
        bigR_th0_rderiv = g0.bigr_r[1]
        gradr_th0 = g0.grad_r[1]
        Btor_th0 = g0.bt[1]
        Bpol_th0 = g0.bp[1]
        Bmag_th0 = g0.b[1]
        Bmag_th0_rderiv = -g0.b[1] / (rmaj * g0.grad_r[1]) * (g0.gcos1[1] + g0.gcos2[1])

        g = geo_interp(mg, theta)
        for it in 1:n_theta
            k_par[it] = 1.0 / (q * rmaj * g.g_theta[it])
            w_theta[it] = g.g_theta[it] / g.b[it]
            sum = sum + w_theta[it]

            bigR[it] = g.bigr[it]
            bigR_rderiv[it] = g.bigr_r[it]
            gradpar_bigR[it] = k_par[it] * g.bigr_t[it]
            Bmag[it] = g.b[it]
            Btor[it] = g.bt[it]
            Bpol[it] = g.bp[it]
            Bmag_rderiv[it] = -g.b[it] / (rmaj * g.grad_r[it]) * (g.gcos1[it] + g.gcos2[it])
            gradpar_Bmag[it] = k_par[it] * g.dbdt[it]
            gradr[it] = g.grad_r[it]
            v_drift_x[it] = -rho / (rmaj * Bmag[it]) * g.grad_r[it] * g.gsin[it]
            v_drift_th[it] = -rho / (rmaj * Bmag[it] * g.l_t[it]) * (g.gcos1[it] * g.bt[it] / Bmag[it] - g.gsin[it] * g.nsin[it] * g.grad_r[it]) * r
            theta_nc[it] = g.theta_nc[it]
        end

        I_div_psip = mg.f * q / r

        for it in 1:n_theta
            acc = zero(T)
            for id in -2:2
                jt = thcyc(it + id, n_theta)
                acc = acc + gradr[jt] * cderiv(id) / (12.0 * d_theta)
            end
            gradpar_gradr[it] = acc * k_par[it]
        end
    else
        # concentric circular geometry
        bigR_th0 = rmaj + r
        bigR_th0_rderiv = one(T)
        gradr_th0 = one(T)
        eps_r = r / rmaj
        if p.equilibrium_model == 1
            # large aspect ratio geometry
            Btor_th0 = 1.0 - eps_r
            Bpol_th0 = r / (q * rmaj) * (1.0 - eps_r)
            Bmag_th0 = (1.0 - eps_r) * sign_bunit
            Bmag_th0_rderiv = -1.0 / rmaj * sign_bunit
        else
            # s-alpha geometry
            Btor_th0 = 1.0 / (1.0 + eps_r)
            Bpol_th0 = r / (q * rmaj) / (1.0 + eps_r)
            Bmag_th0 = (1.0 / (1.0 + eps_r)) * sign_bunit
            Bmag_th0_rderiv = (1.0 / (1.0 + eps_r)^2) * sign_bunit * (-1.0 / rmaj)
        end

        for it in 1:n_theta
            th = theta[it]
            k_par[it] = 1.0 / (q * rmaj) * sign_bunit
            bigR[it] = rmaj + r * cos(th)
            bigR_rderiv[it] = cos(th)
            gradpar_bigR[it] = -r * sin(th) * k_par[it]
            gradr[it] = 1.0
            gradpar_gradr[it] = 0.0
            theta_nc[it] = th

            if p.equilibrium_model == 1
                Bmag[it] = (1.0 - eps_r * cos(th)) * sign_bunit
                Btor[it] = 1.0 - eps_r * cos(th)
                Bpol[it] = r / (q * rmaj) * (1.0 - eps_r * cos(th))
                Bmag_rderiv[it] = -1.0 / rmaj * cos(th) * sign_bunit
                gradpar_Bmag[it] = k_par[it] * eps_r * sin(th) * sign_bunit
                v_drift_x[it] = -rho / rmaj * sin(th) / (1.0 - eps_r * cos(th))^2
                v_drift_th[it] = -rho / (rmaj) * cos(th) / (1.0 - eps_r * cos(th))^2
            else
                Bmag[it] = (1.0 / (1.0 + eps_r * cos(th))) * sign_bunit
                Btor[it] = 1.0 / (1.0 + eps_r * cos(th))
                Bpol[it] = r / (q * rmaj) / (1.0 + eps_r * cos(th))
                Bmag_rderiv[it] = (1.0 / (1.0 + eps_r * cos(th))^2) * sign_bunit * (-1.0 / rmaj * cos(th))
                gradpar_Bmag[it] = k_par[it] * eps_r * sin(th) / (1.0 + eps_r * cos(th))^2 * sign_bunit
                v_drift_x[it] = -rho / rmaj * sin(th)
                v_drift_th[it] = -rho / rmaj * cos(th)
            end

            # flux-surface average weights
            w_theta[it] = 1.0 * sign_bunit / Bmag[it]
            sum = sum + w_theta[it]
        end

        I_div_psip = rmaj * q / r
    end

    for it in 1:n_theta
        w_theta[it] = w_theta[it] / sum
    end

    Bmag2_avg = zero(T)
    Bmag2inv_avg = zero(T)
    gradpar_Bmag2_avg = zero(T)
    Btor2_avg = zero(T)
    bigRinv_avg = zero(T)
    for it in 1:n_theta
        Bmag2_avg = Bmag2_avg + w_theta[it] * Bmag[it]^2
        Bmag2inv_avg = Bmag2inv_avg + w_theta[it] / Bmag[it]^2
        gradpar_Bmag2_avg = gradpar_Bmag2_avg + w_theta[it] * gradpar_Bmag[it]^2
        Btor2_avg = Btor2_avg + w_theta[it] * Btor[it]^2
        bigRinv_avg = bigRinv_avg + w_theta[it] / bigR[it]
    end

    ftrap = compute_fractrap(Bmag, w_theta, sign_bunit, Bmag2_avg)

    return NEOGeometry{T}(n_theta, theta, d_theta, k_par, v_drift_x, v_drift_th, gradr, gradpar_gradr, w_theta,
        Btor, Bpol, Bmag, Bmag_rderiv, gradpar_Bmag, bigR, bigR_rderiv, gradpar_bigR, theta_nc,
        bigR_th0, bigR_th0_rderiv, gradr_th0, Btor_th0, Bpol_th0, Bmag_th0, Bmag_th0_rderiv,
        I_div_psip, Bmag2_avg, Bmag2inv_avg, Btor2_avg, bigRinv_avg, gradpar_Bmag2_avg, ftrap)
end

"""
    compute_fractrap(Bmag, w_theta, sign_bunit, Bmag2_avg)

Trapped-particle fraction (`compute_fractrap`): open Newton-Cotes rule in
lambda (500 points) over the closed flux-surface average in theta.
"""
function compute_fractrap(Bmag::AbstractVector{T}, w_theta::AbstractVector{T}, sign_bunit, Bmag2_avg) where {T<:Real}
    n_theta = length(Bmag)
    nlambda = 500
    Bmax = Bmag[1] * sign_bunit
    for it in 2:n_theta
        if _val(Bmag[it] * sign_bunit) > _val(Bmax)
            Bmax = Bmag[it] * sign_bunit
        end
    end

    dlambda = 1.0 / (nlambda - 1)
    ft = zero(T)
    for i in 2:(nlambda-1)
        lam = (i - 1) * dlambda
        # open integration for lambda
        if i == 2 || i == nlambda - 1
            fac_lambda = 23.0 / 12.0
        elseif i == 3 || i == nlambda - 2
            fac_lambda = 7.0 / 12.0
        else
            fac_lambda = 1.0
        end
        # closed integration for th: <sqrt(1-lambda*B/Bmax)>
        sum_th = zero(T)
        for it in 1:n_theta
            sum_th = sum_th + w_theta[it] * sqrt(1.0 - lam * Bmag[it] * sign_bunit / Bmax)
        end
        ft = ft + fac_lambda * lam / sum_th
    end
    ft = ft * dlambda * 0.75 * Bmag2_avg / Bmax^2
    return 1.0 - ft
end

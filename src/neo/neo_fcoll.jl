# Field-particle collision integrals F(m,n,λ) and their complements
# FB(m,n,λ): a literal port of neo_compute_fcoll.f90.
#
#   G(m,n,λ) = GB(n,m,1/λ),   G(m,n,λ) + GB(m,n,λ) = B((m+1)/2, (n+1)/2)
#
# f and fb are defined where m+n >= -1 and zero elsewhere. Branch selections
# compare on primal values so the port stays ForwardDiff-generic.

const FCOLL_LAMBDA_LARGE = 10.0
const FCOLL_EXACT_BITS = 256

_val(x) = ForwardDiff.value(x)

"""
    compute_fcoll(m0, lambda; exact=false) -> (f, fb)

NEO's `neo_compute_fcoll`: the tables `f[m,n]` and `fb[m,n]` for
`m, n ∈ -m0:m0` (as `OffsetMatrix`), at mass-ratio parameter
`lambda = (vth_a/vth_b)^2`.

NEO's backward recursions for the negative-index entries (stage C) are
unstable in double precision when `0.1 <= lambda <= 10` (unequal-mass
ion pairs and like-species collisions): the entries at large negative `m`
(or `n`) and large `n` (or `m`) that feed the high-Legendre-order
field-particle blocks lose from ~6 digits (`lambda=1`) to everything
(`lambda=0.17`, ix >~ 10). `exact=true` evaluates the same recursions in
`BigFloat` (`FCOLL_EXACT_BITS` bits) and rounds the tables to `T`; this is
what the Fortran would give with enough precision, at ~10x the cost of the
tables (and only for `Float64` input).
"""
function compute_fcoll(m0::Int, lambda_in::T; exact::Bool=false) where {T<:Real}
    if exact
        T === Float64 || error("compute_fcoll: exact=true needs a Float64 lambda (got $T)")
        fB, fbB = setprecision(BigFloat, FCOLL_EXACT_BITS) do
            return compute_fcoll(m0, BigFloat(lambda_in))
        end
        return OffsetArray(Float64.(parent(fB)), -m0:m0, -m0:m0), OffsetArray(Float64.(parent(fbB)), -m0:m0, -m0:m0)
    end
    g = OffsetArray(zeros(T, 2m0 + 1, 2m0 + 1), -m0:m0, -m0:m0)
    gb = OffsetArray(zeros(T, 2m0 + 1, 2m0 + 1), -m0:m0, -m0:m0)
    f = OffsetArray(zeros(T, 2m0 + 1, 2m0 + 1), -m0:m0, -m0:m0)
    fb = OffsetArray(zeros(T, 2m0 + 1, 2m0 + 1), -m0:m0, -m0:m0)

    # gam[i] = Γ(i/2)
    gam = zeros(2m0 + 2)
    gam[2] = 1.0
    for i in 1:m0
        gam[2*(i+1)] = i * gam[2*i]
    end
    gam[1] = sqrt(pi)
    for i in 1:m0
        gam[2*i+1] = (i - 0.5) * gam[2*i-1]
    end

    beta = OffsetArray(zeros(m0 + 1, m0 + 1), 0:m0, 0:m0)
    for m in 0:m0, n in 0:m0
        beta[m, n] = gam[m+1] * gam[n+1] / gam[m+n+2]
    end

    # STAGE A: GB(m,n,lambda)
    lambda = lambda_in
    r = 1.0 / (1.0 + lambda)

    for m in 0:2:m0, n in 1:2:m0
        # CASE 1: m even, n odd
        x = zero(T)
        for i in 0:((n-1)÷2)
            x = x + (lambda / (1.0 + lambda))^i / (1.0 + lambda)^((m + 1) / 2.0) * gam[2*i+m+1] / gam[m+n+2] * gam[n+1] / gam[2*i+2]
        end
        gb[m, n] = x
    end

    for m in 1:2:m0, n in 0:2:m0
        # CASE 2: m odd, n even
        if _val(lambda) < FCOLL_LAMBDA_LARGE
            x = zero(T)
            for i in 0:((m-1)÷2)
                x = x + r^i * gam[2*i+1] / (gam[2*i+2] * gam[1])
            end
            gb[m, n] = beta[m, n] * (1.0 - x * sqrt(lambda / (1.0 + lambda)))
        else
            # complementary large-lambda sum (16 digits of precision)
            c = gam[m] / (gam[m+1] * gam[1]) * r^((m - 1) ÷ 2)
            x = zero(T)
            for i in ((m+1)÷2):((m+1)÷2+trunc(Int, 16 * log(10.0) / log(_val(lambda))))
                c = c * (i - 0.5) / (i) * r
                x = x + c
            end
            gb[m, n] = beta[m, n] * x * sqrt(lambda / (1.0 + lambda))
        end

        x = zero(T)
        for i in 0:(n÷2-1)
            x = x + (lambda / (1 + lambda))^(i + 0.5) / (1 + lambda)^(0.5 + m / 2.0) * gam[2*i+2+m] / gam[m+n+2] * gam[n+1] / gam[2*i+3]
        end
        gb[m, n] = gb[m, n] + x
    end

    # STAGE B: G(n,m,lambda) = GB(m,n,1/lambda)
    lambda = 1 / lambda_in
    r = 1.0 / (1.0 + lambda)

    for m in 0:2:m0, n in 1:2:m0
        x = zero(T)
        for i in 0:((n-1)÷2)
            x = x + (lambda / (1.0 + lambda))^i / (1.0 + lambda)^((m + 1) / 2.0) * gam[2*i+m+1] / gam[m+n+2] * gam[n+1] / gam[2*i+2]
        end
        g[n, m] = x
    end

    for m in 1:2:m0, n in 0:2:m0
        if _val(lambda) < FCOLL_LAMBDA_LARGE
            x = zero(T)
            for i in 0:((m-1)÷2)
                x = x + r^i * gam[2*i+1] / (gam[2*i+2] * gam[1])
            end
            g[n, m] = beta[m, n] * (1.0 - x * sqrt(lambda / (1.0 + lambda)))
        else
            c = gam[m] / (gam[m+1] * gam[1]) * r^((m - 1) ÷ 2)
            x = zero(T)
            for i in ((m+1)÷2):((m+1)÷2+trunc(Int, 16 * log(10.0) / log(_val(lambda))))
                c = c * (i - 0.5) / i * r
                x = x + c
            end
            g[n, m] = beta[m, n] * x * sqrt(lambda / (1.0 + lambda))
        end

        x = zero(T)
        for i in 0:(n÷2-1)
            x = x + (lambda / (1 + lambda))^(i + 0.5) / (1 + lambda)^(0.5 + m / 2.0) * gam[2*i+2+m] / gam[m+n+2] * gam[n+1] / gam[2*i+3]
        end
        g[n, m] = g[n, m] + x
    end

    # STAGE C: special (negative) elements
    lambda = lambda_in

    # Case 1: special elements of g
    if _val(lambda) < 1 / FCOLL_LAMBDA_LARGE
        # asymptotic series for small lambda
        for n in 2:m0
            m = 1 - n
            c = 1.0 / (n + 1)
            gmn = T(c)
            for i in 1:16
                c = -c * (1 + 0.5 / i) * (2 * i + n - 1.0) / (2 * i + n + 1.0) * lambda
                gmn = gmn + c
            end
            g[m, n] = 2 * gmn * lambda^(0.5 * (n + 1))

            for m in (3-n):2:-1
                g[m, n] = 2.0 / (m + n) * (0.5 * (m - 1) * g[m-2, n] + lambda^((n + 1) / 2.0) / (1 + lambda)^((m + n) / 2.0))
            end
        end
    else
        g[-1, 0] = log((sqrt(1 + lambda) + sqrt(lambda)) / (sqrt(1.0 + lambda) - sqrt(lambda)))

        m = -1
        for n in 2:m0
            g[m, n] = 2.0 / (m + n) * (0.5 * (n - 1) * g[m, n-2] - lambda^((n - 1) / 2.0) / (1 + lambda)^((m + n) / 2.0))
        end

        for n in 0:m0
            for m in 0:-1:(-m0+2)
                if m + n >= 1
                    g[m-2, n] = 2.0 / (m - 1.0) * (0.5 * (m + n) * g[m, n] - lambda^((n + 1) / 2.0) / (1 + lambda)^((m + n) / 2.0))
                end
            end
        end
    end

    # Case 2: special elements of gb
    if _val(lambda) > FCOLL_LAMBDA_LARGE
        # small-lambda asymptotic series for g to get gb: gb(n,m,lambda) = g(m,n,1/lambda)
        lambda = 1 / lambda_in

        for n in 2:m0
            m = 1 - n
            c = 1.0 / (n + 1)
            gbnm = T(c)
            for i in 1:16
                c = -c * (1 + 0.5 / i) * (2 * i + n - 1.0) / (2 * i + n + 1.0) * lambda
                gbnm = gbnm + c
            end
            gb[n, m] = 2 * gbnm * lambda^(0.5 * (n + 1))

            for m in (3-n):2:-1
                gb[n, m] = 2.0 / (m + n) * (0.5 * (m - 1) * gb[n, m-2] + lambda^((n + 1) / 2.0) / (1 + lambda)^((m + n) / 2.0))
            end
        end

        lambda = lambda_in
    else
        gb[0, -1] = log((sqrt(1.0 + lambda) + 1.0) / (sqrt(1.0 + lambda) - 1.0))

        n = -1
        for m in 2:m0
            gb[m, n] = 2.0 / (m + n) * (0.5 * (m - 1) * gb[m-2, n] - lambda^((n + 1) / 2.0) / (1 + lambda)^((m + n) / 2.0))
        end

        for m in 0:m0
            for n in 0:-1:(-m0+2)
                if m + n >= 1
                    gb[m, n-2] = 2 / (n - 1.0) * (0.5 * (m + n) * gb[m, n] - lambda^((n - 1) / 2.0) / (1 + lambda)^((m + n) / 2.0))
                end
            end
        end
    end

    # f, fb from g, gb
    for m in -m0:m0, n in -m0:m0
        if m + n >= -1
            f[m, n] = g[m, n] * 0.25 * gam[m+n+2] / lambda^((n + 1) / 2.0)
            fb[m, n] = gb[m, n] * 0.25 * gam[m+n+2] / lambda^((n + 1) / 2.0)
        end
    end

    return f, fb
end

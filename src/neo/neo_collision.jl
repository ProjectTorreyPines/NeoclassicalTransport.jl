# Collision matrices in the Laguerre/Legendre basis: ENERGY_coll_ints of
# neo_energy_grid.f90 for every NEO collision model:
#   1  Connor (reduced, mass-ratio limits)
#   2  reduced Hirshman-Sigmar (HS0)
#   3  full Hirshman-Sigmar
#   4  full linearized Fokker-Planck
#   5  full FP test-particle part + ad-hoc (rank-1 momentum/energy restoring) field part
# The test-particle part (pitch-angle scattering + energy diffusion/drag) and
# the field-particle part (Rosenbluth potentials, in closed form through the
# fcoll tables) are diagonal in the Legendre order.
#
# emat_coll_test[is,js,ie,je,ix] acts on g_is (summed over the background js
# at assembly), emat_coll_field[is,js,ie,je,ix] acts on g_js.

"""
    NEOCollision{T}

`test` and `field`: `(1:ns, 1:ns, 0:n_energy, 0:n_energy, 0:n_xi)` collision
matrices (`emat_coll_test`, `emat_coll_field`), see [`collision_ints`](@ref).
"""
struct NEOCollision{T<:Real}
    test::OffsetArray{T,5,Array{T,5}}
    field::OffsetArray{T,5,Array{T,5}}
end

fcoll_margin(p::NEOParams) = 2 * E_ALPHA * p.n_energy + 2 * p.n_xi + 4

"""
    collision_ints(p::NEOParams, basis::NEOBasis=NEOBasis(p); serial=false, fcoll_exact=false) -> NEOCollision

`ENERGY_coll_ints`: the collision matrices of every species pair for
`p.collision_model`. The pairs are independent and run on threads unless
`serial=true` (or only one thread is available); each pair fills its own
`[is,js,…]` slice. `fcoll_exact=true` evaluates the field-particle integral
tables in extended precision (see [`compute_fcoll`](@ref)) instead of
reproducing NEO's double-precision recursion.
"""
function collision_ints(p::NEOParams{T}, basis::NEOBasis=NEOBasis(p); serial::Bool=false, fcoll_exact::Bool=false) where {T<:Real}
    ns, ne, nxi = p.n_species, p.n_energy, p.n_xi
    if p.collision_model == 3 && any(abs(_val(p.temp[is] - p.temp[1])) > 1e-3 for is in 2:ns)
        @warn "NEONative: full HS collisions (collision_model=3) with unequal temperatures"
    end
    test = OffsetArray(zeros(T, ns, ns, ne + 1, ne + 1, nxi + 1), 1:ns, 1:ns, 0:ne, 0:ne, 0:nxi)
    field = OffsetArray(zeros(T, ns, ns, ne + 1, ne + 1, nxi + 1), 1:ns, 1:ns, 0:ne, 0:ne, 0:nxi)
    pairs = [(is, js) for js in 1:ns for is in 1:ns]
    if serial || Threads.nthreads() == 1
        for (is, js) in pairs
            _coll_pair!(test, field, is, js, p, basis; fcoll_exact)
        end
    else
        Threads.@threads for ip in eachindex(pairs)
            is, js = pairs[ip]
            _coll_pair!(test, field, is, js, p, basis; fcoll_exact)
        end
    end
    if p.collision_model == 5
        _adhoc_field!(test, field, p, basis)
    end
    return NEOCollision{T}(test, field)
end

# one (is,js) block, loop structure and operation order of the Fortran kept
function _coll_pair!(test, field, is::Int, js::Int, p::NEOParams{T}, basis::NEOBasis; fcoll_exact::Bool=false) where {T<:Real}
    ne, nxi = p.n_energy, p.n_xi
    cm = p.collision_model
    Z, mass, dens, temp, vth = p.z, p.mass, p.dens, p.temp, p.vth
    mygamma2 = basis.mygamma2
    beta_l = basis.xi_beta_l
    lag = basis.lag
    fmarg = fcoll_margin(p)

    # (pol part of dens from rotation is added in the coll term of the kinetic equation)
    tauinv_ab = p.nu[is] * (1.0 * Z[js])^2 / (1.0 * Z[is])^2 * dens[js] / dens[is]
    tauinv_ba = p.nu[js] * (1.0 * Z[is])^2 / (1.0 * Z[js])^2 * dens[is] / dens[js]
    lambda = (vth[is] / vth[js])^2
    fcoll, fcoll_bar = compute_fcoll(fmarg, lambda; exact=fcoll_exact)
    if cm in (1, 2, 3)
        fcollinv, fcollinv_bar = compute_fcoll(fmarg, 1.0 / lambda; exact=fcoll_exact)
    else
        fcollinv = fcoll   # unused
    end

    sqlp = sqrt(lambda / pi)
    mrat = mass[is] / mass[js]
    trat = 1.0 - temp[is] / temp[js]
    lam15 = lambda^1.5
    # (m n vth)_b / (m n vth)_a, the field-particle normalization of the reduced models
    mnv = (mass[js] * dens[js] * vth[js]) / (mass[is] * dens[is] * vth[is])
    same_mass = is == js || abs(_val(mass[is] - mass[js])) < eps()

    for ie in 0:ne, je in 0:ne, ix in 0:nxi
        jx = ix
        t = zero(T)
        fld = zero(T)
        for ke in 0:ie, me in 0:je
            zarg0 = (-1.0)^(ke + me)
            zarg1 = lag[ie, ke, ix] * lag[je, me, jx]
            zarg2 = mygamma2[2+2*ke] * mygamma2[2+2*me]
            xab = E_ALPHA * (ke + me) + beta_l[ix] + beta_l[jx]
            xa = E_ALPHA * ke + beta_l[ix]
            yb = E_ALPHA * me + beta_l[jx]

            if cm == 1
                # ---- Connor model
                if same_mass
                    # case 1: ma = mb
                    t = t - tauinv_ab * sqlp * ix * (ix + 1) * zarg0 * zarg1 * (fcoll[xab-1, 0] - fcoll[xab-3, 2]) / zarg2
                    if ix == 1
                        fld = fld + mnv * 2.0 * tauinv_ba / sqrt(lambda * pi) * zarg0 * zarg1 * (fcoll[xa, 0] - fcoll[xa-2, 2]) / zarg2 /
                                    (fcoll[1, 0] - fcoll[-1, 2]) * (fcollinv[yb, 0] - fcollinv[yb-2, 2])
                    end
                elseif _val(mass[is]) < _val(mass[js])
                    # case 2: ma < mb (e-i, i-z)
                    if xab > 0
                        t = t - tauinv_ab * 0.25 * ix * (ix + 1) * zarg0 * zarg1 * mygamma2[xab] / zarg2
                    end
                    if ix == 1
                        fld = fld + mnv * (2.0 / 3.0) * tauinv_ba * temp[js] / temp[is] / sqrt(lambda * pi) * zarg0 * zarg1 *
                                    mygamma2[xa+1] / mygamma2[2] / zarg2 * mygamma2[yb+4]
                    end
                else
                    # case 3: ma > mb (i-e, z-i)
                    t = t - tauinv_ab * sqlp * temp[is] / temp[js] * (1.0 / 3.0) * ix * (ix + 1) * zarg0 * zarg1 * mygamma2[xab+3] / zarg2
                    if ix == 1
                        fld = fld + mnv * 0.5 * tauinv_ba * temp[js] / temp[is] * zarg0 * zarg1 * mygamma2[xa+4] / mygamma2[5] / zarg2 * mygamma2[yb+1]
                    end
                end

            elseif cm == 2
                # ---- HS0
                t = t - tauinv_ab * sqlp * ix * (ix + 1) * zarg0 * zarg1 * (fcoll[xab-1, 0] - fcoll[xab-3, 2]) / zarg2
                if ix == 1
                    rs = mnv * 4.0 * tauinv_ba * temp[js] / temp[is] * (1.0 + mrat) / sqrt(lambda * pi) *
                         zarg0 * zarg1 * fcoll[xa, 2] / zarg2 / fcoll[1, 2] * fcollinv[yb, 2]
                    ru = 2.0 * tauinv_ab * sqlp * zarg0 * zarg1 *
                         (fcoll[xab-1, 0] - fcoll[xab-3, 2] - 2.0 * temp[is] / temp[js] * (1.0 + mass[js] / mass[is]) * fcoll[xab-1, 2]) / zarg2
                    fld = fld + rs
                    t = t + ru
                end

            elseif cm == 3
                # ---- full HS
                t = t - tauinv_ab * sqlp * ix * (ix + 1) * zarg0 * zarg1 * (fcoll[xab-1, 0] - fcoll[xab-3, 2]) / zarg2

                if ix == 1
                    # slowing-down
                    rs = mnv * 4.0 * tauinv_ba * temp[js] / temp[is] * (1.0 + mrat) / sqrt(lambda * pi) *
                         zarg0 * zarg1 * fcoll[xa, 2] / zarg2 / fcoll[1, 2] * fcollinv[yb, 2]
                    # u-restoring
                    ru = 2.0 * tauinv_ab * sqlp * zarg0 * zarg1 *
                         (fcoll[xab-1, 0] - fcoll[xab-3, 2] - 2.0 * temp[is] / temp[js] * (1.0 + mass[js] / mass[is]) * fcoll[xab-1, 2]) / zarg2
                    # h-heating friction
                    r1 = 4.0 * (2.0 * mrat - 1.0) * fcoll[xa, 2] + 8.0 * temp[is] / temp[js] * (1.0 + mass[js] / mass[is]) * fcoll[xa+2, 2] -
                         4.0 * (1.0 + mrat) * (fcoll[xa+2, 0] - fcoll[xa, 2])
                    x3 = 3
                    r2 = 4.0 * (2.0 * mrat - 1.0) * fcoll[x3, 2] + 8.0 * temp[is] / temp[js] * (1.0 + mass[js] / mass[is]) * fcoll[x3+2, 2] -
                         4.0 * (1.0 + mrat) * (fcoll[x3+2, 0] - fcoll[x3, 2])
                    r3 = 4.0 * (2.0 * mass[js] / mass[is] - 1.0) * fcollinv[yb, 2] + 8.0 * temp[js] / temp[is] * (1.0 + mrat) * fcollinv[yb+2, 2] -
                         4.0 * (1.0 + mass[js] / mass[is]) * (fcollinv[yb+2, 0] - fcollinv[yb, 2])
                    rh = (mass[js] * dens[js] * vth[js]^3) / (mass[is] * dens[is] * vth[is]^3) * 1.5 * tauinv_ba / sqrt(lambda * pi) *
                         zarg0 * zarg1 * r1 / zarg2 / r2 * r3
                    # k-heating friction
                    r1 = 12.0 * fcoll[xa, 2] - 8.0 * temp[is] / temp[js] * (1.0 + mass[js] / mass[is]) * fcoll[xa+2, 2] + 4.0 * (fcoll[xa+2, 0] - fcoll[xa, 2])
                    r2 = 12.0 * fcoll[x3, 2] - 8.0 * temp[is] / temp[js] * (1.0 + mass[js] / mass[is]) * fcoll[x3+2, 2] + 4.0 * (fcoll[x3+2, 0] - fcoll[x3, 2])
                    r3 = 12.0 * fcoll[yb, 2] - 8.0 * temp[is] / temp[js] * (1.0 + mass[js] / mass[is]) * fcoll[yb+2, 2] + 4.0 * (fcoll[yb+2, 0] - fcoll[yb, 2])
                    rk = tauinv_ab * sqlp * zarg0 * zarg1 * r1 / zarg2 / r2 * r3
                    fld = fld + rs + rh
                    t = t + ru + rk
                end

                if ix == 2
                    # nup-energy restoring
                    r1 = 8.0 * (1.0 + 1.5 * mrat) * fcoll[xa-1, 2] + 8.0 * temp[is] / temp[js] * (1.0 + mass[js] / mass[is]) * fcoll[xa+1, 2] -
                         4.0 * (1.0 + 1.5 * mrat) * (fcoll[xa+1, 0] - fcoll[xa-1, 2])
                    x2 = 2
                    r2 = 8.0 * (1.0 + 1.5 * mrat) * fcoll[x2-1, 2] + 8.0 * temp[is] / temp[js] * (1.0 + mass[js] / mass[is]) * fcoll[x2+1, 2] -
                         4.0 * (1.0 + 1.5 * mrat) * (fcoll[x2+1, 0] - fcoll[x2-1, 2])
                    r3 = 8.0 * (1.0 + 1.5 * mass[js] / mass[is]) * fcollinv[yb-1, 2] + 8.0 * temp[js] / temp[is] * (1.0 + mrat) * fcollinv[yb+1, 2] -
                         4.0 * (1.0 + 1.5 * mass[js] / mass[is]) * (fcollinv[yb+1, 0] - fcollinv[yb-1, 2])
                    rp = (temp[js] * dens[js]) / (temp[is] * dens[is]) * tauinv_ba / sqrt(lambda * pi) * zarg0 * zarg1 * r1 / zarg2 / r2 * r3
                    fld = fld + rp
                    # pi-energy restoring
                    r1 = -8.0 * temp[is] / temp[js] * (1.0 + mass[js] / mass[is]) * fcoll[xa+yb-1, 2] + 4.0 * fcoll[xa+yb-1, 0]
                    rpi = tauinv_ab * sqlp * zarg0 * zarg1 * r1 / zarg2
                    t = t + rpi
                end

                if ix == 0
                    # energy diffusion
                    r2 = 8.0 * temp[js] / temp[is] * (1.0 + mrat) * fcollinv[yb+1, 2] - 4.0 * fcollinv[yb+1, 0]
                    x2 = 2
                    r3 = 8.0 * temp[is] / temp[js] * (1.0 + mass[js] / mass[is]) * fcoll[x2+1, 2] - 4.0 * fcoll[x2+1, 0]
                    r1 = (temp[js] * dens[js]) / (temp[is] * dens[is]) * tauinv_ba / tauinv_ab / lambda * r2 / r3
                    rd = -4.0 * tauinv_ab * sqlp * xa * zarg0 * zarg1 * 0.5 * yb * fcoll[xa+yb-3, 2] / zarg2
                    rv = 4.0 * tauinv_ab * sqlp * xa * zarg0 * zarg1 * r1 * fcoll[xa-1, 2] / zarg2
                    t = t + rd
                    fld = fld + rv
                end

            else
                # ---- full linearized FP op (models 4 and 5): test-particle part
                t = t - tauinv_ab * sqlp * ix * (ix + 1) * zarg0 * zarg1 * (fcoll[xab-1, 0] - fcoll[xab-3, 2]) / zarg2

                rd = 4.0 * tauinv_ab * sqlp * (E_ALPHA * ke + beta_l[ix]) * zarg0 * zarg1 *
                     (trat * fcoll[xab-1, 2] - 0.5 * (E_ALPHA * me + beta_l[jx]) * fcoll[xab-3, 2]) / zarg2
                t = t + rd

                # field-particle part (Rosenbluth potentials, Legendre order ix)
                r1 = mrat * (1.0 + lambda)^(-0.5 * (xab + 3)) * mygamma2[xab+3]

                r2 = -2.0 / (ix + 0.5) * (mrat - ix * (1.0 - mrat)) * fcoll[xa-ix+1, yb+jx+2]
                r3 = -2.0 / (ix + 0.5) * (1.0 + ix * (1.0 - mrat)) * fcoll_bar[xa+ix+2, yb-jx+1]
                r4 = -ix * (ix - 1.0) / (ix * ix - 0.25) * fcoll[xa-ix+3, yb+jx+2]
                r5 = -ix * (ix - 1.0) / (ix * ix - 0.25) * fcoll_bar[xa+ix+2, yb-jx+3]
                r6 = (ix + 1.0) * (ix + 2.0) / (ix + 1.5) / (ix + 0.5) * fcoll[xa-ix+1, yb+jx+4]
                r7 = (ix + 1.0) * (ix + 2.0) / (ix + 1.5) / (ix + 0.5) * fcoll_bar[xa+ix+4, yb-jx+1]

                fld = fld + tauinv_ab * 2.0 / sqrt(pi) * lam15 * lambda^(0.5 * (E_ALPHA * me + beta_l[ix])) *
                            zarg0 * zarg1 * (r1 + r2 + r3 + r4 + r5 + r6 + r7) / zarg2
            end
        end
        test[is, js, ie, je, ix] = t
        field[is, js, ie, je, ix] = fld
    end
    return nothing
end

# collision_model=5: replace the field-particle operator with the ad-hoc
# momentum- (ix=1) and energy- (ix=0) restoring rank-1 operators built from
# the test-particle matrix
function _adhoc_field!(test, field, p::NEOParams{T}, basis::NEOBasis) where {T<:Real}
    ns, ne = p.n_species, p.n_energy
    mass, dens, temp, vth = p.mass, p.dens, p.temp, p.vth
    beta_l = basis.xi_beta_l
    fill!(field, zero(T))
    for is in 1:ns, js in 1:ns, ie in 0:ne, je in 0:ne
        ix = 1
        k1 = (1 - beta_l[ix]) ÷ E_ALPHA
        if abs(_val(test[is, js, k1, k1, ix])) > eps()
            rs = -(mass[js] * dens[js] * vth[js]) / (mass[is] * dens[is] * vth[is]) *
                 test[is, js, ie, k1, ix] * test[js, is, k1, je, ix] / test[is, js, k1, k1, ix]
        else
            rs = zero(T)
        end
        field[is, js, ie, je, ix] = rs

        ix = 0
        k0 = (2 - beta_l[ix]) ÷ E_ALPHA
        if abs(_val(test[is, js, k0, k0, ix])) > eps()
            ru = -(dens[js] * temp[js]) / (dens[is] * temp[is]) *
                 test[is, js, ie, k0, ix] * test[js, is, k0, je, ix] / test[is, js, k0, k0, ix]
        else
            ru = zero(T)
        end
        field[is, js, ie, je, ix] = ru
    end
    return nothing
end

"""
    collision_ints_mono(p::NEOParams, basis::NEOBasis=NEOBasis(p)) -> (test_mono, field_mono)

The monomial-basis full-FP matrices of `write_fullcoll_mono` (what NEO writes
to `out.neo.diagnostic_coll*` with `WRITE_CMOMENTS_FLAG=1`), divided by
`tauinv_ab` like the file. Diagnostic only.
"""
function collision_ints_mono(p::NEOParams{T}, basis::NEOBasis=NEOBasis(p)) where {T<:Real}
    ns, ne, nxi = p.n_species, p.n_energy, p.n_xi
    Z, mass, dens, temp, vth = p.z, p.mass, p.dens, p.temp, p.vth
    mygamma2 = basis.mygamma2
    beta_l = basis.xi_beta_l
    fmarg = fcoll_margin(p)
    test = OffsetArray(zeros(T, ns, ns, ne + 1, ne + 1, nxi + 1), 1:ns, 1:ns, 0:ne, 0:ne, 0:nxi)
    field = OffsetArray(zeros(T, ns, ns, ne + 1, ne + 1, nxi + 1), 1:ns, 1:ns, 0:ne, 0:ne, 0:nxi)

    for is in 1:ns, js in 1:ns
        tauinv_ab = p.nu[is] * (1.0 * Z[js])^2 / (1.0 * Z[is])^2 * dens[js] / dens[is]
        lambda = (vth[is] / vth[js])^2
        fcoll, fcoll_bar = compute_fcoll(fmarg, lambda)
        mrat = mass[is] / mass[js]
        for ie in 0:ne, je in 0:ne, ix in 0:nxi
            jx = ix
            xab = E_ALPHA * (ie + je) + beta_l[ix] + beta_l[jx]
            t = -tauinv_ab * sqrt(lambda / pi) * ix * (ix + 1) * (fcoll[xab-1, 0] - fcoll[xab-3, 2])
            rd = 4.0 * tauinv_ab * sqrt(lambda / pi) * (E_ALPHA * ie + beta_l[ix]) *
                 ((1.0 - temp[is] / temp[js]) * fcoll[xab-1, 2] - 0.5 * (E_ALPHA * je + beta_l[jx]) * fcoll[xab-3, 2])
            t = t + rd

            r1 = mrat * (1.0 + lambda)^(-0.5 * (xab + 3)) * mygamma2[xab+3]
            xa = E_ALPHA * ie + beta_l[ix]
            yb = E_ALPHA * je + beta_l[jx]
            r2 = -2.0 / (ix + 0.5) * (mrat - ix * (1.0 - mrat)) * fcoll[xa-ix+1, yb+jx+2]
            r3 = -2.0 / (ix + 0.5) * (1.0 + ix * (1.0 - mrat)) * fcoll_bar[xa+ix+2, yb-jx+1]
            r4 = -ix * (ix - 1.0) / (ix * ix - 0.25) * fcoll[xa-ix+3, yb+jx+2]
            r5 = -ix * (ix - 1.0) / (ix * ix - 0.25) * fcoll_bar[xa+ix+2, yb-jx+3]
            r6 = (ix + 1.0) * (ix + 2.0) / (ix + 1.5) / (ix + 0.5) * fcoll[xa-ix+1, yb+jx+4]
            r7 = (ix + 1.0) * (ix + 2.0) / (ix + 1.5) / (ix + 0.5) * fcoll_bar[xa+ix+4, yb-jx+1]
            fld = tauinv_ab * 2.0 / sqrt(pi) * lambda^1.5 * lambda^(0.5 * (E_ALPHA * je + beta_l[ix])) * (r1 + r2 + r3 + r4 + r5 + r6 + r7)

            test[is, js, ie, je, ix] = t / tauinv_ab
            field[is, js, ie, je, ix] = fld / tauinv_ab
        end
    end
    return test, field
end

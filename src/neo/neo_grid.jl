# Energy/pitch-angle basis of NEO (neo_energy_grid.f90): the associated
# Laguerre basis in x = v/vth and Legendre basis in xi, with every energy
# integral done in closed form through Gamma functions of half-integers.
#
# The basis depends only on (n_energy, n_xi, laguerre_method), so it is
# computed once per size and cached.

const E_ALPHA = 2  # neo_energy_grid.f90: e_alpha=2 (should be 1 or 2)

"""
    gamma2(n::Integer) -> Γ(n/2)

Exactly NEO's `gamma2` (`neo_energy_grid.f90`): product recursion from Γ(1)
or Γ(1/2), so it agrees with the Fortran to the last bit.
"""
function gamma2(n::Integer)
    if iseven(n)
        g = 1.0
        i1 = 4
    else
        g = sqrt(pi)
        i1 = 3
    end
    for i in i1:2:n
        g = g * (i / 2.0 - 1.0)
    end
    return g
end

"""
    NEOBasis

Plasma-independent tables of the energy/xi basis (`ENERGY_basis_ints_alloc`,
`ENERGY_basis_ints`). All arrays are 0-based in the energy and xi indices like
the Fortran. `lag[ie,ke,ix]` is the Laguerre expansion coefficient ratio that
NEO recomputes inline as `zarg1` for every matrix element.
"""
struct NEOBasis
    n_energy::Int
    n_xi::Int
    laguerre_method::Int
    e_lag::OffsetVector{Int,Vector{Int}}
    xi_beta_l::OffsetVector{Int,Vector{Int}}
    mygamma2::Vector{Float64}                 # mygamma2[n] = Γ(n/2), 1-based like Fortran
    lag::OffsetArray{Float64,3,Array{Float64,3}}  # (0:ne, 0:ne, 0:nxi) [ie,ke,ix], ke <= ie
    evec_e0::OffsetMatrix{Float64,Matrix{Float64}}
    evec_e1::OffsetMatrix{Float64,Matrix{Float64}}
    evec_e2::OffsetMatrix{Float64,Matrix{Float64}}
    evec_e05::OffsetMatrix{Float64,Matrix{Float64}}
    evec_e105::OffsetMatrix{Float64,Matrix{Float64}}
    emat_e05::OffsetArray{Float64,4,Array{Float64,4}}    # (0:ne,0:ne,0:nxi,1:2)
    emat_en05::OffsetArray{Float64,4,Array{Float64,4}}
    emat_e05de::OffsetArray{Float64,4,Array{Float64,4}}
    emat_e0::OffsetArray{Float64,4,Array{Float64,4}}     # (…,1:1)
    emat_e1::OffsetArray{Float64,4,Array{Float64,4}}     # (…,1:3)
end

const _BASIS_CACHE = Dict{Tuple{Int,Int,Int},NEOBasis}()
const _BASIS_LOCK = ReentrantLock()

"""
    NEOBasis(n_energy, n_xi, laguerre_method=1)

Cached basis tables for the given sizes (see [`NEOBasis`](@ref)).
"""
function NEOBasis(n_energy::Int, n_xi::Int, laguerre_method::Int=1)
    key = (n_energy, n_xi, laguerre_method)
    lock(_BASIS_LOCK) do
        return get!(_BASIS_CACHE, key) do
            return _build_basis(n_energy, n_xi, laguerre_method)
        end
    end
end

NEOBasis(p::NEOParams) = NEOBasis(p.n_energy, p.n_xi, p.laguerre_method)

function _build_basis(ne::Int, nxi::Int, laguerre_method::Int)
    e_lag = OffsetVector(zeros(Int, nxi + 1), 0:nxi)
    xi_beta_l = OffsetVector(zeros(Int, nxi + 1), 0:nxi)
    if laguerre_method == 1
        # Laguerre 1/2+3/2
        e_lag .= 3
        e_lag[0] = 1
        xi_beta_l .= 1
        xi_beta_l[0] = 0
    elseif laguerre_method == 2
        # Sonine
        for ix in 0:nxi
            e_lag[ix] = 2 * ix + 1
            xi_beta_l[ix] = ix
        end
    elseif laguerre_method == 3
        # Laguerre 1/2
        e_lag .= 1
        xi_beta_l .= 0
    elseif laguerre_method == 4
        # Laguerre 3/2
        e_lag .= 3
        xi_beta_l .= 1
    else
        error("NEONative: laguerre_method=$laguerre_method invalid")
    end

    xarg = 4 * E_ALPHA * ne + 4 * nxi + 12
    mygamma2 = [gamma2(n) for n in 1:xarg]

    # zarg1 factor of the Laguerre expansion, exactly as NEO forms it
    lag = OffsetArray(zeros(ne + 1, ne + 1, nxi + 1), 0:ne, 0:ne, 0:nxi)
    for ix in 0:nxi, ie in 0:ne, ke in 0:ie
        lag[ie, ke, ix] = mygamma2[2+2*ie+e_lag[ix]] / mygamma2[2+2*(ie-ke)] / mygamma2[2+2*ke+e_lag[ix]]
    end

    evec_e0 = OffsetArray(zeros(ne + 1, nxi + 1), 0:ne, 0:nxi)
    evec_e1 = similar(evec_e0)
    evec_e2 = similar(evec_e0)
    evec_e05 = similar(evec_e0)
    evec_e105 = similar(evec_e0)
    fill!(evec_e1, 0.0)
    fill!(evec_e2, 0.0)
    fill!(evec_e05, 0.0)
    fill!(evec_e105, 0.0)

    for ie in 0:ne, ix in 0:nxi
        for ke in 0:ie
            zarg0 = (-1.0)^ke
            zarg1 = lag[ie, ke, ix]
            zarg2 = mygamma2[2+2*ke]
            evec_e0[ie, ix] += 0.5 * zarg0 * zarg1 * (mygamma2[E_ALPHA*ke+xi_beta_l[ix]+3] / zarg2)
            evec_e1[ie, ix] += 0.5 * zarg0 * zarg1 * (mygamma2[E_ALPHA*ke+xi_beta_l[ix]+5] / zarg2)
            evec_e2[ie, ix] += 0.5 * zarg0 * zarg1 * (mygamma2[E_ALPHA*ke+xi_beta_l[ix]+7] / zarg2)
            evec_e05[ie, ix] += 0.5 * zarg0 * zarg1 * (mygamma2[E_ALPHA*ke+xi_beta_l[ix]+4] / zarg2)
            evec_e105[ie, ix] += 0.5 * zarg0 * zarg1 * (mygamma2[E_ALPHA*ke+xi_beta_l[ix]+6] / zarg2)
        end
    end

    emat_e05 = OffsetArray(zeros(ne + 1, ne + 1, nxi + 1, 2), 0:ne, 0:ne, 0:nxi, 1:2)
    emat_en05 = OffsetArray(zeros(ne + 1, ne + 1, nxi + 1, 2), 0:ne, 0:ne, 0:nxi, 1:2)
    emat_e05de = OffsetArray(zeros(ne + 1, ne + 1, nxi + 1, 2), 0:ne, 0:ne, 0:nxi, 1:2)
    emat_e0 = OffsetArray(zeros(ne + 1, ne + 1, nxi + 1, 1), 0:ne, 0:ne, 0:nxi, 1:1)
    emat_e1 = OffsetArray(zeros(ne + 1, ne + 1, nxi + 1, 3), 0:ne, 0:ne, 0:nxi, 1:3)

    for ie in 0:ne, je in 0:ne, ix in 0:nxi
        for ke in 0:ie, me in 0:je
            # diagonal in xi
            jx = ix
            zarg0 = (-1.0)^(ke + me)
            zarg1 = lag[ie, ke, ix] * lag[je, me, jx]
            zarg2 = mygamma2[2+2*ke] * mygamma2[2+2*me]
            xarg = E_ALPHA * (ke + me) + xi_beta_l[ix] + xi_beta_l[jx] + 3
            emat_e0[ie, je, ix, 1] += 0.5 * zarg0 * zarg1 * (mygamma2[xarg] / zarg2)

            # xi +/- 1
            for kx in 1:2
                jx = kx == 1 ? ix - 1 : ix + 1
                if 0 <= jx <= nxi
                    zarg1 = lag[ie, ke, ix] * lag[je, me, jx]
                    xarg = E_ALPHA * (ke + me) + xi_beta_l[ix] + xi_beta_l[jx] + 4
                    emat_e05[ie, je, ix, kx] += 0.5 * zarg0 * zarg1 * (mygamma2[xarg] / zarg2)
                    xarg = E_ALPHA * (ke + me) + xi_beta_l[ix] + xi_beta_l[jx] + 2
                    emat_en05[ie, je, ix, kx] += 0.5 * zarg0 * zarg1 * (mygamma2[xarg] / zarg2)
                    emat_e05de[ie, je, ix, kx] += 0.5 * zarg0 * zarg1 * (mygamma2[xarg] / zarg2) * (E_ALPHA * me + xi_beta_l[jx])
                end
            end

            # xi -2, 0, +2
            for kx in 1:3
                jx = kx == 1 ? ix - 2 : (kx == 2 ? ix : ix + 2)
                if 0 <= jx <= nxi
                    zarg1 = lag[ie, ke, ix] * lag[je, me, jx]
                    xarg = E_ALPHA * (ke + me) + xi_beta_l[ix] + xi_beta_l[jx] + 5
                    emat_e1[ie, je, ix, kx] += 0.5 * zarg0 * zarg1 * (mygamma2[xarg] / zarg2)
                end
            end
        end
    end

    return NEOBasis(ne, nxi, laguerre_method, e_lag, xi_beta_l, mygamma2, lag,
        evec_e0, evec_e1, evec_e2, evec_e05, evec_e105,
        emat_e05, emat_en05, emat_e05de, emat_e0, emat_e1)
end

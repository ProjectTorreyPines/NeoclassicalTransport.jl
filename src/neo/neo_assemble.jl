# Assembly of the local drift-kinetic system (neo_do.f90): the row/column
# index map, the sparsity pattern (which depends only on the grid sizes and
# the constraint class and is cached), the matrix values, and the RHS source.
#
# Rows and columns are ordered mindx(is,ie,ix,it), theta fastest. Every row
# holds, in this order: the collision entries (test particle over je, then
# field particle over je and js), the streaming entries (je, theta stencil
# id=-2..2 without 0, jx=ix-1 then ix+1), and the trapping/rotation entries
# (je, jx=ix-1 then ix+1). Constraint rows (<g>=0) replace the kinetic
# equation at ix=0, it=1 for ie in {0,1} (collision models 3-5) or for every
# ie (models 1-2). The triplet layout is fixed per row, so threads can fill
# disjoint slices and the serial and threaded fills are identical.

"""
    NEOPattern

Sparsity pattern of the kinetic matrix for given grid sizes: the index map
`mindx`, the triplet columns `J` with per-row offsets `row_ptr`, the merged
CSC skeleton (`colptr`, `rowval`) and the map `pos` from triplet to CSC slot.
"""
struct NEOPattern
    n_species::Int
    n_energy::Int
    n_xi::Int
    n_theta::Int
    constraint_all_ie::Bool
    n_row::Int
    mindx::OffsetArray{Int,4,Array{Int,4}}   # (1:ns, 0:ne, 0:nxi, 1:nθ)
    row_ptr::Vector{Int}                      # n_row+1
    J::Vector{Int}
    colptr::Vector{Int}
    rowval::Vector{Int}
    pos::Vector{Int}
end

const _PATTERN_CACHE = Dict{NTuple{5,Int},NEOPattern}()
const _PATTERN_LOCK = ReentrantLock()

constraint_all_ie(p::NEOParams) = p.collision_model in (1, 2)

"""
    NEOPattern(p::NEOParams)

Cached pattern for the sizes and constraint class of `p`.
"""
function NEOPattern(p::NEOParams)
    key = (p.n_species, p.n_energy, p.n_xi, p.n_theta, Int(constraint_all_ie(p)))
    lock(_PATTERN_LOCK) do
        return get!(_PATTERN_CACHE, key) do
            return _build_pattern(key...)
        end
    end
end

@inline is_constraint_row(ie::Int, ix::Int, it::Int, all_ie::Bool) = ix == 0 && it == 1 && (all_ie || ie == 0 || ie == 1)

function _build_pattern(ns::Int, ne::Int, nxi::Int, nth::Int, all_ie_int::Int)
    all_ie = all_ie_int == 1
    n = ns * (ne + 1) * (nxi + 1) * nth

    mindx = OffsetArray(zeros(Int, ns, ne + 1, nxi + 1, nth), 1:ns, 0:ne, 0:nxi, 1:nth)
    i = 0
    for is in 1:ns, ie in 0:ne, ix in 0:nxi, it in 1:nth
        i += 1
        mindx[is, ie, ix, it] = i
    end

    # entries per row
    row_ptr = zeros(Int, n + 1)
    row_ptr[1] = 1
    i = 0
    for is in 1:ns, ie in 0:ne, ix in 0:nxi, it in 1:nth
        i += 1
        if is_constraint_row(ie, ix, it, all_ie)
            cnt = nth
        else
            nvalid = (ix - 1 >= 0 ? 1 : 0) + (ix + 1 <= nxi ? 1 : 0)
            cnt = (ne + 1) * (1 + ns) + (ne + 1) * 4 * nvalid + (ne + 1) * nvalid
        end
        row_ptr[i+1] = row_ptr[i] + cnt
    end
    nnz_trip = row_ptr[n+1] - 1

    J = zeros(Int, nnz_trip)
    I = zeros(Int, nnz_trip)
    for is in 1:ns, ie in 0:ne, ix in 0:nxi, it in 1:nth
        i = mindx[is, ie, ix, it]
        k = row_ptr[i] - 1
        if is_constraint_row(ie, ix, it, all_ie)
            for jt in 1:nth
                k += 1
                I[k] = i
                J[k] = mindx[is, ie, 0, jt]
            end
        else
            # collisions: test particle, then field particle
            for je in 0:ne
                k += 1
                I[k] = i
                J[k] = mindx[is, je, ix, it]
            end
            for je in 0:ne, js in 1:ns
                k += 1
                I[k] = i
                J[k] = mindx[js, je, ix, it]
            end
            # streaming
            for je in 0:ne, id in -2:2
                id == 0 && continue
                jt = thcyc(it + id, nth)
                if ix - 1 >= 0
                    k += 1
                    I[k] = i
                    J[k] = mindx[is, je, ix-1, jt]
                end
                if ix + 1 <= nxi
                    k += 1
                    I[k] = i
                    J[k] = mindx[is, je, ix+1, jt]
                end
            end
            # trapping and rotation
            for je in 0:ne
                if ix - 1 >= 0
                    k += 1
                    I[k] = i
                    J[k] = mindx[is, je, ix-1, it]
                end
                if ix + 1 <= nxi
                    k += 1
                    I[k] = i
                    J[k] = mindx[is, je, ix+1, it]
                end
            end
        end
        @assert k == row_ptr[i+1] - 1
    end

    # merged CSC skeleton and triplet -> slot map
    S = sparse(I, J, ones(nnz_trip), n, n)
    colptr = S.colptr
    rowval = S.rowval
    pos = zeros(Int, nnz_trip)
    for k in 1:nnz_trip
        j = J[k]
        lo = colptr[j]
        hi = colptr[j+1] - 1
        idx = searchsortedfirst(view(rowval, lo:hi), I[k]) + lo - 1
        @assert rowval[idx] == I[k]
        pos[k] = idx
    end

    return NEOPattern(ns, ne, nxi, nth, all_ie, n, mindx, row_ptr, J, colptr, rowval, pos)
end

"""
    NEOCoefficients{T}

The per-(species, theta) kinetic-equation coefficients NEO forms in the
assembly loop: `stream`, `trap`, `driftx`, `rotkin`, `driftxrot1..3`.
"""
struct NEOCoefficients{T<:Real}
    stream::Matrix{T}
    trap::Matrix{T}
    driftx::Matrix{T}
    rotkin::Matrix{T}
    driftxrot1::Matrix{T}
    driftxrot2::Matrix{T}
    driftxrot3::Matrix{T}
end

function kinetic_coefficients(p::NEOParams{T}, geo::NEOGeometry, rot::NEORotation) where {T<:Real}
    ns, nth = p.n_species, p.n_theta
    Z, mass, temp, vth, rho, omega_rot = p.z, p.mass, p.temp, p.vth, p.rho, p.omega_rot
    I_div_psip = geo.I_div_psip
    c = NEOCoefficients{T}((zeros(T, ns, nth) for _ in 1:7)...)
    for is in 1:ns, it in 1:nth
        k_par, Bmag, Btor, bigR = geo.k_par[it], geo.Bmag[it], geo.Btor[it], geo.bigR[it]
        # vpar bhat dot grad -> stream * xi * sqrt(ene)
        c.stream[is, it] = sqrt(2.0) * vth[is] * k_par / (12 * geo.d_theta)
        # mu bdot grad B d/dvpar -> trap * (1-xi^2) d/dxi * sqrt(ene)
        c.trap[is, it] = (geo.gradpar_Bmag[it] / Bmag) * sqrt(0.5) * vth[is]
        # vdrift dot grad r -> driftx * (1+xi^2) d/d(r/a) * ene
        c.driftx[is, it] = geo.v_drift_x[it] * mass[is] / (1.0 * Z[is]) * (vth[is])^2
        # rotation
        c.rotkin[is, it] = 0.5 * sqrt(2.0) * vth[is] *
                           (-Z[is] / temp[is] * k_par * rot.phi_rot_deriv[it] + omega_rot^2 * bigR / vth[is]^2 * geo.gradpar_bigR[it])
        c.driftxrot1[is, it] = I_div_psip * mass[is] / (1.0 * Z[is]) * rho / Bmag * (vth[is])^2 *
                               (-Z[is] / temp[is] * k_par * rot.phi_rot_deriv[it] + omega_rot^2 * bigR / vth[is]^2 * geo.gradpar_bigR[it])
        c.driftxrot2[is, it] = I_div_psip / Btor * mass[is] / (1.0 * Z[is]) * rho * vth[is] * 2.0 * sqrt(2.0) * geo.gradpar_bigR[it] * omega_rot
        c.driftxrot3[is, it] = 1.0 / sqrt(2.0) * vth[is]^2 * mass[is] / (1.0 * Z[is]) * rho / Bmag * Btor / (Bmag * I_div_psip) *
                               (2.0 * geo.gradr[it] * geo.gradpar_gradr[it] - geo.gradpar_Bmag[it] / Bmag * geo.gradr[it]^2)
    end
    return c
end

"""
    assemble(p, basis, coll, geo, rot, pattern; serial=false) -> (A, b, coef)

Matrix (as `SparseMatrixCSC` on the cached pattern) and RHS of the kinetic
system, plus the coefficient tables reused by the transport moments. The
triplet values are filled row by row (threaded over rows unless `serial`),
then accumulated into the CSC slots, duplicates summed in a fixed order.
"""
function assemble(p::NEOParams{T}, basis::NEOBasis, coll::NEOCollision, geo::NEOGeometry, rot::NEORotation,
    pattern::NEOPattern; serial::Bool=false) where {T<:Real}
    n = pattern.n_row
    coef = kinetic_coefficients(p, geo, rot)
    V = zeros(T, length(pattern.J))
    b = zeros(T, n)

    if serial || Threads.nthreads() == 1
        for i in 1:n
            _fill_row!(V, b, i, p, basis, coll, geo, rot, coef, pattern)
        end
    else
        Threads.@threads for i in 1:n
            _fill_row!(V, b, i, p, basis, coll, geo, rot, coef, pattern)
        end
    end

    nzval = zeros(T, length(pattern.rowval))
    pos = pattern.pos
    @inbounds for k in eachindex(V)
        nzval[pos[k]] += V[k]
    end
    A = SparseMatrixCSC(n, n, pattern.colptr, pattern.rowval, nzval)
    return A, b, coef
end

# decompose a row index back into (is, ie, ix, it): theta fastest
@inline function row_indices(i::Int, pattern::NEOPattern)
    nth, nxi1, ne1 = pattern.n_theta, pattern.n_xi + 1, pattern.n_energy + 1
    r = i - 1
    it = r % nth + 1
    r ÷= nth
    ix = r % nxi1
    r ÷= nxi1
    ie = r % ne1
    is = r ÷ ne1 + 1
    return is, ie, ix, it
end

function _fill_row!(V, b, i::Int, p::NEOParams{T}, basis::NEOBasis, coll::NEOCollision, geo::NEOGeometry,
    rot::NEORotation, coef::NEOCoefficients, pattern::NEOPattern) where {T<:Real}
    ns, ne, nxi = p.n_species, p.n_energy, p.n_xi
    is, ie, ix, it = row_indices(i, pattern)
    k = pattern.row_ptr[i] - 1

    if is_constraint_row(ie, ix, it, pattern.constraint_all_ie)
        # <f_ie> = 0
        for jt in 1:p.n_theta
            k += 1
            V[k] = geo.w_theta[jt]
        end
        b[i] = zero(T)
        return nothing
    end

    stream = coef.stream[is, it]
    trap = coef.trap[is, it]
    rotkin = coef.rotkin[is, it]
    dens_fac = rot.dens_fac
    test, field = coll.test, coll.field
    emat_e05, emat_e05de, emat_en05 = basis.emat_e05, basis.emat_e05de, basis.emat_en05

    # collisions: test particle
    for je in 0:ne
        k += 1
        a = zero(T)
        for ks in 1:ns
            a = a - test[is, ks, ie, je, ix] * dens_fac[ks, it]
        end
        V[k] = a
    end
    # field particle
    for je in 0:ne, js in 1:ns
        k += 1
        V[k] = -field[is, js, ie, je, ix] * dens_fac[js, it]
    end

    # streaming
    for je in 0:ne, id in -2:2
        id == 0 && continue
        cd = cderiv(id)
        if ix - 1 >= 0
            k += 1
            V[k] = stream * ix / (2 * ix - 1.0) * cd * emat_e05[ie, je, ix, 1]
        end
        if ix + 1 <= nxi
            k += 1
            V[k] = stream * (ix + 1.0) / (2 * ix + 3.0) * cd * emat_e05[ie, je, ix, 2]
        end
    end

    # trapping and rotation
    for je in 0:ne
        if ix - 1 >= 0
            k += 1
            V[k] = trap * ix * (ix - 1.0) / (2 * ix - 1.0) * emat_e05[ie, je, ix, 1] +
                   rotkin * ix / (2 * ix - 1.0) * emat_e05de[ie, je, ix, 1] -
                   rotkin * ix * (ix - 1.0) / (2 * ix - 1.0) * emat_en05[ie, je, ix, 1]
        end
        if ix + 1 <= nxi
            k += 1
            V[k] = -trap * (ix + 1.0) * (ix + 2.0) / (2 * ix + 3.0) * emat_e05[ie, je, ix, 2] +
                   rotkin * (ix + 1.0) / (2 * ix + 3.0) * emat_e05de[ie, je, ix, 2] +
                   rotkin * (ix + 1.0) * (ix + 2.0) / (2 * ix + 3.0) * emat_en05[ie, je, ix, 2]
        end
    end

    b[i] = rhs_source(is, ie, ix, it, p, basis, geo, rot, coef)
    return nothing
end

# set_RHS_source: first-order source -(vdrift dot grad (F0 + Ze/T Phi0)) and E_par
function rhs_source(is::Int, ie::Int, ix::Int, it::Int, p::NEOParams{T}, basis::NEOBasis, geo::NEOGeometry,
    rot::NEORotation, coef::NEOCoefficients) where {T<:Real}
    Z, temp, vth = p.z[is], p.temp[is], p.vth[is]
    omega_rot, omega_rot_deriv = p.omega_rot, p.omega_rot_deriv
    bigR, bigR_th0, Btor, Bmag = geo.bigR[it], geo.bigR_th0, geo.Btor[it], geo.Bmag[it]
    driftx, driftxrot1, driftxrot2, driftxrot3 = coef.driftx[is, it], coef.driftxrot1[is, it], coef.driftxrot2[is, it], coef.driftxrot3[is, it]
    evec_e0, evec_e1, evec_e2, evec_e05, evec_e105 = basis.evec_e0, basis.evec_e1, basis.evec_e2, basis.evec_e05, basis.evec_e105

    # src = (1/F0) dF0/dr + Ze/T dPhi0/dr
    src_F0_Ln = -(p.dlnndr[is] - 1.5 * p.dlntdr[is])   # ene^0 part
    src_F0_Lt = -p.dlntdr[is]                           # ene^1 part
    src_P0 = (1.0 * Z) / temp * p.dphi0dr               # ene^0 Er part

    src_Rot1 = -omega_rot * bigR_th0^2 / vth^2 * omega_rot_deriv -
               p.dlntdr[is] * Z / temp * rot.phi_rot[it] +
               p.dlntdr[is] * (omega_rot / vth)^2 * 0.5 * (bigR^2 - bigR_th0^2) -
               omega_rot^2 * bigR_th0 / vth^2 * geo.bigR_th0_rderiv
    src_Rot2 = omega_rot_deriv * bigR / vth

    if ix == 0
        return -(4.0 / 3.0) * driftx * ((src_F0_Ln + src_P0 + src_Rot1) * evec_e1[ie, ix] + src_F0_Lt * evec_e2[ie, ix]) -
               driftxrot1 * ((src_F0_Ln + src_P0 + src_Rot1) * evec_e0[ie, ix] + src_F0_Lt * evec_e1[ie, ix]) -
               src_Rot2 * 1.0 / 3.0 * driftxrot2 * sqrt(2.0) * Btor / Bmag * evec_e1[ie, ix] -
               src_Rot2 * 4.0 / 3.0 * driftx * omega_rot * bigR / vth * evec_e1[ie, ix] -
               src_Rot2 * driftxrot1 * omega_rot * bigR / vth * evec_e0[ie, ix]
    elseif ix == 1
        return sqrt(2.0) * vth * evec_e05[ie, ix] * (1.0 * Z) / temp * p.epar0 * Bmag / geo.Bmag2_avg -
               driftxrot2 * ((src_F0_Ln + src_P0 + src_Rot1) * evec_e05[ie, ix] + src_F0_Lt * evec_e105[ie, ix]) -
               src_Rot2 * 8.0 / 5.0 * driftx * sqrt(2.0) * Btor / Bmag * evec_e105[ie, ix] -
               src_Rot2 * driftxrot1 * sqrt(2.0) * Btor / Bmag * evec_e05[ie, ix] -
               src_Rot2 * driftxrot2 * omega_rot * bigR / vth * evec_e05[ie, ix] -
               src_Rot2 * 2.0 / 5.0 * driftxrot3 * evec_e105[ie, ix]
    elseif ix == 2
        return -(2.0 / 3.0) * driftx * ((src_F0_Ln + src_P0 + src_Rot1) * evec_e1[ie, ix] + src_F0_Lt * evec_e2[ie, ix]) -
               src_Rot2 * 2.0 / 3.0 * driftxrot2 * sqrt(2.0) * Btor / Bmag * evec_e1[ie, ix] -
               src_Rot2 * 2.0 / 3.0 * driftx * omega_rot * bigR / vth * evec_e1[ie, ix]
    elseif ix == 3
        return -src_Rot2 * 2.0 / 5.0 * driftx * sqrt(2.0) * Btor / Bmag * evec_e105[ie, ix] +
               src_Rot2 * 2.0 / 5.0 * driftxrot3 * evec_e105[ie, ix]
    else
        return zero(T)
    end
end

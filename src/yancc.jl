import Yancc

# Drift kinetic neoclassical fluxes from Yancc.jl, the Julia port of yancc:
# R. Conlin and M. Landreman, Plasma Phys. Control. Fusion 68 (2026) 115001,
# https://doi.org/10.1088/1361-6587/aeb628

"""
    miller_field(; r, a, R0, kappa, delta, zeta=0, shift=0, s_kappa=0, s_delta=0, s_zeta=0, s_zmag=0,
                 q, Bunit, btccw=1, ipccw=1, ntheta=41)

Magnetic field of a Miller local equilibrium (R.L. Miller et al., Phys. Plasmas 5, 973
(1998), with squareness as in NEO) as a `Yancc.Field` on `ntheta` points uniform in
the Miller angle `θ`:

    R = R0 + r cos(θ + asin(δ) sin θ - ζ sin 2θ),    Z = Z0 + κ r sin θ

Lengths in m, `Bunit = (q/r) dψ/dr` in T; the shears are those of NEO
(`s_kappa = r/κ dκ/dr`, `s_delta = r dδ/dr`, `s_zeta = r dζ/dr`, `shift = dR0/dr`,
`s_zmag = dZ0/dr`). `btccw`, `ipccw = ±1` give the directions of the toroidal field
and of the plasma current, positive when counter-clockwise seen from above.

The radial label of the field is `ρ = r/a`: fluxes computed with it are per unit `ρ`
and profile gradients are derivatives in `ρ`.
"""
function miller_field(; r::Real, a::Real, R0::Real, kappa::Real, delta::Real, zeta::Real=0.0, shift::Real=0.0,
    s_kappa::Real=0.0, s_delta::Real=0.0, s_zeta::Real=0.0, s_zmag::Real=0.0, q::Real, Bunit::Real,
    btccw::Real=1, ipccw::Real=1, ntheta::Integer=41)

    theta = range(0, 2π; length=ntheta + 1)[1:ntheta]
    x = asin(delta)
    R, Rr, Jr, dl = zeros(ntheta), zeros(ntheta), zeros(ntheta), zeros(ntheta)
    for (k, t) in enumerate(theta)
        arg = t + x * sin(t) - zeta * sin(2t)
        R[k] = R0 + r * cos(arg)
        R_r = Rr[k] = shift + cos(arg) - sin(arg) * (s_delta / cos(x) * sin(t) - s_zeta * sin(2t))
        R_t = -r * sin(arg) * (1 + x * cos(t) - 2 * zeta * cos(2t))
        Z_r = s_zmag + kappa * (1 + s_kappa) * sin(t)
        Z_t = kappa * r * cos(t)
        Jr[k] = R_r * Z_t - R_t * Z_r   # Jacobian of (r, θ) in the poloidal plane
        dl[k] = hypot(R_t, Z_t)
    end
    @assert all(>(0), Jr) "Miller surface is not nested"

    dpsidr = r * abs(Bunit) / abs(q)             # poloidal flux per radian
    F = 2π * r * abs(Bunit) / (sum(Jr ./ R) * 2π / ntheta)   # from q = F/(2π) ∮ dl/(R² Bp)
    Bp = dpsidr .* (dl ./ Jr) ./ R
    Bt = F ./ R
    Bmag = hypot.(Bp, Bt)

    # (ρ, θ, φ) with φ the cylindrical angle of the right-handed (R, φ, Z): θ runs
    # counter-clockwise in the (R, Z) plane, opposite to the poloidal field of a
    # current along +φ, and the Jacobian is negative.
    st, sz = -sign(ipccw), sign(btccw)
    col(v) = reshape(collect(float.(v)), ntheta, 1)
    # the major radius of the nodes, which a Yancc with toroidal rotation takes
    major_radius = yancc_has_rotation() ? (; R=col(R), dRdrho=col(a .* Rr)) : (;)
    return Yancc.Field(r / a, col(st .* Bp ./ dl), col(sz .* F ./ R .^ 2), col(st .* Bp .* dl), col(sz .* F .+ 0 .* R),
        col(Bmag), col(-a .* R .* Jr); R_major=R0, a_minor=a, iota=st * sz / abs(q), major_radius...)
end

"Whether the installed Yancc.jl solves for plasmas with sonic toroidal rotation (`Yancc.Rotation`)."
yancc_has_rotation() = isdefined(Yancc, :Rotation)

"""
    miller_field(input_neo::InputNEO, a, Bunit; ntheta=41)

Field of the Miller surface described by `input_neo`, with `a` [m] the minor radius
that normalizes its lengths and `Bunit` [T]. The surface is built from the values of
`ForwardDiff.Dual` inputs: it carries no derivatives.
"""
function miller_field(input_neo::InputNEO, a::Real, Bunit::Real; ntheta::Integer=41)
    val(x) = ismissing(x) ? 0.0 : Float64(ForwardDiff.value(x))
    a, Bunit = val(a), val(Bunit)
    return miller_field(; r=val(input_neo.RMIN_OVER_A) * a, a, R0=val(input_neo.RMAJ_OVER_A) * a,
        kappa=val(input_neo.KAPPA), delta=val(input_neo.DELTA), zeta=val(input_neo.ZETA), shift=val(input_neo.SHIFT),
        s_kappa=val(input_neo.S_KAPPA), s_delta=val(input_neo.S_DELTA), s_zeta=val(input_neo.S_ZETA),
        s_zmag=val(input_neo.S_ZMAG), q=val(input_neo.Q), Bunit,
        btccw=ismissing(input_neo.BTCCW) ? 1 : input_neo.BTCCW, ipccw=ismissing(input_neo.IPCCW) ? 1 : input_neo.IPCCW, ntheta)
end

"""
    YanccGeometry

Flux surfaces of one equilibrium for [`yancc`](@ref), valid while the equilibrium is.
"""
struct YanccGeometry
    time::Float64
    gridpoints::Vector{Int}
    fields::Vector{Yancc.Field}
    a::Float64   # minor radius normalizing the radial label ρ = r/a [m]
end

"""
    get_yancc_geometry(eqt, cp1d, gridpoints_cp; ntheta=41)

Miller surfaces at the `cp1d` grid points `gridpoints_cp`, the same that
[`InputNEO`](@ref) describes to NEO.
"""
function get_yancc_geometry(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d,
    gridpoints_cp::AbstractVector{<:Integer}; ntheta::Integer=41)
    eqt1d = eqt.profiles_1d
    rho_cp = cp1d.grid.rho_tor_norm
    a = GACODE.r_min_core_profiles(eqt1d, rho_cp)[end] / IMAS.cgs.m_to_cm
    bunit = IMAS.interp1d(eqt1d.rho_tor_norm, GACODE.bunit(eqt1d))
    fields = [miller_field(InputNEO(eqt, cp1d, i), a, bunit(rho_cp[i]); ntheta) for i in gridpoints_cp]
    return YanccGeometry(eqt.time, collect(Int, gridpoints_cp), fields, Float64(ForwardDiff.value(a)))
end

"""
    yancc_species(cp1d, rmin_over_a, gridpoint_cp)

Maxwellians of the ions (in the order of `cp1d.ion`) and of the electrons at one grid
point, with gradients in `ρ = r/a`, and the indices in `cp1d.ion` of the ions kept
(those with a density). As in `InputNEO` and `hirshmansigmar`, the density of a species
is its total density (thermal + fast) and its temperature is `ion.temperature`. Ion
charges are the average charge states.
"""
function yancc_species(cp1d::IMAS.core_profiles__profiles_1d, rmin_over_a::AbstractVector{<:Real}, gridpoint_cp::Integer)
    i = gridpoint_cp
    mass_unit = IMAS.mks.m_u / Yancc.PROTON_MASS
    ne = cp1d.electrons.density
    S = promote_type(eltype(ne), eltype(cp1d.electrons.temperature), eltype(rmin_over_a))
    plain(x) = Float64(ForwardDiff.value(x))
    function maxwellian(species, T, n)
        dlnT = IMAS.calc_z(rmin_over_a, T, :backward)[i]
        dlnn = IMAS.calc_z(rmin_over_a, n, :backward)[i]
        return Yancc.LocalMaxwellian(species, S(T[i]), S(n[i]), S(T[i] * dlnT), S(n[i] * dlnn))
    end
    species = Yancc.LocalMaxwellian{S}[]
    kept = Int[]
    for (k, ion) in enumerate(cp1d.ion)
        n = ion.density
        n[i] > 1e-10 * ne[i] || continue
        Z = IMAS.avgZ(ion.element[1].z_n, ion.temperature[i])
        push!(species, maxwellian(Yancc.Species(plain(ion.element[1].a) * mass_unit, Z), ion.temperature, n))
        push!(kept, k)
    end
    push!(species, maxwellian(Yancc.Electron, cp1d.electrons.temperature, ne))
    return species, kept
end

"""
    yancc(gridpoint_cp, eqt, cp1d, field; a, nx=6, nalpha=65, p1="4z", rtol=1e-8, tangent_rtol=1e-6,
          max_residual=1e-3, lnlambda=nothing, cache=nothing, rho_s, rmin)

Neoclassical fluxes at one `cp1d` grid point from the drift kinetic equation with the
linearized Fokker-Planck collision operator, solved by Yancc.jl on the flux surface
`field` (see [`get_yancc_geometry`](@ref), whose `a` normalizes the radial label).
Returns a `GACODE.FluxSolution` in gyro-Bohm units; the energy fluxes include
convection and the toroidal stress is not computed.

`nx` speeds, `nalpha` pitch angles and the stencil `p1` set the velocity space
resolution; `lnlambda` is the Coulomb logarithm (one for each pair of species when
`nothing`). A `Yancc.TokamakCache` passed as `cache` carries the preconditioner and
the solution from one call to the next at the same grid point.

The species are those of `InputNEO` and `hirshmansigmar`: every ion of `cp1d` with its
total density (thermal + fast) and its own temperature, and the electrons.

`rtol` is the tolerance of the iterative solve. With a `cache`, the fluxes for small
changes of the profiles (below 1e-5, as in a finite difference Jacobian) are computed
as a correction to the last full solve and vary smoothly; without one they are smooth
only down to about `rtol`.

The solve can stop far above `rtol` where the kinetic equation is ill conditioned, as for
the temperatures of a few eV or of many MeV that a transport solver may try on its way.
Fluxes of such a solution are meaningless: a relative residual above `max_residual` is an
error, which a transport solver can take as a rejected step.

Profiles of `ForwardDiff.Dual` numbers give fluxes with their derivatives: the kinetic
equation is solved in plain numbers and each direction of differentiation by one more
linear solve to `tangent_rtol`, which a `cache` reuses between calls at the same profiles
(the chunks of a Jacobian). The flux surface carries no derivatives.

If you use these fluxes, please cite R. Conlin and M. Landreman, Plasma Phys. Control.
Fusion 68 (2026) 115001.
"""
function yancc(
    gridpoint_cp::Int,
    eqt::IMAS.equilibrium__time_slice,
    cp1d::IMAS.core_profiles__profiles_1d,
    field::Yancc.Field;
    a::Real=field.a_minor,
    nx::Integer=6,
    nalpha::Integer=65,
    p1::AbstractString="4z",
    rtol::Real=1e-8,
    tangent_rtol::Real=1e-6,
    max_residual::Real=1e-3,
    lnlambda::Union{Nothing,Real}=nothing,
    cache::Union{Nothing,Yancc.TokamakCache}=nothing,
    rho_s::AbstractVector{<:Real}=GACODE.rho_s(cp1d, eqt),
    rmin::AbstractVector{<:Real}=GACODE.r_min_core_profiles(eqt.profiles_1d, cp1d.grid.rho_tor_norm)
)
    species, kept = yancc_species(cp1d, rmin ./ (a * IMAS.cgs.m_to_cm), gridpoint_cp)
    op = Yancc.DKEOperator(field, Yancc.UniformPitchAngleGrid(nalpha), Yancc.MaxwellSpeedGrid(nx), species, 0.0; p1, lnlambda)
    sol = cache === nothing ? Yancc.solve_dke_tokamak(op; rtol, tangent_rtol) : Yancc.solve_dke_tokamak!(cache, op; rtol, tangent_rtol)
    if !(sol.residual <= max_residual)
        error("yancc: the drift kinetic solve stopped at the relative residual $(sol.residual) (max_residual = $max_residual) " *
              "at grid point $gridpoint_cp, the equation is ill conditioned for these plasma parameters")
    end
    return yancc_to_GB(sol, kept, eqt, cp1d, gridpoint_cp, a; rho_s)
end

"""
    yancc_to_GB(sol, kept, eqt, cp1d, gridpoint_cp, a; rho_s)

`GACODE.FluxSolution` from a Yancc solution whose species are the ions `kept` of
`cp1d.ion` followed by the electrons, on a surface labelled by `ρ = r/a`: the fluxes
`⟨Γ·∇r⟩ = a ⟨Γ·∇ρ⟩` in the gyro-Bohm units of `GACODE.flux_gacode_to_imas`.
"""
function yancc_to_GB(sol::Yancc.DKESolution, kept::AbstractVector{<:Integer}, eqt::IMAS.equilibrium__time_slice,
    cp1d::IMAS.core_profiles__profiles_1d, gridpoint_cp::Integer, a::Real; rho_s::AbstractVector{<:Real}=GACODE.rho_s(cp1d, eqt))
    ne = cp1d.electrons.density_thermal[gridpoint_cp]
    Te = cp1d.electrons.temperature[gridpoint_cp]
    a_gb = eqt.boundary.minor_radius
    Q_GB = GACODE.gyrobohm_energy_flux(ne, Te, rho_s[gridpoint_cp], a_gb)
    Γ_GB = GACODE.gyrobohm_particle_flux(ne, Te, rho_s[gridpoint_cp], a_gb)
    Qe = a * sol.heat_flux[end] / Q_GB
    particle_flux_i = zeros(typeof(Qe), length(cp1d.ion))
    particle_flux_i[kept] .= a .* sol.particle_flux[1:end-1] ./ Γ_GB
    return GACODE.FluxSolution{typeof(Qe)}(Qe, a * sum(sol.heat_flux[1:end-1]) / Q_GB,
        a * sol.particle_flux[end] / Γ_GB, particle_flux_i, zero(Qe))
end

"Preconditioner and last solution of [`yancc`](@ref) at one grid point, see `Yancc.TokamakCache`."
const YanccCache = Yancc.TokamakCache

"""
    yancc_species(input_neo::InputNEO, ne, Te, a)

Maxwellians of the species of `input_neo`, in its order, for the electron density `ne`
[m⁻³], temperature `Te` [eV] and minor radius `a` [m] that normalize it, and the
Coulomb logarithm implied by its `NU_1`. With [`miller_field`](@ref) this poses to Yancc the problem that
`input_neo` poses to NEO without rotation; [`yancc_operator`](@ref) includes the rotation.
"""
function yancc_species(input_neo::InputNEO, ne::Real, Te::Real, a::Real)
    get(name, k) = Float64(getfield(input_neo, Symbol(name, "_", k)))
    md = IMAS.cgs.md / 1e3   # kg
    species = map(1:input_neo.N_SPECIES) do k
        T, n = get("TEMP", k) * Te, get("DENS", k) * ne
        Yancc.LocalMaxwellian(Yancc.Species(get("MASS", k) * md / Yancc.PROTON_MASS, get("Z", k)), T, n,
            -T * get("DLNTDR", k), -n * get("DLNNDR", k))
    end
    # NU_1 = √2 π n₁ Z₁⁴ e⁴ lnΛ / (√m₁ T₁^{3/2}) a / √(Te/m_D)
    s1 = species[1]
    e2 = Yancc.ELEMENTARY_CHARGE^2 / (4π * Yancc.EPSILON_0)
    T1 = s1.temperature * Yancc.ELEMENTARY_CHARGE
    nu1 = sqrt(2) * π * s1.density * (s1.species.charge / Yancc.ELEMENTARY_CHARGE)^4 * e2^2 / (sqrt(s1.species.mass) * T1^1.5)
    lnlambda = Float64(input_neo.NU_1) * sqrt(Te * Yancc.ELEMENTARY_CHARGE / md) / a / nu1
    return species, lnlambda
end

"""
    yancc_rotation(input_neo::InputNEO, T_norm, a)

Toroidal rotation of `input_neo` as a `Yancc.Rotation`, for the temperature `T_norm` [eV]
and the minor radius `a` [m] that normalize it: the angular velocity
`OMEGA_ROT √(T_norm/m_D)/a` along the `+ζ` of [`miller_field`](@ref) and its derivative in
`ρ = r/a`, `OMEGA_ROT_DERIV √(T_norm/m_D)/a`. The densities and their gradients are those at
the outboard midplane, as in NEO with `ROTATION_MODEL = 2`.

The fluxes of particles and energy depend only on the sign of `OMEGA_ROT OMEGA_ROT_DERIV`,
those of momentum and the flows also on that of `OMEGA_ROT`.
"""
function yancc_rotation(input_neo::InputNEO, T_norm::Real, a::Real)
    yancc_has_rotation() || error("yancc_rotation: the installed Yancc.jl has no toroidal rotation (Yancc.Rotation)")
    md = IMAS.cgs.md / 1e3   # kg
    rate = sqrt(Float64(T_norm) * Yancc.ELEMENTARY_CHARGE / md) / Float64(a)
    val(x) = ismissing(x) ? 0.0 : Float64(ForwardDiff.value(x))
    return Yancc.Rotation(val(input_neo.OMEGA_ROT) * rate, val(input_neo.OMEGA_ROT_DERIV) * rate)
end

"""
    yancc_operator(input_neo::InputNEO, n_norm, T_norm, a; rotation=true, ntheta=41, nalpha=65, nx=10, p1="4z")

`Yancc.DKEOperator` of the local problem that `input_neo` poses to NEO, for the density
`n_norm` [m⁻³], the temperature `T_norm` [eV] and the minor radius `a` [m] that normalize it
(for an `InputNEO` from IMAS data those of the electrons): the Miller surface, with
`Bunit` from `RHO_STAR`, the species, the Coulomb logarithm from `NU_1` and, with
`rotation`, the toroidal rotation of `ROTATION_MODEL = 2` ([`yancc_rotation`](@ref)).
There is no radial electric field besides that of the rotation (`DPHI0DR = 0`).

The solution depends on the normalizing values only through the dimensionless inputs: any
`n_norm`, `T_norm` and `a` give the same [`yancc_neo_outputs`](@ref).
"""
function yancc_operator(input_neo::InputNEO, n_norm::Real, T_norm::Real, a::Real; rotation::Bool=true, ntheta::Integer=41,
    nalpha::Integer=65, nx::Integer=10, p1::AbstractString="4z")
    md = IMAS.cgs.md / 1e3   # kg
    vnorm = sqrt(Float64(T_norm) * Yancc.ELEMENTARY_CHARGE / md)
    Bunit = md * vnorm / (Yancc.ELEMENTARY_CHARGE * a * Float64(input_neo.RHO_STAR))
    field = miller_field(input_neo, a, Bunit; ntheta)
    species, lnlambda = yancc_species(input_neo, n_norm, T_norm, a)
    rotation || return Yancc.DKEOperator(field, Yancc.UniformPitchAngleGrid(nalpha), Yancc.MaxwellSpeedGrid(nx), species, 0.0; p1, lnlambda)
    return Yancc.DKEOperator(field, Yancc.UniformPitchAngleGrid(nalpha), Yancc.MaxwellSpeedGrid(nx), species, 0.0; p1, lnlambda,
        rotation=yancc_rotation(input_neo, T_norm, a))
end

"""
    yancc_neo_outputs(op, sol, input_neo::InputNEO, n_norm, T_norm, a)

Results of the Yancc solution `sol` of `op = yancc_operator(input_neo, n_norm, T_norm, a)` in
the conventions and units of NEO's output, one value for each species of `input_neo`
unless noted:

  - `pflux_dke`, `eflux_dke`, `mflux_dke`: fluxes of particles, energy and toroidal angular
    momentum of the drift kinetic solution, as in `out.neo.transport_flux`: in the gyro-Bohm
    units `n_e c_s (ρ_s/a)²`, `n_e c_s T_e (ρ_s/a)²` and `n_e a T_e (ρ_s/a)²` of the electrons
    of `input_neo`, `c_s = √(T_e/m_D)` and `ρ_s` in `Bunit`. The energy flux is that of the laboratory frame.
  - `pflux_gv`, `eflux_gv`, `mflux_gv`: the fluxes that do not come from collisions
    (gyroviscous), zero for an up-down symmetric surface and without rotation.
  - `pflux_tgyro`, `eflux_tgyro`, `mflux_tgyro`: their sums, with the energy flux of the
    rotating frame (that of the laboratory frame less `ω` times the momentum flux).
  - `uparB`: `⟨u∥ B⟩ / (v_norm Bunit)`, `v_norm = √(T_norm/m_D)`.
  - `vpol`: poloidal flow at the outboard midplane over `v_norm`, positive along the Miller
    angle `θ`.
  - `jpar` (one number): `⟨j∥ B⟩ / (e n_norm v_norm Bunit)`.

Flows and current are of first order in `ρ_s/a`, with the value of `RHO_STAR`, as in
`out.neo.transport`.

Checked against Fortran NEO (`ROTATION_MODEL = 2`, collision model 4) for a counter-clockwise
current and a clockwise toroidal field, on up-down symmetric surfaces; the gyroviscous
fluxes are not.
"""
function yancc_neo_outputs(op::Yancc.DKEOperator, sol::Yancc.DKESolution, input_neo::InputNEO, n_norm::Real, T_norm::Real, a::Real)
    e = Yancc.ELEMENTARY_CHARGE
    md = IMAS.cgs.md / 1e3   # kg
    get(name, k) = Float64(getfield(input_neo, Symbol(name, "_", k)))
    ns = input_neo.N_SPECIES
    n_norm, T_norm, a = Float64(n_norm), Float64(T_norm), Float64(a)
    vnorm = sqrt(T_norm * e / md)
    rho_star = Float64(input_neo.RHO_STAR)
    Bunit = md * vnorm / (e * a * rho_star)

    # gyro-Bohm units of the electrons; Yancc's fluxes are per unit ρ = r/a
    ie = findfirst(k -> get("Z", k) < 0, 1:ns)
    ie === nothing && error("yancc_neo_outputs: input_neo has no electrons to define the gyro-Bohm units")
    ne, Te = get("DENS", ie), get("TEMP", ie)
    Γ_GB = n_norm * vnorm * rho_star^2 / a * ne * Te^1.5
    Q_GB = Γ_GB * T_norm * e * Te
    Π_GB = n_norm * T_norm * e * rho_star^2 * ne * Te^2   # n_e a T_e (ρ_s/a)², per unit ρ

    rotating = yancc_has_rotation() && op.rotation !== nothing
    G, Q, P = sol.particle_flux, sol.heat_flux, rotating ? sol.angular_momentum_flux : zero(sol.heat_flux)
    if rotating
        H = Yancc.nondissipative_fluxes(op)
        lab(G, Q, P) = Yancc.energy_flux(op, G, Q, P; frame=:laboratory)
        rot(G, Q, P) = Yancc.energy_flux(op, G, Q, P; frame=:rotating)
        GH, PH = H.particle_flux, H.angular_momentum_flux
        E, EH = lab(G, Q, P), lab(GH, H.heat_flux, PH)
        Er = rot(G, Q, P) .+ rot(GH, H.heat_flux, PH)
        omega = op.rotation.omega
    else
        GH, PH, EH = zero(G), zero(P), zero(Q)
        E, Er = Q, Q
        omega = 0.0
    end

    # Poloidal flow at θ = 0 (the first node): the flow is K B/n + ω₁ R ζ̂ with the rotation
    # ω₁ = -(dp/dψ - n m ω² R ∂R/∂ψ)/(Z e n) of the force balance there, so the parallel flow
    # gives K B/n = V∥ - ω₁ R B_ζ/B, whose poloidal part is along B_p. ψ grows outwards for a
    # counter-clockwise current and B_p is then along -θ.
    fd = op.field
    ipccw = ismissing(input_neo.IPCCW) ? 1.0 : Float64(sign(input_neo.IPCCW))
    r = Float64(input_neo.RMIN_OVER_A) * a
    R1 = Float64(input_neo.RMAJ_OVER_A) * a + r
    dR1dr = 1 + (ismissing(input_neo.SHIFT) ? 0.0 : Float64(input_neo.SHIFT))
    B1 = fd.Bmag[1, 1]
    Bp1 = sqrt(max(B1^2 - (fd.G / R1)^2, 0.0))
    dpsidr = r * Bunit / abs(Float64(input_neo.Q))
    vpol = map(1:ns) do k
        T, Z, m = get("TEMP", k) * T_norm, get("Z", k), get("MASS", k) * md
        force = -T * (get("DLNNDR", k) + get("DLNTDR", k)) / a - (m / e) * omega^2 * R1 * dR1dr   # per particle [eV/m]
        w1 = force / (Z * dpsidr)   # -ω₁ for a counter-clockwise current
        (-ipccw * sol.Vpar[1, 1, k] - w1 * fd.G / B1) * Bp1 / B1 / vnorm
    end

    return (;
        pflux_dke=G ./ Γ_GB, eflux_dke=E ./ Q_GB, mflux_dke=P ./ Π_GB,
        pflux_gv=GH ./ Γ_GB, eflux_gv=EH ./ Q_GB, mflux_gv=PH ./ Π_GB,
        pflux_tgyro=(G .+ GH) ./ Γ_GB, eflux_tgyro=Er ./ Q_GB, mflux_tgyro=(P .+ PH) ./ Π_GB,
        uparB=sol.VparB ./ (vnorm * Bunit), vpol, jpar=sol.JparB / (e * n_norm * vnorm * Bunit))
end

#= ================ =#
#  NEO-NN inference
#= ================ =#
# Neural-network surrogates for NEO (Fokker-Planck) neoclassical transport,
# trained on NEO databases spanning d3d, mastu and nstx equilibria.
#
# Inference contract (matches the training pipeline exactly):
# - inputs are RAW linear values in NEO bulk-ion-normalized units; the forward
#   pass applies log10 to any `*_log10` feature (only NU_1_log10), then
#   standardizes (x-xm)/xσ; outputs are de-standardized y*yσ+ym.
# - outputs are PHYSICAL, SIGNED, NEO "tgyro"-block (electron-gyroBohm)
#   normalized values — directly comparable to run_neo's parsed tgyro block.
#   (NEO's tgyro block divides by electron-GB factors regardless of which
#   species normalizes the run, so no output conversion is needed.)
# - the training DB is BULK-ION (deuterium) normalized: DENS_1 = TEMP_1 = 1,
#   v_norm = sqrt(k*T_D/m_D). The InputNEO constructor is ELECTRON normalized
#   (n_norm = ne, t_norm = Te) — InputNEONN performs the conversion.

import Flux
import Dates
import BSON
import Memoize
import Measurements
using Statistics: mean, std

const _log_suffix = "_log10"

#= ============ =#
#  Model structs
#= ============ =#

abstract type NEONNfluxmodel end

"""
    NEONNmodel

One member of a NEO-NN ensemble: a Flux chain plus the training-set
standardization constants and feature/output names. Weights are promoted to
Float64 at load time so ForwardDiff duals and Measurements propagate.
"""
struct NEONNmodel <: NEONNfluxmodel
    fluxmodel::Flux.Chain
    name::String
    date::Dates.DateTime
    xnames::Vector{String}
    ynames::Vector{String}
    xm::Vector{Float64}
    xσ::Vector{Float64}
    ym::Vector{Float64}
    yσ::Vector{Float64}
    xbounds::Matrix{Float64}
    ybounds::Matrix{Float64}
end

"""
    NEONNensemble

Ensemble of [`NEONNmodel`](@ref)s (typically 20). Property access other than
`.models` delegates to the first member (`ens.xnames`, `ens.ynames`, ...).
"""
struct NEONNensemble <: NEONNfluxmodel
    models::Vector{NEONNmodel}
end

function Base.getproperty(ensemble::NEONNensemble, field::Symbol)
    if field === :models
        return getfield(ensemble, field)
    elseif field === :fluxmodel
        error("Cannot access .fluxmodel of a NEONNensemble; use .models[k].fluxmodel")
    else
        return getfield(getfield(ensemble, :models)[1], field)
    end
end

function Base.show(io::IO, ::MIME"text/plain", model::NEONNmodel)
    println(io, "NEONNmodel")
    println(io, "name: $(model.name)")
    println(io, "date: $(model.date)")
    println(io, "n inputs: $(length(model.xnames))")
    return println(io, "n outputs: $(length(model.ynames))")
end

function Base.show(io::IO, mime::MIME"text/plain", ensemble::NEONNensemble)
    println(io, "NEONNensemble")
    println(io, "n models: $(length(ensemble.models))")
    return show(io, mime, ensemble.models[1])
end

#= =========== =#
#  Load models
#= =========== =#

function dict2mod(savedict::AbstractDict)
    return NEONNmodel(
        Flux.fmap(Flux.f64, savedict[:fluxmodel]),
        String(savedict[:name]),
        savedict[:date],
        String.(savedict[:xnames]),
        String.(savedict[:ynames]),
        Float64.(vec(savedict[:xm])),
        Float64.(vec(savedict[:xσ])),
        Float64.(vec(savedict[:ym])),
        Float64.(vec(savedict[:yσ])),
        Float64.(savedict[:xbounds]),
        Float64.(savedict[:ybounds]))
end

function dict2ens(savedict::AbstractDict)
    return NEONNensemble([dict2mod(modict) for modict in values(savedict)])
end

"""
    loadmodel(filename::AbstractString)

Load a NEO-NN model or ensemble by name (resolved via [`resolve_model_path`](@ref),
so bare names like `"neonn_d3d+mastu+nstx_flux"` find the shipped models)
or by absolute path.
"""
function loadmodel(filename::AbstractString)
    fullpath = resolve_model_path(filename; extensions=[".bson"])
    savedict = BSON.load(fullpath, @__MODULE__)
    if typeof(first(keys(savedict))) <: Integer
        return dict2ens(savedict)
    else
        return dict2mod(savedict)
    end
end

Memoize.@memoize function loadmodelonce(filename::String)
    return loadmodel(filename)
end

#= ================================== =#
#  InputNEO → NN features (InputNEONN)
#= ================================== =#

"""
    InputNEONN{T<:Real}

The NEO-NN feature vector as a struct: field names are exactly the network's
xnames (with the `_log10` suffix stripped — the forward pass applies log10
itself, so `NU_1` here is the LINEAR collision frequency).

All values are in the training-DB convention: 3 species (1 = bulk hydrogenic
ion, 2 = one lumped impurity, 3 = electrons), normalized to the BULK ION
(`DENS_1 = TEMP_1 = 1`, `v_norm = sqrt(k*T_1/m_1)`, lengths over a).
Constructed from an [`InputNEO`](@ref) via species lumping + normalization
conversion (the InputNEO IMAS constructor is electron-normalized).
"""
Base.@kwdef struct InputNEONN{T<:Real}
    DELTA::T
    DENS_2::T
    DENS_3::T
    DLNNDR_1::T
    DLNNDR_2::T
    DLNNDR_3::T
    DLNTDR_1::T
    DLNTDR_2::T
    DLNTDR_3::T
    KAPPA::T
    NU_1::T
    OMEGA_ROT::T
    OMEGA_ROT_DERIV::T
    Q::T
    RHO_STAR::T
    RMAJ_OVER_A::T
    RMIN_OVER_A::T
    SHEAR::T
    SHIFT::T
    S_DELTA::T
    S_KAPPA::T
    TEMP_2::T
    TEMP_3::T
    ZETA::T
    ZMAG_OVER_A::T
end

Base.eltype(::Type{InputNEONN{T}}) where {T<:Real} = T
Base.eltype(nn::InputNEONN) = eltype(typeof(nn))

function _req(input_neo::InputNEO, field::Symbol)
    value = getfield(input_neo, field)
    ismissing(value) && error("NEO-NN: InputNEO field $field is missing — construct InputNEO(eqt, cp1d, gridpoint) first")
    return value
end

"""
    InputNEONN(input_neo::InputNEO)

Build the NEO-NN feature vector from an `InputNEO` (as produced by
`InputNEO(eqt, cp1d, gridpoint_cp)`):

1. classify species by charge: `Z < 0` → electrons, `Z == 1` → hydrogenic
   (folded into the bulk via quasineutrality, mirroring how the training
   pipeline treated fast/secondary hydrogenic ions; temperatures of folded
   species are discarded), `Z > 1` → impurities, lumped into one species with
   charge-weighted Z and density-weighted T and gradients;
2. convert from the constructor's electron normalization (n_norm = ne,
   t_norm = Te) to the training DB's bulk-ion normalization
   (DENS_1 = TEMP_1 = 1, v_norm = sqrt(k*T_1/m_1)).

The element type follows the `InputNEO`'s, so `ForwardDiff.Dual`s survive the
lumping and the normalization conversion.

The training databases are D + C plasmas; DT plasmas work through the
hydrogenic fold but are outside the training distribution — enable
`warn_nn_train_bounds` to monitor extrapolation.
"""
function InputNEONN(input_neo::InputNEO{S}) where {S<:Real}
    T = float(S)
    n_species = _req(input_neo, :N_SPECIES)

    # species 1 must be the bulk hydrogenic ion (the training normalizer).
    # Z_i in InputNEO is z_n normalized to ion 1, so Z_1 == 1 always; use the
    # mass (in deuterium units) to catch a non-hydrogenic first ion.
    mass_1 = _req(input_neo, :MASS_1)
    mass_1 < 1.6 || error(
        "NEO-NN: first ion has MASS_1 = $mass_1 (deuterium units) — not hydrogenic. " *
        "The NEO-NN training convention requires species 1 to be the bulk hydrogenic ion.")

    Z = [_req(input_neo, Symbol("Z_$i")) for i in 1:n_species]
    dens = [_req(input_neo, Symbol("DENS_$i")) for i in 1:n_species]
    temp = [_req(input_neo, Symbol("TEMP_$i")) for i in 1:n_species]
    dlnndr = [_req(input_neo, Symbol("DLNNDR_$i")) for i in 1:n_species]
    dlntdr = [_req(input_neo, Symbol("DLNTDR_$i")) for i in 1:n_species]

    ielec = findall(<(0), Z)
    length(ielec) == 1 || error("NEO-NN: expected exactly one electron species (Z < 0), found $(length(ielec))")
    ielec = only(ielec)
    ihydro = [i for i in 1:n_species if Z[i] == 1]                 # includes species 1
    iimp = [i for i in 1:n_species if Z[i] > 1]
    isempty(iimp) && error(
        "NEO-NN: no impurity species (Z > 1) found. The NEO-NN models require a " *
        "3-species plasma (bulk ion + lumped impurity + electrons); a pure " *
        "hydrogenic plasma is outside the training space.")

    # --- lump impurities (electron-normalized quantities; ñ ≡ n/ne so electron DENS = 1)
    n_imp = sum(dens[j] for j in iimp)
    Z_imp = sum(Z[j] * dens[j] for j in iimp) / n_imp                       # charge-weighted
    T_imp = sum(dens[j] * temp[j] for j in iimp) / n_imp                    # density-weighted
    aLT_imp = sum(dens[j] * dlntdr[j] for j in iimp) / n_imp
    aLn_imp = sum(dens[j] * dlnndr[j] for j in iimp) / n_imp

    # --- bulk: quasineutrality Σ Z_i·ñ_i = ñ_e = 1 folds all hydrogenic density into species 1
    Z_b = 1.0
    n_b = (dens[ielec] - Z_imp * n_imp) / Z_b
    n_b > 0 || error("NEO-NN: quasineutrality gives non-positive bulk density (n_b = $n_b); check impurity densities")
    aLn_b = (dens[ielec] * dlnndr[ielec] - Z_imp * n_imp * aLn_imp) / (Z_b * n_b)

    # --- electron→bulk-ion normalization conversion
    # τ = T_1/Te; every v_norm-based quantity rescales by v_norm_e/v_norm_D = 1/sqrt(τ)
    τ = temp[1]
    sqrtτ = sqrt(τ)

    return InputNEONN{T}(;
        DELTA=_req(input_neo, :DELTA),
        DENS_2=n_imp / n_b,
        DENS_3=dens[ielec] / n_b,
        DLNNDR_1=aLn_b,
        DLNNDR_2=aLn_imp,
        DLNNDR_3=dlnndr[ielec],
        DLNTDR_1=dlntdr[1],
        DLNTDR_2=aLT_imp,
        DLNTDR_3=dlntdr[ielec],
        KAPPA=_req(input_neo, :KAPPA),
        NU_1=_req(input_neo, :NU_1) * (n_b / dens[1]) / sqrtτ,   # ν₁₁ ∝ n₁ (folded density), then v_norm rescale
        OMEGA_ROT=_req(input_neo, :OMEGA_ROT) / sqrtτ,
        OMEGA_ROT_DERIV=_req(input_neo, :OMEGA_ROT_DERIV) / sqrtτ,
        # the training DBs carry SIGNED q, negative for all three devices
        # (GACODE helicity convention with the devices' IPCCW/BTCCW); the
        # InputNEO constructor stores abs(q) — map into the training convention
        Q=-abs(_req(input_neo, :Q)),
        RHO_STAR=_req(input_neo, :RHO_STAR) * sqrtτ,             # ρ_norm ∝ sqrt(T_norm)
        RMAJ_OVER_A=_req(input_neo, :RMAJ_OVER_A),
        RMIN_OVER_A=_req(input_neo, :RMIN_OVER_A),
        SHEAR=_req(input_neo, :SHEAR),
        SHIFT=_req(input_neo, :SHIFT),
        S_DELTA=_req(input_neo, :S_DELTA),
        S_KAPPA=_req(input_neo, :S_KAPPA),
        TEMP_2=T_imp / τ,
        TEMP_3=temp[ielec] / τ,
        ZETA=_req(input_neo, :ZETA),
        ZMAG_OVER_A=_req(input_neo, :ZMAG_OVER_A))
end

#= ================= =#
#  Feature extraction
#= ================= =#

const _XNAMES_FIELD_SYMBOLS_CACHE = Dict{UInt64,Any}()

# Map the model's xnames (log10 suffix stripped) to InputNEONN field symbols,
# cached as a Val so extraction unrolls at compile time. Keyed on the MODEL's
# xnames: the 24-feature d3d nets simply never read TEMP_2.
function _get_xnames_without_log10_suffix(model::NEONNfluxmodel)
    key = hash(model.xnames)
    return get!(_XNAMES_FIELD_SYMBOLS_CACHE, key) do
        return Val(Tuple(Symbol(endswith(x, _log_suffix) ? x[1:end-length(_log_suffix)] : x) for x in model.xnames))
    end
end

@generated function _extract_fields!(inputs::AbstractVector, obj, ::Val{symbols}) where {symbols}
    exprs = [:(@inbounds inputs[$i] = getfield(obj, $(QuoteNode(s)))) for (i, s) in enumerate(symbols)]
    return Expr(:block, exprs..., :(return inputs))
end

#= ============ =#
#  Forward pass
#= ============ =#

"""
    flux_array(model::NEONNmodel, x::AbstractMatrix{T}; warn_nn_train_bounds=true)

Forward pass of one ensemble member. `x` is `(n_features, n_samples)` of RAW
linear values ordered like `model.xnames`; returns `(n_outputs, n_samples)` of
physical (tgyro-normalized, signed) values. Generic in `T` so ForwardDiff
duals and Measurements propagate.
"""
function flux_array(model::NEONNmodel, x::AbstractMatrix{T}; warn_nn_train_bounds::Bool=true) where {T<:Real}
    N, M = size(x)
    N == length(model.xnames) || error("NEO-NN: input has $N features, model $(model.name) expects $(length(model.xnames))")

    xx = copyto!(similar(x, float(T)), x)
    for i in 1:N
        if endswith(model.xnames[i], _log_suffix)
            for j in 1:M
                xx[i, j] = log10(xx[i, j])
            end
        end
    end

    # bounds are in the log10-transformed, un-normalized space; check the first
    # sample only to avoid warning spam over a radial batch
    if warn_nn_train_bounds
        for i in 1:N
            val = xx[i, 1]
            if isnan(val) || isinf(val)
                error("NEO-NN: $(model.xnames[i]) = $(x[i, 1]) is not allowed")
            elseif val < model.xbounds[i, 1]
                @warn("NEO-NN extrapolation: $(model.xnames[i]) = $val is below training bound $(model.xbounds[i, 1])")
            elseif val > model.xbounds[i, 2]
                @warn("NEO-NN extrapolation: $(model.xnames[i]) = $val is above training bound $(model.xbounds[i, 2])")
            end
        end
    end

    xx .= (xx .- model.xm) ./ model.xσ
    yy = model.fluxmodel(xx)
    return yy .* model.yσ .+ model.ym
end

function flux_array(model::NEONNmodel, x::AbstractVector{T}; warn_nn_train_bounds::Bool=true) where {T<:Real}
    return vec(flux_array(model, reshape(x, length(x), 1); warn_nn_train_bounds))
end

"""
    flux_array(ensemble::NEONNensemble, x; uncertain=false, warn_nn_train_bounds=true)

Ensemble forward pass: mean over members; with `uncertain=true` returns
`Measurements.measurement.(mean, std)` (corrected sample std over members).
"""
function flux_array(ensemble::NEONNensemble, x::AbstractMatrix{T}; uncertain::Bool=false, warn_nn_train_bounds::Bool=true) where {T<:Real}
    nmodels = length(ensemble.models)
    nouts = length(ensemble.ynames)
    nsamples = size(x, 2)
    # Measurements.measurement(::Dual, ::Dual) recurses without bound
    uncertain && T <: ForwardDiff.Dual && error(
        "NEO-NN: uncertain=true is not differentiable — take the gradient of the " *
        "ensemble mean, or evaluate the uncertainty on the primal (Float64) input")

    all_y = Array{float(T),3}(undef, nouts, nsamples, nmodels)
    for (k, model) in enumerate(ensemble.models)
        # warn only for the first member to avoid nmodels-fold warning spam
        all_y[:, :, k] = flux_array(model, x; warn_nn_train_bounds=warn_nn_train_bounds && k == 1)
    end

    mean_y = dropdims(mean(all_y; dims=3); dims=3)
    if uncertain && nmodels > 1 && !(T <: Measurements.Measurement)
        std_y = dropdims(std(all_y; dims=3, mean=reshape(mean_y, nouts, nsamples, 1), corrected=true); dims=3)
        return Measurements.measurement.(mean_y, std_y)
    else
        return mean_y
    end
end

function flux_array(ensemble::NEONNensemble, x::AbstractVector{T}; uncertain::Bool=false, warn_nn_train_bounds::Bool=true) where {T<:Real}
    return vec(flux_array(ensemble, reshape(x, length(x), 1); uncertain, warn_nn_train_bounds))
end

function (model::NEONNfluxmodel)(x::AbstractArray; kw...)
    return flux_array(model, x; kw...)
end

#= ========== =#
#  Public API
#= ========== =#

const _FLUX_YNAMES = ["OUT_$(t)flux_tgyro_$(s)" for t in ("p", "e", "m") for s in ("ion1", "ion2", "elec")]
const _FLOW_YNAMES = ["OUT_jpar", "OUT_vpol_elec", "OUT_vpol_ion1", "OUT_vpol_ion2"]

function _validate_group(ensemble::NEONNfluxmodel, expected::Vector{String}, group::String, model_filename::String)
    Set(ensemble.ynames) == Set(expected) || error(
        "NEO-NN: model '$model_filename' has outputs $(ensemble.ynames), " *
        "not a $group model. Flux models end in _flux, flow models in _flow.")
    return nothing
end

function _build_X(ensemble::NEONNfluxmodel, input_neos::Vector{<:InputNEO})
    nn_inputs = InputNEONN.(input_neos)
    # element type comes from the inputs (Float64 normally, Dual under AD), so
    # the whole forward pass stays differentiable
    T = isempty(nn_inputs) ? Float64 : promote_type(eltype.(nn_inputs)...)
    X = Matrix{T}(undef, length(ensemble.xnames), length(nn_inputs))
    xnames_val = _get_xnames_without_log10_suffix(ensemble)
    for (i, nn_input) in enumerate(nn_inputs)
        _extract_fields!(view(X, :, i), nn_input, xnames_val)
    end
    return X
end

#= =============== =#
#  Radial blending
#= =============== =#
# The DIII-D models come in radial families trained on region-specific DBs
# (core rho 0.10-0.90, near-edge 0.68-0.94, edge 0.80-0.99 — the same radial
# windows as the TGLF-NN DIII-D DBs). For the family bases below,
# `run_neonn`/`run_neonn_flow` evaluate the named (core) model everywhere and
# then overwrite each radial point in the near-edge / edge region with the
# corresponding region net's prediction — the same hard-switch blending, at
# the same boundaries (0.881 / 0.975 in r/a, here `RMIN_OVER_A`; TGLF-NN's
# `RMIN_LOC` is the same coordinate), as TurbulentTransport.jl uses for the
# `sat3_em_d3d_azf-1_withnegD` TGLF-NN model family.
const _RADIAL_BLEND_NEAREDGE_RMIN = 0.881
const _RADIAL_BLEND_EDGE_RMIN = 0.975
# base (model_filename minus _flux/_flow) =>
#   (near-edge base, edge base, near-edge switch r/a, edge switch r/a)
const _RADIAL_BLEND_VARIANTS = Dict(
    "neonn_d3d" => ("neonn_d3dnearedge", "neonn_d3dedge",
        _RADIAL_BLEND_NEAREDGE_RMIN, _RADIAL_BLEND_EDGE_RMIN),
    # *_withnegD: trained on the positive + negative triangularity DBs jointly
    # (TGLF-NN naming convention)
    "neonn_d3d_withnegD" => ("neonn_d3dnearedge_withnegD", "neonn_d3dedge_withnegD",
        _RADIAL_BLEND_NEAREDGE_RMIN, _RADIAL_BLEND_EDGE_RMIN),
    # spherical-tokamak family (MAST-U + NSTX); the ST DBs span the same radial
    # windows as the DIII-D ones (core 0.10-0.90, near-edge 0.68-0.94, edge
    # 0.80-0.99), so the switch points are shared.
    "neonn_mastu+nstx_withnegD" => (
        "neonn_mastunearedge+nstxnearedge_withnegD", "neonn_mastuedge+nstxedge_withnegD",
        _RADIAL_BLEND_NEAREDGE_RMIN, _RADIAL_BLEND_EDGE_RMIN),
)

# Overwrite the columns of Y whose RMIN_OVER_A falls in the near-edge / edge
# region with the region net's prediction. No-op for models outside
# `_RADIAL_BLEND_VARIANTS`. Mutates and returns Y.
function _radial_blend!(Y::AbstractMatrix, X::AbstractMatrix, ensemble::NEONNfluxmodel,
    model_filename::String, expected_ynames::Vector{String}, group::String;
    uncertain::Bool, warn_nn_train_bounds::Bool)

    m = match(r"^(.*)_(flux|flow)$", model_filename)
    m === nothing && return Y
    variants = get(_RADIAL_BLEND_VARIANTS, m.captures[1], nothing)
    variants === nothing && return Y

    k_rmin = findfirst(isequal("RMIN_OVER_A"), ensemble.xnames)
    if k_rmin === nothing
        @warn "RMIN_OVER_A not found in xnames for radial-dependent model blending"
        return Y
    end

    nearedge_base, edge_base, r_nearedge, r_edge = variants
    rmin = view(X, k_rmin, :)
    nearedge_mask = (rmin .>= r_nearedge) .& (rmin .< r_edge)
    edge_mask = rmin .>= r_edge
    for (mask, base) in ((nearedge_mask, nearedge_base), (edge_mask, edge_base))
        any(mask) || continue
        variant = "$(base)_$(m.captures[2])"
        vens = loadmodelonce(variant)
        _validate_group(vens, expected_ynames, group, variant)
        # X columns and Y rows are reused across the family: schemas must match.
        vens.xnames == ensemble.xnames || error(
            "NEO-NN radial blending: '$variant' xnames differ from '$model_filename'")
        vens.ynames == ensemble.ynames || error(
            "NEO-NN radial blending: '$variant' ynames differ from '$model_filename'")
        Y[:, mask] .= flux_array(vens, X[:, mask]; uncertain, warn_nn_train_bounds)
    end
    return Y
end

"""
    run_neonn(input_neos::Vector{<:InputNEO};
              model_filename="neonn_d3d+mastu+nstx_flux",
              uncertain=false, warn_nn_train_bounds=true)

Evaluate the NEO-NN flux surrogate at each radial point and return a
`Vector{GACODE.FluxSolution}` in gyroBohm units — drop-in comparable to
[`run_neo`](@ref) (same tgyro-block normalization, no NEO executable needed).

The default model is the joint `d3d+mastu+nstx` net (recommended); the
single-device nets (`neonn_d3d_flux`, `neonn_mastu+nstx_flux`)
are selectable by name — see [`available_models`](@ref).

Radial families blend automatically: with `model_filename` set to a family's
core net (`neonn_d3d_flux`, the joint ± triangularity
`neonn_d3d_withnegD_flux`, or the spherical-tokamak
`neonn_mastu+nstx_withnegD_flux`), radial points with
`RMIN_OVER_A >= 0.881` use the family's near-edge net and points with
`RMIN_OVER_A >= 0.975` its edge net — the same region switching
TurbulentTransport.jl applies for the `sat3_em_d3d_azf-1_withnegD` TGLF-NN
model (see `_RADIAL_BLEND_VARIANTS`).

`PARTICLE_FLUX_i` has length 2: `[bulk ion, lumped impurity]` (the NN species
reduction), unlike `run_neo` which returns one entry per plasma ion.
`STRESS_TOR_i` sums the two ion momentum-flux channels (electron momentum flux
is discarded, exactly as `run_neo` does). With `uncertain=true` values are
`Measurements.Measurement`s (ensemble mean ± std).

Note: signed outputs follow the training devices' helicity convention
(IPCCW/BTCCW are not NN inputs).
"""
function run_neonn(input_neos::Vector{<:InputNEO};
    model_filename::String="neonn_d3d+mastu+nstx_flux",
    uncertain::Bool=false,
    warn_nn_train_bounds::Bool=true)

    ensemble = loadmodelonce(model_filename)
    _validate_group(ensemble, _FLUX_YNAMES, "flux", model_filename)

    X = _build_X(ensemble, input_neos)
    Y = flux_array(ensemble, X; uncertain, warn_nn_train_bounds)
    _radial_blend!(Y, X, ensemble, model_filename, _FLUX_YNAMES, "flux"; uncertain, warn_nn_train_bounds)

    iy = Dict(name => k for (k, name) in enumerate(ensemble.ynames))
    T = eltype(Y)
    return [
        GACODE.FluxSolution{T}(
            Y[iy["OUT_eflux_tgyro_elec"], i],
            Y[iy["OUT_eflux_tgyro_ion1"], i] + Y[iy["OUT_eflux_tgyro_ion2"], i],
            Y[iy["OUT_pflux_tgyro_elec"], i],
            T[Y[iy["OUT_pflux_tgyro_ion1"], i], Y[iy["OUT_pflux_tgyro_ion2"], i]],
            Y[iy["OUT_mflux_tgyro_ion1"], i] + Y[iy["OUT_mflux_tgyro_ion2"], i])
        for i in eachindex(input_neos)
    ]
end

function run_neonn(input_neo::InputNEO; kw...)
    return only(run_neonn([input_neo]; kw...))
end

"""
    NEOFlowSolution

Neoclassical flow quantities from the NEO-NN flow surrogate, in NEO
bulk-ion-normalized units:

- `vpol_ion1`, `vpol_ion2`, `vpol_elec`: poloidal velocity at θ=0 in units of
  `v_norm = sqrt(k*T_1/m_1)` (bulk-ion thermal speed, post-lumping)
- `jpar`: `⟨j∥·B⟩/B_unit` in units of `e*n_1*v_norm`

Multiply by the dimensional factors (built from the BULK-ION density and
temperature, not electron values) to obtain SI/CGS quantities.
"""
struct NEOFlowSolution{T<:Real}
    vpol_ion1::T
    vpol_ion2::T
    vpol_elec::T
    jpar::T
end

"""
    run_neonn_flow(input_neos::Vector{<:InputNEO};
                   model_filename="neonn_d3d+mastu+nstx_flow",
                   uncertain=false, warn_nn_train_bounds=true)

Evaluate the NEO-NN flow surrogate (poloidal velocities + parallel current)
at each radial point; returns a `Vector{NEOFlowSolution}`. Same conventions
and options as [`run_neonn`](@ref), including the automatic DIII-D radial
family blending (`neonn_d3d_flow`, `neonn_d3d_withnegD_flow`, `neonn_mastu+nstx_withnegD_flow`).
"""
function run_neonn_flow(input_neos::Vector{<:InputNEO};
    model_filename::String="neonn_d3d+mastu+nstx_flow",
    uncertain::Bool=false,
    warn_nn_train_bounds::Bool=true)

    ensemble = loadmodelonce(model_filename)
    _validate_group(ensemble, _FLOW_YNAMES, "flow", model_filename)

    X = _build_X(ensemble, input_neos)
    Y = flux_array(ensemble, X; uncertain, warn_nn_train_bounds)
    _radial_blend!(Y, X, ensemble, model_filename, _FLOW_YNAMES, "flow"; uncertain, warn_nn_train_bounds)

    iy = Dict(name => k for (k, name) in enumerate(ensemble.ynames))
    T = eltype(Y)
    return [
        NEOFlowSolution{T}(
            Y[iy["OUT_vpol_ion1"], i],
            Y[iy["OUT_vpol_ion2"], i],
            Y[iy["OUT_vpol_elec"], i],
            Y[iy["OUT_jpar"], i])
        for i in eachindex(input_neos)
    ]
end

function run_neonn_flow(input_neo::InputNEO; kw...)
    return only(run_neonn_flow([input_neo]; kw...))
end

#= ================================= =#
#  Radial electric field (force balance)
#= ================================= =#

"""
    NeoclassicalErSolution

Outboard-midplane radial electric field from neoclassical radial force balance,
built on the NEO-NN poloidal flow. All fields SI:

- `Er`: E_R at the outboard midplane [V/m], the sum of the three terms below
- `Er_pressure`: `(dp_s/dR)/(Z_s e n_s)` [V/m] (diamagnetic term)
- `Er_vtor`: `-v_φ,s B_Z` [V/m] (toroidal-rotation term)
- `Er_vpol`: `+v_θ,s B_φ` [V/m] (poloidal-flow term, the NEO-NN contribution)
- `vpol`: species poloidal velocity at θ=0 [m/s] (`+Z` at the outboard midplane)
- `B_tor`, `B_pol`: signed `B_φ = F/R` and `B_Z = sign(Ip)|B_p|` at θ=0 [T]
- `grad_r`: `|∇r|` at θ=0, converting `d/dr` (minor radius) to `d/dR`
- `R_omp`: major radius of the θ=0 point [m]

`Er / grad_r` is the same field expressed against the minor-radius coordinate.
"""
struct NeoclassicalErSolution{T<:Real}
    Er::T
    Er_pressure::T
    Er_vtor::T
    Er_vpol::T
    vpol::T
    B_tor::T
    B_pol::T
    grad_r::T
    R_omp::T
end

"""
    neoclassical_Er(eqt, cp1d, gridpoints, v_tor;
                    species=:impurity, model_filename="neonn_d3d+mastu+nstx_flow",
                    uncertain=false, warn_nn_train_bounds=true)

Close radial force balance for one species with the NEO-NN poloidal flow and
return `Vector{NeoclassicalErSolution}`, one per entry of `gridpoints`
(indices into `cp1d.grid.rho_tor_norm`).

Working in cylindrical `(R, φ, Z)` at the outboard midplane, where `R̂` is the
radial and `Ẑ` the poloidal direction, the `R` component of
`Z_s e n_s (E + v_s × B) = ∇p_s` gives

    E_R = (dp_s/dR)/(Z_s e n_s) - v_φ,s B_Z + v_θ,s B_φ

`v_tor` is the **toroidal velocity of that species** in m/s, signed along the
IMAS toroidal angle (a scalar, or one value per gridpoint). It is a required
input because it is not neoclassical: it is set by momentum transport and
torque, and in the standard diagnostic application it comes from CER on the
same impurity whose `v_θ` the network supplies.

`species` selects the force-balance species and the matching network output,
in the nets' 3-species reduction:
`:bulk` (hydrogenic, `vpol_ion1`), `:impurity` (lumped, `vpol_ion2`) or
`:electron` (`vpol_elec`). Densities, temperatures, charge and gradients are
taken from the same [`InputNEONN`](@ref) the network is evaluated on, so the
lumping is identical to the one behind the flow prediction.

Requires the 2D equilibrium: `B_p` at θ=0 is evaluated from the ψ map, since
the 1D `r B_unit / q` shortcut misses `|∇r|` at θ=0 (1.5-2.4x on a shaped
equilibrium).

Differentiable w.r.t. the profiles: a `cp1d` carrying `ForwardDiff.Dual`s
propagates through the network and the force balance (`Er` comes back as a
Dual). `eqt` must stay Float64 — the 2D ψ interpolant is not Dual-capable —
which is the frozen-geometry setting FUSE's AD path uses anyway.

CAVEATS
- The nets carry the training devices' helicity convention (IPCCW/BTCCW are not
  inputs, `q` is forced negative), so the sign of `vpol` — and hence of
  `Er_vpol` — follows that convention. Validate against a known case before
  trusting the sign on a discharge of opposite helicity; the terms are returned
  separately so `Er_vpol` can be negated.
- `OMEGA_ROT`/`OMEGA_ROT_DERIV` ARE network inputs (taken from
  `cp1d.rotation_frequency_tor_sonic`), so if you derive them from an assumed
  E_r you should iterate this closure to a fixed point.

# Example
```julia
gps = [argmin(abs.(cp1d.grid.rho_tor_norm .- r)) for r in (0.90, 0.95, 0.99)]
v_tor = [8.0e4, 6.0e4, 3.0e4]  # m/s, measured impurity toroidal rotation
sols = NeoclassicalTransport.neoclassical_Er(eqt, cp1d, gps, v_tor; species=:impurity)
```
"""
function neoclassical_Er(
    eqt::IMAS.equilibrium__time_slice,
    cp1d::IMAS.core_profiles__profiles_1d,
    gridpoints::AbstractVector{<:Integer},
    v_tor;
    species::Symbol=:impurity,
    model_filename::String="neonn_d3d+mastu+nstx_flow",
    uncertain::Bool=false,
    warn_nn_train_bounds::Bool=true)

    species in (:bulk, :impurity, :electron) ||
        error("neoclassical_Er: species must be :bulk, :impurity or :electron, got :$species")
    isempty(eqt.profiles_2d) &&
        error("neoclassical_Er: needs the 2D equilibrium (eqt.profiles_2d) to evaluate B_p at θ=0")

    v_tors = v_tor isa Real ? fill(float(v_tor), length(gridpoints)) : collect(float.(v_tor))
    length(v_tors) == length(gridpoints) ||
        error("neoclassical_Er: v_tor has $(length(v_tors)) entries but there are $(length(gridpoints)) gridpoints")

    input_neos = [InputNEO(eqt, cp1d, gp) for gp in gridpoints]
    nns = [InputNEONN(ineo) for ineo in input_neos]
    flows = run_neonn_flow(input_neos; model_filename, uncertain, warn_nn_train_bounds)

    eqt1d = eqt.profiles_1d
    rho = cp1d.grid.rho_tor_norm
    rmin = GACODE.r_min_core_profiles(eqt1d, rho) ./ IMAS.cgs.m_to_cm    # cm -> m
    a = rmin[end]
    R_omp = IMAS.interp1d(eqt1d.rho_tor_norm, eqt1d.r_outboard).(rho)     # m
    Fpol = IMAS.interp1d(eqt1d.rho_tor_norm, eqt1d.f).(rho)               # T*m, signed
    qprof = IMAS.interp1d(eqt1d.rho_tor_norm, eqt1d.q).(rho)
    bunit = IMAS.interp1d(eqt1d.rho_tor_norm, GACODE.bunit(eqt1d)).(rho)  # T
    _, _, PSI_interpolant = IMAS.ψ_interpolant(eqt.profiles_2d)
    Zaxis = eqt.global_quantities.magnetic_axis.z
    sIp = sign(eqt.global_quantities.ip)

    # bulk-ion normalizing speed: sqrt(k*T_1/m_D), the training convention that
    # InputNEONN's electron->ion rescale (1/sqrt(τ), deuterium mass) implies
    k = IMAS.cgs.k          # erg/eV
    md = IMAS.cgs.md        # g
    Te = cp1d.electrons.temperature   # eV
    ne = cp1d.electrons.density       # m^-3

    return map(enumerate(gridpoints)) do (i, gp)
        nn = nns[i]
        T1 = Te[gp] / nn.TEMP_3                       # eV, bulk hydrogenic ion
        v_norm = sqrt(k * T1 / md) / IMAS.cgs.m_to_cm # cm/s -> m/s
        n_b = ne[gp] / nn.DENS_3                      # m^-3, folded bulk ion
        n_imp = nn.DENS_2 * n_b                       # m^-3, lumped impurity
        Z_imp = (ne[gp] - n_b) / n_imp                # quasineutrality, Z_bulk = 1

        if species === :bulk
            vpol_n, Zs, Ts, aLn, aLT = flows[i].vpol_ion1, one(Z_imp), T1, nn.DLNNDR_1, nn.DLNTDR_1
        elseif species === :impurity
            vpol_n, Zs, Ts, aLn, aLT = flows[i].vpol_ion2, Z_imp, nn.TEMP_2 * T1, nn.DLNNDR_2, nn.DLNTDR_2
        else
            vpol_n, Zs, Ts, aLn, aLT = flows[i].vpol_elec, -one(Z_imp), Te[gp], nn.DLNNDR_3, nn.DLNTDR_3
        end

        # the ψ map is Float64-only; on a Dual cp1d the grid-interpolated R_omp is a
        # Dual with zero partials (frozen geometry), so evaluating at its primal
        # loses nothing
        Bp_mag = IMAS.Bp(PSI_interpolant, ForwardDiff.value(R_omp[gp]), Zaxis)   # |B_p| at θ=0 [T]
        B_pol = sIp * Bp_mag                                  # +Z at θ=0 for Ip along +φ
        B_tor = Fpol[gp] / R_omp[gp]
        # |∇r| at θ=0 = (∂ψ/∂R)/(dψ/dr) with dψ/dr = r*B_unit/|q| (GACODE)
        # |∇r| is a magnitude: bunit is signed (negative for B0 < 0, e.g. DIII-D)
        grad_r = rmin[gp] > 0 ? R_omp[gp] * Bp_mag * abs(qprof[gp]) / (rmin[gp] * abs(bunit[gp])) : one(Bp_mag)

        # DLNNDR = -a dln n/dr, DLNTDR = -a dln T/dr; T in eV makes this V/m
        Er_pressure = -(Ts / Zs) * (aLn + aLT) / a * grad_r
        vpol = vpol_n * v_norm
        Er_vtor = -v_tors[i] * B_pol
        Er_vpol = vpol * B_tor

        # promote so `uncertain=true` (Measurements in the flow terms only) still
        # gives a concretely-typed solution
        NeoclassicalErSolution(promote(
            Er_pressure + Er_vtor + Er_vpol,
            Er_pressure, Er_vtor, Er_vpol,
            vpol, B_tor, B_pol, grad_r, R_omp[gp])...)
    end
end

function neoclassical_Er(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d,
    gridpoint::Integer, v_tor::Real; kw...)
    return only(neoclassical_Er(eqt, cp1d, [gridpoint], v_tor; kw...))
end

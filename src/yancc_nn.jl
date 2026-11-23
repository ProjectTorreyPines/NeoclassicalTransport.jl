#= ================== =#
#  Yancc-NN surrogates
#= ================== =#
# Neural-network surrogates of the drift kinetic solver (`yancc`) with sonic toroidal
# rotation. Each model is an ensemble of residual multilayer perceptrons, stored in
# `models/` as `<name>.toml` (names, normalisation, training bounds, architecture) and
# `<name>.bin` (little-endian Float32 weights). The forward pass is written out here, so
# no machine-learning package is needed, and it is generic in the element type, so
# `ForwardDiff.Dual`s propagate.
#
# Inputs and outputs are those of a 3-species local problem (1 = bulk hydrogenic ion,
# 2 = one lumped impurity, 3 = electrons) normalised to the bulk ion; the outputs follow
# the normalisation of NEO's transport output (fluxes in electron gyro-Bohm units).

import TOML
import Statistics
import Measurements
import SHA
import Downloads

struct YanccNNDense
    W::Matrix{Float64}
    b::Vector{Float64}
end

"""
    YanccNNEnsemble

Ensemble of residual multilayer perceptrons with shared normalisation: `members[k]` is the
list of dense layers of member `k` in forward order (input layer, 4 per residual block,
output layer). `xnames` ending in `_log10` are passed to the network as `log10` of the
input; `xmin`/`xmax` are the bounds of the training set in that space.
"""
struct YanccNNEnsemble
    name::String
    xnames::Vector{String}
    ynames::Vector{String}
    xm::Vector{Float64}
    xsigma::Vector{Float64}
    ym::Vector{Float64}
    ysigma::Vector{Float64}
    xmin::Vector{Float64}
    xmax::Vector{Float64}
    nblocks::Int
    activation::Symbol
    members::Vector{Vector{YanccNNDense}}
end

function Base.show(io::IO, ::MIME"text/plain", ens::YanccNNEnsemble)
    print(io, "YanccNNEnsemble $(ens.name): $(length(ens.members)) members, $(length(ens.xnames)) inputs, ",
        "$(length(ens.ynames)) outputs ($(join(ens.ynames, ", ")))")
end

_yanccnn_dir() = joinpath(_pkg_root(), "models")

#= ============================= =#
#  Weights stored with Git LFS
#= ============================= =#
# The `.bin` files are tracked with Git LFS. A package installation does not run
# `git lfs pull`, so the file may be a pointer; the weights are then fetched once, in place
# or, when the installation is read-only, into a scratch space of the package version.

const _LFS_POINTER_PREFIX = "version https://git-lfs.github.com/spec/v1"
const _LFS_REPO = "ProjectTorreyPines/NeoclassicalTransport.jl"
const _MODEL_DOWNLOAD_LOCK = ReentrantLock()

"""
    is_lfs_pointer(path) -> Bool

Whether `path` is a Git LFS pointer instead of the file itself.
"""
function is_lfs_pointer(path::AbstractString)
    isfile(path) && sizeof(_LFS_POINTER_PREFIX) <= filesize(path) <= 500 || return false
    return String(read(path, sizeof(_LFS_POINTER_PREFIX))) == _LFS_POINTER_PREFIX
end

function _lfs_pointer_info(path::AbstractString)
    content = read(path, String)
    oid = match(r"oid sha256:([0-9a-fA-F]{64})", content)
    size = match(r"size (\d+)", content)
    (oid === nothing || size === nothing) && error("Malformed Git LFS pointer at '$path'")
    return (oid=lowercase(oid.captures[1]), size=parse(Int, size.captures[1]))
end

_sha256_of_file(path::AbstractString) = open(io -> bytes2hex(SHA.sha256(io)), path)

# Git references to try on media.githubusercontent.com; the SHA-256 of the pointer decides.
function _lfs_candidate_refs()
    refs = String[]
    ref = strip(get(ENV, "NEOCLASSICALTRANSPORT_MODELS_REF", ""))
    isempty(ref) || push!(refs, ref)
    root = _pkg_root()
    if ispath(joinpath(root, ".git"))
        try
            push!(refs, readchomp(`git -C $root rev-parse HEAD`))
        catch
        end
    end
    ver = pkgversion(@__MODULE__)
    ver === nothing || push!(refs, "v$ver")
    push!(refs, "master")
    return unique!(refs)
end

"""
    ensure_model_file!(path) -> path of the file to read

If `path` is a Git LFS pointer, fetch the file it stands for and return where it is: `path`
itself when the installation is writable, a scratch space otherwise. A download is accepted
only if its SHA-256 is the one of the pointer. Returns `path` unchanged for an ordinary file.
"""
function ensure_model_file!(path::AbstractString)
    is_lfs_pointer(path) || return String(path)
    return lock(_MODEL_DOWNLOAD_LOCK) do
        is_lfs_pointer(path) || return String(path)
        spec = _lfs_pointer_info(path)
        root = _pkg_root()
        rel = relpath(path, root)
        in_place = _dir_writable(dirname(path))
        if in_place && ispath(joinpath(root, ".git"))
            # objects already in the local LFS store first, then the remote
            for cmd in (`git -C $root lfs checkout -- $rel`, `git -C $root lfs pull --include $rel`)
                try
                    run(pipeline(setenv(cmd; dir=root); stdout=devnull, stderr=devnull, stdin=devnull))
                    !is_lfs_pointer(path) && _sha256_of_file(path) == spec.oid && return String(path)
                catch
                end
            end
        end
        target = String(path)
        if !in_place
            ver = pkgversion(@__MODULE__)
            target = joinpath(Scratch.get_scratch!(@__MODULE__, "models-v$(ver === nothing ? "dev" : ver)"), basename(path))
            isfile(target) && filesize(target) == spec.size && _sha256_of_file(target) == spec.oid && return target
        end
        tmp = target * ".download"
        tried = String[]
        for ref in _lfs_candidate_refs()
            push!(tried, ref)
            try
                Downloads.download("https://media.githubusercontent.com/media/$(_LFS_REPO)/$(ref)/$(rel)", tmp)
                if filesize(tmp) == spec.size && _sha256_of_file(tmp) == spec.oid
                    mv(tmp, target; force=true)
                    return target
                end
            catch
            end
        end
        rm(tmp; force=true)
        error("Could not fetch the weights for '$path' (Git LFS pointer; tried the references $(join(tried, ", "))). " *
              "In a git checkout run `git lfs pull`; otherwise set NEOCLASSICALTRANSPORT_MODELS_REF to a commit that has the file.")
    end
end

"""
    yanccnn_models()

Names of the Yancc-NN models shipped with the package.
"""
yanccnn_models() = sort!([f[1:end-5] for f in readdir(_yanccnn_dir()) if startswith(f, "yanccnn_") && endswith(f, ".toml")])

"""
    load_yanccnn(name)

Read the Yancc-NN model `name` (one of [`yanccnn_models`](@ref), or the path of a model's
files without extension).
"""
function load_yanccnn(name::AbstractString)
    stem = isfile(name * ".toml") ? String(name) : joinpath(_yanccnn_dir(), name)
    isfile(stem * ".toml") || error("Yancc-NN model '$name' not found. Available: $(join(yanccnn_models(), ", "))")
    meta = TOML.parsefile(stem * ".toml")
    nin, nout, channels, nblocks, nmodels = meta["nin"], meta["nout"], meta["channels"], meta["nblocks"], meta["nmodels"]
    shapes = vcat([(channels, nin)], fill((channels, channels), 4 * nblocks), [(nout, channels)])
    nvalues = nmodels * sum(o * i + o for (o, i) in shapes)
    weights = ensure_model_file!(stem * ".bin")
    filesize(weights) == 4 * nvalues ||
        error("Yancc-NN model '$name': $weights has $(filesize(weights)) bytes, expected $(4 * nvalues)")
    members = open(weights, "r") do io
        return map(1:nmodels) do _
            return map(shapes) do (o, i)
                W = reshape(Float64.(ltoh.(read!(io, Vector{Float32}(undef, o * i)))), o, i)
                b = Float64.(ltoh.(read!(io, Vector{Float32}(undef, o))))
                return YanccNNDense(W, b)
            end
        end
    end
    activation = Symbol(meta["activation"])
    activation in (:gelu_tanh, :elu) || error("Yancc-NN model '$name': activation $activation is not supported")
    return YanccNNEnsemble(meta["name"], meta["xnames"], [replace(y, "OUT_" => "") for y in meta["ynames"]],
        meta["xm"], meta["xsigma"], meta["ym"], meta["ysigma"], meta["xmin"], meta["xmax"], nblocks, activation, members)
end

const _YANCCNN_CACHE = Dict{String,YanccNNEnsemble}()
const _YANCCNN_LOCK = ReentrantLock()

"[`load_yanccnn`](@ref), reading each model once."
function load_yanccnn_once(name::AbstractString)
    return lock(_YANCCNN_LOCK) do
        return get!(() -> load_yanccnn(name), _YANCCNN_CACHE, String(name))
    end
end

_yanccnn_gelu_tanh(x) = x / 2 * (1 + tanh(sqrt(2 / pi) * (x + 0.044715 * x^3)))
_yanccnn_elu(x) = x > 0 ? x : expm1(x)

function _yanccnn_member(layers::Vector{YanccNNDense}, act::F, nblocks::Integer, x::AbstractVector) where {F}
    h = act.(layers[1].W * x .+ layers[1].b)
    k = 1
    for _ in 1:nblocks
        z = h
        for j in 1:3
            z = act.(layers[k+j].W * z .+ layers[k+j].b)
        end
        z = layers[k+4].W * z .+ layers[k+4].b
        h = act.(z .+ h)
        k += 4
    end
    return layers[end].W * h .+ layers[end].b
end

"""
    yanccnn_array(ens::YanccNNEnsemble, x::AbstractVector; warn_nn_train_bounds=true)

Outputs of every member at one point: a `(n outputs, n members)` matrix. `x` holds the
inputs in the order of `ens.xnames`, as plain values (the `log10` of the `_log10` inputs is
taken here).
"""
function yanccnn_array(ens::YanccNNEnsemble, x::AbstractVector{T}; warn_nn_train_bounds::Bool=true) where {T<:Real}
    length(x) == length(ens.xnames) || error("Yancc-NN $(ens.name): $(length(x)) inputs given, $(length(ens.xnames)) expected")
    xx = [endswith(ens.xnames[i], "_log10") ? log10(x[i]) : float(x[i]) for i in eachindex(x)]
    for i in eachindex(xx)
        isfinite(xx[i]) || error("Yancc-NN $(ens.name): $(ens.xnames[i]) = $(x[i]) is not allowed")
        if warn_nn_train_bounds && !(ens.xmin[i] <= xx[i] <= ens.xmax[i])
            @warn "Yancc-NN $(ens.name) extrapolation: $(ens.xnames[i]) = $(ForwardDiff.value(xx[i])) is outside the training range [$(ens.xmin[i]), $(ens.xmax[i])]"
        end
    end
    xn = (xx .- ens.xm) ./ ens.xsigma
    act = ens.activation === :gelu_tanh ? _yanccnn_gelu_tanh : _yanccnn_elu
    y = Matrix{eltype(xn)}(undef, length(ens.ynames), length(ens.members))
    for (k, layers) in enumerate(ens.members)
        y[:, k] = _yanccnn_member(layers, act, ens.nblocks, xn) .* ens.ysigma .+ ens.ym
    end
    return y
end

#= ============================== =#
#  InputNEO → inputs of the networks
#= ============================== =#

"""
    yanccnn_inputs(input_neo::InputNEO)

Inputs of the Yancc-NN models, as a `NamedTuple`, from an `InputNEO` built by
`InputNEO(eqt, cp1d, gridpoint_cp)` (electron normalised, any number of ions):

 1. the hydrogenic ions (`Z == 1`) are folded into one bulk ion whose density follows from
    quasi-neutrality, and the impurities are lumped into one species with charge-weighted
    `Z` and density-weighted temperature and gradients;
 2. densities, temperatures, `NU_1`, `RHO_STAR` and the rotation are renormalised from the
    electrons to the bulk ion (`DENS_1 = TEMP_1 = 1`, `v_norm = sqrt(T_1/m_D)`).

`Q` is given the sign of the training data (negative).
"""
function yanccnn_inputs(input_neo::InputNEO)
    req(field) = (v = getfield(input_neo, field); ismissing(v) ? error("Yancc-NN: InputNEO field $field is missing") : v)
    n_species = req(:N_SPECIES)
    req(:MASS_1) < 1.6 || error("Yancc-NN: the first ion must be the bulk hydrogenic ion (MASS_1 = $(req(:MASS_1)) in deuterium masses)")
    Z = [req(Symbol("Z_$i")) for i in 1:n_species]
    dens = [req(Symbol("DENS_$i")) for i in 1:n_species]
    temp = [req(Symbol("TEMP_$i")) for i in 1:n_species]
    dlnndr = [req(Symbol("DLNNDR_$i")) for i in 1:n_species]
    dlntdr = [req(Symbol("DLNTDR_$i")) for i in 1:n_species]

    ielec = findall(<(0), Z)
    length(ielec) == 1 || error("Yancc-NN: expected one electron species (Z < 0), found $(length(ielec))")
    ie = only(ielec)
    iimp = [i for i in 1:n_species if Z[i] > 1]
    isempty(iimp) && error("Yancc-NN: no impurity (Z > 1): the models are for a bulk ion, one lumped impurity and electrons")

    n_imp = sum(dens[j] for j in iimp)
    Z_imp = sum(Z[j] * dens[j] for j in iimp) / n_imp
    T_imp = sum(dens[j] * temp[j] for j in iimp) / n_imp
    aLT_imp = sum(dens[j] * dlntdr[j] for j in iimp) / n_imp
    aLn_imp = sum(dens[j] * dlnndr[j] for j in iimp) / n_imp

    n_b = dens[ie] - Z_imp * n_imp
    n_b > 0 || error("Yancc-NN: quasi-neutrality gives a bulk ion density of $n_b")
    aLn_b = (dens[ie] * dlnndr[ie] - Z_imp * n_imp * aLn_imp) / n_b

    tau = temp[1]                   # T_1/T_e: speeds normalised to the bulk ion are 1/sqrt(tau) of the electron-normalised ones
    sqrt_tau = sqrt(tau)
    return (
        DELTA=req(:DELTA), DENS_2=n_imp / n_b, DENS_3=dens[ie] / n_b,
        DLNNDR_1=aLn_b, DLNNDR_2=aLn_imp, DLNNDR_3=dlnndr[ie],
        DLNTDR_1=dlntdr[1], DLNTDR_2=aLT_imp, DLNTDR_3=dlntdr[ie],
        KAPPA=req(:KAPPA), NU_1=req(:NU_1) * (n_b / dens[1]) / sqrt_tau,
        OMEGA_ROT=req(:OMEGA_ROT) / sqrt_tau, OMEGA_ROT_DERIV=req(:OMEGA_ROT_DERIV) / sqrt_tau,
        Q=-abs(req(:Q)), RHO_STAR=req(:RHO_STAR) * sqrt_tau,
        RMAJ_OVER_A=req(:RMAJ_OVER_A), RMIN_OVER_A=req(:RMIN_OVER_A), SHEAR=req(:SHEAR), SHIFT=req(:SHIFT),
        S_DELTA=req(:S_DELTA), S_KAPPA=req(:S_KAPPA), TEMP_2=T_imp / tau, TEMP_3=temp[ie] / tau,
        ZETA=req(:ZETA), ZMAG_OVER_A=req(:ZMAG_OVER_A))
end

function _yanccnn_x(ens::YanccNNEnsemble, input_neo::InputNEO)
    inputs = yanccnn_inputs(input_neo)
    return [getfield(inputs, Symbol(replace(name, "_log10" => ""))) for name in ens.xnames]
end

# Ensemble statistics of rows of `y` (outputs × members): the mean, or with `uncertain` the
# mean ± the corrected sample standard deviation over the members as a `Measurement`.
function _yanccnn_stat(v::AbstractVector{T}, uncertain::Bool) where {T<:Real}
    m = Statistics.mean(v)
    uncertain || return m
    T <: ForwardDiff.Dual && error("Yancc-NN: uncertain=true is not differentiable; take derivatives of the ensemble mean")
    return Measurements.measurement(m, length(v) > 1 ? Statistics.std(v; mean=m) : zero(m))
end

function _yanccnn_outputs(ens::YanccNNEnsemble, input_neo::InputNEO, expected, warn_nn_train_bounds::Bool)
    Set(ens.ynames) == Set(expected) || error("Yancc-NN: model $(ens.name) has the outputs $(ens.ynames), not $(expected)")
    y = yanccnn_array(ens, _yanccnn_x(ens, input_neo); warn_nn_train_bounds)
    return Dict(name => view(y, k, :) for (k, name) in enumerate(ens.ynames))
end

const _YANCCNN_FLUX = ["$(f)flux_tgyro_$(s)" for f in ("p", "e", "m") for s in ("ion1", "ion2", "elec")]
const _YANCCNN_FLOW = ["jpar", "vpol_elec", "vpol_ion1", "vpol_ion2"]

"""
    run_yanccnn(input_neo::InputNEO; model="yanccnn_mastu+nstx_flux", uncertain=false, warn_nn_train_bounds=true)
    run_yanccnn(input_neos::Vector{<:InputNEO}; ...)

Neoclassical fluxes from the Yancc-NN surrogate, as a `GACODE.FluxSolution` in gyro-Bohm
units (ensemble mean). `PARTICLE_FLUX_i` has two entries, the bulk ion and the lumped
impurity; `ENERGY_FLUX_i` and `STRESS_TOR_i` sum the two ions.

With `uncertain=true` the values are `Measurements.Measurement`s: ensemble mean ± standard
deviation over the members (for the summed ion channels, of the sum in each member).

`model`: see [`yanccnn_models`](@ref). The shipped models were trained on MAST-U and NSTX
core plasmas (r/a from 0.1 to 0.9, deuterium with carbon) with measured sonic rotation;
`warn_nn_train_bounds` warns when an input is outside the range of the training set. The
signs of the outputs are those of the training devices' field and current directions,
which are not inputs.
"""
function run_yanccnn(input_neo::InputNEO; model::AbstractString="yanccnn_mastu+nstx_flux", uncertain::Bool=false, warn_nn_train_bounds::Bool=true)
    y = _yanccnn_outputs(load_yanccnn_once(model), input_neo, _YANCCNN_FLUX, warn_nn_train_bounds)
    stat(v) = _yanccnn_stat(v, uncertain)
    Qe = stat(y["eflux_tgyro_elec"])
    return GACODE.FluxSolution{typeof(Qe)}(Qe, stat(y["eflux_tgyro_ion1"] .+ y["eflux_tgyro_ion2"]), stat(y["pflux_tgyro_elec"]),
        [stat(y["pflux_tgyro_ion1"]), stat(y["pflux_tgyro_ion2"])], stat(y["mflux_tgyro_ion1"] .+ y["mflux_tgyro_ion2"]))
end

run_yanccnn(input_neos::AbstractVector{<:InputNEO}; kw...) = [run_yanccnn(input_neo; kw...) for input_neo in input_neos]

"""
    YanccNNFlowSolution

Flows from the Yancc-NN flow surrogate, normalised to the bulk ion:

  - `vpol_ion1`, `vpol_ion2`, `vpol_elec`: poloidal velocity at the outboard midplane in
    units of `v_norm = sqrt(T_1/m_D)`;
  - `jpar`: `<j_par B>/B_unit` in units of `e n_1 v_norm`.
"""
struct YanccNNFlowSolution{T<:Real}
    vpol_ion1::T
    vpol_ion2::T
    vpol_elec::T
    jpar::T
end

"""
    run_yanccnn_flow(input_neo::InputNEO; model="yanccnn_mastu+nstx_flow", uncertain=false, warn_nn_train_bounds=true)
    run_yanccnn_flow(input_neos::Vector{<:InputNEO}; ...)

Poloidal flows and parallel current from the Yancc-NN flow surrogate, as a
[`YanccNNFlowSolution`](@ref) (ensemble mean, or mean ± standard deviation over the members
with `uncertain=true`). Options as in [`run_yanccnn`](@ref).
"""
function run_yanccnn_flow(input_neo::InputNEO; model::AbstractString="yanccnn_mastu+nstx_flow", uncertain::Bool=false, warn_nn_train_bounds::Bool=true)
    y = _yanccnn_outputs(load_yanccnn_once(model), input_neo, _YANCCNN_FLOW, warn_nn_train_bounds)
    stat(v) = _yanccnn_stat(v, uncertain)
    return YanccNNFlowSolution(stat(y["vpol_ion1"]), stat(y["vpol_ion2"]), stat(y["vpol_elec"]), stat(y["jpar"]))
end

run_yanccnn_flow(input_neos::AbstractVector{<:InputNEO}; kw...) = [run_yanccnn_flow(input_neo; kw...) for input_neo in input_neos]

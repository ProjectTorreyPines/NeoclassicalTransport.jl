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
import LinearAlgebra
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

# One member for a batch of points (columns of `x`), in place: `y` receives the (normalised)
# outputs; `ws` holds three work matrices (hidden width, points). Every weight matrix is
# read once for the whole batch. Generic in the element type, so dual numbers can be carried
# through.
function _yanccnn_member!(y::AbstractMatrix, ws::NTuple{3,<:AbstractMatrix}, layers::Vector{YanccNNDense}, act::F, nblocks::Integer, x::AbstractMatrix) where {F}
    h, z, t = ws
    LinearAlgebra.mul!(h, layers[1].W, x)
    h .= act.(h .+ layers[1].b)
    k = 1
    for _ in 1:nblocks
        src = h
        for j in 1:3
            dst = isodd(j) ? z : t
            LinearAlgebra.mul!(dst, layers[k+j].W, src)
            dst .= act.(dst .+ layers[k+j].b)
            src = dst
        end
        LinearAlgebra.mul!(t, layers[k+4].W, z)
        h .= act.(t .+ layers[k+4].b .+ h)
        k += 4
    end
    LinearAlgebra.mul!(y, layers[end].W, h)
    y .+= layers[end].b
    return y
end

# Value and derivative of the activations
function _yanccnn_gelu_tanh_d(x)
    c = sqrt(2 / pi)
    t = tanh(c * (x + 0.044715 * x^3))
    return x / 2 * (1 + t), (1 + t) / 2 + x / 2 * (1 - t^2) * c * (1 + 3 * 0.044715 * x^2)
end
_yanccnn_elu_d(x) = x > 0 ? (x, one(x)) : (expm1(x), exp(x))

# One member for a batch of points with the Jacobian of its (normalised) outputs with respect
# to its inputs, in place. `ws` = (h, z, t, d) work matrices (width, points) and (J, Jz, Jt)
# work arrays (width, inputs, points); `Jy` is (outputs, inputs, points).
function _yanccnn_member_jacobian!(y::Matrix{Float64}, Jy::AbstractArray{Float64,3}, ws, layers::Vector{YanccNNDense}, actd::F, nblocks::Integer, x::Matrix{Float64}) where {F}
    h, z, t, d, J, Jz, Jt = ws
    width, nin, np = size(J)
    flat(A) = reshape(A, size(A, 1), nin * np)
    d3 = reshape(d, width, 1, np)
    function activate!(dst, src, b)
        @inbounds for s in axes(dst, 2), i in axes(dst, 1)
            dst[i, s], d[i, s] = actd(src[i, s] + b[i])
        end
    end
    LinearAlgebra.mul!(h, layers[1].W, x)
    activate!(h, h, layers[1].b)
    J .= d3 .* layers[1].W
    k = 1
    for _ in 1:nblocks
        src, Jsrc = h, J
        for j in 1:3
            dst, Jdst = isodd(j) ? (z, Jz) : (t, Jt)
            LinearAlgebra.mul!(dst, layers[k+j].W, src)
            LinearAlgebra.mul!(flat(Jdst), layers[k+j].W, flat(Jsrc))
            activate!(dst, dst, layers[k+j].b)
            Jdst .*= d3
            src, Jsrc = dst, Jdst
        end
        LinearAlgebra.mul!(t, layers[k+4].W, z)
        LinearAlgebra.mul!(flat(Jt), layers[k+4].W, flat(Jz))
        t .+= h
        activate!(h, t, layers[k+4].b)
        J .= d3 .* (Jt .+ J)
        k += 4
    end
    LinearAlgebra.mul!(y, layers[end].W, h)
    y .+= layers[end].b
    LinearAlgebra.mul!(flat(Jy), layers[end].W, flat(J))
    return y, Jy
end

"""
    yanccnn_array(ens::YanccNNEnsemble, x::AbstractVector; warn_nn_train_bounds=false)
    yanccnn_array(ens::YanccNNEnsemble, x::AbstractMatrix; warn_nn_train_bounds=false)

Outputs of every member: a `(n outputs, n members)` matrix for one point, or
`(n outputs, n points, n members)` for a batch of points given as the columns of `x`
(faster than point by point). `x` holds the inputs in the order of `ens.xnames`, as plain
values (the `log10` of the `_log10` inputs is taken here).
"""
function yanccnn_array(ens::YanccNNEnsemble, x::AbstractMatrix{T}; warn_nn_train_bounds::Bool=false) where {T<:Real}
    size(x, 1) == length(ens.xnames) || error("Yancc-NN $(ens.name): $(size(x, 1)) inputs given, $(length(ens.xnames)) expected")
    xn = [endswith(ens.xnames[i], "_log10") ? log10(x[i, s]) : float(x[i, s]) for i in axes(x, 1), s in axes(x, 2)]
    for s in axes(xn, 2), i in axes(xn, 1)
        isfinite(xn[i, s]) || error("Yancc-NN $(ens.name): $(ens.xnames[i]) = $(x[i, s]) is not allowed")
        if warn_nn_train_bounds && !(ens.xmin[i] <= xn[i, s] <= ens.xmax[i])
            @warn "Yancc-NN $(ens.name) extrapolation: $(ens.xnames[i]) = $(ForwardDiff.value(xn[i, s])) is outside the training range [$(ens.xmin[i]), $(ens.xmax[i])]"
        end
    end
    xn .= (xn .- ens.xm) ./ ens.xsigma
    act = ens.activation === :gelu_tanh ? _yanccnn_gelu_tanh : _yanccnn_elu
    return _yanccnn_members(ens, act, xn)
end

function yanccnn_array(ens::YanccNNEnsemble, x::AbstractVector{T}; warn_nn_train_bounds::Bool=false) where {T<:Real}
    return dropdims(yanccnn_array(ens, reshape(x, :, 1); warn_nn_train_bounds); dims=2)
end

function _yanccnn_members(ens::YanccNNEnsemble, act::F, xn::AbstractMatrix{T}) where {F,T}
    width, np = length(ens.members[1][1].b), size(xn, 2)
    ws = (Matrix{T}(undef, width, np), Matrix{T}(undef, width, np), Matrix{T}(undef, width, np))
    y = Array{T,3}(undef, length(ens.ynames), np, length(ens.members))
    for (k, layers) in enumerate(ens.members)
        yk = view(y, :, :, k)
        _yanccnn_member!(yk, ws, layers, act, ens.nblocks, xn)
        yk .= yk .* ens.ysigma .+ ens.ym
    end
    return y
end

# Dual-number inputs. The outputs and their Jacobian with respect to the inputs are computed
# in plain numbers, and the partials of the inputs are carried through that Jacobian. Forward
# mode in chunks evaluates the same point once per chunk with different partials, so the
# result of the last points is kept (exact match of the model and of every input value).
const _YANCCNN_JAC_CAPACITY = 256
const _YANCCNN_JAC_CACHE = Dict{Tuple{String,Vector{Float64}},Tuple{Matrix{Float64},Array{Float64,3}}}()
const _YANCCNN_JAC_ORDER = Tuple{String,Vector{Float64}}[]
const _YANCCNN_JAC_LOCK = ReentrantLock()

# For each point (column of `xv`): outputs (outputs, members) and Jacobian (outputs, inputs, members)
function _yanccnn_members_jacobian(ens::YanccNNEnsemble, xv::Matrix{Float64})
    np = size(xv, 2)
    keys_ = [(ens.name, xv[:, s]) for s in 1:np]
    res = lock(() -> Any[get(_YANCCNN_JAC_CACHE, key, nothing) for key in keys_], _YANCCNN_JAC_LOCK)
    todo = [s for s in 1:np if res[s] === nothing && findfirst(==(keys_[s]), keys_) == s]
    if !isempty(todo)
        actd = ens.activation === :gelu_tanh ? _yanccnn_gelu_tanh_d : _yanccnn_elu_d
        width, nin, nout, nmem, nt = length(ens.members[1][1].b), size(xv, 1), length(ens.ynames), length(ens.members), length(todo)
        ws = (ntuple(_ -> Matrix{Float64}(undef, width, nt), 4)..., ntuple(_ -> Array{Float64,3}(undef, width, nin, nt), 3)...)
        x = xv[:, todo]
        yk = Matrix{Float64}(undef, nout, nt)
        Jk = Array{Float64,3}(undef, nout, nin, nt)
        Y = [Matrix{Float64}(undef, nout, nmem) for _ in todo]
        J = [Array{Float64,3}(undef, nout, nin, nmem) for _ in todo]
        for (k, layers) in enumerate(ens.members)
            _yanccnn_member_jacobian!(yk, Jk, ws, layers, actd, ens.nblocks, x)
            for (q, _) in enumerate(todo)
                Y[q][:, k] .= view(yk, :, q) .* ens.ysigma .+ ens.ym
                J[q][:, :, k] .= view(Jk, :, :, q) .* ens.ysigma
            end
        end
        lock(_YANCCNN_JAC_LOCK) do
            for (q, s) in enumerate(todo)
                res[s] = (Y[q], J[q])
                haskey(_YANCCNN_JAC_CACHE, keys_[s]) && continue
                if length(_YANCCNN_JAC_ORDER) >= _YANCCNN_JAC_CAPACITY
                    delete!(_YANCCNN_JAC_CACHE, popfirst!(_YANCCNN_JAC_ORDER))
                end
                _YANCCNN_JAC_CACHE[keys_[s]] = (Y[q], J[q])
                push!(_YANCCNN_JAC_ORDER, keys_[s])
            end
        end
    end
    for s in 1:np      # repeated points within the batch
        res[s] === nothing && (res[s] = res[findfirst(==(keys_[s]), keys_)])
    end
    return res
end

function _yanccnn_members(ens::YanccNNEnsemble, act::F, xn::AbstractMatrix{ForwardDiff.Dual{Tg,V,N}}) where {F,Tg,V<:AbstractFloat,N}
    res = _yanccnn_members_jacobian(ens, Float64[ForwardDiff.value(v) for v in xn])
    y = Array{ForwardDiff.Dual{Tg,V,N},3}(undef, length(ens.ynames), size(xn, 2), length(ens.members))
    for s in axes(xn, 2)
        Y, J = res[s]::Tuple{Matrix{Float64},Array{Float64,3}}
        for k in axes(Y, 2), i in axes(Y, 1)
            p = zero(ForwardDiff.partials(xn[1, s]))
            for j in axes(xn, 1)
                p += V(J[i, j, k]) * ForwardDiff.partials(xn[j, s])
            end
            y[i, s, k] = ForwardDiff.Dual{Tg}(V(Y[i, k]), p)
        end
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

# Outputs of every member for a batch of InputNEO: name => (points, members) matrix
function _yanccnn_outputs(ens::YanccNNEnsemble, input_neos::AbstractVector{<:InputNEO}, expected, warn_nn_train_bounds::Bool)
    Set(ens.ynames) == Set(expected) || error("Yancc-NN: model $(ens.name) has the outputs $(ens.ynames), not $(expected)")
    y = yanccnn_array(ens, reduce(hcat, [_yanccnn_x(ens, input_neo) for input_neo in input_neos]); warn_nn_train_bounds)
    return Dict(name => view(y, k, :, :) for (k, name) in enumerate(ens.ynames))
end

const _YANCCNN_FLUX = ["$(f)flux_tgyro_$(s)" for f in ("p", "e", "m") for s in ("ion1", "ion2", "elec")]
const _YANCCNN_FLOW = ["jpar", "vpol_elec", "vpol_ion1", "vpol_ion2"]

"""
    run_yanccnn(input_neo::InputNEO; model="yanccnn_mastu+nstx_flux", uncertain=false, warn_nn_train_bounds=false)
    run_yanccnn(input_neos::Vector{<:InputNEO}; ...)    # one batched evaluation, faster than one call per InputNEO

Neoclassical fluxes from the Yancc-NN surrogate, as a `GACODE.FluxSolution` in gyro-Bohm
units (ensemble mean). `PARTICLE_FLUX_i` has two entries, the bulk ion and the lumped
impurity; `ENERGY_FLUX_i` and `STRESS_TOR_i` sum the two ions.

With `uncertain=true` the values are `Measurements.Measurement`s: ensemble mean ± standard
deviation over the members (for the summed ion channels, of the sum in each member).

`model`: see [`yanccnn_models`](@ref). The shipped models were trained on MAST-U and NSTX
(`yanccnn_mastu+nstx_*`, the default) or DIII-D (`yanccnn_d3d_*`) core plasmas (r/a from 0.1 to 0.9,
deuterium with carbon) with measured sonic rotation;
`warn_nn_train_bounds=true` warns when an input is outside the range of the training set. The
signs of the outputs are those of the training devices' field and current directions,
which are not inputs.
"""
function run_yanccnn(input_neos::AbstractVector{<:InputNEO}; model::AbstractString="yanccnn_mastu+nstx_flux", uncertain::Bool=false, warn_nn_train_bounds::Bool=false)
    y = _yanccnn_outputs(load_yanccnn_once(model), input_neos, _YANCCNN_FLUX, warn_nn_train_bounds)
    stat(name, s) = _yanccnn_stat(view(y[name], s, :), uncertain)
    stat(name1, name2, s) = _yanccnn_stat(view(y[name1], s, :) .+ view(y[name2], s, :), uncertain)
    return map(eachindex(input_neos)) do s
        Qe = stat("eflux_tgyro_elec", s)
        GACODE.FluxSolution{typeof(Qe)}(Qe, stat("eflux_tgyro_ion1", "eflux_tgyro_ion2", s), stat("pflux_tgyro_elec", s),
            [stat("pflux_tgyro_ion1", s), stat("pflux_tgyro_ion2", s)], stat("mflux_tgyro_ion1", "mflux_tgyro_ion2", s))
    end
end

run_yanccnn(input_neo::InputNEO; kw...) = only(run_yanccnn([input_neo]; kw...))

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
    run_yanccnn_flow(input_neo::InputNEO; model="yanccnn_mastu+nstx_flow", uncertain=false, warn_nn_train_bounds=false)
    run_yanccnn_flow(input_neos::Vector{<:InputNEO}; ...)

Poloidal flows and parallel current from the Yancc-NN flow surrogate, as a
[`YanccNNFlowSolution`](@ref) (ensemble mean, or mean ± standard deviation over the members
with `uncertain=true`). Options as in [`run_yanccnn`](@ref).
"""
function run_yanccnn_flow(input_neos::AbstractVector{<:InputNEO}; model::AbstractString="yanccnn_mastu+nstx_flow", uncertain::Bool=false, warn_nn_train_bounds::Bool=false)
    y = _yanccnn_outputs(load_yanccnn_once(model), input_neos, _YANCCNN_FLOW, warn_nn_train_bounds)
    stat(name, s) = _yanccnn_stat(view(y[name], s, :), uncertain)
    return [YanccNNFlowSolution(stat("vpol_ion1", s), stat("vpol_ion2", s), stat("vpol_elec", s), stat("jpar", s)) for s in eachindex(input_neos)]
end

run_yanccnn_flow(input_neo::InputNEO; kw...) = only(run_yanccnn_flow([input_neo]; kw...))

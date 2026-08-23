# Export one trained TrainFluxNN NEO-NN ensemble into this repo:
#   1. copy the training BSON to models/neonn_<device>_<group>.bson
#      (dropping the TGLFNNensemble_/tgyro/v1/hyperparameter tokens)
#   2. write models/neonn_<device>_<group>/ with one .onnx per member plus the
#      xnames/ynames/xm/xsigma/ym/ysigma/xbounds_* sidecars (same layout as
#      utilities/convert_nn.ipynb, which this script automates for one model).
#
# Usage (from the repo root, network access required on first run):
#   julia utilities/export_neonn_model.jl <training_ensemble.bson> <device> <group>
# e.g.
#   julia utilities/export_neonn_model.jl \
#     /path/TrainFluxNN/models/gpu/TGLFNNensemble_neonn_tgyro_d3dnegdedge_flux_v1_nmodels_20_..._cuda_true.bson \
#     d3dnegdedge flux

# Pass a pre-built env via --project (with this repo dev'd + the deps below)
# and NCT_EXPORT_SKIP_PKG=1 to skip the temp-env setup.
if get(ENV, "NCT_EXPORT_SKIP_PKG", "0") != "1"
    using Pkg
    Pkg.activate(; temp=true)
    Pkg.develop(; path=dirname(@__DIR__))
    Pkg.add(["ONNXNaiveNASflux", "Flux", "ONNXRunTime"])
end

using NeoclassicalTransport, Flux, ONNXNaiveNASflux, ONNXRunTime
const NCT = NeoclassicalTransport

src_bson, device, group = ARGS[1], ARGS[2], ARGS[3]
group in ("flux", "flow") || error("group must be flux or flow")
isfile(src_bson) || error("not found: $src_bson")

name = "neonn_$(device)_$(group)"
dest = joinpath(dirname(@__DIR__), "models", "$name.bson")
cp(src_bson, dest; force=true)
println("copied -> $dest")

# Memoized loader would hit any previously-loaded copy; load fresh via BSON path
ens = NCT.loadmodel(name)
println("$name: $(length(ens.models)) members, $(length(ens.xnames)) inputs, $(length(ens.ynames)) outputs")

writevec(path, v) = open(path, "w") do io
    for x in v
        println(io, x)
    end
end

outdir = mkpath(joinpath(dirname(@__DIR__), "models", name))
for (i, model) in enumerate(ens.models)
    ONNXNaiveNASflux.save(joinpath(outdir, "$(name)_model_$i.onnx"), Flux.f32(model.fluxmodel))
end
writevec(joinpath(outdir, "xnames.txt"), ens.xnames)
writevec(joinpath(outdir, "ynames.txt"), ens.ynames)
writevec(joinpath(outdir, "xm.txt"), Float32.(ens.xm))
writevec(joinpath(outdir, "xsigma.txt"), Float32.(ens.xσ))
writevec(joinpath(outdir, "ym.txt"), Float32.(ens.ym))
writevec(joinpath(outdir, "ysigma.txt"), Float32.(ens.yσ))
writevec(joinpath(outdir, "xbounds_min.txt"), Float32.(ens.xbounds[:, 1]))
writevec(joinpath(outdir, "xbounds_max.txt"), Float32.(ens.xbounds[:, 2]))
println("exported $(length(ens.models)) members + sidecars -> $outdir")

# Verify ONNX round-trip at the center of the training space (Float32-level)
nx = length(ens.xnames)
xt = Float32.(ens.xm)
xlin = copy(Float64.(xt))
for (i, n) in enumerate(ens.xnames)
    endswith(n, "_log10") && (xlin[i] = 10.0^xlin[i])
end
y_julia = NCT.flux_array(ens, xlin; warn_nn_train_bounds=false)
xn = (xt .- Float32.(ens.xm)) ./ Float32.(ens.xσ)
acc = zeros(Float32, length(ens.ynames))
for i in 1:length(ens.models)
    sess = ONNXRunTime.load_inference(joinpath(outdir, "$(name)_model_$i.onnx"))
    out = sess(Dict(only(ONNXRunTime.input_names(sess)) => permutedims(reshape(xn, nx, 1))))
    global acc
    acc .+= vec(permutedims(out[only(ONNXRunTime.output_names(sess))]))
end
y_onnx = (acc ./ length(ens.models)) .* Float32.(ens.yσ) .+ Float32.(ens.ym)
reldiff = maximum(abs.(y_onnx .- y_julia) ./ (abs.(y_julia) .+ 1e-10))
println("ONNX vs Julia max rel diff = ", round(reldiff; sigdigits=3))
reldiff < 1e-3 || error("ONNX round-trip mismatch ($reldiff) — do not ship")
println("EXPORT OK: $name")

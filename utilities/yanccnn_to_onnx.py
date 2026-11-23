#!/usr/bin/env python3
"""Convert a Yancc-NN model (models/<name>.toml + <name>.bin) to ONNX.

    python utilities/yanccnn_to_onnx.py models/yanccnn_mastu+nstx_flux [outdir] [--check]

Writes into <outdir>/<name>/ (default: next to the model):

  <name>.onnx              the whole ensemble in one graph. Input "x": [batch, n inputs],
                           plain values in the order of xnames.txt (the log10 of the
                           *_log10 inputs is taken in the graph). Outputs "mean" and "std":
                           [batch, n outputs] in the order of ynames.txt, mean and corrected
                           sample standard deviation over the members.
  <name>_model_<k>.onnx    member k alone, normalised input -> normalised output:
                           x_n = (x - xm) / xsigma after the log10, y = y_n * ysigma + ym.
  xnames.txt, ynames.txt, xm.txt, xsigma.txt, ym.txt, ysigma.txt,
  xbounds_min.txt, xbounds_max.txt   one value on each line.

Needs numpy and onnx (Python >= 3.11); --check also needs onnxruntime and compares the
graphs with a NumPy evaluation of the same weights.
"""
import os
import sys
import tomllib

import numpy as np
import onnx
from onnx import TensorProto, helper, numpy_helper

OPSET = 17


def read_model(stem):
    with open(stem + ".toml", "rb") as f:
        meta = tomllib.load(f)
    nin, nout, ch, nb = meta["nin"], meta["nout"], meta["channels"], meta["nblocks"]
    shapes = [(ch, nin)] + [(ch, ch)] * (4 * nb) + [(nout, ch)]
    raw = np.fromfile(stem + ".bin", dtype="<f4")
    if raw.size != meta["nmodels"] * sum(o * i + o for o, i in shapes):
        raise ValueError(f"{stem}.bin: unexpected size {raw.size}")
    members, p = [], 0
    for _ in range(meta["nmodels"]):
        layers = []
        for o, i in shapes:
            W = raw[p:p + o * i].reshape((o, i), order="F").copy()  # column-major on disk
            p += o * i
            b = raw[p:p + o].copy()
            p += o
            layers.append((W, b))
        members.append(layers)
    return meta, members


class Graph:
    def __init__(self):
        self.nodes, self.inits, self.n = [], [], 0

    def const(self, value, name=None):
        name = name or f"c{len(self.inits)}"
        self.inits.append(numpy_helper.from_array(np.asarray(value, dtype=np.float32), name))
        return name

    def op(self, kind, inputs, **attrs):
        self.n += 1
        out = f"t{self.n}"
        self.nodes.append(helper.make_node(kind, inputs, [out], **attrs))
        return out

    def act(self, x, kind):
        if kind == "elu":
            return self.op("Elu", [x], alpha=1.0)
        # gelu_tanh(x) = x/2 (1 + tanh(sqrt(2/pi) (x + 0.044715 x^3)))
        x3 = self.op("Mul", [self.op("Mul", [x, x]), x])
        inner = self.op("Add", [x, self.op("Mul", [x3, self.c044])])
        t = self.op("Tanh", [self.op("Mul", [inner, self.csq])])
        return self.op("Mul", [self.op("Mul", [x, self.chalf]), self.op("Add", [t, self.cone])])

    def scalars(self):
        self.c044, self.csq = self.const(0.044715), self.const(np.sqrt(2 / np.pi))
        self.chalf, self.cone = self.const(0.5), self.const(1.0)

    def member(self, x, layers, kind, nblocks, tag):
        def dense(h, k):
            W, b = layers[k]
            return self.op("Gemm", [h, self.const(W, f"{tag}W{k}"), self.const(b, f"{tag}b{k}")], transB=1)

        h = self.act(dense(x, 0), kind)
        k = 1
        for _ in range(nblocks):
            z = h
            for j in range(3):
                z = self.act(dense(z, k + j), kind)
            h = self.act(self.op("Add", [dense(z, k + 3), h]), kind)
            k += 4
        return dense(h, k)

    def model(self, name, inputs, outputs):
        graph = helper.make_graph(self.nodes, name, inputs, outputs, initializer=self.inits)
        model = helper.make_model(graph, opset_imports=[helper.make_opsetid("", OPSET)], producer_name="NeoclassicalTransport.jl")
        model.ir_version = 8
        onnx.checker.check_model(model)
        return model


def member_model(meta, layers, name):
    g = Graph()
    g.scalars()
    y = g.member("x", layers, meta["activation"], meta["nblocks"], "")
    g.nodes.append(helper.make_node("Identity", [y], ["y"]))
    return g.model(name, [helper.make_tensor_value_info("x", TensorProto.FLOAT, ["batch", meta["nin"]])],
                   [helper.make_tensor_value_info("y", TensorProto.FLOAT, ["batch", meta["nout"]])])


def ensemble_model(meta, members):
    g = Graph()
    g.scalars()
    is_log = np.array([n.endswith("_log10") for n in meta["xnames"]])
    x = "x"
    if is_log.any():
        # log10 on the *_log10 columns only (the other columns may be negative)
        mask = g.const(is_log.astype(np.float32)[None, :])
        keep = g.const((~is_log).astype(np.float32)[None, :])
        safe = g.op("Add", [g.op("Mul", [x, mask]), keep])  # 1 where no log is taken
        logx = g.op("Mul", [g.op("Log", [safe]), g.const(1 / np.log(10.0))])
        x = g.op("Add", [g.op("Mul", [x, keep]), g.op("Mul", [logx, mask])])
    xn = g.op("Div", [g.op("Sub", [x, g.const(np.array(meta["xm"])[None, :])]), g.const(np.array(meta["xsigma"])[None, :])])
    ysig, ym = g.const(np.array(meta["ysigma"])[None, :]), g.const(np.array(meta["ym"])[None, :])
    outs = []
    for k, layers in enumerate(members):
        y = g.op("Add", [g.op("Mul", [g.member(xn, layers, meta["activation"], meta["nblocks"], f"m{k}_"), ysig]), ym])
        outs.append(g.op("Unsqueeze", [y, "axes0"]))
    g.inits.append(numpy_helper.from_array(np.array([0], dtype=np.int64), "axes0"))
    stack = g.op("Concat", outs, axis=0)  # [members, batch, outputs]
    mean = g.op("ReduceMean", [stack], axes=[0], keepdims=0)
    dev = g.op("Sub", [stack, g.op("Unsqueeze", [mean, "axes0"])])
    n = len(members)
    var = g.op("Mul", [g.op("ReduceMean", [g.op("Mul", [dev, dev])], axes=[0], keepdims=0), g.const(n / max(n - 1, 1))])
    g.nodes.append(helper.make_node("Identity", [mean], ["mean"]))
    g.nodes.append(helper.make_node("Sqrt", [var], ["std"]))
    nin, nout = meta["nin"], meta["nout"]
    return g.model(meta["name"], [helper.make_tensor_value_info("x", TensorProto.FLOAT, ["batch", nin])],
                   [helper.make_tensor_value_info(o, TensorProto.FLOAT, ["batch", nout]) for o in ("mean", "std")])


def numpy_member(layers, kind, nblocks, x):
    def act(v):
        if kind == "elu":
            return np.where(v > 0, v, np.expm1(v))
        return 0.5 * v * (1 + np.tanh(np.sqrt(2 / np.pi) * (v + 0.044715 * v**3)))

    def dense(h, k):
        return h @ layers[k][0].T.astype(np.float64) + layers[k][1]

    h = act(dense(x, 0))
    k = 1
    for _ in range(nblocks):
        z = h
        for j in range(3):
            z = act(dense(z, k + j))
        h = act(dense(z, k + 3) + h)
        k += 4
    return dense(h, k)


def check(meta, members, outdir):
    import onnxruntime as ort

    rng = np.random.default_rng(0)
    lo, hi = np.array(meta["xmin"]), np.array(meta["xmax"])
    xt = lo + rng.random((64, meta["nin"])) * (hi - lo)  # in the space of the networks (after log10)
    is_log = np.array([n.endswith("_log10") for n in meta["xnames"]])
    xraw = np.where(is_log, 10.0**xt, xt)
    xn = (xt - np.array(meta["xm"])) / np.array(meta["xsigma"])
    ref = np.stack([numpy_member(l, meta["activation"], meta["nblocks"], xn) for l in members])
    scale = np.abs(ref).max(axis=(0, 1))
    worst = 0.0
    for k in range(len(members)):
        sess = ort.InferenceSession(os.path.join(outdir, f"{meta['name']}_model_{k + 1}.onnx"))
        worst = max(worst, (np.abs(sess.run(None, {"x": xn.astype(np.float32)})[0] - ref[k]) / scale).max())
    phys = ref * np.array(meta["ysigma"]) + np.array(meta["ym"])
    sess = ort.InferenceSession(os.path.join(outdir, meta["name"] + ".onnx"))
    mean, std = sess.run(None, {"x": xraw.astype(np.float32)})
    pscale = np.abs(phys).max(axis=(0, 1))
    dmean = (np.abs(mean - phys.mean(axis=0)) / pscale).max()
    dstd = (np.abs(std - phys.std(axis=0, ddof=1)) / pscale).max()
    print(f"check on 64 points: members {worst:.1e}, ensemble mean {dmean:.1e}, std {dstd:.1e} (relative to the output range)")
    if max(worst, dmean, dstd) > 1e-4:
        raise SystemExit("ONNX graphs do not reproduce the weights")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        raise SystemExit(__doc__)
    stem = args[0][:-5] if args[0].endswith(".toml") else args[0]
    meta, members = read_model(stem)
    outdir = os.path.join(args[1] if len(args) > 1 else os.path.dirname(stem), meta["name"])
    os.makedirs(outdir, exist_ok=True)
    for k, layers in enumerate(members):
        onnx.save(member_model(meta, layers, f"{meta['name']}_model_{k + 1}"), os.path.join(outdir, f"{meta['name']}_model_{k + 1}.onnx"))
    onnx.save(ensemble_model(meta, members), os.path.join(outdir, meta["name"] + ".onnx"))
    for fname, values in (("xnames", meta["xnames"]), ("ynames", meta["ynames"]), ("xm", meta["xm"]), ("xsigma", meta["xsigma"]),
                          ("ym", meta["ym"]), ("ysigma", meta["ysigma"]), ("xbounds_min", meta["xmin"]), ("xbounds_max", meta["xmax"])):
        with open(os.path.join(outdir, fname + ".txt"), "w") as f:
            f.writelines(f"{v}\n" for v in values)
    print(f"{meta['name']}: {len(members)} members, {meta['nin']} inputs, {meta['nout']} outputs -> {outdir}")
    if "--check" in sys.argv:
        check(meta, members, outdir)


if __name__ == "__main__":
    main()

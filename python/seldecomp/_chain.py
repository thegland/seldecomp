"""chain: the decomposition of each pair of adjacent values of by(), and their total."""
import numpy as np
import pandas as pd

from ._estimate import Fit, b_layout, estimate
from ._prepare import fmt_val
from ._results import (CHECK_COLS, MISS_COLS, TITLE, _miss_rows, _stage_frame, common_macros,
                       stage_weight_names)

NAN = float("nan")


def estimate_chain(ctx):
    """One decomposition per pair of adjacent values (the lower value is the pair's focal
    group), on the observations of the two values."""
    vals = ctx.vals
    by = ctx.byvals
    chain = Fit()
    chain.is_chain = True
    chain.fits = []
    chain.messages = []
    for j in range(len(vals) - 1):
        tu = ctx.touse & np.array([(v == vals[j]) or (v == vals[j + 1]) for v in by])
        f = estimate(ctx, tu, vals[j], vals[j + 1], pair=j + 1)
        chain.fits.append(f)
        chain.messages += f.messages
    chain.b_names, chain.b_values = chain_b(chain, ctx)
    return chain


def chain_b(chain, ctx):
    """p<j> equations of every pair, then the total: the sums of the pairs' differences,
    explained and unexplained parts (and of their marginal explained parts by stage)."""
    names, vals = [], []
    for j, f in enumerate(chain.fits):
        n, v = b_layout(f, ctx.ylist, prefix="p%d" % (j + 1))
        names += n
        vals += v
    f0 = chain.fits[0]
    nst = len(f0.stages)
    tot = [0.0, 0.0, 0.0]
    marg = [0.0] * nst
    for f in chain.fits:
        R = f.E.res[0]
        tot[0] += f.E.m1[0] - f.E.m2[0]
        tot[1] += R[nst - 1, 1]
        tot[2] += R[nst - 1, 2]
        for s in range(nst):
            marg[s] += R[s, 4]
    names += ["total:difference", "total:explained", "total:unexplained"]
    vals += tot
    if f0.mode:
        names += ["total_marg_expl:%s" % n for n in f0.stagenames]
        vals += marg
    return names, vals


def merge_pair_notes(fits):
    """The notes of the pairs, each headed "p<j>: "; pairs with the same note share it
    ("p1, p3: ...")."""
    units = []
    for j, f in enumerate(fits):
        cur = None
        for nt in f.notes:
            if nt.startswith("    ") and cur is not None:
                cur[1].append(nt)
            else:
                cur = [nt, []]
                units.append((j + 1, cur))
    merged = []
    for j, (head, cont) in units:
        for m in merged:
            if m[1] == head and m[2] == cont:
                m[0].append(j)
                break
        else:
            merged.append([[j], head, cont])
    out = []
    for pairs, head, cont in merged:
        out.append("%s: %s" % (", ".join("p%d" % p for p in pairs), head))
        out += cont
    return out


def chain_results(ctx, chain, res):
    o = ctx.opt
    fits = chain.fits
    e, M = res.e, res.matrices
    f0, fl = fits[0], fits[-1]
    J = len(fits)
    nst = len(f0.stages)
    M["b"] = pd.DataFrame([chain.b_values], index=["y1"], columns=chain.b_names)
    b = dict(zip(chain.b_names, chain.b_values))
    # the estimation sample: the observations of the pairs (insupport leaves some out)
    n = len(ctx.df)
    sample = np.zeros(n, dtype=bool)
    for f in fits:
        sample[f.idx] = True
    e["N"] = int(sample.sum())
    e["N_pairs"] = J
    e["N_stages"] = nst
    e["N_blocks"] = len(f0.blocks)
    e["N_terms"] = len(f0.terms)
    e["N_zerowt"] = ctx.N_zerowt
    e["mincell"] = o["mincell"]
    e["level"] = ctx.level
    e["mean_focal"] = f0.E.m1[0]
    e["mean_ref"] = fl.E.m2[0]
    e["difference"] = b["total:difference"]
    e["explained"] = b["total:explained"]
    e["unexplained"] = b["total:unexplained"]
    e["pct_explained"] = 100 * e["explained"] / e["difference"] if e["difference"] != 0 else NAN
    e["converged"] = min(float(f.E.res[0][:, 5].min()) for f in fits)
    e["targets_met"] = min(f.targets_met for f in fits)
    staged = f0.mode != ""
    e["k_eq"] = J * (1 + 4 * staged) + 1 + (1 if staged else 0)
    e["diffsign"] = "focal_group - ref_group"
    e["cmd"] = "seldecomp"
    e["cmdline"] = ctx.cmdline
    e["title"] = TITLE
    e["estat_cmd"] = "seldecomp_estat"
    e["chain"] = "chain"
    e["depvar"] = ctx.depvar
    e["by"] = ctx.by
    e["by_levels"] = " ".join(fmt_val(v) for v in ctx.vals)
    # (a bar inside a label would read as the separator)
    e["by_labels"] = " | ".join(ctx.group_label(v).replace("|", "/") for v in ctx.vals)
    e["focal_group"] = fmt_val(ctx.vals[0])
    e["ref_group"] = fmt_val(ctx.vals[-1])
    e["focal_group_label"] = ctx.group_label(ctx.vals[0])
    e["ref_group_label"] = ctx.group_label(ctx.vals[-1])
    common_macros(ctx, f0, e)
    rows = []
    for f in fits:
        R = f.E.res[0]
        last = R[nst - 1]
        row = [f.E.m1[0], last[0], f.E.m2[0], f.E.m1[0] - f.E.m2[0], last[1], last[2], last[3],
               int(f.g1.sum()), int((~f.g1).sum()), f.E.W1, f.E.W2]
        if o["insupport"]:
            row += [f.ins["N_focal_all"], f.ins["mean_focal_all"][0]]
        rows.append(row)
    cols = ["mean_focal", "mean_raked_focal", "mean_ref", "difference", "explained", "unexplained",
            "pct_explained", "N_focal", "N_ref", "W_focal", "W_ref"]
    if o["insupport"]:
        cols += ["N_focal_all", "mean_focal_all"]
    pn = ["p%d" % (j + 1) for j in range(J)]
    M["by_pair"] = pd.DataFrame(rows, index=pn, columns=cols)
    M["pairs"] = pd.DataFrame([[f.focal, f.ref] for f in fits], index=pn, columns=["focal_group", "ref_group"])
    for j, f in enumerate(fits):
        sfx = "" if j == 0 else str(j + 1)
        M["bystage" + sfx] = _stage_frame(f.E.res[0], f.stagenames, False)
        M["raking_checks" + sfx] = pd.DataFrame(f.E.checks, index=f.stagenames, columns=CHECK_COLS)
        miss = _miss_rows(f.E)
        if miss:
            M["raking_misses" + sfx] = pd.DataFrame(miss, index=["r1"] * len(miss), columns=MISS_COLS)
            for k, (s, t, l, tg, ach, d) in enumerate(f.E.misses):
                key = "raking_miss%s%d" % (sfx + "_" if sfx else "", k + 1)
                e[key + "_var"] = f.terms[t].name
                e[key + "_level"] = f.coded[t].labels[l - 1]
        k = 0
        for v, q in ctx.contlist:
            k += 1
            M["cuts%d_p%d" % (k, j + 1)] = pd.DataFrame([f.cuts[v]], index=[v],
                                                       columns=["cut%d" % (i + 1) for i in range(q - 1)])
    if o["ratio"]:
        rr, vv = [], []
        for j, f in enumerate(fits):
            rr += ["p%d:focal_group" % (j + 1), "p%d:raked_focal" % (j + 1), "p%d:ref_group" % (j + 1)]
            vv += [[f.E.n1[0], f.E.d1], [f.E.n1r[0], f.E.d1r], [f.E.n2[0], f.E.d2]]
        M["ratio_parts"] = pd.DataFrame(vv, index=rr, columns=["numerator", "denominator"])
    res.notes = list(ctx.xnotes) + merge_pair_notes(fits)
    for i, nt in enumerate(res.notes):
        e["note%d" % (i + 1)] = nt
    e["N_notes"] = len(res.notes)
    res.messages = list(ctx.messages) + list(chain.messages)
    res.trace = list(ctx.tracelines)
    res.sample = pd.Series(sample, index=ctx.df.index, name="sample")
    wcols = {}
    for j, f in enumerate(fits):
        for s, nm in enumerate(stage_weight_names("rw_p%d" % (j + 1), f)):
            col = np.full(n, NAN)
            col[f.idx] = np.where(f.g1, f.E.rw[s], f.w)
            wcols[nm] = col
    res.weights = pd.DataFrame(wcols, index=ctx.df.index)

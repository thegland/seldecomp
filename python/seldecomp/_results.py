"""The results of a seldecomp run."""
import numpy as np
import pandas as pd

from ._engine import STAGE_COLS
from ._prepare import fmt_val
from .errors import SeldecompError

NAN = float("nan")
DIAG = STAGE_COLS[5:]
CHECK_COLS = ["iterations_converged", "missed_variables", "missed_targets", "targets_checked", "control_scale"]
MISS_COLS = ["stage", "term", "level", "target", "achieved", "discrepancy"]
TITLE = "Raking decomposition (Pylypchuk-Selden)"


class SeldecompResults:
    """What seldecomp returns.

    e           the stored results that are numbers or text, under the names of the Stata
                command's e() results (e["explained"], e["N_focal"], e["stagenames"], ...)
    b           the decomposition results, a Series indexed by "equation:name"
    V           the bootstrap variance-covariance matrix of b (None without bayesboot)
    matrices    every stored matrix as a DataFrame, by name (bystage, raking_checks, ...);
                each is also an attribute: res.bystage
    notes       the note lines
    weights     the weights of each stage, one column per stage: the raked weights of the focal
                group and the base weights of the reference group
    sample      a boolean Series marking the estimation sample
    figure      the matplotlib Figure drawn by the makeplot option (None without it)
    bsdraws     the bootstrap replicates (None without bayesboot)
    bsweights   the bootstrap weights of each replicate (bayesboot(generate()) only)
    """

    def __init__(self):
        self.e = {}
        self.matrices = {}
        self.notes = []
        self.messages = []
        self.trace = []
        self.weights = None
        self.sample = None
        self.bsweights = None
        self.bsdraws = None
        self.figure = None
        self._generated = {}
        self._replay = {}

    def __getattr__(self, name):
        m = self.__dict__.get("matrices", {})
        if name in m:
            return m[name]
        raise AttributeError(name)

    @property
    def b(self):
        return self.matrices["b"].iloc[0]

    @property
    def V(self):
        return self.matrices.get("V")

    def _finish(self, data):
        """Write the generated variables into the data."""
        import warnings
        with warnings.catch_warnings():
            # many replicate-weight columns are added one by one
            warnings.simplefilter("ignore", pd.errors.PerformanceWarning)
            # pandas 2 warns when the data is a slice of another frame; the columns go into
            # the data given, as documented (pandas 3 has no such warning)
            if hasattr(pd.errors, "SettingWithCopyWarning"):
                warnings.simplefilter("ignore", pd.errors.SettingWithCopyWarning)
            for name, values in self._generated.items():
                data[name] = values

    # the output: Display writes the tables
    def display(self, **options):
        from ._display import display
        return display(self, **options)

    def __str__(self):
        from ._display import render
        return render(self)

    __repr__ = __str__

    def estat_bootstrap(self, **options):
        from ._display import estat_bootstrap
        return estat_bootstrap(self, **options)

    def makeplot(self, type=None, **options):
        from ._plot import makeplot
        return makeplot(self, type, **options)


def _stage_frame(R, rows, individually):
    cols = list(STAGE_COLS)
    if individually:
        cols[1], cols[2], cols[4], cols[26] = "explained", "unexplained", "contribution", "pct_contribution"
    return pd.DataFrame(R.copy(), index=list(rows), columns=cols)


def _miss_rows(E, term_of=lambda t: t):
    return [(s + 1, term_of(t) + 1, l, tg, ach, d) for s, t, l, tg, ach, d in E.misses]


def fit_results(ctx, fit, res, sfx=""):
    """The stored results of one decomposition; sfx: the suffix of a chain pair's matrices."""
    o = ctx.opt
    E = fit.E
    e = res.e
    M = res.matrices
    nst = len(fit.stages)
    R = E.res[0]
    variant = fit.mode == "individually"
    nout = len(ctx.ylist)
    g1 = fit.g1
    e["N"] = len(fit.idx)
    e["N_focal"] = int(g1.sum())
    e["N_ref"] = int((~g1).sum())
    e["W_focal"] = E.W1
    e["W_ref"] = E.W2
    e["N_zerowt"] = ctx.N_zerowt
    e["mean_focal"] = E.m1[0]
    e["mean_raked_focal"] = R[nst - 1, 0]
    e["mean_ref"] = E.m2[0]
    e["difference"] = E.m1[0] - E.m2[0]
    e["explained"] = R[nst - 1, 1]
    e["unexplained"] = R[nst - 1, 2]
    e["pct_explained"] = R[nst - 1, 3]
    e["diffsign"] = "focal_group - ref_group"
    e["N_stages"] = nst
    staged = fit.mode != ""
    e["k_eq"] = nout * (1 + (2 if variant else 4) * staged + 4 * (fit.ER is not None) + (fit.ED is not None))
    e["level"] = ctx.level
    for c, name in enumerate(DIAG):
        e[name] = R[nst - 1, 5 + c]
    cmin = min(float(X.res[0][:, 5].min()) for X in (E, fit.ER, fit.ED) if X is not None)
    e["converged"] = cmin
    e["targets_met"] = fit.targets_met
    M["bystage" + sfx] = _stage_frame(R, fit.stagenames, variant)
    for j in range(1, nout):
        M["bystage%d" % (j + 1)] = _stage_frame(E.res[j], fit.stagenames, variant)
    M["raking_checks" + sfx] = pd.DataFrame(E.checks, index=fit.stagenames, columns=CHECK_COLS)
    miss = _miss_rows(E)
    if miss:
        M["raking_misses" + sfx] = pd.DataFrame(miss, index=["r1"] * len(miss), columns=MISS_COLS)
        for k, (s, t, l, tg, ach, d) in enumerate(E.misses):
            key = "raking_miss%s%d" % (sfx + "_" if sfx else "", k + 1)
            e[key + "_var"] = fit.terms[t].name
            e[key + "_level"] = fit.coded[t].labels[l - 1]
    if fit.ER is not None:
        M["bystage_rev"] = _stage_frame(fit.ER.res[0], fit.rnames, False)
        for j in range(1, nout):
            M["bystage_rev%d" % (j + 1)] = _stage_frame(fit.ER.res[j], fit.rnames, False)
        M["raking_checks_rev"] = pd.DataFrame(fit.ER.checks, index=fit.rnames, columns=CHECK_COLS)
        miss = _miss_rows(fit.ER, lambda t: fit.rorder[t])
        if miss:
            M["raking_misses_rev"] = pd.DataFrame(miss, index=["r1"] * len(miss), columns=MISS_COLS)
            for k, (s, t, l, tg, ach, d) in enumerate(fit.ER.misses):
                to = fit.rorder[t]
                e["raking_miss_rev%d_var" % (k + 1)] = fit.terms[to].name
                e["raking_miss_rev%d_level" % (k + 1)] = fit.coded[to].labels[l - 1]
        e["reverse"] = "reverse"
    if fit.ED is not None:
        M["bystage_detail"] = _stage_frame(fit.ED.res[0], fit.dnames, False)
        for j in range(1, nout):
            M["bystage_detail%d" % (j + 1)] = _stage_frame(fit.ED.res[j], fit.dnames, False)
        if fit.detail_same:
            M["raking_checks_detail"] = pd.DataFrame(E.checks, index=fit.dnames, columns=CHECK_COLS)
            if "raking_misses" in M:
                M["raking_misses_detail"] = M["raking_misses"].copy()
        else:
            M["raking_checks_detail"] = pd.DataFrame(fit.ED.checks, index=fit.dnames, columns=CHECK_COLS)
            miss = _miss_rows(fit.ED)
            if miss:
                M["raking_misses_detail"] = pd.DataFrame(miss, index=["r1"] * len(miss), columns=MISS_COLS)
        if "raking_misses_detail" in M:
            src = E if fit.detail_same else fit.ED
            for k, (s, t, l, tg, ach, d) in enumerate(src.misses):
                e["raking_miss_detail%d_var" % (k + 1)] = fit.terms[t].name
                e["raking_miss_detail%d_level" % (k + 1)] = fit.coded[t].labels[l - 1]
        e["detail"] = "detail"
    if e["N_focal_zero"] > 0:
        e["explained_support"] = e["mean_focal"] - e["mean_focal_cs"]
        e["explained_cs"] = e["mean_focal_cs"] - e["mean_raked_focal"]
    wt, wl = E.worst[nst - 1]
    if e["badcontrols"] > 0 and wt is not None:
        e["worstvar"] = fit.terms[wt].name
        e["worstcat"] = fit.coded[wt].labels[wl - 1] if wl is not None else ""
    if o["ratio"]:
        rows, vals = [], []
        for j, y in enumerate(ctx.ylist):
            rows += ["%s:focal_group" % y, "%s:raked_focal" % y, "%s:ref_group" % y]
            vals += [[E.n1[j], E.d1], [E.n1r[j], E.d1r], [E.n2[j], E.d2]]
        M["ratio_parts"] = pd.DataFrame(vals, index=rows, columns=["numerator", "denominator"])
    if fit.ins is not None:
        e["N_focal_all"] = fit.ins["N_focal_all"]
        e["W_focal_all"] = fit.ins["W_focal_all"]
        e["mean_focal_all"] = fit.ins["mean_focal_all"][0]


def build_results(ctx, fit):
    o = ctx.opt
    res = SeldecompResults()
    res._ctx = ctx
    res._fit = fit
    e = res.e
    M = res.matrices
    if getattr(fit, "is_chain", False):
        from ._chain import chain_results
        chain_results(ctx, fit, res)
        return res
    E = fit.E
    terms, blocks = fit.terms, fit.blocks
    nout = len(ctx.ylist)
    fit_results(ctx, fit, res)
    M["b"] = pd.DataFrame([fit.b_values], index=["y1"], columns=fit.b_names)
    e["N_blocks"] = len(blocks)
    e["N_terms"] = len(terms)
    e["mincell"] = o["mincell"]
    e["cmd"] = "seldecomp"
    e["cmdline"] = ctx.cmdline
    e["title"] = TITLE
    e["estat_cmd"] = "seldecomp_estat"
    e["depvar"] = " ".join(ctx.ylist)
    if nout > 1:
        e["outcomes"] = " ".join(ctx.ylist)
        e["N_outcomes"] = nout
    e["by"] = ctx.by
    e["focal_group"] = fmt_val(fit.focal)
    e["ref_group"] = fmt_val(fit.ref)
    e["focal_group_label"] = ctx.group_label(fit.focal)
    e["ref_group_label"] = ctx.group_label(fit.ref)
    common_macros(ctx, fit, e)
    for t, c in enumerate(fit.coded):
        e["levels%d" % (t + 1)] = " | ".join(c.labels)
    k = 0
    for v, q in ctx.contlist:
        k += 1
        M["cuts%d" % k] = pd.DataFrame([fit.cuts[v]], index=[v], columns=["cut%d" % (i + 1) for i in range(q - 1)])
    if ctx.contlist:
        e["cutbins"] = " ".join(str(n) for n in fit.cutbins)
    if nout > 1:
        rows = []
        for j in range(nout):
            R = E.res[j]
            last = R[len(fit.stages) - 1]
            row = [E.m1[j], last[0], E.m2[j], E.m1[j] - E.m2[j], last[1], last[2], last[3], NAN, NAN, NAN]
            if fit.ins is not None:
                row.append(fit.ins["mean_focal_all"][j])
            rows.append(row)
        cols = ["mean_focal", "mean_raked_focal", "mean_ref", "difference", "explained", "unexplained",
                "pct_explained", "pct_explained_se", "pct_explained_ll", "pct_explained_ul"]
        if fit.ins is not None:
            cols.append("mean_focal_all")
        M["by_outcome"] = pd.DataFrame(rows, index=ctx.ylist, columns=cols)
    if o["covariatedetail"]:
        M["covariates"] = covariate_detail(fit)
    res.notes = list(fit.notes)
    for i, nt in enumerate(res.notes):
        e["note%d" % (i + 1)] = nt
    e["N_notes"] = len(res.notes)
    res.messages = list(ctx.messages) + list(fit.messages)
    res.trace = list(ctx.tracelines)
    # the raked weights of each stage, and the sample
    n = len(ctx.df)
    sample = np.zeros(n, dtype=bool)
    sample[fit.idx] = True
    res.sample = pd.Series(sample, index=ctx.df.index, name="sample")
    wcols = {}
    names = stage_weight_names("rw", fit)
    for s, nm in enumerate(names):
        col = np.full(n, NAN)
        col[fit.idx] = np.where(fit.g1, E.rw[s], fit.w)
        wcols[nm] = col
    res.weights = pd.DataFrame(wcols, index=ctx.df.index)
    if ctx.gen:
        gnames = stage_weight_names(ctx.gen, fit)
        for nm in gnames:
            if nm in ctx.df.columns and not ctx.genreplace:
                raise SeldecompError("variable %s already defined" % nm, 110)
        for s, nm in enumerate(gnames):
            res._generated[nm] = res.weights.iloc[:, s].to_numpy()
        e["generate"] = ctx.gen
    return res


def stage_weight_names(stub, fit):
    """generate(): the last stage's weights under the name given, the earlier stages' as
    name_s1, name_s2, ..."""
    nst = len(fit.stages)
    return ["%s_s%d" % (stub, s + 1) for s in range(nst - 1)] + [stub]


def common_macros(ctx, fit, e):
    o = ctx.opt
    terms, blocks = fit.terms, fit.blocks
    if ctx.wtype:
        e["wtype"] = ctx.wtype
        e["wexp"] = "= %s" % ctx.wexp
    e["rakespec"] = ctx.rakespec
    e["stagenames"] = " ".join(fit.stagenames)
    e["blocks"] = " ".join(b.name for b in blocks)
    e["rakevars"] = " ".join(t.name for t in terms)
    for t, term in enumerate(terms):
        e["term%d" % (t + 1)] = term.name
    for b, blk in enumerate(blocks):
        e["block%d" % (b + 1)] = " ".join(terms[t].name for t in blk.terms)
    for s, st in enumerate(fit.stages):
        e["stage%d" % (s + 1)] = " ".join(terms[t].name for t in st)
    e["rakeopts"] = ctx.rakeopts
    if fit.mode:
        e["stagetype"] = fit.mode
    if o["fixothermargins"]:
        e["fixothermargins"] = "fixothermargins"
    if o["ratio"]:
        e["ratio"] = o["ratio"]
    if o["insupport"]:
        e["insupport"] = "insupport"
    if ctx.contlist:
        e["contvars"] = ctx.conttext
        e["contsample"] = ctx.contsample
        e["cutvars"] = " ".join(v for v, q in ctx.contlist)
        e["cutq"] = " ".join(str(q) for v, q in ctx.contlist)
        if ctx.recut:
            e["recut"] = "recut"


def covariate_detail(fit):
    """The weighted share of each level of each term in the focal group, the reference group
    and the raked focal group, with the unweighted counts."""
    g1, w = fit.g1, fit.w
    rw = fit.E.rw[-1]
    W1 = float(np.sum(w[g1]))
    W2 = float(np.sum(w[~g1]))
    WR = float(np.nansum(rw[g1]))
    rows, vals = [], []
    for t, c in enumerate(fit.coded):
        for l in range(1, c.nlev + 1):
            m1 = g1 & (c.codes == l)
            m2 = (~g1) & (c.codes == l)
            vals.append([float(np.sum(w[m1])) / W1 if m1.any() else 0.0,
                         float(np.sum(w[m2])) / W2 if m2.any() else 0.0,
                         float(np.nansum(rw[m1])) / WR if (m1.any() and WR > 0) else 0.0,
                         int(m1.sum()), int(m2.sum())])
            rows.append("%d:%d" % (t + 1, l))
    return pd.DataFrame(vals, index=rows, columns=["focal_group", "ref_group", "raked_focal", "N_focal", "N_ref"])

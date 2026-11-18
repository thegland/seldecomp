"""bayesboot(): Bayesian bootstrap replicates, standard errors and confidence intervals."""
import numpy as np
import pandas as pd
from scipy import stats

from ._estimate import RepFail, rep_values
from ._prepare import bin_codes, cut_points
from ._parse import is_string_series
from ._stata import invnormal, normal, pctile_plain
from ._stata import macro_value as M
from .errors import SeldecompError

NAN = float("nan")
CI_TYPES = ("normal", "t", "percentile", "bc", "bca")


def ci_types(bs):
    """The interval types shown, the first being the coefficient table's: those requested, in
    the order typed, then those that citable adds; without any, t with cluster() and normal
    otherwise."""
    cluster = bool(bs["cluster"])
    dflt = "t" if cluster else "normal"
    if bs["normal"] or bs["tinterval"]:
        dflt = ""
    flag = {"normal": bs["normal"], "t": bs["tinterval"], "percentile": bs["percentile"],
            "bc": bs["bc"], "bca": bs["bca"]}
    types = [t for t in CI_TYPES if flag[t] or (bs["citable"] and t in (dflt, "percentile", "bc"))]
    explicit = [("t" if t == "tinterval" else t) for t in bs.get("_order", []) if flag.get("t" if t == "tinterval" else t)]
    types = explicit + [t for t in types if t not in explicit]
    if not types:
        types = ["t" if cluster else "normal"]
    return types


def recode(fit, ctx, w):
    """recut: the codes of every term with a binned variable, from cut points recomputed with
    the weights w. Returns (codes, nlev) for all terms."""
    df = ctx.df
    g1 = fit.g1
    cuts = {}
    for v, q in ctx.contlist:
        x = df[v].to_numpy(dtype=float)[fit.idx]
        if ctx.contsample == "pooled":
            cond = np.ones(len(x), dtype=bool)
        elif ctx.contsample == "focalgroup":
            cond = g1
        else:
            cond = ~g1
        cond = cond & (w > 0)
        cuts[v] = cut_points(x[cond], w[cond], q, v)
    codes, nlev = [], []
    for t, term in enumerate(fit.terms):
        if "c" not in term.kinds:
            codes.append(fit.coded[t].codes)
            nlev.append(fit.coded[t].nlev)
            continue
        comps = []
        for v, kind in zip(term.vars, term.kinds):
            col = df[v]
            if kind == "c":
                comps.append(bin_codes(col.to_numpy(dtype=float)[fit.idx], cuts[v]).astype(float))
            elif isinstance(col.dtype, pd.CategoricalDtype):
                comps.append(col.cat.codes.to_numpy()[fit.idx].astype(float))
            elif is_string_series(col):
                comps.append(pd.factorize(col.to_numpy()[fit.idx], sort=True)[0].astype(float))
            else:
                comps.append(col.to_numpy(dtype=float)[fit.idx])
        X = np.column_stack(comps)
        uniq, inv = np.unique(X, axis=0, return_inverse=True)
        codes.append(np.asarray(inv).reshape(-1) + 1)
        nlev.append(len(uniq))
    return codes, nlev


class _Runner:
    """Computes the row of a replicate (the coefficient vector and the explained shares) for
    a weight vector over the estimation sample."""

    def __init__(self, ctx, fit, sample_idx):
        self.ctx = ctx
        self.fit = fit
        self.chain = getattr(fit, "is_chain", False)
        self.fits = fit.fits if self.chain else [fit]
        self.pos = [np.searchsorted(sample_idx, f.idx) for f in self.fits]
        self.force = bool(ctx.opt["force"])
        self.nst = len(self.fits[0].stages)
        self.staged = self.fits[0].mode != ""
        self.obs_rows = None

    def pair(self, j, w, state):
        f = self.fits[j]
        wj = w[self.pos[j]]
        st = {"frc": 0, "itl": 0, "itst": state.setdefault("itst", []), "itrev": 0}
        if self.ctx.recut and self.ctx.contlist:
            try:
                st["codes"] = recode(f, self.ctx, wj)
            except SeldecompError:
                raise RepFail("rk")
        try:
            vals, pct = rep_values(f, self.ctx, wj, self.force, st)
        finally:
            state["itrev"] = max(state.get("itrev", 0), st["itrev"])
        state["frc"] = max(state["frc"], st["frc"])
        state["itl"] = max(state["itl"], st["itl"])
        return vals, pct

    def row(self, w, state, only=None):
        if not self.chain:
            return self.pair(0, w, state)
        rows, pcts = [], []
        for j in range(len(self.fits)):
            if only is not None and j not in only:
                v, p = self.obs_rows[j]
            else:
                try:
                    v, p = self.pair(j, w, state)
                except RepFail as ex:
                    ex.pair = j + 1
                    raise
            rows.append(v)
            pcts.append(p[0])
        vals = [x for v in rows for x in v]
        tot = [sum(v[3] for v in rows), sum(v[4] for v in rows), sum(v[5] for v in rows)]
        vals += tot
        if self.staged:
            vals += [sum(v[6 + s] for v in rows) for s in range(self.nst)]
        pcts.append(100 * tot[1] / tot[0] if tot[0] != 0 else NAN)
        return vals, pcts


def _summ(col, ob, level, df, ties, accel):
    """One column of replicates: (mean, bias, se, z0, normal, t, percentile, bc, bca). The
    steps, and the numbers carried from one step to the next as a Stata macro holds them,
    are those of the Stata program, so that a percentile that falls between two replicates
    falls between the same two."""
    v = col[~np.isnan(col)]
    n = len(v)
    alo = M((100 - level) / 200.0)
    ahi = M(1 - alo)
    zz = M(invnormal(1 - (100 - level) / 200.0))
    tt = M(stats.t.ppf(1 - (100 - level) / 200.0, df)) if df > 0 else NAN
    mean = float(v.mean())
    se = M(float(v.std(ddof=1)))
    out = {"b_bs": mean, "bias": mean - ob, "se": se, "z0": NAN,
           "normal": (ob - zz * se, ob + zz * se), "t": (ob - tt * se, ob + tt * se),
           "percentile": tuple(pctile_plain(v, p=[100 * alo, 100 * ahi])),
           "bc": (NAN, NAN), "bca": (NAN, NAN)}
    nle = int(np.sum(v <= ob))
    neq = int(np.sum(v == ob))
    p0 = M((nle - 0.5 * neq) / n) if ties else M(nle / n)
    if 0 < p0 < 1:
        z0 = M(invnormal(p0))
        out["z0"] = z0
        pl, pu = M(normal(2 * z0 - zz)), M(normal(2 * z0 + zz))
        if 0 < pl < 1 and 0 < pu < 1 and pl < pu:
            out["bc"] = tuple(pctile_plain(v, p=[100 * pl, 100 * pu]))
        if accel is not None and not np.isnan(accel):
            a = M(accel)
            pl = M(normal(z0 + (z0 - zz) / (1 - a * (z0 - zz))))
            pu = M(normal(z0 + (z0 + zz) / (1 - a * (z0 + zz))))
            if 0 < pl < 1 and 0 < pu < 1 and pl < pu:
                out["bca"] = tuple(pctile_plain(v, p=[100 * pl, 100 * pu]))
    return out


def acceleration(JK):
    """The BCa acceleration of each column of jackknife values."""
    a = np.full(JK.shape[1], NAN)
    for j in range(JK.shape[1]):
        v = JK[:, j]
        v = v[~np.isnan(v)]
        if len(v) < 2:
            continue
        d = v.mean() - v
        den = 6 * (np.sum(d ** 2)) ** 1.5
        a[j] = np.sum(d ** 3) / den if den > 0 else 0.0
    return a


def summarize(draws, pdraws, ob, pob, level, df, ties, accel):
    """Summary of the completed replicates: a dict of arrays, and the explained shares'
    standard errors and percentile intervals."""
    ok = ~np.isnan(draws).any(axis=1)
    D = draws[ok]
    kb = D.shape[1]
    keys = ("b_bs", "bias", "se", "z0")
    S = {k: np.full(kb, NAN) for k in keys}
    for t in CI_TYPES:
        S["ci_" + t] = np.full((2, kb), NAN)
    for j in range(kb):
        r = _summ(D[:, j], ob[j], level, df, ties, None if accel is None else accel[j])
        for k in keys:
            S[k][j] = r[k]
        for t in CI_TYPES:
            S["ci_" + t][:, j] = r[t]
    S["V"] = np.cov(D, rowvar=False, ddof=1) if kb > 1 else np.array([[D[:, 0].var(ddof=1)]])
    alo = M((100 - level) / 200.0)
    ahi = M(1 - alo)
    pse, pci = [], []
    for k in range(pdraws.shape[1]):
        if pob[k] is None or np.isnan(pob[k]):
            pse.append(NAN)
            pci.append((NAN, NAN))
            continue
        v = pdraws[ok, k]
        pse.append(float(v.std(ddof=1)))
        pci.append(tuple(pctile_plain(v, p=[100 * alo, 100 * ahi])))
    return S, pse, pci


class _Dots:
    """The progress display that dots(#) asks for: one dot every # replicates (an x for a
    failed one), fifty to a line."""

    def __init__(self, title, n, every):
        self.every = max(1, every)
        self.count = 0
        print("")
        print("%s (%d)" % (title, n))
        print("----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5 ")

    def tick(self, i, failed):
        if i % self.every:
            return
        self.count += 1
        end = "%6d\n" % i if self.count % 50 == 0 else ""
        print("x" if failed else ".", end=end, flush=True)

    def close(self):
        if self.count % 50:
            print("")


def bootstrap(ctx, fit, res):
    bs = ctx.bs
    o = ctx.opt
    e, M = res.e, res.matrices
    chain = getattr(fit, "is_chain", False)
    sample = res.sample.to_numpy()
    sample_idx = np.flatnonzero(sample)
    n = len(sample_idx)
    w = ctx.w[sample_idx]
    reps = bs["reps"]
    level = ctx.level
    # sampling units: the observations, or the clusters
    N_clust = None
    if bs["cluster"]:
        cl = ctx.df[bs["cluster"]].to_numpy()[sample_idx]
        units, unit_of = np.unique(cl, return_inverse=True)
        m = len(units)
        N_clust = m
        if m < 3:
            raise SeldecompError("bayesboot(cluster()): fewer than 3 clusters in the estimation sample", 459)
    else:
        unit_of = np.arange(n)
        m = n
        if m < 3:
            raise SeldecompError("bayesboot(): fewer than 3 observations", 459)
    # saving(filename [, replace]): checked before any replication runs
    savefile = None
    if bs["saving"]:
        from ._parse import bs_saving_file, check_bs_saving
        check_bs_saving(bs["saving"])
        savefile = bs_saving_file(bs["saving"])[0]
    seed = bs["rseed"]
    if seed is None:
        seed = int(np.random.SeedSequence().generate_state(1)[0] % 2147483647) + 1
    rng = np.random.default_rng(seed)
    shape = (m - 2) / m
    W = float(w.sum())
    given = bs.get("_weights")
    if given is not None:
        given = np.asarray(given, dtype=float)
        if given.shape != (n, reps):
            raise SeldecompError("the replicate weights given must have one row per observation of the "
                                 "sample and one column per replicate", 198)
    run = _Runner(ctx, fit, sample_idx)
    names = list(M["b"].columns)
    ob = np.asarray(M["b"].iloc[0].to_numpy(dtype=float))
    kb = len(names)
    if chain:
        pob = [f.E.res[0][len(f.stages) - 1, 3] for f in fit.fits] + [e["pct_explained"]]
        pnames = ["p%d:pct_explained" % (j + 1) for j in range(len(fit.fits))] + ["total:pct_explained"]
    else:
        nst = len(fit.stages)
        pob = [fit.E.res[j][nst - 1, 3] for j in range(len(ctx.ylist))]
        pnames = ["pct_explained"] if len(ctx.ylist) == 1 else ["%s:pct_explained" % y for y in ctx.ylist]
    npc = len(pnames)
    draws = np.full((reps, kb), NAN)
    pdraws = np.full((reps, npc), NAN)
    # the replicate weights are kept only when generate() asks for them
    bw_all = np.full((n, reps), NAN) if bs["generate"] else None
    nfail = niter = nforce = 0
    nf = {"it": 0, "tg": 0, "tb": 0, "rk": 0}
    fpairs = []
    itst_all = []
    itrev = 0
    dots = None
    if bs.get("_dots") and bs["dots"] > 0 and not bs["nodots"]:
        dots = _Dots("Bayesian bootstrap replications", reps, bs["dots"])
    for r in range(reps):
        if given is not None:
            bw = given[:, r]
        else:
            g = rng.gamma(shape, 1.0, size=m)[unit_of]
            wg = w * g
            bw = wg * W / float(wg.sum())
        if bw_all is not None:
            bw_all[:, r] = bw
        state = {"frc": 0, "itl": 0, "itst": itst_all, "itrev": 0}
        try:
            vals, pct = run.row(bw, state)
        except RepFail as ex:
            itrev = max(itrev, state["itrev"])
            nfail += 1
            nf[ex.why] += 1
            if chain and getattr(ex, "pair", None) and ex.pair not in fpairs:
                fpairs.append(ex.pair)
            if dots:
                dots.tick(r + 1, True)
            continue
        draws[r, :] = vals
        pdraws[r, :] = pct
        niter += state["itl"]
        nforce += state["frc"]
        itrev = max(itrev, state["itrev"])
        if dots:
            dots.tick(r + 1, False)
    if dots:
        dots.close()
    why = ""
    if nf["it"] > 0 and nf["tg"] == 0:
        why = "raking stopped at iterate() short of its targets"
    elif nf["tg"] > 0:
        why = "raking did not reach its targets"
        if nf["it"] > 0:
            why += " (%d stopped at iterate())" % nf["it"]
    if nf["tb"] > 0:
        why += (", or " if why else "") + "the trimming bounds contradicted each other"
    if nf["rk"] > 0:
        why += (", or " if why else "") + "raking failed"
    if chain and fpairs:
        why += " (%s %s)" % ("pair" if len(fpairs) == 1 else "pairs", ", ".join("p%d" % p for p in sorted(fpairs)))
    N_reps = reps - nfail
    if N_reps < 2:
        raise SeldecompError("bayesboot(): %d of %d replications succeeded; 2 or more are required. "
                             "Failures: %s" % (N_reps, reps, why), 430)

    # the jackknife of the BCa acceleration: one deletion per sampling unit
    accel = None
    njfail = njforce = 0
    if bs["bca"]:
        if chain:
            run.obs_rows = []
            for j in range(len(fit.fits)):
                run.obs_rows.append(run.pair(j, w, {"frc": 0, "itl": 0}))
            inc = [np.unique(unit_of[run.pos[j]]) for j in range(len(fit.fits))]
        JK = np.full((m, kb + npc), NAN)
        jdots = _Dots("Jackknife for BCa acceleration", m, bs["dots"]) if dots else None
        for u in range(m):
            jw = np.where(unit_of == u, 0.0, w)
            state = {"frc": 0, "itl": 0}
            only = None
            if chain:
                only = [j for j in range(len(fit.fits)) if u in inc[j]]
            try:
                vals, pct = run.row(jw, state, only=only)
            except RepFail:
                njfail += 1
                if jdots:
                    jdots.tick(u + 1, True)
                continue
            JK[u, :kb] = vals
            JK[u, kb:] = pct
            njforce += state["frc"]
            if jdots:
                jdots.tick(u + 1, False)
        if jdots:
            jdots.close()
        accel = acceleration(JK)

    df_bs = (N_clust if N_clust is not None else n) - 1
    S, pse, pci = summarize(draws, pdraws, ob, pob, level, df_bs, bs["ties"],
                            None if accel is None else accel[:kb])
    M["V"] = pd.DataFrame(S["V"], index=names, columns=names)
    for k in ("b_bs", "bias", "se", "z0"):
        M[k] = pd.DataFrame([S[k]], index=["r1"], columns=names)
    M["reps"] = pd.DataFrame([[float(N_reps)] * kb], index=["r1"], columns=names)
    for t in CI_TYPES:
        if t == "bca" and not bs["bca"]:
            continue
        M["ci_" + t] = pd.DataFrame(S["ci_" + t], index=["ll", "ul"], columns=names)
    if bs["bca"]:
        M["accel"] = pd.DataFrame([accel[:kb]], index=["r1"], columns=names)
    rn = ["r%d" % (i + 1) for i in range(reps)]
    M["bsdraws"] = pd.DataFrame(draws, index=rn, columns=names)
    if npc > 1 or not np.isnan(pob[0]):
        M["bsdraws_pct"] = pd.DataFrame(pdraws, index=rn, columns=pnames)
    if chain:
        J = len(fit.fits)
        M["pct_explained_ci"] = pd.DataFrame([list(pci[J])], index=["r1"], columns=["ll", "ul"])
        bp = M["by_pair"]
        bp["pct_explained_se"] = pse[:J]
        bp["pct_explained_ll"] = [c[0] for c in pci[:J]]
        bp["pct_explained_ul"] = [c[1] for c in pci[:J]]
        e["pct_explained_se"] = pse[J]
    else:
        M["pct_explained_ci"] = pd.DataFrame([list(pci[0])], index=["r1"], columns=["ll", "ul"])
        e["pct_explained_se"] = pse[0]
        if npc > 1:
            bo = M["by_outcome"]
            bo["pct_explained_se"] = pse
            bo["pct_explained_ll"] = [c[0] for c in pci]
            bo["pct_explained_ul"] = [c[1] for c in pci]
    e["N_reps"] = N_reps
    e["N_misreps"] = nfail
    e["N_iterlimit"] = niter
    if o["force"]:
        e["N_forced"] = nforce
    e["seed"] = seed
    e["bca"] = 1 if bs["bca"] else 0
    e["ties"] = 1 if bs["ties"] else 0
    if bs["bca"]:
        e["N_jkfail"] = njfail
    if N_clust is not None:
        e["N_clust"] = N_clust
        e["clustvar"] = bs["cluster"]
    e["df_bs"] = df_bs
    types = ci_types(bs)
    # the tests of the coefficient table: the first of normal and t among the types shown;
    # with neither, t when the replicates are drawn by cluster
    nt = [t for t in types if t in ("normal", "t")]
    if (nt[0] if nt else ("t" if N_clust is not None else "normal")) == "t":
        e["df_r"] = df_bs
    e["level"] = level
    e["vce"] = "bayesboot"
    e["vcetype"] = "Bootstrap"
    e["bsnames"] = " ".join(names)
    e["citypes"] = " ".join(types)
    if nfail > 0:
        e["missing"] = "missing"
    # replications that stopped at the iterate() limit where the point estimate did: its note
    # counts them (one message per cause), unless the reverse order or a term stage stopped too
    nstop = nf["it"] + niter
    inote = None if chain else getattr(fit, "inote", None)
    if (inote is not None and nfail == nf["it"] and nstop > 0 and not itrev
            and all(st in fit.itstages for st in itst_all)):
        dtx = " (%d short of the targets, dropped)" % nfail if nfail > 0 else ""
        caution = "; bootstrap results may be incorrect." if nfail > 0 else ""
        res.notes[inote] = ("%sraking stopped at the iterate() limit, also in %d of %d replications%s; "
                            "raise iterate()%s" % (fit.istx, nstop, reps, dtx, caution))
        e["note%d" % (inote + 1)] = res.notes[inote]
    if bw_all is not None:
        res.bsweights = pd.DataFrame(bw_all, index=ctx.df.index[sample_idx],
                                     columns=["bw_%d" % (r + 1) for r in range(reps)])
    res.bsdraws = M["bsdraws"]
    if bs["generate"]:
        parts = [x.strip() for x in str(bs["generate"]).split(",")]
        stub = parts[0]
        if len(stub.split()) != 1 or (len(parts) > 1 and parts[1] != "replace") or len(parts) > 2:
            raise SeldecompError("bayesboot(generate()): type generate(stub) or generate(stub, replace)", 198)
        replace = len(parts) > 1
        for r in range(reps):
            nm = "%s_%d" % (stub, r + 1)
            if nm in ctx.df.columns and not replace:
                raise SeldecompError("variable %s already defined" % nm, 110)
            col = np.full(len(ctx.df), NAN)
            col[sample_idx] = bw_all[:, r]
            res._generated[nm] = col
        e["bsgenerate"] = stub
    if bs["saving"]:
        out = M["bsdraws"].copy()
        if "bsdraws_pct" in M:
            for c in M["bsdraws_pct"].columns:
                out[c] = M["bsdraws_pct"][c]
        out.to_csv(savefile, index=False)
        e["bsfile"] = savefile
    res._pct = {"names": pnames, "se": pse, "ci": pci, "draws": pdraws}
    res._bsinfo = {"N_fail": nfail, "why": why, "N_iter": niter, "N_force": nforce, "N_jkforce": njforce,
                   "N_jkfail": njfail, "nf": nf, "reps": reps}

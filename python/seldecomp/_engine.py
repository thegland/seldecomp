"""The decomposition for one weight vector, all stages."""
import numpy as np

from ._rake import rake, TRIM_ONCE, TRIM_SOMETIMES, TRIM_OFTEN
from ._stata import plural

# columns of the stage matrix (e(bystage) in Stata)
STAGE_COLS = [
    "raked_focal", "cum_expl", "cum_unexpl", "pct_explained", "increment",
    "converged", "badcontrols", "maxreldif", "maxctrl",
    "raked_mean", "raked_min", "raked_max", "raked_sd", "raked_cv",
    "factor_mean", "factor_min", "factor_max", "factor_sd", "factor_cv",
    "N_raked", "N_focal_zero", "N_at_bound", "N_eff_base", "N_eff_raked", "ess_ratio",
    "N_focal_mincell", "pct_increment", "pct_W_focal_zero", "mean_focal_zero", "mean_focal_cs",
]
NCOL = len(STAGE_COLS)
NAN = float("nan")


class RakeOptions:
    """The raking options of one run."""

    def __init__(self, iterate=10000, tolerance=1e-8, ctrltolerance=1e-6, trimhirel=None,
                 trimlorel=None, trimhiabs=None, trimloabs=None, trimfrequency="sometimes",
                 trimsupportoverride=False, mincell=10):
        self.iterate = iterate
        self.tolerance = tolerance
        self.ctrltolerance = ctrltolerance
        self.trimhirel = trimhirel
        self.trimlorel = trimlorel
        self.trimhiabs = trimhiabs
        self.trimloabs = trimloabs
        self.trimfrequency = trimfrequency
        self.trimsupportoverride = trimsupportoverride
        self.mincell = mincell

    @property
    def hastrim(self):
        return any(b is not None for b in (self.trimhirel, self.trimlorel, self.trimhiabs, self.trimloabs))

    @property
    def tfreq(self):
        if not self.hastrim:
            return 0
        return {"once": TRIM_ONCE, "often": TRIM_OFTEN}.get(self.trimfrequency, TRIM_SOMETIMES)

    @property
    def incl(self):
        """Focal observations with a level absent from the reference group stay in the raking
        (a lower bound holds them) unless trimsupportoverride sets them aside."""
        return (not self.trimsupportoverride) and (self.trimloabs is not None or self.trimlorel is not None)


class Design:
    """The coded estimation sample of one decomposition.

    g1: True for the focal group. ys: the outcomes (the first is the dependent variable).
    den: the denominator of ratio(), or None. codes[t]: the level (1..nlev[t]) of term t.
    """

    def __init__(self, g1, ys, den, codes, nlev, tnames):
        self.g1 = np.asarray(g1, dtype=bool)
        self.ys = [np.asarray(y, dtype=float) for y in ys]
        self.den = None if den is None else np.asarray(den, dtype=float)
        self.codes = [np.asarray(c, dtype=np.int64) for c in codes]
        self.nlev = list(nlev)
        self.tnames = list(tnames)
        self.n = len(self.g1)


class EngineFail(Exception):
    """The decomposition could not be computed for these weights."""

    def __init__(self, rc, msg="", **info):
        Exception.__init__(self, msg)
        self.rc = rc
        self.msg = msg
        self.info = info


class EngineOut:
    pass


def _msd(x):
    """mean, min, max, sd (n - 1), cv (sd / mean when the mean is positive)"""
    n = len(x)
    m = float(np.sum(x) / n)
    sd = float(np.sqrt(np.sum((x - m) ** 2) / (n - 1))) if n > 1 else NAN
    return [m, float(x.min()), float(x.max()), sd, (sd / m if m > 0 else NAN)]


def _wsum(x, w):
    return float(np.sum(x * w))


def _wmean(x, w, total):
    """The weighted mean of x (total: the sum of the weights). The mean of values that are
    all the same is that value itself, as Stata's summarize gives it."""
    if len(x) and x.min() == x.max():
        return float(x[0])
    return _wsum(x, w) / total


def engine(D, w, stages, opts, holdothers=False, contrib=False, skipdiag=False,
           want_weights=False, want_misses=False, trace=None):
    """Rake the focal group at every stage and decompose the difference of each outcome.

    stages: a list of lists of term numbers (0-based), the terms raked to the reference
    group at each stage. holdothers: the other terms are held at the focal group's own
    shares. contrib: the increment of a stage is its own explained part, and that of the last
    stage the explained part minus their sum (individually).
    """
    w = np.asarray(w, dtype=float)
    g1 = D.g1
    g0 = ~g1
    ny = len(D.ys)
    nterms = len(D.codes)
    nstages = len(stages)
    out = EngineOut()

    W1 = float(np.sum(w[g1]))
    W2 = float(np.sum(w[g0]))
    m1 = [(_wmean(y[g1], w[g1], W1) if W1 > 0 else NAN) for y in D.ys]
    m2 = [(_wmean(y[g0], w[g0], W2) if W2 > 0 else NAN) for y in D.ys]
    if W1 <= 0 or W2 <= 0 or np.isnan(m1[0]) or np.isnan(m2[0]):
        raise EngineFail(2000, "one of the groups has zero total weight")
    n1 = n2 = d1 = d2 = None
    if D.den is not None:
        d1 = _wmean(D.den[g1], w[g1], W1)
        d2 = _wmean(D.den[g0], w[g0], W2)
        if d1 == 0 or d2 == 0 or np.isnan(d1) or np.isnan(d2):
            which = "Focal" if (d1 == 0 or np.isnan(d1)) else "Reference"
            raise EngineFail(459, "ratio(): the mean of %s is zero in the %s Group" % ("{den}", which))
        n1, n2 = list(m1), list(m2)
        m1 = [a / d1 for a in n1]
        m2 = [a / d2 for a in n2]
    diff = [a - b for a, b in zip(m1, m2)]

    rs = g1 & (w > 0)
    idx_rs = np.flatnonzero(rs)
    w_rs = w[idx_rs]
    codes_rs = [c[idx_rs] for c in D.codes]
    # the reference group's share of each level of each term
    s2 = g0 & (w > 0)
    P2 = [np.bincount(D.codes[t][s2], weights=w[s2], minlength=D.nlev[t] + 1)[1:] / W2
          for t in range(nterms)]
    # a level with no reference observation has share exactly 0
    present2 = [np.bincount(D.codes[t][s2], minlength=D.nlev[t] + 1)[1:] > 0 for t in range(nterms)]
    for t in range(nterms):
        P2[t] = np.where(present2[t], P2[t], 0.0)

    res = [np.full((nstages, NCOL), NAN) for _ in range(ny)]
    checks = np.full((nstages, 5), NAN)
    misses = []
    zero_notes, thin_notes, held_notes = [], [], []
    trimchg = [0] * nstages
    supmiss = [0] * nstages
    worst = [(None, None)] * nstages
    rw_stages = []
    prevexpl = [0.0] * ny
    n1r = [NAN] * ny
    d1r = NAN
    tb_user = (opts.trimhirel, opts.trimlorel, opts.trimhiabs, opts.trimloabs)
    incl = opts.incl

    for s, st in enumerate(stages):
        st = list(st)
        ot = [t for t in range(nterms) if t not in st] if holdothers else []
        # focal observations with a level of a stage term that the reference group lacks
        z = np.zeros(len(idx_rs), dtype=bool)
        zerolist = []
        for t in st:
            c = codes_rs[t]
            cnt = np.bincount(c, minlength=D.nlev[t] + 1)[1:]
            for l in np.flatnonzero((P2[t] == 0) & (cnt > 0)):
                z |= (c == l + 1)
                zerolist.append((t, int(l) + 1, int(cnt[l])))
        nz = int(z.sum())
        keep = np.ones(len(idx_rs), dtype=bool) if incl else ~z
        Nrake = int(keep.sum())
        Nin = len(idx_rs) - nz
        nsup = nz
        nzero = 0 if incl else nz
        if Nin == 0:
            # the terms whose absent levels set every focal observation aside
            sole = None
            for t in st:
                c = codes_rs[t]
                if np.all(P2[t][c - 1] == 0):
                    sole = t
                    break
            raise EngineFail(2000, "", stage=s, zerolist=zerolist, nsup=nsup, septerm=sole)
        margins_spec = []
        nthin = 0
        thinlist = []
        fails = []
        for t in st:
            c = codes_rs[t]
            cnt_keep = np.bincount(c[keep], minlength=D.nlev[t] + 1)[1:]
            cnt_all = np.bincount(c, minlength=D.nlev[t] + 1)[1:]
            pres = np.flatnonzero(cnt_keep > 0)
            for l in pres:
                if cnt_keep[l] < opts.mincell and P2[t][l] > 0:
                    nthin += 1
                    thinlist.append((t, int(l) + 1, int(cnt_keep[l])))
            for l in np.flatnonzero((cnt_keep == 0) & (P2[t] > 0)):
                fails.append((t, int(l) + 1, bool(cnt_all[l] > 0)))
            if not fails and len(pres) > 1:
                sh = P2[t][pres]
                psum = 0.0
                for v in sh:
                    psum += v
                margins_spec.append((t, pres + 1, W1 * sh / psum))
        if fails:
            raise EngineFail(459, "", stage=s, fails=fails, nzero=nzero, zerolist=zerolist)
        ir = idx_rs[keep]
        wr = w_rs[keep]
        Wrss = float(np.sum(wr))
        for t in ot:
            c = codes_rs[t][keep]
            tot = np.bincount(c, weights=wr, minlength=D.nlev[t] + 1)[1:]
            cnt = np.bincount(c, minlength=D.nlev[t] + 1)[1:]
            pres = np.flatnonzero(cnt > 0)
            if len(pres) > 1:
                sh = tot[pres] / Wrss
                psum = 0.0
                for v in sh:
                    psum += v
                margins_spec.append((t, pres + 1, W1 * sh / psum))
        sc = Nrake / Wrss
        ntrim = 0
        tchg = 0
        badtargets = 0
        ntargets = 0
        bad0 = 0
        if not margins_spec and opts.hastrim:
            # constant controls still constrain the total weight; with trimming the total goes
            # through the raking so that the bounds are applied and checked
            t0 = st[0]
            lev0 = int(codes_rs[t0][keep].min())
            margins_spec = [(t0, np.array([lev0]), np.array([W1]))]
        if not margins_spec:
            rw_k = wr * W1 / Wrss
            conv, badc, mrd, mctrl = 1, 0, 0.0, 0.0
            wt_, wl_ = None, None
            raw_conv = 1
        else:
            ws = wr * sc
            margins = []
            for t, levs, targets in margins_spec:
                pos = np.full(D.nlev[t] + 1, -1, dtype=np.int64)
                pos[levs] = np.arange(len(levs))
                idx = pos[codes_rs[t][keep]]
                if np.any(idx < 0):
                    # (a constant control: only the observations of its level are adjusted)
                    idx = np.where(idx < 0, 0, idx)
                margins.append((idx, np.asarray(levs), np.asarray(targets, dtype=float) * sc))
            tb = (opts.trimhirel, opts.trimlorel,
                  None if opts.trimhiabs is None else opts.trimhiabs * sc,
                  None if opts.trimloabs is None else opts.trimloabs * sc)
            if trace and nstages > 1:
                trace("")
                trace("Stage %d" % (s + 1))
            r = rake(ws, margins, iterate=opts.iterate, tol=opts.tolerance, ctrltol=opts.ctrltolerance,
                     tb=tb, tfreq=opts.tfreq, trace=trace, want_misses=want_misses)
            if r.tclash[0] + r.tclash[1] > 0:
                parts = []
                if r.tclash[0] > 0:
                    parts.append("trimloabs() is above trimhirel() times the base weight of %d focal group %s"
                                 % (r.tclash[0], plural(r.tclash[0], "observation")))
                if r.tclash[1] > 0:
                    parts.append("trimlorel() times the base weight of %d focal group %s is above trimhiabs()"
                                 % (r.tclash[1], plural(r.tclash[1], "observation")))
                raise EngineFail(198, "stage %d: the trimming bounds contradict each other: %s"
                                 % (s + 1, ", and ".join(parts)), trimclash=True, stage=s)
            conv = raw_conv = r.converged
            badc = r.badcontrols
            bad0 = r.bad0
            mrd = r.maxreldif
            mctrl = r.maxctrl
            badtargets = r.badtargets
            ntargets = r.ntargets
            ntrim = r.n_at_bound
            for k, lev, target, achieved, disc in r.misses:
                misses.append((s, margins_spec[k][0], int(lev), target / sc, achieved / sc, disc))
            wt_ = margins_spec[r.worst_margin][0] if r.worst_margin is not None else None
            wl_ = int(r.worst_level) if r.worst_level is not None else None
            supmiss[s] = bad0
            if badc > 0 and opts.hastrim and bad0 == 0:
                # the bounds are the cause of a miss only when trimming changed the weights:
                # compare with the raking of the stage without the bounds
                ru = rake(ws, margins, iterate=opts.iterate, tol=opts.tolerance, ctrltol=opts.ctrltolerance)
                fd = np.abs(r.w - ru.w) / ws
                tchg = 1 if float(fd.max()) > opts.tolerance else 0
                if trace:
                    from ._display import para
                    for line in para("stage %d: trimming changed the raked weights: %s (largest change of a "
                                     "factor %.2e)" % (s + 1, "yes" if tchg else "no", float(fd.max())), 0, 4):
                        trace(line)
            rw_k = r.w / sc
        trimchg[s] = tchg
        worst[s] = (wt_, wl_)
        checks[s, :] = [raw_conv, badc, badtargets, ntargets, sc]

        # raked means and the decomposition
        if D.den is not None:
            d1r = _wmean(D.den[ir], rw_k, float(np.sum(rw_k)))
            if d1r == 0 or np.isnan(d1r):
                raise EngineFail(459, "ratio(): stage %d: the mean of {den} under the raked focal group weights is zero" % (s + 1))
        srw = float(np.sum(rw_k))
        if not skipdiag:
            dg_rw = _msd(rw_k)
            dg_f = _msd(rw_k / wr)
            rwsq = float(np.sum(rw_k ** 2))
            neff_base = Wrss ** 2 / float(np.sum(wr ** 2))
            neff_raked = srw ** 2 / rwsq if rwsq > 0 else NAN
        zsel = z & ~keep
        for j in range(ny):
            m1r = _wmean(D.ys[j][ir], rw_k, srw)
            if D.den is not None:
                n1r[j] = m1r
                m1r = m1r / d1r
            expl = m1[j] - m1r
            unexpl = m1r - m2[j]
            pct = 100 * expl / diff[j] if diff[j] != 0 else NAN
            if not contrib:
                incr = expl - prevexpl[j]
                prevexpl[j] = expl
            elif s < nstages - 1:
                incr = expl
                prevexpl[j] = prevexpl[j] + expl
            else:
                incr = expl - prevexpl[j]
            R = res[j]
            R[s, 0:5] = [m1r, expl, unexpl, pct, incr]
            R[s, 5:9] = [conv, badc, mrd, mctrl]
            if not skipdiag:
                R[s, 9:14] = dg_rw
                R[s, 14:19] = dg_f
                R[s, 22] = neff_base
                R[s, 23] = neff_raked
                R[s, 24] = neff_base / neff_raked
            R[s, 19] = Nrake
            R[s, 20] = nzero
            R[s, 21] = ntrim
            R[s, 25] = nthin
            R[s, 26] = 100 * incr / diff[j] if diff[j] != 0 else NAN
            if skipdiag:
                pass
            elif nzero > 0:
                R[s, 27] = 100 * (W1 - Wrss) / W1
                wz = w_rs[zsel]
                iz = idx_rs[zsel]
                mz = _wmean(D.ys[j][iz], wz, float(np.sum(wz)))
                mcs = _wmean(D.ys[j][ir], wr, Wrss)
                if D.den is not None:
                    dz = _wmean(D.den[iz], wz, float(np.sum(wz)))
                    dcs = _wmean(D.den[ir], wr, Wrss)
                    mz = mz / dz if dz != 0 else NAN
                    mcs = mcs / dcs if dcs != 0 else NAN
                R[s, 28] = mz
                R[s, 29] = mcs
            else:
                R[s, 27] = 0
                R[s, 29] = m1[j]
        if nthin > 0:
            thin_notes.append((s, nthin, thinlist))
        if want_weights:
            full = np.full(D.n, NAN)
            full[idx_rs] = 0.0
            full[ir] = rw_k
            rw_stages.append(full)
        if nzero > 0:
            zero_notes.append((s, nzero, zerolist))
        elif incl and nsup > 0:
            held_notes.append((s, nsup, zerolist))

    out.res = res
    out.checks = checks
    out.misses = misses
    out.m1, out.m2 = m1, m2
    out.W1, out.W2 = W1, W2
    out.n1, out.n2, out.n1r = n1, n2, n1r
    out.d1, out.d2, out.d1r = d1, d2, d1r
    out.zero, out.thin, out.held = zero_notes, thin_notes, held_notes
    out.trimchg, out.supmiss, out.worst = trimchg, supmiss, worst
    out.rw = rw_stages
    return out

"""Raking (iterative proportional fitting) with trimming."""
import numpy as np

from ._stata import EPSFLOAT, reldif

TRIM_ONCE, TRIM_SOMETIMES, TRIM_OFTEN = 1, 2, 3


class RakeResult:
    """What one raking returns: the weights and the checks made on them."""

    __slots__ = ("w", "converged", "badcontrols", "maxreldif", "maxctrl", "worst_margin",
                 "worst_level", "cycles", "ntrim", "tclash", "bad0", "badtargets", "ntargets",
                 "n_at_bound", "misses")


def trim(w, w0, tb):
    """Bound each weight: above by the smaller of w0 * hirel and hiabs, below by the larger
    of w0 * lorel and loabs. tb = (hirel, lorel, hiabs, loabs), None for an absent bound."""
    hirel, lorel, hiabs, loabs = tb
    if hirel is not None or hiabs is not None:
        hi = np.full(len(w), np.inf)
        if hirel is not None:
            hi = w0 * hirel
        if hiabs is not None:
            hi = np.minimum(hi, hiabs)
        w = np.minimum(w, hi)
    if lorel is not None or loabs is not None:
        lo = np.full(len(w), -np.inf)
        if lorel is not None:
            lo = w0 * lorel
        if loabs is not None:
            lo = np.maximum(lo, loabs)
        w = np.maximum(w, lo)
    return w


def bounds_clash(w0, tb):
    """Observations whose bounds contradict: (loabs above w0 * hirel, w0 * lorel above hiabs)."""
    hirel, lorel, hiabs, loabs = tb
    c1 = c2 = 0
    if loabs is not None and hirel is not None:
        c1 = int(np.sum(loabs > w0 * hirel))
    if lorel is not None and hiabs is not None:
        c2 = int(np.sum(w0 * lorel > hiabs))
    return c1, c2


def rake(w, margins, iterate=10000, tol=1e-8, ctrltol=1e-6, tb=(None, None, None, None),
         tfreq=0, trace=None, want_misses=False):
    """Rake the weights w to the margins.

    margins: a list of (idx, levels, targets); idx gives each observation's position among
    the margin's levels (0-based), levels their codes, targets their target totals.
    A cycle visits the margins in order and multiplies the weights of each level by
    target / current total (a level whose weights are all 0 is left alone). With trimming
    bounds the weights are bounded after each margin (often), after each cycle (sometimes) or
    once after the last cycle (once). The raking has converged when no weight changed by more
    than tol in relative terms over a cycle.
    """
    res = RakeResult()
    w = np.asarray(w, dtype=float).copy()
    w0 = w.copy()
    K = len(margins)
    ntargets = sum(len(m[1]) for m in margins)
    res.misses = []
    res.tclash = bounds_clash(w0, tb) if tfreq else (0, 0)
    if res.tclash[0] + res.tclash[1] > 0:
        res.w = None
        res.converged = 0
        res.badcontrols = res.maxreldif = res.maxctrl = float("nan")
        res.worst_margin = res.worst_level = None
        res.cycles = 0
        res.ntrim = (0, 0, 0, 0)
        res.bad0 = 0
        res.badtargets = float("nan")
        res.ntargets = ntargets
        res.n_at_bound = 0
        return res
    if trace:
        trace("Raking: %d margins, %d levels, %d observations" % (K, ntargets, len(w)))
        if tfreq:
            trace("  trimming after each %s" % ("margin" if tfreq == TRIM_OFTEN else
                                               "cycle" if tfreq == TRIM_SOMETIMES else "raking (once)"))
    conv = False
    obj = float("nan")
    ncyc = 0
    for it in range(1, iterate + 1):
        wprev = w.copy()
        for idx, levels, targets in margins:
            tot = np.bincount(idx, weights=w, minlength=len(levels))
            f = np.ones(len(levels))
            nz = tot != 0
            f[nz] = targets[nz] / tot[nz]
            w = w * f[idx]
            if tfreq == TRIM_OFTEN:
                w = trim(w, w0, tb)
        if tfreq == TRIM_SOMETIMES:
            w = trim(w, w0, tb)
        pos = wprev != 0
        d = np.abs(w[pos] - wprev[pos]) / wprev[pos]
        conv = not bool(np.any(d > tol))
        obj = float(d.max()) if len(d) else float("nan")
        ncyc = it
        if trace:
            trace("  cycle %5d: largest relative change of a weight %10.3e" % (it, obj))
        if conv:
            break
    if tfreq == TRIM_ONCE:
        w = trim(w, w0, tb)
    badc = bad0 = badtargets = 0
    mrdmax = 0.0
    wk = wl = None
    for k, (idx, levels, targets) in enumerate(margins):
        achieved = np.bincount(idx, weights=w, minlength=len(levels))
        rd = reldif(achieved, targets)
        mr = float(rd.max())
        badtargets += int(np.sum(rd > ctrltol))
        if mr > ctrltol:
            badc += 1
            if want_misses:
                worst = int(np.flatnonzero(rd == mr)[-1])
                res.misses.append((k, levels[worst], float(targets[worst]), float(achieved[worst]),
                                   float(rd[worst])))
        bad0 += int(np.sum((targets == 0) & (rd > ctrltol)))
        mrdmax = max(mrdmax, mr)
        hit = np.flatnonzero(reldif(rd, mrdmax) < EPSFLOAT)
        if len(hit):
            wk, wl = k, levels[int(hit[-1])]
    hirel, lorel, hiabs, loabs = tb
    ntr = [0, 0, 0, 0]
    atbound = np.zeros(len(w), dtype=bool)
    with np.errstate(divide="ignore", invalid="ignore"):
        if hiabs is not None:
            hit = (w / hiabs) > 1 - 1e-6
            ntr[0] = int(hit.sum())
            atbound |= hit
        if loabs is not None:
            hit = (w / loabs) < 1 + 1e-6
            ntr[1] = int(hit.sum())
            atbound |= hit
        if hirel is not None:
            hit = (w / w0) > (1 - 1e-6) * hirel
            ntr[2] = int(hit.sum())
            atbound |= hit
        if lorel is not None:
            hit = (w / w0) < (1 + 1e-6) * lorel
            ntr[3] = int(hit.sum())
            atbound |= hit
    if trace:
        trace("  %s after %d cycles; margins missed by more than ctrltolerance(): %d"
              % ("converged" if conv else "not converged", ncyc, badc))
        trace("  largest relative difference from a control total %10.3e" % mrdmax)
        if tfreq:
            trace("  weights at a trimming bound: hiabs %d, loabs %d, hirel %d, lorel %d" % tuple(ntr))
    res.w = w
    res.converged = 1 if conv else 0
    res.badcontrols = badc
    res.maxreldif = float("nan") if tfreq == TRIM_ONCE else obj
    res.maxctrl = mrdmax
    res.worst_margin = wk
    res.worst_level = wl
    res.cycles = ncyc
    res.ntrim = tuple(ntr)
    res.bad0 = bad0
    res.badtargets = badtargets
    res.ntargets = ntargets
    res.n_at_bound = int(atbound.sum())
    return res

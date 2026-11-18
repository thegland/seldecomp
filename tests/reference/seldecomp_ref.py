#!/usr/bin/env python3
"""Reference implementation of the seldecomp raking decomposition.

Generates tests/data/seldecomp_test.csv and tests/expected_values.do, a do-file
that stores the expected values as scalars so tests/test_reference.do can assert
that seldecomp reproduces them, and tests/expected_support.do (cases K and L: a
level found only in Group 1 under sequentially fixothermargins).

The staged cases follow seldecomp's two modes of the option sequentially: case F
(sequentially alone) rakes the blocks added so far to Group 2 and leaves the
blocks not yet added free; case G (sequentially fixothermargins, decompose's
hold=True) also rakes the blocks not yet added back to Group 1's own shares.

Only numpy is required.  Weighted quantiles follow Stata's _pctile definition
(no altdef): with sorted values x_(1..n) and cumulative weights W_(i),
percentile p is x_(i) for the first i with W_(i) >= p*W, except that when
W_(i) == p*W exactly the average of x_(i) and x_(i+1) is used.
"""
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data", "seldecomp_test.csv")
EXPECTED = os.path.join(ROOT, "expected_values.do")
EXPECTED_JSON = os.path.join(ROOT, "data", "expected_values.json")
EXPECTED_SUPPORT = os.path.join(ROOT, "expected_support.do")


# ----------------------------------------------------------------------------
# raking (iterative proportional fitting), the algorithm of seldecomp's Mata raker
# ----------------------------------------------------------------------------
def rake(w, cats, targets, iters=2000, tol=1e-10):
    """w: base weights; cats: list of integer-coded arrays (levels 0..K-1);
    targets: list of arrays of target totals per level.  Returns raked w."""
    w = w.astype(float).copy()
    for _ in range(iters):
        wprev = w.copy()
        for c, t in zip(cats, targets):
            cur = np.bincount(c, weights=w, minlength=len(t))
            f = np.ones(len(t))
            nz = cur > 0
            f[nz] = t[nz] / cur[nz]
            w = w * f[c]
        rel = np.abs(w - wprev) / np.where(wprev > 0, wprev, 1.0)
        if rel.max() < tol:
            return w, True
    return w, False


def wmean(y, w):
    return float(np.sum(w * y) / np.sum(w))


def codes(x):
    lv, c = np.unique(x, return_inverse=True)
    return lv, c


def stata_pctile(x, w, q):
    """Cut points for q quantile groups per Stata's _pctile (weighted)."""
    order = np.argsort(x, kind="mergesort")
    xs, ws = x[order], w[order]
    W = ws.sum()
    cw = np.cumsum(ws)
    cuts = []
    for j in range(1, q):
        P = j / q * W
        tol = 1e-9 * W
        # first i with cw[i] >= P up to rounding: the float cumulative sum can land just
        # below an exact tie (0.2 + 0.7 < 0.9), so search with a tolerance (a run on Stata 19
        # confirmed the tie rule: x = 1..4, w = .2 .7 .7 1.1, nq(3) gives 2.5 and 4)
        i = int(np.searchsorted(cw, P - tol, side="left"))
        if abs(cw[i] - P) <= tol and i + 1 < len(xs):
            cuts.append((xs[i] + xs[i + 1]) / 2)
        else:
            cuts.append(xs[i])
    return cuts


def binned(x, cuts):
    # bin j (1-based) if cuts[j-2] < x <= cuts[j-1]; last bin above last cut
    return np.searchsorted(np.array(cuts), x, side="left")


# ----------------------------------------------------------------------------
# decomposition
# ----------------------------------------------------------------------------
def decompose(y, g, w, varlist, focal, hold=False, stages=None):
    """varlist: list of integer-coded arrays (each a raking variable, levels
    0..K-1 over the pooled sample).  stages: list of lists of indices into
    varlist (cumulative blocks); None = one stage with everything.
    hold=True holds the variables not raked at a stage at Group 1's own shares
    (seldecomp's fixothermargins); hold=False leaves them free.
    Returns dict of results per stage."""
    A = g == focal
    B = ~A
    W1 = w[A].sum()
    m1 = wmean(y[A], w[A])
    m2 = wmean(y[B], w[B])
    if stages is None:
        stages = [list(range(len(varlist)))]
    out = []
    prev_expl = 0.0
    for s in stages:
        cats, targets = [], []
        # the Group 1 observations that keep a positive raked weight at this stage: none of
        # their levels of the variables raked to Group 2 here lacks Group 2 weight (a level
        # absent from Group 2 has target 0, so its observations are zeroed); the held shares
        # (hold=True, fixothermargins) are taken over these observations, as seldecomp's
        # Engine does
        keep = np.ones(A.sum(), dtype=bool)
        for j in s:
            c = varlist[j]
            shareB = np.bincount(c[B], weights=w[B], minlength=c.max() + 1)
            keep &= shareB[c[A]] > 0
        for j, c in enumerate(varlist):
            K = c.max() + 1
            if j in s:
                shareB = np.bincount(c[B], weights=w[B], minlength=K) / w[B].sum()
                # support: levels with positive share in B but absent in A
                presentA = np.bincount(c[A], minlength=K) > 0
                if np.any((shareB > 0) & ~presentA):
                    raise ValueError("support failure")
                cats.append(c[A])
                targets.append(W1 * shareB)
            elif hold:
                wk = w[A][keep]
                shareA = np.bincount(c[A][keep], weights=wk, minlength=K) / wk.sum()
                cats.append(c[A])
                targets.append(W1 * shareA)
        wr, conv = rake(w[A], cats, targets)
        m1r = wmean(y[A], wr)
        expl = m1 - m1r
        res = dict(mean_1=m1, mean_2=m2, raked_1=m1r, difference=m1 - m2,
                   explained=expl, unexplained=m1r - m2,
                   pct_explained=100 * expl / (m1 - m2), increment=expl - prev_expl,
                   converged=conv, raked_sum=float(wr.sum()), W1=float(W1),
                   raked_min=float(wr.min()), raked_max=float(wr.max()))
        prev_expl = expl
        out.append(res)
    return out


def ols_ob_explained(y, g, w, X_builder, focal):
    """Two-fold OB explained part with focal-group coefficients."""
    A = g == focal
    B = ~A
    XA, XB = X_builder(A), X_builder(B)
    wa = w[A]
    beta = np.linalg.lstsq((XA * wa[:, None]).T @ XA, (XA * wa[:, None]).T @ y[A], rcond=None)[0]
    xbarA = (XA.T @ wa) / wa.sum()
    xbarB = (XB.T @ w[B]) / w[B].sum()
    return float((xbarA - xbarB) @ beta)


# ----------------------------------------------------------------------------
# data
# ----------------------------------------------------------------------------
def make_data(n=2000, seed=20260924):
    rng = np.random.default_rng(seed)
    group = rng.choice([0, 1], n, p=[0.45, 0.55])
    # x1: 3-level categorical, distribution differs by group
    x1 = np.where(group == 0, rng.choice([1, 2, 3], n, p=[0.5, 0.3, 0.2]),
                  rng.choice([1, 2, 3], n, p=[0.2, 0.3, 0.5]))
    # x2: binary, correlated with x1 within group
    p2 = 0.3 + 0.15 * (x1 - 1) + 0.15 * group
    x2 = (rng.uniform(size=n) < p2).astype(int)
    # x3: continuous
    x3 = rng.normal(0, 1, n) + 0.4 * (x1 - 2) + 0.5 * group
    # x4: string categorical with 4 levels (coded as letters)
    x4c = rng.choice([0, 1, 2, 3], n, p=[0.4, 0.3, 0.2, 0.1])
    x4c = np.where((group == 1) & (rng.uniform(size=n) < 0.3), 3, x4c)
    x4 = np.array(["north", "south", "east", "west"])[x4c]
    # continuous outcome with nonlinearity and interaction
    y = (1.0 * x1 + 2.0 * x2 + 1.5 * (x1 == 3) * x2 + 0.8 * x3 + 0.3 * x3 ** 2
         + 0.5 * x4c + 0.7 * group + 0.4 * group * (x1 == 2) + rng.normal(0, 1, n))
    ybin = (y + rng.logistic(0, 1, n) > np.median(y)).astype(int)
    wt = np.round(rng.uniform(0.5, 3.0, n), 4)
    cluster = rng.integers(1, 61, n)
    strata = 1 + (cluster - 1) // 20
    return dict(id=np.arange(1, n + 1), group=group, y=y, ybin=ybin, x1=x1, x2=x2,
                x3=x3, x4=x4, wt=wt, cluster=cluster, strata=strata)


def main():
    d = make_data()
    n = len(d["id"])
    with open(DATA, "w", newline="") as f:
        wr = csv.writer(f)
        keys = list(d.keys())
        wr.writerow(keys)
        for i in range(n):
            row = []
            for k in keys:
                v = d[k][i]
                if k in ("y", "x3"):
                    row.append(f"{v:.10g}")
                else:
                    row.append(v)
            wr.writerow(row)

    y, g, wt = d["y"], d["group"], d["wt"].astype(float)
    ones = np.ones(n)
    _, c1 = codes(d["x1"])
    _, c2 = codes(d["x2"])
    _, c4 = codes(d["x4"])
    _, c12 = codes(d["x1"] * 10 + d["x2"])
    E = {}

    # A. one categorical, unweighted, focal = 0 (default)
    r = decompose(y, g, ones, [c1], focal=0)[0]
    E["A"] = r
    E["A"]["ob_explained"] = ols_ob_explained(
        y, g, ones, lambda m: np.column_stack([np.ones(m.sum()), c1[m] == 1, c1[m] == 2]), 0)

    # B. one categorical, pweights
    r = decompose(y, g, wt, [c1], focal=0)[0]
    E["B"] = r
    E["B"]["ob_explained"] = ols_ob_explained(
        y, g, wt, lambda m: np.column_stack([np.ones(m.sum()), c1[m] == 1, c1[m] == 2]), 0)

    # C. two categoricals additive (approximate OB)
    r = decompose(y, g, wt, [c1, c2], focal=0)[0]
    E["C"] = r
    E["C"]["ob_explained"] = ols_ob_explained(
        y, g, wt, lambda m: np.column_stack([np.ones(m.sum()), c1[m] == 1, c1[m] == 2, c2[m]]), 0)

    # D. full interaction x1##x2 (exact OB with interacted dummies)
    r = decompose(y, g, wt, [c12], focal=0)[0]
    E["D"] = r
    lv = np.unique(c12)
    E["D"]["ob_explained"] = ols_ob_explained(
        y, g, wt, lambda m: (c12[m][:, None] == lv[None, :]).astype(float), 0)

    # E. focal = 1, string variable x4 + x1, weighted
    r = decompose(y, g, wt, [c4, c1], focal=1)[0]
    E["E"] = r

    # F. sequentially x1, x2, x4 (the blocks not yet added are free: sequentially alone)
    st = decompose(y, g, wt, [c1, c2, c4], focal=0, stages=[[0], [0, 1], [0, 1, 2]])
    for k, r in enumerate(st, 1):
        E[f"F{k}"] = r
    # G. sequentially fixothermargins (the blocks not yet added held at Group 1's shares)
    st = decompose(y, g, wt, [c1, c2, c4], focal=0, hold=True, stages=[[0], [0, 1], [0, 1, 2]])
    for k, r in enumerate(st, 1):
        E[f"G{k}"] = r

    # H. continuous x3 in quartiles (pooled, weighted cut points) + x2
    cuts = stata_pctile(d["x3"], wt, 4)
    b3 = binned(d["x3"], cuts)
    r = decompose(y, g, wt, [b3, c2], focal=0)[0]
    E["H"] = r
    E["H"]["cuts"] = [float(c) for c in cuts]
    # H2: quartiles computed in focal group only (focalgroup)
    A = g == 0
    cutsA = stata_pctile(d["x3"][A], wt[A], 4)
    b3A = binned(d["x3"], cutsA)
    r = decompose(y, g, wt, [b3A, c2], focal=0)[0]
    E["H2"] = r
    E["H2"]["cuts"] = [float(c) for c in cutsA]
    # H3: quartiles in reference group only (the refgroup suboption of contvars())
    B = g == 1
    cutsB = stata_pctile(d["x3"][B], wt[B], 4)
    b3B = binned(d["x3"], cutsB)
    r = decompose(y, g, wt, [b3B, c2], focal=0)[0]
    E["H3"] = r
    E["H3"]["cuts"] = [float(c) for c in cutsB]

    # I. binary outcome, unweighted, x1 x2 x4
    r = decompose(d["ybin"].astype(float), g, ones, [c1, c2, c4], focal=0)[0]
    E["I"] = r

    with open(EXPECTED_JSON, "w") as f:
        json.dump(E, f, indent=1, default=float)

    with open(EXPECTED, "w") as f:
        f.write("* Expected values generated by tests/reference/seldecomp_ref.py -- do not edit\n")
        f.write("* Each scalar is named exp_<case>_<quantity>\n")
        for case, r in E.items():
            for k, v in r.items():
                if k == "cuts":
                    for j, c in enumerate(v, 1):
                        f.write(f"scalar exp_{case}_cut{j} = {c:.15g}\n")
                elif isinstance(v, (bool, np.bool_)):
                    f.write(f"scalar exp_{case}_{k} = {int(v)}\n")
                else:
                    f.write(f"scalar exp_{case}_{k} = {float(v):.15g}\n")
    # K and L: a level found only in Group 1 (x1 = 9 for the Group 1 observations with id <= 60,
    # made in Stata by -replace-), sequentially fixothermargins with the held shares over the
    # observations that
    # keep a positive weight at the stage (Engine's rule); K adds x1z first, L second (so that
    # stage 1 holds x1z at the whole Group 1's shares and stage 2 zeroes level 9). Written to
    # expected_support.do, so that the files above are not rewritten by these cases
    x1z = np.where((g == 0) & (d["id"] <= 60), 9, d["x1"])
    _, c1z = codes(x1z)
    S = {}
    for case, vl in (("K", [c1z, c2, c4]), ("L", [c2, c1z, c4])):
        st = decompose(y, g, wt, vl, focal=0, hold=True, stages=[[0], [0, 1], [0, 1, 2]])
        for k, r in enumerate(st, 1):
            S[f"{case}{k}"] = r
    with open(EXPECTED_SUPPORT, "w") as f:
        f.write("* Expected values generated by tests/reference/seldecomp_ref.py -- do not edit\n")
        f.write("* Cases K and L (a level found only in Group 1, sequentially fixothermargins); scalars exp_<case>_<quantity>\n")
        for case, r in S.items():
            for k, v in r.items():
                if isinstance(v, (bool, np.bool_)):
                    f.write(f"scalar exp_{case}_{k} = {int(v)}\n")
                else:
                    f.write(f"scalar exp_{case}_{k} = {float(v):.15g}\n")
    for case in ("A", "B", "C", "D"):
        print(case, "raking explained", f"{E[case]['explained']:.8f}",
              "OB explained", f"{E[case]['ob_explained']:.8f}",
              "diff", f"{abs(E[case]['explained'] - E[case]['ob_explained']):.2e}")
    for case in ("F1", "F2", "F3", "G1", "G2", "G3", "H", "I"):
        print(case, {k: round(v, 6) for k, v in E[case].items() if isinstance(v, float)})
    return 0


if __name__ == "__main__":
    sys.exit(main())

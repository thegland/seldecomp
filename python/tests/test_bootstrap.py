"""bayesboot(): the draws of the Python version (the cases with Stata's own replicate weights
are in test_results.py)."""
import numpy as np

import common
import seldecomp as sd


def test_seed_reproduces():
    df = common.load_data()
    df = df[df["id"] <= 400].copy()
    a = sd.seldecomp(df, "y x1 x2", by="group", weights="wt", bayesboot="reps(30) rseed(7)")
    b = sd.seldecomp(df, "y x1 x2", by="group", weights="wt", bayesboot="reps(30) rseed(7)")
    c = sd.seldecomp(df, "y x1 x2", by="group", weights="wt", bayesboot="reps(30) rseed(8)")
    assert np.array_equal(a.bsdraws.to_numpy(), b.bsdraws.to_numpy())
    assert not np.array_equal(a.bsdraws.to_numpy(), c.bsdraws.to_numpy())
    assert a.e["seed"] == 7 and a.e["N_reps"] == 30


def test_replicate_weights():
    """Each replicate multiplies the weights by Gamma((m - 2)/m, 1) draws, one per observation
    or per cluster, and keeps their total."""
    df = common.load_data()
    df = df[df["id"] <= 600].copy()
    res = sd.seldecomp(df, "y x1 x2", by="group", weights="wt", bayesboot="reps(200) rseed(3) generate(bw)")
    w = df["wt"].to_numpy()
    bw = res.bsweights.to_numpy()
    assert bw.shape == (600, 200) and "bw_200" in df.columns
    assert np.allclose(bw.sum(axis=0), w.sum())
    assert (bw > 0).all()
    g = bw / w[:, None]
    g = g / g.mean(axis=0)
    m = 600
    shape = (m - 2) / m
    # a Gamma(shape, 1) draw divided by its mean has variance 1/shape
    assert abs(g.var() - 1 / shape) < 0.08
    res = sd.seldecomp(df, "y x1 x2", by="group", weights="wt",
                       bayesboot="reps(50) rseed(4) cluster(cluster) generate(bc)")
    f = res.bsweights.to_numpy() / w[:, None]
    for c in np.unique(df["cluster"].to_numpy())[:10]:
        rows = (df["cluster"] == c).to_numpy()
        assert np.allclose(f[rows], f[rows][0])
    assert res.e["N_clust"] == df["cluster"].nunique() and res.e["df_bs"] == res.e["N_clust"] - 1
    # the replicate weights are kept only when generate() asks for them
    res = sd.seldecomp(df, "y x1 x2", by="group", weights="wt", bayesboot="reps(20) rseed(4)")
    assert res.bsweights is None


def test_standard_errors_near_stata():
    """With 1,000 replications of its own, the standard errors are those of Stata's 1,000
    replications up to simulation error (each has a relative error of about 2 percent)."""
    for name in common.case_names("long"):
        case = common.load_case(name)
        df = common.load_data()
        res = sd.seldecomp(df, case["cmd"].replace("rseed(77)", "rseed(1001)").replace("rseed(78)", "rseed(1002)"))
        se = res.matrices["se"].iloc[0].to_numpy(float)
        stata = np.array(case["matrices"]["se"]["data"][0])
        assert list(res.matrices["se"].columns) == case["matrices"]["se"]["cols"]
        ratio = se / stata
        assert ratio.min() > 0.87 and ratio.max() < 1.13, (name, ratio)
        assert np.allclose(res.b.to_numpy(float), np.array(case["matrices"]["b"]["data"][0]), rtol=1e-6)

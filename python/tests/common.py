"""What the tests share: the test data with the derived variables of the fixture cases, and the
comparison with the results frozen from the Stata program."""
import json
import math
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURES = os.path.join(HERE, "fixtures")
DATA = os.path.join(HERE, "..", "..", "tests", "data", "seldecomp_test.csv")
VALUE_LABELS = {"group": {0: "Group zero", 1: "Group one"}}
LABELS = {"y": "Outcome", "x1": "First factor", "x2": "Second factor"}


def load_data():
    """tests/data/seldecomp_test.csv with the variables that make_fixtures.do adds."""
    df = pd.read_csv(DATA, dtype={"x4": str})
    i = df["id"].to_numpy()
    df["g3"] = i % 3
    df["hrs"] = 20.0 + (i % 25)
    df["earn"] = df["y"] * df["hrs"]
    df["xs"] = np.where((df["group"] == 0) & (i % 40 == 0), 9, df["x1"])
    df["xr"] = np.where((df["group"] == 1) & (i % 50 == 0), 8, df["x1"])
    wz = df["wt"].to_numpy(dtype=float).copy()
    wz[i % 37 == 0] = 0.0
    wz[i % 91 == 0] = np.nan
    df["wz"] = wz
    ym = df["y"].to_numpy(dtype=float).copy()
    ym[i % 53 == 0] = np.nan
    df["ym"] = ym
    df["ta"] = df["x2"]
    df["tb"] = 1 - df["x2"]
    df["xc"] = np.where((df["g3"] == 0) & (i % 30 == 0), 7, df["x1"])
    df["gs"] = np.where(df["group"] == 0, "alpha", "beta")
    df["xj"] = np.where((i == 9) | ((df["group"] == 1) & (i % 20 == 0)), 5, df["x1"])
    hm = df["hrs"].to_numpy(dtype=float).copy()
    hm[i % 61 == 0] = np.nan
    df["hm"] = hm
    cm = df["cluster"].to_numpy(dtype=float).copy()
    cm[i % 71 == 0] = np.nan
    df["cm"] = cm
    # (Stata's round() rounds halves away from zero)
    x3 = df["x3"].to_numpy(dtype=float)
    df["xt"] = np.sign(x3) * np.floor(np.abs(x3) + 0.5)
    return df


def load_case(name, folder=None):
    with open(os.path.join(folder or FIXTURES, name + ".json")) as f:
        d = json.load(f)

    def num(x):
        return float("nan") if x is None else float(x)

    if d["rc"] == 0 and "scalars" in d:
        d["scalars"] = {k: num(v) for k, v in d["scalars"].items()}
        for m in d["matrices"].values():
            m["data"] = [[num(x) for x in row] for row in m["data"]]
    return d


def close(a, b, rtol, atol):
    if isinstance(a, float) and math.isnan(a):
        return isinstance(b, float) and math.isnan(b)
    if isinstance(b, float) and math.isnan(b):
        return False
    return abs(a - b) <= atol + rtol * abs(b)


def compare(res, case, rtol=1e-6, atol=1e-8, skip_macros=(), skip_scalars=(), skip_matrices=()):
    """The differences between a Python result and the Stata fixture, as a list of text lines
    (empty when they agree)."""
    out = []
    for k, v in case["scalars"].items():
        if k in skip_scalars:
            continue
        if k not in res.e:
            out.append("scalar %s missing (Stata %r)" % (k, v))
        elif not close(float(res.e[k]), v, rtol, atol):
            out.append("scalar %s: %r, Stata %r" % (k, res.e[k], v))
    for k, v in case["macros"].items():
        if k in skip_macros or k in ("predict", "properties", "marginsok", "marginsnotok", "cmd2"):
            continue
        if k not in res.e:
            out.append("macro %s missing (Stata %r)" % (k, v[:60]))
        elif str(res.e[k]).strip() != v.strip():
            out.append("macro %s: %r, Stata %r" % (k, str(res.e[k])[:70], v[:70]))
    for k, m in case["matrices"].items():
        if k in skip_matrices:
            continue
        # a chain without missed targets stores a 1 x 1 missing matrix under this name in Stata
        if len(m["data"]) == 1 and len(m["data"][0]) == 1 and math.isnan(m["data"][0][0]) and k not in res.matrices:
            continue
        if k not in res.matrices:
            out.append("matrix %s missing" % k)
            continue
        df = res.matrices[k]
        if list(df.columns) != m["cols"]:
            out.append("matrix %s columns: %r, Stata %r" % (k, list(df.columns)[:6], m["cols"][:6]))
            continue
        if [str(r) for r in df.index] != m["rows"]:
            out.append("matrix %s rows: %r, Stata %r" % (k, list(df.index)[:6], m["rows"][:6]))
            continue
        A = df.to_numpy(dtype=float)
        for i, row in enumerate(m["data"]):
            for j, v in enumerate(row):
                if not close(float(A[i, j]), v, rtol, atol):
                    out.append("matrix %s[%s, %s]: %r, Stata %r" % (k, m["rows"][i], m["cols"][j], float(A[i, j]), v))
                    break
            else:
                continue
            break
    extra_s = [k for k in res.e if k not in case["scalars"] and k not in case["macros"] and k != "cmdline"]
    if extra_s:
        out.append("extra in Python e: %s" % " ".join(sorted(extra_s)[:12]))
    extra_m = [k for k in res.matrices if k not in case["matrices"]]
    if extra_m:
        out.append("extra matrices in Python: %s" % " ".join(sorted(extra_m)))
    return out


def replicate_weights(name, folder=None):
    """The replicate weights that Stata drew for a bayesboot case (one column per replicate,
    the rows of the estimation sample in order), or None."""
    path = os.path.join(folder or FIXTURES, name + "_bsweights.csv.gz")
    if not os.path.exists(path):
        return None
    bw = pd.read_csv(path)
    return bw[[c for c in bw.columns if c != "id"]].to_numpy(dtype=float)


def run_case(name, **kw):
    """Run a fixture case in Python, with Stata's replicate weights when it has them."""
    case = load_case(name)
    df = load_data()
    bw = replicate_weights(name)
    if bw is not None:
        kw["_bsweights"] = bw
    import contextlib
    import io
    import warnings
    import seldecomp as sd
    # (the trace of noisily, printed as the raking runs, is kept in res.trace; a warning is
    # also a line of res.messages)
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore", sd.SeldecompWarning)
        res = sd.seldecomp(df, case["cmd"], value_labels=VALUE_LABELS, labels=LABELS, **kw)
    return case, df, res


def case_names(kind=None):
    """The names of the fixture cases: estimation cases (no underscore in the name), or the
    follow-ups (kind="follow")."""
    names = sorted(f[:-5] for f in os.listdir(FIXTURES) if f.endswith(".json"))
    if kind == "follow":
        return [n for n in names if "_" in n]
    if kind == "long":
        return [n for n in names if n.startswith("s")]
    return [n for n in names if "_" not in n and not n.startswith("s")]


# ---- fixtures/more: more cases, on the test data and on the datasets of fixtures/more/data ----
MORE = os.path.join(FIXTURES, "more")
_more = {}


def more_cases():
    """The cases of fixtures/more/cases.txt, in order: (name, dataset, command, follow-ups),
    the follow-ups as (name, command)."""
    if "cases" not in _more:
        cases = []
        with open(os.path.join(MORE, "cases.txt")) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("*"):
                    continue
                head, cmd = line.lstrip("+").split(":", 1)
                if line.startswith("+"):
                    cases[-1][3].append((head.strip(), cmd.strip()))
                else:
                    name, ds = [x.strip() for x in head.split("@")]
                    cases.append((name, ds, cmd.strip(), []))
        _more["cases"] = cases
    return _more["cases"]


def more_stata(name):
    """What Stata gave for a case or a follow-up of fixtures/more (None when it was not run: a
    follow-up of a case that Stata refused)."""
    if "stata" not in _more:
        import gzip
        with gzip.open(os.path.join(MORE, "stata.json.gz"), "rt", encoding="utf-8") as f:
            _more["stata"] = json.load(f)
        _more["weights"] = np.load(os.path.join(MORE, "weights.npz"))
    d = _more["stata"].get(name)
    if d is not None and "scalars" in d and not d.get("_read"):
        d["scalars"] = {k: (float("nan") if v is None else float(v)) for k, v in d["scalars"].items()}
        for m in d["matrices"].values():
            m["data"] = [[float("nan") if x is None else float(x) for x in row] for row in m["data"]]
        d["_read"] = True
    return d


def more_data(ds):
    """A dataset of fixtures/more with its variable and value labels: t is the test data, the
    others are the .dta files of fixtures/more/data."""
    if ds == "t":
        return load_data(), LABELS, VALUE_LABELS
    if ds not in _more:
        df = pd.read_stata(os.path.join(MORE, "data", ds + ".dta"), convert_categoricals=False)
        with open(os.path.join(MORE, "data", ds + ".labels.json")) as f:
            lab = json.load(f)
        vl = {v: {(int(float(k)) if float(k) == int(float(k)) else float(k)): t for k, t in d.items()}
              for v, d in lab["value_labels"].items()}
        _more[ds] = (df, lab["labels"], vl)
    df, labels, vl = _more[ds]
    return df.copy(), labels, vl


def run_more(name, ds, cmd):
    """A case of fixtures/more in Python, with Stata's replicate weights when it has them:
    (Stata's results, the data, the results or None, return code, message, warnings)."""
    import contextlib
    import io
    import warnings
    import seldecomp as sd
    case = more_stata(name)
    df, labels, vl = more_data(ds)
    kw = {}
    if name + ":bs" in _more["weights"].files:
        kw["_bsweights"] = _more["weights"][name + ":bs"]
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            res = sd.seldecomp(df, cmd, labels=labels, value_labels=vl, **kw)
            rc, msg = 0, ""
        except sd.SeldecompError as ex:
            res, rc, msg = None, ex.rc, ex.message
    warned = [str(w.message) for w in caught if issubclass(w.category, sd.SeldecompWarning)]
    return case, df, res, rc, msg, warned


def more_generated(name):
    """The weights that generate() made in Stata for a case: {column: values, one per row}."""
    more_stata(name)
    pre = name + ":gen:"
    return {k[len(pre):]: _more["weights"][k] for k in _more["weights"].files if k.startswith(pre)}


def squeeze(lines):
    """Lines of output as one text, blanks squeezed (a wrapped paragraph matches)."""
    return " ".join(" ".join(l.split()) for l in lines if l.strip())


def trimmed(lines):
    """Lines without trailing blanks and without blank lines at either end."""
    out = [l.rstrip() for l in lines]
    while out and out[0] == "":
        out.pop(0)
    while out and out[-1] == "":
        out.pop()
    return out

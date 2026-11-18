"""The estimation sample: weights, groups, and the coded raking variables."""
import numpy as np
import pandas as pd

from ._parse import is_string_series, unquote
from ._stata import pctile, plural
from .errors import SeldecompError


def fmt_val(v, single=False):
    """A value as Stata shows it in text: an integer without decimals, otherwise up to 16
    significant digits without a leading zero. single: the value comes from a 4-byte float
    column; it is shown in the fewest digits (7 to 9) that identify it (.1, not
    .1000000014901161)."""
    if isinstance(v, str):
        return v
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if f == int(f) and abs(f) < 1e15:
        return str(int(f))
    s = "%.16g" % f
    if single:
        for d in (7, 8, 9):
            t = "%.*g" % (d, f)
            if float(np.float32(float(t))) == f:
                s = t
                break
    if s.startswith("0."):
        s = s[1:]
    elif s.startswith("-0."):
        s = "-" + s[2:]
    return s


def fmt_count(n):
    return "{:,}".format(int(n))


def group_values_text(vals, labels):
    """The values of the by() variable as the messages list them: "0 (Single) and 1
    (Married)"; "1, 2, and 3"; at most 6, then "..."."""
    items = []
    for v in vals[:6]:
        t = '"%s"' % v if isinstance(v, str) else fmt_val(v)
        if labels and v in labels:
            t = "%s (%s)" % (t, labels[v])
        items.append(t)
    if len(vals) > 6:
        return ", ".join(items) + ", ..."
    if len(items) == 2:
        return items[0] + " and " + items[1]
    if len(items) > 2:
        return ", ".join(items[:-1]) + ", and " + items[-1]
    return items[0] if items else ""


def bin_label(b, cuts):
    """The text of quantile group b (1-based) given the formatted cut points."""
    q = len(cuts) + 1
    if b == 1:
        return "<= %s" % cuts[0]
    if b == q:
        return "> %s" % cuts[-1]
    return "(%s, %s]" % (cuts[b - 2], cuts[b - 1])


def bin_codes(x, cuts):
    """Quantile group of each value: 1 for x <= c1, ..., q for x above the last cut point."""
    return np.searchsorted(np.asarray(cuts, dtype=float), x, side="left") + 1


def cut_points(x, w, q, what):
    if len(x) < q:
        raise SeldecompError("contvars(): too few observations to form %d quantile groups of %s"
                             % (q, what), 2001)
    return [float(c) for c in pctile(x, w, nq=q)]


class Coded:
    """One term coded on a sample: codes 1..nlev, the labels of the levels, and for each
    level the values of the term's components."""

    def __init__(self, codes, nlev, labels, values):
        self.codes = codes
        self.nlev = nlev
        self.labels = labels
        self.values = values


def code_term(term, df, touse, cutinfo, value_labels):
    """Number the levels of a term 1, 2, ... in the ascending order of its components on the
    sample touse. cutinfo: {variable: (cut points, formatted cut points)} for c. variables."""
    comps = []
    kinds = []
    for v, kind in zip(term.vars, term.kinds):
        col = df[v]
        if isinstance(col.dtype, pd.CategoricalDtype):
            # a pandas categorical: its categories in their order, shown by name
            comps.append(col.cat.codes.to_numpy()[touse].astype(float))
            kinds.append(("g", list(col.cat.categories)))
        elif is_string_series(col):
            vals = col.to_numpy()[touse]
            kinds.append("s")
            comps.append(vals.astype(str))
        elif kind == "c":
            cuts = cutinfo[v][0]
            comps.append(bin_codes(col.to_numpy(dtype=float)[touse], cuts))
            kinds.append("c")
        else:
            comps.append(col.to_numpy(dtype=float)[touse])
            kinds.append("f" if col.dtype == np.float32 else "i")
    frame = pd.DataFrame({i: c for i, c in enumerate(comps)})
    uniq = frame.drop_duplicates().sort_values(list(frame.columns), kind="mergesort").reset_index(drop=True)
    uniq["_code"] = np.arange(1, len(uniq) + 1)
    merged = frame.merge(uniq, how="left", on=list(frame.columns), sort=False)
    codes = merged["_code"].to_numpy(dtype=np.int64)
    labels = []
    values = []
    nv = len(term.vars)
    for row in uniq.itertuples(index=False):
        pieces = []
        vals = []
        for j, (v, kind) in enumerate(zip(term.vars, kinds)):
            val = row[j]
            vals.append(val)
            if isinstance(kind, tuple):
                piece = str(kind[1][int(val)])
            elif kind == "s":
                piece = str(val)
            elif kind == "c":
                piece = bin_label(int(val), cutinfo[v][1])
            else:
                piece = fmt_val(val, single=(kind == "f"))
                vl = (value_labels or {}).get(v)
                if vl is not None and float(val) == int(val):
                    lab = vl.get(int(val), vl.get(val))
                    if lab not in (None, ""):
                        piece = "%s (%s)" % (piece, lab)
            if nv > 1:
                piece = "%s=%s" % (v, piece)
            pieces.append(piece)
        labels.append(" ".join(pieces).replace("|", "/"))
        values.append(vals)
    return Coded(codes, len(uniq), labels, values)


def recode_present(coded, keep):
    """The codes renumbered 1, 2, ... over the levels still present where keep is true (in
    their order); returns the new Coded for the kept observations."""
    c = coded.codes[keep]
    present = np.unique(c)
    newcode = np.zeros(coded.nlev + 1, dtype=np.int64)
    newcode[present] = np.arange(1, len(present) + 1)
    return Coded(newcode[c], len(present), [coded.labels[p - 1] for p in present],
                 [coded.values[p - 1] for p in present])


def stage_text(nstages, stages):
    """"stage 2", "stages 1-3, 5": the stages (1-based) a note refers to; empty with one stage."""
    stages = list(stages)
    if nstages <= 1 or not stages:
        return ""
    parts = []
    a = b = stages[0]
    for s in stages[1:]:
        if s == b + 1:
            b = s
        else:
            parts.append("%d-%d" % (a, b) if b > a else "%d" % a)
            a = b = s
    parts.append("%d-%d" % (a, b) if b > a else "%d" % a)
    return ("stages " if len(stages) > 1 else "stage ") + ", ".join(parts)


def exclusion_note(n, n1, n2, what, by, tail=""):
    other = n - n1 - n2
    xo = ", %s with other values of %s" % (fmt_count(other), by) if other > 0 else ""
    return ("%s %s with %s excluded (%s in the focal group, %s in the reference group%s)%s"
            % (fmt_count(n), plural(n, "observation"), what, fmt_count(n1), fmt_count(n2), xo, tail))

"""makeplot(): the graphs of a staged decomposition (sequentially or individually), drawn with
matplotlib from the stored results.

    stages           the cumulative explained and unexplained parts at each raking stage, with
                     the total difference as a dashed line
    marginal         the total difference, the explained part, the marginal explained part of
                     each stage and the unexplained part (the only type after individually)
    waterfall        each stage's marginal explained part as a step from 0, then the explained
                     part, the unexplained part up to the total difference, and the total
                     difference (the default for sequentially)
    waterfalllevels  the same steps in levels of the mean: from the focal group mean through the
                     stages to the raked focal group mean, then to the reference group mean
"""
import math
import re
import textwrap
from decimal import Decimal, ROUND_HALF_UP

from ._parse import split_top, unquote
from ._stata import fmt_g, invnormal
from .errors import SeldecompError

TYPES = ("stages", "marginal", "waterfall", "waterfalllevels")
CI_FLAGS = ("normal", "tinterval", "percentile", "bc", "bca")
CI_NAME = {"normal": "normal", "tinterval": "t", "t": "t", "percentile": "percentile", "bc": "bc",
           "bca": "bca"}
CI_LABEL = {"normal": "normal", "t": "t", "percentile": "percentile", "bc": "bias-corrected", "bca": "BCa"}
FMT_RE = re.compile(r"^%(-?)(0?)([0-9]+)([.,])([0-9]+)([efg])(c?)$")
MARGINAL_TITLES = ("Marginal explanatory contribution at each raking stage",
                   "Marginal explanatory contribution by raking stage",
                   "Marginal explanatory contribution at each stage",
                   "Marginal explanatory contribution by stage")
# the kinds of rows: a mean or the total difference, the explained part, a block, the unexplained
# part, the remainder of individually
MEAN, EXPL, BLOCK, UNEXPL, REMAINDER = 1, 2, 3, 4, 5
GREY_FILL, GREY_EDGE, GREY_TEXT, LINE = "0.75", "0.5", "0.375", "0.5"


def _err(message, rc=198):
    return SeldecompError(message, rc)


# ---- numbers and texts ------------------------------------------------------------------------

def _round_half_up(x, dec):
    q = Decimal(x).quantize(Decimal(1).scaleb(-dec), rounding=ROUND_HALF_UP)
    return "{:.{}f}".format(q, dec)


def _efmt(x, digits, room):
    s = "%.*e" % (digits, x)
    while digits > 0 and len(s.lstrip("-")) > room:
        digits -= 1
        s = "%.*e" % (digits, x)
    return s


def _thousands(s):
    sign = "-" if s.startswith("-") else ""
    ip, dot, fp = s.lstrip("-").partition(".")
    return sign + ("{:,}".format(int(ip)) if ip else "") + dot + fp


def sfmt(x, fmt):
    """x as Stata's string(x, fmt) shows it, fmt a numeric format %w.df, %w.de or %w.dg (with c
    for thousands separators, a comma for the decimal comma, a leading 0 for zero padding)."""
    left, zero, w, mark, d, kind, comma = FMT_RE.match(fmt).groups()
    w, d, comma = int(w), int(d), comma == "c"
    if x is None or not math.isfinite(x):
        return "."
    x = float(x)
    if kind == "e":
        s = "%.*e" % (d, x)
    elif kind == "f":
        s = ("{:,.%df}" % d).format(x) if comma else "%.*f" % (d, x)
        if len(s) > w:
            s = _efmt(x, d if d > 0 else 2, w - 1)
    elif d == 0:
        s = fmt_g(x, w)
        if comma and "e" not in s:
            s = _thousands(s)
    elif x == 0:
        s = "0"
    else:
        # d significant digits (at least 2); the integer part whole when it has more digits;
        # e-notation for small numbers and when the number does not fit the width
        sig = max(d, 2)
        k = math.floor(math.log10(abs(x)))
        s = None
        if k >= -5:
            if k + 1 > sig:
                s = _round_half_up(x, 0)
            else:
                t = _round_half_up(x, sig - 1 - k)
                if len(t.lstrip("-").split(".")[0]) <= sig:
                    s = t.rstrip("0").rstrip(".") if "." in t else t
            if s is not None:
                if comma:
                    s = _thousands(s)
                if s.startswith("0."):
                    s = s[1:]
                elif s.startswith("-0."):
                    s = "-" + s[2:]
                if len(s.lstrip("-")) > w - 1:
                    s = None
        if s is None:
            s = _efmt(x, sig - 1, w - 1)
    if zero and not left and len(s) < w and "e" not in s:
        s = ("-" + s[1:].rjust(w - 1, "0")) if s.startswith("-") else s.rjust(w, "0")
    if mark == ",":
        s = s.replace(",", "\0").replace(".", ",").replace("\0", ".")
    return s


def _valid_format(fmt):
    m = FMT_RE.match(fmt) if isinstance(fmt, str) else None
    return m is not None and int(m.group(3)) >= 1 and int(m.group(5)) < int(m.group(3))


def _num(v, fmt):
    """A value labelled on a bar: in the format, a value that reads 0 without its sign."""
    s = sfmt(v, fmt).strip()
    if re.match(r"^-[0.,]+$", s):
        s = s[1:]
    return s


def _pct(p):
    """A share in percent, rounded as Stata's round() does (half up), never -0."""
    return "0" if abs(p) < 0.5 else "%d" % math.floor(p + 0.5)


def _pfit(s, n):
    """s shortened to at most n characters, ending with "..." when it is cut: after the last
    whole word when that keeps more than half of the room."""
    if len(s) <= n:
        return s
    t = s[:n - 3]
    if s[n - 3] != " ":
        p = t.rfind(" ") + 1
        if p > (n - 3) / 2:
            t = t[:p - 1]
    return re.sub(r"[ ,;:]+$", "", t) + "..."


def _wrap(s, n, tail=""):
    """s on at most two lines of n characters, then tail (Stata's seldecomp_pwrap): one line when
    s has at most n characters; else split at a space with both lines at most n, after a comma
    when that fits, else where the longer line is shortest; when no split fits, the words that
    fit in n on the first line and the rest on the second, shortened by _pfit (a first word
    longer than n shortened likewise; a text without spaces is one shortened line)."""
    if len(s) <= n:
        return s + tail
    s = s.strip()
    m = len(s)
    if m <= n:
        return s + tail
    if " " not in s:
        return _pfit(s, n) + tail
    best, bk, bc = 0, None, 0
    for i in range(1, m - 1):
        if s[i] != " ":
            continue
        a, b = s[:i].rstrip(), s[i + 1:].lstrip()
        k = max(len(a), len(b))
        c = 1 if a.endswith(",") else 0
        if k <= n and (c > bc or (c == bc and (bk is None or k < bk))):
            best, bk, bc = i, k, c
    if not best:
        for i in range(min(n, m - 1), 0, -1):
            if s[i] == " ":
                best = i
                break
    if best:
        a, b = s[:best].rstrip(), s[best + 1:].lstrip()
    else:
        i = s.index(" ")
        a, b = _pfit(s[:i], n), s[i + 1:].lstrip()
    return a + "\n" + _pfit(b, n) + tail


def _axis_title(a, b):
    """The title of the value axis, "a: b", on two lines beyond 40 characters."""
    if len(a + ": " + b) <= 40:
        return a + ": " + b
    return _pfit(a, 40) + ":\n" + b


def _tx(s):
    """A text of the data or of blocklabels as matplotlib shows it as typed: every $ escaped,
    since two of them would start mathtext."""
    return s.replace("$", r"\$")


def _gfmt(v):
    return fmt_g(v, 12)


def _tick(lo, hi):
    """A round step (1, 2, 2.5 or 5 times a power of 10) that gives about five ticks over lo to
    hi."""
    raw = (hi - lo) / 5
    mag = 10.0 ** math.floor(math.log10(raw))
    f = raw / mag
    return mag * (1 if f <= 1 else 2 if f <= 2 else 2.5 if f <= 2.5 else 5 if f <= 5 else 10)


# ---- the plotspec -------------------------------------------------------------------------------

def _interval_type(ci, citype, flags):
    """The interval type asked for (None: the default of the results), with the refusals of
    the plotspec."""
    asked = []
    if isinstance(ci, str):
        asked.append(ci)
        ci = True
    if citype not in (None, False, ""):
        asked.append(str(citype))
    asked += [f for f in CI_FLAGS if flags.get(f)]
    names = []
    for a in asked:
        if a.lower() not in CI_NAME:
            raise _err("makeplot(): %s is not an interval type; the types are normal, tinterval, "
                       "percentile, bc, and bca" % a)
        if CI_NAME[a.lower()] not in names:
            names.append(CI_NAME[a.lower()])
    if len(names) > 1:
        raise _err("makeplot(): specify only one interval type: normal, tinterval, percentile, bc, or bca")
    if names and not ci:
        raise _err("makeplot(): specify ci with an interval type")
    return bool(ci), (names[0] if names else None)


def _check_spec(type, chainstep, outcome, format):
    if type is not None:
        hits = [p for p in TYPES if p == str(type).strip().lower()]
        if not hits:
            raise _err("makeplot(): plot type %s not allowed; the types are stages, marginal, "
                       "waterfall, and waterfalllevels" % type)
        type = hits[0]
    if chainstep is not None:
        try:
            ok = not isinstance(chainstep, bool) and float(chainstep) == int(float(chainstep))
        except (TypeError, ValueError):
            ok = False
        if not ok or int(float(chainstep)) < 1:
            raise _err("makeplot(chainstep()) must be a positive integer")
        chainstep = int(float(chainstep))
    if outcome is not None:
        words = str(outcome).split()
        if len(words) > 1:
            raise _err("makeplot(outcome()): %s names %d outcomes; name one" % (outcome, len(words)))
        outcome = words[0] if words else None
    if format is not None and not _valid_format(format):
        raise _err("makeplot(format()): %s is not a numeric display format" % format)
    return type, chainstep, outcome, format


def _block_texts(blocks, given):
    """blocklabels: the text of each block (its name when none is given), with the refusals of
    Stata's Plot_blabels."""
    pairs = []
    if isinstance(given, dict):
        pairs = [(str(k), "" if v is None else str(v)) for k, v in given.items()]
    elif isinstance(given, (list, tuple)) and all(isinstance(q, (list, tuple)) and len(q) == 2 for q in given):
        pairs = [(str(k), "" if v is None else str(v)) for k, v in given]
    elif given not in (None, False, ""):
        toks = split_top(str(given))
        if len(toks) % 2:
            raise _err('makeplot(blocklabels()): type pairs of a block name and its text, as in '
                       'blocklabels(demo "Demographics")')
        pairs = [(toks[i], unquote(toks[i + 1])) for i in range(0, len(toks), 2)]
    texts = dict((bn, bn) for bn in blocks)
    seen = set()
    for name, text in pairs:
        if name not in texts:
            raise _err("makeplot(blocklabels()): %s is not a block; the blocks are %s"
                       % (name, " ".join(blocks)))
        if name in seen:
            raise _err("makeplot(blocklabels()): block %s is given twice" % name)
        seen.add(name)
        if text == "":
            raise _err('makeplot(blocklabels()): the text of %s is empty; for a blank label, type " "'
                       % name)
        texts[name] = text
    return [texts[bn] for bn in blocks], seen


# ---- what is drawn ------------------------------------------------------------------------------

class _Plot:
    """One graph: what is read from the results, and the texts."""


def _setup(res, type, swapaxis, ci, citype, reverse, chainstep, outcome, fmt, blocklabels, varlabels,
           level=None):
    e = res.e
    ctx = getattr(res, "_ctx", None)
    labels = dict(getattr(ctx, "labels", None) or {})
    varlabels = bool(varlabels)
    P = _Plot()
    P.e = e
    P.b = res.b
    P.ns = int(e.get("N_stages", 0))
    eblocks = str(e.get("blocks", "")).split()
    P.rev = bool(reverse)
    P.rp = ""
    blocks = list(eblocks)
    if P.rev:
        if not e.get("reverse"):
            raise _err("makeplot(reverse) requires results estimated with sequentially reverse")
        P.rp = "rev_"
        blocks = blocks[::-1]
    texts, given = _block_texts(blocks, blocklabels)

    def vlab(s):
        if not varlabels:
            return s
        return ", ".join(" # ".join(labels.get(v) or v for v in t.split("#")) for t in s.split())

    # varlabels: a block named after its terms and without a text of blocklabels() is shown by
    # their variable labels
    for k, bn in enumerate(blocks):
        if not varlabels or bn in given or bn not in eblocks:
            continue
        terms = str(e.get("block%d" % (eblocks.index(bn) + 1), ""))
        if bn in (terms, terms.replace("#", "_")):
            texts[k] = vlab(terms)
    P.variant = e.get("stagetype") == "individually"
    if P.variant:
        if type not in (None, "marginal"):
            raise _err("makeplot(%s) not allowed with individually; use marginal" % type)
        type = "marginal"
        blocks.append("_remainder")
        texts.append("Remainder")
    P.blocks, P.texts = blocks, texts
    P.fmt = fmt or "%9.3g"
    P.user_fmt = fmt is not None
    P.ci = None
    if ci and res.V is None:
        raise _err("makeplot(ci): these results have no confidence intervals")
    # outcomes(): the first outcome unless outcome= names another; chain: the pair of
    # chainstep=, else the whole chain
    P.dv = str(e.get("depvar", ""))
    P.eq0, P.pfx = "overall", ""
    P.D = float(e.get("difference", float("nan")))
    oi = 1
    outs = str(e.get("outcomes", "")).split()
    if outs:
        if outcome is None:
            outcome = outs[0]
            print("(makeplot: outcome %s, the first of the outcomes; outcome=name draws another)" % outcome)
        elif outcome not in outs:
            raise _err("makeplot(outcome()): %s is not one of the outcomes (%s)" % (outcome, " ".join(outs)))
        P.eq0, P.pfx, P.dv = outcome, outcome + "_", outcome
        P.D = _bval(P.b, outcome + ":difference")
        oi = outs.index(outcome) + 1
    elif outcome is not None and outcome != P.dv:
        raise _err("makeplot(outcome()): %s is not the outcome (%s)" % (outcome, P.dv))
    P.chain = bool(e.get("chain"))
    P.pi = 0
    if chainstep is not None and not P.chain:
        raise _err("makeplot(chainstep()) requires results from a chain run")
    if P.chain:
        npairs = int(e.get("N_pairs", 0))
        if chainstep is not None:
            if chainstep > npairs:
                raise _err("makeplot(chainstep()) must be 1 to %d" % npairs)
            P.pi = chainstep
            P.eq0, P.pfx = "p%d" % P.pi, "p%d_" % P.pi
        else:
            if type in ("stages", "waterfalllevels"):
                raise _err("makeplot(%s) requires chainstep(#) with chain results" % type)
            P.eq0, P.pfx = "total", "total_"
        P.D = _bval(P.b, P.eq0 + ":difference")
        oi = P.pi
    P.type = type or "waterfall"
    P.horiz = (P.type == "marginal") != bool(swapaxis)
    if ci:
        P.ci, P.cilabel, P.cilevel = _intervals(res, citype, level)
    # the rounding threshold: 1e-13 times the largest of the outcome's (the pair's) means
    if P.chain:
        eqs = ["p%d" % oi] if oi > 0 else ["p%d" % (j + 1) for j in range(int(e.get("N_pairs", 0)))]
    else:
        eqs = [outs[oi - 1]] if outs else ["overall"]
    s = 0.0
    for eq in eqs:
        for r in ("focal_group", "raked_focal", "ref_group"):
            if eq + ":" + r in P.b.index:
                s = max(s, abs(float(P.b[eq + ":" + r])))
    P.eps = 1e-13 * s
    P.Dz = abs(P.D) < P.eps
    P.Dtxt = _num(0.0 if P.Dz else P.D, P.fmt)
    # the groups: their value labels, else "by = value"
    gr, sr = [], []
    for g in (1, 2):
        gv = str(e.get("focal_group" if g == 1 else "ref_group", ""))
        gl = str(e.get("focal_group_label" if g == 1 else "ref_group_label", ""))
        if P.chain:
            kk = max(P.pi, 1) if g == 1 else (P.pi + 1 if P.pi > 0 else int(e["N_pairs"]) + 1)
            gv = str(e.get("by_levels", "")).split()[kk - 1]
            gl = str(e.get("by_labels", "")).split("|")[kk - 1].strip()
        numeric = _is_number(gv)
        if gl != gv:
            name = gl[len(gv) + 2:-1]
        elif numeric:
            name = vlab(str(e.get("by", ""))) + " = " + gv
        else:
            name = gv
        gr.append(name)
        sr.append("(" + name + ")" if gl == gv and numeric else name)
    P.gr, P.sr = gr, sr
    ssfx = ""
    if P.chain and P.pi == 0:
        ssfx = " (sum of %d pairs)" % int(e["N_pairs"])
    if P.rev:
        ssfx += ", reverse order"
    if len(sr[0] + " - " + sr[1] + ssfx) <= 72:
        P.sub = sr[0] + " - " + sr[1] + ssfx
    else:
        P.sub = _pfit(sr[0], 70) + " -\n" + _pfit(sr[1], 70) + ssfx
    # the value axis: the outcome's variable label (its name without one); ratio(): that of the
    # denominator too
    rv = str(e.get("ratio", "") or "")
    P.ymeas = "ratio" if rv else "mean"
    P.yname = (labels.get(P.dv) or P.dv) + ((" / " + (labels.get(rv) or rv)) if rv else "")
    return P


def _is_number(s):
    try:
        return math.isfinite(float(s))
    except ValueError:
        return False


def _bval(b, name):
    if name not in b.index:
        raise _err("makeplot: %s not found in e(b)" % name, 111)
    return float(b[name])


def _intervals(res, citype, level=None):
    """The interval of each name of e(b) as {name: (ll, ul)}, the label of its type and the
    level: the stored intervals of the type asked for, else the first type of the results; at
    another level than the results', the intervals recomputed from the bootstrap replicates."""
    e = res.e
    if res.V is None:
        raise _err("makeplot(ci): these results have no confidence intervals")
    stored = res.matrices
    lev = float(e.get("level", 95))
    if level is not None:
        from ._display import _Display, _level
        if _level(level) != lev:
            lev = _level(level)
            if "bsdraws" not in stored:
                raise _err("makeplot(ci): these results store intervals only at level(%s)" % _gfmt(float(e["level"])))
            import pandas as pd
            ci = _Display(res, False).display_ci(lev)[0]
            names = str(e.get("bsnames", "")).split()
            stored = {k: pd.DataFrame(v, index=["ll", "ul"], columns=names) for k, v in ci.items()
                      if k.startswith("ci_") and k != "pct_explained_ci" and getattr(v, "shape", (0, 0))[1:] == (len(names),)}
    t = citype
    if t is None:
        words = str(e.get("citypes", "") or "").split()
        if words:
            t = CI_NAME.get(words[0], words[0])
        else:
            dfr, ncl = e.get("df_r", float("nan")), e.get("N_clust", float("nan"))
            t = "t" if (0 < dfr < float("inf")) or (0 < ncl < float("inf")) else "normal"
    C = stored.get("ci_" + t)
    out = {}
    if C is not None:
        for name in C.columns:
            out[str(name)] = (float(C.loc["ll", name]), float(C.loc["ul", name]))
    elif t in ("normal", "t"):
        from scipy import stats
        level = lev
        if t == "t":
            df = e.get("df_bs", e.get("df_r"))
            if not (df is not None and 0 < float(df) < float("inf")):
                raise _err("makeplot(tinterval): these results have no usable degrees of freedom")
            crit = float(stats.t.ppf(1 - (100 - level) / 200, float(df)))
        else:
            crit = invnormal(1 - (100 - level) / 200)
        V = res.V
        for name, v in res.b.items():
            var = float(V.loc[name, name]) if name in V.index else float("nan")
            if var >= 0:
                out[str(name)] = (float(v) - crit * math.sqrt(var), float(v) + crit * math.sqrt(var))
    else:
        raise _err("makeplot(%s): these results have no %s intervals; request them when estimating "
                   "inference" % (t, t))
    return out, CI_LABEL[t], lev


def _ci_of(P, name):
    if P.ci is None:
        return None
    lu = P.ci.get(name)
    if lu is None or not (math.isfinite(lu[0]) and math.isfinite(lu[1])):
        return None
    return lu


class _Row:
    def __init__(self, **kw):
        self.__dict__.update(kw)


def _label_value(P, v, kind, plus):
    """The label of a bar: its value (a value that reads 0 without a sign, + before a positive
    step) and its share of the difference."""
    vz = 0.0 if abs(v) < P.eps else v
    lab = _num(vz, P.fmt)
    if plus and vz > 0 and not re.match(r"^[0.,]+$", lab):
        lab = "+" + lab
    return vz, lab


def _rows_waterfall(P):
    nb = len(P.blocks)
    names = [P.pfx + P.rp + "marg_expl:" + bn for bn in P.blocks]
    names += [P.eq0 + ":explained", P.eq0 + ":unexplained", P.eq0 + ":difference"]
    vals = [_bval(P.b, n) for n in names]
    rows = []
    run = 0.0
    for r, (name, val) in enumerate(zip(names, vals)):
        if r < nb:
            kind = REMAINDER if (P.variant and r == nb - 1) else BLOCK
            w0, w1 = run, run + val
            run = w1
            tick = P.texts[r]
        elif r == nb:
            kind, w0, w1, tick = EXPL, 0.0, val, "Expl."
        elif r == nb + 1:
            kind, w0, w1, tick = UNEXPL, vals[nb], vals[nb + 2], "Unexpl."
        else:
            kind, w0, w1, tick = MEAN, 0.0, val, "Total diff"
        vz, lab = _label_value(P, val, kind, kind in (BLOCK, UNEXPL, REMAINDER))
        if not P.Dz and kind != MEAN:
            lab += " (%s%%)" % _pct(100 * vz / P.D)
        lu = _ci_of(P, name)
        ll = ul = None
        if lu is not None:
            ll, ul = w1 - val + lu[0], w1 - val + lu[1]
        up = val >= 0
        anc = w1 if lu is None else (max(w1, ul) if up else min(w1, ll))
        rows.append(_Row(name=name, kind=kind, pos=float(r + 1), w0=w0, w1=w1, val=val, ll=ll, ul=ul,
                         lab=lab, anc=anc, up=up, tick=tick, conn=w1 if r < nb + 2 else None))
    return rows


def _rows_marginal(P):
    nb = len(P.blocks)
    names = [P.eq0 + ":difference", P.eq0 + ":explained"]
    names += [P.pfx + P.rp + "marg_expl:" + bn for bn in P.blocks] + [P.eq0 + ":unexplained"]
    nr = nb + 3
    rows = []
    for r, name in enumerate(names, 1):
        val = _bval(P.b, name)
        if r == 1:
            kind, pos, tick = MEAN, 1.0, "Total diff"
        elif r == 2:
            kind, pos, tick = EXPL, r + 0.2, "Expl."
        elif r < nr:
            kind = REMAINDER if (P.variant and r == nr - 1) else BLOCK
            pos, tick = r + 0.2, P.texts[r - 3]
        else:
            kind, pos, tick = UNEXPL, r + 0.4, "Unexpl."
        vz, lab = _label_value(P, val, kind, False)
        if not P.Dz and kind != MEAN:
            lab += " (%s%%)" % _pct(100 * vz / P.D)
        lu = _ci_of(P, name)
        ll, ul = (None, None) if lu is None else lu
        up = val >= 0
        anc = val if lu is None else (max(val, ul) if up else min(val, ll))
        rows.append(_Row(name=name, kind=kind, pos=pos, w0=0.0, w1=val, val=val, ll=ll, ul=ul,
                         lab=lab, anc=anc, up=up, tick=tick, conn=None))
    return rows


def _mean_labels(vals, vzs, eps):
    """The labels of the three means of the waterfall in levels: 4 significant digits, and at
    least the decimals that show the smallest gap between them, sharing the most decimals any
    of them keeps (up to their own); thousands separated from 10,000."""
    mx = max(abs(v) for v in vals)
    gaps = [abs(vals[1] - vals[0]), abs(vals[2] - vals[1]), abs(vals[2] - vals[0])]
    gaps = [g for g in gaps if g >= eps]
    dg = 0
    if gaps:
        dg = max(0, -math.floor(math.log10(min(gaps))))
        if mx > 0 and math.floor(math.log10(mx)) + 1 + dg > 10:
            dg = 0
    sep = mx >= 9999.5

    def fixed(v, d):
        s = ("{:,.%df}" % d).format(v) if sep else "%.*f" % (d, v)
        if s.startswith("0."):
            s = s[1:]
        elif s.startswith("-0."):
            s = "-" + s[2:]
        return s

    labs, dm = [], []
    for vz in vzs:
        d = max(0, dg) if vz == 0 else max(0, 3 - math.floor(math.log10(abs(vz))), dg)
        dm.append(d)
        if vz == 0:
            s = "0"
        elif mx >= 1e15 or d > 20:
            s = sfmt(vz, "%9.4g").strip()
        else:
            s = ("{:,.%df}" % d).format(vz) if sep else "%.*f" % (d, vz)
            k = d - dg
            while k > 0 and s.endswith("0"):
                s = s[:-1]
                k -= 1
            if s.endswith("."):
                s = s[:-1]
            if s.startswith("0."):
                s = s[1:]
            elif s.startswith("-0."):
                s = "-" + s[2:]
        labs.append(s)
    if mx < 1e15:
        dcol = max([len(s) - s.index(".") - 1 for s in labs if "." in s and "e" not in s] + [0])
        for i, vz in enumerate(vzs):
            s = labs[i]
            k = len(s) - s.index(".") - 1 if "." in s else 0
            t = min(dcol, dm[i])
            if vz != 0 and k < t and "e" not in s:
                labs[i] = fixed(vz, t)
    return labs


def _rows_levels(P):
    nb = len(P.blocks)
    names = [P.eq0 + ":focal_group"] + [P.pfx + P.rp + "marg_expl:" + bn for bn in P.blocks]
    names += [P.eq0 + ":raked_focal", P.eq0 + ":unexplained", P.eq0 + ":ref_group"]
    nr = nb + 4
    vals = [_bval(P.b, n) for n in names]
    # the row names of the groups: on two lines from 20 characters in the horizontal view; the
    # vertical view angles them beyond seven bars or with long names, and shortens them
    longest = max([len(t) for t in P.texts] + [0])
    glen = max(len(P.gr[0]), len(P.gr[1]))
    P.angled = not P.horiz and (nr > 7 or longest > 10 or glen > 7)
    # (beyond nine rows two lines would run into the neighbouring rows: one line, shortened)
    if P.horiz and nr <= 9:
        rw1, rw2, rwr = _wrap(P.gr[0], 20), _wrap(P.gr[1], 20), _wrap(P.gr[0], 20, ", raked")
    else:
        wl = 20 if P.horiz else 18 if nr <= 5 else 15 if nr <= 8 else 13 if nr <= 12 else 12
        rw1, rw2 = _pfit(P.gr[0], wl), _pfit(P.gr[1], wl)
        rwr = rw1 + ", raked" if (P.horiz or P.angled) else P.gr[0] + "\nraked"
    rows = []
    level = vals[0]
    for r, (name, val) in enumerate(zip(names, vals), 1):
        if r == 1:
            kind, w0, w1, tick = MEAN, None, val, rw1
        elif r <= nb + 1:
            kind = REMAINDER if (P.variant and r == nb + 1) else BLOCK
            w0, w1, tick = level, level - val, P.texts[r - 2]
            level = w1
        elif r == nb + 2:
            kind, w0, w1, tick = MEAN, None, val, rwr
        elif r == nb + 3:
            kind, w0, w1, tick = UNEXPL, vals[nb + 1], vals[nr - 1], "Unexpl."
        else:
            kind, w0, w1, tick = MEAN, None, val, rw2
        lu = _ci_of(P, name)
        ll = ul = None
        if lu is not None:
            ll, ul = (lu[0], lu[1]) if kind == MEAN else (w1 + val - lu[1], w1 + val - lu[0])
        rows.append(_Row(name=name, kind=kind, pos=float(r), w0=w0, w1=w1, val=val, ll=ll, ul=ul,
                         tick=tick, conn=w1 if r < nr else None))
    # the labels: the steps as the waterfall's (a step is minus its value in e(b); its share is
    # its value over the difference), the means by _mean_labels or in format=
    means = [rows[0], rows[nb + 1], rows[nr - 1]]
    mvz = [0.0 if abs(m.val) < P.eps else m.val for m in means]
    if P.user_fmt:
        mlabs = [_num(v, P.fmt) for v in mvz]
    else:
        mlabs = _mean_labels([m.val for m in means], mvz, P.eps)
    for m, lab in zip(means, mlabs):
        m.lab = lab
    for row in rows:
        if row.kind == MEAN:
            continue
        stp = -row.val
        vz, lab = _label_value(P, stp, row.kind, True)
        if not P.Dz:
            share = 100 * row.val / P.D
            lab += " (%s%%)" % _pct(0.0 if (vz == 0 or abs(share) < 0.5) else share)
        row.lab = lab
    # the value axis: ticks at a round step; the means' bars start at 0 when 0 is within the
    # levels and intervals, else at the end of the axis nearest 0
    ends = [x for row in rows for x in (row.w0, row.w1, row.ll, row.ul) if x is not None]
    lo, hi = min(ends), max(ends)
    rg = hi - lo
    if rg == 0:
        rg = max(abs(hi), 1.0)
    pad = 0.08 * rg
    P.lo0, P.hi0 = lo, hi
    lo, hi = lo - pad, hi + pad
    if P.lo0 >= 0:
        lo = max(lo, 0.0)
    if P.hi0 <= 0:
        hi = min(hi, 0.0)
    if hi <= lo:
        hi = lo + 1
    tk = _tick(lo, hi)
    P.tk = tk
    P.t0 = math.floor(lo / tk + 1e-9) * tk
    P.t1 = math.ceil(hi / tk - 1e-9) * tk
    P.base_at = "zero" if P.lo0 <= 0 <= P.hi0 else ("low" if P.lo0 > 0 else "high")
    for row in rows:
        if row.kind == MEAN:
            d = row.w1 if P.base_at == "zero" else (1.0 if P.base_at == "low" else -1.0)
        else:
            d = -row.val
        row.up = d >= 0
        row.anc = row.w1 if row.ll is None else (max(row.w1, row.ul) if row.up else min(row.w1, row.ll))
    return rows


def _levels_title(P):
    """From the <focal group> mean to the <reference group> mean, on one line up to 56
    characters, with ", stage by stage" when that fits too; else on two lines."""
    m = P.ymeas
    full = "From the %s %s to the %s %s" % (P.sr[0], m, P.sr[1], m)
    if len(full) + 16 <= 56:
        return full + ", stage by stage"
    if len(full) <= 56:
        return full
    return "From the %s %s\nto the %s %s" % (_pfit(P.sr[0], 42), m, _pfit(P.sr[1], 42), m)


# ---- drawing ------------------------------------------------------------------------------------

def _bar(ax, horiz, pos, w0, w1, width, gid, **style):
    if horiz:
        bars = ax.barh([pos], [w1 - w0], left=[w0], height=width, **style)
    else:
        bars = ax.bar([pos], [w1 - w0], bottom=[w0], width=width, **style)
    bars.patches[0].set_gid(gid)
    return bars.patches[0]


def _spike(ax, horiz, pos, ll, ul, gid):
    """A capped interval from ll to ul at the row position pos."""
    cap = 0.08
    if horiz:
        ax.hlines([pos], [ll], [ul], colors="black", linewidth=0.8, gid=gid)
        ax.vlines([ll, ul], pos - cap, pos + cap, colors="black", linewidth=0.8, gid=gid + ":caps")
    else:
        ax.vlines([pos], [ll], [ul], colors="black", linewidth=0.8, gid=gid)
        ax.hlines([ll, ul], pos - cap, pos + cap, colors="black", linewidth=0.8, gid=gid + ":caps")


def _label(ax, horiz, pos, anc, up, text, gid, size="small", color="black"):
    """A label beyond a bar's end: above or below (right or left in the horizontal view)."""
    if horiz:
        xy, off = (anc, pos), ((3, 0) if up else (-3, 0))
        ha, va = ("left" if up else "right"), "center"
    else:
        xy, off = (pos, anc), ((0, 2) if up else (0, -2))
        ha, va = "center", ("bottom" if up else "top")
    return ax.annotate(text, xy, xytext=off, textcoords="offset points", ha=ha, va=va, fontsize=size,
                       color=color, gid=gid, annotation_clip=False)


def _line(ax, horiz, a, b, at, **style):
    """A line across the rows a to b at the value at."""
    if horiz:
        return ax.plot([at, at], [a, b], **style)[0]
    return ax.plot([a, b], [at, at], **style)[0]


def _style(kind, ptype):
    if kind in (BLOCK, EXPL):
        return dict(color="C0")
    if kind == REMAINDER:
        return dict(color="C0", alpha=0.35)
    if kind == UNEXPL:
        return dict(facecolor=GREY_FILL, edgecolor=GREY_EDGE, linewidth=0.8)
    if ptype == "waterfalllevels":
        return dict(facecolor="none", edgecolor="black", linewidth=1.2)
    return dict(facecolor="none", edgecolor="black", linewidth=1.0, linestyle="--")


def _rows_axis(ax, horiz, pos, ticks, lo, hi, angled=False):
    """The rows: named, without ticks; from the top in the horizontal view."""
    ticks = [_tx(t) for t in ticks]
    if horiz:
        ax.set_yticks(pos)
        ax.set_yticklabels(ticks)
        ax.set_ylim(hi, lo)
        ax.tick_params(axis="y", length=0)
    else:
        ax.set_xticks(pos)
        if angled:
            ax.set_xticklabels(ticks, rotation=30, ha="right", rotation_mode="anchor")
        else:
            ax.set_xticklabels(ticks)
        ax.set_xlim(lo, hi)
        ax.tick_params(axis="x", length=0)


def _set_value_lim(ax, horiz, lo, hi):
    if horiz:
        ax.set_xlim(lo, hi)
    else:
        ax.set_ylim(lo, hi)


def _render(fig):
    try:
        fig.draw_without_rendering()
    except AttributeError:
        fig.canvas.draw()


def _span(ax, horiz, texts):
    """The extent of the texts along the value axis, in data units."""
    inv = ax.transData.inverted()
    lo = hi = None
    for t in texts:
        bb = t.get_window_extent()
        (x0, y0), (x1, y1) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
        a, b = (x0, x1) if horiz else (y0, y1)
        a, b = min(a, b), max(a, b)
        lo = a if lo is None else min(lo, a)
        hi = b if hi is None else max(hi, b)
    return lo, hi


def _fit(fig, ax, horiz, texts, lo, hi, widen):
    """Widen the value axis until the labels lie within it: widen(lo, hi, need_lo, need_hi)
    gives the new ends."""
    for _ in range(8):
        _set_value_lim(ax, horiz, lo, hi)
        if not texts:
            break
        try:
            _render(fig)
            a, b = _span(ax, horiz, texts)
        except Exception:
            break
        nlo, nhi = widen(lo, hi, a, b)
        if (nlo, nhi) == (lo, hi):
            break
        lo, hi = nlo, nhi
    _set_value_lim(ax, horiz, lo, hi)
    return lo, hi


def _widen_margin(lo, hi, a, b):
    m = 0.03 * (hi - lo)
    return (a - m if a < lo else lo), (b + m if b > hi else hi)


def _draw_rows(P, fig, ax, rows):
    """The waterfalls and the marginal plot: a bar per row, the intervals, the labels, the
    connectors from each bar's end to the next bar."""
    horiz = P.horiz
    nr = len(rows)
    small = "x-small" if (not horiz and nr > (7 if P.type == "waterfalllevels" else 6)) else "small"
    texts = []
    for row in rows:
        if not (P.type == "waterfalllevels" and row.kind == MEAN):
            _bar(ax, horiz, row.pos, row.w0, row.w1, 0.6, row.name, **_style(row.kind, P.type))
        if row.ll is not None:
            _spike(ax, horiz, row.pos, row.ll, row.ul, "ci:" + row.name)
        texts.append(_label(ax, horiz, row.pos, row.anc, row.up, row.lab, "label:" + row.name, small))
        if row.conn is not None:
            _line(ax, horiz, row.pos + 0.3, row.pos + 0.7, row.conn, color=LINE, linewidth=0.8,
                  linestyle=(0, (2, 2)), gid="connector:" + row.name)
    top = nr + (0.9 if P.type == "marginal" else 0.5)
    # the row names are angled in the vertical view beyond six bars or with a name longer than
    # 10 characters (the waterfall in levels sets its own rule)
    angled = getattr(P, "angled", False)
    if not horiz and P.type != "waterfalllevels":
        angled = nr > 6 or max([len(t) for t in P.texts] + [0]) > 10
    _rows_axis(ax, horiz, [r.pos for r in rows], [r.tick for r in rows], 0.5, top, angled)
    if P.type == "waterfalllevels":
        _levels_axis(P, fig, ax, rows, texts)
    else:
        ends = [x for row in rows for x in (row.w0, row.w1, row.anc, row.ll, row.ul) if x is not None]
        if P.type == "waterfall":
            ends += [0.0, P.D]
        lo, hi = min(ends + [0.0]), max(ends + [0.0])
        span = hi - lo if hi > lo else 1.0
        lo, hi = lo - 0.03 * span, hi + 0.03 * span
        _fit(fig, ax, horiz, texts, lo, hi, _widen_margin)
        (ax.axvline if horiz else ax.axhline)(0, color=LINE, linewidth=0.8, zorder=0.5, gid="zero")


def _levels_axis(P, fig, ax, rows, texts):
    """The value axis of the waterfall in levels: whole ticks, widened by whole ticks until the
    labels fit (an end at 0 stays); then the means' bars from the base."""
    tk = P.tk

    def widen(t0, t1, a, b):
        n0, n1 = t0, t1
        if a < t0 and t0 != 0:
            n0 = max(t0 - tk, 0.0) if t0 > 0 else t0 - tk
        if b > t1 and t1 != 0:
            n1 = min(t1 + tk, 0.0) if t1 < 0 else t1 + tk
        return n0, n1

    t0, t1 = _fit(fig, ax, P.horiz, texts, P.t0, P.t1, widen)
    if round((t1 - t0) / tk) + 1 > 8:
        tk = _tick(t0, t1)
        t0 = math.floor(t0 / tk + 1e-9) * tk
        t1 = math.ceil(t1 / tk - 1e-9) * tk
    ticks = [t0 + i * tk for i in range(int(round((t1 - t0) / tk)) + 1)]
    if P.horiz:
        ax.set_xticks(ticks)
    else:
        ax.set_yticks(ticks)
    _set_value_lim(ax, P.horiz, t0, t1)
    base = 0.0 if P.base_at == "zero" else (t0 if P.base_at == "low" else t1)
    for row in rows:
        if row.kind == MEAN:
            row.w0 = base
            _bar(ax, P.horiz, row.pos, base, row.w1, 0.6, row.name, **_style(MEAN, P.type))
    (ax.axvline if P.horiz else ax.axhline)(base, color=LINE, linewidth=0.8, zorder=0.5, gid="zero")


def _draw_stages(P, fig, ax):
    """The cumulative explained and unexplained parts at each stage, side by side; the total
    difference as a dashed line; the labels give each part's share of the difference."""
    horiz = P.horiz
    D = P.D
    ns = P.ns
    rows = []
    for s in range(ns):
        bn = P.blocks[s]
        cum = _bval(P.b, P.pfx + P.rp + "cum_expl:" + bn)
        unx = _bval(P.b, P.pfx + P.rp + "cum_unexpl:" + bn)
        for v, nm, off in ((cum, "cum_expl", -0.2), (unx, "cum_unexpl", 0.2)):
            name = P.pfx + P.rp + nm + ":" + bn
            lu = _ci_of(P, name)
            ll, ul = (None, None) if lu is None else lu
            up = v >= 0
            anc = v if lu is None else (max(v, ul) if up else min(v, ll))
            lab = "" if P.Dz else _pct(100 * v / D) + "%"
            rows.append(_Row(name=name, kind=EXPL if nm == "cum_expl" else UNEXPL, pos=s + 1 + off,
                             w0=0.0, w1=v, val=v, ll=ll, ul=ul, up=up, anc=anc, lab=lab))
    ends = [0.0, D] + [r.anc for r in rows]
    lo, hi = min(ends), max(ends)
    # a label that the dashed line of the total would strike through is moved beyond the line:
    # the line lies within the label's reach beyond its anchor, a share of the span of 0, the
    # difference and the anchors with a margin of 8% (10% horizontal) on the sides with bars
    if not P.Dz:
        pad = (hi - lo) * (0.1 if horiz else 0.08)
        rg = hi - lo + pad * (any(r.val < 0 for r in rows) + any(r.val >= 0 for r in rows))
        for r in rows:
            sg = 1.0 if r.val >= 0 else -1.0
            ext = (0.02 + 0.018 * len(r.lab)) if horiz else 0.05
            if 0 < sg * (D - r.anc) < ext * rg:
                r.anc = D
    # the value axis: also the ends of the intervals that reach across 0
    ends += [x for r in rows for x in (r.ll, r.ul) if x is not None]
    lo, hi = min(ends), max(ends)
    span = hi - lo if hi > lo else 1.0
    lo, hi = lo - 0.03 * span, hi + 0.03 * span
    texts = []
    for r in rows:
        _bar(ax, horiz, r.pos, 0.0, r.w1, 0.34, r.name, **_style(r.kind, "stages"))
        if r.ll is not None:
            _spike(ax, horiz, r.pos, r.ll, r.ul, "ci:" + r.name)
        texts.append(_label(ax, horiz, r.pos, r.anc, r.up, r.lab, "label:" + r.name, "small",
                            "black" if r.kind == EXPL else GREY_TEXT))
    _line(ax, horiz, 0.5, ns + 0.5, D, color="black", linewidth=1.0, linestyle="--", gid="total")
    ticks = [P.texts[0]] + ["+ " + t for t in P.texts[1:ns]]
    angled = not horiz and (ns > 4 or max([len(t) for t in P.texts] + [0]) > 10)
    _rows_axis(ax, horiz, [float(s + 1) for s in range(ns)], ticks, 0.5, ns + 0.5, angled)
    _fit(fig, ax, horiz, texts, lo, hi, _widen_margin)
    (ax.axvline if horiz else ax.axhline)(0, color=LINE, linewidth=0.8, zorder=0.5, gid="zero")


def _below(fig, ax, gap):
    """The distance in points from the bottom of the Axes to below its tick labels, axis title
    and legend, plus gap."""
    try:
        _render(fig)
        return (ax.get_window_extent().y0 - ax.get_tightbbox().y0) * 72.0 / fig.dpi + gap
    except Exception:
        return 30.0


def _stages_legend(P, fig, ax):
    """The legend of the stages plot, below the plot: the explained and unexplained parts and
    the total difference with its value."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.transforms import offset_copy
    handles = [Patch(**_style(EXPL, "stages")), Patch(**_style(UNEXPL, "stages")),
               Line2D([], [], color="black", linewidth=1.0, linestyle="--")]
    names = ["Expl.", "Unexpl.", "Total diff (%s)" % P.Dtxt]
    under = offset_copy(ax.transAxes, fig=fig, y=-_below(fig, ax, 4.0), units="points")
    ax.legend(handles, names, loc="upper center", bbox_to_anchor=(0.5, 0.0), bbox_transform=under, ncol=3,
              fontsize="small")


def _fit_text(t, width):
    """The lines of the text t wrapped at spaces to fit width pixels."""
    w = t.get_window_extent().width
    text = t.get_text()
    if not text or w <= width:
        return
    lines = text.split("\n")
    n = max(12, int(max(len(x) for x in lines) * width / w * 0.95))
    t.set_text("\n".join(y for x in lines for y in (textwrap.wrap(x, n) or [""])))


def _marginal_title(fig):
    """The longest of the titles that fits the width of the figure on one line."""
    width = fig.get_figwidth() * 72.0
    for t in MARGINAL_TITLES:
        if len(t) * 0.55 * 14.4 <= width * 0.95:
            return t
    return MARGINAL_TITLES[-1]


# ---- the call -----------------------------------------------------------------------------------

def makeplot(res, type=None, *, swapaxis=False, ci=False, citype=None, normal=False, tinterval=False,
             percentile=False, bc=False, bca=False, reverse=False, chainstep=None, outcome=None,
             format=None, blocklabels=None, varlabels=False, level=None, twoway=None, ax=None,
             title=None, subtitle=None, figsize=None, **axes_set):
    """Draw the graph of a staged decomposition and return the matplotlib Figure.

    type         "stages", "marginal", "waterfall" or "waterfalllevels"; by default waterfall,
                 and marginal (the only type) after individually
    swapaxis     swap the horizontal and vertical axes (marginal draws horizontal bars, the
                 others vertical ones)
    ci           draw the confidence intervals of the bootstrap results; ci=True with the
                 interval type citype= ("normal", "tinterval", "percentile", "bc" or "bca"), one
                 of the flags normal=, tinterval=, percentile=, bc=, bca=, or ci="bca"; by
                 default the first interval type of the results (res.e["citypes"]). In the
                 waterfalls an interval is drawn around the end of the bar whose length it
                 concerns
    reverse      the stages of the reverse order (results of sequentially reverse)
    chainstep    with chain results, the pair drawn (1, 2, ...); by default the whole chain
                 (the waterfall and the marginal plot only)
    outcome      with outcomes(), the outcome drawn (by default the first)
    format       the number format of the values labelled, a Stata format such as "%9.2f"
                 (by default 3 significant digits)
    blocklabels  the texts of the rows of the blocks: {"demo": "Demographics"}, a list of
                 (block, text) pairs, or the Stata text 'demo "Demographics" geo Region'
    varlabels    show the variable labels of the blocks named after their variables, and of by()
    level        with ci, the confidence level of the intervals (by default that of the results)
    twoway       a dict of the matplotlib arguments below (Stata's twoway options are not
                 available)

    The matplotlib arguments: ax= draws in these Axes (the title and subtitle are then the
    Axes title, the interval note is placed below the Axes), else a new figure is made with
    pyplot; figsize= its size in inches; title= and subtitle= replace the title and the
    subtitle ("" leaves them out); any other keyword is passed to Axes.set (xlabel=, ylabel=,
    xlim=, ylim=, ...), after the graph is drawn.
    """
    e = res.e
    if not e.get("stagetype"):
        raise _err("makeplot requires results estimated with sequentially or individually")
    if twoway not in (None, False):
        if not isinstance(twoway, dict):
            raise _err("makeplot(twoway()) is not available; give the matplotlib arguments ax=, title=, "
                       "subtitle=, figsize= or Axes.set keywords to res.makeplot()")
        tw = dict(twoway)
        ax = tw.pop("ax", ax)
        title = tw.pop("title", title)
        subtitle = tw.pop("subtitle", subtitle)
        figsize = tw.pop("figsize", figsize)
        tw.update(axes_set)
        axes_set = tw
    ci, citype = _interval_type(ci, citype, dict(normal=normal, tinterval=tinterval, percentile=percentile,
                                                 bc=bc, bca=bca))
    type, chainstep, outcome, format = _check_spec(type, chainstep, outcome, format)
    import matplotlib.axes
    for k in axes_set:
        if not hasattr(matplotlib.axes.Axes, "set_" + k):
            raise _err("makeplot(): option %s not allowed" % k)
    if level is not None and not ci:
        raise _err("makeplot(): level= sets the level of the intervals; it requires ci")
    P = _setup(res, type, swapaxis, ci, citype, reverse, chainstep, outcome, format, blocklabels, varlabels,
               level)

    import matplotlib.pyplot as plt
    own = ax is None
    if own:
        try:
            fig, ax = plt.subplots(figsize=figsize, layout="constrained")
        except TypeError:
            fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure
    try:
        _draw(P, fig, ax, own, title, subtitle, axes_set)
    except Exception:
        if own:
            plt.close(fig)
        raise
    return fig


def _draw(P, fig, ax, own, title, subtitle, axes_set):
    """The graph in the Axes ax of the figure fig (own: made by makeplot)."""
    e = P.e
    if P.type == "waterfall":
        rows = _rows_waterfall(P)
        main = "Decomposition of the difference by raking stage"
        sub = P.sub
    elif P.type == "marginal":
        rows = _rows_marginal(P)
        main = _marginal_title(fig)
        sub = P.sub
    elif P.type == "waterfalllevels":
        rows = _rows_levels(P)
        main = _levels_title(P)
        sub = "Reverse order" if P.rev else None
    else:
        rows = None
        main = "Explained and unexplained difference by raking stage"
        sub = P.sub
    what = P.ymeas if P.type == "waterfalllevels" else "difference in %ss" % P.ymeas
    axis_title = _tx(_axis_title(P.yname, what))
    main = _tx(main)
    sub = None if sub is None else _tx(sub)
    if P.horiz:
        ax.set_xlabel(axis_title)
        ax.set_ylabel("")
    else:
        ax.set_ylabel(axis_title)
        ax.set_xlabel("")
    if title is not None:
        main = title
    if subtitle is not None:
        sub = subtitle
    heading = None
    if own:
        if main:
            heading = fig.suptitle(main)
        if sub:
            ax.set_title(sub, fontsize="medium")
    else:
        ax.set_title("\n".join(t for t in (main, sub) if t))
    note = None
    if P.ci is not None:
        note = "%s%% %s intervals" % (_gfmt(P.cilevel), P.cilabel)
    if note and own:
        fig.supxlabel(note, x=0.01, ha="left", fontsize="x-small")
    if rows is None:
        _draw_stages(P, fig, ax)
        _stages_legend(P, fig, ax)
    else:
        _draw_rows(P, fig, ax, rows)
    if note and not own:
        # in Axes of the user's: below the Axes, its tick labels, axis title and legend
        ax.annotate(note, (0, 0), xycoords="axes fraction", xytext=(0, -_below(fig, ax, 2.0)),
                    textcoords="offset points", ha="left", va="top", fontsize="x-small", gid="note",
                    annotation_clip=False)
    if axes_set:
        try:
            ax.set(**axes_set)
        except (AttributeError, TypeError, ValueError) as exc:
            raise _err("makeplot(): %s" % exc)
    # the title over the middle of the plot, as the subtitle; both on lines that fit the figure
    # when centred there
    try:
        _render(fig)
        if heading is not None:
            box = ax.get_position()
            mid = (box.x0 + box.x1) / 2
            room = 2 * min(mid, 1 - mid) * fig.bbox.width * 0.98
            # the title spans the figure, as Stata's title(, span): on lines that fit the
            # figure's width, centred over the plot when it fits there, else moved in from the
            # edge of the figure
            _fit_text(heading, fig.bbox.width * 0.98)
            _fit_text(ax.title, room)
            _render(fig)
            box = ax.get_position()
            half = heading.get_window_extent().width / fig.bbox.width / 2
            heading.set_x(min(max((box.x0 + box.x1) / 2, half + 0.01), 0.99 - half))
        else:
            _fit_text(ax.title, ax.get_window_extent().width * 1.1)
    except Exception:
        pass

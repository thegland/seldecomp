"""seldecomp(): the command."""
import re

import numpy as np
import pandas as pd

from . import _parse as P
from ._engine import RakeOptions
from ._estimate import RESERVED_BLOCKS, estimate
from ._prepare import exclusion_note, fmt_val, group_values_text
from ._results import build_results
from ._stata import fmt_g, plural, reldif
from .errors import SeldecompError

DEFAULTS = {
    "by": None, "focalgroup": None, "chain": False, "ratio": None, "outcomes": None,
    "contvars": None, "insupport": False, "individually": False, "sequentially": False,
    "fixothermargins": False, "reverse": False, "detail": False, "makeplot": None,
    "bayesboot": None, "iterate": 10000, "tolerance": 1e-8, "ctrltolerance": 1e-6,
    "trimhiabs": None, "trimloabs": None, "trimhirel": None, "trimlorel": None,
    "trimsupportoverride": False, "trimfrequency": None, "force": False, "generate": None,
    "covariatedetail": False, "diagnostics": False, "mincell": 10, "noisily": False,
    "cformat": None, "pformat": None, "sformat": None, "varlabels": False, "level": None,
}
# the formats that _display.fmts reads


class Context:
    """What every part of a run needs: the data, the options, the sample."""

    def group_label(self, v):
        t = fmt_val(v)
        lab = (self.by_labels or {}).get(v)
        return "%s (%s)" % (t, lab) if lab not in (None, "") else t


def _text(val):
    """A number as the user would type it in an option."""
    if isinstance(val, str):
        return val.strip()
    f = float(val)
    if f == int(f) and abs(f) < 1e15:
        return str(int(f))
    return ("%.15g" % f)


def _suboption_text(name, val):
    """A suboption list given as a dict, in the Stata form."""
    if name == "contvars":
        d = dict(val)
        tail = []
        if d.pop("recut", False):
            tail.append("recut")
        sample = d.pop("sample", None)
        for k in ("pooled", "focalgroup", "refgroup"):
            if d.pop(k, False):
                sample = k
        if sample and sample != "pooled":
            tail.append(sample)
        body = " ".join("%s %s" % (k, v) for k, v in d.items())
        return body + (", " + " ".join(tail) if tail else "")
    parts = []
    head = ""
    for k, v in val.items():
        if k.startswith("_") or v is None or v is False:
            continue
        if name == "makeplot" and k == "type":
            head = str(v)
        elif v is True:
            parts.append(k)
        elif isinstance(v, dict):
            parts.append("%s(%s)" % (k, " ".join('%s "%s"' % kv for kv in v.items())))
        else:
            parts.append("%s(%s)" % (k, v))
    body = " ".join(parts)
    if name == "makeplot":
        return head + (", " + body if body else "")
    return body


def _contvars_text(d):
    """contvars= as a dict, in the Stata form: the pairs of a variable and its number of
    quantile groups, then the suboptions (the keys given True, and the value of sample)."""
    pairs, words = [], []
    for k, v in d.items():
        if k == "sample":
            if v not in (None, False):
                words.append(str(v))
        elif isinstance(v, (bool, np.bool_)):
            if v:
                words.append(k)
        else:
            pairs.append("%s %s" % (k, v))
    return " ".join(pairs) + (", " + " ".join(words) if words else "")


def _command_text(varlist, wtype, wexp, items):
    """The call written as the Stata command line (stored as e["cmdline"])."""
    text = "seldecomp " + " ".join(varlist.split())
    if wexp is not None:
        text += " [%s=%s]" % ({"pweight": "pw", "aweight": "aw", "iweight": "iw"}.get(wtype, "pw"), wexp)
    opts = _options_text(items)
    return text + (", " + opts if opts else "")


def _options_text(items):
    """Keyword options written as Stata options."""
    opts = []
    for name, val in items:
        if val is None or val is False:
            continue
        if val is True:
            opts.append(name)
        elif isinstance(val, dict):
            opts.append("%s(%s)" % (name, _suboption_text(name, val)))
        elif isinstance(val, (list, tuple)):
            opts.append("%s(%s)" % (name, " ".join(str(v) for v in val)))
        else:
            opts.append("%s(%s)" % (name, val))
    return " ".join(opts)


def _keyword_opts(items):
    """The options given as keywords, as options typed (P.Opt): True is the word alone, a
    number or text the text in its parentheses, a list its words; a dict keeps its value."""
    out = []
    for k, v in items:
        if isinstance(v, (bool, np.bool_)) or (P.spec_name(k, P.MAIN_SPEC) in P.FLAG_OPTIONS
                                                and isinstance(v, (int, float, np.integer, np.floating))):
            v = bool(v)
        if v is None or v is False:
            continue
        if v is True or (isinstance(v, dict) and not v):
            out.append(P.Opt(k))
        elif isinstance(v, dict):
            out.append(P.Opt(k, True, _suboption_text(P.spec_name(k, P.MAIN_SPEC) or k, v), v))
        elif isinstance(v, (list, tuple)):
            out.append(P.Opt(k, True, " ".join(str(x) for x in v)))
        elif isinstance(v, str):
            out.append(P.Opt(k, True, v))
        else:
            out.append(P.Opt(k, True, _text(v)))
    return out


def _is_string(s):
    """Whether a column holds text (a categorical: text categories)."""
    if isinstance(s.dtype, pd.CategoricalDtype):
        return not pd.api.types.is_numeric_dtype(s.cat.categories)
    return P.is_string_series(s)


def _uniq(values, order=None):
    u = pd.unique(np.asarray(values, dtype=object)).tolist()
    return sorted(u, key=order.get) if order is not None else sorted(u)


def _group_pick(t, by, g1, g2, labels, gvl, zwhy):
    """Group_pick: which of the two values (g1 < g2) of a numeric by() the text of
    focalgroup() names: a value (exactly, else the one within a relative 1e-7), else a value
    label -> (1 or 2, a note or None). gvl() and zwhy() give the texts of the messages."""
    x = P.stata_number(t)
    isnum = x == x
    kv = 0
    if isnum:
        if x == g1:
            kv = 1
        elif x == g2:
            kv = 2
        else:
            d1, d2 = reldif(x, g1) < 1e-7, reldif(x, g2) < 1e-7
            kv = 1 if d1 and not d2 else (2 if d2 and not d1 else 0)
    m = [bool(labels) and float(g) == int(g) and labels.get(int(g)) == t for g in (g1, g2)]
    nl = sum(m)
    kl = 1 if m[0] else (2 if m[1] else 0)
    if nl == 2 or (kv and nl == 1 and kl != kv):
        why = "the value label of both values" if nl == 2 else \
            "a value of %s and the value label of the other value" % by
        raise SeldecompError("focalgroup(%s) is ambiguous: %s; its values are %s" % (t, why, gvl()), 198)
    if kv:
        return kv, None
    if kl:
        note = None
        if isnum:
            note = ("note: focalgroup(%s): %s is not a value of %s; it is taken as the value label of %s"
                    % (t, t, by, fmt_g((g1, g2)[kl - 1], 21)))
        return kl, note
    if labels:
        raise SeldecompError("focalgroup(%s) is neither a value nor a value label of %s%s; its values are %s"
                             % (t, by, zwhy(), gvl()), 459)
    if not isnum:
        raise SeldecompError("focalgroup(%s) must be a number; the values of %s are %s" % (t, by, gvl()), 198)
    raise SeldecompError("focalgroup(%s) is not a value of %s%s; its values are %s"
                         % (t, by, zwhy(), gvl()), 459)


def _zero_named(t, by, byvals, zw, labels, bystr):
    """Zero_named: when the text of focalgroup() names a value whose observations were all
    left out for a zero weight, ": " and that reason, else ""."""
    if not zw.any():
        return ""
    zn = np.zeros(len(byvals), dtype=bool)
    if bystr:
        zn = zw & np.array([v == t for v in byvals])
    else:
        vals = _uniq(byvals[zw])
        x = P.stata_number(t)
        if x == x:
            e = [v for v in vals if v == x] or [v for v in vals if reldif(v, x) < 1e-7]
            if len(e) == 1:
                zn |= zw & (byvals == e[0])
        for v in vals:
            if labels and v == int(v) and labels.get(int(v)) == t:
                zn |= zw & (byvals == v)
    nz = int(zn.sum())
    if nz == 0:
        return ""
    gvz = group_values_text(_uniq(byvals[zn]), labels)
    if nz == 1:
        return ": the one observation with %s = %s and nonmissing data has zero weight" % (by, gvz)
    return ": all %d observations with %s = %s and nonmissing data have zero weight" % (nz, by, gvz)


def _zero_group(fg, by, byvals, touse, zw, labels, bystr):
    """Zero_group: by() takes one value; when the observations of the other value all have
    zero weight, the text of that reason (naming the group when it can be told), else None."""
    k = byvals[touse][0] if bystr else min(byvals[touse])
    zo = zw & np.array([v != k for v in byvals])
    nz = int(zo.sum())
    if nz == 0:
        return None
    zv = _uniq(byvals[zo])
    gvz = group_values_text(zv, labels)
    if len(zv) > 1:
        return ("the observations with the other values of %s (%s) and nonmissing data all have zero weight"
                % (by, gvz))
    z = zv[0]
    role = 0
    if fg is None:
        role = 1 if z < k else 2
    elif bystr:
        t = P.unquote(fg)
        zhit = 2 if t == z else (1 if t == k else 0)
        role = 3 - zhit if zhit else 0
    else:
        try:
            kz, _ = _group_pick(P.unquote(fg), by, k, z, labels, lambda: "", lambda: "")
            role = 3 - kz
        except SeldecompError:
            role = 0
    lead = "the one" if nz == 1 else "all %d" % nz
    obs = plural(nz, "observation")
    have = "has" if nz == 1 else "have"
    if role:
        return "%s %s %s (%s = %s) with nonmissing data %s zero weight" % (
            lead, "focal group" if role == 1 else "reference group", obs, by, gvz, have)
    return "%s %s of the other group (%s = %s) with nonmissing data %s zero weight" % (
        lead, obs, by, gvz, have)


def _syntax(data, spec, by, weights, kwitems, value_labels, sub):
    """The command line as the Stata command reads it before its own checks: the options
    refused by -syntax- or by the command's first checks (Reject_nodiagnostics, the replay,
    Comma_check, Top_error), each with Stata's message, in Stata's order -> (line, option
    values, rows selected by if and in, weights or None)."""
    line = P.split_command(spec)
    kwopts = _keyword_opts(([("by", by)] if by is not None else []) + list(kwitems))
    tok_err = None
    try:
        text_opts = P.option_tokens(line.options) if line.options is not None else []
    except SeldecompError as ex:
        tok_err, text_opts = ex, []
    # the same option in the text and as a keyword
    if by is not None and any(P.spec_name(o.word, P.MAIN_SPEC) == "by" for o in text_opts):
        raise SeldecompError("option by() given twice", 198)
    for o in kwopts[1 if by is not None else 0:]:
        name = P.spec_name(o.word, P.MAIN_SPEC)
        if name is not None and any(P.spec_name(t.word, P.MAIN_SPEC) == name for t in text_opts):
            raise SeldecompError("option %s given twice: in the command text and as a keyword" % name, 198)
    if weights is not None and line.qual("weight") is not None:
        raise SeldecompError("weights given twice: in the command text and as weights=", 198)
    opts = text_opts + kwopts
    try:
        mask, w = P.check_qualifiers(line, data)
        ok = tok_err is None
    except SeldecompError as ex:
        mask, w, ok, qual_err = None, None, False, ex
    # Reject_nodiagnostics, then the replay that seldecomp without a variable list is
    if ok:
        P.reject_nodiagnostics(opts)
    if line.zero == "" or line.zero.startswith(","):
        raise SeldecompError("last estimates not found", 301)
    # bayesboot() and makeplot() with empty parentheses are the words alone
    for i, o in enumerate(opts):
        if o.paren and o.arg == "" and o.value is None and re.match(
                r"^(bayes|bayesb|bayesbo|bayesboo|bayesboot|make|makep|makepl|makeplo|makeplot)$", o.word):
            opts[i] = P.Opt(o.word)
    if kwopts:
        kwtext = " ".join(o.text() for o in kwopts)
        line.options = (line.options + " " if line.options is not None else "") + kwtext
        line.zero = line.zero + (" " if P.find_top(line.zero, ",") >= 0 else ", ") + kwtext
    P.comma_check(line, ok, list(data.columns))
    if mask is None:
        raise qual_err
    if line.anything == "":
        raise SeldecompError("varlist required", 100)
    if tok_err is not None:
        raise tok_err
    try:
        vals = P.run_syntax(opts, P.MAIN_SPEC, list(data.columns), required=("by",))
    except SeldecompError as err:
        P.top_error(err, opts, data, mask if sub is None else mask & sub, value_labels)
        raise
    if weights is not None:
        w = P.weight_values(weights, data)
    return line, vals, mask, w


def seldecomp(data, spec, by=None, *, weights=None, subset=None, labels=None, value_labels=None,
              **options):
    """Raking-based reweighting decomposition of group differences in means.

    Decomposes the difference in the mean of an outcome between a focal group and a reference
    group into explained and unexplained components, by reweighting (raking) the focal group
    to match the reference group with respect to the distribution of the raking variables.

    Parameters
    ----------
    data : pandas.DataFrame
        The data.
    spec : str
        ``"y rakespec"``, or the Stata command line (``"y rakespec [if] [in] [weight],
        by(groupvar) [options]"``, with or without the leading word ``seldecomp``).
        *rakespec* lists the raking variables: ``x`` or ``i.x`` (categorical), ``c.x``
        (continuous, with ``contvars``), ``x1##x2`` (interaction), ``(x1 x2)`` or
        ``(name: x1 x2)`` (a set of variables).
    by : str
        The variable defining the focal group and the reference group (lowest value becomes
        the focal group).
    weights : str, optional
        The column of sample weights.
    subset : array-like of bool, optional
        The rows to use.
    labels : dict, optional
        Variable labels: ``{"var": "variable label"}``.
    value_labels : dict, optional
        Value labels: ``{"var": {0: "No", 1: "Yes"}}``.
    **options
        The options of the Stata command, under the same names (abbreviations accepted): a
        flag is ``True`` (``sequentially=True``), an option with an argument takes a number
        or text (``focalgroup=1``, ``ratio="hours"``), an option with suboptions takes the
        Stata text or a dict (``contvars="age 5"``, ``bayesboot={"reps": 200, "rseed": 1}``).
        The options are documented at https://github.com/thegland/seldecomp/tree/main/python

    Returns
    -------
    SeldecompResults
        ``print(res)`` shows the tables; ``res.e``, ``res.b``, ``res.V`` and the matrices
        (``res.bystage``, ...) hold the stored results.

    Raises
    ------
    SeldecompError
        With the message and the return code (``.rc``) of the Stata command.

    Warns
    -----
    SeldecompWarning
        The warning of the Stata command; also shown by ``print(res)``.

    Examples
    --------
    >>> import seldecomp as sd
    >>> res = sd.seldecomp(df, "wage i.collgrad i.south", by="married")
    >>> res = sd.seldecomp(df, "wage (demo: race collgrad) (geo: south smsa) c.age",
    ...                    by="married", contvars="age 5", sequentially=True,
    ...                    bayesboot="reps(200) rseed(123) citable")
    >>> res = sd.seldecomp(df, "wage collgrad south [pw=wt] if age < 40, by(married)")
    """
    if not isinstance(data, pd.DataFrame):
        raise SeldecompError("the first argument must be a pandas DataFrame", 198)
    spec = str(spec)
    bsweights = options.pop("_bsweights", None)
    value_labels = value_labels or {}
    cols = list(data.columns)
    sub = None
    if subset is not None:
        sub = np.asarray(subset)
        if sub.dtype != bool:
            sub = sub != 0
        if len(sub) != len(data):
            raise SeldecompError("subset= must have one value per row of the data", 198)
    line, vals, mask, w = _syntax(data, spec, by, weights, options.items(), value_labels, sub)
    varlist = line.anything
    wq = line.qual("weight")
    if wq is not None:
        wtype, wexp = wq[0], wq[1]
    else:
        wtype, wexp = ("pweight", weights) if weights is not None else (None, None)
    items = [] if by is None else [("by", by)]
    items += [(P.spec_name(k, P.MAIN_SPEC) or k, v) for k, v in options.items()]
    if P.looks_like_command(spec):
        cmdline = spec.strip()
        if not cmdline.startswith("seldecomp"):
            cmdline = "seldecomp " + cmdline
        # keywords given with a command line are written into it
        k = P.find_top(cmdline, ",")
        head, tail = (cmdline, "") if k < 0 else (cmdline[:k].rstrip(), cmdline[k + 1:].strip())
        if weights is not None and wq is None:
            head += " [pw=%s]" % weights
        tail = (tail + " " + _options_text(items)).strip()
        cmdline = head + (", " + tail if tail else "")
    else:
        cmdline = _command_text(varlist, wtype, wexp, items)
    if sub is not None:
        mask = mask & sub

    # ---- the options, as the command reads them
    opt = dict(DEFAULTS)
    for k in P.FLAG_OPTIONS:
        opt[k] = bool(vals.get(k, False))
    for k in ("by", "focalgroup", "contvars", "tolerance", "ctrltolerance", "trimhirel", "trimhiabs",
              "trimlorel", "trimloabs", "trimfrequency", "generate", "cformat", "pformat", "sformat",
              "ratio"):
        opt[k] = vals.get(k)
    opt["iterate"] = vals.get("iterate", 10000)
    opt["mincell"] = vals.get("mincell", 10)
    opt["outcomes"] = list(vals.get("outcomes") or [])
    for k in ("bayesboot", "makeplot"):
        opt[k] = vals[k + "2"] if k + "2" in vals else (True if vals.get(k) else None)

    ctx = Context()
    ctx.df = data
    ctx.cmdline = cmdline
    ctx.value_labels = value_labels
    ctx.labels = labels or {}
    ctx.messages = []
    ctx.tracelines = []

    def trace(text):
        # noisily: each line is shown as the raking produces it
        print(text)
        ctx.tracelines.append(text)

    ctx.trace = trace

    # ---- the checks of the options, in the order of the Stata command ------------
    # the display formats
    for k in ("cformat", "pformat", "sformat"):
        if opt[k] is not None:
            opt[k] = P.check_display_format(k, opt[k])
    seq, ind = opt["sequentially"], opt["individually"]
    if seq and ind:
        raise SeldecompError("sequentially and individually may not be combined", 198)
    if ind:
        for k in ("reverse", "detail"):
            if opt[k]:
                raise SeldecompError("%s may not be combined with individually" % k, 198)
    if opt["fixothermargins"] and not (seq or ind):
        raise SeldecompError("fixothermargins requires sequentially or individually", 198)
    gen, genreplace = None, False
    if opt["generate"] is not None:
        g = str(opt["generate"])
        k = g.find(",")
        gname, grest = (g, "") if k < 0 else (g[:k], g[k:])
        gopt = grest.strip()[1:].strip()
        if len(gname.split()) != 1 or (grest.strip()[:1] not in ("", ",")) or gopt not in ("", "replace"):
            raise SeldecompError("generate(): type generate(newvar) or generate(newvar, replace)", 198)
        gen, genreplace = gname.strip(), gopt == "replace"
    ctx.gen, ctx.genreplace = gen, genreplace
    if opt["chain"]:
        if opt["focalgroup"] is not None:
            raise SeldecompError("focalgroup() may not be combined with chain", 198)
        for bad, what in ((ind, "individually"), (opt["reverse"], "reverse"), (opt["detail"], "detail"),
                          (opt["outcomes"], "outcomes()"), (gen, "generate()"),
                          (opt["covariatedetail"], "covariatedetail")):
            if bad:
                raise SeldecompError("%s may not be combined with chain" % what, 198)
    if wtype == "fweight":
        raise SeldecompError("fweights not allowed\nexpand the data first: expand %s" % wexp, 101)
    if opt["makeplot"] is not None and not (seq or ind):
        raise SeldecompError("makeplot requires sequentially or individually", 198)
    if opt["reverse"] and not (seq or ind):
        raise SeldecompError("reverse requires sequentially", 198)
    if opt["detail"] and not (seq or ind):
        raise SeldecompError("detail requires sequentially", 198)
    mp = None
    if opt["makeplot"] is not None:
        mp = P.parse_makeplot(_suboption_text("makeplot", opt["makeplot"])
                              if isinstance(opt["makeplot"], dict) else opt["makeplot"])
        if mp["chainstep"] is not None and not opt["chain"]:
            raise SeldecompError("makeplot(chainstep()) requires chain", 198)
        if ind and mp["type"] not in (None, "marginal"):
            raise SeldecompError("makeplot(%s) not allowed with individually; use marginal" % mp["type"], 198)
        if mp["reverse"] and not opt["reverse"]:
            raise SeldecompError("makeplot(reverse) requires sequentially reverse", 198)
        if mp["ci"] and opt["bayesboot"] is None:
            raise SeldecompError("makeplot(ci) requires bayesboot()", 198)
        if isinstance(mp["twoway"], str):
            raise SeldecompError("makeplot(): twoway() takes Stata graph options, which are not read "
                                 "here; give matplotlib arguments instead, as in "
                                 "res.makeplot(title=..., figsize=...)", 198)
    if opt["mincell"] < 0:
        raise SeldecompError("mincell() must be a nonnegative integer", 198)
    # the raking options, with the text they are recorded as
    rtext = {}
    for k in ("tolerance", "ctrltolerance", "trimhirel", "trimhiabs", "trimlorel", "trimloabs"):
        if opt[k] is not None:
            t = str(opt[k]).strip()
            f = P.stata_number(t) if len(t.split()) == 1 else float("nan")
            if f != f:
                raise SeldecompError("%s() must be a number" % k, 198)
            if f <= 0:
                raise SeldecompError("%s() must be a positive number" % k, 198)
            # (a plus sign is dropped from the text that is recorded)
            rtext[k] = t[1:] if t.startswith("+") else t
            opt[k] = f
    if opt["iterate"] < 1:
        raise SeldecompError("iterate() must be a positive integer", 198)
    for k, dflt in (("tolerance", "1e-8"), ("ctrltolerance", "1e-6")):
        if opt[k] is None:
            rtext[k] = dflt
            opt[k] = float(dflt)
    anytrim = any(opt[k] is not None for k in ("trimhirel", "trimhiabs", "trimlorel", "trimloabs"))
    if opt["trimhiabs"] is not None and opt["trimloabs"] is not None and opt["trimhiabs"] <= opt["trimloabs"]:
        raise SeldecompError("trimhiabs() must be greater than trimloabs()", 198)
    if opt["trimhirel"] is not None and opt["trimlorel"] is not None and opt["trimhirel"] <= opt["trimlorel"]:
        raise SeldecompError("trimhirel() must be greater than trimlorel()", 198)
    if opt["trimfrequency"] is not None:
        tf = str(opt["trimfrequency"])
        if tf not in ("once", "sometimes", "often"):
            raise SeldecompError("trimfrequency() must be once, sometimes, or often", 198)
        if not anytrim:
            raise SeldecompError("trimfrequency() requires a trimming option", 198)
        opt["trimfrequency"] = tf
    if opt["trimsupportoverride"] and opt["trimloabs"] is None and opt["trimlorel"] is None:
        raise SeldecompError("trimsupportoverride requires trimloabs() or trimlorel()", 198)
    ctx.ropt = RakeOptions(iterate=opt["iterate"], tolerance=opt["tolerance"],
                           ctrltolerance=opt["ctrltolerance"],
                           trimhirel=opt["trimhirel"], trimlorel=opt["trimlorel"], trimhiabs=opt["trimhiabs"],
                           trimloabs=opt["trimloabs"], trimfrequency=opt["trimfrequency"] or "sometimes",
                           trimsupportoverride=opt["trimsupportoverride"], mincell=opt["mincell"])
    # e(rakeopts): the relative bounds and trimfrequency(), then the absolute bounds
    ro = "iterate(%d) tolerance(%s) ctrltolerance(%s)" % (opt["iterate"], rtext["tolerance"],
                                                          rtext["ctrltolerance"])
    for k in ("trimhirel", "trimlorel"):
        if opt[k] is not None:
            ro += " %s(%s)" % (k, rtext[k])
    if opt["trimfrequency"] is not None:
        ro += " trimfrequency(%s)" % opt["trimfrequency"]
    for k in ("trimhiabs", "trimloabs"):
        if opt[k] is not None:
            ro += " %s(%s)" % (k, rtext[k])
    if opt["trimsupportoverride"]:
        ro += " trimsupportoverride"
    ctx.rakeopts = ro
    # bayesboot()
    bs = None
    bsgen, bsgenreplace = None, False
    if opt["bayesboot"] is not None:
        given = opt["bayesboot"]
        if isinstance(given, dict):
            # a dict: read as the text it stands for
            given = _suboption_text("bayesboot", given)
        bs = P.parse_bayesboot(given, cols)
        if isinstance(opt["bayesboot"], dict):
            bs.update({k: v for k, v in opt["bayesboot"].items() if k.startswith("_")})
        if mp is not None and mp["citype"] == "bca" and not bs["bca"]:
            raise SeldecompError("makeplot(bca): BCa intervals require bayesboot(bca)", 198)
        if bs["generate"] not in (None, ""):
            bsgen, bsgenreplace = P.check_bs_generate(bs["generate"])
        if bs["saving"] not in (None, ""):
            P.check_bs_saving(bs["saving"])
        if bsweights is not None:
            bs["_weights"] = bsweights
    ctx.bs = bs
    # contvars()
    ctx.contlist, ctx.contsample, ctx.recut, ctx.conttext = [], "pooled", False, ""
    if opt["contvars"] is not None:
        given = opt["contvars"]
        if isinstance(given, dict):
            given = _contvars_text(given)
        ctx.contlist, ctx.contsample, ctx.recut, ctx.conttext = P.parse_contvars(given, data)
    if ctx.recut and bs is not None and opt["insupport"]:
        raise SeldecompError("contvars(..., recut) may not be combined with insupport and bayesboot()", 198)
    if ctx.recut and bs is not None and opt["chain"]:
        raise SeldecompError("contvars(..., recut) may not be combined with chain and bayesboot()", 198)
    ctx.makeplot = mp
    ctx.level = float(bs["level"]) if (bs is not None and bs["level"] is not None) else 95.0
    ctx.opt = opt

    # ---- depvar and the estimation sample --------------------------------------
    toks = varlist.split(None, 1)
    depvar = toks[0]
    rakespec = toks[1] if len(toks) > 1 else ""
    if rakespec.strip() == "":
        raise SeldecompError("at least one raking variable is required", 102)
    if "." in depvar or "#" in depvar:
        raise SeldecompError("dependent variable %s: factor-variable and time-series operators are "
                             "not allowed" % depvar, 101)
    try:
        dl = P.unab(depvar, cols)
    except SeldecompError:
        raise SeldecompError("dependent variable %s not found" % depvar, 111)
    if len(dl) != 1:
        raise SeldecompError("one dependent variable allowed", 103)
    depvar = dl[0]
    if _is_string(data[depvar]):
        raise SeldecompError("dependent variable %s must be numeric" % depvar, 108)
    ylist = [depvar]
    for v in opt["outcomes"]:
        if _is_string(data[v]):
            raise SeldecompError("outcomes(): %s must be numeric" % v, 108)
        if v in ylist:
            if v == depvar:
                raise SeldecompError("outcomes(): %s is the dependent variable; list only the other "
                                     "outcomes" % v, 198)
            raise SeldecompError("outcomes(): %s listed more than once" % v, 198)
        ylist.append(v)
    if mp is not None and mp["outcome"] is not None and mp["outcome"] not in ylist:
        if len(ylist) > 1:
            raise SeldecompError("makeplot(outcome()): %s is not one of the outcomes (%s)"
                                 % (mp["outcome"], " ".join(ylist)), 198)
        raise SeldecompError("makeplot(outcome()): %s is not the outcome (%s)" % (mp["outcome"], depvar), 198)
    if opt["ratio"] and _is_string(data[opt["ratio"]]):
        raise SeldecompError("ratio(): %s must be a numeric variable" % opt["ratio"], 108)
    byname = opt["by"]
    ctx.by = byname
    ctx.ylist = ylist
    ctx.depvar = depvar
    ctx.rakespec = rakespec

    def notmissing(col):
        s = data[col]
        if _is_string(s):
            return P._column_value(data, col)[1] != ""
        return s.notna().to_numpy()

    touse = mask.copy()
    touse &= notmissing(depvar)
    touse &= notmissing(byname)
    if w is None:
        w = np.ones(len(data))
    ctx.wtype, ctx.wexp = wtype, (None if wexp is None else str(wexp).strip())

    blocks, terms = P.parse_rakespec(rakespec, cols)
    for t in terms:
        if depvar in t.vars:
            raise SeldecompError("the dependent variable %s cannot be a raking variable" % depvar, 198)
        for v in opt["outcomes"]:
            if v in t.vars:
                raise SeldecompError("outcomes(): %s cannot be a raking variable" % v, 198)
    contnames = [v for v, q in ctx.contlist]
    cvars = []
    for t in terms:
        for v, kind in zip(t.vars, t.kinds):
            if kind == "c":
                if v not in cvars:
                    cvars.append(v)
                if v not in contnames:
                    raise SeldecompError("c.%s requires an entry in contvars(), e.g. contvars(%s 4)"
                                         % (v, v), 198)
            elif v in contnames:
                raise SeldecompError("%s is in contvars() but is not marked c.%s" % (v, v), 198)
            touse &= notmissing(v)
    for v in contnames:
        if v not in cvars:
            raise SeldecompError("%s is in contvars() but is not a raking variable" % v, 198)
    for t in terms:
        if byname in t.vars:
            raise SeldecompError("the by() variable %s cannot be a raking variable" % byname, 198)
    # the continuous variables in the order of the raking list
    ctx.contlist = sorted(ctx.contlist, key=lambda vq: cvars.index(vq[0]))
    ctx.terms, ctx.blocks = terms, blocks

    byraw = data[byname]
    byorder = None
    if isinstance(byraw.dtype, pd.CategoricalDtype):
        # a pandas categorical: its categories in their order
        byorder = {c: i for i, c in enumerate(byraw.cat.categories)}
        byvals = byraw.astype(object).to_numpy()
        bystr = not pd.api.types.is_numeric_dtype(byraw.cat.categories)
    elif _is_string(byraw):
        byvals = P._column_value(data, byname)[1]
        bystr = True
    else:
        byvals = byraw.to_numpy()
        bystr = False
    ctx.byvals = byvals
    ctx.by_labels = ctx.value_labels.get(byname, {})
    bylab = ctx.by_labels

    # the other exclusions, counted by group: a missing outcomes() or ratio() value, a missing
    # weight, a missing bayesboot() cluster id; then the zero weights
    xl = touse.copy()
    for v in ylist[1:] + ([opt["ratio"]] if opt["ratio"] else []):
        xl &= notmissing(v)
    xlw = touse & ~xl
    touse = touse & xl
    xmw = touse & np.isnan(w)
    touse = touse & ~xmw
    xmc = np.zeros(len(data), dtype=bool)
    if bs is not None and bs["cluster"]:
        xmc = touse & ~notmissing(bs["cluster"])
        touse = touse & ~xmc
    zw = np.zeros(len(data), dtype=bool)
    if wexp is not None:
        if np.any(touse & (w < 0)):
            raise SeldecompError("negative weights encountered", 402)
        zw = touse & (w == 0)
        touse = touse & ~zw
    ctx.w = w
    ctx.N_zerowt = int(zw.sum())

    # the block names (with sequentially, the equation names of e(b))
    staged = seq or ind
    bnames = [b.name if P.valid_name(b.name) else "block%d" % (i + 1) for i, b in enumerate(blocks)]
    if staged:
        for i, nm in enumerate(bnames):
            if nm in bnames[:i]:
                raise SeldecompError("duplicate block name %s" % nm, 198)
            if nm in RESERVED_BLOCKS or (ind and nm == "_remainder"):
                raise SeldecompError("%s cannot be used as a block name" % nm, 198)
        if opt["fixothermargins"]:
            bv = [[v for t in b.terms for v in terms[t].vars] for b in blocks]
            for i in range(len(blocks)):
                for j in range(i + 1, len(blocks)):
                    common = [v for v in bv[i] if v in bv[j]]
                    if common:
                        raise SeldecompError("with fixothermargins, a variable may be in one block only: "
                                             "%s is in blocks %s and %s" % (common[0], bnames[i], bnames[j]),
                                             198)
            if opt["detail"]:
                for b in blocks:
                    for i, t1 in enumerate(b.terms):
                        for t2 in b.terms[i + 1:]:
                            common = [v for v in terms[t1].vars if v in terms[t2].vars]
                            if common:
                                raise SeldecompError("with detail and fixothermargins, a variable may be in "
                                                     "one term only: %s is in %s and %s"
                                                     % (common[0], terms[t1].name, terms[t2].name), 198)
    if mp is not None and mp["blocklabels"] is not None:
        mp["blocklabels"] = P.check_blocklabels(mp["blocklabels"], bnames)

    # ---- the two groups of by() --------------------------------------------------
    if not touse.any():
        nzw = ctx.N_zerowt
        if nzw > 0:
            raise SeldecompError("no observations: %s zero weight"
                                 % ("the one observation with nonmissing data has" if nzw == 1 else
                                    "all %d observations with nonmissing data have" % nzw), 2000)
        raise SeldecompError("no observations", 2000)
    if opt["chain"] and bystr:
        raise SeldecompError("chain requires a numeric by() variable", 198)
    vals_in = _uniq(byvals[touse], byorder)
    ng = len(vals_in)
    gvl = lambda: group_values_text(vals_in, bylab)  # noqa: E731
    focal = ref = None
    if opt["chain"]:
        if ng < 2:
            raise SeldecompError("by(%s) takes one value (%s); chain requires two or more"
                                 % (byname, " ".join(fmt_val(v) for v in vals_in)), 459)
        if mp is not None:
            if mp["chainstep"] is not None:
                if mp["chainstep"] > ng - 1:
                    raise SeldecompError("makeplot(chainstep()) must be 1 to %d" % (ng - 1), 198)
            elif mp["type"] in ("stages", "waterfalllevels"):
                raise SeldecompError("makeplot(%s) requires chainstep(#) with chain" % mp["type"], 198)
    elif ng != 2:
        if ng == 1:
            zwhy = None
            if ctx.N_zerowt > 0:
                zwhy = _zero_group(opt["focalgroup"], byname, byvals, touse, zw, bylab, bystr)
            txt = gvl() if bystr else fmt_val(vals_in[0])
            zwhy = zwhy or "the other group has no observations with nonmissing data and a positive weight"
            raise SeldecompError("by(%s) takes one value (%s): %s" % (byname, txt, zwhy), 459)
        raise SeldecompError("by(%s) takes %d values (%s); it must take two" % (byname, ng, gvl()), 459)
    else:
        focal, ref = vals_in
        fg = opt["focalgroup"]
        if fg is not None:
            t = P.unquote(fg)
            zwhy = lambda: _zero_named(t, byname, byvals, zw, bylab, bystr)  # noqa: E731
            if bystr:
                if t not in vals_in:
                    raise SeldecompError("focalgroup(%s) is not a value of %s%s; its values are %s"
                                         % (t, byname, zwhy(), gvl()), 459)
                kfoc = vals_in.index(t) + 1
            else:
                kfoc, note = _group_pick(t, byname, vals_in[0], vals_in[1], bylab, gvl, zwhy)
                if note:
                    from ._display import para
                    ctx.messages += para(note, 0, 6)
            focal, ref = (vals_in[0], vals_in[1]) if kfoc == 1 else (vals_in[1], vals_in[0])
    vals = vals_in
    ctx.vals = vals

    def by_count(m, v):
        return int(np.sum(m & np.array([x == v for x in byvals])))

    ctx.xnotes = []
    if opt["chain"]:
        # chain: the observations left out, counted at each value of by()
        def by_text(m):
            u = sorted(pd.unique(byvals[m]).tolist())
            if len(u) > 6:
                return "in %d values of %s" % (len(u), byname)
            return "%s = " % byname + "; ".join(
                "%s: %s" % (fmt_val(v), "{:,}".format(int(np.sum(m & np.array([x == v for x in byvals])))))
                for v in u)

        for m, what in ((xlw, "a missing value of %s (the ratio() variable)" % opt["ratio"]),
                        (xmw, "a missing weight"), (xmc, "a missing bayesboot() cluster id"),
                        (zw, "zero weight")):
            n = int(m.sum())
            if n:
                ctx.xnotes.append("%s %s with %s excluded (%s)" % (
                    "{:,}".format(n), "observation" if n == 1 else "observations", what, by_text(m)))
    else:
        for m, what, tail in ((xlw, None, ""), (xmw, "a missing weight", ""),
                              (xmc, "a missing bayesboot() cluster id", "")):
            n = int(m.sum())
            if n == 0:
                continue
            if what is None:
                if opt["outcomes"] and opt["ratio"]:
                    what = "a missing value of an outcomes() or the ratio() variable"
                elif opt["outcomes"]:
                    what = "a missing value of an outcomes() variable"
                else:
                    what = "a missing value of %s (the ratio() variable)" % opt["ratio"]
                if opt["outcomes"]:
                    tail = "; every outcome uses the same observations"
            ctx.xnotes.append(exclusion_note(n, by_count(m, focal), by_count(m, ref), what, byname, tail))
        if ctx.N_zerowt > 0:
            ctx.xnotes.append(exclusion_note(ctx.N_zerowt, by_count(zw, focal), by_count(zw, ref),
                                             "zero weight", byname))
    ctx.fixed_cuts = None
    ctx.touse = touse

    # contvars(): the observations the quantile groups are formed on (a chain: in each pair)
    if not opt["chain"] and ctx.contlist:
        cond = touse
        if ctx.contsample != "pooled":
            isfocal = np.array([x == focal for x in byvals])
            cond = touse & (isfocal if ctx.contsample == "focalgroup" else ~isfocal)
        nbin = int(cond.sum())
        for v, q in ctx.contlist:
            if nbin < q:
                raise SeldecompError("contvars(): too few observations to form %d quantile groups of %s"
                                     % (q, v), 2001)

    # ---- the names the results take: outcomes() equations, generate() variables. They are
    # checked when the sample is settled (after insupport), before the first raking
    def names_check():
        if len(ylist) > 1:
            eqs = []
            for v in ylist:
                veqs = [v]
                if ind:
                    veqs += [v + "_marg_expl", v + "_raked_focal"]
                elif seq:
                    veqs += [v + "_marg_expl", v + "_cum_expl", v + "_cum_unexpl", v + "_raked_focal"]
                if opt["detail"]:
                    veqs.append(v + "_detail")
                if opt["reverse"]:
                    veqs += [v + "_rev_marg_expl", v + "_rev_cum_expl", v + "_rev_cum_unexpl",
                             v + "_rev_raked_focal"]
                for q in veqs:
                    if q in eqs:
                        raise SeldecompError("outcomes(): two outcomes give the equation name %s; rename one"
                                             % q, 198)
                eqs += veqs
        wvars = [t for t in re.sub(r"[^A-Za-z0-9_]", " ", str(wexp or "")).split() if t in cols]
        rakevars = [v for t in terms for v in t.vars]
        clust = bs["cluster"] if bs is not None else None
        inputs = (("the dependent variable", [depvar]), ("an outcomes() variable", opt["outcomes"]),
                  ("the ratio() variable", [opt["ratio"]] if opt["ratio"] else []),
                  ("the by() variable", [byname]),
                  ("a contvars() variable", contnames), ("a raking variable", rakevars),
                  ("the weight variable", wvars), ("the bayesboot(cluster()) variable", [clust] if clust else []))
        gennames = []
        if gen is not None:
            nstages = len(blocks) + 1 if ind else (len(blocks) if seq else 1)
            gennames = [gen] + ["%s_s%d" % (gen, s) for s in range(1, nstages)]
            for gn in gennames:
                if not P.valid_name(gn):
                    raise SeldecompError("generate(): %s is not a valid variable name" % gn, 198)
                role = next((r for r, vl in inputs if gn in vl), None)
                if role is not None:
                    raise SeldecompError("generate(): %s is %s; choose another name" % (gn, role), 198)
                if not genreplace and gn in cols:
                    raise SeldecompError("variable %s already defined" % gn, 110)
        if bsgen is not None:
            for r in range(1, bs["reps"] + 1):
                nm = "%s_%d" % (bsgen, r)
                if not bsgenreplace:
                    if not P.valid_name(nm):
                        raise SeldecompError("%s invalid name" % nm, 198)
                    if nm in cols:
                        raise SeldecompError("variable %s already defined" % nm, 110)
                else:
                    role = next((rr for rr, vl in inputs if nm in vl), None)
                    if role is not None:
                        raise SeldecompError("bayesboot(generate(%s, replace)) would overwrite %s, %s; choose "
                                             "another name" % (bsgen, nm, role), 198)
                if nm in gennames:
                    raise SeldecompError("generate() and bayesboot(generate()) would both create %s" % nm, 198)

    ctx.before_raking = names_check

    if opt["chain"]:
        from ._chain import estimate_chain
        fitobj = estimate_chain(ctx)
    else:
        fitobj = estimate(ctx, touse, focal, ref)
    res = build_results(ctx, fitobj)
    if bs is not None:
        from ._boot import bootstrap
        bootstrap(ctx, fitobj, res)
    res._finish(data)
    if mp is not None:
        # the results exist whatever happens to the graph: an error carries them along
        try:
            res.figure = res.makeplot(varlabels=opt["varlabels"],
                                      **{k: v for k, v in mp.items() if v not in (None, False)})
        except SeldecompError as ex:
            ex.results = res
            raise
    return res


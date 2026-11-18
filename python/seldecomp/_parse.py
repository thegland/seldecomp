"""Parsing of the Stata-style command text as the Stata command reads it: the command line
(the variable list, the if and in qualifiers, the weight, the options, with the messages of
Stata's -syntax-), Stata expressions, the list of raking variables, and the suboptions of
contvars(), bayesboot() and makeplot()."""
import math
import os
import re
import warnings

import numpy as np
import pandas as pd

from ._stata import numeric_format
from .errors import SeldecompError

# ---- the option lists ----------------------------------------------------------------------------
# Each list as the Stata command declares it to -syntax-, in that order: (name, shortest
# abbreviation, kind). A flag and a passthru of the same name take the word typed without and
# with parentheses (bayesboot, bayesboot(...)).
MAIN_SPEC = (
    ("by", 2, "varname"), ("focalgroup", 5, "string"), ("contvars", 4, "asis"),
    ("sequentially", 3, "flag"), ("individually", 3, "flag"), ("fixothermargins", 3, "flag"),
    ("insupport", 5, "flag"), ("bayesboot", 5, "flag"), ("bayesboot", 5, "passthru"),
    ("iterate", 4, "integer"), ("tolerance", 3, "string"), ("ctrltolerance", 7, "string"),
    ("trimhirel", 9, "string"), ("trimhiabs", 9, "string"), ("trimlorel", 9, "string"),
    ("trimloabs", 9, "string"), ("trimfrequency", 8, "string"), ("trimsupportoverride", 19, "flag"),
    ("generate", 3, "string"), ("covariatedetail", 3, "flag"), ("makeplot", 4, "flag"),
    ("makeplot", 4, "passthru"), ("mincell", 3, "integer"), ("noisily", 3, "flag"),
    ("diagnostics", 4, "flag"), ("force", 5, "flag"), ("cformat", 7, "passthru"),
    ("pformat", 7, "passthru"), ("sformat", 7, "passthru"), ("varlabels", 4, "flag"),
    ("ratio", 3, "varname"), ("outcomes", 3, "varlist"), ("reverse", 3, "flag"),
    ("detail", 3, "flag"), ("chain", 5, "flag"),
)
BAYESBOOT_SPEC = (
    ("reps", 1, "integer"), ("cluster", 2, "varname"), ("rseed", 5, "string"), ("bca", 3, "flag"),
    ("ties", 4, "flag"), ("saving", 3, "asis"), ("generate", 3, "string"), ("level", 1, "cilevel"),
    ("citable", 3, "flag"), ("dots", 4, "integer"), ("nodots", 6, "flag"), ("normal", 4, "flag"),
    ("tinterval", 4, "flag"), ("percentile", 4, "flag"), ("bc", 2, "flag"),
)
CONTVARS_SPEC = (("pooled", 6, "flag"), ("focalgroup", 10, "flag"), ("refgroup", 8, "flag"),
                 ("recut", 5, "flag"))
MAKEPLOT_SPEC = (
    ("swapaxis", 4, "flag"), ("ci", 2, "flag"), ("reverse", 3, "flag"), ("chainstep", 5, "string"),
    ("outcome", 3, "string"), ("format", 3, "string"), ("blocklabels", 5, "asis"),
    ("twoway", 2, "asis"), ("normal", 4, "flag"), ("tinterval", 4, "flag"),
    ("percentile", 4, "flag"), ("bc", 2, "flag"), ("bca", 3, "flag"),
)
# option name -> shortest accepted abbreviation
MAIN_OPTIONS = {n: m for n, m, k in MAIN_SPEC}
FLAG_OPTIONS = {n for n, m, k in MAIN_SPEC if k == "flag" and n not in ("bayesboot", "makeplot")}
BAYESBOOT_OPTIONS = {n: m for n, m, k in BAYESBOOT_SPEC}
BAYESBOOT_FLAGS = {n for n, m, k in BAYESBOOT_SPEC if k == "flag"}
MAKEPLOT_OPTIONS = {n: m for n, m, k in MAKEPLOT_SPEC}
MAKEPLOT_FLAGS = {n for n, m, k in MAKEPLOT_SPEC if k == "flag"}
PLOT_TYPES = ("stages", "marginal", "waterfall", "waterfalllevels")
CONTVARS_FLAGS = {n: m for n, m, k in CONTVARS_SPEC}
CI_TYPES = ("normal", "tinterval", "percentile", "bc", "bca")

# a name: a letter (Unicode letters too, as in Stata 14 and later) or _, then letters, digits, _
_NAME = re.compile(r"^[^\W\d]\w*$")
# a number as Stata reads one (confirm number)
_NUMBER = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")
_INTEGER = re.compile(r"^-?\d+$")


# ---- text --------------------------------------------------------------------------------------
def split_top(text, seps=" \t"):
    """Split text at the separators that are outside parentheses, brackets and double
    quotes. Empty pieces are dropped."""
    out, cur, depth, quote = [], [], 0, False
    for ch in text:
        if quote:
            cur.append(ch)
            if ch == '"':
                quote = False
            continue
        if ch == '"':
            quote = True
            cur.append(ch)
        elif ch in "([":
            depth += 1
            cur.append(ch)
        elif ch in ")]":
            depth -= 1
            cur.append(ch)
        elif depth == 0 and ch in seps:
            if cur:
                out.append("".join(cur))
                cur = []
        else:
            cur.append(ch)
    if cur:
        out.append("".join(cur))
    return out


def find_top(text, ch):
    """Position of the first ch outside parentheses, brackets and quotes, or -1."""
    depth, quote = 0, False
    for i, c in enumerate(text):
        if quote:
            if c == '"':
                quote = False
            continue
        if c == '"':
            quote = True
        elif c in "([":
            depth += 1
        elif c in ")]":
            depth -= 1
        elif depth == 0 and c == ch:
            return i
    return -1


def match_name(word, table, what="option"):
    """The full name in table that word abbreviates (at least the table's minimum length; lower
    case, as Stata reads option names)."""
    w = word
    for full, minlen in table.items():
        if len(w) >= minlen and full.startswith(w):
            return full
    raise SeldecompError("%s %s not allowed" % (what, word), 198)


def unquote(s):
    s = str(s).strip()
    if len(s) >= 2 and s[0] == '"' and s[-1] == '"':
        return s[1:-1]
    return s


def is_string_series(s):
    return pd.api.types.is_string_dtype(s) or pd.api.types.is_object_dtype(s)


def and_list(words):
    """"a", "a and b", "a, b and c"."""
    words = list(words)
    if len(words) <= 1:
        return "".join(words)
    return ", ".join(words[:-1]) + " and " + words[-1]


def bind_tokens(text):
    """The words of text as gettoken with bind reads them: a parenthesized group or a quoted
    string stays in the word it is part of."""
    out, cur, depth, quote = [], [], 0, False
    for ch in text:
        if quote:
            cur.append(ch)
            quote = ch != '"'
            continue
        if ch == '"':
            quote = True
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif depth <= 0 and ch.isspace():
            if cur:
                out.append("".join(cur))
                cur = []
            continue
        cur.append(ch)
    if cur:
        out.append("".join(cur))
    return out


def stata_number(text):
    """The number text holds as confirm number reads it, else nan."""
    t = str(text).strip()
    if not _NUMBER.match(t):
        return float("nan")
    f = float(t)
    return f if math.isfinite(f) and abs(f) < 8.988e307 else float("nan")


# ---- variable lists ------------------------------------------------------------------------------
def valid_name(s):
    return _NAME.match(s) is not None and len(s) <= 32


def unab(text, columns):
    """The variables a Stata varlist names, as unab expands it: names, the wildcards * ? ~ and
    ranges a-b (in the order of the columns); no abbreviations. Raises "x invalid name" (198)
    or "variable x not found" (111)."""
    cols = [c for c in columns if isinstance(c, str)]
    if str(text).strip() in cols:
        return [str(text).strip()]
    out = []
    for tok in str(text).split():
        if tok in cols:
            out.append(tok)
            continue
        if "-" in tok and not tok.startswith("-") and not tok.endswith("-"):
            a, b = tok.split("-", 1)
            if _NAME.match(a) and _NAME.match(b):
                for x in (a, b):
                    if x not in cols:
                        raise SeldecompError("variable %s not found" % x, 111)
                i, j = cols.index(a), cols.index(b)
                if i > j:
                    raise SeldecompError("%s invalid name" % tok, 198)
                out += cols[i:j + 1]
                continue
        if re.search(r"[*?~]", tok) and re.match(r"^[\w*?~]+$", tok) and not tok[0].isdigit():
            rx = re.compile("^" + re.escape(tok).replace(r"\*", ".*").replace(r"\?", ".")
                            .replace(r"\~", ".*") + "$")
            hits = [c for c in cols if rx.match(c)]
            if not hits:
                raise SeldecompError("variable %s not found" % tok, 111)
            out += hits
            continue
        if not valid_name(tok):
            raise SeldecompError("%s invalid name" % tok, 198)
        if tok not in cols:
            raise SeldecompError("variable %s not found" % tok, 111)
        out.append(tok)
    return out


# ---- -syntax- --------------------------------------------------------------------------------------
class Opt:
    """An option as typed (or given as a keyword): its word, whether it has parentheses, the
    text inside them, and the value of a keyword."""

    def __init__(self, word, paren=False, arg=None, value=None):
        self.word = word
        self.paren = paren
        self.arg = arg
        self.value = value

    def shown(self):
        """The option as -syntax- names it in its messages."""
        return self.word + ("()" if self.paren else "")

    def text(self):
        return self.word + ("(%s)" % self.arg if self.paren else "")


def option_tokens(text):
    """The options of a list as -syntax- reads them (a word, and the text in the parentheses
    that follow it). A comma followed by more options is refused as -syntax- refuses it."""
    out, i, n = [], 0, len(text)
    while i < n:
        c = text[i]
        if c.isspace():
            i += 1
            continue
        if c == ",":
            rest = text[i + 1:].strip()
            if rest == "":
                break
            raise SeldecompError("invalid '%s'" % re.match(r"\w+|\S", rest).group(0), 198)
        if c == '"':
            j = text.find('"', i + 1)
            j = n if j < 0 else j + 1
            out.append(Opt(text[i:j][:-1] if text[j - 1] == '"' and j - i > 1 else text[i:j]))
            i = j
            continue
        m = re.match(r"\w+", text[i:])
        if not m:
            j = i
            while j < n and not text[j].isspace():
                j += 1
            out.append(Opt(text[i:j]))
            i = j
            continue
        word = m.group(0)
        i += len(word)
        if i < n and text[i] == "(":
            depth, quote, j = 0, False, i
            while j < n:
                ch = text[j]
                if quote:
                    quote = ch != '"'
                elif ch == '"':
                    quote = True
                elif ch == "(":
                    depth += 1
                elif ch == ")":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            out.append(Opt(word, True, text[i + 1:j]))
            i = j + 1
        else:
            out.append(Opt(word))
    return out


def _abbrev(word, name, minlen):
    return len(word) >= minlen and name.startswith(word)


def _fits(o, kind):
    return (not o.paren) if kind == "flag" else o.paren


def _claims(opts, spec, loose):
    """Which entry of spec each option takes: {option index: (spec index, state)}, state "ok",
    "wrong" (the other form: a flag with parentheses, an option without them) or "off" (the
    no-form of a flag). An option takes the first word of its name, whatever its form; a flag
    and a passthru of the same name (bayesboot, bayesboot()) take the first word of their form
    first; the no-forms of the flags are matched last. With loose (an option list ending in *,
    as Top_error and Subopt_error read it), a word of the other form is left over."""
    claim, taken = {}, set()
    twice = {e[0] for e in spec if sum(1 for f in spec if f[0] == e[0]) > 1}
    passes = ("ok", "wrong") if not loose else ("ok",)
    order = [si for si, e in enumerate(spec) if not loose and e[2] in ("varname", "varlist")]
    order += [si for si in range(len(spec)) if si not in order]
    for state in passes:
        for si in order:
            name, mn, kind = spec[si]
            if si in taken:
                continue
            for oi, o in enumerate(opts):
                if oi in claim or not _abbrev(o.word, name, mn):
                    continue
                fits = _fits(o, kind)
                if state == "ok" and not fits and (loose or name in twice):
                    continue
                claim[oi] = (si, "ok" if fits else "wrong")
                taken.add(si)
                break
    for si, (name, mn, kind) in enumerate(spec):
        if si in taken or kind != "flag":
            continue
        for oi, o in enumerate(opts):
            if oi not in claim and o.word.startswith("no") and _abbrev(o.word[2:], name, mn):
                if loose and o.paren:
                    continue
                claim[oi] = (si, "wrong" if o.paren else "off")
                taken.add(si)
                break
    return claim


def loose_leftovers(opts, spec):
    """The options that an option list ending in * leaves in its options macro (every
    argument read as text): the words it does not know, those given twice, and those of the
    wrong form, in the order typed."""
    claim = _claims(opts, spec, loose=True)
    return [o for oi, o in enumerate(opts) if oi not in claim]


def loose_values(opts, spec):
    """The arguments as an option list ending in * takes them (text, as typed)."""
    claim = _claims(opts, spec, loose=True)
    out = {}
    for oi, (si, state) in claim.items():
        name, mn, kind = spec[si]
        o = opts[oi]
        if state == "off":
            out.setdefault(name, False)
        elif kind == "flag":
            out[name] = True
        elif o.arg is not None and o.arg.strip() != "":
            out[name + ("2" if kind == "passthru" and (name, mn, "flag") in spec else "")] = o.arg
    return out


def _vars_value(o, name, kind, columns):
    arg = o.arg.strip()
    if columns is None:
        return arg.split() if kind == "varlist" else arg
    try:
        vl = unab(arg, columns)
    except SeldecompError as ex:
        raise SeldecompError("%s\n(error in option %s())" % (ex.message, name), ex.rc)
    if kind == "varname":
        if len(vl) > 1:
            raise SeldecompError("%s(): too many variables specified" % name, 103)
        return vl[0]
    return vl


def run_syntax(opts, spec, columns, required=()):
    """The options as -syntax- checks them -> {name: value}. The varname and varlist options
    are checked first (a missing variable: 111), then each option in the order typed: a word
    not in spec or given twice is "not allowed", a flag with parentheses or an option without
    them and an integer option without an integer "incorrectly specified" (198). Values: True
    or False (a no-form) for a flag; the text in the parentheses (or the keyword's value) for
    the others, an int for integer, a variable name or list for varname and varlist; options
    with empty parentheses are left out."""
    claim = _claims(opts, spec, loose=False)
    values = {}
    for si, (name, mn, kind) in enumerate(spec):
        if kind not in ("varname", "varlist"):
            continue
        hit = [oi for oi, (s, st) in claim.items() if s == si and st == "ok"]
        if hit:
            o = opts[hit[0]]
            if o.arg is not None and o.arg.strip() != "":
                values[name] = _vars_value(o, name, kind, columns)
        if name in required and name not in values:
            raise SeldecompError("option %s() required" % name, 198)
    for oi, o in enumerate(opts):
        if oi not in claim:
            raise SeldecompError("option %s not allowed" % o.shown(), 198)
        si, state = claim[oi]
        name, mn, kind = spec[si]
        if state == "wrong":
            raise SeldecompError("option %s incorrectly specified" % o.shown(), 198)
        if kind in ("varname", "varlist"):
            continue
        if state == "off":
            values[name] = False
            continue
        if kind == "flag":
            values[name] = True
            continue
        arg = o.arg if o.arg is not None else ""
        if kind == "integer":
            f = stata_number(arg)
            if f != f or f != int(f):
                raise SeldecompError("option %s incorrectly specified" % o.shown(), 198)
            values[name] = int(f)
        elif kind == "cilevel":
            t = arg.strip()
            if t == "":
                continue
            f = stata_number(t)
            if f != f:
                raise SeldecompError("%s found where number expected" % t, 198)
            if not 10 <= f <= 99.99:
                raise SeldecompError("level() must be between 10 and 99.99 inclusive", 198)
            if abs(f * 100 - round(f * 100)) > 1e-7:
                raise SeldecompError("level() can have at most two digits after the decimal point", 198)
            values[name] = f
        else:
            if arg.strip() == "" and not (kind == "passthru" and arg != ""):
                continue
            key = name + "2" if kind == "passthru" and (name, mn, "flag") in spec else name
            values[key] = o.value if o.value is not None else arg
    return values


def spec_name(word, spec):
    """The name in spec that word (typed positively or as a no-form) stands for, or None."""
    for name, mn, kind in spec:
        if _abbrev(word, name, mn):
            return name
    if word.startswith("no"):
        for name, mn, kind in spec:
            if kind == "flag" and _abbrev(word[2:], name, mn):
                return name
    return None


def subopt_syntax(opt, text, spec, columns):
    """A suboption list (bayesboot(), contvars() after its comma) as Stata's Subopt_error
    reads it: a suboption the list does not know (or given twice, or of the wrong form) is
    "opt(): option x not allowed" (198, with a hint for seed() and strata() in bayesboot());
    any other refusal of -syntax- keeps its message and code, followed by
    "(error in option opt())"."""
    try:
        opts = option_tokens(text)
    except SeldecompError as ex:
        raise SeldecompError("%s\n(error in option %s())" % (ex.message, opt), ex.rc)
    try:
        return run_syntax(opts, spec, columns), opts
    except SeldecompError as ex:
        left = loose_leftovers(opts, spec)
        if left:
            nm = re.sub(r"\(.*$", "()", left[0].text(), flags=re.S)
            hint = ""
            if opt == "bayesboot" and nm in ("seed", "seed()"):
                hint = "; use rseed()"
            if opt == "bayesboot" and re.match(r"^str(a(t(a)?)?)?(\(\))?$", nm):
                hint = "; for a stratified design use svy jackknife: or svy brr:"
            raise SeldecompError("%s(): option %s not allowed%s" % (opt, nm, hint), 198)
        raise SeldecompError("%s\n(error in option %s())" % (ex.message, opt), ex.rc)


# ---- the command line ----------------------------------------------------------------------------
_WEIGHT_WORDS = (("aweight", 2, "aweight"), ("pweight", 2, "pweight"), ("iweight", 2, "iweight"),
                 ("fweight", 2, "fweight"), ("frequency", 4, "fweight"), ("cellsize", 4, "aweight"),
                 ("weight", 1, "aweight"))


def weight_type(word):
    """The weight type a weight word stands for ([pw=...], [pweigh=...], [w=...]), or None."""
    for full, mn, typ in _WEIGHT_WORDS:
        if len(word) >= mn and full.startswith(word):
            return typ
    return None


def _qualifier_at(s, i):
    """The qualifier that starts at s[i] ("if", "in", "weight") or None."""
    if s[i] == "[":
        m = re.match(r"\[\s*([A-Za-z]+)\s*=", s[i:])
        if m and weight_type(m.group(1)):
            return "weight"
        return None
    if (i == 0 or s[i - 1].isspace()) and s[i:i + 2] in ("if", "in"):
        after = s[i + 2:i + 3]
        if after == "" or after.isspace() or after == "(":
            return s[i:i + 2]
    return None


def _next_qualifier(s, start, kinds):
    """The position of the next qualifier of kinds outside parentheses and quotes, or len(s)."""
    depth, quote, i, n = 0, False, start, len(s)
    while i < n:
        c = s[i]
        if quote:
            quote = c != '"'
        elif c == '"':
            quote = True
        elif c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        elif depth == 0 and _qualifier_at(s, i) in kinds:
            return i
        i += 1
    return n


class Line:
    """A command line split as -syntax- splits it: zero (the line without the command name),
    anything (the variable list), the qualifiers in the order typed (("if", expression),
    ("in", range), ("weight", (type, expression, text as typed))), options (the text after
    the comma, or None), and error (a refusal of the line itself, raised after the
    qualifiers before it are checked)."""

    def __init__(self, zero):
        self.zero = zero
        self.anything = ""
        self.quals = []
        self.options = None
        self.error = None

    def qual(self, kind):
        for k, v in self.quals:
            if k == kind:
                return v
        return None

    def qual_text(self):
        """The qualifiers as -syntax- returns them, in its order: if, in, weight."""
        out = []
        if self.qual("if") is not None:
            out.append("if " + self.qual("if"))
        if self.qual("in") is not None:
            out.append("in " + self.qual("in"))
        w = self.qual("weight")
        if w is not None:
            out.append("[%s= %s]" % (w[0], w[1]))
        return " ".join(out)


def _depth(text):
    """The parentheses left open at the end of text (outside quotes)."""
    depth, quote = 0, False
    for c in text:
        if quote:
            quote = c != '"'
        elif c == '"':
            quote = True
        elif c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
    return depth


def split_command(cmd):
    """The command line -> Line."""
    zero = str(cmd).strip()
    if re.match(r"^seldecomp\b", zero):
        zero = zero[len("seldecomp"):].strip()
    line = Line(zero)
    k = find_top(zero, ",")
    if k < 0 and _depth(zero) > 0:
        # a parenthesis left open does not hide the comma of the options
        q = False
        for j, c in enumerate(zero):
            q = (not q) if c == '"' else q
            if c == "," and not q:
                k = j
                break
    left = zero if k < 0 else zero[:k]
    if k >= 0:
        line.options = zero[k + 1:]
    i = _next_qualifier(left, 0, ("if", "in", "weight"))
    line.anything = left[:i].strip()
    if _depth(line.anything) > 0:
        line.error = SeldecompError("invalid varlist: unmatched open parenthesis or bracket", 198)
    n = len(left)
    while i < n and line.error is None:
        if left[i].isspace():
            i += 1
            continue
        kind = _qualifier_at(left, i)
        if kind == "if" and line.qual("if") is None:
            j = _next_qualifier(left, i + 2, ("in", "weight"))
            line.quals.append(("if", left[i + 2:j].strip()))
            i = j
        elif kind == "in" and line.qual("in") is None:
            m = re.match(r"in\s*([^\s/\[]*)\s*(/\s*([^\s/\[]*))?", left[i:])
            line.quals.append(("in", " ".join(m.group(0)[2:].split())))
            i += m.end()
        elif kind == "weight":
            j = left.find("]", i)
            j = n if j < 0 else j
            inner = left[i + 1:j]
            word, exp = inner.split("=", 1)
            if line.qual("weight") is not None:
                line.error = SeldecompError("too many weights", 198)
                break
            line.quals.append(("weight", (weight_type(word.strip()), exp.strip(), left[i:j + 1])))
            i = j + 1
        elif kind in ("if", "in"):
            line.error = SeldecompError("invalid syntax", 198)
        else:
            w = re.match(r"\w+|\S", left[i:]).group(0)
            line.error = SeldecompError("invalid '%s'" % w, 198)
    return line


def check_qualifiers(line, df):
    """The qualifiers checked in the order typed, as -syntax- checks them -> (rows selected by
    if and in, the weights or None)."""
    n = len(df)
    mask = np.ones(n, dtype=bool)
    w = None
    for kind, val in line.quals:
        if kind == "if":
            mask &= if_mask(val, df)
        elif kind == "in":
            mask &= in_range(val, n)
        else:
            w = weight_values(val[1], df)
    if line.error is not None:
        raise line.error
    return mask, w


def looks_like_command(spec):
    """Whether the text is a full Stata command line (it has options after a comma, a weight
    in brackets, or an if or in qualifier) and not just "y rakespec"."""
    if find_top(spec, ",") >= 0:
        return True
    if re.search(r"\[\s*[A-Za-z]+\s*=", spec):
        return True
    toks = split_top(spec)
    return "if" in toks or "in" in toks


def parse_command(cmd):
    """A Stata command line -> dict(varlist, if_, in_, wtype, wexp, options), options as
    {full option name: True, False (a no-form) or the argument text}."""
    line = split_command(cmd)
    w = line.qual("weight")
    options = {}
    if line.options is not None:
        for o in option_tokens(line.options):
            if re.match(r"^(bayes|bayesb|bayesbo|bayesboo|bayesboot|make|makep|makepl|makeplo|makeplot)$",
                        o.word) and o.paren and o.arg == "":
                o = Opt(o.word)
            name = spec_name(o.word, MAIN_SPEC)
            if name is None:
                continue
            if o.paren:
                options[name] = o.arg.strip()
            elif isinstance(options.get(name), str):
                continue
            else:
                options[name] = not (o.word.startswith("no")
                                     and not _abbrev(o.word, name, MAIN_OPTIONS[name]))
    return {"varlist": line.anything, "if_": line.qual("if"), "in_": line.qual("in"),
            "wtype": w[0] if w else None, "wexp": w[1] if w else None, "options": options}


def in_range(text, n):
    """A Stata in range ("101/1500", "f/100", "-10/l", "5") -> boolean mask over n rows, with
    Stata's messages for a range it refuses."""
    t = str(text).strip()
    parts = [p.strip() for p in t.split("/")]
    if len(parts) > 2:
        raise SeldecompError("'/' invalid observation number", 198)

    def pos(s):
        if s in ("f", "F"):
            return 1
        if s in ("l", "L"):
            return n
        if re.match(r"^[-+]?\d+$", s):
            v = int(s)
            if v == 0 or v < -n:
                raise SeldecompError("'%s' invalid observation number" % s, 198)
            return n + v + 1 if v < 0 else v
        raise SeldecompError("'%s' invalid observation number" % (s if s else ","), 198)

    a = pos(parts[0])
    b = pos(parts[1]) if len(parts) > 1 else a
    if a > n or b > n or a > b:
        raise SeldecompError("observation numbers out of range", 198)
    mask = np.zeros(n, dtype=bool)
    mask[a - 1:b] = True
    return mask


# ---- Stata expressions -------------------------------------------------------------------------
class ExprUnsupported(Exception):
    """An expression (a function, a subscript) that the translator does not read."""


_EXPR_TOKEN = re.compile(
    r"\s*(?:(?P<num>(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)"
    r"|(?P<miss>\.[a-z]?(?![\w.]))"
    r"|(?P<str>`\"[^\"]*\"'|\"[^\"]*\")"
    r"|(?P<name>[^\W\d]\w*)"
    r"|(?P<op>==|!=|~=|>=|<=|[-+*/^<>&|!~(),=\[\]]))")


def _expr_tokens(text):
    toks, i, n = [], 0, len(text)
    while i < n:
        if text[i:].strip() == "":
            break
        m = _EXPR_TOKEN.match(text, i)
        if not m:
            raise SeldecompError("invalid syntax", 198)
        kind = m.lastgroup
        toks.append((kind, m.group(kind)))
        i = m.end()
    return toks


def _column_value(df, name):
    """A column as an expression reads it: ("n", float64 with nan for missing) or ("s", str)."""
    s = df[name]
    if isinstance(s.dtype, pd.CategoricalDtype):
        if pd.api.types.is_numeric_dtype(s.cat.categories):
            return "n", s.astype("float64").to_numpy()
        s = s.astype(object)
    if pd.api.types.is_bool_dtype(s) or pd.api.types.is_numeric_dtype(s):
        try:
            return "n", s.astype("float64").to_numpy()
        except (TypeError, ValueError):
            return "n", pd.to_numeric(s, errors="coerce").astype("float64").to_numpy()
    vals = [x if isinstance(x, str) else ("" if x is None or x is pd.NA or (isinstance(x, float) and x != x)
                                          else str(x)) for x in s]
    return "s", np.array(vals, dtype=object)


def _clean(x):
    x = np.asarray(x, dtype=float)
    return np.where(np.isfinite(x), x, np.nan)


class _Eval:
    """A recursive-descent reading of a Stata expression, evaluated on the rows of a frame in
    Stata's order of operations (!, ^, negation, * and /, + and -, the comparisons, &, |).
    A missing number is larger than every number in a comparison; arithmetic with it is
    missing; a 4-byte float column is compared in double precision, as Stata compares it."""

    def __init__(self, toks, df):
        self.t = toks
        self.k = 0
        self.df = df
        self.n = len(df)
        self.cache = {}

    def peek(self):
        return self.t[self.k] if self.k < len(self.t) else (None, None)

    def operand_follows(self, op):
        if self.k >= len(self.t):
            prev = self.t[self.k - 2][1] if self.k >= 2 else ""
            raise SeldecompError("%s%s invalid name" % (prev, op), 198)

    def full(self, x):
        return np.full(self.n, x, dtype=object if isinstance(x, str) else float)

    def num(self, v):
        if v[0] != "n":
            raise SeldecompError("type mismatch", 109)
        return v[1]

    def expr(self):
        a = self.and_()
        while self.peek() == ("op", "|"):
            self.k += 1
            self.operand_follows("|")
            b = self.and_()
            a = ("n", ((self.num(a) != 0) | (self.num(b) != 0)).astype(float))
        return a

    def and_(self):
        a = self.cmp()
        while self.peek() == ("op", "&"):
            self.k += 1
            self.operand_follows("&")
            b = self.cmp()
            a = ("n", ((self.num(a) != 0) & (self.num(b) != 0)).astype(float))
        return a

    def cmp(self):
        a = self.add()
        while self.peek()[0] == "op" and self.peek()[1] in ("==", "!=", "~=", ">", "<", ">=", "<="):
            op = self.t[self.k][1]
            self.k += 1
            self.operand_follows(op)
            b = self.add()
            a = ("n", _compare(op, a, b).astype(float))
        return a

    def add(self):
        a = self.mul()
        while self.peek()[0] == "op" and self.peek()[1] in ("+", "-"):
            op = self.t[self.k][1]
            self.k += 1
            self.operand_follows(op)
            b = self.mul()
            if op == "+" and a[0] == "s" and b[0] == "s":
                a = ("s", np.array([x + y for x, y in zip(a[1], b[1])], dtype=object))
            elif op == "+":
                a = ("n", _clean(self.num(a) + self.num(b)))
            else:
                a = ("n", _clean(self.num(a) - self.num(b)))
        return a

    def mul(self):
        a = self.neg()
        while self.peek()[0] == "op" and self.peek()[1] in ("*", "/"):
            op = self.t[self.k][1]
            self.k += 1
            self.operand_follows(op)
            b = self.neg()
            with np.errstate(all="ignore"):
                x = self.num(a) * self.num(b) if op == "*" else self.num(a) / self.num(b)
            a = ("n", _clean(x))
        return a

    def neg(self):
        if self.peek() == ("op", "-"):
            self.k += 1
            self.operand_follows("-")
            return "n", -self.num(self.neg())
        if self.peek() == ("op", "+"):
            self.k += 1
            self.operand_follows("+")
            return "n", self.num(self.neg())
        return self.pow()

    def pow(self):
        a = self.not_()
        while self.peek() == ("op", "^"):
            self.k += 1
            self.operand_follows("^")
            b = self.neg()
            with np.errstate(all="ignore"):
                a = ("n", _clean(np.power(self.num(a), self.num(b))))
        return a

    def not_(self):
        if self.peek()[0] == "op" and self.peek()[1] in ("!", "~"):
            op = self.t[self.k][1]
            self.k += 1
            self.operand_follows(op)
            return "n", (self.num(self.not_()) == 0).astype(float)
        return self.atom()

    def atom(self):
        kind, text = self.peek()
        if kind is None:
            raise SeldecompError("invalid syntax", 198)
        if (kind, text) in (("op", "-"), ("op", "+")):
            return self.neg()
        self.k += 1
        if kind == "num":
            return "n", self.full(float(text))
        if kind == "miss":
            return "n", self.full(np.nan)
        if kind == "str":
            return "s", self.full(text[2:-2] if text.startswith("`") else text[1:-1])
        if kind == "name":
            if self.peek() == ("op", "("):
                self.k += 1
                args = []
                if self.peek() != ("op", ")"):
                    while True:
                        args.append(self.expr())
                        if self.peek() == ("op", ","):
                            self.k += 1
                            continue
                        break
                if self.peek() != ("op", ")"):
                    raise SeldecompError("too few ')' or ']'", 132)
                self.k += 1
                return _function(text, args, self)
            if self.peek() == ("op", "["):
                raise ExprUnsupported()
            if text == "_n":
                return "n", np.arange(1, self.n + 1, dtype=float)
            if text == "_N":
                return "n", self.full(float(self.n))
            if text == "_pi":
                return "n", self.full(math.pi)
            if text in self.df.columns:
                if text not in self.cache:
                    self.cache[text] = _column_value(self.df, text)
                return self.cache[text]
            raise SeldecompError("%s not found" % text, 111)
        if text == "(":
            v = self.expr()
            if self.peek() != ("op", ")"):
                if self.peek()[0] is None:
                    raise SeldecompError("too few ')' or ']'", 132)
                raise SeldecompError("invalid syntax", 198)
            self.k += 1
            return v
        raise SeldecompError("invalid syntax", 198)


def _compare(op, a, b):
    if a[0] != b[0]:
        raise SeldecompError("type mismatch", 109)
    x, y = a[1], b[1]
    if a[0] == "n":
        x = np.where(np.isnan(x), np.inf, x)
        y = np.where(np.isnan(y), np.inf, y)
    if op == "==":
        r = x == y
    elif op in ("!=", "~="):
        r = x != y
    elif op == ">":
        r = x > y
    elif op == "<":
        r = x < y
    elif op == ">=":
        r = x >= y
    else:
        r = x <= y
    return np.asarray(r, dtype=bool)


def _missing_of(v):
    return v[1] == "" if v[0] == "s" else np.isnan(v[1])


def _function(name, args, ev):
    """The functions the translator reads."""
    nargs = len(args)

    def need(lo, hi):
        if not lo <= nargs <= hi:
            raise SeldecompError("invalid syntax", 198)

    def strs(i):
        if args[i][0] != "s":
            raise SeldecompError("type mismatch", 109)
        return args[i][1]

    if name in ("missing", "mi"):
        need(1, 10 ** 6)
        r = np.zeros(ev.n, dtype=bool)
        for a in args:
            r |= _missing_of(a)
        return "n", r.astype(float)
    if name == "inlist":
        need(2, 10 ** 6)
        r = np.zeros(ev.n, dtype=bool)
        for a in args[1:]:
            r |= _compare("==", args[0], a)
        return "n", r.astype(float)
    if name == "inrange":
        need(3, 3)
        z, a, b = args
        if not z[0] == a[0] == b[0]:
            raise SeldecompError("type mismatch", 109)
        if z[0] == "s":
            return "n", ((a[1] <= z[1]) & (z[1] <= b[1])).astype(float)
        zv, av, bv = z[1], a[1], b[1]
        lo = np.where(np.isnan(av), -np.inf, av)
        hi = np.where(np.isnan(bv), np.inf, bv)
        r = (zv >= lo) & (zv <= hi)
        r = np.where(np.isnan(zv), np.isnan(av) & np.isnan(bv), r)
        return "n", r.astype(float)
    one = {"abs": np.abs, "sqrt": np.sqrt, "exp": np.exp, "ln": np.log, "log": np.log,
           "log10": np.log10, "int": np.trunc, "trunc": np.trunc, "floor": np.floor, "ceil": np.ceil}
    if name in one:
        need(1, 1)
        with np.errstate(all="ignore"):
            return "n", _clean(one[name](ev.num(args[0])))
    if name == "float":
        need(1, 1)
        return "n", ev.num(args[0]).astype(np.float32).astype(float)
    if name == "round":
        need(1, 2)
        x = ev.num(args[0])
        y = ev.num(args[1]) if nargs == 2 else np.ones(ev.n)
        with np.errstate(all="ignore"):
            return "n", _clean(y * np.floor(x / y + 0.5))
    if name == "mod":
        need(2, 2)
        x, y = ev.num(args[0]), ev.num(args[1])
        with np.errstate(all="ignore"):
            return "n", _clean(x - y * np.floor(x / y))
    if name in ("min", "max"):
        need(1, 10 ** 6)
        m = np.vstack([ev.num(a) for a in args])
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = np.nanmin(m, axis=0) if name == "min" else np.nanmax(m, axis=0)
        return "n", r
    if name == "cond":
        need(3, 4)
        c = ev.num(args[0])
        if args[1][0] != args[2][0]:
            raise SeldecompError("type mismatch", 109)
        r = np.where(c != 0, args[1][1], args[2][1])
        if nargs == 4:
            r = np.where(np.isnan(c), args[3][1], r)
        return args[1][0], r
    if name == "real":
        need(1, 1)
        return "n", np.array([stata_number(x) for x in strs(0)], dtype=float)
    if name in ("strlen", "length", "ustrlen"):
        need(1, 1)
        return "n", np.array([len(x) for x in strs(0)], dtype=float)
    if name in ("strpos", "ustrpos"):
        need(2, 2)
        return "n", np.array([x.find(y) + 1 for x, y in zip(strs(0), strs(1))], dtype=float)
    if name in ("substr", "usubstr"):
        need(3, 3)
        s, p, m = strs(0), ev.num(args[1]), ev.num(args[2])
        out = []
        for x, a, b in zip(s, p, m):
            if a != a or a == 0:
                out.append("")
                continue
            a = int(a)
            start = a - 1 if a > 0 else len(x) + a
            if start < 0 or start >= len(x):
                out.append("")
                continue
            out.append(x[start:] if b != b else x[start:start + max(int(b), 0)])
        return "s", np.array(out, dtype=object)
    strfun = {"lower": str.lower, "strlower": str.lower, "upper": str.upper, "strupper": str.upper,
              "trim": str.strip, "strtrim": str.strip, "ltrim": str.lstrip, "strltrim": str.lstrip,
              "rtrim": str.rstrip, "strrtrim": str.rstrip}
    if name in strfun:
        need(1, 1)
        f = strfun[name]
        return "s", np.array([f(x) if "trim" not in name else f(x, " ") for x in strs(0)], dtype=object)
    raise ExprUnsupported()


_RESERVED = ("if", "in", "using", "with", "_all", "_b", "_coef", "_cons", "byte", "double", "float", "int",
             "long", "_pred", "_rc", "_skip", "strL")


def eval_expr(text, df, words=False):
    """A Stata expression evaluated on the rows of df -> ("n", float array) or ("s", str array).
    Raises SeldecompError with Stata's message (a variable not found: 111; invalid syntax: 198;
    type mismatch: 109) or ExprUnsupported. words: a word left over after the expression is
    "invalid 'word'" (an if qualifier), else invalid syntax."""
    toks = _expr_tokens(text)
    if not toks:
        raise SeldecompError("invalid syntax", 198)
    ev = _Eval(toks, df)
    v = ev.expr()
    if ev.k < len(toks):
        if toks[ev.k] == ("op", ")"):
            raise SeldecompError("too many ')' or ']'", 132)
        if words and toks[ev.k][0] in ("name", "num") and toks[ev.k][1] not in _RESERVED:
            raise SeldecompError("invalid '%s'" % toks[ev.k][1], 198)
        raise SeldecompError("invalid syntax", 198)
    return v


def if_mask(expr, df):
    """A Stata if expression -> boolean mask: comparison, logical and arithmetic operators,
    missing(), inlist(), inrange() and the common math and string functions are read as Stata
    reads them (anything else should be passed as subset=). As in Stata, a missing number is
    larger than every number (x > 5 is true where x is missing; x < . selects the nonmissing
    values), and an expression that is missing counts as true."""
    try:
        kind, v = eval_expr(expr, df, words=True)
    except ExprUnsupported:
        raise SeldecompError("if %s: the expression could not be evaluated; pass a boolean mask "
                             "as subset= instead" % expr, 198)
    if kind != "n":
        raise SeldecompError("type mismatch", 109)
    return np.asarray(v != 0, dtype=bool)


def weight_values(expr, df):
    """The weight expression (or a column's name) -> float array (nan where missing)."""
    if str(expr).strip() == "":
        raise SeldecompError("invalid syntax", 198)
    if str(expr).strip() in df.columns:
        kind, v = _column_value(df, str(expr).strip())
        if kind != "n":
            raise SeldecompError("type mismatch", 109)
        return np.asarray(v, dtype=float)
    try:
        kind, v = eval_expr(expr, df)
    except ExprUnsupported:
        raise SeldecompError("[weight=%s]: the expression could not be evaluated; pass the weights as "
                             "a column instead" % expr, 198)
    if kind != "n":
        raise SeldecompError("type mismatch", 109)
    return np.asarray(v, dtype=float)


# ---- the command's own checks of the line ----------------------------------------------------------
_BARE = [(n, m) for n, m, k in MAIN_SPEC if k == "flag"]


def _isopt(w):
    """Whether the word is one of the command's options without argument (seldecomp_isopt)."""
    for name, mn in _BARE:
        if mn <= len(w) <= len(name) and name[:len(w)] == w:
            return True
    return False


def _commaby(s):
    """(text before by, by() and the rest) when the text before the first comma holds the word
    by( outside parentheses and quotes (seldecomp_commaby), else None."""
    n = len(s)
    d, q, c, b = 0, False, n + 1, 0
    for i in range(1, n + 1):
        ch = s[i - 1]
        if ch == '"':
            q = not q
        elif not q and ch == "(":
            if d == 0 and not b and i >= 3 and s[i - 3:i - 1] == "by":
                if i == 3:
                    b = 1
                elif s[i - 4] == " ":
                    b = i - 2
            d += 1
        elif not q and ch == ")":
            d -= 1
        elif not q and d == 0 and ch == ",":
            c = i
            break
    if not b:
        return None
    pre = s[:c - 1]
    post = pre[b - 1:] + (" " + s[c:] if c <= n else "")
    return pre[:b - 1], post


def _squeeze(s):
    return " ".join(s.split())


def _strip_seldby(s):
    return _squeeze(re.sub(r"seldby(var|idx|vars)[(][^)]*[)]", "", s))


def comma_check(line, ok, columns):
    """Comma_check: an option typed before the comma (by(married) among the raking variables,
    or by() after an if or in qualifier), or an option word that is no variable typed there
    while options follow the comma, is refused (198) with the command written as it should
    be. ok: whether the line parses (qualifiers and options)."""
    first, moved, pre, post = None, [], "", ""
    options = line.options or ""
    if ok:
        keep = []
        for k, tok in enumerate(bind_tokens(line.anything), 1):
            m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)[(]", tok)
            nm = m.group(1) if m else ""
            if re.match(r"^[A-Za-z_][A-Za-z0-9_]*[(][^()]*[)][.]", tok):
                nm = ""
            if nm in ("seldbyvar", "seldbyidx", "seldbyvars"):
                continue
            mv = False
            if nm:
                first = first or nm
                mv = True
            elif k > 1 and _NAME.match(tok) and tok not in columns:
                mv = first is not None or _isopt(tok)
            (moved if mv else keep).append(tok)
        pre = " ".join(keep) + " " + line.qual_text()
        post = " ".join(moved) + " " + options
    also = ""
    if first is None:
        cb = _commaby(line.zero)
        if cb is None:
            opts = _strip_seldby(options)
            if not moved or opts == "":
                return
            post = _strip_seldby(options + " " + " ".join(moved))
            raise SeldecompError("%s must follow the comma:\n    seldecomp %s, %s"
                                 % (and_list(moved), _squeeze(pre), post), 198)
        first = "by"
        pre, post = cb
        sub = split_command(cb[0])
        if sub.error is None and sub.options is None:
            keep, moved = [], []
            for k, tok in enumerate(bind_tokens(sub.anything), 1):
                mv = k > 1 and _NAME.match(tok) is not None and tok not in columns and _isopt(tok)
                (moved if mv else keep).append(tok)
            if moved:
                pre = " ".join(keep) + " " + sub.qual_text()
                post = cb[1] + " " + " ".join(moved)
                also = " (and %s)" % and_list(moved)
    raise SeldecompError("option %s() must follow the comma%s:\n    seldecomp %s, %s"
                         % (first, also, _squeeze(pre), _strip_seldby(post)), 198)


def reject_nodiagnostics(opts):
    """-syntax- would take nodiagnostics as diagnostics switched off: refused (198)."""
    for o in opts:
        w = o.text()
        if len(w) >= 6 and "nodiagnostics"[:len(w)] == w:
            raise SeldecompError("option nodiagnostics not allowed", 198)


_HABITS = (("bb", "reps", "r", "#"), ("bb", "rseed", "rseed", "#"), ("bb", "seed", "seed", "#"),
           ("bb", "cluster", "cl", "varname"), ("bb", "saving", "sav", "filename"),
           ("bb", "dots", "dots", "#"), ("bb", "nodots", "nodots", ""), ("bb", "bca", "bca", ""),
           ("bb", "ties", "ties", ""), ("bb", "level", "l", "#"), ("bb", "citable", "cit", ""),
           ("bb", "normal", "norm", ""), ("bb", "tinterval", "tint", ""),
           ("bb", "percentile", "perc", ""), ("bb", "bc", "bc", ""), ("oax", "swap", "swap", ""),
           ("oax", "weight", "w", ""), ("oax", "pooled", "pool", ""), ("oax", "vce", "vce", ""))


def _habit(w):
    """(class, full name, kind of argument) of a bayesboot() suboption or an oaxaca option
    the word abbreviates, or (None, None, None)."""
    for c, f, m, a in _HABITS:
        if len(w) >= len(m) and f[:len(w)] == w:
            return c, f, a
    return None, None, None


def _word_name(tok):
    return re.sub(r"[(].*$", "", tok, flags=re.S).lower()


def top_error(err, opts, df, sample, value_labels=None):
    """Top_error: -syntax- refused the options (err). When the first option it does not know
    is recut or replace typed outside contvars() or generate(), a suboption of bayesboot(),
    or an option of oaxaca (swap, weight(), pooled, vce()), the message names seldecomp's
    equivalent (198); otherwise err. sample: the rows of if and in."""
    left = loose_leftovers(opts, MAIN_SPEC)
    if not left:
        raise err
    vals = loose_values(opts, MAIN_SPEC)
    vals["_labels"] = value_labels or {}
    w1 = left[0].text()
    for word, opt in (("recut", "contvars"), ("replace", "generate")):
        if w1.lower() != word:
            continue
        cur = vals.get(opt)
        if cur is None:
            hold = "varname #" if opt == "contvars" else "newvar"
            raise SeldecompError("%s is a suboption of %s(): type %s(%s, %s)"
                                 % (word, opt, opt, hold, word), 198)
        if "," in cur:
            raise SeldecompError("%s is a suboption of %s(): type %s(%s %s)"
                                 % (word, opt, opt, cur, word), 198)
        raise SeldecompError("%s is a suboption of %s(): type %s(%s, %s)" % (word, opt, opt, cur, word), 198)
    cls, oword, harg = _habit(_word_name(w1))
    if cls is None:
        raise err
    if cls == "oax":
        raise SeldecompError(_oaxaca_message(left[0], oword, vals, df, sample), 198)
    bb_in = vals.get("bayesboot2", "")
    gname, gtok = [], {}
    for t in bind_tokens(bb_in):
        c, f, a = _habit(_word_name(t))
        if c == "bb":
            f = "rseed" if f == "seed" else f
            if f not in gname:
                gname.append(f)
                gtok[f] = t
    names, toks, dtoks, dins = [], [], [], []
    seed = False
    for o in left:
        tok = o.text()
        c, f, a = _habit(_word_name(tok))
        if c != "bb":
            continue
        ttok = tok
        par = "()" if "(" in tok else ""
        if f == "seed":
            seed = True
            f = "rseed"
            tok = "r" + tok
        if par == "" and a:
            par = "()"
            tok = "dots(1)" if f == "dots" else "%s(%s)" % (f, a)
        hn = f + par
        if hn not in names:
            names.append(hn)
        if f in gname:
            dtoks.append(ttok)
            if gtok[f] not in dins:
                dins.append(gtok[f])
            continue
        gname.append(f)
        gtok[f] = tok
        toks.append(tok)
    sug = _squeeze(bb_in + " " + " ".join(toks))
    sh = "; use rseed(), not seed()" if seed else ""
    if "level()" in names:
        sh += "; with a replication prefix, give level() to the prefix"
    dh = ""
    if dtoks:
        dh = "; %s %s left out, as bayesboot() has %s already" % (
            and_list(dtoks), "is" if len(dtoks) == 1 else "are", and_list(dins))
    raise SeldecompError("%s %s of bayesboot(): type bayesboot(%s)%s%s"
                         % (and_list(names), "is a suboption" if len(names) == 1 else "are suboptions",
                            sug, dh, sh), 198)


def _oaxaca_message(o, oword, vals, df, sample):
    """The message for an option of oaxaca: the other group named in focalgroup()."""
    from ._stata import fmt_g, reldif
    other = "<the other group>"
    fg = vals.get("focalgroup")
    byv = vals.get("by")
    oarg = o.arg.strip() if o.paren and o.arg is not None else ""
    if byv is not None and df is not None:
        try:
            names = unab(byv, df.columns)
        except SeldecompError:
            names = []
        if len(names) == 1:
            kind, col = _column_value(df, names[0])
            if kind == "n":
                x = col[sample]
                x = x[~np.isnan(x)]
                if len(x) and x.min() < x.max():
                    vmin, vmax = float(x.min()), float(x.max())
                    kf = 0
                    if fg is not None:
                        f = stata_number(fg)
                        if f == f:
                            d1, d2 = reldif(f, vmin) < 1e-7, reldif(f, vmax) < 1e-7
                            kf = 1 if f == vmin else (2 if f == vmax else (1 if d1 and not d2 else
                                                                              (2 if d2 and not d1 else 0)))
                        labels = vals["_labels"].get(names[0]) or {}
                        if kf == 0 and labels:
                            m = [v == int(v) and labels.get(int(v)) == fg for v in (vmin, vmax)]
                            kf = 1 if m[0] and not m[1] else (2 if m[1] and not m[0] else 0)
                    if kf == 1 or fg is None:
                        other = fmt_g(vmax, 12)
                    elif kf == 2:
                        other = fmt_g(vmin, 12)
            elif fg is not None:
                u = sorted(set(x for x in col[sample] if x != ""))
                if len(u) == 2:
                    if fg == u[0]:
                        other = u[1]
                    elif fg == u[1]:
                        other = u[0]
    sw = "focalgroup(%s)" % other
    if fg is not None:
        sw = "focalgroup(%s) in place of focalgroup(%s)" % (other, fg)
    if oword == "swap":
        return "swap not allowed; to rake the other group, use %s" % sw
    if oword == "weight" and oarg in ("0", "0.0"):
        return "weight(0) not allowed; to rake the other group, use %s" % sw
    if oword == "weight" and oarg in ("1", "1.0"):
        return "weight(1) not allowed; it is the default"
    if oword == "weight":
        return "weight() not allowed; the default corresponds to weight(1), and %s to weight(0)" % sw
    if oword == "pooled":
        return "pooled not allowed; pooled is a suboption of contvars()"
    vt = re.sub(r"[ ,].*$", "", oarg).lower()
    m = re.match(r"^[^ ,]+[ ]+([^ ,]+)", oarg)
    cv = m.group(1) if m else ""
    if vt and "bootstrap".startswith(vt):
        return "vce(bootstrap) not allowed; use bayesboot() or the bootstrap: prefix"
    if vt and "jackknife".startswith(vt):
        return "vce(jackknife) not allowed; use the jackknife: prefix"
    if vt and "cluster".startswith(vt) and cv:
        return "vce(cluster) not allowed; use bayesboot(cluster(%s))" % cv
    return "vce() not allowed; use bayesboot() or a replication prefix"


# ---- formats -------------------------------------------------------------------------------------
def confirm_numeric_format(fmt):
    """Stata's confirm numeric format: a numeric display format (%w.d with e, f or g, and c for
    f and g; the hexadecimal, binary and date formats), else Stata's message (7)."""
    toks = str(fmt).split()
    if len(toks) > 1:
        raise SeldecompError("'%s' found where nothing expected" % toks[1][0], 7)
    f = toks[0] if toks else ""
    if not (numeric_format(f) or re.match(r"^%-?[1-9]\d*[xHL]$", f) or re.match(r"^%-?t[cCdwmqhybg]", f)
            or re.match(r"^%-?d", f)):
        raise SeldecompError("'%s' found where numeric format expected" % f, 7)


def other_format_width(fmt):
    """The width of a hexadecimal or binary format (%21x, %16H, %8L); None for a date format."""
    m = re.match(r"^%-?([1-9]\d*)[xHL]$", str(fmt).strip())
    return int(m.group(1)) if m else None


FORMAT_WIDTHS = {"cformat": 9, "pformat": 5, "sformat": 8}


def check_display_format(name, fmt):
    """A display format option: Stata's test first; a format Stata reads and this package does
    not show (hexadecimal, binary, date) is refused, unless it is wider than the column, when
    the default is used as for any format too wide."""
    confirm_numeric_format(fmt)
    f = str(fmt).strip()
    if not numeric_format(f):
        w = other_format_width(f)
        if w is None or w <= FORMAT_WIDTHS[name]:
            raise SeldecompError("%s(%s): date, hexadecimal and binary formats are not shown; give a %%w.d "
                                 "format with e, f or g" % (name, f), 198)
    return f


# ---- the list of raking variables ----------------------------------------------------------------
class Term:
    def __init__(self, vars_, kinds, name):
        self.vars = vars_
        self.kinds = kinds
        self.name = name


class Block:
    def __init__(self, name, terms):
        self.name = name
        self.terms = terms


def _valid_name(s):
    return valid_name(s)


def parse_term(tok, columns):
    """One token of the raking list -> a list of terms (one per variable, or one interaction)."""
    if "(" in tok or ")" in tok or "," in tok:
        raise SeldecompError("invalid term %s in the raking variable list: i.(varlist), nested "
                             "parentheses, and commas are not allowed" % tok, 198)
    t = tok.replace("##", "#")
    if t.startswith("#") or t.endswith("#") or "##" in t:
        raise SeldecompError("invalid interaction %s in the raking variable list" % tok, 198)
    inter = "#" in t
    vars_, kinds = [], []
    for piece in t.split("#"):
        if piece.startswith("c."):
            kind, v = "c", piece[2:]
        elif piece.startswith("i."):
            kind, v = "i", piece[2:]
        elif "." in piece:
            raise SeldecompError("%s: only the i. and c. operators are allowed in the raking "
                                 "variable list" % piece, 198)
        else:
            kind, v = "i", piece
        if v == "":
            raise SeldecompError("invalid term %s in the raking variable list" % tok, 198)
        try:
            vl = unab(v, columns)
        except SeldecompError:
            raise SeldecompError("variable %s not found" % v, 111)
        if inter and len(vl) > 1:
            raise SeldecompError("%s: wildcards and variable lists are not allowed inside "
                                 "interaction terms" % tok, 198)
        for vv in vl:
            vars_.append(vv)
            kinds.append(kind)
    if not vars_:
        raise SeldecompError("invalid term %s in the raking variable list" % tok, 198)
    if inter:
        if len(set(vars_)) < len(vars_):
            raise SeldecompError("%s: a variable may appear only once in an interaction" % tok, 198)
        return [Term(vars_, kinds, "#".join(vars_))]
    return [Term([v], [k], v) for v, k in zip(vars_, kinds)]


def rake_tokens(text):
    """The words of the raking list as gettoken with match(paren) reads them: a word that
    begins with a parenthesis runs to the parenthesis that closes it."""
    out, i, n = [], 0, len(text)
    while i < n:
        if text[i].isspace():
            i += 1
            continue
        j = i
        if text[i] == "(":
            depth = 0
            while j < n:
                if text[j] == "(":
                    depth += 1
                elif text[j] == ")":
                    depth -= 1
                    if depth == 0:
                        j += 1
                        break
                j += 1
        else:
            while j < n and not text[j].isspace():
                j += 1
        out.append(text[i:j])
        i = j
    return out


def parse_rakespec(text, columns):
    """The list of raking variables -> (blocks, terms)."""
    m = re.search(r"[A-Za-z_][A-Za-z0-9_]*[(][^()]*[)][.][A-Za-z_][A-Za-z0-9_]*", text)
    if m:
        raise SeldecompError("%s: only the i. and c. operators are allowed in the raking "
                             "variable list" % m.group(0), 198)
    blocks, terms = [], []
    for tok in rake_tokens(text):
        if tok.startswith("("):
            inner = tok[1:-1] if tok.endswith(")") else tok[1:]
            bname = ""
            if ":" in inner:
                c = inner.index(":")
                bname, inner = inner[:c], inner[c + 1:]
                if bname.strip() == "":
                    bname, inner = ":", inner
                bname = bname.strip()
                if len(bname.split()) != 1:
                    raise SeldecompError("invalid block name in (%s: ...)" % bname, 198)
                if not _valid_name(bname):
                    raise SeldecompError("invalid block name %s; block names must be valid Stata "
                                         "names" % bname, 198)
            inner = inner.strip()
            if inner == "":
                raise SeldecompError("empty block () in the raking variable list", 198)
            ids = []
            for t1 in inner.split():
                for term in parse_term(t1, columns):
                    terms.append(term)
                    ids.append(len(terms) - 1)
            if bname == "":
                bname = terms[ids[0]].name.replace("#", "_")
            blocks.append(Block(bname, ids))
        else:
            for term in parse_term(tok, columns):
                terms.append(term)
                blocks.append(Block(term.name.replace("#", "_"), [len(terms) - 1]))
    for i, t in enumerate(terms):
        for t2 in terms[:i]:
            if t.name == t2.name or sorted(t.vars) == sorted(t2.vars):
                raise SeldecompError("raking variable %s specified more than once" % t.name, 198)
    if not terms:
        raise SeldecompError("at least one raking variable is required", 102)
    return blocks, terms


# ---- contvars() ----------------------------------------------------------------------------------
def parse_contvars(arg, df=None):
    """contvars: "x1 4 x2 5, recut focalgroup" or {"x1": 4, ...} (with optional keys recut,
    sample) -> (list of (variable, bins), sample, recut, text), checked in Stata's order."""
    sample, recut = "pooled", False
    columns = None if df is None else list(df.columns)
    if isinstance(arg, dict):
        d = dict(arg)
        recut = bool(d.pop("recut", False))
        chosen = [k for k in ("pooled", "focalgroup", "refgroup") if d.pop(k, False)]
        if "sample" in d:
            chosen.append(d.pop("sample"))
        if len(chosen) > 1:
            raise SeldecompError("contvars(): only one of pooled, focalgroup, and refgroup is allowed", 198)
        if chosen:
            sample = chosen[0]
        words = []
        for k, v in d.items():
            words += [str(k), str(v)]
        text = " ".join(words)
    else:
        text = str(arg).strip()
        k = find_top(text, ",")
        left, right = (text, "") if k < 0 else (text[:k], text[k + 1:])
        if left.strip() == "":
            raise SeldecompError("contvars() required\n(error in option contvars())", 100)
        vals, toks = subopt_syntax("contvars", right, CONTVARS_SPEC, columns)
        chosen = [n for n in ("pooled", "focalgroup", "refgroup") if vals.get(n)]
        if len(chosen) > 1:
            raise SeldecompError("contvars(): only one of pooled, focalgroup, and refgroup is allowed", 198)
        if chosen:
            sample = chosen[0]
        recut = bool(vals.get("recut"))
        words = left.split()
    out = []
    i = 0
    while i < len(words):
        cv = words[i]
        cq = words[i + 1] if i + 1 < len(words) else ""
        i += 2
        if cv.startswith("c."):
            cv = cv[2:]
        if columns is not None:
            try:
                vl = unab(cv, columns)
            except SeldecompError:
                if cv in ("pooled", "focalgroup", "refgroup", "recut"):
                    fix = " ".join("%s %s" % p for p in out) or "age 4"
                    raise SeldecompError("contvars(): put a comma before %s, as in contvars(%s, %s)"
                                         % (cv, fix, cv), 198)
                raise SeldecompError("contvars(): variable %s not found" % cv, 111)
            if len(vl) != 1:
                raise SeldecompError("contvars(): one variable name at a time (found %s)" % " ".join(vl), 198)
            cv = vl[0]
        if cv in [o[0] for o in out]:
            raise SeldecompError("contvars(): %s listed more than once" % cv, 198)
        if not _INTEGER.match(cq.lstrip("+")):
            raise SeldecompError("contvars(): %s must be followed by the number of quantile groups" % cv, 198)
        if int(cq) < 2:
            raise SeldecompError("contvars(): number of quantile groups for %s must be at least 2" % cv, 198)
        if df is not None and is_string_series(df[cv]):
            raise SeldecompError("contvars(): %s is not numeric" % cv, 108)
        out.append((cv, int(cq)))
    if not out:
        raise SeldecompError("contvars(): nothing specified", 198)
    return out, sample, recut, text


# ---- bayesboot() ---------------------------------------------------------------------------------
def parse_bayesboot(arg, columns=None):
    """bayesboot: True, "reps(200) rseed(1) citable" or a dict -> dict of suboptions, checked in
    Stata's order (the suboption list, reps(), dots(), rseed(), generate(), saving())."""
    d = {"reps": 100, "cluster": None, "rseed": None, "bca": False, "ties": False, "saving": None,
         "generate": None, "level": None, "citable": False, "dots": 1, "nodots": False,
         "normal": False, "tinterval": False, "percentile": False, "bc": False}
    d["_order"] = []
    if arg is True:
        return d
    if isinstance(arg, dict):
        for k, v in arg.items():
            if k.startswith("_"):
                d[k] = v
                continue
            name = spec_name(k, BAYESBOOT_SPEC)
            if name is None:
                raise SeldecompError("bayesboot(): option %s not allowed" % k, 198)
            d[name] = v
            if name == "dots":
                d["_dots"] = True
            if name in CI_TYPES and v:
                d["_order"].append(name)
        if d["cluster"] is not None and columns is not None:
            d["cluster"] = _vars_value(Opt("cluster", True, str(d["cluster"])), "cluster", "varname", columns)
        try:
            d["reps"] = int(d["reps"])
            d["dots"] = int(d["dots"])
        except (TypeError, ValueError):
            raise SeldecompError("option reps() incorrectly specified\n(error in option bayesboot())", 198)
        if d["level"] is not None:
            f = stata_number(d["level"])
            if f != f:
                raise SeldecompError("%s found where number expected\n(error in option bayesboot())"
                                     % d["level"], 198)
            if not 10 <= f <= 99.99:
                raise SeldecompError("level() must be between 10 and 99.99 inclusive\n"
                                     "(error in option bayesboot())", 198)
            if abs(f * 100 - round(f * 100)) > 1e-7:
                raise SeldecompError("level() can have at most two digits after the decimal point\n"
                                     "(error in option bayesboot())", 198)
            d["level"] = f
        if d["rseed"] is not None:
            d["rseed"] = str(d["rseed"])
    else:
        text = str(arg).strip()
        m = re.match(r"^([0-9]+)([ ,]|$)", text)
        if m:
            raise SeldecompError("bayesboot(%s) not allowed; type bayesboot(reps(%s))"
                                 % (m.group(1), m.group(1)), 198)
        vals, toks = subopt_syntax("bayesboot", text, BAYESBOOT_SPEC, columns)
        for k, v in vals.items():
            d[k] = v
        # the interval types in the order typed
        for t in bind_tokens(text):
            for rx, ty in (("^norm[a-z]*$", "normal"), ("^tint[a-z]*$", "tinterval"),
                           ("^perc[a-z]*$", "percentile"), ("^bc$", "bc"), ("^bca$", "bca")):
                if re.match(rx, t) and ty not in d["_order"]:
                    d["_order"].append(ty)
        if "dots" in vals:
            d["_dots"] = True
        d["nodots"] = bool(d["nodots"])
    if d["reps"] < 2:
        raise SeldecompError("bayesboot(reps()) must be at least 2", 198)
    if d["dots"] == 0:
        d["nodots"] = True
    elif d["dots"] < 0:
        raise SeldecompError("bayesboot(dots()) must be a nonnegative integer", 198)
    if d["dots"] == 0:
        d["dots"] = 1
    if d["rseed"] is not None:
        t = str(d["rseed"]).strip()
        if t.startswith("+"):
            t = t[1:]
        if not _INTEGER.match(t):
            raise SeldecompError("bayesboot(rseed()) must be an integer", 198)
        if not 0 <= int(t) <= 2147483647:
            raise SeldecompError("bayesboot(rseed()) must be an integer between 0 and 2147483647", 198)
        d["rseed"] = int(t)
    return d


def check_bs_generate(text):
    """bayesboot(generate(stub [, replace])) -> (stub, replace)."""
    g = str(text)
    k = find_top(g, ",")
    stub, rest = (g, "") if k < 0 else (g[:k], g[k + 1:])
    words = stub.split()
    ok = len(words) == 1 and valid_name(words[0])
    if ok:
        try:
            r = run_syntax(option_tokens(rest), (("replace", 7, "flag"),), None)
        except SeldecompError:
            ok = False
    if not ok:
        raise SeldecompError("bayesboot(generate()): type generate(stub) or generate(stub, replace)", 198)
    return words[0], bool(r.get("replace"))


def bs_saving_file(text):
    """bayesboot(saving(filename [, replace])) -> (the file, replace). The file is a CSV file;
    a name without an extension gets .csv."""
    s = str(text)
    k = find_top(s, ",")
    name, rest = (s, "") if k < 0 else (s[:k], s[k + 1:])
    if name.strip() == "":
        raise SeldecompError("filename required", 100)
    r = run_syntax(option_tokens(rest), (("replace", 7, "flag"),), None)
    toks = split_top(name.strip())
    if len(toks) > 1:
        raise SeldecompError("bayesboot(saving()): invalid file name; put names with spaces in double "
                             "quotes", 198)
    f = unquote(toks[0])
    if os.path.splitext(f)[1] == "":
        f += ".csv"
    return f, bool(r.get("replace"))


def check_bs_saving(text):
    """bayesboot(saving()), checked before any replication runs, with Stata's messages."""
    f, replace = bs_saving_file(text)
    if not replace and os.path.exists(f):
        raise SeldecompError("file %s already exists" % f, 602)
    if replace:
        folder = os.path.dirname(os.path.abspath(f))
        if (os.path.exists(f) and not os.access(f, os.W_OK)) or not os.access(folder, os.W_OK):
            raise SeldecompError("bayesboot(saving()): file %s cannot be written" % f, 603)


# ---- makeplot() ----------------------------------------------------------------------------------
_PLOT_OPTIONS_TEXT = ("the options are swapaxis, ci, normal, tinterval, percentile, bc, bca, reverse, "
                      "chainstep(), outcome(), format(), blocklabels(), and twoway()")


def parse_makeplot(arg):
    """makeplot: True, "waterfalllevels, blocklabels(demo Demographics) format(%9.2f)" or a
    dict -> dict(type, swapaxis, ci, citype, reverse, chainstep, outcome, format, blocklabels,
    twoway), checked as Stata's Plot_spec checks it. blocklabels: a list of (block, text)
    pairs, checked against the blocks by check_blocklabels."""
    d = {"type": None, "swapaxis": False, "ci": False, "citype": None, "reverse": False,
         "chainstep": None, "outcome": None, "format": None, "blocklabels": None, "twoway": None}
    if arg is True or arg is None:
        return d
    citypes = []
    if isinstance(arg, dict):
        for k, v in arg.items():
            if k == "type":
                d["type"] = v
            elif k in CI_TYPES:
                if v:
                    citypes.append(k)
            else:
                name = spec_name(k, MAKEPLOT_SPEC)
                if name is None:
                    raise SeldecompError("makeplot(): option %s not allowed; %s"
                                         % (k, _PLOT_OPTIONS_TEXT), 198)
                d[name] = v
        if d["type"] is not None and d["type"] not in PLOT_TYPES:
            t = str(d["type"])
            tw = "%s is not a plot type; " % t if _NAME.match(t) else ""
            raise SeldecompError("makeplot(): %stype the plot type first (stages, marginal, waterfall, or "
                                 "waterfalllevels), then the options after a comma" % tw, 198)
    else:
        s = str(arg).lstrip()
        opts = ""
        if s.startswith(","):
            opts = s[1:]
        elif s:
            m = re.match(r'"[^"]*"?|[^\s,]+', s)
            typ, rest = m.group(0), s[m.end():]
            if typ not in PLOT_TYPES:
                tw = "%s is not a plot type; " % typ if _NAME.match(typ) else ""
                raise SeldecompError("makeplot(): %stype the plot type first (stages, marginal, waterfall, or "
                                     "waterfalllevels), then the options after a comma" % tw, 198)
            d["type"] = typ
            rest = rest.lstrip()
            if rest.startswith(","):
                opts = rest[1:]
            elif rest:
                raise SeldecompError("makeplot(): put a comma between the plot type and the options, as in "
                                     "makeplot(%s, ci)" % typ, 198)
        try:
            toks = option_tokens(opts)
        except SeldecompError:
            raise SeldecompError("makeplot(): option not allowed; %s" % _PLOT_OPTIONS_TEXT, 198)
        left = loose_leftovers(toks, MAKEPLOT_SPEC)
        if left:
            bad = re.sub(r"[(].*$", "", left[0].text(), flags=re.S)
            bad = bad + " " if _NAME.match(bad) else ""
            raise SeldecompError("makeplot(): option %snot allowed; %s" % (bad, _PLOT_OPTIONS_TEXT), 198)
        for k, v in loose_values(toks, MAKEPLOT_SPEC).items():
            if k in CI_TYPES:
                if v:
                    citypes.append(k)
            elif v is not False:
                d[k] = v if v is True or k in ("blocklabels", "twoway") else unquote(v)
    if len(citypes) > 1:
        raise SeldecompError("makeplot(): specify only one interval type: normal, tinterval, percentile, bc, "
                             "or bca", 198)
    if citypes:
        if not d["ci"]:
            raise SeldecompError("makeplot(): specify ci with an interval type", 198)
        d["citype"] = citypes[0]
    if isinstance(d["twoway"], str):
        try:
            tw = loose_values(option_tokens(d["twoway"]), (("by", 2, "asis"),))
        except SeldecompError:
            tw = {}
        if tw.get("by"):
            raise SeldecompError("makeplot(): by() is not allowed inside twoway()", 198)
    if d["chainstep"] is not None:
        t = str(d["chainstep"]).strip()
        if not _INTEGER.match(t) or int(t) < 1:
            raise SeldecompError("makeplot(chainstep()) must be a positive integer", 198)
        d["chainstep"] = int(t)
    if d["outcome"] is not None:
        names = str(d["outcome"]).split()
        if len(names) > 1:
            raise SeldecompError("makeplot(outcome()): %s names %d outcomes; name one"
                                 % (d["outcome"], len(names)), 198)
        d["outcome"] = names[0] if names else None
    if d["format"] is not None:
        f = str(d["format"]).strip()
        ok = re.match(r"^%-?0?[0-9]+[.,][0-9]+[efg]c?$", f) is not None
        if ok:
            try:
                confirm_numeric_format(f)
            except SeldecompError:
                ok = False
        if not ok:
            raise SeldecompError("makeplot(format()): %s is not a numeric display format" % d["format"], 198)
    if isinstance(d["blocklabels"], str):
        toks = split_top(d["blocklabels"])
        if len(toks) % 2:
            raise SeldecompError('makeplot(blocklabels()): type pairs of a block name and its text, as in '
                                 'blocklabels(demo "Demographics")', 198)
        d["blocklabels"] = [(toks[i], unquote(toks[i + 1]) if toks[i + 1].startswith('"') else toks[i + 1])
                            for i in range(0, len(toks), 2)]
    elif isinstance(d["blocklabels"], dict):
        d["blocklabels"] = list(d["blocklabels"].items())
    return d


def check_blocklabels(pairs, blocks):
    """makeplot(blocklabels()): every name a block, no block twice, no empty text (Plot_blabels)
    -> {block: text}."""
    out = {}
    for nm, lab in pairs:
        if nm not in blocks:
            raise SeldecompError("makeplot(blocklabels()): %s is not a block; the blocks are %s"
                                 % (nm, " ".join(blocks)), 198)
        if nm in out:
            raise SeldecompError("makeplot(blocklabels()): block %s is given twice" % nm, 198)
        if lab == "":
            raise SeldecompError('makeplot(blocklabels()): the text of %s is empty; for a blank label, '
                                 'type " "' % nm, 198)
        out[nm] = lab
    return out


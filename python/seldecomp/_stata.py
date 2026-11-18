"""Numerical and formatting conventions shared with the Stata program."""
import math
import re
from decimal import ROUND_HALF_UP, Decimal, localcontext

import numpy as np

EPSFLOAT = 2.0 ** -23


def reldif(a, b):
    """|a - b| / (|b| + 1), as Stata's reldif()."""
    return np.abs(a - b) / (np.abs(b) + 1.0)


def pctile(x, w=None, nq=None, p=None):
    """Quantiles as Stata's _pctile computes them (its default definition).

    With sorted values x_(1..n) and cumulative weights W_(i), the quantile at fraction f is
    x_(i) for the first i with W_(i) > f * W, except that when W_(i-1) equals f * W exactly it
    is the average of x_(i-1) and x_(i). nq=q gives the q - 1 cut points of q quantile
    groups; p gives the percentiles listed (0 < p < 100).
    """
    x = np.asarray(x, dtype=float)
    if w is None:
        w = np.ones(len(x))
    w = np.asarray(w, dtype=float)
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    cw = np.cumsum(w[order])
    total = cw[-1]
    if nq is not None:
        fracs = [j / nq for j in range(1, nq)]
    else:
        fracs = [pp / 100.0 for pp in p]
    out = []
    n = len(xs)
    # the float cumulative sum can land just below an exact tie (0.2 + 0.7 < 0.9), so the
    # search and the tie test use a tolerance
    tolw = 1e-9 * total
    for f in fracs:
        target = f * total
        i = int(np.searchsorted(cw, target - tolw, side="left"))
        if i >= n:
            i = n - 1
        if abs(cw[i] - target) <= tolw and i + 1 < n:
            out.append((xs[i] + xs[i + 1]) / 2.0)
        else:
            out.append(xs[i])
    return out


def pctile_plain(x, p):
    """Unweighted percentiles as Stata's _pctile gives them: with the n values sorted, the
    average of x_(i) and x_(i+1) when n * p equals 100 * i (to a relative 1e-15), else
    x_(ceil(n * p / 100)). Each percentile is taken as a Stata macro holds it."""
    xs = np.sort(np.asarray(x, dtype=float))
    n = len(xs)
    out = []
    for pp in p:
        t = n * macro_value(pp)
        i = int(round(t / 100.0))
        if i >= 1 and abs(t - 100.0 * i) <= 1e-15 * (100.0 * i):
            out.append(xs[n - 1] if i >= n else (xs[i - 1] + xs[i]) / 2.0)
        else:
            i = int(math.ceil(t / 100.0))
            out.append(xs[min(max(i, 1), n) - 1])
    return [float(v) for v in out]


def macro_value(x):
    """A number as a Stata macro holds it (the text of format %18.0g, read back)."""
    x = float(x)
    if x != x or math.isinf(x):
        return x
    return float(fmt_g(x, 18))


_NA = (2.2352520354606839287, 161.02823106855587881, 1067.6894854603709582, 18154.981253343561249,
       0.065682337918207449113)
_NB = (47.20258190468824187, 976.09855173777669322, 10260.932208618978205, 45507.789335026729956)
_NC = (0.39894151208813466764, 8.8831497943883759412, 93.506656132177855979, 597.27027639480026226,
       2494.5375852903726711, 6848.1904505362823326, 11602.651437647350124, 9842.7148383839780218,
       1.0765576773720192317e-8)
_ND = (22.266688044328115691, 235.38790178262499861, 1519.377599407554805, 6485.558298266760755,
       18615.571640885098091, 34900.952721145977266, 38912.003286093271411, 19685.429676859990727)
_NP = (0.21589853405795699, 0.1274011611602473639, 0.022235277870649807, 0.001421619193227893466,
       2.9112874951168792e-5, 0.02307344176494017303)
_NQ = (1.28426009614491121, 0.468238212480865118, 0.0659881378689285515, 0.00378239633202758244,
       7.29751555083966205e-5)


def normal(x):
    """The standard normal distribution function, computed as Stata's normal() computes it
    (Cody's rational approximations), so that the two agree to the last bit."""
    x = float(x)
    if x != x:
        return x
    y = abs(x)
    if y <= 0.67448975:
        xnum = xden = 0.0
        if y > 1.11e-16:
            xsq = x * x
            xnum = _NA[4] * xsq
            xden = xsq
            for i in range(3):
                xnum = (xnum + _NA[i]) * xsq
                xden = (xden + _NB[i]) * xsq
        return 0.5 + x * (xnum + _NA[3]) / (xden + _NB[3])
    if y > 38.0:
        return 0.0 if x < 0 else 1.0
    if y <= 5.656854249492380195206754896838:
        xnum = _NC[8] * y
        xden = y
        for i in range(7):
            xnum = (xnum + _NC[i]) * y
            xden = (xden + _ND[i]) * y
        temp = (xnum + _NC[7]) / (xden + _ND[7])
    else:
        xsq = 1.0 / (x * x)
        xnum = _NP[5] * xsq
        xden = xsq
        for i in range(4):
            xnum = (xnum + _NP[i]) * xsq
            xden = (xden + _NQ[i]) * xsq
        temp = xsq * (xnum + _NP[4]) / (xden + _NQ[4])
        temp = (0.398942280401432677939946059934 - temp) / y
    xsq = math.trunc(y * 16) / 16
    cum = math.exp(-xsq * xsq * 0.5) * math.exp(-(y - xsq) * (y + xsq) * 0.5) * temp
    return cum if x <= 0 else 1.0 - cum


_IA = (3.387132872796366608, 133.14166789178437745, 1971.5909503065514427, 13731.693765509461125,
       45921.953931549871457, 67265.770927008700853, 33430.575583588128105, 2509.0809287301226727)
_IB = (1.0, 42.313330701600911252, 687.1870074920579083, 5394.1960214247511077, 21213.794301586595867,
       39307.89580009271061, 28729.085735721942674, 5226.495278852854561)
_IC = (1.42343711074968357734, 4.6303378461565452959, 5.7694972214606914055, 3.64784832476320460504,
       1.27045825245236838258, 0.24178072517745061177, 0.0227238449892691845833, 7.7454501427834140764e-4)
_ID = (1.0, 2.05319162663775882187, 1.6763848301838038494, 0.68976733498510000455, 0.14810397642748007459,
       0.0151986665636164571966, 5.475938084995344946e-4, 1.05075007164441684324e-9)
_IE = (6.6579046435011037772, 5.4637849111641143699, 1.7848265399172913358, 0.29656057182850489123,
       0.026532189526576123093, 0.0012426609473880784386, 2.71155556874348757815e-5, 2.01033439929228813265e-7)
_IF = (1.0, 0.59983220655588793769, 0.13692988092273580531, 0.0148753612908506148525, 7.868691311456132591e-4,
       1.8463183175100546818e-5, 1.4215117583164458887e-7, 2.04426310338993978564e-15)


def _poly(co, r):
    v = co[7]
    for k in range(6, -1, -1):
        v = v * r + co[k]
    return v


def invnormal(p):
    """The inverse of the standard normal distribution function, computed as Stata's
    invnormal() computes it (Wichura's algorithm AS 241), so that the two agree to the last bit."""
    p = float(p)
    if p != p or p <= 0.0 or p >= 1.0:
        return float("nan")
    q = p - 0.5
    if abs(q) <= 0.425:
        r = 0.180625 - q * q
        return q * _poly(_IA, r) / _poly(_IB, r)
    r = math.sqrt(-math.log(p if q < 0 else 1.0 - p))
    if r <= 5.0:
        v = _poly(_IC, r - 1.6) / _poly(_ID, r - 1.6)
    else:
        v = _poly(_IE, r - 5.0) / _poly(_IF, r - 5.0)
    return -v if q < 0 else v


def _missing(x):
    return x is None or (isinstance(x, float) and math.isnan(x))


def _fixed(ax, ndec):
    """ax (>= 0) with ndec decimals, an exact half rounded up (as the g formats round it)."""
    with localcontext() as c:
        c.prec = 400
        q = Decimal(ax).quantize(Decimal(1).scaleb(-ndec), rounding=ROUND_HALF_UP)
    return "{:.{}f}".format(q, ndec)


def _commas(s):
    """Thousands separators in the integer part of a number's text."""
    sign = "-" if s.startswith("-") else ""
    body = s[len(sign):]
    ip, dot, fp = body.partition(".")
    if len(ip) > 3:
        ip = "{:,}".format(int(ip))
    return sign + ip + dot + fp


def _g_exp(ax, width, dec):
    """The exponential form of %<width>.<dec>g (without the sign): as many mantissa decimals as
    fit in width - 1 columns, at least one (none with a three-digit exponent that leaves no
    room), at most dec - 1 when dec is given."""
    nexp = 3 if ax >= 1e100 or ax < 1e-99 else 2
    edec = max(1, width - 7) - (1 if nexp == 3 else 0)
    if dec > 0:
        edec = min(edec, dec - 1)
    edec = max(edec, 0)
    s = "%.*e" % (edec, ax)
    if edec == 0:
        m, e = s.split("e")
        s = m + ".e" + e
    return s


def _g_text(x, width, dec, comma):
    """(text, exponential?) of x in %<width>.<dec>g[c], without the padding."""
    ax = abs(x)
    sign = "-" if x < 0 else ""
    X = math.floor(math.log10(ax))
    if 10.0 ** X > ax:
        X -= 1
    elif 10.0 ** (X + 1) <= ax:
        X += 1
    A = width - 1
    P = width - 2
    for _ in range(2):
        if X >= P or X < -5 or A + X < 2:
            return sign + _g_exp(ax, width, dec), True
        usecomma = False
        if X >= 0:
            digits = P if dec == 0 else max(min(dec, P), X + 1)
            if comma and X >= 3:
                nc = X // 3
                if min(digits, P - nc) >= X + 1:
                    digits = min(digits, P - nc)
                    usecomma = True
        else:
            digits = min(P, A + X)
            if dec > 0:
                digits = min(digits, dec)
        ndec = max(digits - X - 1, 0)
        s = _fixed(ax, ndec)
        # rounding that adds an integer digit (99999999.5) can call for the exponential form
        if X >= 0 and len(s.split(".")[0]) > X + 1:
            if dec > 0:
                return sign + _g_exp(ax, width, dec), True
            X += 1
            continue
        break
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    if s.startswith("0."):
        s = s[1:]
    if usecomma:
        s = _commas(s)
    return sign + s, False


def fmt_gs(x, width=9, dec=0, comma=False):
    """A number as Stata's %<width>.<dec>g[c] shows it, without the padding (as string() gives
    it). Fixed notation shows width - 2 significant digits at most (dec at most when dec > 0,
    though never fewer than the integer digits), within width - 1 columns, without trailing
    zeros or a leading zero; exponential notation when the integer part has width - 2 digits
    or more, or the first significant digit lies beyond the fifth decimal or leaves room for
    fewer than two digits. With comma, thousands separators when they fit without dropping
    integer digits."""
    if _missing(x) or (isinstance(x, float) and math.isinf(x)):
        return "."
    x = float(x)
    if x == 0:
        return "0"
    return _g_text(x, width, dec, comma)[0]


def fmt_g(x, width=9):
    """A number as Stata's %<width>.0g shows it, without the padding."""
    return fmt_gs(x, width)


NUMERIC_FMT = re.compile(r"^%(-?)(0?)([1-9][0-9]*)([.,])([0-9]+)([fg]c?|e)$")


def numeric_format(fmt):
    """Whether fmt is a format fmt_num shows: %[-][0]w.d with f, g or e (c after f or g for
    commas; a comma in place of the period for a decimal comma), d smaller than w."""
    m = NUMERIC_FMT.match(str(fmt).strip())
    return m is not None and int(m.group(5)) < int(m.group(3))


def parse_fmt(fmt):
    """A Stata numeric format -> (left-aligned, width, decimals, kind, comma)."""
    return _parse_fmt(fmt)[:5]


def _parse_fmt(fmt):
    """-> (left-aligned, width, decimals, kind, comma, zero-padded, decimal comma)"""
    m = NUMERIC_FMT.match(fmt.strip())
    if m is None:
        raise ValueError("format %s" % fmt)
    kind = m.group(6)
    return (m.group(1) == "-", int(m.group(3)), int(m.group(5)), kind[0], kind.endswith("c"),
            m.group(2) == "0", m.group(4) == ",")


def fmt_width(fmt):
    """The width of a Stata format (fmtwidth())."""
    return parse_fmt(fmt)[1]


def _e_text(x, width, dec, left):
    """%w.de: a column for the sign, and as many decimals as asked for (.0: as many as fit)
    that fit in the width (at least 8 columns are used)."""
    ax = abs(x)
    room = max(width, 8) - 7
    d = min(dec, room) if dec > 0 else room
    while True:
        t = "%.*e" % (max(d, 0), ax)
        extra = len(t.split("e")[1]) - 3
        fit = min(dec, room - extra) if dec > 0 else room - extra
        if fit >= d or d <= 0:
            break
        d = fit
    if d <= 0:
        t = t.replace("e", ".e", 1) if "." not in t else t
    sign = "-" if x < 0 else ("" if left else " ")
    return sign + t


def fmt_num(x, fmt):
    """A number in a Stata numeric format (%w.df, %w.dg, %w.de, with optional comma: %9.0fc),
    padded to its width as display shows it (a missing value is "."), or, when fmt does not
    start with %, a Python format specification. A fixed text too wide for the format is shown
    whole up to max(width, 8) columns, beyond that in exponential notation, which keeps a
    column for the sign."""
    if not fmt.startswith("%"):
        return format(x, fmt)
    left, width, dec, kind, comma, zero, euro = _parse_fmt(fmt)
    expo = False
    if _missing(x) or (isinstance(x, float) and math.isinf(x)):
        s = "."
    elif kind == "f":
        x = float(x)
        s = "%.*f" % (dec, x)
        if comma:
            sc = _commas(s)
            if len(sc) <= width:
                s = sc
        if len(s) > max(width, 8):
            s = ("-" if x < 0 else "") + _g_exp(abs(x), width, dec + 1 if dec > 0 else 0)
            expo = True
    elif kind == "e":
        s = _e_text(float(x), width, dec, left)
    elif float(x) == 0:
        s = "0"
    else:
        s, expo = _g_text(float(x), width, dec, comma)
    if expo and not s.startswith("-"):
        width = max(width, len(s) + 1)
    if zero and kind == "f" and not expo and s != "." and len(s) < width:
        # %0w.df: zeros between the sign and the digits
        s = ("-" if s.startswith("-") else "") + s.lstrip("-").rjust(width - (1 if s.startswith("-") else 0), "0")
    if euro and s != ".":
        s = s.translate({46: 44, 44: 46})
    return s.ljust(width) if left else s.rjust(width)


def fmt_str(x, fmt):
    """string(x, fmt): the number in the format, without the padding."""
    return fmt_num(x, fmt).strip()


def stata_round(x, y=1.0):
    """Stata's round(x, y): the multiple of y nearest x, halves away from zero."""
    if _missing(x):
        return x
    q = x / y
    r = math.floor(abs(q) + 0.5)
    return math.copysign(r, q) * y


def plural(n, word):
    """word, with an s unless n is 1 (Stata's plural())."""
    return word if n == 1 else word + "s"

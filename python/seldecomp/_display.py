"""The printed output: the text of the Stata command's display, line for line (the programs
Display ... Display_covdetail of seldecomp.ado, whose names the methods below keep), at 78
columns. Plain text: what SMCL renders, with {p a b c 78} paragraphs wrapped here."""
import math
import re

import numpy as np
from scipy import stats

from ._stata import fmt_num, fmt_str, fmt_width, invnormal, normal, plural, stata_round

NAN = float("nan")
RULE = "-" * 78
TOP = "-" * 17 + "-" + "-" * 60          # {hline 17}{c TT}{hline 60}
MID = "-" * 17 + "+" + "-" * 60          # {hline 17}{c +}{hline 60}
BOT = TOP                                # {hline 17}{c BT}{hline 60}
CI_KEYS = ("ci_normal", "ci_t", "ci_percentile", "ci_bc", "ci_bca", "pct_explained_ci")
H_COEF = "   Estimate"
H_SEL = "std. err."
H_CI = "conf. interval"


# ---- the primitives: display's columns and formats, SMCL paragraphs ----------------------

class _Col(int):
    """_col(#) in a display line."""


def col(n):
    return _Col(n)


def di(*parts):
    """One display line: text, col(#) (move to column #, never back), and (format, value)
    pairs (a number in a Stata format, or a string in %#s / %-#s)."""
    s = ""
    for p in parts:
        if isinstance(p, _Col):
            if len(s) < p - 1:
                s += " " * (p - 1 - len(s))
        elif isinstance(p, tuple):
            f, v = p
            if f.endswith("s"):
                w = int(f[1:-1].lstrip("-") or 0)
                v = str(v)
                s += v.ljust(w) if f[1] == "-" else v.rjust(w)
            else:
                s += fmt_num(v, f)
        else:
            s += p
    return s


def para(text, first, rest, right=2, width=78):
    """{p first rest right width}text{p_end}: the words (runs of blanks collapse) filled into
    lines that end by column width - right, the first indented first columns, the others
    rest."""
    words = text.split()
    lim = width - right
    lines = []
    cur = " " * first
    has = False
    for w in words:
        if has and len(cur) + 1 + len(w) > lim:
            lines.append(cur)
            cur = " " * rest + w
        else:
            cur = cur + (" " if has else "") + w
        has = True
    lines.append(cur)
    return lines


def ralign(text, n):
    """{ralign n:text}"""
    return text.rjust(n)


def _miss(x):
    return x is None or (isinstance(x, float) and math.isnan(x))


def _num(x):
    """A value as a float, missing as nan."""
    if _miss(x):
        return NAN
    try:
        return float(x)
    except (TypeError, ValueError):
        return NAN


def _big(x):
    """A value for a comparison as Stata makes it: missing is larger than any number."""
    x = _num(x)
    return math.inf if math.isnan(x) else x


def _g12(x):
    """string(x, "%12.0g")"""
    return fmt_str(_num(x), "%12.0g")


def abbrev(s, n):
    """Stata's abbrev(): a name longer than n characters (its first 32 are read) cut to n,
    with ~ before its last character and the digits that come just before it, or before the
    digits it ends with."""
    if n < 5:
        n = 5
    s = s[:32]
    if len(s) <= n:
        return s
    k = len(s) - 1
    if not s[k].isdigit():
        k -= 1
    while k >= 0 and s[k].isdigit():
        k -= 1
    tail = s[k + 1:]
    if not tail or len(tail) > n - 2:
        tail = s[-1]
    return s[:n - 1 - len(tail)] + "~" + tail


def vltrim(s, n):
    """seldecomp_vltrim: a text cut to n characters, ending with ~ when it is cut."""
    return s if len(s) <= n else s[:n - 1] + "~"


def _inv_normal(level):
    return invnormal(1 - (100 - level) / 200.0)


def _inv_t(df, level):
    return float(stats.t.isf((100 - level) / 200.0, df))


def ixlab(lab, v):
    """seldecomp_ixlab: the label of a cell of an interaction (v1=a v2=b) as its components'
    value labels (or values), joined by " # "."""
    if not lab.startswith(v[0] + "="):
        return lab
    rest = lab[len(v[0]) + 1:]
    pc = [""] * len(v)
    for j in range(1, len(v)):
        p = rest.find(" " + v[j] + "=")
        if p < 0:
            return lab
        pc[j - 1] = rest[:p]
        rest = rest[p + len(v[j]) + 2:]
    pc[-1] = rest
    out = []
    for piece in pc:
        m = re.match(r"^-?[0-9]+ [(](.*)[)]$", piece)
        out.append(m.group(1) if m else piece)
    return " # ".join(out)


# ---- the display ----------------------------------------------------------------------

class _Display:
    """One display of the results: the e() results, the varlabels switch, the lines."""

    def __init__(self, res, varlabels=False):
        self.res = res
        self.E = res.e
        self.M = res.matrices
        self.ctx = getattr(res, "_ctx", None)
        self.vl = bool(varlabels)
        self.out = []

    # e(), as Stata reads it: an absent scalar is missing, an absent macro is ""
    def ev(self, k):
        return _num(self.E.get(k))

    def es(self, k):
        v = self.E.get(k)
        return "" if v is None else str(v)

    def has(self, name):
        return name in self.M and self.M[name] is not None

    def mat(self, name):
        return self.M[name].to_numpy(dtype=float)

    def bvec(self):
        b = self.M["b"]
        return list(b.columns), b.iloc[0].to_numpy(dtype=float)

    def add(self, *lines):
        self.out.extend(lines)

    def p(self, text, first=0, rest=4, right=2, width=78):
        self.out.extend(para(text, first, rest, right, width))

    # ---- labels (Vlab, Block_lab, Row_lab, Eq_lab, Eq_title) ----
    def vl1(self, v):
        labs = getattr(self.ctx, "labels", None) or {}
        lab = labs.get(v)
        if lab is None or str(lab).strip() == "":
            return v
        return str(lab)

    def vlab(self, s):
        """seldecomp_vlab: with varlabels, each name of a list replaced by its variable label,
        the parts of a term joined by " # ", the terms by ", "."""
        if not self.vl:
            return s
        out = []
        for t in s.split():
            out.append(" # ".join(self.vl1(q) for q in t.split("#")))
        return ", ".join(out)

    def Vlab(self, names, w=None):
        t = self.vlab(names)
        if w is not None:
            t = abbrev(names, w) if t == names else vltrim(t, w)
        return t

    def Block_lab(self, bn, w=None):
        blks = self.es("blocks").split()
        t = bn
        if bn in blks and self.vl:
            terms = self.es("block%d" % (blks.index(bn) + 1))
            if bn == terms or bn == terms.replace("#", "_"):
                t = terms
        return self.Vlab(t, w)

    def Row_lab(self, eq, nm):
        if eq in ("overall", ""):
            mtxt = "ratio" if self.es("ratio") else "mean"
            t = abbrev(nm, 16)
            if nm == "focal_group":
                t = "Focal " + mtxt
            elif nm == "raked_focal":
                t = "Raked ratio" if self.es("ratio") else "Raked focal mean"
            elif nm == "ref_group":
                t = "Reference " + mtxt
            elif nm == "difference":
                t = "Difference"
            elif nm == "explained":
                t = "Explained"
            elif nm == "unexplained":
                t = "Unexplained"
            return t
        if eq in ("detail", "rev_detail"):
            return self.Vlab(nm, 16)
        return self.Block_lab(nm, 16)

    def Eq_lab(self, eq):
        if eq.startswith("rev_"):
            eq = eq[4:]
        if eq in ("overall", ""):
            return "Overall"
        if eq == "detail":
            return "By term"
        if self.es("stagetype") == "individually":
            return "By block"
        return "By stage"

    def Eq_title(self, eq):
        rev = eq.startswith("rev_")
        if rev:
            eq = eq[4:]
        ind = self.es("stagetype") == "individually"
        mtxt = "ratio" if self.es("ratio") else "mean"
        t = ""
        if eq == "marg_expl":
            t = "Contribution of each block raked alone, and the remainder" if ind else "Marginal explained"
        elif eq == "cum_expl":
            t = "Cumulative explained"
        elif eq == "cum_unexpl":
            t = "Cumulative unexplained"
        elif eq == "raked_focal":
            t = ("Raked focal group %s, each block raked alone" % mtxt) if ind else "Raked focal group " + mtxt
        elif eq == "detail":
            t = "Marginal explained, by term"
        if rev and t:
            t += ", reverse order"
        return " " + t if t else ""

    def vlnote(self, s):
        """seldecomp_vlnote: with varlabels, the variable names of a note's level listings (a
        name followed by " = " or "#", or after "#") replaced by their labels."""
        if not self.vl:
            return s
        out = ""
        rest = s
        while True:
            m = re.search(r"[A-Za-z_][A-Za-z0-9_]*", rest)
            if not m:
                break
            w = m.group(0)
            p = rest.find(w)
            out += rest[:p]
            rest = rest[p + len(w):]
            if rest[:3] == " = " or rest[:1] == "#" or out[-1:] == "#":
                out += self.vl1(w)
            else:
                out += w
        return out + rest

    # ---- Display_eps, Display_pct, Display_pct6, Display_fmts, Display_vce ----
    def eps(self, o=1):
        if self.es("outcomes"):
            eq0 = [self.es("outcomes").split()[o - 1]]
        else:
            eq0 = ["overall"]
        if self.es("chain"):
            eq0 = ["p%d" % o] if o != 0 else ["p%d" % j for j in range(1, int(self.ev("N_pairs")) + 1)]
        names, vals = self.bvec()
        pos = {n: i for i, n in enumerate(names)}
        s = 0.0
        for eq in eq0:
            for r in ("focal_group", "raked_focal", "ref_group"):
                k = pos.get("%s:%s" % (eq, r))
                if k is not None and not math.isnan(vals[k]):
                    s = max(s, abs(vals[k]))
        return 1e-13 * s

    @staticmethod
    def pct(x):
        s = fmt_str(x, "%9.1f")
        try:
            if float(s) == 0:
                s = "0.0"
        except ValueError:
            pass
        return s

    @staticmethod
    def pct6(x):
        if _miss(x):
            return "."
        s = fmt_str(x, "%6.1f")
        if abs(x) >= 1000 or len(s) > 6:
            s = fmt_str(x, "%6.0f")
        if abs(x) >= 1e5 or len(s) > 6:
            sg = "-" if x < 0 else ""
            k = math.floor(math.log10(abs(x)))
            m = stata_round(abs(x) / 10.0 ** k, .1)
            if m >= 10:
                m = m / 10
                k += 1
            s = sg + fmt_str(m, "%3.1f") + "e" + str(k)
            if len(s) > 6:
                m = stata_round(abs(x) / 10.0 ** k)
                if m >= 10:
                    m = 1
                    k += 1
                s = sg + fmt_str(m, "%12.0g") + "e" + str(k)
        return s

    def fmts(self, cformat=None, pformat=None, sformat=None, check=False):
        out = {}
        for f, given, maxw, dflt in (("c", cformat, 9, "%9.0g"), ("p", pformat, 5, "%5.3f"),
                                     ("s", sformat, 8, "%8.2f")):
            fmt = "" if given in (None, False) else str(given).strip()
            if fmt:
                from ._parse import check_display_format, other_format_width
                fmt = check_display_format(f + "format", fmt)
                if (other_format_width(fmt) or fmt_width(fmt)) > maxw:
                    if not check:
                        self.add("note: invalid %sformat(), using default" % f)
                    fmt = ""
            out[f] = fmt or dflt
        return out["c"], out["p"], out["s"]

    def vce(self):
        lab = self.es("vcetype")
        head = lab
        bs = self.es("vce") in ("bootstrap", "bayesboot")
        if self.es("bsnames"):
            lab, head, bs = "Bayesian bootstrap", "Bootstrap", True
        elif self.es("vce") == "bayesboot":
            lab = head = "Bayesian bootstrap"
        return lab, head, bs

    def sep(self):
        self.add("", "=" * 78)

    def jkfail(self):
        n = self.ev("N_jkfail")
        if n > 0 and not math.isnan(n):
            self.p("Caution: %s jackknife replicates failed; BCa intervals may be incorrect." % _g12(n))

    def fold_has(self, txt):
        return any(txt in str(nt) for nt in self.res.notes)

    # ---- Display_ci: the intervals at a level -------------------------------------------
    def display_ci(self, level):
        """The bootstrap intervals at level: e()'s at the estimation level, else recomputed
        from the replicates (e(bsdraws)) as the estimation computes them. Returns (dict of
        arrays: ci_<type> 2 x k over e(bsnames), pct_explained_ci rows x 2; the level text)."""
        from ._boot import summarize
        bca = self.ev("bca") == 1
        nout = max(1, len(self.es("outcomes").split()))
        ch = bool(self.es("chain"))
        if ch:
            nout = int(self.ev("N_pairs")) + 1
        ci = {}
        if level == self.ev("level") or not self.has("bsdraws"):
            for m in ("ci_normal", "ci_t", "ci_percentile", "ci_bc") + (("ci_bca",) if bca else ()):
                if self.has(m):
                    ci[m] = self.mat(m)
            if ch:
                BP = self.M["by_pair"]
                top = BP[["pct_explained_ll", "pct_explained_ul"]].to_numpy(dtype=float)
                ci["pct_explained_ci"] = np.vstack([top, self.mat("pct_explained_ci")[:1, :2]])
                ci["_pct_rows"] = list(BP.index) + ["total"]
            elif nout > 1:
                ci["pct_explained_ci"] = self.mat("by_outcome")[:, 8:10]
                ci["_pct_rows"] = list(self.M["by_outcome"].index)
            elif self.has("pct_explained_ci"):
                ci["pct_explained_ci"] = self.mat("pct_explained_ci")[:1, :2]
                ci["_pct_rows"] = ["r1"]
            return ci, _g12(self.ev("level"))
        names, vals = self.bvec()
        pos = {n: i for i, n in enumerate(names)}
        bsn = self.es("bsnames").split()
        ob = np.array([vals[pos[n]] for n in bsn])
        draws = self.mat("bsdraws")
        haspct = self.has("bsdraws_pct")
        if haspct:
            pdraws = self.mat("bsdraws_pct")
            if ch:
                pob = list(self.M["by_pair"].iloc[:, 6].to_numpy(dtype=float)) + [self.ev("pct_explained")]
            elif nout > 1:
                pob = list(self.mat("by_outcome")[:, 6])
            else:
                pob = [self.ev("pct_explained")]
        else:
            pdraws = np.zeros((draws.shape[0], 0))
            pob = []
        accel = self.mat("accel")[0] if (bca and self.has("accel")) else None
        S, pse, pci = summarize(draws, pdraws, ob, pob, float(level), self.ev("df_bs"),
                                self.ev("ties") == 1, accel)
        for m in ("ci_normal", "ci_t", "ci_percentile", "ci_bc") + (("ci_bca",) if bca else ()):
            ci[m] = S[m]
        if haspct:
            ci["pct_explained_ci"] = np.array([list(c) for c in pci], dtype=float).reshape(-1, 2)
            if ch:
                ci["_pct_rows"] = ["p%d" % k for k in range(1, nout)] + ["total"]
            elif nout > 1:
                ci["_pct_rows"] = self.es("outcomes").split()
            else:
                ci["_pct_rows"] = ["r1"]
        return ci, _g12(level)

    # ---- Display_header -----------------------------------------------------------------
    def header(self):
        """The title and counts, the groups, weights, inference, the mode lines, the chain's pairs
        and the notes."""
        by = self.Vlab(self.es("by"))
        ch = bool(self.es("chain"))
        g = [self.es("focal_group_label"), self.es("ref_group_label")]
        gw = [24 + len(by) + len(x) for x in g]
        gfit = all(w <= 47 for w in gw)
        self.add("", di(self.es("title"), col(49), "Number of obs", col(67), "= ", ("%10.0fc", self.ev("N"))))
        if ch:
            t = "Chain of adjacent pairs of " + by
            if len(t) <= 47:
                self.add(di(t, col(49), "Pairs", col(67), "= ", ("%10.0fc", self.ev("N_pairs"))))
            else:
                self.add(di(col(49), "Pairs", col(67), "= ", ("%10.0fc", self.ev("N_pairs"))), t)
        elif gfit:
            self.add(di("Focal Group (raked): %s = %s" % (by, g[0]), col(49), "Focal obs", col(67), "= ",
                        ("%10.0fc", self.ev("N_focal"))))
            self.add(di("Reference Group:     %s = %s" % (by, g[1]), col(49), "Reference obs", col(67), "= ",
                        ("%10.0fc", self.ev("N_ref"))))
        else:
            self.add(di(col(49), "Focal obs", col(67), "= ", ("%10.0fc", self.ev("N_focal"))))
            self.add(di(col(49), "Reference obs", col(67), "= ", ("%10.0fc", self.ev("N_ref"))))
        wfit = True
        wtp = self.es("wtype")
        wtxt = ""
        if wtp and ch:
            wtxt = "Weights: %s %s" % (self.es("wtype"), self.es("wexp"))
            wfit = False
        elif wtp:
            wtxt = "Weights: %s %s" % (self.es("wtype"), self.es("wexp"))
            wfit = len(wtxt) <= 47
            wcol, wf = 69, "%10.0gc"
            if max(_big(self.ev("W_focal")), _big(self.ev("W_ref"))) >= 1e7:
                wcol, wf = 68, "%11.0fc"
            if wfit:
                self.add(di(wtxt, col(49), "Focal wgt sum", col(67), "=", col(wcol), (wf, self.ev("W_focal"))))
            else:
                self.add(di(col(49), "Focal wgt sum", col(67), "=", col(wcol), (wf, self.ev("W_focal"))))
            self.add(di(col(49), "Reference wgt sum", col(67), "=", col(wcol), (wf, self.ev("W_ref"))))
        if not gfit and not ch:
            for k, gt in ((0, "Focal Group (raked): "), (1, "Reference Group:     ")):
                t = "%s%s = %s" % (gt, by, g[k])
                if gw[k] > 78:
                    self.p(t, 0, 21)
                else:
                    self.add(t)
        if not wfit:
            if len(wtxt) > 78:
                self.p(wtxt, 0, 9)
            else:
                self.add(wtxt)
        vcelab, vcehead, vcebs = self.vce()
        if vcelab:
            itxt = "Inference: " + vcelab
            iw = len(itxt)
            nreps = self.ev("N_reps")
            if math.isnan(nreps):
                self.p(itxt, 0, 11) if iw > 78 else self.add(itxt)
            elif iw <= 47:
                self.add(di(itxt, col(49), "Replications", col(67), "= ", ("%10.0fc", nreps)))
            else:
                self.add(di(col(49), "Replications", col(67), "= ", ("%10.0fc", nreps)))
            if self.es("bsnames") and not math.isnan(self.ev("N_clust")):
                cvl = self.Vlab(self.es("clustvar"))
                self.add(ralign("(Replications based on %s clusters in %s)"
                                % (fmt_str(self.ev("N_clust"), "%12.0fc"), cvl), 78))
            if not math.isnan(nreps) and iw > 47:
                self.p(itxt, 0, 11) if iw > 78 else self.add(itxt)
        if self.es("outcomes"):
            self.p("Outcomes (one raking): " + self.Vlab(self.es("outcomes")))
        if self.es("ratio"):
            num = self.Vlab(self.es("depvar"))
            if self.es("outcomes"):
                num = "each outcome"
            self.p("Ratio of means: %s / %s" % (num, self.Vlab(self.es("ratio"))))
        how = "others held at focal group margins" if self.es("fixothermargins") else "others unconstrained"
        blocks = self.es("blocks").split()
        st = self.es("stagetype")
        if st == "individually":
            self.p("Blocks raked individually (%s):" % how)
            for s in range(1, int(self.ev("N_blocks")) + 1):
                bn = blocks[s - 1]
                bd = self.es("block%d" % s)
                bnl, bdl = self.Block_lab(bn), self.Vlab(bd)
                lead = di("block ", ("%2.0f", s), ": ") + bnl
                self.p(lead if (bn == bd or bnl == bdl) else lead + " (%s)" % bdl, 2, 14)
        elif st:
            ns = int(self.ev("N_stages"))
            self.p("Blocks added sequentially (%s):" % how)
            for s in range(1, ns + 1):
                bn = blocks[s - 1]
                bd = self.es("block%d" % s)
                bnl, bdl = self.Block_lab(bn), self.Vlab(bd)
                lead = di("stage ", ("%2.0f", s), ": + ") + bnl
                self.p(lead if (bn == bd or bnl == bdl) else lead + " (%s)" % bdl, 2, 14)
            if self.es("reverse"):
                rb = " ".join("+ " + self.Block_lab(blocks[s - 1]) for s in range(ns, 0, -1))
                self.p("and in reverse order: " + rb, 2, 4)
            if self.es("detail"):
                self.p("and term by term within each block (detail)", 2, 4)
        else:
            self.p("Raking variables: " + self.Vlab(self.es("rakevars")), 0, 18)
        if self.es("insupport") and ch:
            BP = self.M["by_pair"]
            nx = float((BP["N_focal_all"] - BP["N_focal"]).sum())
            if nx > 0:
                self.p("Support (insupport): %s focal group obs excluded over the pairs (levels absent "
                       "from the pair's reference group)" % fmt_str(nx, "%12.0fc"))
            else:
                self.p("Support (insupport): no observation excluded")
        elif self.es("insupport"):
            nx = self.ev("N_focal_all") - self.ev("N_focal")
            if nx > 0:
                dvl = self.Vlab(self.es("depvar").split()[0])
                mtxt = ("ratio" if self.es("ratio") else "mean") + " of " + dvl
                self.p("Support (insupport): %s of %s focal group obs excluded (levels absent from the "
                       "reference group); focal group %s before the exclusion: %s"
                       % (fmt_str(nx, "%12.0fc"), fmt_str(self.ev("N_focal_all"), "%12.0fc"), mtxt,
                          fmt_str(self.ev("mean_focal_all"), "%9.0g")))
            else:
                self.p("Support (insupport): no observation excluded")
        if self.es("cutvars"):
            parts = []
            cq = self.es("cutq").split()
            cb = self.es("cutbins").split()
            for c, cv in enumerate(self.es("cutvars").split()):
                q = cq[c]
                qtxt = "%s quantile bins" % q
                if c < len(cb) and int(cb[c]) < int(q):
                    qtxt = "%s nonempty of %s quantile bins (ties)" % (cb[c], q)
                parts.append("%s: %s" % (self.Vlab(cv), qtxt))
            sepc = ";" if len(parts) > 1 else ","
            csx = ("each pair's %s sample" if ch else "the %s sample") % self.es("contsample")
            self.p("(%s%s cut points from %s)" % (", ".join(parts), sepc, csx), 2, 3)
        if ch:
            BP = self.M["by_pair"]
            labs = [x.strip() for x in self.es("by_labels").split("|")]
            bypfx = by + " = " if len(by) + 3 + max(len(x) for x in labs) <= 26 else ""
            self.add("", di("Pair  Focal Group (raked)", col(34), "Reference Group", col(61), "Focal obs",
                            col(71), "Ref. obs"))
            nf = BP["N_focal"].to_numpy(dtype=float)
            nr = BP["N_ref"].to_numpy(dtype=float)
            for j in range(len(BP)):
                t = []
                for q in (0, 1):
                    x = bypfx + labs[j + q]
                    t.append(x[:23] + "..." if len(x) > 26 else x)
                self.add(di(("%4s", "p%d" % (j + 1)), "  ", t[0], col(34), t[1], col(61), ("%9.0fc", nf[j]),
                            col(70), ("%9.0fc", nr[j])))
            self.p("total: the sum of the pairs, from %s = %s to %s = %s" % (by, labs[0], by, labs[-1]), 0, 6)
        # the notes; those of thin target levels and of missed targets (with the levels listed
        # under them) are left to the raking diagnostics
        thin = missnote = False
        shown = False
        for nt0 in self.res.notes:
            nt = self.vlnote(str(nt0))
            cont = nt[:4] == "    "
            if not cont:
                thin = bool(re.search(r"^(p[0-9]+(, p[0-9]+)*: )?([0-9,]+ )?thin target levels?[,:]", nt))
                missnote = bool(re.search(r"(^|: )(trimming left [0-9,]+ raking margins? unmet|[0-9,]+ raking "
                                          r"margins? not met|raking stopped .* with [0-9,]+ raking margins? "
                                          r"not met)", nt))
            if thin or missnote:
                continue
            if not shown:
                self.add("")
            shown = True
            if cont:
                self.p(nt, 6, 8)
            else:
                self.p("note: " + nt, 0, 6)

    # ---- Display_table ------------------------------------------------------------------
    def table(self, level, cf, pf, sf, ci, cilevel, outcome=1, atest=True, citypes=()):
        """The coefficient table of an outcome, its footnotes and the explained share."""
        hasV = self.has("V")
        vcelab, vcehead, vcebs = self.vce()
        names, b = self.bvec()
        k = len(names)
        V = self.mat("V") if hasV else None
        depvar = self.Vlab(self.es("depvar"), 16)
        outs = self.es("outcomes").split()
        nout = max(1, len(outs))
        j0 = (outcome - 1) * k // nout
        j1 = outcome * k // nout
        P = self.ev("pct_explained")
        PSE = self.ev("pct_explained_se")
        MZ = self.ev("mean_focal_zero")
        ES = self.ev("explained_support")
        ECS = self.ev("explained_cs")
        ylab = self.es("depvar")
        eq0 = "overall"
        yo = None
        pos = {n: i for i, n in enumerate(names)}
        if nout > 1:
            yo = outs[outcome - 1]
            ylab = yo
            eq0 = yo
            depvar = self.Vlab(yo, 16)
            dfn = pos["%s:difference" % yo]
            P = 100 * b[dfn + 1] / b[dfn] if b[dfn] != 0 else NAN
            PSE = NAN
            if self.has("by_outcome"):
                BO = self.mat("by_outcome")
                P = BO[outcome - 1, 6]
                PSE = BO[outcome - 1, 7]
            if outcome > 1 and self.ev("N_focal_zero") > 0:
                S = self.mat("bystage%d" % outcome)
                ns = int(self.ev("N_stages"))
                MZ = S[ns - 1, 28]
                mcs = S[ns - 1, 29]
                ES = b[dfn - 3] - mcs
                ECS = mcs - b[dfn - 2]
        if self.es("ratio"):
            ylab = ylab + "/" + self.es("ratio")
        eps = self.eps(outcome)
        df = self.ev("df_r")
        st = "t" if not math.isnan(df) else "z"
        crit = _inv_t(df, level) if not math.isnan(df) else _inv_normal(level)
        ownbs = bool(self.es("bsnames"))
        mt = "normal"
        ccrit = crit
        CM = None
        if ownbs:
            mt = citypes[0]
            if mt == "normal":
                ccrit = _inv_normal(level)
            elif mt == "t":
                ccrit = _inv_t(self.ev("df_bs"), level)
            else:
                CM = ci["ci_" + mt]
        bsn = self.es("bsnames").split()
        bpos = {n: i for i, n in enumerate(bsn)}
        cw, pw, sw = fmt_width(cf), fmt_width(pf), fmt_width(sf)
        self.add("", TOP)
        if hasV:
            vw = vcehead.split()
            vlast = vw[-1] if vw else ""
            vfirst = " ".join(vw[:-1])
            cilab = "[%s%% %s]" % (_g12(level), H_CI)
            if mt not in ("normal", "t"):
                cilab = "[%s%% %s]" % (cilevel, H_CI)
            cicol = 79 - len(cilab)
            if vfirst:
                self.add(di(" " * 17 + "|", col(32 + (9 - len(vfirst)) // 2), vfirst))
            vcol = 32 + (9 - len(vlast)) // 2
            if vcebs:
                cihd = self.ci_label(mt)
                self.add(di(" " * 17 + "|", col(vcol), vlast, col(69 - len(cihd) // 2), cihd))
            else:
                self.add(di(" " * 17 + "|", col(vcol), vlast))
            self.add(di(ralign(depvar, 16) + " |", col(19), H_COEF, col(32), H_SEL, col(49), st, col(51),
                        "P>|%s|" % st, col(cicol), cilab))
        else:
            self.add(di(ralign(depvar, 16) + " |", col(19), H_COEF))
        self.add(MID)
        preveq = None
        for j in range(j0, j1):
            fn = names[j]
            eq, nm = fn.split(":", 1)
            if nout > 1:
                eq = "overall" if eq == yo else eq[len(yo) + 1:]
            if not hasV and eq != "overall":
                continue
            if eq != preveq:
                if preveq is not None:
                    self.add(MID)
                self.add(di(("%-17s", self.Eq_lab(eq)), "|" + self.Eq_title(eq)))
                preveq = eq
            coef = b[j]
            if abs(coef) < eps:
                coef = 0.0
            line = di(ralign(self.Row_lab(eq, nm), 16) + " |", col(30 - cw), (cf, coef))
            if hasV:
                se = math.sqrt(V[j, j]) if V[j, j] >= 0 else NAN
                if se < eps:
                    se = 0.0
                if not math.isnan(se) and se > 0:
                    z = coef / se
                    if not math.isnan(df):
                        pv = 2 * float(stats.t.sf(abs(z), df))
                    else:
                        pv = 2 * normal(-abs(z))
                    if mt in ("normal", "t"):
                        ll, ul = coef - ccrit * se, coef + ccrit * se
                    else:
                        cc = bpos[fn]
                        ll, ul = CM[0, cc], CM[1, cc]
                    if abs(ll) < eps:
                        ll = 0.0
                    if abs(ul) < eps:
                        ul = 0.0
                    line = di(line, col(41 - cw), (cf, se), col(50 - sw), (sf, z), col(56 - pw), (pf, pv),
                              col(67 - cw), (cf, ll), col(79 - cw), (cf, ul))
            self.add(line)
        self.add(BOT)
        self.add("Difference = focal - reference. Explained = focal - raked focal.")
        if ownbs and hasV and (mt == "t" or not math.isnan(df)):
            dfw = "clusters" if self.es("clustvar") else "observations"
            self.p("(t: %s degrees of freedom, the number of %s less 1)" % (_g12(self.ev("df_bs")), dfw), 0, 1)
        cd = pos.get("%s:difference" % eq0)
        if cd is not None and abs(b[cd]) < eps:
            P = NAN
        if not math.isnan(P):
            sh = self.pct(P)
            if hasV and not math.isnan(PSE):
                pci = ci.get("pct_explained_ci") if ci else None
                self.p("Explained share of the difference: %s%%, bootstrap std. err. %s" % (sh, self.pct(PSE)))
                if pci is not None:
                    self.p("%s%% percentile CI: %s%% to %s%%"
                           % (cilevel, self.pct(pci[outcome - 1, 0]), self.pct(pci[outcome - 1, 1])), 4, 4)
            else:
                self.p("Explained share of the difference: %s%%" % sh)
        else:
            self.add("Explained share of the difference: n/a (the difference is 0)")
        if self.ev("N_focal_zero") > 0:
            mtxt = "ratio" if self.es("ratio") else "mean"
            ES, ECS, MZ = [0.0 if abs(q) < eps else q for q in (ES, ECS, MZ)]
            self.p("explained = %s (change of support: %s focal group obs with %s%% of its weight, %s %s %s, "
                   "get zero weight) + %s (reweighting within the common support)"
                   % (fmt_num(ES, "%9.0g"), _g12(self.ev("N_focal_zero")),
                      fmt_num(self.ev("pct_W_focal_zero"), "%4.1f"), mtxt, ylab, fmt_num(MZ, "%9.0g"),
                      fmt_num(ECS, "%9.0g")))
        if hasV and outcome == 1:
            self.rep_footnotes(atest, chain=False)

    def rep_footnotes(self, atest, chain):
        """At replay, the failed replications, those stopped at the iterate() limit, and those
        force kept (at estimation, the messages before the output say so)."""
        nf = self.ev("N_misreps")
        nf = 0 if math.isnan(nf) else nf
        ni = self.ev("N_iterlimit")
        ni = 0 if math.isnan(ni) else ni
        nr = self.ev("N_reps") + nf
        if not atest and self.es("bsnames") and (nf > 0 or ni > 0):
            fin = (not chain) and self.fold_has(", also in %s of %s replications" % (_g12(nf + ni), _g12(nr)))
            if fin:
                pass
            elif nf > 0:
                self.p("(%s of %s replications failed and %s dropped; bootstrap results may be incorrect)"
                       % (_g12(nf), _g12(nr), "was" if nf == 1 else "were"), 0, 1)
            else:
                self.p("(%s of %s replications stopped before the weights converged)"
                       % (_g12(ni), _g12(nr)), 0, 1)
        nfo = self.ev("N_forced")
        if not atest and nfo > 0 and not math.isnan(nfo):
            self.p("(force kept %s of %s replications that missed raking targets)"
                   % (_g12(nfo), _g12(nr)), 0, 1)

    @staticmethod
    def ci_label(ty):
        return {"t": "t-based", "percentile": "Percentile", "bc": "Bias-corrected",
                "bca": "BCa"}.get(ty, "Normal-based")

    @staticmethod
    def ci_code(ty):
        return {"t": "   (T)", "percentile": "   (P)", "bc": "  (BC)", "bca": " (BCa)"}.get(ty, "   (N)")

    def ci_key(self, types):
        keys = []
        for ty in types:
            if ty == "normal":
                keys.append("N: normal")
            elif ty == "t":
                keys.append("T: t with %s degrees of freedom" % _g12(self.ev("df_bs")))
            elif ty == "percentile":
                keys.append("P: percentile")
            elif ty == "bc":
                keys.append("BC: bias-corrected")
            else:
                keys.append("BCa: bias-corrected and accelerated")
        key = ", ".join(keys)
        if self.ev("ties") == 1 and ("bc" in types or "bca" in types):
            key += " (adjusted for ties)"
        return "Key: " + key

    # ---- Display_stages, Display_contrib ------------------------------------------------
    def pair_title(self, base, j):
        lv = self.es("by_levels").split()
        v1, v2 = lv[j - 1], lv[j]
        t = "%s: p%d, %s %s to %s" % (base, j, self.Vlab(self.es("by")), v1, v2)
        if len(t) > 78:
            t = "%s: p%d, %s to %s" % (base, j, v1, v2)
        return t

    def title(self, t):
        """A table title after a blank line, wrapped beyond 78 columns."""
        self.add("")
        if len(t) > 78:
            self.p(t)
        else:
            self.add(t)

    def stages(self, cf, outcome=1, reverse=False, pair=0):
        """The stage table (sequentially) of an outcome, of the reverse order or of a chain's
        pair; with individually, the table of the block contributions."""
        if self.es("stagetype") == "individually":
            self.contrib(cf, outcome)
            return
        cw = fmt_width(cf)
        ns = int(self.ev("N_stages"))
        mname = "bystage" + ("_rev" if reverse else "") + (str(outcome) if outcome > 1 else "")
        blocks = self.es("blocks").split()
        if reverse:
            blocks = blocks[::-1]
        fom = ", fixothermargins" if self.es("fixothermargins") else ""
        ttl = "Decomposition by raking stage (sequentially%s)" % fom
        if reverse:
            ttl += ", reverse order"
        outs = self.es("outcomes").split()
        if outs:
            ttl += ": " + self.Vlab(outs[outcome - 1], 16)
        if pair > 0:
            mname = "bystage" + (str(pair) if pair > 1 else "")
            ttl = self.pair_title(ttl, pair)
        S = self.mat(mname)
        self.title(ttl)
        self.add(RULE)
        self.add(di(col(33), "-" * 8 + " Cumulative " + "-" * 8, col(62), "-" * 4 + " Marginal " + "-" * 3))
        rlab = "Raked ratio" if self.es("ratio") else "Raked mean"
        self.add(di("Stage Block", col(20), rlab, col(33), "Explained", col(43), "Unexplained", col(60), "%",
                    col(63), "Explained", col(78), "%"))
        self.add(RULE)
        eps = self.eps(pair if pair > 0 else outcome)
        eq0 = outs[outcome - 1] if outs else "overall"
        if pair > 0:
            eq0 = "p%d" % pair
        names, b = self.bvec()
        pos = {n: i for i, n in enumerate(names)}
        dz = abs(b[pos["%s:difference" % eq0]]) < eps
        det = bool(self.es("detail")) and not reverse
        dmulti = False
        if det:
            dmat = "bystage_detail" + (str(outcome) if outcome > 1 else "")
            hasD = self.has(dmat)
            D = self.mat(dmat) if hasD else None
            deq = (outs[outcome - 1] + "_" if outs else "") + "detail"
            dcol, dnm = [], []
            for jj, fn in enumerate(names):
                eq, nm = fn.split(":", 1)
                if eq == deq:
                    dcol.append(jj)
                    dnm.append(nm)
            g1 = b[pos["%s:focal_group" % eq0]]
            dd = b[pos["%s:difference" % eq0]]
        tt = 0
        for s in range(1, ns + 1):
            nm = blocks[s - 1]
            v = {c: (0.0 if abs(S[s - 1, c - 1]) < eps else S[s - 1, c - 1]) for c in (1, 2, 3, 5)}
            p4 = NAN if dz else S[s - 1, 3]
            p27 = NAN if dz else S[s - 1, 26]
            if abs(p4) < .05:
                p4 = 0.0
            if abs(p27) < .05:
                p27 = 0.0
            self.add(di(("%5.0f", s), " ", ("%-12s", self.Block_lab(nm, 12)), col(31 - cw), (cf, v[1]),
                        col(42 - cw), (cf, v[2]), col(54 - cw), (cf, v[3]), col(55), ("%6s", self.pct6(p4)),
                        col(72 - cw), (cf, v[5]), col(73), ("%6s", self.pct6(p27))))
            if not det:
                continue
            nbt = len(self.es("block%d" % s).split())
            cprev = S[s - 2, 1] if s > 1 else 0.0
            cw0 = 0.0
            for u in range(tt, tt + nbt):
                if nbt < 2:
                    continue
                dmulti = True
                inc = b[dcol[u]]
                if hasD:
                    r1 = D[u, 0]
                    cwb = D[u, 1] - cprev
                    q = D[u, 26]
                else:
                    cwb = cw0 + inc
                    r1 = g1 - (cprev + cwb)
                    q = 100 * inc / dd if dd != 0 else NAN
                cw0 = cwb
                r1, cwb, inc = [0.0 if abs(x) < eps else x for x in (r1, cwb, inc)]
                if dz:
                    q = NAN
                if abs(q) < .05:
                    q = 0.0
                self.add(di(col(9), ("%-12s", self.Vlab(dnm[u], 12)), col(31 - cw), (cf, r1), col(42 - cw),
                            (cf, cwb), col(72 - cw), (cf, inc), col(73), ("%6s", self.pct6(q))))
            tt += nbt
        self.add(RULE)
        if outcome > 1:
            return
        if 0 < pair < self.ev("N_pairs"):
            return
        if reverse:
            self.p("Reverse order: blocks added from last to first.")
            return
        self.p("%: percent of the difference.")
        if dmulti:
            self.p("Indented: the block's terms, added in turn (Cumulative: within the block).")

    def contrib(self, cf, outcome=1):
        """individually: each block's contribution raked alone, their sum, all blocks and the
        remainder."""
        cw = fmt_width(cf)
        ns = int(self.ev("N_stages"))
        nb = ns - 1
        blocks = self.es("blocks").split()
        S = self.mat("bystage" + (str(outcome) if outcome > 1 else ""))
        fom = ", fixothermargins" if self.es("fixothermargins") else ""
        ttl = "Contribution of each block, raked alone (individually%s)" % fom
        outs = self.es("outcomes").split()
        if outs:
            ttl += ": " + self.Vlab(outs[outcome - 1], 16)
        self.title(ttl)
        self.add(RULE)
        rlab = "Raked ratio" if self.es("ratio") else "Raked mean"
        self.add(di("Block", col(20), rlab, col(37), "Explained", col(53), "%"))
        self.add(RULE)
        eps = self.eps(outcome)
        eq0 = outs[outcome - 1] if outs else "overall"
        names, b = self.bvec()
        pos = {n: i for i, n in enumerate(names)}
        dz = abs(b[pos["%s:difference" % eq0]]) < eps
        tot = 0.0
        psum = NAN if dz else 0.0

        def pfmt(x):
            return "%6.1f" if (math.isnan(x) or abs(x) < 1000) else ("%6.0f" if abs(x) < 1e5 else "%6.0g")

        pp, fp, v = {}, {}, {}
        for s in range(1, ns + 1):
            p = NAN if dz else (S[s - 1, 26] if s < ns else S[s - 1, 3])
            if abs(p) < .05:
                p = 0.0
            fp[s] = pfmt(p)
            pp[s] = p
            if s < ns:
                tot += S[s - 1, 4]
                if not dz:
                    psum += S[s - 1, 26]
            for c in (1, 2, 5):
                x = S[s - 1, c - 1]
                v[c, s] = 0.0 if abs(x) < eps else x
        if abs(tot) < eps:
            tot = 0.0
        for s in range(1, nb + 1):
            self.add(di(("%-21s", self.Block_lab(blocks[s - 1], 21)), col(31 - cw), (cf, v[1, s]), col(46 - cw),
                        (cf, v[5, s]), col(48), (fp[s], pp[s])))
        self.add(RULE)
        pr = NAN if dz else S[ns - 1, 26]
        if abs(psum) < .05:
            psum = 0.0
        if abs(pr) < .05:
            pr = 0.0
        self.add(di("Sum of the blocks", col(46 - cw), (cf, tot), col(48), (pfmt(psum), psum)))
        self.add(di("All blocks (overall)", col(31 - cw), (cf, v[1, ns]), col(46 - cw), (cf, v[2, ns]), col(48),
                    (fp[ns], pp[ns])))
        self.add(di("Remainder", col(46 - cw), (cf, v[5, ns]), col(48), (pfmt(pr), pr)))
        self.add(RULE)
        if outcome > 1:
            return
        how = "are unconstrained"
        if self.es("fixothermargins"):
            how = "are held at their focal group distributions"
        self.p("%%: percent of the difference. Other blocks %s. Remainder: explained part of all blocks less "
               "the sum of the blocks' parts." % how)

    # ---- the raking diagnostics ---------------------------------------------------------
    def diag_sets(self):
        sets = []
        if self.has("bystage"):
            sets.append("base")
        if self.es("chain"):
            for j in range(2, int(self.ev("N_pairs")) + 1):
                if self.has("bystage%d" % j):
                    sets.append(str(j))
        else:
            for q in ("rev", "detail"):
                if self.has("bystage_" + q):
                    sets.append("_" + q)
        return sets

    def diag_get(self, ds):
        sx = "" if ds == "base" else ds
        S = self.mat("bystage" + sx)
        saved = self.has("raking_checks" + sx)
        if saved:
            C = self.mat("raking_checks" + sx)
        else:
            C = np.full((S.shape[0], 5), NAN)
            for j in range(S.shape[0]):
                C[j, 0] = S[j, 5]
                if _big(S[j, 6]) > 0 and S[j, 5] == 0:
                    C[j, 0] = NAN
                C[j, 1] = S[j, 6]
        ctxt = caption = ""
        lp = "raking_miss" + sx
        if ds == "_rev":
            ctxt, caption = "reverse order", "Reverse order"
        elif ds == "_detail":
            ctxt, caption = "detail", "Term-by-term raking"
        elif self.es("chain"):
            j = 1 if ds == "base" else int(ds)
            ctxt, caption = "p%d" % j, "Pair p%d" % j
            if j > 1:
                lp = "raking_miss%d_" % j
        elif self.has("bystage_rev"):
            ctxt, caption = "forward order", "Forward order"
        return dict(sx=sx, S=S, C=C, context=ctxt, caption=caption, labelprefix=lp, saved=saved)

    def diag_scan(self):
        issue = False
        for ds in self.diag_sets():
            C = self.diag_get(ds)["C"]
            for j in range(C.shape[0]):
                if _big(C[j, 0]) != 1 or _big(C[j, 1]) > 0:
                    issue = True
        return issue

    def raking_status(self):
        """The one-line raking status: convergence and targets, the stages that failed."""
        sets = self.diag_sets()
        if not sets:
            return
        baditer = badtarget = unknown = multi = False
        iw, tw = "", ""
        for ds in sets:
            g = self.diag_get(ds)
            C = g["C"]
            ns = C.shape[0]
            if ns > 1 or ds != "base":
                multi = True
            il, tl = [], []
            for j in range(ns):
                if _big(C[j, 0]) < 1:
                    baditer = True
                    il.append(j + 1)
                if _big(C[j, 0]) == math.inf:
                    unknown = True
                if _big(C[j, 1]) > 0:
                    badtarget = True
                    tl.append(j + 1)
            for k, lst in (("i", il), ("t", tl)):
                if not lst:
                    continue
                loc = ""
                if ns > 1:
                    loc = ("stage " if len(lst) == 1 else "stages ") + ", ".join(str(x) for x in lst)
                if g["context"]:
                    loc = g["context"] + (": " + loc if loc else "")
                if loc:
                    if k == "i":
                        iw = iw + ("; " if iw else "") + loc
                    else:
                        tw = tw + ("; " if tw else "") + loc
        allw = "all pairs and stages" if self.es("chain") else "all stages"
        if not baditer and not badtarget and not unknown:
            pfx = "Raking (%s)" % allw if multi else "Raking"
            self.p("%s: weight iterations converged; all targets met." % pfx)
            return
        if baditer:
            it = "Convergence failure"
        elif unknown:
            it = "Iteration convergence unavailable"
        else:
            it = "all weight iterations converged" if multi else "weight iterations converged"
        if baditer and iw:
            it += " (%s)" % iw
        tt = "targets *not* met" if badtarget else "all targets met"
        if badtarget and tw:
            tt += " (%s)" % tw
        self.p("Raking: %s; %s." % (it, tt))

    def raking_footer(self, shown, issue):
        """The raking status after the tables and, without the diagnostics, the range of the
        raking factors and the thin target levels."""
        sets = self.diag_sets()
        if not sets:
            return
        if not shown or issue:
            self.add("")
            self.raking_status()
        if shown:
            return
        lo, hi = math.inf, 0.0
        thin = ""
        nthin = 0.0
        multi = False
        for ds in sets:
            g = self.diag_get(ds)
            S = g["S"]
            ns = S.shape[0]
            if ns > 1 or ds != "base":
                multi = True
            for j in range(ns):
                lo = min(lo, _big(S[j, 15]))
                hi = max(hi, _big(S[j, 16]))
            if ds == "base" or self.es("chain"):
                if not math.isnan(S[ns - 1, 25]):
                    nthin += S[ns - 1, 25]
                n = fmt_str(S[ns - 1, 25], "%12.0fc")
                if self.es("chain"):
                    thin = thin + ("; " if thin else "") + "%s = %s" % (g["context"], n)
                else:
                    thin = n
        lx = fmt_str(lo if lo != math.inf else NAN, "%9.4g")
        hx = fmt_str(hi if hi != math.inf else NAN, "%9.4g")
        allw = "all pairs and stages" if self.es("chain") else "all stages"
        span = " (%s)" % allw if multi else ""
        self.p("Raking factors%s: %s-%s." % (span, lx, hx))
        if nthin > 0:
            self.p("Thin target levels (<%s focal obs): %s." % (_g12(self.ev("mincell")), thin))

    def diag(self):
        """The raking diagnostics of each set of stages: raking factors and effective sizes, and
        when a check fails the failure counts and the worst missed targets."""
        sets = self.diag_sets()
        if not sets:
            self.add("", "(raking diagnostics are unavailable for these saved results)")
            return
        anyissue = self.diag_scan()
        mincell = _g12(self.ev("mincell"))
        self.add("", "Raking diagnostics")
        for ds in sets:
            g = self.diag_get(ds)
            S, C = g["S"], g["C"]
            ns = S.shape[0]
            staged = ns > 1
            dx = 0 if staged else -5
            issue = any(_big(C[j, 0]) != 1 or _big(C[j, 1]) > 0 for j in range(ns))
            if g["caption"]:
                self.add("", g["caption"])
            self.add(RULE)
            self.add(di(col(9 + dx), "-- Raking factor --", col(35 + dx), "Kish effective N", col(55 + dx),
                        "ESS"))
            head = "Stage" if staged else ""
            self.add(di(head, col(12 + dx), "Min", col(21 + dx), "Max", col(30 + dx), "CV", col(38 + dx),
                        "Base", col(47 + dx), "Raked", col(55 + dx), "ratio", col(61 + dx), "Zero wt",
                        col(70 + dx), "At bound"))
            self.add(RULE)
            for j in range(ns):
                cv = S[j, 18]
                if abs(cv) < 1e-13:
                    cv = 0.0
                nf = "%10.1fc" if _big(S[j, 22]) < 999999.5 else "%10.0fc"
                de = S[j, 24]
                dfm = "%7.3f" if _big(de) < 100 else "%7.0g"
                sl = fmt_str(j + 1, "%5.0f") if staged else ""
                self.add(di(sl, col(6 + dx), ("%9.4g", S[j, 15]), ("%9.4g", S[j, 16]), ("%8.4g", cv),
                            (nf, S[j, 22]), (nf, S[j, 23]), (dfm, de), ("%9.0fc", S[j, 20]),
                            ("%9.0fc", S[j, 21])))
            self.add(RULE)
            self.p("Raking factor = raked / base weight. ESS ratio = base / raked effective N.")
            self.p("Zero wt: zero raked weights. At bound: distinct weights at a trimming bound.")
            self.p("Thin target levels (<%s focal obs): %s." % (mincell, fmt_str(S[ns - 1, 25], "%12.0fc")))
            if not issue:
                continue
            self.add("", "Raking failure counts")
            self.add(RULE)
            self.add(di(head, col(8 + dx), "Converged", col(21 + dx), "Targets met", col(43 + dx), "# Missed",
                        col(60 + dx), "# Missed"))
            self.add(di(col(42 + dx), "Variables", col(60 + dx), "Targets"))
            self.add(RULE)
            for j in range(ns):
                c1 = C[j, 0]
                cv = "." if math.isnan(c1) else ("yes" if c1 == 1 else "no")
                tm = "yes" if C[j, 1] == 0 else "no"
                sl = fmt_str(j + 1, "%5.0f") if staged else ""
                self.add(di(sl, col(8 + dx), cv, col(21 + dx), tm, col(42 + dx), ("%9.0f", C[j, 1]),
                            col(58 + dx), ("%9.0f", C[j, 2])))
            self.add(RULE)
            self.p("Converged: weight iterations. Targets met: within ctrltolerance().")
            self.p("Missed variables: raking variables with one or more missed level targets.")
            self.p("Missed targets: individual level totals outside ctrltolerance().")
            if not g["saved"]:
                self.p("Cell-level checks are unavailable for these saved results; re-estimate to display "
                       "them.")
            # without a matrix of missed targets (a stage that failed by convergence alone), the
            # table holds one row of missing values: Stata's matrix of an absent e() result is a
            # missing scalar
            mname = "raking_misses" + g["sx"]
            Mm = self.mat(mname) if self.has(mname) else np.full((1, 6), NAN)
            self.add("", "Worst missed target for each affected variable")
            self.add(RULE)
            self.add(di(head, col(7 + dx), "Variable", col(20 + dx), "Level", col(44 + dx), "Target",
                        col(54 + dx), "Achieved", col(66 + dx), "Discrepancy"))
            self.add(RULE)
            for j in range(Mm.shape[0]):
                mv = self.es("%s%d_var" % (g["labelprefix"], j + 1))
                ml = self.es("%s%d_level" % (g["labelprefix"], j + 1))
                sl = fmt_str(Mm[j, 0], "%5.0f") if staged else ""
                self.add(di(sl, col(7 + dx), self.Vlab(mv, 12), col(20 + dx), vltrim(ml, 18), col(39 + dx),
                            ("%10.4g", Mm[j, 3]), col(51 + dx), ("%11.4g", Mm[j, 4]), col(65 + dx),
                            ("%11.4g", Mm[j, 5])))
            self.add(RULE)
            ctl = 1e-6
            m = re.search(r"ctrltolerance\(([^)]+)\)", self.es("rakeopts"))
            if m:
                ctl = _num(m.group(1))
            self.p("Discrepancy = |s*achieved - s*target| / (|s*target| + 1), where s scales base weights to "
                   "the raking sample size.")
            self.p("Targets are met when discrepancy <= %s." % fmt_str(ctl, "%9.0g"))
        if not anyissue:
            self.add("")
            self.raking_status()

    # ---- Display_bootci -----------------------------------------------------------------
    def bootci(self, level, cf, ci, cilevel, outcome=1, citypes=()):
        """The bootstrap interval table of an outcome: a row per interval type shown."""
        cw = fmt_width(cf)
        names, b = self.bvec()
        pos = {n: i for i, n in enumerate(names)}
        bsn = self.es("bsnames").split()
        kb = len(bsn)
        depvar = self.Vlab(self.es("depvar"), 16)
        outs = self.es("outcomes").split()
        nout = max(1, len(outs))
        j0 = (outcome - 1) * kb // nout
        j1 = outcome * kb // nout
        yo = None
        if nout > 1:
            yo = outs[outcome - 1]
            depvar = self.Vlab(yo, 16)
        eps = self.eps(outcome)
        lv = cilevel
        self.add("", di("Bayesian bootstrap confidence intervals", col(49), "Replications", col(67), "= ",
                        ("%10.0fc", self.ev("N_reps"))))
        self.boot_head(depvar, lv)
        bias = self.mat("bias")[0]
        se = self.mat("se")[0]
        preveq = None
        for j in range(j0, j1):
            fn = bsn[j]
            eq, nm = fn.split(":", 1)
            if nout > 1:
                eq = "overall" if eq == yo else eq[len(yo) + 1:]
            if eq != preveq:
                if preveq is not None:
                    self.add(MID)
                self.add(di(("%-17s", self.Eq_lab(eq)), "|" + self.Eq_title(eq)))
                preveq = eq
            self.boot_rows(self.Row_lab(eq, nm), b[pos[fn]], bias[j], se[j], j, ci, citypes, cw, cf, eps)
        self.add(BOT)
        self.p(self.ci_key(citypes))
        if cilevel != _g12(level) and outcome == 1:
            self.p("(intervals at the estimation level, %s%%: level() cannot change them for these "
                   "results)" % lv, 0, 1)

    def boot_head(self, depvar, lv):
        self.add(TOP)
        self.add(di(("%16s", ""), " |", col(42), "Bootstrap"))
        self.add(di(ralign(depvar, 16) + " |", col(19), H_COEF, col(37), "Bias", col(42), H_SEL, col(54),
                    "[%s%% %s]" % (lv, H_CI)))
        self.add(MID)

    def boot_rows(self, label, x0, x1, x2, j, ci, citypes, cw, cf, eps):
        """The rows of an element in the bootstrap interval table: the estimate, bias, std.
        err. and the first interval type, then a row for each other type."""
        x = [x0, x1, x2]
        for M in citypes:
            x += [ci["ci_" + M][0, j], ci["ci_" + M][1, j]]
        x = [0.0 if abs(v) < eps else v for v in x]
        c1, c2, c3, c4, c5 = 30 - cw, 41 - cw, 51 - cw, 62 - cw, 73 - cw
        for r, M in enumerate(citypes):
            i1, i2 = 3 + 2 * r, 4 + 2 * r
            cc = self.ci_code(M)
            if r == 0:
                self.add(di(ralign(label, 16) + " |", col(c1), (cf, x[0]), col(c2), (cf, x[1]), col(c3),
                            (cf, x[2]), col(c4), (cf, x[i1]), col(c5), (cf, x[i2]), cc))
            else:
                self.add(di(("%16s", ""), " |", col(c4), (cf, x[i1]), col(c5), (cf, x[i2]), cc))

    # ---- Display_chain, Display_chain_share, Pair_title ----------------------------------
    def chain(self, level, cf, pf, sf, ci, cilevel, atest=True, bootci=False, citypes=()):
        """chain: the coefficient table (bootci: the bootstrap interval table) in one piece, the
        pairs and the total, each with its equations; then the explained shares."""
        hasV = self.has("V")
        vcelab, vcehead, vcebs = self.vce()
        names, b = self.bvec()
        pos = {n: i for i, n in enumerate(names)}
        V = self.mat("V") if hasV else None
        depvar = self.Vlab(self.es("depvar").split()[0], 16)
        J = int(self.ev("N_pairs"))
        df = self.ev("df_r")
        st = "t" if not math.isnan(df) else "z"
        crit = _inv_t(df, level) if not math.isnan(df) else _inv_normal(level)
        ownbs = bool(self.es("bsnames"))
        mt = "normal"
        ccrit = crit
        CM = None
        if ownbs:
            mt = citypes[0]
            if mt == "normal":
                ccrit = _inv_normal(level)
            elif mt == "t":
                ccrit = _inv_t(self.ev("df_bs"), level)
            else:
                CM = ci["ci_" + mt]
        bsn = self.es("bsnames").split()
        bpos = {n: i for i, n in enumerate(bsn)}
        cw, pw, sw = fmt_width(cf), fmt_width(pf), fmt_width(sf)
        lv = cilevel
        if bootci:
            cols = bsn
            self.add("", di("Bayesian bootstrap confidence intervals", col(49), "Replications", col(67), "= ",
                            ("%10.0fc", self.ev("N_reps"))))
            self.add(TOP)
            self.add(di(("%16s", ""), " |", col(42), "Bootstrap"))
            self.add(di(ralign(depvar, 16) + " |", col(19), H_COEF, col(37), "Bias", col(42), H_SEL, col(54),
                        "[%s%% %s]" % (lv, H_CI)))
        else:
            cols = names
            self.add("", TOP)
            if hasV:
                vw = vcehead.split()
                vlast = vw[-1] if vw else ""
                vfirst = " ".join(vw[:-1])
                cilab = "[%s%% %s]" % (_g12(level), H_CI)
                if mt not in ("normal", "t"):
                    cilab = "[%s%% %s]" % (cilevel, H_CI)
                cicol = 79 - len(cilab)
                if vfirst:
                    self.add(di(" " * 17 + "|", col(32 + (9 - len(vfirst)) // 2), vfirst))
                vcol = 32 + (9 - len(vlast)) // 2
                if vcebs:
                    cihd = self.ci_label(mt)
                    self.add(di(" " * 17 + "|", col(vcol), vlast, col(69 - len(cihd) // 2), cihd))
                else:
                    self.add(di(" " * 17 + "|", col(vcol), vlast))
                self.add(di(ralign(depvar, 16) + " |", col(19), H_COEF, col(32), H_SEL, col(49), st, col(51),
                            "P>|%s|" % st, col(cicol), cilab))
            else:
                self.add(di(ralign(depvar, 16) + " |", col(19), H_COEF))
        self.add(MID)
        labs = [x.strip() for x in self.es("by_labels").split("|")]
        prevpan = preveq = None
        eps = 0.0
        bias = self.mat("bias")[0] if bootci else None
        se_bs = self.mat("se")[0] if bootci else None
        for jb, fn in enumerate(cols):
            eq, nm = fn.split(":", 1)
            m = re.match(r"^(p[0-9]+|total)", eq)
            if not m:
                continue
            pan = m.group(1)
            m2 = re.match(r"^(p[0-9]+|total)_(.+)$", eq)
            sub = m2.group(2) if m2 else "overall"
            if sub in ("cum_expl", "cum_unexpl", "raked_focal") and not hasV:
                continue
            if pan != prevpan:
                if prevpan is not None:
                    self.add(MID)
                pi = 0 if pan == "total" else int(pan[1:])
                k1 = 1 if pi == 0 else pi
                k2 = J + 1 if pi == 0 else pi + 1
                ptx = "%s %s to %s%s" % (self.vlab(self.es("by")), labs[k1 - 1], labs[k2 - 1],
                                         ", the sum of the pairs" if pi == 0 else "")
                if len(ptx) > 59:
                    ptx = ptx[:56] + "..."
                self.add(di(("%-17s", pan), "| " + ptx))
                eps = self.eps(pi)
                prevpan = pan
                preveq = None
            if sub != preveq:
                self.add(di("  ", ("%-15s", self.Eq_lab(sub)), "|" + self.Eq_title(sub)))
                preveq = sub
            j = pos[fn]
            coef = b[j]
            if abs(coef) < eps:
                coef = 0.0
            if bootci:
                self.boot_rows(self.Row_lab(sub, nm), coef, bias[jb], se_bs[jb], jb, ci, citypes, cw, cf, eps)
                continue
            line = di(ralign(self.Row_lab(sub, nm), 16) + " |", col(30 - cw), (cf, coef))
            if hasV:
                se = math.sqrt(V[j, j]) if V[j, j] >= 0 else NAN
                if se < eps:
                    se = 0.0
                if not math.isnan(se) and se > 0:
                    z = coef / se
                    if not math.isnan(df):
                        pv = 2 * float(stats.t.sf(abs(z), df))
                    else:
                        pv = 2 * normal(-abs(z))
                    if mt in ("normal", "t"):
                        ll, ul = coef - ccrit * se, coef + ccrit * se
                    else:
                        cc = bpos[fn]
                        ll, ul = CM[0, cc], CM[1, cc]
                    if abs(ll) < eps:
                        ll = 0.0
                    if abs(ul) < eps:
                        ul = 0.0
                    line = di(line, col(41 - cw), (cf, se), col(50 - sw), (sf, z), col(56 - pw), (pf, pv),
                              col(67 - cw), (cf, ll), col(79 - cw), (cf, ul))
            self.add(line)
        self.add(BOT)
        if not bootci:
            self.add("Difference = focal - reference. Explained = focal - raked focal.")
        if not bootci and ownbs and hasV and (mt == "t" or not math.isnan(df)):
            dfw = "clusters" if self.es("clustvar") else "observations"
            self.p("(t: %s degrees of freedom, the number of %s less 1)" % (_g12(self.ev("df_bs")), dfw), 0, 1)
        if bootci:
            self.p(self.ci_key(citypes))
            if cilevel != _g12(level):
                self.p("(intervals at the estimation level, %s%%: level() cannot change them for these results)"
                       % lv, 0, 1)
            return
        self.chain_share(ci, cilevel)
        if hasV:
            self.rep_footnotes(atest, chain=True)

    def chain_share(self, ci, cilevel):
        """chain: the explained share of each pair and of the total."""
        names, b = self.bvec()
        pos = {n: i for i, n in enumerate(names)}
        BP = self.M["by_pair"]
        bp = BP.to_numpy(dtype=float)
        J = int(self.ev("N_pairs"))
        own = bool(self.es("bsnames"))
        pci = ci.get("pct_explained_ci") if (own and ci) else None
        cse = list(BP.columns).index("pct_explained_se") if "pct_explained_se" in BP.columns else None
        self.add("", "Explained share of the difference (%)", TOP)
        if own:
            self.add(di(("%16s", ""), " |", col(25), "Share", col(32), "Std. err.", col(47),
                        "[%s%% percentile CI]" % cilevel))
        else:
            self.add(di(("%16s", ""), " |", col(25), "Share"))
        self.add(MID)
        for j in range(1, J + 2):
            eq = "p%d" % j if j <= J else "total"
            eps = self.eps(j if j <= J else 0)
            dd = b[pos["%s:difference" % eq]]
            sh = bp[j - 1, 6] if j <= J else self.ev("pct_explained")
            if abs(dd) < eps:
                sh = NAN
            se = NAN
            if own:
                se = (bp[j - 1, cse] if cse is not None else NAN) if j <= J else self.ev("pct_explained_se")
            if math.isnan(sh):
                self.add(di(("%16s", eq), " |", col(21), "n/a (the difference is 0)"))
                continue
            sh, se = [0.0 if abs(q) < .05 else q for q in (sh, se)]
            line = di(("%16s", eq), " |", col(21), ("%9.1f", sh))
            if own:
                line = di(line, col(32), ("%9.1f", se))
                if pci is not None:
                    pl, pu = [0.0 if abs(q) < .05 else q for q in (pci[j - 1, 0], pci[j - 1, 1])]
                    line = di(line, col(46), ("%9.1f", pl), col(57), ("%9.1f", pu))
            self.add(line)
        self.add(BOT)

    # ---- Display_covdetail --------------------------------------------------------------
    def covdetail(self):
        """covariatedetail: the observations and weighted shares of each level of each term."""
        Cdf = self.M["covariates"]
        C = Cdf.to_numpy(dtype=float)
        rn = list(Cdf.index)
        self.add("", "Raking variables: weighted shares by level", RULE)
        self.add(di(col(30), "-" * 5 + " Obs. " + "-" * 4, col(47), "-" * 9 + " Weighted share " + "-" * 7))
        self.add(di("Level", col(32), "Focal", col(41), "Ref.", col(49), "Focal", col(58), "Ref.", col(65),
                    "Raked", col(76), "Gap"))
        self.add(RULE)
        prevt = None
        labs = []
        for r in range(len(rn)):
            t, l = rn[r].split(":", 1)
            if t != prevt:
                term = self.es("term%s" % t)
                self.add(self.Vlab(term))
                prevt = t
                tvn = term.replace("#", " ").split()
                labs = [x.strip() for x in self.es("levels%s" % t).split("|")]
                if len(tvn) > 1:
                    labs = [ixlab(x, tvn) for x in labs]
            gap = C[r, 0] - C[r, 1]
            if abs(gap) < .00005:
                gap = 0.0
            lab = labs[int(l) - 1]
            nums = (col(30), ("%7.0fc", C[r, 3]), col(38), ("%7.0fc", C[r, 4]), col(47), ("%7.4f", C[r, 0]),
                    col(55), ("%7.4f", C[r, 1]), col(63), ("%7.4f", C[r, 2]), col(71), ("%8.4f", gap))
            if len(lab) > 26:
                self.p(lab, 2, 4)
                self.add(di(*nums))
            else:
                self.add(di("  " + lab, *nums))
        self.add(RULE)
        self.p("Focal: focal group. Ref.: reference group. Raked: focal group share under the final raked "
               "weights. Gap: Focal minus Reference share (the order of the difference).")

    # ---- Display ------------------------------------------------------------------------
    def display(self, level, diagnostics=False, covariatedetail=False, cformat=None, pformat=None,
                sformat=None, citable=False, replay=False):
        """The whole display; returns the intervals shown and their level."""
        cf, pf, sf = self.fmts(cformat, pformat, sformat)
        atest = not replay
        ownbs = bool(self.es("bsnames"))
        ci, cilevel = {}, _g12(level)
        if ownbs:
            ci, cilevel = self.display_ci(level)
        citypes = []
        if ownbs:
            citypes = self.es("citypes").split() or ["normal"]
            if citable:
                add = ["percentile", "bc"]
                if not [t for t in citypes if t in ("normal", "t")]:
                    add.append("t" if self.es("clustvar") else "normal")
                if self.ev("bca") == 1:
                    add.append("bca")
                citypes = citypes + [a for a in add if a not in citypes]
        citab = len(citypes) > 1 and not (replay and not citable)
        chain = bool(self.es("chain"))
        self.header()
        if ownbs:
            self.jkfail()
        issue = self.diag_scan()
        showdiag = bool(diagnostics) or issue
        if showdiag:
            self.sep()
            self.diag()
        if covariatedetail:
            self.sep()
            if not self.has("covariates"):
                self.add("", "(covariatedetail was not specified at estimation; nothing to display)")
            else:
                self.covdetail()
        self.sep()
        ns = int(self.ev("N_stages"))
        detail = bool(self.es("detail"))
        multi_terms = self.ev("N_terms") > self.ev("N_blocks")
        stg = ns > 1 or (not chain and detail and multi_terms)
        if not stg:
            self.add("", "Decomposition")
        if chain:
            if ns > 1:
                for p in range(1, int(self.ev("N_pairs")) + 1):
                    self.stages(cf, pair=p)
            self.chain(level, cf, pf, sf, ci, cilevel, atest=atest, citypes=citypes)
            if citab:
                self.chain(level, cf, pf, sf, ci, cilevel, bootci=True, citypes=citypes)
        else:
            nout = max(1, len(self.es("outcomes").split()))
            for o in range(1, nout + 1):
                if ns > 1 or (detail and multi_terms):
                    self.stages(cf, outcome=o)
                if ns > 1 and self.es("reverse"):
                    self.stages(cf, outcome=o, reverse=True)
                self.table(level, cf, pf, sf, ci, cilevel, outcome=o, atest=atest, citypes=citypes)
                if citab:
                    self.bootci(level, cf, ci, cilevel, outcome=o, citypes=citypes)
        if citable and not ownbs:
            if not self.has("V"):
                self.add("", "(citable: no confidence intervals without bayesboot() or a prefix)")
            else:
                self.add("", "(citable: no bayesboot() replicates; the CIs above are %s)" % self.vce()[0])
        self.raking_footer(showdiag, issue)
        return ci, cilevel

    # ---- Display_estat ------------------------------------------------------------------
    def estat(self, all=False, normal=False, percentile=False, bc=False, bca=False, nolegend=False,
              verbose=False):
        """estat bootstrap: every element of e(bsnames) with the intervals asked for, at the level of
        the last display."""
        if all:
            normal = percentile = bc = True
            if self.ev("bca") == 1:
                bca = True
        if not (normal or percentile or bc or bca):
            bc = True
        kinds = [k for k, on in (("normal", normal), ("percentile", percentile), ("bc", bc), ("bca", bca))
                 if on]
        lv = self.ev("level")
        rl = (getattr(self.res, "_replay", None) or {}).get("level")
        if rl is not None:
            lv = rl
        ci, cilevel = self.display_ci(lv)
        cf = self.fmts(check=True)[0]
        cw = fmt_width(cf)
        names, b = self.bvec()
        pos = {n: i for i, n in enumerate(names)}
        bsn = self.es("bsnames").split()
        kb = len(bsn)
        outs = self.es("outcomes").split()
        nout = max(1, len(outs))
        self.add("", di(self.es("title"), col(49), "Number of obs", col(67), "= ", ("%10.0fc", self.ev("N"))))
        self.add(di(col(49), "Replications", col(67), "= ", ("%10.0f", self.ev("N_reps"))))
        if not nolegend:
            self.add("")
            if verbose:
                self.out.extend(_p2col("Command:", self.es("cmdline")))
        self.add("")
        if self.es("clustvar") and not math.isnan(self.ev("N_clust")):
            self.add(ralign("(Replications based on %s clusters in %s)"
                            % (fmt_str(self.ev("N_clust"), "%12.0fc"), self.es("clustvar")), 78))
        dv = self.es("depvar").split()
        depvar = abbrev(dv[0], 16) if len(dv) == 1 else ""
        lvs = fmt_str(stata_round(float(cilevel), .01), "%9.0g")
        if len(lvs) > 2:
            sp = " " if len(lvs) == 4 else ""
            cilab = "   [%s%s%% %s]" % (sp, lvs, H_CI.replace("nterval", "nt.", 1))
        else:
            cilab = "  [%s%% %s]" % (lvs, H_CI)
        self.add(TOP)
        self.add(di(col(18), "|", col(42), "Bootstrap"))
        self.add(di(ralign(depvar, 16) + " |", col(30 - len(H_COEF)), H_COEF, col(37), "Bias", col(42), H_SEL,
                    col(51), cilab))
        c1, c2, c3, c4, c5 = 30 - cw, 41 - cw, 51 - cw, 62 - cw, 73 - cw
        nper = kb // nout
        bias = self.mat("bias")[0]
        se = self.mat("se")[0]
        preveq = prevpan = None
        prevo = -1
        eps = 0.0
        ch = bool(self.es("chain"))
        for j in range(kb):
            fn = bsn[j]
            eq, nm = fn.split(":", 1)
            o = j // nper + 1
            if ch:
                m = re.match(r"^p([0-9]+)", eq)
                o = int(m.group(1)) if m else 0
            if o != prevo:
                eps = self.eps(o)
                prevo = o
            deq = eq
            pan = None
            if ch:
                m = re.match(r"^(p[0-9]+|total)", eq)
                pan = m.group(1) if m else ""
                m2 = re.match(r"^(p[0-9]+|total)_(.+)$", eq)
                deq = m2.group(2) if m2 else "overall"
            elif nout > 1:
                pan = outs[o - 1]
                deq = "overall" if eq == pan else eq[len(pan) + 1:]
            if eq != preveq:
                self.add(MID)
                if pan != prevpan:
                    if ch:
                        self.add(di(("%-17s", pan), "|"))
                    else:
                        self.add(di(("%-17s", "Outcome"), "| " + pan))
                    prevpan = pan
                self.add(di(("%-17s", self.Eq_lab(deq)), "|" + self.Eq_title(deq)))
                preveq = eq
            x0, x1, x2 = [0.0 if abs(v) < eps else v for v in (b[pos[fn]], bias[j], se[j])]
            first = True
            for k in kinds:
                key = {"normal": "N", "percentile": "P", "bc": "BC"}.get(k, "BCa")
                ll, ul = ci["ci_" + k][0, j], ci["ci_" + k][1, j]
                if abs(ll) < eps:
                    ll = 0.0
                if abs(ul) < eps:
                    ul = 0.0
                if first:
                    self.add(di(ralign(self.Row_lab(deq, nm), 16) + " |", col(c1), (cf, x0), col(c2), (cf, x1),
                                col(c3), (cf, x2), col(c4), (cf, ll), col(c5), (cf, ul), col(73),
                                ("%6s", "(%s)" % key)))
                else:
                    self.add(di(("%16s", ""), " |", col(c4), (cf, ll), col(c5), (cf, ul), col(73),
                                ("%6s", "(%s)" % key)))
                first = False
        self.add(BOT)
        if not nolegend:
            self.add("Difference = focal - reference. Explained = focal - raked focal.")
        ties = " (adjusted for ties) " if self.ev("ties") == 1 else ""
        cp = 9 if bca else (8 if bc else 7)
        keys = []
        if normal:
            keys.append((cp - 1, "N: Normal"))
        if percentile:
            keys.append((cp - 1, "P: Percentile"))
        if bc:
            keys.append((cp - 2, "BC: Bias-corrected" + ties))
        if bca:
            keys.append((cp - 3, "BCa: Bias-corrected and accelerated" + ties))
        for i, (c, t) in enumerate(keys):
            self.add(di("Key:" if i == 0 else "", col(c), t))
        nf = self.ev("N_misreps")
        if nf > 0 and not math.isnan(nf):
            nr = self.ev("N_reps") + nf
            self.p("(%s of %s replications failed and %s dropped; bootstrap results may be incorrect)"
                   % (_g12(nf), _g12(nr), "was" if nf == 1 else "were"), 0, 1)
        self.jkfail()
        return ci, cilevel


def _p2col(label, text):
    """{p2colset 7 16 20 2}{p2col :label}text{p_end}: the label from column 7, the text from
    column 16, its continuation lines from column 20, ending by column 78 (80 columns less the
    right margin of 2)."""
    lines = para(text, 15, 19, 2, 80)
    return [(" " * 6 + label).ljust(15) + lines[0][15:]] + lines[1:]


# ---- the estimation's messages ----------------------------------------------------------

def _message_lines(msg):
    """A message of the estimation as printed: warnings and the detail note as paragraphs,
    the rest (the raking trace, short notes) as they are."""
    if msg.startswith("warning: "):
        return para(msg, 0, 4)
    if msg.startswith("note: term ") or msg.startswith("note: block "):
        return para(msg, 0, 6)
    return [msg]


def _boot_saved(res):
    """bayesboot(saving()): the file of the replicates (the progress dots of the replicate loop
    are printed by the loop itself, as it runs)."""
    if res.e.get("bsfile"):
        return ["(replicate results saved in %s)" % res.e["bsfile"]]
    return []


def _boot_notes(res):
    """The notes of bayesboot() after its replicate loop: the failed replications (what
    failed them), those that stopped at the iterate() limit and were kept, those force kept."""
    info = getattr(res, "_bsinfo", None)
    if not info:
        return []
    e = res.e
    chain = bool(e.get("chain"))
    reps = int(info.get("reps", 0))
    nfail = int(info.get("N_fail", 0))
    niter = int(info.get("N_iter", 0))
    nfailit = int((info.get("nf") or {}).get("it", 0))
    nforce = int(info.get("N_force", 0))
    njforce = int(info.get("N_jkforce", 0))
    why = info.get("why", "")
    stoptxt = "at the iterate() limit"
    nstop = nfailit + niter
    out = []
    folded = (not chain) and nstop > 0 and any(
        (", also in %d of %d replications" % (nstop, reps)) in str(nt) for nt in res.notes)
    if folded:
        pass
    elif nfail > 0 and nfail == nfailit and niter > 0:
        out += para("note: raking stopped %s in %d of %d replications; %d of them missed their targets and %s "
                    "dropped; raise iterate(); bootstrap results may be incorrect."
                    % (stoptxt, nstop, reps, nfail, "was" if nfail == 1 else "were"), 0, 6)
    elif nfail > 0:
        out += para("note: %d of %d replications failed and %s dropped: in %s, %s; bootstrap results may be "
                    "incorrect." % (nfail, reps, "was" if nfail == 1 else "were",
                                    "it" if nfail == 1 else "them", why), 0, 6)
    elif niter > 0:
        out += para("note: %d of %d replications stopped %s with the raking targets met and %s kept; raise "
                    "iterate()" % (niter, reps, stoptxt, "was" if niter == 1 else "were"), 0, 6)
    if nforce + njforce > 0:
        jkt = (" and %d jackknife %s" % (njforce, plural(njforce, "deletion"))) if njforce > 0 else ""
        out += para("note: force kept %d of %d %s%s that missed raking targets"
                    % (nforce, reps, plural(reps, "replication"), jkt), 0, 6)
    return out


def _estimation_options(res):
    """The display options of the estimation command (Estimate's dispopts)."""
    ctx = res._ctx
    o = ctx.opt
    bs = getattr(ctx, "bs", None) or {}
    chain = bool(res.e.get("chain"))
    return dict(diagnostics=bool(o.get("diagnostics")),
                covariatedetail=bool(o.get("covariatedetail")) and not chain,
                varlabels=bool(o.get("varlabels")),
                cformat=o.get("cformat"), pformat=o.get("pformat"), sformat=o.get("sformat"),
                citable=bool(bs.get("citable")))


def _ci_frames(res, ci, cilevel):
    """The intervals shown, as DataFrames named as Stata's r() matrices."""
    import pandas as pd
    out = {}
    bsn = str(res.e.get("bsnames", "")).split()
    for m in CI_KEYS[:-1]:
        if m in ci:
            out[m] = pd.DataFrame(np.asarray(ci[m]), index=["ll", "ul"], columns=bsn)
    if "pct_explained_ci" in ci:
        out["pct_explained_ci"] = pd.DataFrame(np.asarray(ci["pct_explained_ci"]),
                                               index=ci.get("_pct_rows"), columns=["ll", "ul"])
    out["level"] = float(cilevel)
    return out


def _level(x):
    """level(): a confidence level between 10 and 99.99."""
    from .errors import SeldecompError
    try:
        v = float(x)
    except (TypeError, ValueError):
        raise SeldecompError("level() must be a number", 198)
    if not 10 <= v <= 99.99:
        raise SeldecompError("level() must be between 10 and 99.99 inclusive", 198)
    return v


# ---- the interface ----------------------------------------------------------------------

def render(res):
    """The output of the command at estimation, as text: the messages before the tables, the
    notes of bayesboot(), then the display with the display options of the command."""
    lines = []
    for m in res.messages:
        lines += _message_lines(str(m))
    lines += _boot_saved(res)
    lines += _boot_notes(res)
    opts = _estimation_options(res)
    d = _Display(res, opts.pop("varlabels"))
    d.display(_num(res.e.get("level")), replay=False, **opts)
    return "\n".join(lines + d.out)


def _replay(res, level=None, diagnostics=False, covariatedetail=False, citable=False, varlabels=False,
            cformat=None, pformat=None, sformat=None):
    lv = _num(res.e.get("level")) if level in (None, False, "") else _level(level)
    d = _Display(res, varlabels)
    ci, cilevel = d.display(lv, diagnostics=bool(diagnostics), covariatedetail=bool(covariatedetail),
                            cformat=cformat or None, pformat=pformat or None, sformat=sformat or None,
                            citable=bool(citable), replay=True)
    return d, ci, cilevel


def replay_text(res, **options):
    """The output of a replay (the command typed without arguments, with display options:
    level, diagnostics, covariatedetail, citable, varlabels, cformat, pformat, sformat). With
    bayesboot() results, a replay at another level than the estimation's makes that the level
    of estat bootstrap, as the replays after it keep it current."""
    d, ci, cilevel = _replay(res, **_replay_options(options))
    if res.e.get("bsnames"):
        rec = getattr(res, "_replay", None)
        if rec is None:
            rec = res._replay = {}
        lv = float(cilevel)
        if lv != _num(res.e.get("level")) or "level" in rec:
            rec["level"] = lv
    return "\n".join(d.out)


def _replay_options(options):
    """The replay's options under their full names (Stata's abbreviations accepted)."""
    from ._parse import match_name
    from .errors import SeldecompError
    out = {}
    for k, v in options.items():
        try:
            nm = match_name(k, REPLAY_OPTIONS)
        except SeldecompError:
            raise SeldecompError("option %s not allowed at replay; replay takes level(), diagnostics, "
                                 "covariatedetail, citable, varlabels, cformat(), pformat(), sformat(), and "
                                 "makeplot()" % (k if v is True else k + "()"), 198)
        if nm == "makeplot":
            raise SeldecompError("makeplot() is not drawn by res.display(); use res.makeplot()", 198)
        out[nm] = v
    return out


REPLAY_OPTIONS = {"level": 1, "diagnostics": 4, "covariatedetail": 3, "makeplot": 4, "cformat": 7,
                  "pformat": 7, "sformat": 7, "citable": 3, "varlabels": 4}
ESTAT_OPTIONS = {"all": 3, "normal": 3, "percentile": 1, "bc": 2, "bca": 3, "nolegend": 3, "verbose": 1}


def _estat_options(options):
    """estat bootstrap's options under their full names (legend=False is nolegend), refused as
    seldecomp_estat.ado refuses them."""
    from ._parse import match_name
    from .errors import SeldecompError
    out = {}
    for k, v in options.items():
        if k == "legend":
            out["nolegend"] = not v
            continue
        if k in ("l", "le", "lev", "leve", "level"):
            raise SeldecompError("option level() not allowed: estat bootstrap uses the level of the last "
                                 "display; type seldecomp, level(#) first", 198)
        out[match_name(k, ESTAT_OPTIONS, "estat bootstrap: option")] = bool(v)
    return out


def estat_text(res, **options):
    """The output of estat bootstrap (options all, normal, percentile, bc, bca, nolegend,
    verbose): the intervals at the level of the last display."""
    return "\n".join(_estat(res, **options)[0].out)


def _estat(res, **options):
    from .errors import SeldecompError
    if not res.e.get("bsnames"):
        raise SeldecompError("estat bootstrap requires results from bayesboot() or the bootstrap prefix", 321)
    opts = _estat_options(options)
    if opts.get("bca") and _num(res.e.get("bca")) != 1:
        raise SeldecompError("option bca requires results from bayesboot() with bca", 198)
    d = _Display(res, False)
    ci, cilevel = d.estat(**opts)
    return d, ci, cilevel


def display(res, **options):
    """Print the replay's output."""
    print(replay_text(res, **options))


def estat_bootstrap(res, **options):
    """Print estat bootstrap's output; return the intervals shown (Stata's r(): level, ci_normal,
    ci_t, ci_percentile, ci_bc, ci_bca, pct_explained_ci)."""
    d, ci, cilevel = _estat(res, **options)
    print("\n".join(d.out))
    return _ci_frames(res, ci, cilevel)

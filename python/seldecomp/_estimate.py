"""One decomposition on a sample: coding, stages, the point estimates, and the notes."""
import numpy as np

from ._engine import Design, EngineFail, engine
from ._parse import unquote
from ._prepare import (bin_codes, code_term, cut_points, fmt_count, fmt_val,
                       recode_present, stage_text)
from ._stata import fmt_g, plural
from .errors import SeldecompError, warn

NAN = float("nan")
RESERVED_BLOCKS = ("overall", "explained", "marg_expl", "cum_expl", "cum_unexpl", "raked_focal")


class Fit:
    """What one decomposition (the whole run, or one pair of a chain) leaves behind."""


def block_names(blocks, terms, staged, variant, messages, quiet):
    """Default and checked block names; a message for a block shown as block<b>."""
    import re
    for b, blk in enumerate(blocks):
        ok = re.match(r"^[A-Za-z_][A-Za-z_0-9]*$", blk.name) is not None and len(blk.name) <= 32
        if not ok:
            if staged and not quiet:
                messages.append("note: block %d is shown as block%d; name it with (name: ...)" % (b + 1, b + 1))
            blk.name = "block%d" % (b + 1)
    if staged:
        seen = []
        for blk in blocks:
            if blk.name in seen:
                raise SeldecompError("duplicate block name %s" % blk.name, 198)
            seen.append(blk.name)
            if blk.name in RESERVED_BLOCKS or (variant and blk.name == "_remainder"):
                raise SeldecompError("%s cannot be used as a block name" % blk.name, 198)


def check_fixed(blocks, terms, detail):
    """With fixothermargins a variable may be in one block only (with detail: one term only)."""
    bv = [[v for t in blk.terms for v in terms[t].vars] for blk in blocks]
    for b in range(len(blocks)):
        for b2 in range(b + 1, len(blocks)):
            common = [v for v in bv[b] if v in bv[b2]]
            if common:
                raise SeldecompError("with fixothermargins, a variable may be in one block only: "
                                     "%s is in blocks %s and %s" % (common[0], blocks[b].name, blocks[b2].name), 198)
    if detail:
        for blk in blocks:
            for i, t1 in enumerate(blk.terms):
                for t2 in blk.terms[i + 1:]:
                    common = [v for v in terms[t1].vars if v in terms[t2].vars]
                    if common:
                        raise SeldecompError("with detail and fixothermargins, a variable may be in one "
                                             "term only: %s is in %s and %s"
                                             % (common[0], terms[t1].name, terms[t2].name), 198)


def stage_defs(blocks, mode):
    """(stages, names): the terms raked to the reference group at each stage."""
    allterms = [t for blk in blocks for t in blk.terms]
    if mode == "individually":
        stages = [list(blk.terms) for blk in blocks] + [allterms]
        names = [blk.name for blk in blocks] + ["overall"]
    elif mode == "sequentially":
        stages, names, cum = [], [], []
        for blk in blocks:
            cum = cum + list(blk.terms)
            stages.append(list(cum))
            names.append(blk.name)
    else:
        stages, names = [allterms], ["overall"]
    return stages, names


def miss_message(pre, nst, s, nbad, itstop, wtxt, sup=False, trim=False):
    """The error of a raking that did not reach its targets."""
    head = "%sraking did not reach its targets%s: %d %s missed by more than ctrltolerance()" % (
        pre, " at stage %d" % s if nst > 1 else "", nbad, plural(nbad, "margin"))
    if sup:
        return (head + wtxt + "; the lower trimming bounds keep focal group observations with levels "
                "absent from the reference group above 0; specify trimsupportoverride, insupport, or force")
    if trim:
        return head + wtxt + "; the trimming bounds are binding: relax them, or specify force to keep the results"
    return head + (" (raking stopped at iterate())" if itstop else "") + wtxt + "; specify force to keep the results"


def worst_share(rw, codes, level, w, g1, W2):
    """"raked share 32.92%, target 49.71%": the raked focal share of a level and the
    reference group's share."""
    ok = ~np.isnan(rw)
    a = 100 * float(np.sum(rw[ok & (codes == level)])) / float(np.sum(rw[ok]))
    b = 100 * float(np.sum(w[(~g1) & (codes == level)])) / W2
    fa = "%.2f%%" % a
    fb = "%.2f%%" % b
    if fa == fb:
        fa, fb = "%.6f%%" % a, "%.6f%%" % b
    return "raked share %s, target %s" % (fa, fb)


def estimate(ctx, touse, focal, ref, pair=0, quiet=False):
    """The decomposition of the difference between the focal and the reference value of by()
    on the observations touse (a boolean mask over the rows of the data). Returns a Fit."""
    o = ctx.opt
    df = ctx.df
    chp = "pair %d: " % pair if pair else ""
    fit = Fit()
    notes = []
    terms, blocks = ctx.terms, [type(b)(b.name, list(b.terms)) for b in ctx.blocks]
    nterms, nblocks = len(terms), len(blocks)
    mode = "sequentially" if o["sequentially"] else ("individually" if o["individually"] else "")
    staged = mode != ""
    fit.messages = []
    # (a chain says it once, with its first pair)
    block_names(blocks, terms, staged, mode == "individually", fit.messages, quiet or pair > 1)
    if staged and o["fixothermargins"]:
        check_fixed(blocks, terms, o["detail"])

    byv = ctx.byvals
    g1_full = np.zeros(len(df), dtype=bool)
    g1_full[touse] = np.array([v == focal for v in byv[touse]])
    w_full = ctx.w

    # cut points of the binned continuous variables
    cutinfo = {}
    cutbins = []
    for v, q in ctx.contlist:
        if ctx.fixed_cuts is not None and (pair, v) in ctx.fixed_cuts:
            cuts = ctx.fixed_cuts[(pair, v)]
        else:
            if ctx.contsample == "pooled":
                cond = touse
            elif ctx.contsample == "focalgroup":
                cond = touse & g1_full
            else:
                cond = touse & ~g1_full
            x = df[v].to_numpy(dtype=float)
            cuts = cut_points(x[cond], w_full[cond], q, v)
        cutinfo[v] = (cuts, [fmt_g(c, 9) for c in cuts])
        x = df[v].to_numpy(dtype=float)
        cutbins.append(len(np.unique(bin_codes(x[touse], cuts))))
    fit.cuts = {v: cutinfo[v][0] for v, q in ctx.contlist}
    fit.cutbins = cutbins

    coded = [code_term(t, df, touse, cutinfo, ctx.value_labels) for t in terms]
    for t, c in enumerate(coded):
        if c.nlev == 1 and not quiet:
            fit.messages.append("note: %s%s is constant in the estimation sample" % (chp, terms[t].name))

    idx = np.flatnonzero(touse)
    g1 = g1_full[idx]
    w = w_full[idx]
    ys = [df[y].to_numpy(dtype=float)[idx] for y in ctx.ylist]
    den = df[o["ratio"]].to_numpy(dtype=float)[idx] if o["ratio"] else None

    # insupport: focal observations with a level that the reference group lacks leave the sample
    fit.ins = None
    insnotes = []
    if o["insupport"]:
        N1all = int(g1.sum())
        W1all = float(np.sum(w[g1]))
        m1all = []
        d1all = float(np.sum(den[g1] * w[g1]) / W1all) if den is not None else None
        for y in ys:
            m = float(np.sum(y[g1] * w[g1]) / W1all)
            if den is not None:
                m = m / d1all if d1all != 0 else NAN
            m1all.append(m)
        out = np.zeros(len(idx), dtype=bool)
        inslist = []
        in2 = []
        for t, c in enumerate(coded):
            has2 = np.bincount(c.codes[~g1], minlength=c.nlev + 1)[1:] > 0
            cnt1 = np.bincount(c.codes[g1], minlength=c.nlev + 1)[1:]
            in2.append(has2)
            for l in np.flatnonzero(~has2 & (cnt1 > 0)):
                out |= g1 & (c.codes == l + 1)
                inslist.append((t, int(l) + 1, int(cnt1[l])))
        insfail = []
        for t, c in enumerate(coded):
            for l in np.flatnonzero(in2[t]):
                m = g1 & (c.codes == l + 1)
                if m.any() and not (m & ~out).any():
                    insfail.append((t, int(l) + 1))
        nexcl = int(out.sum())
        if nexcl == N1all:
            tot = {}
            for t, l, cnt in inslist:
                tot[t] = tot.get(t, 0) + cnt
            sole = [t for t in tot if tot[t] == N1all]
            if sole:
                raise SeldecompError("insupport excludes every focal group observation (levels of %s "
                                     "absent from the reference group)" % terms[sole[0]].name, 2000)
            txt = "; ".join("%s: %s obs" % (terms[t].name, fmt_count(n)) for t, n in tot.items())
            raise SeldecompError("insupport excludes every focal group observation (levels absent "
                                 "from the reference group: %s)" % txt, 2000)
        if insfail:
            lines = ["insupport: after the exclusion, these levels of the reference group (%s = %s) "
                     "have no focal group observations:" % (ctx.by, ctx.group_label(ref))]
            for k, (t, l) in enumerate(insfail):
                if k >= 10:
                    lines.append("    ... and %d more level(s)" % (len(insfail) - 10))
                    break
                lines.append("    %s = %s" % (terms[t].name, coded[t].labels[l - 1]))
            raise SeldecompError("\n".join(lines), 459)
        if nexcl > 0:
            if len(inslist) == 1:
                t, l, cnt = inslist[0]
                insnotes.append("insupport excluded %s = %s (%s obs)"
                                % (terms[t].name, coded[t].labels[l - 1], fmt_count(cnt)))
            else:
                insnotes.append("insupport excluded:")
                for k, (t, l, cnt) in enumerate(inslist):
                    if k >= 10:
                        insnotes.append("    ... (%d more)" % (len(inslist) - 10))
                        break
                    insnotes.append("    %s = %s (%s obs)" % (terms[t].name, coded[t].labels[l - 1], fmt_count(cnt)))
            keep = ~out
            before = [c.nlev for c in coded]
            coded = [recode_present(c, keep) for c in coded]
            for t, c in enumerate(coded):
                if c.nlev == 1 and before[t] > 1 and not quiet:
                    fit.messages.append("note: %s%s is constant in the estimation sample" % (chp, terms[t].name))
            idx = idx[keep]
            g1, w = g1[keep], w[keep]
            ys = [y[keep] for y in ys]
            den = den[keep] if den is not None else None
        fit.ins = {"N_focal_all": N1all, "W_focal_all": W1all, "mean_focal_all": m1all, "nexcl": nexcl}

    fit.idx = idx
    fit.g1 = g1
    fit.w = w
    fit.coded = coded
    fit.terms = terms
    fit.blocks = blocks
    fit.focal, fit.ref = focal, ref
    D = Design(g1, ys, den, [c.codes for c in coded], [c.nlev for c in coded], [t.name for t in terms])
    fit.design = D
    stages, snames = stage_defs(blocks, mode)
    nstages = len(stages)
    fit.stages, fit.stagenames, fit.mode = stages, snames, mode
    hold = bool(o["fixothermargins"])
    contrib = mode == "individually"
    fit.hold, fit.contrib = hold, contrib
    trace = ctx.trace if o["noisily"] else None

    def lev(t, l):
        return coded[t].labels[l - 1]

    def run_engine(design, stg, prefix, hold_=hold, contrib_=contrib, stage_label=None):
        try:
            return engine(design, w, stg, ctx.ropt, holdothers=hold_, contrib=contrib_,
                          want_weights=True, want_misses=True, trace=trace)
        except EngineFail as f:
            raise support_error(f, prefix, len(stg), stage_label)

    def support_error(f, prefix, nst, stage_label=None):
        info = f.info
        if f.rc == 459 and "fails" in info:
            fails = info["fails"]
            lines = ["%sraking support failure%s: these levels of the reference group (%s = %s) have no "
                     "focal group observations (%s = %s):"
                     % (prefix, " at stage %d" % (info["stage"] + 1) if nst > 1 else "",
                        ctx.by, ctx.group_label(ref), ctx.by, ctx.group_label(focal))]
            for k, (t, l, z) in enumerate(fails):
                if k >= 10:
                    lines.append("    ... and %d more %s" % (len(fails) - 10, plural(len(fails) - 10, "level")))
                    break
                lines.append("    %s = %s%s" % (terms[t].name, lev(t, l),
                                                " (only in set-aside focal group observations)" if z else ""))
            if info.get("nzero", 0) > 0:
                nz = info["nzero"]
                lines.append("%d focal group %s set aside (levels absent from the reference group):"
                             % (nz, plural(nz, "observation")))
                for k, (t, l, cnt) in enumerate(info["zerolist"]):
                    if k >= 10:
                        lines.append("      ...")
                        break
                    lines.append("      %s = %s (%d obs)" % (terms[t].name, lev(t, l), cnt))
            rem = "collapse or drop these levels"
            if not (o["insupport"] and fit.ins and fit.ins["nexcl"] > 0) and pair == 0 and other_ok():
                oin = " in place of focalgroup(%s)" % unquote(str(o["focalgroup"]).strip()) \
                    if o["focalgroup"] is not None else ""
                rem += ", or rake the other group with focalgroup(%s)%s" % (fmt_val(ref), oin)
            lines.append(rem)
            return SeldecompError("\n".join(lines), 459)
        if f.rc == 2000 and "zerolist" in info:
            s = info["stage"] + 1
            if info.get("septerm") is not None:
                msg = "stage %d: every focal group observation has a level of %s absent from the reference group" % (
                    s, terms[info["septerm"]].name)
            else:
                tot = {}
                for t, l, cnt in info["zerolist"]:
                    tot[t] = tot.get(t, 0) + cnt
                txt = "; ".join("%s: %s obs" % (terms[t].name, fmt_count(n)) for t, n in tot.items())
                msg = "stage %d: every focal group observation has a level absent from the reference group (%s)" % (s, txt)
            return SeldecompError(chp + prefix + msg, 2000)
        msg = f.msg.replace("{den}", o["ratio"] or "")
        return SeldecompError(chp + prefix + msg, f.rc)

    def other_ok():
        # the support check with the roles of the groups swapped
        cod = coded
        if ctx.contsample != "pooled" and ctx.contlist:
            # cut points from one group's sample: raking the other group cuts at the other
            # group's quantiles, so the binned terms are coded as that run would code them
            if ctx.fixed_cuts is not None:
                return False
            other = touse & (g1_full if ctx.contsample == "refgroup" else ~g1_full)
            cutinfo2 = {}
            for v, q in ctx.contlist:
                x = df[v].to_numpy(dtype=float)
                try:
                    cuts = cut_points(x[other], w_full[other], q, v)
                except SeldecompError:
                    return False
                cutinfo2[v] = (cuts, [fmt_g(c, 9) for c in cuts])
            cod = [code_term(t, df, touse, cutinfo2, ctx.value_labels) for t in terms]
        keep = ~g1
        for c in cod:
            in1 = np.bincount(c.codes[g1], minlength=c.nlev + 1) > 0
            keep = keep & in1[c.codes]
        if not keep.any():
            return False
        for c in cod:
            in1 = np.bincount(c.codes[g1], minlength=c.nlev + 1) > 0
            in2 = np.bincount(c.codes[keep], minlength=c.nlev + 1) > 0
            if np.any(in1 & ~in2):
                return False
        return True

    check = getattr(ctx, "before_raking", None)
    if check is not None:
        ctx.before_raking = None
        check()
    E = run_engine(D, stages, "")
    fit.E = E
    RES = E.res
    W2 = E.W2

    # focal observations held above zero weight by a lower bound: a warning
    if E.held and not o["force"] and not quiet:
        hn, hl = 0, []
        for s, n, zl in E.held:
            if n > hn:
                hn, hl = n, zl
        htx = "; ".join("%s = %s: %s obs" % (terms[t].name, lev(t, l), fmt_count(cnt)) for t, l, cnt in hl[:5])
        if len(hl) > 5:
            htx += "; and %d more" % (len(hl) - 5)
        lb = [nm for nm, val in (("trimloabs()", ctx.ropt.trimloabs), ("trimlorel()", ctx.ropt.trimlorel)) if val is not None]
        lbtxt = " and ".join(lb)
        text = ("%s%s focal group %s %s a level absent from the reference group (%s). %s %s %s raked %s "
                "above 0, so raking will likely fail. trimsupportoverride sets %s to 0."
                % (chp, fmt_count(hn), plural(hn, "observation"), "has" if hn == 1 else "have", htx, lbtxt,
                   "keep" if len(lb) > 1 else "keeps", "its" if hn == 1 else "their", plural(hn, "weight"),
                   "it" if hn == 1 else "them"))
        warn(text)
        fit.messages.append("warning: " + text)

    def check_misses(Eo, RESo, nst, prefix, term_of=lambda t: t, stage_filter=None, rwcodes=None, raked_terms=None):
        """A stage with a missed margin: error 430 unless force. Returns the last such stage."""
        bad = None
        for s in range(nst):
            if stage_filter is not None and not stage_filter(s):
                continue
            if RESo[0][s, 6] > 0:
                bad = s
        if bad is None:
            return None
        if not o["force"]:
            nbad = int(RESo[0][bad, 6])
            wt, wl = Eo.worst[bad]
            wtxt = ""
            if wt is not None and wl is not None:
                to = term_of(wt)
                wtxt = "; worst: %s = %s" % (terms[to].name, lev(to, wl))
                if raked_terms is None or wt in raked_terms(bad):
                    wtxt += " (%s)" % worst_share(Eo.rw[bad], coded[to].codes, wl, w, g1, W2)
            raise SeldecompError(miss_message(chp + prefix, nst, bad + 1, nbad, RESo[0][bad, 5] == 0, wtxt,
                                              sup=Eo.supmiss[bad] > 0, trim=bool(Eo.trimchg[bad])), 430)
        return bad

    badstage = check_misses(E, RES, nstages, "", raked_terms=lambda s: stages[s])
    rconv = [RES[0][s, 5] for s in range(nstages)]
    tmet = 1
    for s in range(nstages):
        if RES[0][s, 6] > 0:
            tmet = 0
            for R in RES:
                R[s, 5] = 0
    fit.targets_met = tmet

    # the notes: exclusions, levels absent from the reference group, thin levels, stage notes
    for nt in ctx.xnotes if pair == 0 else []:
        notes.append(nt)
    notes.extend(insnotes)
    for kind, data in (("zero", E.zero), ("thin", E.thin)):
        keys, st_of = [], {}
        for s, n, cells in data:
            for cell in cells:
                if cell not in keys:
                    keys.append(cell)
                st_of.setdefault(cell, []).append(s + 1)
        groups = []
        for cell in keys:
            txt = stage_text(nstages, st_of[cell])
            for g in groups:
                if g[0] == txt:
                    g[1].append(cell)
                    break
            else:
                groups.append((txt, [cell]))
        for txt, cells in groups:
            nc = len(cells)
            stx = ", " + txt if txt else ""
            nct = "%d " % nc if nc > 1 else ""
            if kind == "zero":
                notes.append("%s%s absent from the reference group (zero raked weight)%s:" % (nct, plural(nc, "level"), stx))
            else:
                notes.append("%sthin %s%s:" % (nct, plural(nc, "target level"), stx))
            for k, (t, l, cnt) in enumerate(cells):
                if kind == "thin" and k >= 8:
                    notes.append("    ... (%d more)" % (nc - 8))
                    break
                notes.append("    %s = %s (%s obs)" % (terms[t].name, lev(t, l), fmt_count(cnt)))
    stoptxt = "at the iterate() limit"
    mt, mta, mtb = {}, {}, {}
    itstages = []
    worstvar, worstcat = {}, {}
    for s in range(nstages):
        wt, wl = E.worst[s]
        worstvar[s] = terms[wt].name if wt is not None else ""
        worstcat[s] = lev(wt, wl) if wt is not None and wl is not None else ""
        nbad = int(RES[0][s, 6])
        plm = plural(nbad, "margin")
        if nbad > 0 and E.trimchg[s]:
            ntr = int(RES[0][s, 21])
            ntx = "; %s %s at a bound" % (fmt_count(ntr), plural(ntr, "weight")) if ntr > 0 else ""
            mta[s] = "trimming left %d raking %s unmet (worst: %s = %s%s" % (nbad, plm, worstvar[s], worstcat[s], ntx)
            mtb[s] = "); results kept because force was specified"
            mt[s] = mta[s] + mtb[s]
        elif nbad > 0:
            stp = "raking stopped %s with " % stoptxt if rconv[s] == 0 else ""
            mta[s] = "%s%d raking %s not met (worst: %s = %s" % (stp, nbad, plm, worstvar[s], worstcat[s])
            mtb[s] = "), kept by force" if rconv[s] == 0 else "); results kept because force was specified"
            mt[s] = mta[s] + mtb[s]
        elif rconv[s] == 0:
            mt[s] = "raking stopped %s before the weights converged; raise iterate()" % stoptxt
            itstages.append(s)
    fnote, fnotea, fnoteb = None, "", ""
    fit.inote = None
    done = set()
    for s in range(nstages):
        if s in mt and s not in done:
            sl = [s]
            for s2 in range(s + 1, nstages):
                if mt.get(s2) == mt[s]:
                    sl.append(s2)
                    done.add(s2)
            sttxt = stage_text(nstages, [x + 1 for x in sl])
            stx = sttxt + ": " if sttxt else ""
            notes.append(stx + mt[s])
            if s in itstages:
                fit.inote = len(notes) - 1
                fit.istx = stx
            if s in mta:
                fnote = len(notes) - 1
                fnotea = stx + mta[s]
                fnoteb = mtb[s]
    fit.itstages = itstages
    fit.fnotes = [fnote] if fnote is not None else []

    # reverse order: the blocks added last to first, the terms renumbered block-wise
    fit.ER = None
    if o["reverse"]:
        rorder = [t for blk in reversed(blocks) for t in blk.terms]
        Dr = Design(g1, ys, den, [coded[t].codes for t in rorder], [coded[t].nlev for t in rorder],
                    [terms[t].name for t in rorder])
        rstages, rnames, cum = [], [], 0
        for blk in reversed(blocks):
            cum += len(blk.terms)
            rstages.append(list(range(cum)))
            rnames.append(blk.name)
        fit.rorder, fit.rdesign, fit.rstages, fit.rnames = rorder, Dr, rstages, rnames
        ER = run_engine(Dr, rstages, "reverse order: ")
        fit.ER = ER
        rbad = check_misses(ER, ER.res, nstages, "reverse order: ", term_of=lambda t: rorder[t],
                            raked_terms=lambda s: rstages[s])
        if rbad is not None:
            if fnote is not None:
                notes[fnote] = fnotea + "; also in reverse order" + fnoteb
                fnotea = fnotea + "; also in reverse order"
            else:
                nrb = int(ER.res[0][rbad, 6])
                wt, wl = ER.worst[rbad]
                wtxt = ""
                if wt is not None and wl is not None:
                    wo = rorder[wt]
                    wtxt = "worst: %s = %s" % (terms[wo].name, lev(wo, wl))
                rstem = "%d raking %s not met" % (nrb, plural(nrb, "margin"))
                if ER.trimchg[rbad]:
                    rstem = "trimming left %d raking %s unmet" % (nrb, plural(nrb, "margin"))
                    ntr = int(ER.res[0][rbad, 21])
                    if ntr > 0:
                        wtxt = (wtxt + "; " if wtxt else "") + "%s %s at a bound" % (fmt_count(ntr), plural(ntr, "weight"))
                if wtxt:
                    wtxt = " (%s)" % wtxt
                notes.append("reverse order: stage %d: %s%s; results kept because force was specified" % (rbad + 1, rstem, wtxt))
                fit.fnotes.append(len(notes) - 1)
        for s in range(nstages):
            if ER.res[0][s, 6] > 0:
                for R in ER.res:
                    R[s, 5] = 0
        if not ctx.ropt.hastrim and badstage is None and rbad is None:
            fd = 0.0
            for j in range(len(ys)):
                a, b = RES[j][nstages - 1, 0], ER.res[j][nstages - 1, 0]
                fd = max(fd, abs(a - b) / (abs(b) + 1))
            if fd > 1000 * ctx.ropt.tolerance:
                notes.append("the final stages of the two orders differ by %s (relative); tighten tolerance()"
                             % ("%.2g" % fd))

    # detail: the terms added one at a time
    fit.ED = None
    if o["detail"]:
        dends = [blk.terms[-1] for blk in blocks]
        dblk = {}
        for b, blk in enumerate(blocks):
            for t in blk.terms:
                dblk[t] = b
        import re
        dnames = []
        for t, term in enumerate(terms):
            dn = term.name.replace("#", "_")
            if re.match(r"^[A-Za-z_][A-Za-z_0-9]*$", dn) is None or len(dn) > 32 or dn in dnames:
                if not quiet:
                    fit.messages.append("note: term %d (%s) is shown as term%d in detail" % (t + 1, term.name, t + 1))
                dn = "term%d" % (t + 1)
            dnames.append(dn)
        fit.dends, fit.dnames = dends, dnames
        dstages = [list(range(t + 1)) for t in range(nterms)]
        fit.dstages = dstages
        if nterms == nblocks:
            ED = E
            fit.detail_same = True
        else:
            fit.detail_same = False
            try:
                ED = engine(D, w, dstages, ctx.ropt, holdothers=hold, want_weights=True, want_misses=True, trace=trace)
            except EngineFail as f:
                info = f.info
                if f.rc == 459 and "fails" in info:
                    raise support_error(f, "detail: ", nterms)
                err = support_error(f, "", nterms)
                msg = err.message[len(chp):]
                m = re.search(r"stage ([0-9]+)", msg)
                if m:
                    msg = msg.replace("stage %s" % m.group(1), "adding %s" % terms[int(m.group(1)) - 1].name, 1)
                raise SeldecompError("detail: " + msg, err.rc)
            dbad = None
            for s in range(nterms):
                if s not in dends and ED.res[0][s, 6] > 0:
                    dbad = s
            if dbad is not None:
                wt, wl = ED.worst[dbad]
                nbd = int(ED.res[0][dbad, 6])
                wtxt = ""
                if wt is not None and wl is not None:
                    wtxt = "worst: %s = %s" % (terms[wt].name, lev(wt, wl))
                if not o["force"]:
                    stx = " (raking stopped at iterate())" if ED.res[0][dbad, 5] == 0 else ""
                    dfix = "specify force to keep the results"
                    if ED.trimchg[dbad]:
                        stx = ""
                        dfix = "the trimming bounds are binding: relax them, or " + dfix
                    raise SeldecompError("%sdetail: adding %s (block %s) misses %d raking %s%s%s; %s"
                                         % (chp, terms[dbad].name, blocks[dblk[dbad]].name, nbd, plural(nbd, "margin"),
                                            " (%s)" % wtxt if wtxt else "", stx, dfix), 430)
                if fnote is not None:
                    notes[fnote] = fnotea + "; also in the detail, adding %s" % terms[dbad].name + fnoteb
                else:
                    dstem = "%d raking %s not met" % (nbd, plural(nbd, "margin"))
                    if ED.trimchg[dbad]:
                        dstem = "trimming left %d raking %s unmet" % (nbd, plural(nbd, "margin"))
                        ntr = int(ED.res[0][dbad, 21])
                        if ntr > 0:
                            wtxt = (wtxt + "; " if wtxt else "") + "%s %s at a bound" % (fmt_count(ntr), plural(ntr, "weight"))
                    notes.append("detail, adding %s (block %s): %s%s; results kept because force was specified"
                                 % (terms[dbad].name, blocks[dblk[dbad]].name, dstem, " (%s)" % wtxt if wtxt else ""))
                    fit.fnotes.append(len(notes) - 1)
            for s in range(nterms):
                if ED.res[0][s, 6] > 0:
                    for R in ED.res:
                        R[s, 5] = 0
        fit.ED = ED

    fit.notes = notes
    fit.worstvar, fit.worstcat = worstvar, worstcat
    fit.b_names, fit.b_values = b_layout(fit, ctx.ylist)
    return fit


def b_layout(fit, ylist, prefix=""):
    """The names and values of the coefficient vector of one decomposition."""
    E = fit.E
    nout = len(ylist)
    nst = len(fit.stages)
    names, vals = [], []
    staged = fit.mode != ""
    variant = fit.mode == "individually"
    for j in range(nout):
        eq0 = "overall"
        pfx = ""
        if nout > 1:
            eq0 = ylist[j]
            pfx = eq0 + "_"
        if prefix:
            eq0 = prefix
            pfx = prefix + "_"
        R = E.res[j]
        last = R[nst - 1]
        vals += [E.m1[j], last[0], E.m2[j], E.m1[j] - E.m2[j], last[1], last[2]]
        names += ["%s:%s" % (eq0, n) for n in ("focal_group", "raked_focal", "ref_group", "difference",
                                               "explained", "unexplained")]
        if staged:
            orders = [(R, "", fit.stagenames)]
            if fit.ER is not None:
                orders.append((fit.ER.res[j], "rev_", fit.rnames))
            for k, (M, rp, sn) in enumerate(orders):
                cols = [(4, "marg_expl"), (0, "raked_focal")] if variant else \
                    [(4, "marg_expl"), (1, "cum_expl"), (2, "cum_unexpl"), (0, "raked_focal")]
                for c, eqn in cols:
                    for s in range(nst):
                        if variant and c == 0 and s == nst - 1:
                            continue
                        sname = "_remainder" if (variant and s == nst - 1) else sn[s]
                        vals.append(M[s, c])
                        names.append("%s%s%s:%s" % (pfx, rp, eqn, sname))
                if k == 0 and fit.ED is not None:
                    Dm = fit.ED.res[j]
                    for t in range(len(fit.terms)):
                        vals.append(Dm[t, 4])
                        names.append("%sdetail:%s" % (pfx, fit.dnames[t]))
    return names, vals


class RepFail(Exception):
    """A replicate (or a jackknife deletion) that failed: why is "it" (stopped at iterate()
    short of the targets), "tg" (targets missed), "tb" (the trimming bounds contradict), or
    "rk" (the raking failed)."""

    def __init__(self, why):
        Exception.__init__(self, why)
        self.why = why


def rep_values(fit, ctx, w, force, state):
    """The coefficient vector recomputed on the weights w (those of fit's observations), and
    the explained share of each outcome. state collects itl (a raking stopped at iterate()
    with its margins met), frc (a missed margin kept by force), itst (the point estimate's
    stages at which a raking stopped at iterate()) and itrev (one did in the reverse order or
    at a term stage of detail)."""
    D = fit.design
    codes_override = state.get("codes")
    if codes_override is not None:
        D = Design(D.g1, D.ys, D.den, codes_override[0], codes_override[1], D.tnames)
    nst = len(fit.stages)
    det = fit.ED is not None
    stg = fit.dstages if det else fit.stages
    hold, contrib = fit.hold, fit.contrib

    def run(design, stages, contrib_, reverse=False):
        try:
            Eo = engine(design, w, stages, ctx.ropt, holdothers=hold, contrib=contrib_, skipdiag=True)
        except EngineFail as f:
            raise RepFail("tb" if f.info.get("trimclash") else "rk")
        fail = None
        R = Eo.res[0]
        for s in range(len(stages)):
            # st: the point estimate's stage that this raking belongs to (with detail, a term
            # stage that is not the last of its block belongs to none)
            if reverse:
                st = None
            elif det:
                st = fit.dends.index(s) if s in fit.dends else None
                if st is None and R[s, 5] == 0:
                    state["itrev"] = 1
            else:
                st = nst - 1 if len(stages) == 1 else s
            if R[s, 6] > 0 and force:
                state["frc"] = 1
            elif R[s, 6] > 0:
                itb = (R[s, 5] == 0) and not Eo.trimchg[s]
                if fail is None:
                    fail = "it" if itb else "tg"
                if itb and reverse:
                    state["itrev"] = 1
                elif itb and st is not None and st not in state["itst"]:
                    state["itst"].append(st)
            elif R[s, 5] == 0 and R[s, 6] == 0:
                state["itl"] = 1
                if reverse:
                    state["itrev"] = 1
                elif st is not None and st not in state["itst"]:
                    state["itst"].append(st)
        if fail:
            raise RepFail(fail)
        return Eo

    state.setdefault("itst", [])
    state.setdefault("itrev", 0)
    Eo = run(D, stg, contrib and not det)
    ER = None
    if fit.ER is not None:
        Dr = fit.rdesign
        if codes_override is not None:
            Dr = Design(D.g1, D.ys, D.den, [codes_override[0][t] for t in fit.rorder],
                        [codes_override[1][t] for t in fit.rorder], Dr.tnames)
        ER = run(Dr, fit.rstages, False, reverse=True)
    vals = []
    nsb = len(stg)
    variant = fit.mode == "individually"
    staged = fit.mode != ""
    for j in range(len(D.ys)):
        R = Eo.res[j]
        vals += [Eo.m1[j], R[nsb - 1, 0], Eo.m2[j], Eo.m1[j] - Eo.m2[j], R[nsb - 1, 1], R[nsb - 1, 2]]
        if staged:
            cols = [4, 0] if variant else [4, 1, 2, 0]
            if det:
                # a block's values are those of the term stage at its last term; its increment
                # the difference of the explained parts at the last terms of it and of the block before
                for c in cols:
                    e0 = None
                    for e in fit.dends:
                        if c != 4:
                            vals.append(R[e, c])
                        elif e0 is None:
                            vals.append(R[e, 1])
                        else:
                            vals.append(R[e, 1] - R[e0, 1])
                        e0 = e
                for s in range(nsb):
                    vals.append(R[s, 4])
            else:
                for c in cols:
                    for s in range(nst):
                        if variant and c == 0 and s == nst - 1:
                            continue
                        vals.append(R[s, c])
            if ER is not None:
                M = ER.res[j]
                for c in cols:
                    for s in range(nst):
                        vals.append(M[s, c])
    pct = [Eo.res[j][nsb - 1, 3] for j in range(len(D.ys))]
    return vals, pct

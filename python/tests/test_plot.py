"""makeplot(): the graphs, read back through the Figure's artists. The bars' starts and lengths,
the intervals' ends and the row names are compared with the stored results (res.b and the
ci_<type> matrices); the labels of the bars with the texts that Stata draws for the same
commands; each refusal raises SeldecompError with rc 198."""
import contextlib
import io
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.text import Text  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (HERE, os.path.join(HERE, "..")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import seldecomp as sd  # noqa: E402
from common import LABELS, VALUE_LABELS, load_data, run_case  # noqa: E402

TOL = 1e-9
SPECS = {
    "seq": "y (a: x1) (b: x2 x4) c.x3 [pw=wt], by(group) contvars(x3 5) sequentially",
    "fixrev": "y (a: x1) (b: x2 x4) c.x3 [pw=wt], by(group) contvars(x3 5) sequentially "
              "fixothermargins reverse",
    "ind": "y (a: x1) (b: x2 x4) c.x3 [pw=wt], by(group) contvars(x3 5) individually",
    "outcomes": "y (a: x1) (b: x2) [pw=wt], by(group) outcomes(ybin earn) sequentially",
    "chain": "y (a: x1) (b: x2) [pw=wt], by(g3) chain sequentially",
    "terms": "y x1 x2 x4 [pw=wt], by(group) sequentially",
    "plain": "y x1 x2 [pw=wt], by(group)",
}
_DATA = []
_FITS = {}


def data():
    if not _DATA:
        _DATA.append(load_data())
    return _DATA[0].copy()


def fit(name):
    """The results of a case, estimated once: a spec above or a bootstrap fixture case."""
    if name not in _FITS:
        if name in SPECS:
            _FITS[name] = sd.seldecomp(data(), "seldecomp " + SPECS[name], labels=LABELS,
                                       value_labels=VALUE_LABELS)
        else:
            _FITS[name] = run_case(name)[2]
    return _FITS[name]


def read(fig, horiz):
    """What a graph shows, by the gid of its artists: bars {b name: (row position, start,
    length)} along the value axis, intervals {b name: (row position, low, high)}, labels {b
    name: text}, the row names, the value-axis limits, the lines and every text of the
    figure. The figure is closed."""
    ax = fig.axes[0]
    out = {"bars": {}, "ci": {}, "labels": {}, "lines": {}}
    for p in ax.patches:
        g = p.get_gid()
        if g is None:
            continue
        if horiz:
            out["bars"][g] = (p.get_y() + p.get_height() / 2, p.get_x(), p.get_width())
        else:
            out["bars"][g] = (p.get_x() + p.get_width() / 2, p.get_y(), p.get_height())
    for c in ax.collections:
        g = c.get_gid() or ""
        if g.startswith("ci:") and not g.endswith(":caps"):
            (x0, y0), (x1, y1) = c.get_segments()[0]
            out["ci"][g[3:]] = (y0, min(x0, x1), max(x0, x1)) if horiz else (x0, min(y0, y1), max(y0, y1))
    for t in ax.texts:
        g = t.get_gid() or ""
        if g.startswith("label:"):
            out["labels"][g[6:]] = t.get_text()
    for ln in ax.lines:
        if ln.get_gid():
            out["lines"][ln.get_gid()] = (list(ln.get_xdata()), list(ln.get_ydata()))
    ticks = ax.get_yticklabels() if horiz else ax.get_xticklabels()
    out["ticks"] = [t.get_text() for t in ticks]
    out["vlim"] = ax.get_xlim() if horiz else ax.get_ylim()
    out["texts"] = [t.get_text() for t in fig.findobj(Text) if t.get_text()]
    out["inside"] = labels_inside(fig, ax, horiz)
    plt.close(fig)
    return out


def labels_inside(fig, ax, horiz):
    """Whether every label of a bar lies within the plot along the value axis."""
    fig.canvas.draw()
    box = ax.get_window_extent()
    for t in ax.texts:
        if not (t.get_gid() or "").startswith("label:"):
            continue
        bb = t.get_window_extent()
        a, b, lo, hi = (bb.x0, bb.x1, box.x0, box.x1) if horiz else (bb.y0, bb.y1, box.y0, box.y1)
        if a < lo - 1 or b > hi + 1:
            return False
    return True


def near(a, b, tol=TOL):
    return abs(a - b) <= tol * max(1.0, abs(b))


def check_bars(shown, rows):
    """rows: (b name, row position, start, length); start None for a mean's bar (from the base)."""
    for name, pos, start, length in rows:
        assert name in shown["bars"], name
        p, s, n = shown["bars"][name]
        assert near(p, pos), (name, p, pos)
        if start is None:
            assert near(s + n, length), (name, s + n, length)
        else:
            assert near(s, start) and near(n, length), (name, s, n, start, length)


def ci_of(res, kind, name):
    C = res.matrices["ci_" + kind]
    return float(C.loc["ll", name]), float(C.loc["ul", name])


def waterfall_rows(res, blocks, pfx="", eq0="overall"):
    b = res.b
    rows, run = [], 0.0
    for k, bn in enumerate(blocks):
        name = pfx + "marg_expl:" + bn
        rows.append((name, k + 1.0, run, b[name]))
        run += b[name]
    E, U, D = b[eq0 + ":explained"], b[eq0 + ":unexplained"], b[eq0 + ":difference"]
    assert near(D - E, U, 1e-9)
    nb = len(blocks)
    rows += [(eq0 + ":explained", nb + 1.0, 0.0, E), (eq0 + ":unexplained", nb + 2.0, E, D - E),
             (eq0 + ":difference", nb + 3.0, 0.0, D)]
    return rows


def levels_rows(res, blocks, pfx="", eq0="overall"):
    b = res.b
    m1, rk, m2 = b[eq0 + ":focal_group"], b[eq0 + ":raked_focal"], b[eq0 + ":ref_group"]
    rows, level = [(eq0 + ":focal_group", 1.0, None, m1)], m1
    for k, bn in enumerate(blocks):
        name = pfx + "marg_expl:" + bn
        rows.append((name, k + 2.0, level, -b[name]))
        level -= b[name]
    assert near(level, rk, 1e-9)
    nb = len(blocks)
    rows += [(eq0 + ":raked_focal", nb + 2.0, None, rk), (eq0 + ":unexplained", nb + 3.0, rk, m2 - rk),
             (eq0 + ":ref_group", nb + 4.0, None, m2)]
    return rows


def marginal_rows(res, blocks, pfx="", eq0="overall"):
    b = res.b
    nb = len(blocks)
    rows = [(eq0 + ":difference", 1.0, 0.0, b[eq0 + ":difference"]),
            (eq0 + ":explained", 2.2, 0.0, b[eq0 + ":explained"])]
    rows += [(pfx + "marg_expl:" + bn, k + 3.2, 0.0, b[pfx + "marg_expl:" + bn])
             for k, bn in enumerate(blocks)]
    rows.append((eq0 + ":unexplained", nb + 3.4, 0.0, b[eq0 + ":unexplained"]))
    return rows


def stages_rows(res, blocks, pfx=""):
    b = res.b
    rows = []
    for s, bn in enumerate(blocks):
        rows.append((pfx + "cum_expl:" + bn, s + 0.8, 0.0, b[pfx + "cum_expl:" + bn]))
        rows.append((pfx + "cum_unexpl:" + bn, s + 1.2, 0.0, b[pfx + "cum_unexpl:" + bn]))
    return rows


def refused(fn, text=None, rc=198):
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            fn()
    except sd.SeldecompError as exc:
        assert exc.rc == rc, (exc.rc, exc.message)
        if text is not None:
            assert text in exc.message, exc.message
        return exc.message
    finally:
        plt.close("all")
    raise AssertionError("not refused: %s" % text)


# ---- the four types after sequentially ---------------------------------------------------------

def test_waterfall_sequentially():
    res = fit("seq")
    shown = read(res.makeplot(), horiz=False)
    check_bars(shown, waterfall_rows(res, ["a", "b", "x3"]))
    assert shown["ticks"] == ["a", "b", "x3", "Expl.", "Unexpl.", "Total diff"]
    # connectors from each bar's end to the next bar
    b = res.b
    xs, ys = shown["lines"]["connector:marg_expl:a"]
    assert near(xs[0], 1.3) and near(xs[1], 1.7) and near(ys[0], b["marg_expl:a"])
    assert "connector:overall:difference" not in shown["lines"]
    assert "Decomposition of the difference by raking stage" in shown["texts"]
    assert "Group zero - Group one" in shown["texts"]
    assert "Outcome: difference in means" in shown["texts"]
    assert shown["inside"]


def test_labels_as_stata():
    """The labels of the bars as Stata draws them for the same command (fixothermargins)."""
    res = fit("fixrev")
    shown = read(res.makeplot("waterfall"), horiz=False)
    assert [shown["labels"][n] for n in ("marg_expl:a", "marg_expl:b", "marg_expl:x3", "overall:explained",
                                         "overall:unexplained", "overall:difference")] == \
        ["-.778 (24%)", "-.93 (28%)", "-.575 (17%)", "-2.28 (69%)", "-1.01 (31%)", "-3.29"]
    shown = read(res.makeplot("waterfalllevels"), horiz=False)
    assert [shown["labels"][n] for n in ("overall:focal_group", "marg_expl:a", "marg_expl:b",
                                         "marg_expl:x3", "overall:raked_focal", "overall:unexplained",
                                         "overall:ref_group")] == \
        ["3.465", "+.778 (24%)", "+.93 (28%)", "+.575 (17%)", "5.748", "+1.01 (31%)", "6.758"]
    assert shown["ticks"] == ["Group zero", "a", "b", "x3", "Group zero, raked", "Unexpl.", "Group one"]
    assert "From the Group zero mean to the Group one mean" in shown["texts"]
    assert "Outcome: mean" in shown["texts"]
    # the axis of Stata's graph: ticks 3 to 8, the means' bars from 3
    assert shown["vlim"] == (3.0, 8.0)
    assert all(near(shown["bars"][n][1], 3.0) for n in ("overall:focal_group", "overall:raked_focal",
                                                         "overall:ref_group"))
    shown = read(res.makeplot("stages"), horiz=False)
    assert [shown["labels"][n] for n in ("cum_expl:a", "cum_unexpl:a", "cum_expl:b", "cum_unexpl:b",
                                         "cum_expl:x3", "cum_unexpl:x3")] == \
        ["24%", "76%", "52%", "48%", "69%", "31%"]
    assert "Total diff (-3.29)" in shown["texts"]
    shown = read(res.makeplot("waterfall", reverse=True, swapaxis=True,
                              blocklabels={"a": "First block", "b": "Second block"}), horiz=True)
    assert [shown["labels"][n] for n in ("rev_marg_expl:x3", "rev_marg_expl:b", "rev_marg_expl:a")] == \
        ["-.548 (17%)", "-.77 (23%)", "-.965 (29%)"]
    assert shown["ticks"] == ["x3", "Second block", "First block", "Expl.", "Unexpl.", "Total diff"]
    assert "Group zero - Group one, reverse order" in shown["texts"]


def test_levels_sequentially():
    res = fit("seq")
    for swap in (False, True):
        shown = read(res.makeplot("waterfalllevels", swapaxis=swap), horiz=swap)
        check_bars(shown, levels_rows(res, ["a", "b", "x3"]))
        # the three means start at the same base: the end of the axis nearest 0 (the means are
        # positive), else 0
        starts = [shown["bars"][n][1]
                  for n in ("overall:focal_group", "overall:raked_focal", "overall:ref_group")]
        assert starts[0] == starts[1] == starts[2]
        assert near(starts[0], min(shown["vlim"])) or starts[0] == 0
        assert shown["inside"]
    assert shown["ticks"][0] == "Group zero" and shown["ticks"][4] == "Group zero, raked"


def test_marginal_and_stages_sequentially():
    res = fit("seq")
    shown = read(res.makeplot("marginal"), horiz=True)
    check_bars(shown, marginal_rows(res, ["a", "b", "x3"]))
    assert shown["ticks"] == ["Total diff", "Expl.", "a", "b", "x3", "Unexpl."]
    assert "Marginal explanatory contribution at each raking stage" in shown["texts"]
    assert shown["inside"]
    shown = read(res.makeplot("marginal", swapaxis=True), horiz=False)
    check_bars(shown, marginal_rows(res, ["a", "b", "x3"]))
    for swap in (False, True):
        shown = read(res.makeplot("stages", swapaxis=swap), horiz=swap)
        check_bars(shown, stages_rows(res, ["a", "b", "x3"]))
        assert shown["ticks"] == ["a", "+ b", "+ x3"]
        xs, ys = shown["lines"]["total"]
        at = xs if swap else ys
        assert near(at[0], res.b["overall:difference"]) and near(at[1], res.b["overall:difference"])
        assert "Explained and unexplained difference by raking stage" in shown["texts"]
        assert shown["inside"]
    # the cumulative parts in the stage table are those drawn
    S = res.matrices["bystage"]
    for bn in ("a", "b", "x3"):
        assert near(S.loc[bn, "cum_expl"], res.b["cum_expl:" + bn])


def test_reverse_fixothermargins():
    res = fit("fixrev")
    rb = ["x3", "b", "a"]
    shown = read(res.makeplot("waterfall", reverse=True), horiz=False)
    check_bars(shown, waterfall_rows(res, rb, pfx="rev_"))
    shown = read(res.makeplot("waterfalllevels", reverse=True), horiz=False)
    check_bars(shown, levels_rows(res, rb, pfx="rev_"))
    assert "Reverse order" in shown["texts"]
    shown = read(res.makeplot("marginal", reverse=True), horiz=True)
    check_bars(shown, marginal_rows(res, rb, pfx="rev_"))
    shown = read(res.makeplot("stages", reverse=True), horiz=False)
    check_bars(shown, stages_rows(res, rb, pfx="rev_"))
    assert shown["ticks"] == ["x3", "+ b", "+ a"]
    # without reverse the forward order
    shown = read(res.makeplot("stages"), horiz=False)
    check_bars(shown, stages_rows(res, ["a", "b", "x3"]))


def test_individually():
    res = fit("ind")
    blocks = ["a", "b", "x3", "_remainder"]
    shown = read(res.makeplot(), horiz=True)
    check_bars(shown, marginal_rows(res, blocks))
    assert shown["ticks"] == ["Total diff", "Expl.", "a", "b", "x3", "Remainder", "Unexpl."]
    # Stata's graph of the same command
    assert [shown["labels"][n] for n in ("overall:difference", "overall:explained", "marg_expl:a",
                                         "marg_expl:b", "marg_expl:x3", "marg_expl:_remainder",
                                         "overall:unexplained")] == \
        ["-3.29", "-2.28 (69%)", "-1.2 (36%)", "-1.05 (32%)", "-.85 (26%)", ".816 (-25%)", "-1.01 (31%)"]
    shown = read(res.makeplot("marginal", swapaxis=True, blocklabels='b "Second block"'), horiz=False)
    check_bars(shown, marginal_rows(res, blocks))
    assert shown["ticks"][3] == "Second block"


def test_outcomes():
    res = fit("outcomes")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        shown = read(res.makeplot("waterfall"), horiz=False)
    assert "outcome y, the first of the outcomes" in buf.getvalue()
    check_bars(shown, waterfall_rows(res, ["a", "b"], pfx="y_", eq0="y"))
    shown = read(res.makeplot("waterfalllevels", outcome="earn", format="%9.2f"), horiz=False)
    check_bars(shown, levels_rows(res, ["a", "b"], pfx="earn_", eq0="earn"))
    b = res.b
    assert shown["labels"]["earn:focal_group"] == "%.2f" % b["earn:focal_group"]
    share = 100 * b["earn_marg_expl:a"] / b["earn:difference"]
    assert shown["labels"]["earn_marg_expl:a"] == "+%.2f (%d%%)" % (-b["earn_marg_expl:a"], int(share + 0.5))
    assert "earn: mean" in shown["texts"]
    shown = read(res.makeplot("stages", outcome="ybin"), horiz=False)
    check_bars(shown, stages_rows(res, ["a", "b"], pfx="ybin_"))
    shown = read(res.makeplot("marginal", outcome="earn"), horiz=True)
    check_bars(shown, marginal_rows(res, ["a", "b"], pfx="earn_", eq0="earn"))


def test_chain():
    res = fit("chain")
    shown = read(res.makeplot(), horiz=False)
    check_bars(shown, waterfall_rows(res, ["a", "b"], pfx="total_", eq0="total"))
    # Stata's graph of the same command
    assert [shown["labels"][n] for n in ("total_marg_expl:a", "total_marg_expl:b", "total:explained",
                                         "total:unexplained", "total:difference")] == \
        ["-.173 (-213%)", "+.135 (167%)", "-.0377 (-47%)", "+.119 (147%)", ".0809"]
    assert "(g3 = 0) - (g3 = 2) (sum of 2 pairs)" in shown["texts"]
    shown = read(res.makeplot("marginal"), horiz=True)
    check_bars(shown, marginal_rows(res, ["a", "b"], pfx="total_", eq0="total"))
    shown = read(res.makeplot("waterfalllevels", chainstep=2), horiz=False)
    check_bars(shown, levels_rows(res, ["a", "b"], pfx="p2_", eq0="p2"))
    assert shown["ticks"] == ["g3 = 1", "a", "b", "g3 = 1\nraked", "Unexpl.", "g3 = 2"]
    assert [shown["labels"][n] for n in ("p2:focal_group", "p2_marg_expl:a", "p2_marg_expl:b",
                                         "p2:raked_focal", "p2:unexplained", "p2:ref_group")] == \
        ["5.365", "+.0637 (-44%)", "-.144 (100%)", "5.285", "-.0636 (44%)", "5.222"]
    assert shown["vlim"] == (5.2, 5.45)
    assert "From the (g3 = 1) mean to the (g3 = 2) mean" in shown["texts"]
    shown = read(res.makeplot("stages", chainstep=1, swapaxis=True), horiz=True)
    check_bars(shown, stages_rows(res, ["a", "b"], pfx="p1_"))
    assert [shown["labels"][n] for n in ("p1_cum_expl:a", "p1_cum_unexpl:a", "p1_cum_expl:b",
                                         "p1_cum_unexpl:b")] == ["173%", "-73%", "188%", "-88%"]
    assert "Total diff (-.0628)" in shown["texts"]
    assert "(g3 = 0) - (g3 = 1)" in shown["texts"]
    shown = read(res.makeplot("waterfall", chainstep=1), horiz=False)
    check_bars(shown, waterfall_rows(res, ["a", "b"], pfx="p1_", eq0="p1"))


# ---- intervals -----------------------------------------------------------------------------------

def test_intervals():
    res = fit("b03")
    rb = ["b", "a"]
    b = res.b
    # the default type: the first of the results' types (percentile bc)
    shown = read(res.makeplot("marginal", ci=True), horiz=True)
    assert "95% percentile intervals" in shown["texts"]
    for name, pos, start, length in marginal_rows(res, ["a", "b"]):
        p, lo, hi = shown["ci"][name]
        ll, ul = ci_of(res, "percentile", name)
        assert near(p, pos) and near(lo, ll) and near(hi, ul), name
    # stages: the intervals of the cumulative parts
    shown = read(res.makeplot("stages", ci=True, citype="bc", reverse=True), horiz=False)
    assert "95% bias-corrected intervals" in shown["texts"]
    for name, pos, start, length in stages_rows(res, rb, pfx="rev_"):
        p, lo, hi = shown["ci"][name]
        ll, ul = ci_of(res, "bc", name)
        assert near(p, pos) and near(lo, ll) and near(hi, ul), name
    # waterfall: the interval of a bar's length, around its end
    shown = read(res.makeplot("waterfall", ci=True, normal=True), horiz=False)
    assert "95% normal intervals" in shown["texts"]
    for name, pos, start, length in waterfall_rows(res, ["a", "b"]):
        p, lo, hi = shown["ci"][name]
        ll, ul = ci_of(res, "normal", name)
        end = start + length
        assert near(lo, end - b[name] + ll) and near(hi, end - b[name] + ul), name
    # waterfalllevels: the means' intervals around them; a step is minus the value in e(b),
    # its interval reversed, around the step's end
    shown = read(res.makeplot("waterfalllevels", ci="tinterval", swapaxis=True), horiz=True)
    assert "95% t intervals" in shown["texts"]
    for name, pos, start, length in levels_rows(res, ["a", "b"]):
        p, lo, hi = shown["ci"][name]
        ll, ul = ci_of(res, "t", name)
        if start is None:
            assert near(lo, ll) and near(hi, ul), name
        else:
            end = start + length
            assert near(lo, end + b[name] - ul) and near(hi, end + b[name] - ll), name
        assert near(p, pos)
    assert shown["inside"]


def test_intervals_at_another_level():
    # level=: the intervals recomputed from the replicates, as res.estat_bootstrap() gives them
    # after res.display(level=90)
    import contextlib
    import io
    res = fit("b03")
    shown = read(res.makeplot("marginal", ci=True, level=90), horiz=True)
    assert "90% percentile intervals" in shown["texts"]
    with contextlib.redirect_stdout(io.StringIO()):
        res.display(level=90)
        at90 = res.estat_bootstrap()["ci_percentile"]
    for name, pos, start, length in marginal_rows(res, ["a", "b"]):
        p, lo, hi = shown["ci"][name]
        assert near(lo, float(at90.loc["ll", name])) and near(hi, float(at90.loc["ul", name])), name
        assert lo > ci_of(res, "percentile", name)[0] - 1e-12
    refused(lambda: res.makeplot("marginal", level=90))
    refused(lambda: res.makeplot("marginal", ci=True, level=5))


def test_intervals_bca():
    res = fit("b10")
    shown = read(res.makeplot(ci=True, bca=True), horiz=True)
    assert "95% BCa intervals" in shown["texts"]
    for name, pos, start, length in marginal_rows(res, ["a", "b", "_remainder"]):
        p, lo, hi = shown["ci"][name]
        ll, ul = ci_of(res, "bca", name)
        assert near(lo, ll) and near(hi, ul), name
    res = fit("b05")
    shown = read(res.makeplot("waterfalllevels", chainstep=1, ci=True), horiz=False)
    assert "95% BCa intervals" in shown["texts"]
    for name, pos, start, length in levels_rows(res, ["a", "b"], pfx="p1_", eq0="p1"):
        assert name in shown["ci"], name
    # without ci, no intervals and no note
    shown = read(res.makeplot("waterfall"), horiz=False)
    assert not shown["ci"] and not any("intervals" in t for t in shown["texts"])


# ---- options -------------------------------------------------------------------------------------

def test_levels_base():
    """The means' bars of the waterfall in levels start at 0 when the levels straddle it, else
    at the end of the axis nearest 0."""
    df = data()
    df["yneg"] = -df["y"]
    df["ymid"] = df["y"] - 5
    for v, where in (("yneg", "high"), ("ymid", "zero")):
        res = sd.seldecomp(df, "seldecomp %s (a: x1) (b: x2 x4) [pw=wt], by(group) sequentially" % v,
                           labels=LABELS, value_labels=VALUE_LABELS)
        for swap in (False, True):
            shown = read(res.makeplot("waterfalllevels", swapaxis=swap), horiz=swap)
            check_bars(shown, levels_rows(res, ["a", "b"]))
            starts = [shown["bars"][n][1] for n in ("overall:focal_group", "overall:raked_focal",
                                                    "overall:ref_group")]
            assert starts[0] == starts[1] == starts[2]
            assert starts[0] == (max(shown["vlim"]) if where == "high" else 0.0)
            assert shown["inside"]


def test_blocklabels_format_varlabels():
    res = fit("seq")
    shown = read(res.makeplot("waterfall", blocklabels='a "Group A" x3 Age'), horiz=False)
    assert shown["ticks"][:3] == ["Group A", "b", "Age"]
    shown = read(res.makeplot("stages", blocklabels={"b": "Second"}), horiz=False)
    assert shown["ticks"] == ["a", "+ Second", "+ x3"]
    shown = read(res.makeplot("marginal", format="%9.2f"), horiz=True)
    b = res.b
    D = b["overall:difference"]
    assert shown["labels"]["overall:difference"] == "%.2f" % D
    v = b["marg_expl:b"]
    assert shown["labels"]["marg_expl:b"] == "%.2f (%d%%)" % (v, int(100 * v / D + 0.5))
    shown = read(res.makeplot("waterfall", format="%9,3f"), horiz=False)
    assert shown["labels"]["overall:difference"] == ("%.3f" % D).replace(".", ",")
    res = fit("terms")
    shown = read(res.makeplot("waterfall"), horiz=False)
    assert shown["ticks"][:3] == ["x1", "x2", "x4"]
    shown = read(res.makeplot("waterfall", varlabels=True), horiz=False)
    assert shown["ticks"][:3] == ["First factor", "Second factor", "x4"]
    shown = read(res.makeplot("waterfall", varlabels=True, blocklabels={"x1": "One"}), horiz=False)
    assert shown["ticks"][:2] == ["One", "Second factor"]
    # an interaction block: its name without varlabels, its variables' labels with them
    res = sd.seldecomp(data(), "seldecomp y x1#x4 x2 [pw=wt], by(group) sequentially", labels=LABELS,
                       value_labels=VALUE_LABELS)
    assert read(res.makeplot(), horiz=False)["ticks"][:2] == ["x1_x4", "x2"]
    ticks = read(res.makeplot(varlabels=True), horiz=False)["ticks"]
    assert ticks[:2] == ["First factor # x4", "Second factor"]


def test_number_formats():
    """The labels' numbers as Stata's string(x, fmt) shows them (the values Stata gives)."""
    from seldecomp._plot import sfmt
    for fmt, x, want in [
            ("%9.3g", -0.778393, "-.778"), ("%9.3g", -0.929566, "-.93"), ("%9.3g", 0.0809, ".0809"),
            ("%9.3g", 1234.5, "1235"), ("%9.3g", 0.0012345, ".00123"), ("%9.3g", 99.95, "100"),
            ("%9.2g", 99.95, "1.0e+02"), ("%9.1g", 1.2345, "1.2"), ("%9.3g", 1e-5, ".00001"),
            ("%9.3g", -2.5e-7, "-2.50e-07"), ("%10.5g", -2.5e-7, "-2.500e-07"), ("%8.4g", 1e15, "1.0e+15"),
            ("%9.3g", 3.14159e9, "3.14e+09"), ("%9.3g", 0.9995, "1"), ("%9.3gc", 12345.6, "12,346"),
            ("%9.2f", 0.5, "0.50"), ("%9.2f", -0.0004, "-0.00"), ("%9.0f", 0.5, "0"),
            ("%9.0f", 1234.5, "1234"),
            ("%9.3f", 123456.7, "1.23e+05"), ("%9.1fc", 1234.5, "1,234.5"), ("%9,2f", 1.2345, "1,23"),
            ("%9.2e", -0.0012345, "-1.23e-03"), ("%9.2f", 3.14159e9, "3.14e+09")]:
        assert sfmt(x, fmt).strip() == want, (fmt, x, sfmt(x, fmt), want)


def test_matplotlib_arguments():
    res = fit("seq")
    fig, ax = plt.subplots()
    out = res.makeplot("marginal", ax=ax, title="My title", xlim=(-5, 1))
    assert out is fig
    assert ax.get_title().startswith("My title\n") and "Group zero - Group one" in ax.get_title()
    assert ax.get_xlim() == (-5.0, 1.0)
    plt.close(fig)
    fig = res.makeplot("stages", figsize=(9, 4), title="", subtitle="Sub")
    assert tuple(fig.get_size_inches()) == (9.0, 4.0)
    texts = [t.get_text() for t in fig.findobj(Text)]
    assert "Sub" in texts and "Explained and unexplained difference by raking stage" not in texts
    plt.close(fig)
    fig = res.makeplot("waterfall", twoway={"ylabel": "Gap", "figsize": (6, 4)})
    assert fig.axes[0].get_ylabel() == "Gap"
    plt.close(fig)
    # in Axes of the user's, with intervals: the note below the Axes
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    boot = fit("b03")
    boot.makeplot("waterfalllevels", ax=axes[1], ci=True)
    fig.canvas.draw()
    note = [t for t in axes[1].texts if t.get_gid() == "note"][0]
    assert note.get_text() == "95% percentile intervals"
    assert note.get_window_extent().y1 <= axes[1].get_window_extent().y0
    assert any((c.get_gid() or "").startswith("ci:") for c in axes[1].collections)
    assert not axes[0].patches
    plt.close(fig)


def test_makeplot_option_of_seldecomp():
    plt.close("all")
    res = sd.seldecomp(data(), "seldecomp " + SPECS["seq"].replace(
        "sequentially", 'sequentially makeplot(waterfalllevels, blocklabels(a "First") format(%9.2f))'),
        labels=LABELS, value_labels=VALUE_LABELS)
    fig = plt.gcf()
    shown = read(fig, horiz=False)
    check_bars(shown, levels_rows(res, ["a", "b", "x3"]))
    assert shown["ticks"][1] == "First"
    assert shown["labels"]["overall:ref_group"] == "%.2f" % res.b["overall:ref_group"]
    res = run_case("b03", makeplot="marginal, ci tinterval")[2]
    shown = read(plt.gcf(), horiz=True)
    assert "95% t intervals" in shown["texts"]
    p, lo, hi = shown["ci"]["overall:explained"]
    assert (lo, hi) == ci_of(res, "t", "overall:explained")


# ---- refusals ------------------------------------------------------------------------------------

def test_refusals():
    seq, ind, chain, outs, boot = fit("seq"), fit("ind"), fit("chain"), fit("outcomes"), fit("b03")
    refused(lambda: fit("plain").makeplot(),
            "makeplot requires results estimated with sequentially or individually")
    refused(lambda: sd.seldecomp(data(), "seldecomp " + SPECS["plain"] + " makeplot"),
            "makeplot requires sequentially or individually")
    refused(lambda: seq.makeplot("pie"), "plot type pie not allowed")
    refused(lambda: ind.makeplot("waterfall"),
            "makeplot(waterfall) not allowed with individually; use marginal")
    refused(lambda: ind.makeplot("stages"), "makeplot(stages) not allowed with individually")
    refused(lambda: seq.makeplot("waterfall", reverse=True),
            "makeplot(reverse) requires results estimated with sequentially reverse")
    refused(lambda: seq.makeplot("waterfall", ci=True),
            "makeplot(ci): these results have no confidence intervals")
    refused(lambda: boot.makeplot("waterfall", ci=True, bca=True),
            "makeplot(bca): these results have no bca intervals")
    refused(lambda: boot.makeplot("waterfall", ci=True, normal=True, bc=True),
            "specify only one interval type")
    refused(lambda: boot.makeplot("waterfall", percentile=True), "specify ci with an interval type")
    refused(lambda: boot.makeplot("waterfall", ci="wide"), "wide is not an interval type")
    refused(lambda: seq.makeplot("waterfall", chainstep=1),
            "makeplot(chainstep()) requires results from a chain run")
    refused(lambda: chain.makeplot("waterfall", chainstep=3), "makeplot(chainstep()) must be 1 to 2")
    refused(lambda: chain.makeplot("waterfall", chainstep=0), "must be a positive integer")
    refused(lambda: chain.makeplot("waterfall", chainstep="x"), "must be a positive integer")
    refused(lambda: chain.makeplot("stages"),
            "makeplot(stages) requires chainstep(#) with chain results")
    refused(lambda: chain.makeplot("waterfalllevels"), "makeplot(waterfalllevels) requires chainstep(#)")
    refused(lambda: outs.makeplot("waterfall", outcome="hrs"),
            "makeplot(outcome()): hrs is not one of the outcomes (y ybin earn)")
    refused(lambda: outs.makeplot("waterfall", outcome="y earn"), "y earn names 2 outcomes; name one")
    refused(lambda: seq.makeplot("waterfall", outcome="earn"),
            "makeplot(outcome()): earn is not the outcome (y)")
    refused(lambda: seq.makeplot("waterfall", format="%9.2x"), "%9.2x is not a numeric display format")
    refused(lambda: seq.makeplot("waterfall", format="%3.3g"), "is not a numeric display format")
    refused(lambda: seq.makeplot("waterfall", blocklabels={"z": "Zed"}),
            "makeplot(blocklabels()): z is not a block; the blocks are a b x3")
    refused(lambda: seq.makeplot("waterfall", blocklabels={"a": ""}), "the text of a is empty")
    refused(lambda: seq.makeplot("waterfall", blocklabels='a "A" a "B"'), "block a is given twice")
    refused(lambda: seq.makeplot("waterfall", blocklabels='a "A" b'),
            "type pairs of a block name and its text")
    refused(lambda: seq.makeplot("waterfall", swapaxes=True), "option swapaxes not allowed")
    refused(lambda: seq.makeplot("waterfall", twoway="title(x)"), "makeplot(twoway()) is not available")
    # the refusals of makeplot() at estimation, made before anything is estimated
    cmd = "seldecomp " + SPECS["seq"]
    for spec, text in (("makeplot(waterfall, ci)", "makeplot(ci) requires bayesboot()"),
                       ("makeplot(waterfall, reverse)", "makeplot(reverse) requires sequentially reverse"),
                       ("makeplot(waterfall, chainstep(1))", "makeplot(chainstep()) requires chain"),
                       ("makeplot(stages, percentile)", "makeplot(): specify ci with an interval type"),
                       ("makeplot(stages, ci bc normal)", "makeplot(): specify only one interval type"),
                       ("makeplot(waterfall, outcome(x1))", "makeplot(outcome()): x1 is not the outcome (y)"),
                       ("makeplot(waterfall, blocklabels(zz Text))", "makeplot(blocklabels()): zz is not a block"),
                       ("makeplot(waterfall, format(abc))", "is not a numeric display format"),
                       ("makeplot(waterfall, twoway(title(x)))", "twoway() takes Stata graph options")):
        refused(lambda spec=spec: sd.seldecomp(data(), cmd + " " + spec), text)
    refused(lambda: sd.seldecomp(data(), cmd + " bayesboot(reps(5) rseed(1)) makeplot(waterfall, ci bca)"),
            "makeplot(bca): BCa intervals require bayesboot(bca)")
    refused(lambda: sd.seldecomp(data(), "seldecomp " + SPECS["ind"] + " makeplot(stages)"),
            "makeplot(stages) not allowed with individually; use marginal")


if __name__ == "__main__":
    for _name in sorted(dir()):
        if _name.startswith("test_"):
            globals()[_name]()
    print("ok")

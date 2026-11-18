"""The printed output against Stata's: what the command prints at estimation, at a replay
with display options, and what estat bootstrap prints, for every case of fixtures/cases.txt."""
import contextlib
import io
import os
import re

import common
import seldecomp as sd
from seldecomp import _display, _parse


def norm(lines):
    """The lines without trailing blanks and without the blank lines at either end."""
    out = [l.rstrip() for l in lines]
    while out and out[0] == "":
        out.pop(0)
    while out and out[-1] == "":
        out.pop()
    return out


def cases_in_order():
    """The case names in the order of cases.txt: a follow-up (replay, estat bootstrap) acts on
    the results as the follow-ups before it left them."""
    names = []
    with open(os.path.join(common.FIXTURES, "cases.txt")) as f:
        for line in f:
            line = line.strip().lstrip("+").strip()
            if line and not line.startswith("*"):
                names.append(line.split(":")[0].strip())
    return names


def follow_text(res, cmd):
    """The Stata command typed after the estimation -> the text of the matching call."""
    cmd = cmd.strip()
    k = cmd.find(",")
    head, opts = (cmd, "") if k < 0 else (cmd[:k].strip(), cmd[k + 1:])
    kw = {}
    for tok in _parse.split_top(opts):
        m = re.match(r"^([A-Za-z_]+)(\((.*)\))?$", tok)
        kw[m.group(1)] = True if m.group(2) is None else m.group(3)
    if head == "seldecomp":
        return _display.replay_text(res, **kw)
    return _display.estat_text(res, **kw)


def test_output_as_stata():
    results = {}
    checked = 0
    for name in cases_in_order():
        case = common.load_case(name)
        if "_" in name:
            text = follow_text(results[name.split("_")[0]], case["cmd"])
        else:
            if case["rc"] != 0:
                continue
            with contextlib.redirect_stdout(io.StringIO()):
                res = common.run_case(name)[2]
            results[name] = res
            # the trace of noisily is printed as the raking runs, ahead of the output
            text = "\n".join(res.trace + [str(res)])
        assert norm(text.split("\n")) == norm(case["output"]), (name, case["cmd"])
        checked += 1
    assert checked >= 110


def test_display_and_estat_return():
    """display() prints the replay and returns nothing; estat_bootstrap() prints its table and
    returns the intervals it shows."""
    with contextlib.redirect_stdout(io.StringIO()) as out:
        res = common.run_case("b02")[2]
        assert res.display(level=90) is None
        r = res.estat_bootstrap(all=True)
    assert "[90% conf. interval]" in out.getvalue()
    assert r["level"] == 90
    for key in ("ci_normal", "ci_percentile", "ci_bc", "ci_bca", "pct_explained_ci"):
        assert key in r, key
    assert list(r["ci_bca"].index) == ["ll", "ul"]
    assert res.e["level"] == 95

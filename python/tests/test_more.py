"""More cases against the Stata program (fixtures/more): commands that Stata refuses, with
their return code and message; stored results and generated weights; the printed output at
estimation, at a replay and of estat bootstrap."""
import re

import numpy as np

import common
import seldecomp as sd
from seldecomp import _display, _parse

_runs = {}


def run(entry):
    name, ds, cmd, _ = entry
    if name not in _runs:
        _runs[name] = common.run_more(name, ds, cmd)
    return _runs[name]


def follow_text(res, cmd):
    """The Stata follow-up command (a replay, estat bootstrap) as the matching Python call."""
    k = cmd.find(",")
    head, opts = (cmd, "") if k < 0 else (cmd[:k].strip(), cmd[k + 1:])
    kw = {}
    for tok in _parse.split_top(opts):
        m = re.match(r"^([A-Za-z_]+)(\((.*)\))?$", tok)
        kw[m.group(1)] = True if m.group(2) is None else m.group(3)
    if head == "seldecomp":
        return _display.replay_text(res, **kw)
    return _display.estat_text(res, **kw)


def test_more_refusals():
    """A command that Stata refuses is refused with the same return code, and the message is
    the end of Stata's output, after notes and warnings only."""
    bad = []
    n = 0
    for entry in common.more_cases():
        case, df, res, rc, msg, warned = run(entry)
        if rc != case["rc"]:
            bad.append("%s: rc %s, Stata %s" % (entry[0], rc, case["rc"]))
            continue
        if rc == 0:
            continue
        n += 1
        stata = common.squeeze(case["output"])
        text = " ".join(msg.split())
        before = stata[:-len(text)].strip() if stata.endswith(text) else None
        if before is None or not (before == "" or before.startswith("note:") or before.startswith("warning:")):
            bad.append("%s: message %r, Stata %r" % (entry[0], text[:120], stata[-160:]))
        elif ("warning: " in stata) != bool(warned):
            bad.append("%s: the warning" % entry[0])
    assert n > 400 and not bad, bad[:10]


def test_more_results():
    """Every stored number, name and matrix, and the weights of generate()."""
    bad = []
    n = 0
    for entry in common.more_cases():
        case, df, res, rc, msg, warned = run(entry)
        if res is None or case["rc"] != 0:
            continue
        n += 1
        diffs = common.compare(res, case)
        if diffs:
            bad.append("%s: %s" % (entry[0], diffs[0]))
        for col, want in common.more_generated(entry[0]).items():
            got = df[col].to_numpy(dtype=float) if col in df.columns else None
            if got is None or not np.allclose(got, want, rtol=1e-9, atol=0, equal_nan=True):
                bad.append("%s: generated %s" % (entry[0], col))
    assert n > 300 and not bad, bad[:10]


def test_more_output():
    """The printed output, line by line: at estimation, and of each follow-up."""
    bad = []
    n = 0
    for entry in common.more_cases():
        case, df, res, rc, msg, warned = run(entry)
        if res is None or case["rc"] != 0:
            continue
        n += 1
        text = "\n".join(res.trace + [_display.render(res)])
        if common.trimmed(case["output"]) != common.trimmed(text.split("\n")):
            bad.append(entry[0])
        for fname, fcmd in entry[3]:
            want = common.more_stata(fname)
            if want is None:
                continue
            n += 1
            try:
                got, frc = follow_text(res, fcmd), 0
            except sd.SeldecompError as ex:
                got, frc = ex.message, ex.rc
            if frc != want["rc"]:
                bad.append("%s: rc %s, Stata %s" % (fname, frc, want["rc"]))
            elif frc:
                if not common.squeeze(want["output"]).endswith(" ".join(got.split())):
                    bad.append(fname)
            elif common.trimmed(want["output"]) != common.trimmed(got.split("\n")):
                bad.append(fname)
    assert n > 450 and not bad, bad[:10]

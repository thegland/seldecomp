"""Every stored result against the Stata program: the cases of fixtures/cases.txt, run by
make_fixtures.do in Stata and frozen as JSON."""
import numpy as np
import pandas as pd

import common
import seldecomp as sd


def test_cases_agree_with_stata():
    """Scalars, text results, matrices and notes of every case that runs."""
    checked = 0
    for name in common.case_names():
        case = common.load_case(name)
        if case["rc"] != 0:
            continue
        case, df, res = common.run_case(name)
        diffs = common.compare(res, case)
        assert not diffs, (name, case["cmd"], diffs[:5])
        checked += 1
    assert checked >= 85


def test_errors_agree_with_stata():
    """The return code and the message of every case that Stata refuses."""
    checked = 0
    for name in common.case_names():
        case = common.load_case(name)
        if case["rc"] == 0:
            continue
        df = common.load_data()
        try:
            sd.seldecomp(df, case["cmd"], value_labels=common.VALUE_LABELS)
        except sd.SeldecompError as ex:
            assert ex.rc == case["rc"], (name, ex.rc, case["rc"], ex.message)
            stata = " ".join(" ".join(l.split()) for l in case["output"] if l.strip())
            assert " ".join(ex.message.split()) in stata, (name, ex.message, stata[-300:])
            checked += 1
        else:
            raise AssertionError("%s: no error, Stata gives %d" % (name, case["rc"]))
    assert checked >= 20


def test_generated_weights():
    """generate(): the weights of every stage, observation by observation."""
    import os
    checked = 0
    for name in common.case_names():
        path = os.path.join(common.FIXTURES, name + "_weights.csv")
        if not os.path.exists(path):
            continue
        case, df, res = common.run_case(name)
        stata = pd.read_csv(path)
        for col in [c for c in stata.columns if c != "id"]:
            a = df[col].to_numpy(dtype=float)
            b = stata[col].to_numpy(dtype=float)
            assert (np.isnan(a) == np.isnan(b)).all(), (name, col)
            ok = ~np.isnan(a)
            assert np.allclose(a[ok], b[ok], rtol=1e-9, atol=1e-12), (name, col)
            assert np.allclose(res.weights[col.replace("rk", "rw")].to_numpy()[ok], b[ok], rtol=1e-9, atol=1e-12)
        checked += 1
    assert checked >= 10

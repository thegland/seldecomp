"""The two forms of the call, the Python-only arguments, and the generated columns."""
import numpy as np
import pandas as pd

import common
import seldecomp as sd


def same(a, b):
    return list(a.b.index) == list(b.b.index) and np.allclose(a.b.to_numpy(float), b.b.to_numpy(float), rtol=1e-12, atol=0)


def test_keywords_equal_command_text():
    df = common.load_data()
    a = sd.seldecomp(df, "seldecomp y (a: x1) (b: x2 x4) c.x3 [pw=wt], by(group) contvars(x3 5) "
                         "sequentially fixothermargins reverse detail")
    b = sd.seldecomp(df, "y (a: x1) (b: x2 x4) c.x3", by="group", weights="wt", contvars="x3 5",
                     sequentially=True, fixothermargins=True, reverse=True, detail=True)
    c = sd.seldecomp(df, "y (a: x1) (b: x2 x4) c.x3", by="group", weights="wt", contvars={"x3": 5},
                     seq=True, fix=True, rev=True, det=True)
    assert same(a, b) and same(a, c)
    assert a.e["stagetype"] == "sequentially" and b.e["contvars"] == "x3 5"


def test_suboptions_as_text_or_dict():
    df = common.load_data()
    a = sd.seldecomp(df, "y x1 c.x3", by="group", weights="wt", contvars="x3 4, focalgroup")
    b = sd.seldecomp(df, "y x1 c.x3", by="group", weights="wt", contvars={"x3": 4, "sample": "focalgroup"})
    assert same(a, b)
    small = df[df["id"] <= 300].copy()
    c = sd.seldecomp(small, "y x1 x2", by="group", weights="wt", bayesboot="reps(20) rseed(5) cluster(cluster) bca")
    d = sd.seldecomp(small, "y x1 x2", by="group", weights="wt",
                     bayesboot={"reps": 20, "rseed": 5, "cluster": "cluster", "bca": True})
    assert np.allclose(c.V.to_numpy(), d.V.to_numpy()) and c.e["citypes"] == d.e["citypes"] == "bca"


def test_subset_and_if():
    df = common.load_data()
    a = sd.seldecomp(df, 'y x1 x2 [pw=wt] if x4 != "east" in 101/1800, by(group)')
    mask = (df["x4"] != "east").to_numpy() & (np.arange(len(df)) >= 100) & (np.arange(len(df)) < 1800)
    b = sd.seldecomp(df, "y x1 x2", by="group", weights="wt", subset=mask)
    assert same(a, b) and a.e["N"] == b.e["N"] == int(mask.sum())
    assert bool((b.sample.to_numpy() == mask).all())


def test_categorical_columns():
    df = common.load_data()
    a = sd.seldecomp(df, "y x1 x2 x4", by="group", weights="wt")
    d2 = df.copy()
    d2["group"] = pd.Categorical(np.where(df["group"] == 0, "zero", "one"), categories=["zero", "one"])
    d2["x4"] = pd.Categorical(d2["x4"], categories=["west", "south", "north", "east"])
    b = sd.seldecomp(d2, "y x1 x2 x4", by="group", weights="wt")
    assert np.allclose(a.b.to_numpy(float), b.b.to_numpy(float), rtol=1e-9)
    assert b.e["focal_group"] == "zero" and b.e["levels3"] == "west | south | north | east"
    c = sd.seldecomp(d2, "y x1 x2 x4", by="group", weights="wt", focalgroup="one")
    assert c.e["focal_group"] == "one" and abs(c.e["difference"] + b.e["difference"]) < 1e-12


def test_warning_is_issued():
    # Stata prints "warning: ..." before the error; here it is a SeldecompWarning
    import warnings
    case = common.load_case("e16")
    df = common.load_data()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            sd.seldecomp(df, case["cmd"])
            raise AssertionError("no error")
        except sd.SeldecompError as ex:
            assert ex.rc == case["rc"] == 430
    assert len(caught) == 1 and caught[0].category is sd.SeldecompWarning
    assert caught[0].filename == __file__
    text = " ".join(" ".join(case["output"]).split())
    assert "warning: " + str(caught[0].message) in text


def test_saving_writes_csv():
    # bayesboot(saving(file)): a CSV file, .csv added to a name without an extension
    import os
    import tempfile
    df = common.load_data()
    small = df[df["id"] <= 300].copy()
    with tempfile.TemporaryDirectory() as folder:
        stem = os.path.join(folder, "draws")
        text = "reps(5) rseed(1) saving(%s%s)"
        res = sd.seldecomp(small, "y x1 x2", by="group", weights="wt", bayesboot=text % (stem, ""))
        assert res.e["bsfile"] == stem + ".csv" and os.path.exists(stem + ".csv")
        saved = pd.read_csv(stem + ".csv")
        assert saved.shape[0] == 5 and "overall:explained" in saved.columns
        try:
            sd.seldecomp(small, "y x1 x2", by="group", weights="wt", bayesboot=text % (stem, ""))
            raise AssertionError("no error")
        except sd.SeldecompError as ex:
            assert ex.rc == 602
        sd.seldecomp(small, "y x1 x2", by="group", weights="wt", bayesboot=text % (stem, ", replace"))


def test_version():
    assert sd.__version__ == "1.0.0"


def test_generate_adds_columns():
    df = common.load_data()
    res = sd.seldecomp(df, "y (a: x1) (b: x2)", by="group", weights="wt", sequentially=True, generate="rk")
    assert "rk" in df.columns and "rk_s1" in df.columns
    ref = (df["group"] == 1).to_numpy()
    assert np.allclose(df["rk"].to_numpy()[ref], df["wt"].to_numpy()[ref])
    focal = ~ref
    assert abs(df["rk"].to_numpy()[focal].sum() - df["wt"].to_numpy()[focal].sum()) < 1e-6
    try:
        sd.seldecomp(df, "y x1 x2", by="group", weights="wt", generate="rk")
    except sd.SeldecompError as ex:
        assert ex.rc == 110
    else:
        raise AssertionError("an existing variable was overwritten without replace")
    sd.seldecomp(df, "y x1 x2", by="group", weights="wt", generate="rk, replace")
    assert res.e["generate"] == "rk"


def test_generate_into_a_slice():
    # the columns go into the frame given, not into the frame it was taken from, without a warning
    import warnings
    df = common.load_data()
    part = df[df["x1"] > 1]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        sd.seldecomp(part, "y x1 x2", by="group", generate="rk", bayesboot="reps(3) rseed(1) generate(bw)")
    assert not [w for w in caught if w.category.__name__ == "SettingWithCopyWarning"]
    assert {"rk", "bw_1", "bw_3"} <= set(part.columns)
    assert not ({"rk", "bw_1"} & set(df.columns))


def test_invalid_calls():
    df = common.load_data()
    for kw, rc in (({"nosuchoption": True}, 198), ({"by": "group", "level": 90}, 198),
                   ({"by": "group", "trimfrequency": "often"}, 198),
                   ({"by": "group", "bayesboot": "200"}, 198)):
        try:
            sd.seldecomp(df, "y x1 x2", **kw)
        except sd.SeldecompError as ex:
            assert ex.rc == rc, (kw, ex.rc, ex.message)
        else:
            raise AssertionError("no error for %r" % (kw,))
    try:
        sd.seldecomp(df, "y x1 x2, by(group)", by="group")
    except sd.SeldecompError as ex:
        assert ex.rc == 198
    else:
        raise AssertionError("by() given twice was accepted")

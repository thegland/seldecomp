"""pack.py <folder> : gathers what make_fixtures.do wrote to <folder> into the two files that the
tests read: stata.json.gz (every case's return code, printed output and stored results, by
name) and weights.npz (the generated weights, "<case>:gen:<column>", and the replicate weights
of the bayesboot() cases, "<case>:bs", one row per observation of the sample, in order)."""
import csv
import gzip
import json
import os
import sys

import numpy as np

src = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
cases, arrays = {}, {}
for f in sorted(os.listdir(src)):
    path = os.path.join(src, f)
    if f.endswith(".json"):
        with open(path) as fh:
            cases[f[:-5]] = json.load(fh)
    elif f.endswith("_weights.csv") or f.endswith("_bsweights.csv"):
        with open(path) as fh:
            rows = list(csv.reader(fh))
        cols = rows[0]
        data = np.array([[float(x) if x not in ("", ".") else np.nan for x in r] for r in rows[1:]], dtype=float)
        if f.endswith("_bsweights.csv"):
            arrays[f[:-len("_bsweights.csv")] + ":bs"] = data[:, 1:]
        else:
            for j, c in enumerate(cols[1:]):
                arrays[f[:-len("_weights.csv")] + ":gen:" + c] = data[:, j + 1]
with gzip.GzipFile(os.path.join(here, "stata.json.gz"), "wb", mtime=0) as fh:
    fh.write(json.dumps(cases, sort_keys=True, separators=(",", ":")).encode("utf-8"))
np.savez_compressed(os.path.join(here, "weights.npz"), **arrays)
print("%d cases, %d arrays" % (len(cases), len(arrays)))

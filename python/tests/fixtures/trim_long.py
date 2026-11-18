"""trim_long.py <folder>: keep, of the long bootstrap runs written by make_fixtures.do for
cases_long.txt (s01.json, s02.json), what the tests read: the command, the scalars, and the
matrices b and se. The replicates themselves are dropped."""
import json
import os
import sys

folder = sys.argv[1]
for name in ("s01", "s02"):
    path = os.path.join(folder, name + ".json")
    with open(path) as f:
        d = json.load(f)
    out = {"case": d["case"], "cmd": d["cmd"], "rc": d["rc"], "scalars": d["scalars"],
           "matrices": {k: d["matrices"][k] for k in ("b", "se")}}
    with open(path, "w") as f:
        json.dump(out, f)
        f.write("\n")
    print(name, os.path.getsize(path), "bytes")

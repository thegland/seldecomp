"""run_tests.py [name ...]: run the tests without pytest. Every function test_* of every
file test_*.py in this folder is called; a name limits the run to the files or functions that
contain it. Exit status 1 if a test fails."""
import importlib
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))


def main(argv):
    want = argv[1:]
    failed = ran = 0
    for fn in sorted(os.listdir(HERE)):
        if not (fn.startswith("test_") and fn.endswith(".py")):
            continue
        mod = importlib.import_module(fn[:-3])
        for name in sorted(n for n in dir(mod) if n.startswith("test_")):
            if want and not any(w in fn or w in name for w in want):
                continue
            t0 = time.time()
            try:
                getattr(mod, name)()
                print("ok    %-22s %-42s %5.1fs" % (fn, name, time.time() - t0))
            except Exception:
                failed += 1
                print("FAIL  %-22s %s" % (fn, name))
                traceback.print_exc()
            ran += 1
    print("%d tests, %d failed" % (ran, failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

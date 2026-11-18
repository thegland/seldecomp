# Tests of the Python version

```
python3 python/tests/run_tests.py
```

`run_tests.py` needs no test framework and reads the package from `python/` (the functions
`test_*` also follow pytest's conventions). The data is `tests/data/seldecomp_test.csv`, the
data of the Stata program's tests.

The tests compare the Python version with the Stata program. `fixtures/` holds what Stata
produced for each case of `fixtures/cases.txt`: the stored results and the printed output
(`<case>.json`), the weights of `generate()` (`<case>_weights.csv`) and the replicate weights
that `bayesboot()` drew (`<case>_bsweights.csv.gz`), so the tests need no Stata.

| File | What it checks |
|---|---|
| `test_results.py` | Every stored number, name and note of every case; the return code and message of every case that is refused; the generated weights, observation by observation. |
| `test_display.py` | The printed output at estimation, at a replay with display options and of estat bootstrap, line by line. |
| `test_bootstrap.py` | The package's own bootstrap draws: reproducible with `rseed`, the replicate weights' distribution, standard errors against two long Stata runs. |
| `test_plot.py` | `makeplot`: the bars, intervals, labels and titles of every plot type against the results; the refusals. |
| `test_api.py` | Keyword form against the Stata command text, suboptions as text or dict, `subset=`, pandas categorical columns, `generate`, `saving`. |
| `test_more.py` | The cases of `fixtures/more`: commands that Stata refuses (return code and message), stored results and generated weights, printed output. |

For the bootstrap cases the Python version is given Stata's replicate weights (the random
number generators differ), so standard errors and intervals are compared exactly.

`fixtures/more` holds 858 more cases (512 of them commands that Stata refuses) on the test
data and on three small datasets, `fixtures/more/data` (`s1`: 600 rows with unusual columns;
`s2`: 40 rows; `s4`: labels and text with special characters). What Stata produced is in two
files: `stata.json.gz` (return codes, printed output, stored results) and `weights.npz`
(generated weights and replicate weights).

## Regenerating the fixtures (needs Stata)

`fixtures/make_fixtures.do` runs the cases and writes the files; its header says how to call
it (from a copy of the folder, never inside the repository). Afterwards compress the replicate
weights (`gzip -n *_bsweights.csv`). The two long runs of `fixtures/cases_long.txt` are run
the same way and trimmed with `fixtures/trim_long.py`.

`fixtures/more/make_fixtures.do` does the same for `fixtures/more/cases.txt` (run it from a
copy of the fixtures folder, in `more`); then `python3 pack.py <output folder>` writes the two
files.

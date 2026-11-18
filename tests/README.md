# seldecomp test suite

Stata tests for `seldecomp.ado`. They need Stata 16.1 or later; `seldecomp` itself needs no
other package (the raking is its own Mata code). Two SSC packages are optional: `oaxaca`
(Jann), needed only by `test_oaxaca.do`, and `ebalance` (Hainmueller and Xu), needed only by
`test_ebalance.do`. `run_tests.do` reports each of those two files as SKIPPED when its package
is missing, so a full pass needs both:

```stata
ssc install oaxaca      // test_oaxaca.do only
ssc install ebalance    // test_ebalance.do only
```

## Running

With the working directory set to this `tests` directory:

```stata
cd /path/to/seldecomp/tests
do run_tests.do
```

`run_tests.do` adds the parent of the working directory to the ado-path, so the `seldecomp.ado`
of the repository root is the one tested. To test a copy elsewhere, set the global
`SELDECOMP_ROOT` to its directory first; the data files are always read relative to the tests
directory, so the working directory must still be `tests`:

```stata
global SELDECOMP_ROOT "/path/to/seldecomp"
```

The tests write scratch files (logs, `.dta` files, graphs) into the working directory. To keep
them out of the repository, copy the repository to a scratch directory and run the suite
there, in batch mode if preferred (`stata-mp -b do run_tests.do`; the result is at the end of
`run_tests.log`). The suite takes about five minutes.

### The result

The runner executes every test file with `capture noisily`, so all files run even if one
fails. It ends with one of:

- `ALL SELDECOMP TESTS PASSED (26 files)`;
- `SELDECOMP TESTS PASSED (<n> of 26 files; <k> skipped)`, when `oaxaca` or `ebalance` is not
  installed (each skipped file is listed as `SKIPPED: <file> (<package> not installed)`);
- `FAILED FILES: <list>`, with return code 9.

Each test file can also be run on its own from the tests directory (`do test_options.do`):
each loads its data, defines its helper programs and ends with a line saying that it passed.
The comments in each file describe its checks one by one.

## Files

| File | What it checks |
|---|---|
| `run_tests.do` | The runner: looks for `oaxaca` and `ebalance`, sets the ado-path, runs the files below, prints the summary. |
| `test_reference.do` | Point estimates against the values of the Python reference (`expected_values.do`, `expected_support.do`): one and two categorical variables, interactions, weights, string variables, the staged decompositions, a level found in one group only, binned continuous variables. |
| `test_analytic.do` | The worked example of Pylypchuk and Selden (2008), footnote 18, against its closed-form solution, stage by stage. |
| `test_oracle.do` | 40 random configurations against an independent raking routine written in Mata from the description of the method: explained parts, generated weights and support failures. |
| `test_oaxaca.do` | Equality with the Oaxaca–Blinder decomposition (`oaxaca`, SSC) where the two methods coincide: one categorical variable, and a full interaction entered as one dummy per cell. For two additive categorical variables the methods differ; the test prints both. |
| `test_ebalance.do` | Cross-check against entropy balancing (`ebalance`, SSC) on the level indicators of the raking variables: the same reweighted mean and weights. |
| `test_options.do` | Options and messages: return codes and message texts of invalid commands, `focalgroup()`, weights (zero, missing, negative), `contvars()`, `generate()`, the raking options (`tolerance()`, `iterate()`, the trimming bounds), `force`, `insupport`, the `by` prefix, replay. |
| `test_bootstrap.do` | `bayesboot()`: stored results, reproducibility with a seed, `cluster()`, `saving()`, the interval types (normal, t, percentile, BC, BCa), failed replicates, `generate()` of replicate weights, the explained share, `estat bootstrap`. |
| `test_svy.do` | The replication prefixes (`svy brr`, `svy jackknife`, `svy sdr`, `svy bootstrap`, `bootstrap`, `jackknife`) on a synthetic stratified design: variances against hand computation, what the prefixes refuse, the check of replicates that will fail and its messages, replay. |
| `test_cutpoints.do` | The `contvars()` cut points in replicates: fixed at the full-sample values under every prefix and `bayesboot()` by default, recomputed with `recut`. |
| `test_plot.do` | `makeplot()`: every plot type, checked through the data of the drawn graph (bars, spikes, labels, axes), its options and its refusals. |
| `test_plot_intervals.do` | The interval types and levels drawn by `makeplot()`: graph coordinates against the coefficient intervals. |
| `test_display.do` | The printed output: formats, the header and its notes, the tables at estimation and replay and after each prefix, the reference transcripts (`reference/transcripts/`), one message per cause of failure. |
| `test_label_characters.do` | Labels with quotes and special characters in tables, replay and graphs; the explained share's standard error and interval. |
| `test_diagnostics_display.do` | `diagnostics` and the diagnostics shown automatically after a failure, their order, replay from saved results, the prefixes. |
| `test_table_labels.do` | Row labels and equation headings of the tables (ratios, several outcomes, `chain`, `estat bootstrap`) within 78 columns. |
| `test_diagnostic_counts.do` | The diagnostics' counts of targets met and of raking iterations, the scale of the weights, support exclusions; by stage and with `reverse`, `detail` and `chain`. |
| `test_diagnostic_bounds.do` | Counts of weights at the trimming bounds and the effective-sample-size ratio. |
| `test_failure_messages.do` | The cautions after failed bootstrap replicates and BCa jackknife deletions, at estimation, at replay and with `chain`; inference unchanged by them. |
| `test_invariance.do` | Invariance of the results to rescaled weights, a linear transformation of the outcome, repeated observations, the order of the rows, and recoded levels and groups; the data left as they were by the command; cut points and labels of binned variables; aweights and iweights. |
| `test_features.do` | `ratio()`, `outcomes()` and `reverse`, each against the plain runs it must reproduce. |
| `test_variants.do` | `individually`, with and without `fixothermargins`, and `insupport`, each against the runs it must reproduce. |
| `test_detail.do` | `detail`: the terms of each block added one at a time, against runs with one block per term. |
| `test_chain.do` | `chain`: each pair against the run on its two values, the totals, joint inference under `bayesboot()` and the prefixes, the tables and the plots. |
| `test_bootstrap_cleanup.do` | Nothing left behind (variables, files, results, seeds) by an interrupted bootstrap or a refused `bayesboot()`, also under `chain`. |
| `test_single_level.do` | Trimming bounds and diagnostics when every raking variable keeps a single level. |
| `test_edge_cases.do` | Edge cases: `focalgroup()` with nearly equal values of `by()`, the count of replicates in the svy warning, labels with SMCL in messages, an unwritable `saving()` file, Break inside the raking, and others. |
| `display_cases.do` | The runs whose output `test_display.do` compares with `reference/transcripts/`. |
| `reference/transcripts/*.txt` | Reference transcripts of `display_cases.do` (Stata 19 and 18 print them identically). |
| `expected_values.do`, `expected_support.do` | Generated scalars `exp_<case>_<quantity>`; do not edit. |
| `data/seldecomp_test.csv` | Generated synthetic data (2,000 observations: `id group y ybin x1 x2 x3 x4 wt cluster strata`). |
| `data/expected_values.json` | The same expected values in JSON. |
| `reference/seldecomp_ref.py` | Python (numpy only) reference implementation that generates the data and the expected values. |

## Tolerances

Means, differences and percentages are compared with `reldif() < 1e-6`; individual raked-weight
extremes (`e(raked_min)`, `e(raked_max)`) and the raked shares of the levels with `1e-5`,
because a single weight or margin is more sensitive to the raking's convergence tolerance
(`tolerance()`, default 1e-8) than a weighted mean is. The tolerance macros sit at
the top of every test file.

## Regenerating the reference transcripts

After a deliberate change to what seldecomp prints, rewrite `reference/transcripts/*.txt` from
the tests directory with the global `SELDECOMP_REGEN_TRANSCRIPTS` set to 1 (`test_display.do`
then writes the logs there instead of comparing), and review the diff:

```stata
global SELDECOMP_REGEN_TRANSCRIPTS 1
do test_display.do
global SELDECOMP_REGEN_TRANSCRIPTS
```

## Regenerating the data and expected values

```sh
python3 reference/seldecomp_ref.py
```

This rewrites `data/seldecomp_test.csv`, `data/expected_values.json`, `expected_values.do` and
`expected_support.do`. A rerun reproduces those files except the last digit of two
`ob_explained` values, which vary with the linear-algebra library. The reference implements
the raking, both modes of `sequentially` (with and without `fixothermargins`), weighted
quantile bins following Stata's `_pctile`, and the two-fold Oaxaca–Blinder explained part used
by `test_oaxaca.do`. Change the seed or the cases in `main()` and rerun; the Stata tests read
whatever the script wrote.

## Manual check: Break

Break cannot be sent from a do-file, so this is checked by hand, in an interactive Stata,
after a change to the raking or the bootstrap. Set up data on which one run takes several
seconds, in a random order, and keep a record of it and of the current results:

```stata
cd /path/to/seldecomp/tests
import delimited using data/seldecomp_test.csv, clear varnames(1) case(preserve) stringcols(8) asdouble
expand 50
gen double u = runiform()
sort u
drop u
gen long ord = _n
datasignature
seldecomp y x1 x2 [pw=wt], by(group)
```

1. Break during the raking: run
   `seldecomp y (a: x1) (b: x2 x4) [pw=wt], by(group) sequentially fixothermargins noisily generate(rk) tolerance(1e-300) iterate(1000000)`
   and press Break while the cycle lines scroll.
2. Break during the bootstrap: run
   `seldecomp y x1 x2 x4 [pw=wt], by(group) bayesboot(reps(500) generate(bw) saving(bsbreak, replace))`
   and press Break while the dots are printed.
3. Break during the BCa jackknife: run
   `seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(5) bca cluster(cluster))`
   and press Break after "Jackknife for BCa acceleration" is printed.

Expected each time: `--Break--` and `r(1);`, and then

- `datasignature` prints the value recorded above, `assert ord == _n` holds and
  `describe, short` shows the data unsorted: no variable was added (no `rk`, `rk_s1` or
  `bw_#`), none changed, the order is the same;
- `seldecomp` replays the results of the last complete run (`e(cmdline)` is
  `seldecomp y x1 x2 [pw=wt], by(group)`);
- no `bsbreak.dta` was written (`dir bsbreak*`), and `mata: mata describe` lists no
  `seldecomp_JK`.

Break almost always lands inside the raking, where the time goes, and is then handled as
above. If it lands between two replicates or deletions (while the weights are drawn), the
command stops at once, with the same results except that the replicate-weight variables
created so far (`bw_1`, ...) and, in the jackknife, the Mata matrix `seldecomp_JK` are left
behind; drop them by hand (`drop bw_*`, `mata: mata drop seldecomp_JK`).

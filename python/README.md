# seldecomp -- Raking-based reweighting decomposition of group differences in means

## Description

`seldecomp` decomposes the difference in means in a variable *y* between two groups, a focal
group and a reference group, using the reweighting method of Pylypchuk and Selden (2008). This
approach relies on reweighting the focal group to match the reference group with respect to the
distribution of one or more binary or categorical raking variables, implemented via raking
(iterative proportional fitting).

Given an initial set of sample weights (or unweighted data) and the reweighted data for the focal
group, `seldecomp` presents three primary objects. First, the full `difference` in means (focal
group mean *y* - reference group mean *y*). Second, the `explained` component of the difference
in means (focal group mean *y* - reweighted focal group mean *y*). Finally, the `unexplained`
component of the difference in means (reweighted focal group mean *y* - reference group mean
*y*). These figures are computed comparing the focal and reference groups after reweighting based
on the full set of raking variables specified and, optionally, on user-specified subsets of the
reweighting variables.

Inference is available using the Bayesian bootstrap. A variety of convenience tools and methods
for visualizing decomposition results are also available.

## Installation

From GitHub:

```
pip install "seldecomp @ git+https://github.com/thegland/seldecomp.git#subdirectory=python"
```

or from a clone of the repository:

```
pip install ./python
```

The package needs Python 3.9 or later, numpy, pandas, and scipy; the plots need matplotlib
(`seldecomp[plots]` in place of `seldecomp` above, or `pip install "./python[plots]"`).

## Syntax

`sd.seldecomp()` (after `import seldecomp as sd`) takes a pandas DataFrame *df* and either the
options as keyword arguments with the names of the Stata options (first form) or the Stata
command line as text, with or without the leading word `seldecomp` (second form). Option names
may be abbreviated as in Stata; variable names may not. It returns the results, *res*
(see [Stored results](#stored-results)).

```
res = sd.seldecomp(df, "y rakespec", by="groupvar", [weights="wt", subset=mask,
                   individually=True, sequentially=True, fixothermargins=True, main_options,
                   multistage_decomposition_options, inference_options, raking_options,
                   reporting_options])

res = sd.seldecomp(df, "y rakespec [if] [in] [weight], by(groupvar) [individually
                   sequentially fixothermargins main_options multistage_decomposition_options
                   inference_options raking_options reporting_options]")
```

where *groupvar* specifies the two groups (lowest value becomes the focal group) and *rakespec*
is a list of raking variables, possibly containing the following elements:

- binary or categorical variables, listed as *x* or `i.`*x* (string variables accepted)
- continuous variables, tagged as `c.`*x*; requires option `contvars` to specify how to convert
  to categorical
- interactions between variables, listed as *x1*`##`*x2* or `c.`*x1*`##i.`*x2*
- sets of variables in parentheses (optionally with a name), such as `(`*x1 x2 x3*`##`*x4*`)` or
  `(`*name*`:` *x1 x2 x3*`)`

By default, `seldecomp` computes differences comparing baseline means to those obtained after
raking to all variables in *rakespec*. `seldecomp` can also produce intermediary decompositions
including just subsets of included variables:

- Specifying `individually` causes `seldecomp` to present decompositions where the focal group
  is raked to match the reference group with respect to each variable, or set of variables in
  parentheses, appearing in *rakespec*, ignoring all other variables
- Specifying `sequentially` causes `seldecomp` to present decompositions raking the focal group
  with respect to the first variable (or variable set) in *rakespec*, then the first two
  variables (or variable sets), then the first three, and so on in sequence until all variables
  are included
- Specifying `fixothermargins` modifies `seldecomp`'s behavior for intermediary decompositions,
  such that the focal group's values for all raking variables not actively being used are held
  constant. This is achieved by specifying the focal group's own initial (within sample) values
  as raking targets for these variables, while using the reference group's values as targets
  for the active variables

Beyond the above, options for implementing decompositions of differences in the ratio of means,
chained decomposition of categorical group variables, multi-outcome decompositions, and automated
handling of continuous raking variables are available under [Main Options](#main). Options
controlling multi-stage decompositions, including production of plots visualizing the
decompositions, are available under
[Multi-stage Decomposition Options](#multi-stage-decomposition). Tools for implementing Bayesian
bootstrap-based inference may be found under [Inference Options](#inference). Options
controlling technical aspects of how raking is implemented may be found under
[Raking Options](#raking). Finally, [Reporting Options](#reporting) controls various aspects of
how output is formatted and whether bonus tables of summary statistics and the like are shown.

## Stata and Python

| Stata | Python |
|---|---|
| `e(b)`, `e(V)`, `e(bystage)`, ... | `res.b`, `res.V`, `res.bystage`, ... |
| `e(explained)`, `e(stagenames)`, ... | `res.e["explained"]`, `res.e["stagenames"]`, ... |
| output tables | `print(res)`, same layout |
| `seldecomp,` *display_options* (replay) | `res.display(diagnostics=True, level=90)` |
| `estat bootstrap, all` | `res.estat_bootstrap(all=True)` |
| `makeplot(`*type*`,` *options*`)` | `makeplot=` at estimation; the Figure is `res.figure` |
| `seldecomp, makeplot(`*type*`,` *options*`)` | `res.makeplot(type, ...)`: a matplotlib Figure |
| `twoway(title() name() ...)` | matplotlib arguments of `res.makeplot()`: `title=`, `ax=`, ... |
| `generate(`*newvar*`)` | `generate="newvar"` adds the column to *df*; also `res.weights` |
| `bayesboot(generate(`*stub*`))` | adds the columns `stub_1`, ... to *df*; also `res.bsweights` |
| `bayesboot(saving(`*file*`))` | `saving(file)`: a CSV file, `file.csv` (a name with an extension is used as typed); also `res.bsdraws` |
| `[if] [in]` | `subset=` a boolean mask (an `if` that cannot be read from the text is refused); `res.sample` |
| `[pw=wt]` | `weights="wt"` |
| variable labels | `labels={"var": "variable label"}` |
| value labels | `value_labels={"var": {0: "No", 1: "Yes"}}` |
| | pandas categorical columns are used in the order of their categories and shown by name |
| `cformat()`, `pformat()`, `sformat()` | the same names: `cformat="%9.3f"` |
| `level()` | none at estimation (as in Stata): `bayesboot` `level(90)`; `res.display(level=90)`; `res.makeplot(..., ci=True, level=90)` |
| `by:`, `statsby:`, `bootstrap:`, `jackknife:`, `svy:` | not provided; inference: `bayesboot` |
| error messages and return codes | `sd.SeldecompError`, with `.message` and `.rc` |
| warnings | `sd.SeldecompWarning`, issued through Python's `warnings`; also shown by `print(res)` |

## Options

### Main

- `by="groupvar"` - specifies the binary variable defining the focal group and the reference
  group
- `focalgroup=value` - specifies the value defining the focal group; default is the lowest value
  (or first in sort order if a string)
- `insupport=True` - drop focal group observations with raking variable levels not represented
  in the reference group. These observations should normally get raked to a weight of 0.
  Specifying this option ensures they are dropped from initial mean calculations as well.
- `chain=True` - allows specifying a categorical *groupvar* in `by`. The lowest value is taken as
  the focal group, and the highest value as the final reference group. Decompositions are then
  shown going from the focal group to the next highest (second lowest overall) value in
  *groupvar*, from the second lowest to the next highest, and so on chaining to the final
  reference group value. A total summing all of the individual chained decompositions (i.e.,
  summing each step from focal group to final reference group) is also given. When inference is
  requested, joint inference is used (calculating all decompositions within a given replicate).
- `outcomes="varlist"` (or a list) - specifies additional *y* variables to decompose
  simultaneously (using the same raking, and with joint inference).
- `ratio="ratiovar"` - decompose the difference in the ratio of means between *y* and *ratiovar*
  (i.e., mean(*y*)/mean(*ratiovar*)).
- `contvars="x1 # x2 # ... xN # [, recut pooled|focalgroup|refgroup]"`, or a dict
  (`{"x1": 4, "x2": 5, "recut": True, "sample": "focalgroup"}`) - for each continuous variable
  *x1*-*xN* in *rakespec*, specifies the number of bins to divide the given continuous variable
  into when converting it into a categorical variable. Bins are formed by dividing the sample
  based on (sample-weighted) quantiles, so *x1* `2` divides *x1* at the sample-weighted median.
  By default (`pooled`), bins are formed based on the distribution of the continuous variable
  within the full estimation sample; options `focalgroup` and `refgroup` allow forming bins based
  on the distribution within just the focal or reference groups respectively. When bootstrap or
  other replication-based inference is used, the same main sample binning scheme is used in each
  replicate. The option `recut` switches to reforming bins based on the quantiles observed within
  each replicate.

### Multi-stage Decomposition

- `individually=True` - compute decompositions raking to each variable, or set of variables,
  specified in *rakespec*
- `sequentially=True` - compute decompositions for variable (or variable set) 1 in *rakespec*,
  then variables 1 and 2, ..., then variables 1, 2, ..., and N.
- `fixothermargins=True` - when computing decompositions using subsets of the variables in
  *rakespec*, hold constant all raking variables not actively being raked to reference group
  targets (i.e., ensure raking does not modify these margins incidentally). Achieved by
  including all of these variables as raking variables, but with targets set to the focal
  group's own initial (within sample) values.
- `detail=True` - with `sequentially`, also compute decompositions adding variables sequentially
  within variable sets (equivalent to not using variable sets, except for the variable sets
  serving a grouping/organizing role).
- `reverse=True` - with `sequentially`, also compute decompositions adding variables in reverse
  order (from right to left).
- `makeplot="[plottype, plotoptions]"`, a dict (`{"type": "plottype", ...}`), or `True` - graph
  the multi-step decomposition results. When `res.makeplot(plottype, plotoptions)` is called,
  graph the multi-step decomposition results from the `seldecomp` run that returned *res*; the
  plot options are then keywords (`ci=True`, `format="%9.2f"`,
  `blocklabels={"demo": "Demographics"}`). Plot types include:
  - `stages` - presents a bar chart of the cumulative total explained and unexplained amounts
    observed at each raking stage.
  - `marginal` - the only option for `individually`. Presents a bar chart of the marginal
    explained quantity observed at each raking stage.
  - `waterfall` - the default option for `sequentially`. Presents a waterfall chart of each
    raking stage's marginal contribution to the explained difference, building to the overall
    explained and total difference.
  - `waterfalllevels` - presents a waterfall chart as above, but with initial values starting at
    the focal group mean rather than zero.

  Plot options include:
  - matplotlib arguments of `res.makeplot()` - for passing through standard matplotlib plot
    options, including `title=`, `ax=`, and so on
  - `swapaxis` - to swap bar charts between horizontal/vertical axis alignments
  - `blocklabels(var1 "text" varsetname1 "text" ...)` - to relabel, for display in the plot,
    independent variables and variable set names used in *rakespec*
  - `ci` - to display confidence intervals on the plot (when inference is used)
  - `normal`, `tinterval`, `percentile`, `bc`, or `bca` - specify the type of confidence interval
    used; see [Inference](#inference) section for more details on these interval types. Default
    is first interval type specified as an inference option or, if nothing is specified,
    `normal` (without clustering) or `tinterval` (with clustering).
  - `outcome(varname)` - to specify the outcome for plot production, when multiple outcomes are
    specified
  - `reverse` - to produce the plot for the reverse-order staging, when `sequentially` `reverse`
    is specified
  - `chainstep(#)` - when `chain` is specified, which comparison step within chain to produce
    plots for; default is from focal group to final reference group
  - `format(%fmt)` - to set the display format for values in bar labels

### Inference

Inference is available via `seldecomp`'s built-in `bayesboot` option.

The option `bayesboot` implements inference using the Bayesian bootstrap. Under `bayesboot`,
bootstrap replicates are formed by reweighting the sample rather than by literal resampling,
after which bootstrap inference proceeds as normal. Reweighting is done in each replicate by
multiplying the initial sample weights (taken to be all 1 if no weights are specified) by random
draws from the Gamma((m-2)/m,1) distribution, where m is either the number of observations within
the sample or the number of clusters within the sample (when clustering). Since all observations
with positive weight in the main sample will have at least some weight in a Bayesian bootstrap
replicate, this approach to inference reduces the likelihood of running into raking feasibility
and support problems while bootstrapping.

`bayesboot` takes the following suboptions, as text
(`bayesboot="reps(200) rseed(123) citable"`), as a dict
(`bayesboot={"reps": 200, "rseed": 123, "citable": True}`), or `bayesboot=True` for the
defaults:

- `reps(#)` - number of replicates; default 100
- `cluster(varname)` - specify a clustering variable; in each replicate, only one weight
  multiplier is drawn per cluster
- `normal`, `tinterval`, `percentile`, `bc`, `bca` - show one or more of these confidence
  intervals:
  - `normal`: normal approximation-based intervals (default without clustering)
  - `tinterval`: t-distribution-based intervals with N-1 degrees of freedom, where N is the
    number of observations or clusters (default when clustering)
  - `percentile`: percentile-based intervals
  - `bc`: bias-corrected percentile intervals
  - `bca`: bias-corrected accelerated percentile intervals
- `citable` - shortcut to show all of normal/t-distribution-based intervals, percentile
  intervals, bias-corrected percentile intervals, and (if explicitly added) bias-corrected
  accelerated percentile intervals
- `ties` - adjust the requested bc/bca confidence intervals for ties
- `level(#)` - specify the confidence level
- `rseed(#)` - randomization seed for bootstrapping
- `dots(#)` - print one dot every # replicates (no dots are printed unless `dots()` is given)
- `nodots` - do not display dots
- `generate(stub [, replace])` - save the bootstrap weights for each replicate *i* to variable
  *stub*`_`*i*
- `saving(filename [, replace])` - save the bootstrap replicate results to the specified file

Note: when raking fails within a given replicate (due to support issues, failure to converge to a
solution, etc.), that replicate is dropped (unless the option `force` is specified) and results
are reported just within the set of replicates that successfully completed.

### Raking

- `iterate=#` - set the maximum number of iterations (default = 10,000). If the cap is hit when
  the focal group's raking margins are too distant from the raking targets, raking is considered
  to have failed and results will not be presented. If the cap is hit because the weights have
  not converged to a final set of weights, but the focal group's raking margins match the raking
  targets, raking is still considered to have succeeded.
- `tolerance=#` - sets the convergence tolerance (default 1e-8). Raking stops when the largest
  relative change in the weights over a full raking cycle is below this tolerance.
- `ctrltolerance=#` - when the focal group value for a given raking margin differs from the
  target value by more than this tolerance (default 1e-6) in relative terms, the raking run is
  considered to have failed.
- `trimhiabs=#`, `trimloabs=#` - set upper and lower bounds on the raked weights. `trimhiabs`
  sets the maximum value raked weights can attain; `trimloabs` sets the minimum value.
- `trimhirel=#`, `trimlorel=#` - set upper and lower bounds on the relative change in weights
  induced by raking. `trimhirel` sets the maximum raking factor (raked weight/initial weight);
  `trimlorel` sets the minimum raking factor.
- `trimsupportoverride=True` - assign raked weights of 0 to observations in the focal group with
  raking variable levels not found in the reference group, even if trimming lower bounds would
  otherwise prevent this.
- `trimfrequency="often"|"sometimes"|"once"` - determines when trimming rules are applied to
  weights. `sometimes`, the default, applies trimming at the end of each raking iteration (after
  all raking margins have been looped through, but before the convergence check). `often` does
  so after each raking margin within an iteration, and can perform better when `sometimes` fails
  to converge. `once` trims the final set of raked weights, without regard for whether the raked
  weights ensure the focal group hits the reference group raking targets.
- `force=True` - presents results even when raking fails to bring focal group raking variable
  means within the specified (or default) tolerance range of the raking targets. Also prevents
  bootstrap replicates where raking fails from being dropped.
- `generate="newvar"` (or `"newvar, replace"`) - save the raked weights to the specified
  variable.

### Reporting

- `covariatedetail=True` - displays a summary statistics table showing observation counts and
  weighted shares for each level of each raking variable in the focal group, raked focal group
  (shares only), and reference group.
- `diagnostics=True` - always show the raking diagnostics which otherwise only appear when a
  raking issue is detected.
- `mincell=#` - defines a thin target level, a raking variable level with fewer than the given
  number of observations (default: 10) in the focal group. Monitored for diagnostic purposes.
- `noisily=True` - display information on each raking cycle (i.e., iteration) as it occurs.
- `varlabels=True` - display variable labels (from `labels=`) in place of variable names, in the
  output and in plots.
- `cformat="%fmt"`, `pformat="%fmt"`, `sformat="%fmt"` - set display formats for the main
  estimates, standard errors, and CI bounds (`cformat`), p-values (`pformat`), and test
  statistics (`sformat`).

## Examples

Setup

```python
import pandas as pd
import seldecomp as sd

# nlsw88.dta: in Stata, sysuse nlsw88 and save nlsw88; it is also at
# https://www.stata-press.com/data/r18/nlsw88.dta
df = pd.read_stata("nlsw88.dta", convert_categoricals=False)
```

Decompose wages between unmarried (focal group, `married` = 0) and married women, raked on two
categorical variables

```python
res = sd.seldecomp(df, "wage i.collgrad i.south", by="married")
print(res)
```

Rake on the joint distribution of two variables and on a continuous variable cut into quartiles

```python
res = sd.seldecomp(df, "wage collgrad##south c.ttl_exp", by="married", contvars="ttl_exp 4")
print(res)
```

Decompose a ratio of means and save the raked weights

```python
df["earnings"] = df["wage"] * df["hours"]
res = sd.seldecomp(df, "earnings i.collgrad i.south", by="married", ratio="hours",
                   generate="rkdwts")
print(res)
```

Multi-stage decomposition adding named sets of variables one at a time, with not-yet-added set
held at the focal group's initial (within sample) values. Plot the marginal contribution of each
block (a waterfall plot by default)

```python
res = sd.seldecomp(df, "wage (demo: race collgrad) (geo: south smsa) c.age", by="married",
                   contvars="age 5", sequentially=True, fixothermargins=True, makeplot=True)
print(res)
```

Redraw the plot from the stored results in *res*, but recentering the waterfall plot around levels
of the dependent variable, and with some relabeling of variable sets

```python
fig = res.makeplot("waterfalllevels", blocklabels={"demo": "Demographics", "geo": "Region"},
                   format="%9.2f")
```

Multi-stage decomposition showing the contribution of each variable set in isolation, allowing
all other raking variable values to vary, with a plot showing these contributions:

```python
res = sd.seldecomp(df, "wage (demo: race collgrad) (geo: south smsa) c.age", by="married",
                   contvars="age 5", individually=True, makeplot=True)
print(res)
```

Bayesian bootstrap inference for a multi-stage decomposition. Show normal, percentile, and
bias-corrected confidence intervals; then the same intervals in the layout of `estat bootstrap`,
and a plot of the stages with their confidence intervals

```python
res = sd.seldecomp(df, "wage (demo: race collgrad) (geo: south smsa)", by="married",
                   sequentially=True, bayesboot="reps(200) rseed(123) citable")
print(res)
res.estat_bootstrap(all=True)
fig = res.makeplot("stages", ci=True)
```

A chained decomposition of three education groups (fewer than 12 years of schooling, 12 to 15, 16
or more), as in Selden, Abdus, and Miller (2019): compare each pair of adjacent groups and their
total, with joint bootstrap inference. The output shows the stage table of each pair; the plot
shows the total

```python
df["edgrp"] = (1 + (df["grade"] >= 12) + (df["grade"] >= 16)).where(df["grade"].notna())
res = sd.seldecomp(df, "wage married (geo: south smsa) (exp: c.ttl_exp)", by="edgrp",
                   chain=True, sequentially=True, contvars="ttl_exp 4",
                   bayesboot="reps(100) rseed(3)")
print(res)
fig = res.makeplot()
```

## Stored results

`seldecomp` stores its results in the object it returns, *res*: scalars and macros in the dict
`res.e`, matrices as attributes. The most used are listed here; see [RESULTS.md](RESULTS.md) for
the full index.

Scalars

- `res.e["N"]` - number of observations
- `res.e["N_focal"]` - number of observations, focal group
- `res.e["N_ref"]` - number of observations, reference group
- `res.e["mean_focal"]` - mean of *depvar (y)*, focal group
- `res.e["mean_raked_focal"]` - mean of *depvar (y)*, focal group, raked weights
- `res.e["mean_ref"]` - mean of *depvar (y)*, reference group
- `res.e["difference"]` - `res.e["mean_focal"]` - `res.e["mean_ref"]`
- `res.e["explained"]` - `res.e["mean_focal"]` - `res.e["mean_raked_focal"]`
- `res.e["unexplained"]` - `res.e["mean_raked_focal"]` - `res.e["mean_ref"]`
- `res.e["pct_explained"]` - 100 * `res.e["explained"]`/`res.e["difference"]`
- `res.e["N_stages"]` - number of stages: 1 without `sequentially` or `individually`; the number
  of blocks with `sequentially`; the number of blocks + 1 with `individually`
- `res.e["N_eff_raked"]` - Kish effective sample size of the raked weights
- `res.e["converged"]` - 1 if the raking of every stage converged and met its targets; 0
  otherwise
- `res.e["targets_met"]` - 1 if every raking margin of every stage was met within
  `ctrltolerance`; 0 otherwise
- `res.e["N_reps"]` - number of completed bootstrap replicates

Macros

- `res.e["cmdline"]` - command as typed (with keyword options, the call written as the Stata
  command line)
- `res.e["depvar"]` - name of dependent variable (*y*); with `outcomes`, all outcomes
- `res.e["by"]` - name of group variable
- `res.e["focal_group"]` - value of group variable, focal group
- `res.e["ref_group"]` - value of group variable, reference group
- `res.e["stagetype"]` - `sequentially` or `individually` (staged results only)

Matrices

- `res.b` - decomposition results (a pandas Series indexed by "equation:name"). Equation
  `overall` has the elements `focal_group`, `raked_focal`, `ref_group`, `difference`,
  `explained`, and `unexplained`; a staged decomposition adds equations by stage
- `res.V` - bootstrap variance-covariance matrix of `res.b`
- `res.bystage` - one row per stage: the raked focal group mean, the explained and unexplained
  parts, and the raking diagnostics

Other

- `res.sample` - marks estimation sample

## References

Pylypchuk, Yuriy, and Thomas M. Selden (2008). A discrete choice decomposition analysis of racial
and ethnic differences in children's health insurance coverage. Journal of Health Economics
27(4): 1109-1128.
DOI: [10.1016/j.jhealeco.2007.12.001](https://doi.org/10.1016/j.jhealeco.2007.12.001)

Selden, Thomas M., Salam Abdus, and G. Edward Miller (2019). Decomposing changes in the growth of
U.S. prescription drug use and expenditures, 1999-2016. Health Services Research 54(4): 752-763.
DOI: [10.1111/1475-6773.13164](https://doi.org/10.1111/1475-6773.13164)

## Authorship

Tom Selden & Thomas A. Hegland (thomas.selden@ahrq.hhs.gov, thomas.hegland@ahrq.hhs.gov)

Source code: https://github.com/thegland/seldecomp

Suggested citation:

> Selden, Thomas and Thomas A. Hegland (2026). "seldecomp: Stata module to compute raking
> (reweighting) decompositions of group differences in means." Statistical Software Components,
> Boston College Department of Economics.

Please also cite Pylypchuk and Selden (2008) for the method.

## Disclaimer

This document reflects only the views of the authors and does not necessarily represent the
views of the Agency for Healthcare Research and Quality, the Department of Health and Human
Services, or the United States government. `seldecomp` does not come with any endorsement by or
warranty from the Agency for Healthcare Research and Quality, the Department of Health and Human
Services, or the United States government.

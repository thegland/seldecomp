# seldecomp -- Full index of stored results (see [README.md](README.md) for the command)

## Full Index of Stored results

`seldecomp` stores the following in the object it returns, *res*: scalars and macros in the dict
`res.e`; matrices as pandas DataFrames with the row and column names of Stata, as attributes
(`res.bystage`) and in the dict `res.matrices`:

### Scalars

- `res.e["N"]` - number of observations
- `res.e["N_focal"]` - number of observations, focal group
- `res.e["N_ref"]` - number of observations, reference group
- `res.e["W_focal"]` - sum of weights, focal group
- `res.e["W_ref"]` - sum of weights, reference group
- `res.e["N_zerowt"]` - number of observations with zero weight; these are not in the sample (0
  without weights)
- `res.e["mean_focal"]` - mean of *depvar (y)*, focal group
- `res.e["mean_raked_focal"]` - mean of *depvar (y)*, focal group, raked weights
- `res.e["mean_ref"]` - mean of *depvar (y)*, reference group
- `res.e["difference"]` - `res.e["mean_focal"]` - `res.e["mean_ref"]`
- `res.e["explained"]` - `res.e["mean_focal"]` - `res.e["mean_raked_focal"]`
- `res.e["unexplained"]` - `res.e["mean_raked_focal"]` - `res.e["mean_ref"]`
- `res.e["pct_explained"]` - 100 * `res.e["explained"]`/`res.e["difference"]`
- `res.e["N_stages"]` - number of stages: 1 without `sequentially` or `individually`; the number
  of blocks with `sequentially`; the number of blocks + 1 with `individually`
- `res.e["N_blocks"]` - number of blocks in *rakespec*
- `res.e["N_terms"]` - number of raking variables (terms)
- `res.e["k_eq"]` - number of equations in `res.b`
- `res.e["N_raked"]` - number of focal group observations in the final raking sample
- `res.e["N_focal_zero"]` - number of focal group observations with zero raked weight: their
  level does not occur in the reference group. With `trimloabs` or `trimlorel` and `force`,
  these observations keep the lower bound and are not counted, unless `trimsupportoverride` is
  specified
- `res.e["N_at_bound"]` - number of distinct weights at one or more trimming bounds
- `res.e["N_eff_base"]` - Kish effective sample size of the base weights, raking sample
- `res.e["N_eff_raked"]` - Kish effective sample size of the raked weights
- `res.e["ess_ratio"]` - `res.e["N_eff_base"]`/`res.e["N_eff_raked"]`
- `res.e["N_focal_mincell"]` - number of focal group target levels below `mincell`
- `res.e["pct_increment"]` - 100 * increment/`res.e["difference"]`, last stage. It equals
  `res.e["pct_explained"]` without `sequentially` or `individually`. With `individually`, it is
  the percent of the remainder
- `res.e["mincell"]` - value of `mincell`
- `res.e["pct_W_focal_zero"]` - percent of focal group weight with zero raked weight
- `res.e["mean_focal_zero"]` - mean of *depvar (y)* for these observations
- `res.e["mean_focal_cs"]` - focal group mean on the raking sample, base weights
- `res.e["explained_support"]` - `res.e["mean_focal"]` - `res.e["mean_focal_cs"]`, the change of
  support; only when `res.e["N_focal_zero"]` > 0
- `res.e["explained_cs"]` - `res.e["mean_focal_cs"]` - `res.e["mean_raked_focal"]`; only when
  `res.e["N_focal_zero"]` > 0
- `res.e["N_notes"]` - number of note lines displayed (`res.e["note1"]`, ...)
- `res.e["converged"]` - 1 if the raking of every stage converged and met its targets; 0
  otherwise. Every stage includes the stages of `reverse` and `detail`. A missed target that
  `force` keeps gives 0. A stop at the `iterate` limit gives 0. Column `converged` of
  `res.bystage`, `res.bystage_rev`, and `res.bystage_detail` gives the value of each stage
- `res.e["targets_met"]` - 1 if every raking margin of every stage was met within
  `ctrltolerance`; 0 otherwise (a miss that `force` keeps). Column `badcontrols` of
  `res.bystage` counts the misses by stage
- `res.e["badcontrols"]` - number of raking variables (margins) with a control total not met
  within `ctrltolerance`
- `res.e["maxreldif"]` - largest relative change of a weight at the last iteration; missing with
  `trimfrequency="once"`
- `res.e["maxctrl"]` - largest relative difference between an achieved and a target control
  total
- `res.e["raked_mean"]` - mean of the raked weights
- `res.e["raked_min"]` - minimum of the raked weights
- `res.e["raked_max"]` - maximum of the raked weights
- `res.e["raked_sd"]` - standard deviation of the raked weights
- `res.e["raked_cv"]` - coefficient of variation of the raked weights
- `res.e["factor_mean"]` - mean of the adjustment factors (raked/base weight)
- `res.e["factor_min"]` - minimum of the adjustment factors
- `res.e["factor_max"]` - maximum of the adjustment factors
- `res.e["factor_sd"]` - standard deviation of the adjustment factors
- `res.e["factor_cv"]` - coefficient of variation of the adjustment factors
- `res.e["level"]` - confidence level
- `res.e["N_reps"]` - number of completed bootstrap replicates
- `res.e["N_misreps"]` - number of replicates dropped because raking failed or missed its
  targets (with `chain`, also because a value of the group variable had no weight)
- `res.e["N_jkfail"]` - number of jackknife replicates that failed while assessing BCa
  acceleration (`bayesboot` `bca` only)
- `res.e["N_iterlimit"]` - number of completed replicates whose raking stopped at the `iterate`
  limit, with every margin met
- `res.e["N_forced"]` - number of replicates that missed raking targets and were kept because
  `force` was specified (`bayesboot` with `force` only)
- `res.e["pct_explained_se"]` - bootstrap standard error of the explained share
- `res.e["seed"]` - random-number seed of the bootstrap
- `res.e["bca"]` - 1 if `bca` was specified, 0 otherwise
- `res.e["ties"]` - 1 if `ties` was specified, 0 otherwise
- `res.e["N_clust"]` - number of clusters (`bayesboot` `cluster()`)
- `res.e["df_bs"]` - degrees of freedom of the t intervals: the number of clusters less 1
  (observations without `cluster()`)
- `res.e["df_r"]` - the same as `res.e["df_bs"]`; only when the coefficient table uses t
  (`tinterval` without `normal` or before it, or neither with `cluster()`)
- `res.e["N_outcomes"]` - number of outcomes (`outcomes` only)
- `res.e["N_pairs"]` - number of pairs (`chain` only)
- `res.e["N_focal_all"]` - number of focal group observations before the exclusion (`insupport`
  only; `res.e["N_focal"]` is the number after it)
- `res.e["W_focal_all"]` - their sum of weights (`insupport` only)
- `res.e["mean_focal_all"]` - their mean of *depvar (y)*, a ratio with `ratio` (`insupport`
  only)

### Macros

- `res.e["cmdline"]` - command as typed (with keyword options, the call written as the Stata
  command line)
- `res.e["title"]` - title in estimation output
- `res.e["depvar"]` - name of dependent variable (*y*); with `outcomes`, all outcomes
- `res.e["outcomes"]` - all outcomes (`outcomes` only)
- `res.e["chain"]` - `chain` or empty
- `res.e["by_levels"]` - (`chain`) values of the group variable, in ascending order
- `res.e["by_labels"]` - (`chain`) these values with their value labels, as displayed, separated
  by `|`
- `res.e["ratio"]` - *ratiovar* of `ratio`
- `res.e["reverse"]` - `reverse` or empty
- `res.e["detail"]` - `detail` or empty
- `res.e["diffsign"]` - `focal_group - ref_group`: the order of the difference
- `res.e["insupport"]` - `insupport` or empty
- `res.e["by"]` - name of group variable
- `res.e["focal_group"]` - value of group variable, focal group
- `res.e["ref_group"]` - value of group variable, reference group
- `res.e["focal_group_label"]` - focal group value with its value label, as displayed
- `res.e["ref_group_label"]` - reference group value with its value label, as displayed
- `res.e["rakespec"]` - *rakespec* as typed
- `res.e["rakevars"]` - names of the raking variables
- `res.e["blocks"]` - names of the blocks
- `res.e["stagenames"]` - names of the stages (`overall` without `sequentially` or
  `individually`)
- `res.e["stage#"]` - raking variables of stage *#*
- `res.e["block#"]` - raking variables of block *#*
- `res.e["term#"]` - name of raking variable *#*
- `res.e["levels#"]` - level labels of raking variable *#*, separated by `|`. A `|` in a label
  shows as `/`. A non-integer value of a `float` variable shows in the fewest digits that
  identify it, as .1
- `res.e["note#"]` - text of note line *#*; a line that lists cells under a note starts with
  four blanks
- `res.e["contvars"]` - contents of `contvars`
- `res.e["stagetype"]` - `sequentially` or `individually` (staged results only)
- `res.e["fixothermargins"]` - `fixothermargins` or empty
- `res.e["wtype"]` - weight type
- `res.e["wexp"]` - weight expression
- `res.e["rakeopts"]` - raking options in effect, with the defaults: `iterate()`,
  `tolerance()`, `ctrltolerance()`, and the trimming options given
- `res.e["generate"]` - name given in `generate`
- `res.e["bsfile"]` - file given in `bayesboot` `saving()`, with `.csv` added to a name
  without an extension
- `res.e["bsgenerate"]` - stub given in `bayesboot` `generate()`
- `res.e["bsnames"]` - names of the bootstrapped elements of `res.b`
- `res.e["citypes"]` - confidence interval types shown, from `normal`, `t`, `percentile`, `bc`,
  and `bca`; the coefficient table shows the first
- `res.e["vce"]` - `bayesboot` or empty
- `res.e["missing"]` - `missing` if some replicates were dropped
- `res.e["worstvar"]` - raking variable with the largest miss of a control total; only when
  `res.e["badcontrols"]` > 0
- `res.e["worstcat"]` - its level
- `res.e["contsample"]` - sample of the quantile cut points
- `res.e["cutvars"]` - binned continuous variables
- `res.e["cutq"]` - number of quantile groups of each binned variable
- `res.e["cutbins"]` - number of nonempty quantile groups of each (fewer than `res.e["cutq"]`
  with ties)
- `res.e["recut"]` - `recut`, if specified
- `res.e["vcetype"]` - title used to label Std. err.
- `res.e["clustvar"]` - name of cluster variable

### Matrices

- `res.b` - decomposition results (a pandas Series indexed by "equation:name";
  `res.matrices["b"]` is the one-row DataFrame). Equation `overall` has the elements
  `focal_group`, `raked_focal`, `ref_group`, `difference`, `explained`, and `unexplained`.
  `sequentially` adds the equations `marg_expl` (the marginal explained part of each block),
  `cum_expl` and `cum_unexpl` (the cumulative parts), and `raked_focal` (the raked focal group
  mean), each by stage. `detail` then adds `detail` (the contribution of each term). `reverse`
  then adds `rev_marg_expl`, `rev_cum_expl`, `rev_cum_unexpl`, and `rev_raked_focal`.
  `individually` adds `marg_expl` (the contribution of each block, and `_remainder`) and
  `raked_focal` (the raked mean of the raking of each block alone). With `outcomes`, each
  outcome *y* has these equations in turn, named *y* (for `overall`), *y*`_marg_expl`, and so on
- `res.V` - bootstrap variance-covariance matrix of `res.b`
- `res.bystage` - one row per stage (named as in `res.e["stagenames"]`), 30 columns:
  `raked_focal`, `cum_expl`, `cum_unexpl`, `pct_explained`, `increment`, then the raking
  diagnostics, the effective sample sizes, `N_focal_mincell`, `pct_increment`, and the
  out-of-support values. These columns have the names of the scalars above (`converged` ...
  `mean_focal_cs`). The scalars hold the values of the last row (`res.e["converged"]`: the
  minimum of column `converged`). With `individually`: one row per block, then a row `overall`
  for the raking of all blocks. Columns 2, 3, 5, and 27 are then `explained`, `unexplained` (the
  parts of the raking of that row), `contribution` (of the block; on the `overall` row, the
  remainder), and `pct_contribution`
- `res.bystage_rev` - (`reverse`) the same for the reverse order, one row per stage (named by
  the block added, last block first)
- `res.raking_checks` - one row per raking stage; columns `iterations_converged`,
  `missed_variables`, `missed_targets`, `targets_checked`, and `control_scale`. The iteration
  flag is measured before `force` marks an unmet-target result as not converged. Missed
  variables count terms with at least one missed level target; missed targets count individual
  level totals. The scale converts original base-weight units to the scale used for
  `ctrltolerance`
- `res.raking_misses` - one row per affected variable and stage, only when targets are missed;
  columns `stage`, `term`, `level`, `target`, `achieved`, and `discrepancy`. Each row identifies
  that variable's worst missed level. Target and achieved totals use original weight units.
  Discrepancy is `abs(s*achieved-s*target)/(abs(s*target)+1)`, where *s* is the stage's
  `control_scale`. Macros `res.e["raking_miss#_var"]` and `res.e["raking_miss#_level"]`
  identify the term and level
- `res.raking_checks_rev` - (`reverse`) checks for the reverse stages; `res.raking_misses_rev`
  gives their missed targets. The label macros are `res.e["raking_miss_rev#_var"]` and
  `res.e["raking_miss_rev#_level"]`. The corresponding `_detail` matrices and macros describe
  term stages under `detail`. Chain pairs 2 and later use numeric suffixes on the matrices and
  labels `res.e["raking_miss`*pair*`_`*row*`_var"]` and
  `res.e["raking_miss`*pair*`_`*row*`_level"]`. Outcomes share the same raking checks and missed
  targets
- `res.bystage_detail` - (`detail`) the same for the stage of each term, one row per term (named
  as in equation `detail`). Column `increment` is the contribution of the term
- `res.bystage#` - (`outcomes`) the same for outcome *#* = 2, 3, ...; `res.bystage` is that of
  the first outcome. The raking diagnostics are the same for all outcomes. `res.bystage_rev#`
  with `reverse`, `res.bystage_detail#` with `detail`
- `res.by_outcome` - (`outcomes`) one row per outcome. Columns `mean_focal`,
  `mean_raked_focal`, `mean_ref`, `difference`, `explained`, `unexplained`, `pct_explained`,
  `pct_explained_se`, `pct_explained_ll`, and `pct_explained_ul` (missing without `bayesboot`);
  with `insupport`, also `mean_focal_all`. The scalars of the same names refer to the first
  outcome
- `res.ratio_parts` - (`ratio`) the numerator and denominator means (columns) of each ratio.
  Rows *y*`:focal_group`, *y*`:raked_focal`, and *y*`:ref_group` for each outcome *y*
- `res.covariates` - (`covariatedetail`) one row per level (*t*`:`*l* = level *l* of
  `res.e["term`*t*`"]`). Columns `focal_group`, `ref_group`, and `raked_focal` (weighted
  shares), then `N_focal` and `N_ref` (unweighted counts)
- `res.cuts#` - cut points of binned continuous variable *#* in `res.e["cutvars"]`
- `res.pairs` - (`chain`) one row per pair (`p1`, ...). Columns `focal_group` and `ref_group`:
  the two values that the pair compares
- `res.by_pair` - (`chain`) one row per pair. Columns `mean_focal`, `mean_raked_focal`,
  `mean_ref`, `difference`, `explained`, `unexplained`, `pct_explained`, `N_focal`, `N_ref`,
  `W_focal`, and `W_ref`; with `insupport`, also `N_focal_all` and `mean_focal_all`; with
  `bayesboot`, also `pct_explained_se`, `pct_explained_ll`, and `pct_explained_ul`
- `res.cuts#_p#` - (`chain`) cut points of binned continuous variable *#*, pair *#*
- `res.reps` - number of completed replicates per bootstrapped element
- `res.pct_explained_ci` - percentile confidence interval of the explained share
- `res.b_bs` - means of the bootstrap replicates
- `res.bias` - bootstrap bias estimates
- `res.se` - bootstrap standard errors
- `res.ci_normal` - normal-based confidence intervals
- `res.ci_t` - t-based confidence intervals (`res.e["df_bs"]` degrees of freedom)
- `res.ci_percentile` - percentile confidence intervals
- `res.ci_bc` - bias-corrected confidence intervals
- `res.ci_bca` - bias-corrected and accelerated confidence intervals (`bca`)
- `res.z0` - bias-correction constants
- `res.accel` - acceleration constants (`bca`)
- `res.bsdraws` - bootstrap replicates: one row per replicate (missing values for a failed one),
  one column per bootstrapped element of `res.b`, named as in `res.e["bsnames"]`
- `res.bsdraws_pct` - replicates of the explained share, one column (`pct_explained`); not
  stored when the difference is 0. With `outcomes`, one column per outcome *y*
  (*y*`:pct_explained`), always stored

### Other

- `res.sample` - marks estimation sample (a boolean Series)
- `res.weights` - the weights of each stage, one column per stage: the raked weights of the
  focal group and the base weights of the reference group (the columns that `generate` adds to
  the DataFrame)
- `res.bsweights` - (`bayesboot` `generate()`) the bootstrap weights of each replicate, one
  column per replicate, one row per observation of the estimation sample (the columns that
  `bayesboot` `generate()` adds)
- `res.figure` - the matplotlib Figure drawn by the `makeplot` option (`None` without it)
- `res.trace` - the lines that `noisily` printed
- `res.notes` - the note lines (`res.e["note1"]`, ...)
- `res.messages` - the lines that `print(res)` shows before the tables, except those of
  `bayesboot`

With `ratio`, `res.e["mean_focal"]`, `res.e["mean_raked_focal"]`, `res.e["mean_ref"]`,
`res.e["mean_focal_zero"]`, `res.e["mean_focal_cs"]`, and the stage values are ratios of means.
With `outcomes`, the scalars (`res.e["mean_focal"]` ... `res.e["pct_explained"]`,
`res.e["pct_explained_se"]`, `res.e["explained_support"]`, `res.e["explained_cs"]`),
`res.pct_explained_ci`, and `res.bystage` refer to the first outcome. The results of the other
outcomes are in `res.by_outcome`, `res.bystage#`, and `res.b`.

With `chain`, `res.b` holds the equations for each pair *j* (identified by numbering the
decompositions pairings from lowest value to next highest value), with the prefix `p`*j*
(`p`*j*, `p`*j*`_marg_expl`, ...). The related items without the prefix refer to the total
across the chained decompositions. `res.e["mean_focal"]` and `res.e["mean_ref"]` are the means
(ratios with `ratio`) of the first and the last values of the chained group variable. With
`insupport`, `res.e["mean_focal"]` is the focal group mean of pair 1 after the exclusion; the
mean of the first value is `mean_focal_all` in `res.by_pair`. `res.e["difference"]` is the sum
of the differences of the pairs. It equals `res.e["mean_focal"]` - `res.e["mean_ref"]`, unless
`insupport` excludes observations of a value other than the first. `res.e["focal_group"]`,
`res.e["ref_group"]`, and their labels are the first and the last values. `res.bystage` is the
stage matrix of pair 1, and `res.bystage#` that of pair *#*. `res.e["converged"]` and
`res.e["targets_met"]` are the minima over the pairs. `res.ratio_parts` has the rows
`p`*j*`:focal_group`, `p`*j*`:raked_focal`, and `p`*j*`:ref_group`.

`res.V` is created by `seldecomp` only with `bayesboot` (`None` otherwise).

`res.estat_bootstrap()` returns the intervals of the results of `bayesboot` as a dict, at the
confidence level of the last `res.display()`, else `res.e["level"]` (`res.display()`, the
replay, returns nothing). These calls do not change `res.e`. The items returned include:

Scalars

- `"level"` - confidence level of the intervals below

Matrices

- `"ci_normal"` - normal-based confidence intervals
- `"ci_t"` - t-based confidence intervals
- `"ci_percentile"` - percentile confidence intervals
- `"ci_bc"` - bias-corrected confidence intervals
- `"ci_bca"` - bias-corrected and accelerated confidence intervals (`bca`)
- `"pct_explained_ci"` - percentile confidence interval of the explained share; with
  `outcomes`, one row per outcome; with `chain`, one row per pair and one for the total

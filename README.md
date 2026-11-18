# seldecomp: Raking-Based Reweighting Decomposition of Group Differences in Means

_seldecomp_ is a Stata program that decomposes the difference in means in a variable y between two groups, a _focal group_ and a _reference group_, using the reweighting method of Pylypchuk and Selden (2008) - an alternative to Oaxaca-Blinder-Kitagawa decomposition. This approach relies on reweighting the focal group to match the reference group with respect to the distribution of one or more binary or categorical raking variables, implemented via nonparametric raking (iterative proportional fitting).

Given an initial set of sample weights (or unweighted data) and the reweighted data for the focal group, seldecomp presents three primary objects. First, the _full difference_ in means (focal group mean y - reference group mean y). Second, the _explained_ component of the difference in means (focal group mean y - reweighted focal group mean y). Finally, the _unexplained_ component of the difference in means (reweighted focal group mean y - reference group mean y). These figures are computed comparing the focal and reference groups after reweighting based on the full set of raking variables specified and, optionally, on user-specified subsets of the reweighting variables.

Other features include:
- Inference via the Bayesian bootstrap or via Stata's built in bootstrap:, jackknife:, and svy: features.
- Multi-stage decompositions, showing how results vary as each raking variable or set of raking variables is added (either individually, or in sequence).
- Production of plots visualizing multi-stage decomposition results.
- Decomposition of ratios of means.
- Chained decompositions of categorical group variables.
- Multi-outcome decompositions.
- Automated handling of continuous raking variables. 

Note that in addition to the main Stata program, a (less extensively tested) implementation in Python is also available.  Its options follow the same patterns of use and naming scheme as the Stata version of the program; it is detailed in the Python subfolder of this project. 

<br>

### How to Install
```
net install seldecomp, from(https://raw.githubusercontent.com/thegland/seldecomp/main/) replace
```
or
```
ssc install seldecomp
```
Note: _seldecomp_ requires Stata version 16.1 or later.

<br>

The Python version of the program may be installed as follows:
```
pip install "seldecomp @ git+https://github.com/thegland/seldecomp.git#subdirectory=python"
```
It requires Python 3.9+, numpy, pandas, scipy, and matplotlib. 

<br>

### Basic Syntax

_For more detail on syntax, please refer to seldecomp's Stata help file._

```
seldecomp y rakespec [if] [in] [aw/pw/iw=weights], by(groupvar) [individually sequentially fixothermargins other_options]
```

- `y` is the outcome variable
- `groupvar` specifies the two groups (lowest value becomes the focal group)
- `rakespec` is a list of raking variables, possibly containing the following elements:
  - binary or categorical variables, listed as `x` or `i.x` (string variables accepted)
  - continuous variables, tagged as `c.x`; requires option `contvars()` to specify how to convert to categorical
  - interactions between variables, listed as `x1##x2` or `c.x1##i.x2`
  - sets of variables in parentheses (optionally with a name), such as `(x1 x2 x3##x4)` or `(name: x1 x2 x3)`

By default, _seldecomp_ computes differences comparing baseline means to those obtained after raking to all variables in `rakespec`. _seldecomp_ can also produce intermediary decompositions including just subsets of included variables:

- Specifying `individually` causes _seldecomp_ to present decompositions where the focal group is raked to match the reference group with respect to each variable, or set of variables in parentheses, appearing in `rakespec`, ignoring all other variables
- Specifying `sequentially` causes _seldecomp_ to present decompositions raking the focal group with respect to the first variable (or variable set) in `rakespec`, then the first two variables (or variable sets), then the first three, and so on in sequence until all variables are included
- Specifying `fixothermargins` modifies _seldecomp's_ behavior for intermediary decompositions, such that the focal group's values for all raking variables not actively being used are held constant (achieved by specifying the focal group's own initial values as raking targets for these variables, while using the reference group's values as targets for the active variables)

Note: _seldecomp_ is compatible with estimation prefixes like `by:`, `statsby:`, `bootstrap:`, `jackknife:`, and `svy:`. It is also compatible with standard post-estimation commands, including `estat bootstrap`, `esttab`, and the use of `seldecomp,` or `seldecomp, makeplot()` for estimate replay.

For more detail on all of _seldecomp's_ features, including for how to implement decompositions of differences in the ratio of means, chained decomposition of categorical group variables, multi-outcome decompositions, automated handling of continuous raking variables, production of plots visualizing multi-stage decompositions, Bayesian bootstrap-based inference, and more - please refer to the Stata help file.

<br>

### Author Information and Contact

*Authors:* Tom Selden & Thomas A. Hegland (thomas.selden@ahrq.hhs.gov, thomas.hegland@ahrq.hhs.gov)

*Citation:* Selden, Thomas and Thomas A. Hegland (2026). "seldecomp: Stata module to compute raking (reweighting) decompositions of group differences in means." Statistical Software Components, Boston College Department of Economics.

Please also cite Pylypchuk and Selden (2008) for the method.

*Disclaimer:* This document reflects only the views of the authors and does not necessarily represent the views of the Agency for Healthcare Research and Quality, the Department of Health and Human Services, or the United States government. _seldecomp_ does not come with any endorsement by or warranty from the Agency for Healthcare Research and Quality, the Department of Health and Human Services, or the United States government.

<br>

### Citations

Pylypchuk, Yuriy, and Thomas M. Selden (2008). A discrete choice decomposition analysis of racial and ethnic differences in children's health insurance coverage. Journal of Health Economics 27(4): 1109-1128. https://doi.org/10.1016/j.jhealeco.2007.12.001

Selden, Thomas M., Salam Abdus, and G. Edward Miller (2019). Decomposing changes in the growth of U.S. prescription drug use and expenditures, 1999-2016. Health Services Research 54(4): 752-763. https://doi.org/10.1111/1475-6773.13164

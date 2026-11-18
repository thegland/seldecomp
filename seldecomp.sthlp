{smcl}
{* *! version 1.0.0  18nov2026}{...}
{vieweralsosee "oaxaca" "help oaxaca"}{...}
{vieweralsosee "fairlie" "help fairlie"}{...}
{vieweralsosee "[R] bootstrap" "help bootstrap"}{...}
{vieweralsosee "[R] bayesboot" "help bayesboot"}{...}
{viewerjumpto "Description" "seldecomp##description"}{...}
{viewerjumpto "Syntax" "seldecomp##syntax"}{...}
{viewerjumpto "Options" "seldecomp##options"}{...}
{viewerjumpto "Examples" "seldecomp##examples"}{...}
{viewerjumpto "Stored results" "seldecomp##results"}{...}
{viewerjumpto "References" "seldecomp##references"}{...}
{viewerjumpto "Author" "seldecomp##author"}{...}
{viewerjumpto "Disclaimer" "seldecomp##disclaimer"}{...}
{cmd:help seldecomp}{right:{browse "https://github.com/thegland/seldecomp"}}
{hline}

{title:Title}

{phang}{bf:seldecomp} {hline 2} Raking-based reweighting decomposition of group differences in means


{marker description}{...}
{title:Description}

{pstd}
{cmd:seldecomp} decomposes the difference in means in a variable 
{it:y} between two groups, a focal group and a reference group, 
using the reweighting method of Pylypchuk and Selden (2008). This 
approach relies on reweighting the focal group to match the reference 
group with respect to the distribution of one or more binary or categorical 
raking variables, implemented via raking (iterative proportional fitting).

{pstd}
Given an initial set of sample weights (or unweighted data) and 
the reweighted data for the focal group, {cmd:seldecomp} presents 
three primary objects. First, the full {cmd:difference} in means 
(focal group mean {it:y} - reference group mean {it:y}). Second, 
the {cmd:explained} component of the difference in means (focal 
group mean {it:y} - reweighted focal group mean {it:y}). Finally, 
the {cmd:unexplained} component of the difference in means (reweighted 
focal group mean {it:y} - reference group mean {it:y}). These 
figures are computed comparing the focal and reference groups after 
reweighting based on the full set of raking variables specified 
and, optionally, on user-specified subsets of the reweighting variables.

{pstd}
Inference is available using a variety of methods, including the 
Bayesian bootstrap, bootstrap, and balanced repeated replication (BRR). 
A variety of convenience tools and methods for visualizing 
decomposition results are also available.


{marker syntax}{...}
{title:Syntax}

{p 4 13 2}
{cmd:seldecomp}
{it:y}
{it:rakespec}
{ifin}
{weight}{cmd:,}
{opt by(groupvar)}
[{opt ind:ividually}
{opt seq:uentially}
{opt fix:othermargins}
{help seldecomp##main:{it:main_options}}
{help seldecomp##multistage:{it:multistage_decomposition_options}}
{help seldecomp##inference:{it:inference_options}}
{help seldecomp##raking:{it:raking_options}}
{help seldecomp##reporting:{it:reporting_options}}]

{pstd}
where {it:groupvar} specifies the two groups (lowest value becomes 
the focal group) and {it:rakespec} is a list of raking variables, 
possibly containing the following elements:{p_end}
{p 8 10 2}- binary or categorical variables, listed as {it:x} or {cmd:i.}{it:x} (string variables accepted){p_end}
{p 8 10 2}- continuous variables, tagged as {cmd:c.}{it:x}; requires 
option {opt contvars()} to specify how to convert to categorical{p_end}
{p 8 10 2}- interactions between variables, listed as {it:x1}{cmd:##}{it:x2} or {cmd:c.}{it:x1}{cmd:##i.}{it:x2}{p_end}
{p 8 10 2}- sets of variables in parentheses (optionally with a name), such as {cmd:(}{it:x1 x2 x3}{cmd:##}{it:x4}{cmd:)} or 
{cmd:(}{it:name}{cmd::} {it:x1 x2 x3}{cmd:)}{p_end}

{pstd}
By default, {cmd:seldecomp} computes differences comparing baseline means to those obtained after raking to all variables 
in {it:rakespec}. {cmd:seldecomp} can also produce intermediary 
decompositions including just subsets of included variables:{p_end}
{p 8 10 2}- Specifying {opt individually} causes {cmd:seldecomp} to present decompositions where the focal group is raked 
to match the reference group with respect to each variable, or set of variables in parentheses, appearing in {it:rakespec}, 
ignoring all other variables{p_end}
{p 8 10 2}- Specifying {opt sequentially} causes {cmd:seldecomp} 
to present decompositions raking the focal group with respect
 to the first variable (or variable set) in {it:rakespec}, then 
the first two variables (or variable sets), then the first three, 
 and so on in sequence until all variables are included{p_end}
{p 8 10 2}- Specifying {opt fixothermargins} modifies {cmd:seldecomp}'s 
behavior for intermediary decompositions, such that the 
focal group's values for all raking variables not actively being used are held constant. This is achieved 
	by specifying the focal group's own initial (within sample) values as raking targets for these variables, 
	while using the reference group's values as targets for the active variables{p_end}

{pstd}
Beyond the above, options for implementing decompositions of 
differences in the ratio of means, chained decomposition of categorical 
group variables, multi-outcome decompositions, and automated handling of continuous raking variables are available under 
{help seldecomp##main:Main Options}. Options controlling multi-stage 
decompositions, including production of plots visualizing 
the decompositions, are available under 
{help seldecomp##multistage:Multi-stage Decomposition Options}. Tools for implementing 
Bayesian bootstrap-based inference, as well as other types of inference, may be found under 
{help seldecomp##inference:Inference Options}. Options controlling technical aspects of how raking is implemented may be 
found under {help seldecomp##raking:Raking Options}. Finally, {help seldecomp##reporting:Reporting Options} controls various 
aspects of how output is formatted and whether bonus tables of summary statistics and the like are shown.

{pstd}
Note: {cmd:seldecomp} is compatible with estimation prefixes like {cmd:by:}, {cmd:statsby:}, {cmd:bootstrap:}, 
{cmd:jackknife:}, and {cmd:svy:} (replicates type inference 
only). It is also compatible with standard post-estimation commands, 
including {cmd:estat bootstrap}, {cmd:esttab}, and the use of 
'{cmd:seldecomp},'  or '{cmd:seldecomp, makeplot()'} for estimate replay.


{marker options}{...}
{title:Options}

{marker main}{...}
{dlgtab:Main}

{phang}
{opt by(groupvar)} - specifies the binary variable defining the focal group and the reference group

{phang}
{opt focal:group(value)} - specifies the value defining the focal group; default is the lowest value (or first in 
	sort order if a string)

{phang}
{opt insup:port} - drop focal group observations with raking variable levels not represented in the reference group. These 
observations should normally get raked to a weight of 0. Specifying this option ensures they are dropped from initial mean
 calculations as well.

{phang}
{opt chain} - allows specifying a categorical {it:groupvar} in {opt by()}. The lowest value is taken as the focal 
group, and the highest value as the final reference group. Decompositions 
are then shown going from the focal group to the next highest (second lowest overall) 
value in {it:groupvar}, from the second lowest to the next highest, and so 
on chaining to the final reference group value. A total summing all of the 
individual chained decompositions (i.e., summing each step from focal group 
to final reference group) is also given. When inference is requested, joint 
inference is used (calculating all decompositions within a given replicate).

{phang}
{opt out:comes(varlist)} - specifies additional {it:y} variables to decompose simultaneously (using the same 
	raking, and with joint inference).

{phang}
{opt rat:io(ratiovar)} - decompose the difference in the ratio of means between {it:y} and {it:ratiovar} (i.e., 
	mean({it:y})/mean({it:ratiovar})).

{phang}
{cmdab:cont:vars}{cmd:(}{it:x1 # x2 # ... xN #} [{cmd:,} {cmd:recut} {cmd:pooled}|{cmd:focalgroup}|{cmd:refgroup}]{cmd:)} - 
for each continuous variable {it:x1}-{it:xN} in {it:rakespec}, specifies the number of bins to divide the given continuous 
variable into when converting it into a categorical variable. Bins are formed by dividing the sample based on 
(sample-weighted) quantiles, so {it:x1} {cmd:2} divides {it:x1} at the sample-weighted median. By default ({cmd:pooled}), 
bins are formed based on the distribution of the continuous 
variable within the full estimation sample; options {cmd:focalgroup} 
and {cmd:refgroup} allow forming bins based on the distribution within just the focal or reference groups respectively. When 
bootstrap or other replication-based inference is used, the same main sample binning scheme is used in each replicate. The 
option {cmd:recut} switches to reforming bins based on the quantiles observed within each replicate.


{marker multistage}{...}
{dlgtab:Multi-stage Decomposition}

{phang}
{opt ind:ividually} - compute decompositions raking to each variable, or set of variables, specified in {it:rakespec}

{phang}
{opt seq:uentially} - compute decompositions for variable (or variable set) 1 in {it:rakespec}, then variables 1 and 2, 
..., then variables 1, 2, ..., and N.

{phang}
{opt fix:othermargins} - when computing decompositions using subsets of the variables in {it:rakespec}, hold constant 
all raking variables not actively being raked to reference group targets (i.e., ensure raking does not modify these 
	margins incidentally). Achieved by including all of these variables as raking variables, but with targets set 
to the focal group's own initial (within sample) values.

{phang}
{opt det:ail} - with {opt sequentially}, also compute decompositions adding variables sequentially within variable sets 
(equivalent to not using variable sets, except for the variable sets serving a grouping/organizing role).

{phang}
{opt rev:erse} - with {opt sequentially}, also compute decompositions adding variables in reverse order (from right to left).

{phang}
{cmdab:make:plot}{cmd:([}{it:plottype} {cmd:,} {it:plotoptions}]{cmd:)} - graph the multi-step decomposition results. 
When {cmd:seldecomp, makeplot(}{it:spec}{cmd:)} is specified, graph the multi-step decomposition results from the prior 
{cmd:seldecomp} run. Plot types include:{p_end}
{p 8 10 2}- {opt stages} - presents a bar chart of the cumulative total explained and unexplained amounts observed at each 
raking stage.{p_end}
{p 8 10 2}- {opt marginal} - the only option for {opt individually}. Presents a bar chart of the marginal explained quantity 
observed at each raking stage.{p_end}
{p 8 10 2}- {opt waterfall} - the default option for {opt sequentially}. Presents a waterfall chart of each raking stage's 
marginal contribution to the explained difference, building to the overall explained and total difference.{p_end}
{p 8 10 2}- {opt waterfalllevels} - presents a waterfall chart as above, but with initial values starting at the focal 
group mean rather than zero.{p_end}

{pmore}
Plot options include:{p_end}
{p 12 14 2}- {cmdab:tw:oway}{cmd:(}{it:{help twoway_options}}{cmd:)} - for passing through standard twoway plot options, 
including {opt title()}, {opt name()}, and so on{p_end}
{p 12 14 2}- {opt swap:axis} - to swap bar charts between horizontal/vertical axis alignments{p_end}
{p 12 14 2}- {cmdab:block:labels}{cmd:(}{it:var1} {cmd:"}{it:text}{cmd:"} 
{it:varsetname1} {cmd:"}{it:text}{cmd:"} ...{cmd:)} 
- to relabel, for display in the plot, independent variables and variable set names used in {it:rakespec}{p_end}
{p 12 14 2}- {opt ci} - to display confidence intervals on the plot (when inference is used){p_end}
{p 12 14 2}- {opt normal}, {opt tinterval}, {opt percentile}, {opt bc}, or {opt bca} - specify the 
type of confidence interval used; see {help seldecomp##inference:Inference} section for more details 
on these interval types. Default is first interval type specified as an inference option or, if 
nothing is specified, {opt normal} (without clustering) or {opt tinterval} (with clustering).{p_end}
{p 12 14 2}- {opt out:come(varname)} - to specify the outcome 
for plot production, when multiple outcomes are specified{p_end}
{p 12 14 2}- {opt rev:erse} - to produce the plot for the reverse-order 
staging, when {opt sequentially} {opt reverse} is specified{p_end}
{p 12 14 2}- {opt chain:step(#)} - when {opt chain} is specified, which comparison step within chain to produce plots for; 
default is from focal group to final reference group{p_end}
{p 12 14 2}- {opt for:mat(%fmt)} - to set the display format for values in bar labels{p_end}


{marker inference}{...}
{dlgtab:Inference}

{pstd}
Inference is available via {cmd:seldecomp}'s built-in {opt bayesboot()} option, and also via the Stata default 
bootstrap, jackknife, and survey-based inference features. The Stata default methods are implemented using the 
prefixes below:{p_end}
{p 8 10 2}- {helpb bootstrap:bootstrap:} {cmd:seldecomp ...}{p_end}
{p 8 10 2}- {helpb jackknife:jackknife:} {cmd:seldecomp ...}{p_end}
{p 8 10 2}- {helpb svy_brr:svy brr:} {cmd:seldecomp ...}{p_end}
{p 8 10 2}- {helpb svy_bootstrap:svy bootstrap:} {cmd:seldecomp ...}{p_end}
{p 8 10 2}- {helpb svy_jackknife:svy jackknife:} {cmd:seldecomp ...}{p_end}
{p 8 10 2}- {helpb svy_sdr:svy sdr:} {cmd:seldecomp ...}{p_end}

{pstd}
The option {opt bayesboot()} implements inference using the Bayesian bootstrap. Under {opt bayesboot()}, bootstrap 
replicates are formed by reweighting the sample rather than by literal resampling, after which bootstrap inference 
proceeds as normal. Reweighting is done in each replicate by multiplying the initial sample weights (taken to be all 
	1 if no weights are specified) by random draws from the Gamma((m-2)/m,1) distribution, where m is either the number 
of observations within the sample or the number of clusters within the sample (when clustering). Since all 
observations with positive weight in the main sample will have at least some 
weight in a Bayesian bootstrap replicate, this approach to inference reduces the likelihood of running into raking 
feasibility and support problems while bootstrapping.

{pstd}
{opt bayesboot()} takes the following suboptions:{p_end}
{p 8 10 2}- {opt r:eps(#)} - number of replicates; default 100{p_end}
{p 8 10 2}- {opt cl:uster(varname)} - specify a clustering variable; in each replicate, only one weight multiplier 
is drawn per cluster{p_end}
{p 8 10 2}- {opt norm:al}, {opt tint:erval}, {opt perc:entile}, {opt bc}, 
{opt bca} - show one or more of these confidence intervals:{p_end}
{p 12 14 2}- {opt normal}: normal approximation-based intervals (default without clustering){p_end}
{p 12 14 2}- {opt tinterval}: t-distribution-based intervals with N-1 degrees of freedom, where N is the number 
of observations or clusters (default when clustering){p_end}
{p 12 14 2}- {opt percentile}: percentile-based intervals{p_end}
{p 12 14 2}- {opt bc}: bias-corrected percentile intervals{p_end}
{p 12 14 2}- {opt bca}: bias-corrected accelerated percentile intervals{p_end}
{p 8 10 2}- {opt cit:able} - shortcut to show all of normal/t-distribution-based intervals, percentile intervals, 
bias-corrected percentile intervals, and (if explicitly added) bias-corrected accelerated percentile intervals{p_end}
{p 8 10 2}- {opt ties} - adjust the requested bc/bca confidence intervals for ties{p_end}
{p 8 10 2}- {opt l:evel(#)} - specify the confidence level    {p_end}
{p 8 10 2}- {opt rseed(#)} - randomization seed for bootstrapping{p_end}
{p 8 10 2}- {opt dots(#)} - print one dot every # replicates; default 1{p_end}
{p 8 10 2}- {opt nodots} - do not display dots{p_end}
{p 8 10 2}- {cmdab:gen:erate}{cmd:(}{it:stub} [{cmd:, replace}]{cmd:)} - save the bootstrap weights for each replicate 
{it:i} to variable {it:stub}{cmd:_}{it:i}{p_end}
{p 8 10 2}- {cmdab:sav:ing}{cmd:(}{it:filename} [{cmd:, replace}]{cmd:)} - save the bootstrap replicate results to the 
specified file{p_end}

{pstd}
Note: across all inference types, when raking fails within a given replicate (due to support issues, failure to converge 
	to a solution, etc.), that replicate is dropped (unless the option{cmd: force} is specified) and results are reported 
just within the set of replicates that successfully completed.


{marker raking}{...}
{dlgtab:Raking}

{phang}
{opt iter:ate(#)} - set the maximum number of iterations  (default = 10,000). If the cap is hit when the focal group's 
raking margins are too distant from the raking targets, raking is considered to have failed and results will not be 
presented. If the cap is hit because the weights have not converged to a final set of weights, but the focal group's 
raking margins match the raking targets, raking is still considered to have succeeded.  

{phang}
{opt tol:erance(#)} - sets the convergence tolerance (default 1e-8). Raking stops when the largest relative change in 
the weights over a full raking cycle is below this tolerance.

{phang}
{opt ctrltol:erance(#)} - when the focal group value for a given raking margin differs from the target value by more 
than this tolerance (default 1e-6) in relative terms, the raking run is considered to have failed.

{phang}
{opt trimhiabs(#)}, {opt trimloabs(#)} - set upper and lower bounds on the raked weights. {opt trimhiabs()} sets the 
maximum value raked weights can attain; {opt trimloabs()} sets the minimum value. 

{phang}
{opt trimhirel(#)}, {opt trimlorel(#)} - set upper and lower bounds on the relative change in weights induced by 
raking. {opt trimhirel()} sets the maximum raking factor (raked 
weight/initial weight); {opt trimlorel()} sets the minimum raking factor.

{phang}
{opt trimsupportoverride} - assign raked weights of 0 to observations in the focal group with raking variable levels 
not found in the reference group, even if trimming lower bounds would otherwise prevent this.    

{phang}
{cmdab:trimfreq:uency}{cmd:(often}|{cmd:sometimes}|{cmd:once)} - determines when trimming rules are applied to weights. 
{cmd:sometimes}, the default, applies trimming at the end of each raking iteration (after all raking margins have been 
	looped through, but before the convergence check). {cmd:often} does so after each raking margin within an iteration, 
and can perform better when {cmd:sometimes} fails to converge. {cmd:once} trims the final set of raked weights, without 
regard for whether the raked weights ensure the focal group hits the reference group raking targets.

{phang}
{opt force} - presents results even when raking fails to bring focal group raking variable means within the specified 
(or default) tolerance range of the raking targets. Also prevents bootstrap replicates
 where raking fails from being dropped.      

{phang}
{cmdab:gen:erate}{cmd:(}{it:newvar} [{cmd:, replace}]{cmd:)} - save the raked weights to the specified variable.


{marker reporting}{...}
{dlgtab:Reporting}

{phang}
{opt cov:ariatedetail} - displays a summary statistics table showing observation counts and weighted shares for each 
level of each raking variable in the focal group, raked focal group (shares only), and reference group.

{phang}
{opt diag:nostics} - always show the raking diagnostics which otherwise only appear when a raking issue is detected.

{phang}
{opt min:cell(#)} - defines a thin target level, a raking variable 
level with fewer than the given number of observations (default: 10) in the 
focal group. Monitored for diagnostic purposes.

{phang}
{opt noi:sily} - display information on each raking cycle (i.e., iteration) as it occurs.

{phang}
{opt varl:abels} - display variable labels in place of variable names, in the output and in plots.

{phang}
{opt cformat(%fmt)}, {opt pformat(%fmt)}, {opt sformat(%fmt)} - set display formats for the main estimates, standard errors, 
and CI bounds ({opt cformat()}), p-values ({opt pformat()}), and test statistics ({opt sformat()}).


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{stata "sysuse nlsw88, clear"}{p_end}

{pstd}Decompose wages between unmarried (focal group, {cmd:married} = 0) and married
women, raked on two categorical variables{p_end}
{phang2}{stata "seldecomp wage i.collgrad i.south, by(married)"}{p_end}

{pstd}Rake on the joint distribution of two variables and on a continuous
variable cut into quartiles{p_end}
{phang2}{stata "seldecomp wage collgrad##south c.ttl_exp, by(married) contvars(ttl_exp 4)"}{p_end}

{pstd}Decompose a ratio of means and save the raked weights {p_end}
{phang2}{stata "generate earnings = wage * hours"}{p_end}
{phang2}{stata "seldecomp earnings i.collgrad i.south, by(married) ratio(hours) generate(rkdwts)"}{p_end}

{pstd}Multi-stage decomposition adding named sets of variables one at a time, with not-yet-added set held at the
focal group's initial (within sample) values. Plot the marginal contribution of each block (a
waterfall plot by default){p_end}
{phang2}{stata "seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) sequentially fixothermargins makeplot"}{p_end}

{pstd}Redraw the plot from the stored results in e(), but recentering the waterfall plot around levels of the dependent 
variable, and with some relabeling of variable sets {p_end}
{phang2}{stata "seldecomp, makeplot(waterfalllevels, blocklabels(demo Demographics geo Region) format(%9.2f))"}{p_end}

{pstd}Multi-stage decomposition showing the contribution of each variable set in isolation, allowing 
all other raking variable values to vary, with a plot showing  these contributions:{p_end}
{phang2}{stata "seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) individually makeplot"}{p_end}

{pstd}Bayesian bootstrap inference for a multi-stage decomposition. Show normal, percentile, and
bias-corrected confidence intervals; then the same intervals in the layout of
{cmd:estat bootstrap}, and a plot of the stages with their confidence
intervals{p_end}
{phang2}{stata "seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially bayesboot(reps(200) rseed(123) citable)"}{p_end}
{phang2}{stata "estat bootstrap, all"}{p_end}
{phang2}{stata "seldecomp, makeplot(stages, ci)"}{p_end}

{pstd}Traditional bootstrap-based inference{p_end}
{phang2}{stata "bootstrap: seldecomp wage race collgrad c.ttl_exp, by(married) contvars(ttl_exp 4)"}{p_end}

{pstd}Survey design-based inference with BRR replicate weights (the variable names are
placeholders){p_end}
{phang2}{stata "svyset psu [pw=wt], strata(stratum) brrweight(brr*) fay(.5) vce(brr)"}{p_end}
{phang2}{stata "svy brr: seldecomp y (demo: race educ) (inc: c.income), by(group) contvars(income 4) sequentially fixothermargins"}{p_end}

{pstd}A chained decomposition of three education groups (fewer than 12 years of schooling, 12 to
15, 16 or more), as in Selden, Abdus, and Miller (2019): compare each pair of adjacent
groups and their total, with joint bootstrap inference. The output shows the stage
table of each pair; the plot shows the total{p_end}
{phang2}{stata "generate edgrp = 1 + (grade >= 12) + (grade >= 16) if grade < ."}{p_end}
{phang2}{stata "seldecomp wage married (geo: south smsa) (exp: c.ttl_exp), by(edgrp) chain sequentially contvars(ttl_exp 4) bayesboot(reps(100) rseed(3))"}{p_end}
{phang2}{stata "seldecomp, makeplot"}{p_end}


{marker references}{...}
{title:References}

{phang}
Pylypchuk, Yuriy, and Thomas M. Selden (2008). A discrete choice decomposition
analysis of racial and ethnic differences in children's health insurance
coverage. Journal of Health Economics 27(4): 1109-1128.
DOI: {browse "https://doi.org/10.1016/j.jhealeco.2007.12.001":10.1016/j.jhealeco.2007.12.001}
{p_end}

{phang}
Selden, Thomas M., Salam Abdus, and G. Edward Miller (2019). Decomposing changes
in the growth of U.S. prescription drug use and expenditures, 1999-2016. Health
Services Research 54(4): 752-763.
DOI: {browse "https://doi.org/10.1111/1475-6773.13164":10.1111/1475-6773.13164}
{p_end}


{marker author}{...}
{title:Authorship}

{pstd}
Tom Selden & Thomas A. Hegland (thomas.selden@ahrq.hhs.gov, thomas.hegland@ahrq.hhs.gov)

{pstd}
Source code: {browse "https://github.com/thegland/seldecomp"}

{pstd}
Suggested citation:{p_end}
{pmore}
Selden, Thomas and Thomas A. Hegland (2026). "seldecomp: Stata module to compute raking (reweighting)
decompositions of group differences in means." Statistical Software Components, Boston College Department of Economics.

{pstd}
Please also cite Pylypchuk and Selden (2008) for the method.


{marker disclaimer}{...}
{title:Disclaimer}

{pstd}This document reflects only the views of the authors and does not necessarily represent the views of the Agency
for Healthcare Research and Quality, the Department of Health and Human Services, or the United States government. The
program with which this help file is associated is the work of the authors and does not come with any endorsement by or
warranty from the Agency for Healthcare Research and Quality, the Department of Health and Human Services, or the United
States government.  {p_end}



{marker results}{...}
{title:Stored results}

{pstd}
{cmd:seldecomp} stores its results in {cmd:e()}. The most used are listed here; see 
{helpb seldecomp_results} for the full index.

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_focal)}}number of observations, focal group{p_end}
{synopt:{cmd:e(N_ref)}}number of observations, reference group{p_end}
{synopt:{cmd:e(mean_focal)}}mean of {it:depvar (y)}, focal group{p_end}
{synopt:{cmd:e(mean_raked_focal)}}mean of {it:depvar (y)}, focal group, raked
    weights{p_end}
{synopt:{cmd:e(mean_ref)}}mean of {it:depvar (y)}, reference group{p_end}
{synopt:{cmd:e(difference)}}{cmd:e(mean_focal)} - {cmd:e(mean_ref)}{p_end}
{synopt:{cmd:e(explained)}}{cmd:e(mean_focal)} - {cmd:e(mean_raked_focal)}{p_end}
{synopt:{cmd:e(unexplained)}}{cmd:e(mean_raked_focal)} - {cmd:e(mean_ref)}{p_end}
{synopt:{cmd:e(pct_explained)}}100 * {cmd:e(explained)}/{cmd:e(difference)}{p_end}
{synopt:{cmd:e(N_stages)}}number of stages: 1 without {cmd:sequentially} or
    {cmd:individually}; the number of blocks with {cmd:sequentially}; the
    number of blocks + 1 with {cmd:individually}{p_end}
{synopt:{cmd:e(N_eff_raked)}}Kish effective sample size of the raked
    weights{p_end}
{synopt:{cmd:e(converged)}}1 if the raking of every stage converged and met
    its targets; 0 otherwise{p_end}
{synopt:{cmd:e(targets_met)}}1 if every raking margin of every stage was met
    within {cmd:ctrltolerance()}; 0 otherwise{p_end}
{synopt:{cmd:e(N_reps)}}number of completed bootstrap replicates{p_end}

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable ({it:y}); with {cmd:outcomes()},
    all outcomes{p_end}
{synopt:{cmd:e(by)}}name of group variable{p_end}
{synopt:{cmd:e(focal_group)}}value of group variable, focal group{p_end}
{synopt:{cmd:e(ref_group)}}value of group variable, reference group{p_end}
{synopt:{cmd:e(stagetype)}}{cmd:sequentially} or {cmd:individually} (staged
    results only){p_end}

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}decomposition results. Equation {cmd:overall} has the
    elements {cmd:focal_group}, {cmd:raked_focal}, {cmd:ref_group},
    {cmd:difference}, {cmd:explained}, and {cmd:unexplained}; a staged
    decomposition adds equations by stage{p_end}
{synopt:{cmd:e(V)}}bootstrap variance-covariance matrix of {cmd:e(b)}{p_end}
{synopt:{cmd:e(bystage)}}one row per stage: the raked focal group mean, the
    explained and unexplained parts, and the raking diagnostics{p_end}

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

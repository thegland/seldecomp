{smcl}
{* *! version 1.0.0  18nov2026}{...}
{vieweralsosee "seldecomp" "help seldecomp"}{...}
{cmd:help seldecomp_results}{right:{browse "https://github.com/thegland/seldecomp"}}
{hline}

{title:Title}

{phang}{bf:seldecomp} {hline 2} Full index of stored results (see {helpb seldecomp} for the command)


{marker results}{...}
{title:Full Index of Stored results}

{pstd}
{cmd:seldecomp} stores the following in {cmd:e()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_focal)}}number of observations, focal group{p_end}
{synopt:{cmd:e(N_ref)}}number of observations, reference group{p_end}
{synopt:{cmd:e(W_focal)}}sum of weights, focal group{p_end}
{synopt:{cmd:e(W_ref)}}sum of weights, reference group{p_end}
{synopt:{cmd:e(N_zerowt)}}number of observations with zero weight; these are
    not in the sample (0 without weights and under {cmd:svy jackknife}){p_end}
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
{synopt:{cmd:e(N_blocks)}}number of blocks in {it:rakespec}{p_end}
{synopt:{cmd:e(N_terms)}}number of raking variables (terms){p_end}
{synopt:{cmd:e(k_eq)}}number of equations in {cmd:e(b)}{p_end}
{synopt:{cmd:e(N_raked)}}number of focal group observations in the final
    raking sample{p_end}
{synopt:{cmd:e(N_focal_zero)}}number of focal group observations with zero
    raked weight: their level does not occur in the reference group. With
    {cmd:trimloabs()} or {cmd:trimlorel()} and  {cmd:force}, these observations keep the lower
    bound and are not counted, unless {cmd:trimsupportoverride} is
    specified{p_end}
{synopt:{cmd:e(N_at_bound)}}number of distinct weights at one or more trimming bounds{p_end}
{synopt:{cmd:e(N_eff_base)}}Kish effective sample size of the base weights,
    raking sample{p_end}
{synopt:{cmd:e(N_eff_raked)}}Kish effective sample size of the raked
    weights{p_end}
{synopt:{cmd:e(ess_ratio)}}{cmd:e(N_eff_base)}/{cmd:e(N_eff_raked)}{p_end}
{synopt:{cmd:e(N_focal_mincell)}}number of focal group target levels below
    {cmd:mincell()}{p_end}
{synopt:{cmd:e(pct_increment)}}100 * increment/{cmd:e(difference)}, last stage.
    It equals {cmd:e(pct_explained)} without {cmd:sequentially} or
    {cmd:individually}. With {cmd:individually}, it is the percent of the
    remainder{p_end}
{synopt:{cmd:e(mincell)}}value of {cmd:mincell()}{p_end}
{synopt:{cmd:e(pct_W_focal_zero)}}percent of focal group weight with zero raked
    weight{p_end}
{synopt:{cmd:e(mean_focal_zero)}}mean of {it:depvar (y)} for these
    observations{p_end}
{synopt:{cmd:e(mean_focal_cs)}}focal group mean on the raking sample, base
    weights{p_end}
{synopt:{cmd:e(explained_support)}}{cmd:e(mean_focal)} - {cmd:e(mean_focal_cs)},
    the change of support; only when {cmd:e(N_focal_zero)} > 0{p_end}
{synopt:{cmd:e(explained_cs)}}{cmd:e(mean_focal_cs)} - {cmd:e(mean_raked_focal)};
    only when {cmd:e(N_focal_zero)} > 0{p_end}
{synopt:{cmd:e(N_notes)}}number of note lines displayed ({cmd:e(note1)},
    ...){p_end}
{synopt:{cmd:e(converged)}}1 if the raking of every stage converged and met
    its targets; 0 otherwise. Every stage includes the stages of
    {cmd:reverse} and {cmd:detail}. A missed target that {cmd:force} keeps
    gives 0. A stop at the {cmd:iterate()} limit gives 0. Column
    {cmd:converged} of {cmd:e(bystage)}, {cmd:e(bystage_rev)}, and
    {cmd:e(bystage_detail)} gives the value of each stage{p_end}
{synopt:{cmd:e(targets_met)}}1 if every raking margin of every stage was met
    within {cmd:ctrltolerance()}; 0 otherwise (a miss that {cmd:force}
    keeps). Column {cmd:badcontrols} of {cmd:e(bystage)}
    counts the misses by stage{p_end}
{synopt:{cmd:e(badcontrols)}}number of raking variables (margins) with a
    control total not met within {cmd:ctrltolerance()}{p_end}
{synopt:{cmd:e(maxreldif)}}largest relative change of a weight at the last
    iteration; missing with {cmd:trimfrequency(once)}{p_end}
{synopt:{cmd:e(maxctrl)}}largest relative difference between an achieved and
    a target control total{p_end}
{synopt:{cmd:e(raked_mean)}}mean of the raked weights{p_end}
{synopt:{cmd:e(raked_min)}}minimum of the raked weights{p_end}
{synopt:{cmd:e(raked_max)}}maximum of the raked weights{p_end}
{synopt:{cmd:e(raked_sd)}}standard deviation of the raked weights{p_end}
{synopt:{cmd:e(raked_cv)}}coefficient of variation of the raked weights{p_end}
{synopt:{cmd:e(factor_mean)}}mean of the adjustment factors (raked/base
    weight){p_end}
{synopt:{cmd:e(factor_min)}}minimum of the adjustment factors{p_end}
{synopt:{cmd:e(factor_max)}}maximum of the adjustment factors{p_end}
{synopt:{cmd:e(factor_sd)}}standard deviation of the adjustment factors{p_end}
{synopt:{cmd:e(factor_cv)}}coefficient of variation of the adjustment
    factors{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}
{synopt:{cmd:e(level_replay)}}(hidden; {cmd:bayesboot()}) level of the last
    replay, after a replay at a level other than {cmd:e(level)}.
    {cmd:estat bootstrap} uses this level{p_end}
{synopt:{cmd:e(N_reps)}}number of completed bootstrap replicates{p_end}
{synopt:{cmd:e(N_misreps)}}number of replicates dropped because raking failed
    or missed its targets (with {cmd:chain}, also because a value of the group
    variable had no weight){p_end}
{synopt:{cmd:e(N_jkfail)}}number of jackknife replicates that failed while assessing
    BCa acceleration ({cmd:bayesboot(bca)} only){p_end}
{synopt:{cmd:e(N_prefail)}}number of replicates that will fail for lack of
    focal group support, found before they run (with {cmd:chain}, also for a
    value of the group variable without weight). Only under {cmd:svy} with
    {cmd:svyset} replicate weights (but not {cmd:svy brr} with
    {cmd:hadamard()}) and under {cmd:svy jackknife} by PSU{p_end}
{synopt:{cmd:e(N_iterlimit)}}number of completed replicates whose raking
    stopped at the {cmd:iterate()} limit, with every margin met{p_end}
{synopt:{cmd:e(N_forced)}}number of replicates that missed raking targets and
    were kept because {cmd:force} was specified ({cmd:bayesboot()} with
    {cmd:force} only){p_end}
{synopt:{cmd:e(pct_explained_se)}}bootstrap standard error of the explained
    share{p_end}
{synopt:{cmd:e(seed)}}random-number seed of the bootstrap{p_end}
{synopt:{cmd:e(bca)}}1 if {cmd:bca} was specified, 0 otherwise{p_end}
{synopt:{cmd:e(ties)}}1 if {cmd:ties} was specified, 0 otherwise{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters ({cmd:bayesboot(cluster())}){p_end}
{synopt:{cmd:e(df_bs)}}degrees of freedom of the t intervals: the number of
    clusters less 1 (observations without {cmd:cluster()}){p_end}
{synopt:{cmd:e(df_r)}}the same as {cmd:e(df_bs)}; only when the coefficient
    table uses t ({cmd:tinterval} without {cmd:normal} or before it, or neither
    with {cmd:cluster()}). {cmd:test} and {cmd:lincom} use it{p_end}
{synopt:{cmd:e(N_outcomes)}}number of outcomes ({cmd:outcomes()} only){p_end}
{synopt:{cmd:e(N_pairs)}}number of pairs ({cmd:chain} only){p_end}
{synopt:{cmd:e(N_focal_all)}}number of focal group observations before the
    exclusion ({cmd:insupport} only; {cmd:e(N_focal)} is the number after
    it){p_end}
{synopt:{cmd:e(W_focal_all)}}their sum of weights ({cmd:insupport}
    only){p_end}
{synopt:{cmd:e(mean_focal_all)}}their mean of {it:depvar (y)}, a ratio with
    {cmd:ratio()} ({cmd:insupport} only){p_end}

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:seldecomp}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(estat_cmd)}}{cmd:seldecomp_estat}, the program for
    {cmd:estat}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable ({it:y}); with {cmd:outcomes()},
    all outcomes{p_end}
{synopt:{cmd:e(outcomes)}}all outcomes ({cmd:outcomes()} only){p_end}
{synopt:{cmd:e(chain)}}{cmd:chain} or empty{p_end}
{synopt:{cmd:e(by_levels)}}({cmd:chain}) values of the group variable, in
    ascending order{p_end}
{synopt:{cmd:e(by_labels)}}({cmd:chain}) these values with their value labels,
    as displayed, separated by {cmd:|}{p_end}
{synopt:{cmd:e(ratio)}}{it:ratiovar} of {cmd:ratio()}{p_end}
{synopt:{cmd:e(reverse)}}{cmd:reverse} or empty{p_end}
{synopt:{cmd:e(detail)}}{cmd:detail} or empty{p_end}
{synopt:{cmd:e(diffsign)}}{cmd:focal_group - ref_group}: the order of the
    difference{p_end}
{synopt:{cmd:e(insupport)}}{cmd:insupport} or empty{p_end}
{synopt:{cmd:e(by)}}name of group variable{p_end}
{synopt:{cmd:e(focal_group)}}value of group variable, focal group{p_end}
{synopt:{cmd:e(ref_group)}}value of group variable, reference group{p_end}
{synopt:{cmd:e(focal_group_label)}}focal group value with its value label, as
    displayed{p_end}
{synopt:{cmd:e(ref_group_label)}}reference group value with its value label,
    as displayed{p_end}
{synopt:{cmd:e(rakespec)}}{it:rakespec} as typed{p_end}
{synopt:{cmd:e(rakevars)}}names of the raking variables{p_end}
{synopt:{cmd:e(blocks)}}names of the blocks{p_end}
{synopt:{cmd:e(stagenames)}}names of the stages ({cmd:overall} without
    {cmd:sequentially} or {cmd:individually}){p_end}
{synopt:{cmd:e(stage}{it:#}{cmd:)}}raking variables of stage {it:#}{p_end}
{synopt:{cmd:e(block}{it:#}{cmd:)}}raking variables of block {it:#}{p_end}
{synopt:{cmd:e(term}{it:#}{cmd:)}}name of raking variable {it:#}{p_end}
{synopt:{cmd:e(levels}{it:#}{cmd:)}}level labels of raking variable {it:#},
    separated by {cmd:|}. A {cmd:|} in a label shows as {cmd:/}. A non-integer
    value of a {cmd:float} variable shows in the fewest digits that identify
    it, as .1{p_end}
{synopt:{cmd:e(note}{it:#}{cmd:)}}text of note line {it:#}; a line that lists
    cells under a note starts with four blanks{p_end}
{synopt:{cmd:e(contvars)}}contents of {cmd:contvars()}{p_end}
{synopt:{cmd:e(stagetype)}}{cmd:sequentially} or {cmd:individually} (staged
    results only){p_end}
{synopt:{cmd:e(fixothermargins)}}{cmd:fixothermargins} or empty{p_end}
{synopt:{cmd:e(wtype)}}weight type{p_end}
{synopt:{cmd:e(wexp)}}weight expression{p_end}
{synopt:{cmd:e(rakeopts)}}raking options in effect, with the defaults:
    {cmd:iterate()}, {cmd:tolerance()}, {cmd:ctrltolerance()}, and the
    trimming options given{p_end}
{synopt:{cmd:e(generate)}}name given in {cmd:generate()}{p_end}
{synopt:{cmd:e(bsfile)}}file given in {cmd:bayesboot(saving())}{p_end}
{synopt:{cmd:e(bsgenerate)}}stub given in {cmd:bayesboot(generate())}{p_end}
{synopt:{cmd:e(bsnames)}}names of the bootstrapped elements of
    {cmd:e(b)}{p_end}
{synopt:{cmd:e(citypes)}}confidence interval types shown, from {cmd:normal},
    {cmd:t}, {cmd:percentile}, {cmd:bc}, and {cmd:bca}; the coefficient table
    shows the first{p_end}
{synopt:{cmd:e(vce)}}{cmd:bayesboot} or empty{p_end}
{synopt:{cmd:e(missing)}}{cmd:missing} if some replicates were dropped{p_end}
{synopt:{cmd:e(worstvar)}}raking variable with the largest miss of a control
    total; only when {cmd:e(badcontrols)} > 0{p_end}
{synopt:{cmd:e(worstcat)}}its level{p_end}
{synopt:{cmd:e(contsample)}}sample of the quantile cut points{p_end}
{synopt:{cmd:e(cutvars)}}binned continuous variables{p_end}
{synopt:{cmd:e(cutq)}}number of quantile groups of each binned
    variable{p_end}
{synopt:{cmd:e(cutbins)}}number of nonempty quantile groups of each (fewer than
    {cmd:e(cutq)} with ties){p_end}
{synopt:{cmd:e(recut)}}{cmd:recut}, if specified{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. err.{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(rngstate)}}random-number state at the start of the
    bootstrap{p_end}

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}decomposition results. Equation {cmd:overall} has the
    elements {cmd:focal_group}, {cmd:raked_focal}, {cmd:ref_group},
    {cmd:difference}, {cmd:explained}, and {cmd:unexplained}.
    {cmd:sequentially} adds the equations {cmd:marg_expl} (the marginal
    explained part of each block), {cmd:cum_expl} and {cmd:cum_unexpl} (the
    cumulative parts), and {cmd:raked_focal} (the raked focal group mean), each
    by stage. {cmd:detail} then adds {cmd:detail} (the contribution of each
    term). {cmd:reverse} then adds {cmd:rev_marg_expl}, {cmd:rev_cum_expl},
    {cmd:rev_cum_unexpl}, and {cmd:rev_raked_focal}. {cmd:individually} adds
    {cmd:marg_expl} (the contribution of each block, and {cmd:_remainder}) and
    {cmd:raked_focal} (the raked mean of the raking of each block alone). With
    {cmd:outcomes()}, each outcome {it:y} has these equations in turn, named
    {it:y} (for {cmd:overall}), {it:y}{cmd:_marg_expl}, and so on{p_end}
{synopt:{cmd:e(V)}}bootstrap variance-covariance matrix of {cmd:e(b)}{p_end}
{synopt:{cmd:e(bystage)}}one row per stage (named as in {cmd:e(stagenames)}),
    30 columns: {cmd:raked_focal}, {cmd:cum_expl}, {cmd:cum_unexpl},
    {cmd:pct_explained}, {cmd:increment}, then the raking diagnostics, the
    effective sample sizes, {cmd:N_focal_mincell}, {cmd:pct_increment}, and the
    out-of-support values. These columns have the names of the scalars above
    ({cmd:converged} ... {cmd:mean_focal_cs}). The scalars hold the values of
    the last row ({cmd:e(converged)}: the minimum of column {cmd:converged}).
    With {cmd:individually}: one row per block, then a row {cmd:overall} for
    the raking of all blocks. Columns 2, 3, 5, and 27 are then
    {cmd:explained}, {cmd:unexplained} (the parts of the raking of that row),
    {cmd:contribution} (of the block; on the {cmd:overall} row, the
    remainder), and {cmd:pct_contribution}{p_end}
{synopt:{cmd:e(bystage_rev)}}({cmd:reverse}) the same for the reverse order,
    one row per stage (named by the block added, last block first){p_end}
{synopt:{cmd:e(raking_checks)}}one row per raking stage; columns
    {cmd:iterations_converged}, {cmd:missed_variables}, {cmd:missed_targets},
    {cmd:targets_checked}, and {cmd:control_scale}. The iteration flag is
    measured before {cmd:force} marks an unmet-target result as not converged.
    Missed variables count terms with at least one missed level target;
    missed targets count individual level totals. The scale converts original
    base-weight units to the scale used for {cmd:ctrltolerance()}{p_end}
{synopt:{cmd:e(raking_misses)}}one row per affected variable and stage, only
    when targets are missed; columns {cmd:stage}, {cmd:term}, {cmd:level},
    {cmd:target}, {cmd:achieved}, and {cmd:discrepancy}. Each row identifies
    that variable's worst missed level. Target and achieved totals use original
    weight units. Discrepancy is {cmd:abs(s*achieved-s*target)/(abs(s*target)+1)},
    where {it:s} is the stage's {cmd:control_scale}. Macros
    {cmd:e(raking_miss}{it:#}{cmd:_var)} and
    {cmd:e(raking_miss}{it:#}{cmd:_level)} identify the term and level{p_end}
{synopt:{cmd:e(raking_checks_rev)}}({cmd:reverse}) checks for the reverse
    stages; {cmd:e(raking_misses_rev)} gives their missed targets. The label
    macros are {cmd:e(raking_miss_rev}{it:#}{cmd:_var)} and
    {cmd:e(raking_miss_rev}{it:#}{cmd:_level)}. The corresponding
    {cmd:_detail} matrices and macros describe term stages under {cmd:detail}.
    Chain pairs 2 and later use numeric suffixes on the matrices and labels
    {cmd:e(raking_miss}{it:pair}{cmd:_}{it:row}{cmd:_var)} and
    {cmd:e(raking_miss}{it:pair}{cmd:_}{it:row}{cmd:_level)}. Outcomes share
    the same raking checks and missed targets{p_end}
{synopt:{cmd:e(bystage_detail)}}({cmd:detail}) the same for the stage of each
    term, one row per term (named as in equation {cmd:detail}). Column
    {cmd:increment} is the contribution of the term{p_end}
{synopt:{cmd:e(bystage}{it:#}{cmd:)}}({cmd:outcomes()}) the same for outcome
    {it:#} = 2, 3, ...; {cmd:e(bystage)} is that of the first outcome. The
    raking diagnostics are the same for all outcomes.
    {cmd:e(bystage_rev}{it:#}{cmd:)} with {cmd:reverse},
    {cmd:e(bystage_detail}{it:#}{cmd:)} with {cmd:detail}{p_end}
{synopt:{cmd:e(by_outcome)}}({cmd:outcomes()}) one row per outcome. Columns
    {cmd:mean_focal}, {cmd:mean_raked_focal}, {cmd:mean_ref},
    {cmd:difference}, {cmd:explained}, {cmd:unexplained},
    {cmd:pct_explained}, {cmd:pct_explained_se}, {cmd:pct_explained_ll}, and
    {cmd:pct_explained_ul} (missing without {cmd:bayesboot()}); with
    {cmd:insupport}, also {cmd:mean_focal_all}. The scalars of the same names
    refer to the first outcome{p_end}
{synopt:{cmd:e(ratio_parts)}}({cmd:ratio()}) the numerator and denominator
    means (columns) of each ratio. Rows {it:y}{cmd::focal_group},
    {it:y}{cmd::raked_focal}, and {it:y}{cmd::ref_group} for each outcome
    {it:y}{p_end}
{synopt:{cmd:e(covariates)}}({cmd:covariatedetail}) one row per level
    ({it:t}{cmd::}{it:l} = level {it:l} of {cmd:e(term}{it:t}{cmd:)}).
    Columns {cmd:focal_group}, {cmd:ref_group}, and {cmd:raked_focal}
    (weighted shares), then {cmd:N_focal} and {cmd:N_ref} (unweighted
    counts){p_end}
{synopt:{cmd:e(cuts}{it:#}{cmd:)}}cut points of binned continuous variable
    {it:#} in {cmd:e(cutvars)}{p_end}
{synopt:{cmd:e(pairs)}}({cmd:chain}) one row per pair ({cmd:p1}, ...).
    Columns {cmd:focal_group} and {cmd:ref_group}: the two values that the
    pair compares{p_end}
{synopt:{cmd:e(by_pair)}}({cmd:chain}) one row per pair. Columns
    {cmd:mean_focal}, {cmd:mean_raked_focal}, {cmd:mean_ref},
    {cmd:difference}, {cmd:explained}, {cmd:unexplained},
    {cmd:pct_explained}, {cmd:N_focal}, {cmd:N_ref}, {cmd:W_focal}, and
    {cmd:W_ref}; with {cmd:insupport}, also {cmd:N_focal_all} and
    {cmd:mean_focal_all}; with {cmd:bayesboot()}, also
    {cmd:pct_explained_se}, {cmd:pct_explained_ll}, and
    {cmd:pct_explained_ul}{p_end}
{synopt:{cmd:e(cuts}{it:#}{cmd:_p}{it:#}{cmd:)}}({cmd:chain}) cut points of
    binned continuous variable {it:#}, pair {it:#}{p_end}
{synopt:{cmd:e(reps)}}number of completed replicates per bootstrapped
    element{p_end}
{synopt:{cmd:e(pct_explained_ci)}}percentile confidence interval of the
    explained share{p_end}
{synopt:{cmd:e(b_bs)}}means of the bootstrap replicates{p_end}
{synopt:{cmd:e(bias)}}bootstrap bias estimates{p_end}
{synopt:{cmd:e(se)}}bootstrap standard errors{p_end}
{synopt:{cmd:e(ci_normal)}}normal-based confidence intervals{p_end}
{synopt:{cmd:e(ci_t)}}t-based confidence intervals ({cmd:e(df_bs)} degrees of
    freedom){p_end}
{synopt:{cmd:e(ci_percentile)}}percentile confidence intervals{p_end}
{synopt:{cmd:e(ci_bc)}}bias-corrected confidence intervals{p_end}
{synopt:{cmd:e(ci_bca)}}bias-corrected and accelerated confidence intervals
    ({cmd:bca}){p_end}
{synopt:{cmd:e(z0)}}bias-correction constants{p_end}
{synopt:{cmd:e(accel)}}acceleration constants ({cmd:bca}){p_end}
{synopt:{cmd:e(bsdraws)}}bootstrap replicates: one row per replicate (missing
    values for a failed one), one column per bootstrapped element of
    {cmd:e(b)}, named as in {cmd:e(bsnames)}. Not stored with more than
    {cmd:c(max_matdim)} replicates{p_end}
{synopt:{cmd:e(bsdraws_pct)}}replicates of the explained share, one column
    ({cmd:pct_explained}); not stored when the difference is 0. With
    {cmd:outcomes()}, one column per outcome {it:y}
    ({it:y}{cmd::pct_explained}), always stored. Like {cmd:e(bsdraws)}, not
    stored with more than {cmd:c(max_matdim)} replicates{p_end}

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{pstd}
With {cmd:ratio()}, {cmd:e(mean_focal)}, {cmd:e(mean_raked_focal)},
{cmd:e(mean_ref)}, {cmd:e(mean_focal_zero)}, {cmd:e(mean_focal_cs)}, and the
stage values are ratios of means. With {cmd:outcomes()}, the scalars
({cmd:e(mean_focal)} ... {cmd:e(pct_explained)}, {cmd:e(pct_explained_se)},
{cmd:e(explained_support)}, {cmd:e(explained_cs)}), {cmd:e(pct_explained_ci)},
and {cmd:e(bystage)} refer to the first outcome. The results of the other
outcomes are in {cmd:e(by_outcome)}, {cmd:e(bystage}{it:#}{cmd:)}, and
{cmd:e(b)}.

{pstd}
With {cmd:chain}, {cmd:e(b)} holds the equations for each pair {it:j} (identified by numbering the 
	decompositions pairings from lowest value to next highest value), with the
prefix {cmd:p}{it:j} ({cmd:p}{it:j}, {cmd:p}{it:j}{cmd:_marg_expl}, ...). The related items without 
the prefix refer to the total across the chained decompositions. {cmd:e(mean_focal)} and
{cmd:e(mean_ref)} are the means (ratios with {cmd:ratio()}) of the first and
the last values of the chained group variable. With {cmd:insupport},
{cmd:e(mean_focal)} is the focal group mean of pair 1 after the exclusion; the
mean of the first value is {cmd:mean_focal_all} in {cmd:e(by_pair)}.
{cmd:e(difference)} is the sum of the differences of the pairs. It equals
{cmd:e(mean_focal)} - {cmd:e(mean_ref)}, unless {cmd:insupport} excludes
observations of a value other than the first. {cmd:e(focal_group)},
{cmd:e(ref_group)}, and their labels are the first and the last values.
{cmd:e(bystage)} is the stage matrix of pair 1, and {cmd:e(bystage}{it:#}{cmd:)}
that of pair {it:#}. {cmd:e(converged)} and {cmd:e(targets_met)} are the
minima over the pairs. {cmd:e(ratio_parts)} has the rows
{cmd:p}{it:j}{cmd::focal_group}, {cmd:p}{it:j}{cmd::raked_focal}, and
{cmd:p}{it:j}{cmd::ref_group}.  

{pstd}
{cmd:e(V)} is created by {cmd:seldecomp} only with {cmd:bayesboot()}, though use of a prefix like 
{cmd:bootstrap:}, {cmd:jackknife:}, etc. will create an {cmd:e(V)} matrix of their own.


{pstd}
A replay of the results of {cmd:bayesboot()} ({cmd:seldecomp} typed without
arguments, with or without {cmd:level()}) leaves the intervals that it
displays in {cmd:r()}, at the requested confidence level. {cmd:estat bootstrap} does
the same. These commands do not change
{cmd:e()}, except {cmd:e(level_replay)} The replay items stored in {cmd:r()} include:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:r(level)}}confidence level of the intervals below{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:r(ci_normal)}}normal-based confidence intervals{p_end}
{synopt:{cmd:r(ci_t)}}t-based confidence intervals (replay only){p_end}
{synopt:{cmd:r(ci_percentile)}}percentile confidence intervals{p_end}
{synopt:{cmd:r(ci_bc)}}bias-corrected confidence intervals{p_end}
{synopt:{cmd:r(ci_bca)}}bias-corrected and accelerated confidence intervals
    ({cmd:bca}){p_end}
{synopt:{cmd:r(pct_explained_ci)}}percentile confidence interval of the
    explained share; with {cmd:outcomes()}, one row per outcome; with
    {cmd:chain}, one row per pair and one for the total{p_end}
{p2colreset}{...}

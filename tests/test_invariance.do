* test_invariance.do -- paths not asserted elsewhere, and invariance properties:
*   1. e(cuts1), e(cuts2) with two binned continuous variables (and a term joining two)
*   2. the level labels of c.x#i.z and i.z#c.x
*   3. aweights and iweights with generate() and with bayesboot()
*   4. invariance: scaling the weights, an affine transform of the outcome, duplicated
*      observations (iweights as frequency weights), the order of the rows, relabelled and
*      renumbered levels and groups
*   5. the command leaves the data and their sort order as they were
*   6. the table redisplayed after svy equals the table at estimation
* The staged runs use sequentially fixothermargins, with
* bayesboot() replicating every stage; the raking options (trimming
* bounds, iterate()). The manual Break checklist is in tests/README.md.
*
* Run from the tests directory:   . do test_invariance.do

version 16
clear
set more off
set varabbrev off
discard             // drop any seldecomp already in memory; the repository copy is what gets tested

* ---- tolerances ---------------------------------------------------------------
local tol   1e-6    // runs whose raking differs in the order of margins (the raking's tolerance)
local itol  1e-10   // the same computation up to the order of floating-point sums

* ---- helpers ------------------------------------------------------------------
capture program drop keep_e
program define keep_e
    // keep_e stub  --  copies e(b), e(bystage), e(V) (when posted) to stub_b, stub_G, stub_V,
    // and every e() scalar to the 1 x k matrix stub_S (column names = scalar names)
    args stub
    matrix `stub'_b = e(b)
    matrix `stub'_G = e(bystage)
    capture matrix drop `stub'_V
    capture matrix `stub'_V = e(V)
    local sn : e(scalars)
    matrix `stub'_S = J(1, `: word count `sn'', .)
    matrix colnames `stub'_S = `sn'
    local j 0
    foreach s of local sn {
        matrix `stub'_S[1, `++j'] = e(`s')
    }
end

capture program drop cmp_scalars
program define cmp_scalars
    // cmp_scalars stub [, mult(names) lin(names) a(#) by(#) skip(names) only(names) tol(#)]
    // every e() scalar (only those in only() when given) equals the one kept by keep_e stub,
    // except: those in mult() equal by() times it, those in lin() equal a() + by() times it,
    // those in skip() are not compared; missing equals missing
    syntax name(name=stub) [, MULT(namelist) LIN(namelist) A(real 0) BY(real 1) ///
        SKIP(namelist) ONLY(namelist) TOL(real 1e-10)]
    local sn : e(scalars)
    local kn : colnames `stub'_S
    local sn : list sn - skip
    local kn2 : list kn - skip
    if "`only'" != "" {
        local sn : list sn & only
        local kn2 : list kn2 & only
        local miss : list only - sn
        if "`miss'" != "" {
            di as err "cmp_scalars: e(`miss') not posted"
            exit 9
        }
    }
    local s1 : list sort sn
    local s2 : list sort kn2
    if "`s1'" != "`s2'" {
        di as err "cmp_scalars: e() scalars differ from those kept in `stub'_S"
        di as err "  now:  `s1'"
        di as err "  kept: `s2'"
        exit 9
    }
    foreach s of local sn {
        local j : list posof "`s'" in kn
        local k = el(`stub'_S, 1, `j')
        local now = e(`s')
        if missing(`now') | missing(`k') {
            if missing(`now') != missing(`k') {
                di as err "e(`s') = `now', kept `k'"
                exit 9
            }
            continue
        }
        local e = `k'
        if `: list s in mult' local e = `by' * `k'
        if `: list s in lin'  local e = `a' + `by' * `k'
        if reldif(`now', `e') > `tol' {
            di as err "e(`s') = " %20.0g `now' ", expected " %20.0g `e' ///
                "  (reldif " %9.3g reldif(`now', `e') " > `tol')"
            exit 9
        }
    }
end

capture program drop cmp_bystage
program define cmp_bystage
    // cmp_bystage K [, mult() lin() a() by() skip() only() tol()]  --  e(bystage) against the
    // kept matrix K column by column name, with the rules of cmp_scalars
    syntax name(name=K) [, MULT(namelist) LIN(namelist) A(real 0) BY(real 1) ///
        SKIP(namelist) ONLY(namelist) TOL(real 1e-10)]
    tempname G
    matrix `G' = e(bystage)
    local cn : colnames `G'
    local kn : colnames `K'
    if "`cn'" != "`kn'" | rowsof(`G') != rowsof(`K') {
        di as err "cmp_bystage: e(bystage) has other rows or columns than `K'"
        exit 9
    }
    local j 0
    foreach n of local cn {
        local ++j
        if `: list n in skip' continue
        if "`only'" != "" & !`: list n in only' continue
        forvalues i = 1/`=rowsof(`G')' {
            local now = `G'[`i', `j']
            local k = `K'[`i', `j']
            if missing(`now') | missing(`k') {
                if missing(`now') != missing(`k') {
                    di as err "e(bystage)[`i', `n'] = `now', kept `k'"
                    exit 9
                }
                continue
            }
            local e = `k'
            if `: list n in mult' local e = `by' * `k'
            if `: list n in lin'  local e = `a' + `by' * `k'
            if reldif(`now', `e') > `tol' {
                di as err "e(bystage)[`i', `n'] = " %20.0g `now' ", expected " %20.0g `e'
                exit 9
            }
        }
    }
end

capture program drop cmp_vars
program define cmp_vars
    // cmp_vars v1 v2 [#]  --  v1 equals # times v2 (default 1) in every observation, to 1e-10
    args v1 v2 f
    if "`f'" == "" local f 1
    tempvar rd
    qui gen double `rd' = reldif(`v1', `f' * `v2')
    qui count if missing(`v1') != missing(`v2')
    if r(N) {
        di as err "`v1' and `v2' are missing in different observations"
        exit 9
    }
    qui sum `rd', meanonly
    if r(max) > 1e-10 & r(max) < . {
        di as err "`v1' vs `f' * `v2': max reldif " r(max)
        exit 9
    }
end

capture program drop check_data
program define check_data
    // check_data "signature" "sortedby"  --  the data (datasignature), the order of the
    // observations (ord == _n) and the sort order are those recorded before the command
    args sig sb
    qui datasignature
    if "`r(datasignature)'" != "`sig'" {
        di as err "the data changed: `r(datasignature)' (before: `sig')"
        exit 9
    }
    qui count if ord != _n
    if r(N) {
        di as err "the order of the observations changed"
        exit 9
    }
    local now : sortedby
    if "`now'" != "`sb'" {
        di as err "the data are sorted by {`now'}, before by {`sb'}"
        exit 9
    }
end

* _sd_same_tables(logfile): the coefficient rows ("name | number ...") printed after the
* markers SDW_T_A, SDW_T_B and SDW_T_C of the log; sets the locals nrows (the rows after
* SDW_T_A) and same (1 if all numeric tokens match in order; display labels/spacing differ)
capture mata: mata drop _sd_same_tables()
mata:
void _sd_same_tables(string scalar fn)
{
    string colvector t, a, b, c
    real colvector sec, isrow
    real scalar i, same

    t = cat(fn)
    sec = runningsum((strtrim(t) :== "SDW_T_A") :| (strtrim(t) :== "SDW_T_B") :|
        (strtrim(t) :== "SDW_T_C") :| (strtrim(t) :== "SDW_T_D"))
    isrow = J(rows(t), 1, 0)
    for (i = 1; i <= rows(t); i++) {
        isrow[i] = regexm(t[i], "^ *[A-Za-z_][A-Za-z0-9_ ]* [|] +-?[0-9]*[.]?[0-9]")
    }
    // Native svy labels differ from descriptive seldecomp labels. Compare every
    // printed coefficient, SE, test statistic, p-value and CI limit in row order.
    for (i = 1; i <= rows(t); i++) {
        if (isrow[i]) t[i] = strtrim(stritrim(substr(t[i], strpos(t[i], "|") + 1, .)))
    }
    a = select(t, isrow :& (sec :== 1))
    b = select(t, isrow :& (sec :== 2))
    c = select(t, isrow :& (sec :== 3))
    same = (rows(a) == rows(b) & rows(a) == rows(c))
    if (same) same = (all(a :== b) & all(a :== c))
    if (!same) {
        printf("{err}rows at estimation:\n")
        a
        printf("{err}rows after svy:\n")
        b
        printf("{err}rows after seldecomp:\n")
        c
    }
    st_local("nrows", strofreal(rows(a)))
    st_local("same", strofreal(same))
}
end

* ---- locate seldecomp.ado (see run_tests.do) ----------------------------------
if `"$SELDECOMP_ROOT"' != "" adopath ++ `"$SELDECOMP_ROOT"'
else adopath ++ `"`c(pwd)'/.."'
capture which seldecomp
if _rc {
    di as err "seldecomp.ado not found: run this do-file from the tests directory " ///
        "or set the global SELDECOMP_ROOT to the repository root"
    exit 601
}
which seldecomp        // record in the log which copy is being tested

* ---- data -----------------------------------------------------------------------
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
recast double y x3 wt
assert _N == 2000
sort id
* a second continuous variable, correlated with x3 (deterministic: no random numbers)
gen double x6 = mod(id * 7919, 1009) / 1009 + .5 * x3
set seed 20260925


* ==============================================================================
* 1. two binned continuous variables: e(cuts1), e(cuts2)
* ==============================================================================
di as txt _n "---- 1. e(cuts1) and e(cuts2)"
seldecomp y c.x3 c.x6 [pw=wt], by(group) contvars(x3 4 x6 3)
assert "`e(cutvars)'" == "x3 x6"
assert "`e(cutq)'" == "4 3"
assert "`e(contsample)'" == "pooled"
matrix b_c = e(b)
matrix C1 = e(cuts1)
matrix C2 = e(cuts2)
capture confirm matrix e(cuts3)
assert _rc != 0
assert rowsof(C1) == 1 & colsof(C1) == 3 & rowsof(C2) == 1 & colsof(C2) == 2
assert "`: rownames C1'" == "x3" & "`: rownames C2'" == "x6"
assert "`: colnames C2'" == "cut1 cut2"
* the cut points are the weighted quantiles of the pooled sample
_pctile x3 [aw=wt], nq(4)
forvalues j = 1/3 {
    assert reldif(C1[1, `j'], r(r`j')) < `itol'
}
_pctile x6 [aw=wt], nq(3)
forvalues j = 1/2 {
    assert reldif(C2[1, `j'], r(r`j')) < `itol'
}
* the level labels of the second variable are made from e(cuts2)
local c1 = strtrim("`: di %9.0g C2[1, 1]'")
local c2 = strtrim("`: di %9.0g C2[1, 2]'")
assert `"`e(levels2)'"' == "<= `c1' | (`c1', `c2'] | > `c2'"
assert "`e(term2)'" == "x6"
* the same raking on bins made by hand from the stored cut points
gen byte b3 = 4
forvalues j = 3(-1)1 {
    replace b3 = `j' if x3 <= C1[1, `j']
}
gen byte b6 = 3
forvalues j = 2(-1)1 {
    replace b6 = `j' if x6 <= C2[1, `j']
}
seldecomp y i.b3 i.b6 [pw=wt], by(group)
assert mreldif(e(b), b_c) < `itol'
* e(cutvars), e(cutq) and e(cuts#) follow the order of the terms, not of contvars()
seldecomp y c.x6 c.x3 [pw=wt], by(group) contvars(x3 4 x6 3)
assert "`e(cutvars)'" == "x6 x3" & "`e(cutq)'" == "3 4"
matrix T1 = e(cuts1)
matrix T2 = e(cuts2)
assert "`: rownames T1'" == "x6" & "`: rownames T2'" == "x3"
assert mreldif(T1, C2) < `itol' & mreldif(T2, C1) < `itol'
assert mreldif(e(b), b_c) < `tol'
* one term joining two binned variables, cut points from the Focal Group (focalgroup)
seldecomp y c.x3#c.x6 x2 [pw=wt], by(group) contvars(x3 2 x6 2, focalgroup)
assert "`e(cutvars)'" == "x3 x6" & "`e(cutq)'" == "2 2" & "`e(contsample)'" == "focalgroup"
assert "`e(term1)'" == "x3#x6"
matrix b_j = e(b)
_pctile x3 [aw=wt] if group == 0, nq(2)
scalar q3 = r(r1)
assert reldif(el(e(cuts1), 1, 1), q3) < `itol'
_pctile x6 [aw=wt] if group == 0, nq(2)
scalar q6 = r(r1)
assert reldif(el(e(cuts2), 1, 1), q6) < `itol'
local a = strtrim("`: di %9.0g q3'")
local b = strtrim("`: di %9.0g q6'")
assert `"`e(levels1)'"' == "x3=<= `a' x6=<= `b' | x3=<= `a' x6=> `b' | x3=> `a' x6=<= `b' | x3=> `a' x6=> `b'"
gen byte j3 = 1 + (x3 > q3)
gen byte j6 = 1 + (x6 > q6)
egen jc = group(j3 j6)
seldecomp y i.jc x2 [pw=wt], by(group)
assert mreldif(e(b), b_j) < `itol'
drop b3 b6 j3 j6 jc


* ==============================================================================
* 2. level labels of c.x#i.z and i.z#c.x
* ==============================================================================
di as txt _n "---- 2. labels of c.x#i.z"
label define x1lab 1 "one" 2 "two" 3 "three"
label values x1 x1lab
seldecomp y c.x3#i.x1 [pw=wt], by(group) contvars(x3 2)
matrix b_cz = e(b)
assert "`e(term1)'" == "x3#x1" & "`e(cutvars)'" == "x3"
scalar q = el(e(cuts1), 1, 1)
local c = strtrim("`: di %9.0g q'")
local L "x3=<= `c' x1=1 (one) | x3=<= `c' x1=2 (two) | x3=<= `c' x1=3 (three)"
local L "`L' | x3=> `c' x1=1 (one) | x3=> `c' x1=2 (two) | x3=> `c' x1=3 (three)"
assert `"`e(levels1)'"' == "`L'"
* the other order of the operands: the same cells, labelled and ordered by x1 first
seldecomp y i.x1#c.x3 [pw=wt], by(group) contvars(x3 2)
assert "`e(term1)'" == "x1#x3"
local L "x1=1 (one) x3=<= `c' | x1=1 (one) x3=> `c' | x1=2 (two) x3=<= `c'"
local L "`L' | x1=2 (two) x3=> `c' | x1=3 (three) x3=<= `c' | x1=3 (three) x3=> `c'"
assert `"`e(levels1)'"' == "`L'"
assert mreldif(e(b), b_cz) < `itol'
* = the cells made by hand
gen byte bx3 = 1 + (x3 > q)
egen cz = group(bx3 x1)
seldecomp y i.cz [pw=wt], by(group)
assert mreldif(e(b), b_cz) < `itol'
* inside a block of a staged run (sequentially fixothermargins), and a value label containing
* the separator |
label define x1lab 2 "two|2", modify
seldecomp y (a: c.x3#i.x1) (b: x2) [pw=wt], by(group) sequentially fixothermargins contvars(x3 2)
assert strpos(`"`e(levels1)'"', "x1=2 (two/2)") > 0
assert `"`e(block1)'"' == "x3#x1"
label values x1
drop bx3 cz


* ==============================================================================
* 3. aweights and iweights with generate() and bayesboot()
* ==============================================================================
di as txt _n "---- 3. aweights and iweights with generate() and bayesboot()"
* the weight type is a label only: every result, the generate() weights, the replicate
* weights, the saving() file and e(V) equal those with pweights
tempfile bsf
foreach t in pw aw iw {
    local bs_`t' `"`bsf'_`t'.dta"'
    seldecomp y (a: x1) (b: x2 x4) [`t'=wt], by(group) sequentially fixothermargins generate(g_`t') ///
        bayesboot(reps(5) rseed(4) nodots generate(bw_`t') saving(`"`bs_`t''"'))
    local wt_`t' `e(wtype)'
    keep_e w_`t'
}
assert "`wt_pw'" == "pweight" & "`wt_aw'" == "aweight" & "`wt_iw'" == "iweight"
foreach t in aw iw {
    assert mreldif(w_`t'_b, w_pw_b) < `itol'
    assert mreldif(w_`t'_G, w_pw_G) < `itol'
    assert mreldif(w_`t'_V, w_pw_V) < `itol'
    assert mreldif(w_`t'_S, w_pw_S) < `itol'
    cmp_vars g_`t' g_pw
    cmp_vars g_`t'_s1 g_pw_s1
    forvalues r = 1/5 {
        cmp_vars bw_`t'_`r' bw_pw_`r'
    }
    preserve
    use `"`bs_pw'"', clear
    mkmat _all, matrix(F_pw)
    use `"`bs_`t''"', clear
    mkmat _all, matrix(F_`t')
    restore
    assert rowsof(F_`t') == 5
    assert mreldif(F_`t', F_pw) < `itol'
}
foreach t in pw aw iw {
    erase `"`bs_`t''"'
}
drop g_* bw_*
* clusters
foreach t in pw aw iw {
    seldecomp y x1 x2 [`t'=wt], by(group) ///
        bayesboot(reps(5) rseed(5) cluster(cluster) nodots)
    keep_e c_`t'
}
foreach t in aw iw {
    assert mreldif(c_`t'_V, c_pw_V) < `itol'
    assert mreldif(c_`t'_S, c_pw_S) < `itol'
}


* ==============================================================================
* 4. invariance
* ==============================================================================
local spec y (a: x1) (b: x2 x4)
local opts by(group) sequentially fixothermargins bayesboot(reps(5) rseed(21) nodots)
* the scale of the weights: quantities measured in weight units scale with it
local wunits W_focal W_ref raked_mean raked_min raked_max raked_sd
* convergence measures of the raking, which depend on the last iteration's rounding
local conv maxreldif maxctrl

di as txt _n "---- 4a. scaling every weight by a constant"
seldecomp `spec' [pw=wt], `opts' generate(sw0)
keep_e s0
foreach c in 1000 .001 7.25 {
    capture drop wc
    gen double wc = `c' * wt
    capture drop swc swc_s1
    seldecomp `spec' [pw=wc], `opts' generate(swc)
    assert mreldif(e(b), s0_b) < `itol'
    assert mreldif(e(V), s0_V) < `itol'          // the bootstrap draws too (same seed)
    cmp_scalars s0, mult(`wunits') by(`c') skip(`conv')
    cmp_bystage s0_G, mult(`wunits') by(`c') skip(`conv')
    cmp_vars swc sw0 `c'
    cmp_vars swc_s1 sw0_s1 `c'
}
* aweights and iweights, and a constant weight = no weights
foreach t in aw iw {
    seldecomp `spec' [`t'=wt], by(group) sequentially fixothermargins
    keep_e t0
    seldecomp `spec' [`t'=wc], by(group) sequentially fixothermargins
    assert mreldif(e(b), t0_b) < `itol'
    cmp_scalars t0, mult(`wunits') by(7.25) skip(`conv')
}
seldecomp `spec', by(group) sequentially fixothermargins
keep_e n0
gen double w3 = 3
seldecomp `spec' [pw=w3], by(group) sequentially fixothermargins
assert mreldif(e(b), n0_b) < `itol'
cmp_scalars n0, mult(`wunits') by(3) skip(`conv')
* absolute trimming bounds scale with the weights (the bounds leave margins unmet: force)
seldecomp y x1 x2 x4 [pw=wt], by(group) trimhiabs(3) trimloabs(.5) force
keep_e a0
assert e(N_at_bound) > 0
seldecomp y x1 x2 x4 [pw=wc], by(group) trimhiabs(21.75) trimloabs(3.625) force
assert mreldif(e(b), a0_b) < `itol'
cmp_scalars a0, mult(`wunits') by(7.25) skip(`conv')
drop sw0 sw0_s1 swc swc_s1 wc w3

di as txt _n "---- 4b. an affine transform of the outcome"
* y -> A + B*y with B < 0: the means go to A + B*mean, the differences to B*difference,
* the shares are unchanged, V to B^2 V, and the percentile interval is reflected
local A 3
local B -2.5
gen double ya = `A' + `B' * y
seldecomp `spec' [pw=wt], `opts'
keep_e y0
matrix y0_P = e(ci_percentile)
seldecomp ya (a: x1) (b: x2 x4) [pw=wt], `opts'
matrix bb = e(b)
matrix PP = e(ci_percentile)
* (the means: the overall equation's three and the raked means by stage, raked_focal:)
local cn : colnames bb
local ce : coleq bb
local j 0
foreach n of local cn {
    local ++j
    local e : word `j' of `ce'
    local f = cond(inlist("`n'", "focal_group", "raked_focal", "ref_group") | "`e'" == "raked_focal", `A', 0)
    assert reldif(bb[1, `j'], `f' + `B' * y0_b[1, `j']) < `itol'
    assert reldif(PP[1, `j'], `f' + `B' * y0_P[2, `j']) < `itol'
    assert reldif(PP[2, `j'], `f' + `B' * y0_P[1, `j']) < `itol'
}
local B2 = (`B')^2
matrix VV = `B2' * y0_V
assert mreldif(e(V), VV) < `itol'
local means mean_focal mean_ref mean_raked_focal mean_focal_cs mean_focal_zero
local diffs difference explained unexplained
local same pct_explained pct_increment pct_explained_se N N_focal N_ref W_focal W_ref converged ///
    badcontrols N_at_bound N_raked N_focal_zero raked_mean raked_min raked_max raked_sd raked_cv ///
    factor_mean factor_min factor_max factor_sd factor_cv N_eff_base N_eff_raked ess_ratio
cmp_scalars y0, lin(`means') mult(`diffs') a(`A') by(`B') only(`means' `diffs' `same')
cmp_bystage y0_G, lin(raked_focal mean_focal_zero mean_focal_cs) mult(cum_expl cum_unexpl increment) ///
    a(`A') by(`B') skip(`conv')
drop ya

di as txt _n "---- 4c. duplicated observations: iweights as frequency weights"
* each observation repeated k = 1, 2 or 3 times with its weight gives the results of the
* original data with k times the weight; with cluster(id) the Bayesian bootstrap draws one
* multiplier per original observation, so e(V) agrees as well
gen byte k = 1 + mod(id, 3)
gen double kw = k * wt
seldecomp `spec' [iw=kw], by(group) sequentially fixothermargins generate(dg0) ///
    bayesboot(reps(5) rseed(23) cluster(id) nodots)
keep_e d0
preserve
expand k
sort id
seldecomp `spec' [iw=wt], by(group) sequentially fixothermargins generate(dg1) ///
    bayesboot(reps(5) rseed(23) cluster(id) nodots)
assert _N == 4001 & e(N) == 4001 & e(N_clust) == 2000
assert mreldif(e(b), d0_b) < `itol'
assert mreldif(e(V), d0_V) < `itol'
cmp_scalars d0, only(W_focal W_ref mean_focal mean_ref mean_raked_focal difference explained unexplained ///
    pct_explained pct_explained_se converged N_at_bound N_clust)
* the raked weights of the copies add up to the raked weight of the original
collapse (sum) dg1 dg1_s1 (first) dg0 dg0_s1, by(id)
cmp_vars dg1 dg0
cmp_vars dg1_s1 dg0_s1
restore
* the same for pweights and aweights (point estimates), and every row twice with iweights
foreach t in pw aw {
    seldecomp `spec' [`t'=kw], by(group) sequentially fixothermargins
    matrix k0 = e(b)
    preserve
    expand k
    seldecomp `spec' [`t'=wt], by(group) sequentially fixothermargins
    assert mreldif(e(b), k0) < `itol'
    restore
}
seldecomp `spec' [iw=wt], by(group) sequentially fixothermargins
matrix k0 = e(b)
preserve
expand 2
seldecomp `spec' [iw=wt], by(group) sequentially fixothermargins
assert mreldif(e(b), k0) < `itol' & e(N) == 4000
restore
drop k kw dg0 dg0_s1

di as txt _n "---- 4d. the order of the rows"
local spec2 y (a: x1) (b: c.x3 x4)
local opts2 by(group) sequentially fixothermargins contvars(x3 4) bayesboot(reps(5) rseed(24) nodots)
sort id
seldecomp `spec2' [pw=wt], `opts2' generate(ro0)
keep_e r0
matrix r0_C = e(cuts1)
gen double u = runiform()
sort u
seldecomp `spec2' [pw=wt], `opts2' generate(ro1)
* point estimates, stages, cut points and raked weights do not depend on the order
assert mreldif(e(b), r0_b) < `itol'
cmp_bystage r0_G, skip(`conv')
cmp_scalars r0, skip(`conv' pct_explained_se)
matrix T1 = e(cuts1)
assert mreldif(T1, r0_C) < `itol'
sort id
cmp_vars ro1 ro0
cmp_vars ro1_s1 ro0_s1
* the Bayesian bootstrap draws follow the data order: in the original order the same seed
* gives the same draws
seldecomp `spec2' [pw=wt], `opts2'
assert mreldif(e(V), r0_V) < `itol'
cmp_scalars r0, skip(`conv')
drop u ro0 ro0_s1 ro1 ro1_s1

di as txt _n "---- 4e. relabelled and renumbered levels and groups"
seldecomp `spec' [pw=wt], by(group) sequentially fixothermargins
keep_e l0
* x1: codes 1 2 3 -> 30 10 20 with value labels; x2: 0/1 -> the strings "no"/"yes"; x4 (a
* string): numeric codes in the reverse order; each under another name
gen int x1r = cond(x1 == 1, 30, cond(x1 == 2, 10, 20))
label define x1rl 10 "ten" 20 "twenty" 30 "thirty"
label values x1r x1rl
gen str3 x2r = cond(x2 == 1, "yes", "no")
encode x4, gen(x4e)
gen byte x4r = 5 - x4e
seldecomp y (a: x1r) (b: x2r x4r) [pw=wt], by(group) sequentially fixothermargins
assert mreldif(e(b), l0_b) < `itol'
cmp_bystage l0_G, skip(`conv')
cmp_scalars l0, skip(`conv')
assert `"`e(levels1)'"' == "10 (ten) | 20 (twenty) | 30 (thirty)"
assert `"`e(levels2)'"' == "no | yes"
assert `"`e(levels3)'"' == "1 | 2 | 3 | 4"
* the groups under other codes (named by focalgroup()), as strings, labelled
gen int grp = cond(group == 0, 7, 3)
seldecomp `spec' [pw=wt], by(grp) focalgroup(7) sequentially fixothermargins
assert mreldif(e(b), l0_b) < `itol'
assert "`e(focal_group)'" == "7" & "`e(ref_group)'" == "3"
cmp_scalars l0, skip(`conv')
gen str5 grs = cond(group == 0, "focal", "ref")
seldecomp `spec' [pw=wt], by(grs) sequentially fixothermargins
assert mreldif(e(b), l0_b) < `itol'
assert "`e(focal_group)'" == "focal"
label define grl 0 "zero" 1 "one"
label values group grl
seldecomp `spec' [pw=wt], by(group) sequentially fixothermargins
assert mreldif(e(b), l0_b) < `itol'
label values group
drop x1r x2r x4e x4r grp grs


* ==============================================================================
* 5. the data and their sort order are left as they were
* ==============================================================================
di as txt _n "---- 5. data and sort order unchanged"
gen x1c = x1                     // a level found only in the Reference Group: support failure (459)
replace x1c = 9 if group == 1 & mod(id, 11) == 0
tempfile bsf5
local k5 0
local cmd`++k5' `"seldecomp y x1 x2 x4 [pw=wt], by(group)"'
local cmd`++k5' `"seldecomp y (a: x1) (b: c.x3 x4) [aw=wt], by(group) sequentially fixothermargins contvars(x3 4) covariatedetail generate(tg)"'
local cmd`++k5' `"seldecomp y x1 x2 [iw=wt], by(group) bayesboot(reps(5) rseed(1) cluster(cluster) nodots generate(tb) saving("`bsf5'.dta", replace))"'
local cmd`++k5' `"seldecomp y x1 [pw=wt], by(group) bayesboot(reps(5) rseed(1) bca cluster(cluster) nodots)"'
local cmd`++k5' `"seldecomp y (a: x1) (b: x2) [pw=wt], by(group) sequentially fixothermargins makeplot(, twoway(name(tplot, replace)))"'
local cmd`++k5' `"seldecomp y x1c x2 [pw=wt], by(group)"'
local cmd`++k5' `"seldecomp y x1 x2 x4 [pw=wt], by(group) iterate(3)"'
local cmd`++k5' `"seldecomp, covariatedetail diagnostics"'
local cmd`++k5' `"by x2: seldecomp y x1 [pw=wt], by(group)"'
local rc5 0 0 0 0 0 459 430 0 0
* first in a random order that is not a sort order, then sorted by x2 id
gen double u = runiform()
sort u
drop u
gen long ord = _n
foreach order in none "x2 id" {
    if "`order'" != "none" {
        sort `order'
        replace ord = _n
    }
    local sb = cond("`order'" == "none", "", "`order'")
    forvalues c = 1/`k5' {
        if `c' == `k5' & "`order'" == "none" continue      // the by prefix needs sorted data
        qui datasignature
        local sig `r(datasignature)'
        capture noisily `cmd`c''
        local rc = _rc
        if `rc' != `: word `c' of `rc5'' {
            di as err `"`cmd`c'': return code `rc'"'
            exit 9
        }
        capture drop tg*
        capture drop tb_*
        check_data "`sig'" "`sb'"
    }
}
capture graph drop tplot
capture erase `"`bsf5'.dta"'
drop x1c ord


* ==============================================================================
* 6. the table redisplayed after svy equals the table at estimation
* ==============================================================================
di as txt _n "---- 6. redisplay after svy"
* 30 strata of 2 PSUs (pairs of the 60 clusters), BRR half-samples from a Hadamard matrix
gen int hs = ceil(cluster / 2)
mata: Hd = (1); for (i = 1; i <= 5; i++) Hd = (Hd, Hd \ Hd, -Hd); st_matrix("H32", Hd)
svyset cluster [pw=wt], strata(hs)
tempfile lg
local lgf `"`lg'.log"'
capture log close ilog
* (test_display.do compares the rows after svy jackknife)
foreach call in "svy brr, hadamard(H32): seldecomp y (a: x1) (b: x2 x4), by(group) sequentially fixothermargins" ///
    "svy brr, hadamard(H32) fay(.5): seldecomp y x1 x2, by(group)" {
    log using `"`lgf'"', text replace name(ilog)
    di "SDW_T_A"
    `call'
    di "SDW_T_B"
    svy
    di "SDW_T_C"
    * with inference seldecomp's own table shows every equation
    seldecomp
    di "SDW_T_D"
    log close ilog
    * the coefficient rows at estimation, after svy (replay) and after seldecomp (replay)
    mata: _sd_same_tables(st_local("lgf"))
    di as txt `"`call': `nrows' coefficient rows, identical: `same'"'
    assert `nrows' == colsof(e(b)) & `same' == 1
    * each renderer centres the variance type over its standard-error column.
    mata: st_local("nbrr", strofreal(sum(strrtrim(cat(st_local("lgf"))) :== "             |                 BRR")))
    assert `nbrr' == 2
    mata: st_local("nbrr", strofreal(sum(strrtrim(cat(st_local("lgf"))) :== "                 |                BRR")))
    assert `nbrr' == 1
}
svyset, clear
capture erase `"`lgf'"'
drop hs

di as res _n "PASSED: test_invariance.do"

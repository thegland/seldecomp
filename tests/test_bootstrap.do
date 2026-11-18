* test_bootstrap.do -- bayesboot() Bayesian bootstrap of seldecomp: e(V), the intervals and
* the replicate statistics; cluster(), strata() refused; saving() and generate() against plain
* runs on the replicate weights; bca and its jackknife; failed replicates, and those kept at
* the iterate() limit; a replay at another level; estat bootstrap. The staged runs use
* sequentially (the blocks not yet added free, the reference's case F) and sequentially
* fixothermargins (held at the Focal Group's margins, case G), with bayesboot() replicating
* every stage; the words varbyvar and sequential inside bayesboot(), and the options notable,
* noheader and nodefinitions (and noheader of estat bootstrap) are refused.
*
* Run from the tests directory:   . do test_bootstrap.do

version 16
clear
set more off
set varabbrev off
discard             // drop any seldecomp already in memory; the repository copy is what gets tested

* ---- tolerances ---------------------------------------------------------------
local tol  1e-6     // point estimates vs reference
local itol 1e-8     // identities among the bootstrap results
local rtol 1e-12    // identical reruns with the same seed

* ---- helpers ------------------------------------------------------------------
capture program drop assert_close
program define assert_close
    // assert_close exp1 exp2 [tol] ["label"]  --  asserts reldif(exp1, exp2) <= tol
    args a b tol name
    if "`tol'" == "" local tol 1e-6
    if "`name'" == "" local name "`a' vs `b'"
    if missing(`a') | missing(`b') {
        di as err "`name': missing value (`a' = " `a' ", `b' = " `b' ")"
        exit 9
    }
    if reldif(`a', `b') > `tol' {
        di as err "`name': " %18.0g `a' " vs " %18.0g `b' ///
            "  (reldif = " %9.3g reldif(`a', `b') " > `tol')"
        exit 9
    }
end

capture program drop row_vs_e
program define row_vs_e
    // row_vs_e M r tol  --  row r of matrix M, whose columns carry the variable names of a
    // bayesboot(saving()) file, equals the results of the last (plain) seldecomp run
    args M r tol
    tempname S
    matrix `S' = e(bystage)
    local names : colnames `M'
    local j 0
    foreach v of local names {
        local ++j
        local col
        if regexm("`v'", "^(marg|cumexpl|cumunexpl|raked)_s([0-9]+)$") {
            local s = regexs(2)
            local col = cond(regexs(1) == "marg", "increment", cond(regexs(1) == "cumexpl", "cum_expl", ///
                cond(regexs(1) == "cumunexpl", "cum_unexpl", "raked_focal")))
        }
        if "`col'" != "" {
            local x = `S'[`s', colnumb(`S', "`col'")]
        }
        else if "`v'" == "focal_group" local x = e(mean_focal)
        else if "`v'" == "raked_focal" local x = e(mean_raked_focal)
        else if "`v'" == "ref_group" local x = e(mean_ref)
        else if inlist("`v'", "difference", "explained", "unexplained", "pct_explained") local x = e(`v')
        else {
            di as err "row_vs_e: unexpected column `v'"
            exit 9
        }
        assert_close `M'[`r',`j'] `x' `tol' "saved replicate `r', `v', vs a plain run"
    }
end

capture program drop jk_accel
program define jk_accel
    // jk_accel JK A  --  A = jackknife acceleration of each column of JK (rows = deletions;
    // rows with missing values, the failed deletions, are left out):
    // a = sum(d^3) / (6 * sum(d^2)^1.5), d = mean(theta_(-i)) - theta_(-i)
    args JK A
    mata: st_matrix("`A'", tb_jk_accel(st_matrix("`JK'")))
end
capture mata: mata drop tb_jk_accel()
mata:
real rowvector tb_jk_accel(real matrix JK)
{
    real rowvector a
    real colvector v, d
    real scalar j
    a = J(1, cols(JK), .)
    for (j = 1; j <= cols(JK); j++) {
        v = select(JK[., j], JK[., j] :< .)
        d = mean(v) :- v
        a[j] = sum(d :^ 3) / (6 * sum(d :^ 2)^1.5)
    }
    return(a)
}
end

capture program drop log_has
program define log_has, rclass
    // log_has "file" "text"  --  r(found) = 1 if the text log contains the text once its
    // lines are joined and runs of blanks squeezed (so that a wrapped note still matches);
    // r(count): how many times
    args file text
    mata: _t1s = stritrim(invtokens(cat(st_local("file"))')); _t1t = st_local("text")
    mata: st_local("n", strofreal((strlen(_t1s) - strlen(subinstr(_t1s, _t1t, ""))) / strlen(_t1t)))
    capture mata: mata drop _t1s _t1t
    return scalar found = (`n' > 0)
    return scalar count = `n'
end

capture program drop mat_same
program define mat_same
    // mat_same A B tol  --  A and B have the same dimensions, missing values in the same
    // cells, and reldif() <= tol in the others
    args A B tol
    mata: st_local("ok", strofreal(tb_mat_same(st_matrix("`A'"), st_matrix("`B'"), `tol')))
    if !`ok' {
        di as err "mat_same: `A' and `B' differ (tolerance `tol')"
        matrix list `A', format(%20.0g)
        matrix list `B', format(%20.0g)
        exit 9
    }
end
capture mata: mata drop tb_mat_same()
mata:
real scalar tb_mat_same(real matrix A, real matrix B, real scalar tol)
{
    real colvector a, b
    if (rows(A) != rows(B) | cols(A) != cols(B)) return(0)
    a = vec(A)
    b = vec(B)
    if (sum((a :>= .) :!= (b :>= .)) > 0) return(0)
    if (sum(a :< .) == 0) return(1)
    return(max(reldif(select(a, a :< .), select(b, a :< .))) <= tol)
}
end

* ---- locate seldecomp.ado -----------------------------------------------------
* the repository copy of seldecomp.ado always goes first on the ado-path (see the -discard-
* at the top of this file)
if `"$SELDECOMP_ROOT"' != "" adopath ++ `"$SELDECOMP_ROOT"'
else adopath ++ `"`c(pwd)'/.."'
capture which seldecomp
if _rc {
    di as err "seldecomp.ado not found: run this do-file from the tests directory " ///
        "or set the global SELDECOMP_ROOT to the repository root"
    exit 601
}
which seldecomp        // record in the log which copy is being tested

* ---- data and expected values -------------------------------------------------
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
recast double y x3 wt
confirm string variable x4
assert _N == 2000
do expected_values.do
sum wt, meanonly
local W = r(sum)


* ==============================================================================
* 1. basic run: reps(30) rseed(12345)
* ==============================================================================
di as txt _n "---- 1. bayesboot(reps(30) rseed(12345))"
seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(30) rseed(12345))

assert "`e(vce)'"        == "bayesboot"
assert "`e(vcetype)'"    == "Bootstrap"
assert "`e(properties)'" == "b V"
assert "`e(clustvar)'"   == ""
assert "`e(stratvar)'"   == ""
assert "`e(rngstate)'"   != ""
assert e(N_reps)        == 30                 // completed replicates, as in -bootstrap-
assert e(N_misreps)     == 0                  // failed replicates
assert e(N_reps) + e(N_misreps) == 30         // together: the requested number
assert e(N_iterlimit)   == 0                  // kept although stopped at iterate()
assert "`e(missing)'"   == ""
assert e(seed)          == 12345
assert e(bca)           == 0
assert e(ties)          == 0
assert e(level)         == 95
assert wordcount("`e(bsnames)'") == 6
assert e(N) == 2000

* point estimates are unaffected by the bootstrap
assert_close e(explained)    exp_C_explained `tol' "explained with bayesboot"
assert_close e(mean_raked_focal) exp_C_raked_1   `tol' "raked mean with bayesboot"

matrix b  = e(b)
matrix V1 = e(V)
assert colsof(b) == 6
assert rowsof(V1) == 6 & colsof(V1) == 6
forvalues j = 1/6 {
    assert el(V1, `j', `j') > 0 & el(V1, `j', `j') < .
}
matrix V1t = V1'
assert mreldif(V1, V1t) < `itol'                       // symmetric
* stripes of V match e(b)
local bn : colfullnames b
local vn : colfullnames V1
assert "`bn'" == "`vn'"

foreach m in ci_normal ci_percentile ci_bc {
    matrix M = e(`m')
    assert rowsof(M) == 2 & colsof(M) == 6
    forvalues j = 1/6 {
        assert el(M, 1, `j') < el(M, 2, `j')             // lower < upper
        assert el(M, 1, `j') < . & el(M, 2, `j') < .
    }
}
foreach m in b_bs bias se z0 reps {
    matrix M = e(`m')
    assert rowsof(M) == 1 & colsof(M) == 6
}
* e(reps): completed replicates per bootstrapped element (= e(N_reps), which
* counts completed replicates; e(N_misreps) counts the failed ones)
matrix R = e(reps)
forvalues j = 1/6 {
    assert el(R, 1, `j') == e(N_reps)
}
local rn : colfullnames R
assert "`rn'" == "`e(bsnames)'"
* explained share: bootstrap SE and percentile CI (not part of e(b))
assert e(pct_explained_se) > 0 & e(pct_explained_se) < .
matrix PC = e(pct_explained_ci)
assert rowsof(PC) == 1 & colsof(PC) == 2
assert el(PC, 1, 1) < el(PC, 1, 2)
local pcn : colnames PC
assert "`pcn'" == "ll ul"
* no BCa results without bca
capture confirm matrix e(ci_bca)
assert _rc != 0
capture confirm matrix e(accel)
assert _rc != 0

* SE = sqrt(diag(V)); normal CI = b -/+ z * SE; bias = mean of replicates - observed
matrix SE = e(se)
matrix CN = e(ci_normal)
matrix BB = e(b_bs)
matrix BI = e(bias)
local z = invnormal(0.975)
forvalues j = 1/6 {
    assert el(SE, 1, `j') > 0
    assert reldif(el(SE, 1, `j'), sqrt(el(V1, `j', `j'))) < `itol'
    assert reldif(el(CN, 1, `j'), el(b, 1, `j') - `z' * el(SE, 1, `j')) < `itol'
    assert reldif(el(CN, 2, `j'), el(b, 1, `j') + `z' * el(SE, 1, `j')) < `itol'
    assert reldif(el(BI, 1, `j'), el(BB, 1, `j') - el(b, 1, `j')) < `itol'
}
* replicate-level identities: difference = focal_group - ref_group = explained + unexplained
assert reldif(el(V1,4,4), el(V1,1,1) + el(V1,3,3) - 2 * el(V1,1,3)) < `itol'
assert reldif(el(V1,4,4), el(V1,5,5) + el(V1,6,6) + 2 * el(V1,5,6)) < `itol'
* the percentile interval contains the mean of the replicates
matrix CP = e(ci_percentile)
forvalues j = 1/6 {
    assert el(CP, 1, `j') <= el(BB, 1, `j') & el(BB, 1, `j') <= el(CP, 2, `j')
}
* replay with bootstrap results
seldecomp
seldecomp, level(90)
* the words notable, noheader and nodefinitions are refused (198) at replay
foreach o in notable noheader nodefinitions {
    rcof "noi seldecomp, `o'" == 198
}
ereturn list


* ==============================================================================
* 2. reproducibility: the same seed gives identical results, another seed does not
* ==============================================================================
di as txt _n "---- 2. seed"
seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(30) rseed(12345))
matrix V2 = e(V)
matrix b2 = e(b)
assert mreldif(V1, V2) < `rtol'
assert mreldif(b, b2) < `rtol'
matrix P2 = e(ci_percentile)
matrix P1 = CP
assert mreldif(P1, P2) < `rtol'

seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(30) rseed(54321))
assert e(seed) == 54321
matrix V3 = e(V)
assert mreldif(V1, V3) > 1e-6

* no rseed(): a seed is drawn and stored, so the run can be reproduced from e(seed)
seldecomp y x1 [pw=wt], by(group) bayesboot(reps(10) nodots)
assert e(seed) < . & e(seed) >= 1
local seed = e(seed)
matrix V4 = e(V)
seldecomp y x1 [pw=wt], by(group) bayesboot(reps(10) nodots rseed(`seed'))
matrix V5 = e(V)
assert mreldif(V4, V5) < `rtol'


* ==============================================================================
* 3. cluster(); strata() is refused
* ==============================================================================
di as txt _n "---- 3. cluster / strata() refused"
seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(30) rseed(12345) cluster(cluster))
assert e(N_clust) == 60
assert "`e(clustvar)'" == "cluster"
assert e(N_reps) == 30
assert missing(e(N_strata))                            // never stored by bayesboot()
assert "`e(stratvar)'" == ""                           // nor is e(stratvar)
matrix Vc = e(V)
forvalues j = 1/6 {
    assert el(Vc, `j', `j') > 0 & el(Vc, `j', `j') < .
}

* bayesboot(strata()), in full or abbreviated (the name is echoed as typed), is refused (198)
* with a pointer to svy jackknife and svy brr
tempfile stlog
log using `"`stlog'"', text replace name(bsstrata)
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(30) rseed(1) strata(strata) cluster(cluster))" == 198
rcof "noi seldecomp y x1 [pw=wt], by(group) bayesboot(reps(2) nodots str(strata))" == 198
log close bsstrata
log_has `"`stlog'"' "bayesboot(): option strata() not allowed; for a stratified design use svy jackknife: or svy brr:"
assert r(count) == 1
log_has `"`stlog'"' "bayesboot(): option str() not allowed; for a stratified design use svy jackknife: or svy brr:"
assert r(count) == 1
log_has `"`stlog'"' "not allowed; for a stratified design use svy jackknife: or svy brr:"
assert r(count) == 2

* replicate weights: within a cluster the multiplier is constant; the total is kept
seldecomp y x1 [pw=wt], by(group) bayesboot(reps(2) rseed(9) cluster(cluster) generate(bc) nodots)
confirm variable bc_1 bc_2
gen double fac = bc_1 / wt
bysort cluster (fac): gen double facrng = reldif(fac[_N], fac[1])
assert facrng < `itol'
sum bc_1, meanonly
local sbc = r(sum)
assert_close `sbc' `W' `itol' "sum of clustered replicate weights"
drop fac facrng bc_1 bc_2


* ==============================================================================
* 4. saving() and generate()
* ==============================================================================
di as txt _n "---- 4. saving() / generate()"
capture erase bsreps.dta
seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(20) rseed(7) saving(bsreps.dta, replace) generate(bw))
assert e(N_reps) == 20
assert "`e(bsfile)'" == "bsreps.dta"
assert "`e(bsgenerate)'" == "bw"
forvalues r = 1/20 {
    confirm variable bw_`r'
    assert "`: type bw_`r''" == "double"
}
capture confirm variable bw_21
assert _rc == 111
* replicate weights are positive on the whole estimation sample and keep the total weight
count if !(bw_1 > 0 & bw_1 < .)
assert r(N) == 0
sum bw_1, meanonly
local sbw = r(sum)
assert_close `sbw' `W' `itol' "sum of replicate weights = sum of wt"
sum bw_20, meanonly
local sbw = r(sum)
assert_close `sbw' `W' `itol' "sum of replicate weights = sum of wt (rep 20)"
* different replicates have different weights
count if reldif(bw_1, bw_2) > 1e-6
assert r(N) > 1000
* the replicate file: 20 rows, one variable per bootstrapped statistic plus the explained share
confirm file bsreps.dta
describe using bsreps.dta
assert r(N) == 20
assert r(k) == 7
matrix BB = e(b_bs)
matrix SE = e(se)
matrix sbb = e(b)
local pct_se = e(pct_explained_se)
preserve
use bsreps.dta, clear
confirm variable focal_group raked_focal ref_group difference explained unexplained pct_explained
assert _N == 20
sum explained
local rmean = r(mean)
local rsd   = r(sd)
assert_close `rmean' el(BB,1,5) `itol' "mean of saved replicates = e(b_bs)"
assert_close `rsd'   el(SE,1,5) `itol' "sd of saved replicates = e(se)"
assert reldif(explained + unexplained, difference) < 1e-10
sum pct_explained
local psd = r(sd)
assert_close `psd' `pct_se' `itol' "sd of saved explained shares = e(pct_explained_se)"
assert reldif(pct_explained, 100 * explained / difference) < 1e-10
* the file is a version-3 bootstrap dataset: characteristics that -bstat- reads
assert "`: char _dta[bs_version]'" == "3"
assert "`: char _dta[N]'" == "2000"
assert "`: char _dta[N_reps]'" == "20"
assert "`: char _dta[command]'" == "seldecomp"
assert "`: char _dta[seed]'" == "7"
local obs : char explained[observed]
assert reldif(`obs', el(sbb, 1, 5)) < `itol'
assert "`: char explained[coleq]'" == "overall"
assert "`: char explained[colname]'" == "explained"
assert "`: char explained[expression]'" != ""
local obs : char pct_explained[observed]
assert "`obs'" != ""
assert reldif(`obs', 100 * el(sbb, 1, 5) / el(sbb, 1, 4)) < `itol'
restore
* -bstat using- accepts the file and reproduces the bootstrap results
bstat using bsreps.dta, notable noheader
assert "`e(cmd)'" == "bstat"
assert e(N_reps) == 20
assert e(N) == 2000
matrix bb = e(b)
matrix bse = e(se)
local c = colnumb(bb, "overall:explained")
assert `c' < .
assert reldif(el(bb, 1, `c'), el(sbb, 1, 5)) < `itol'
assert reldif(el(bse, 1, `c'), el(SE, 1, 5)) < `itol'
* saving() without replace refuses to overwrite
rcof "noi seldecomp y x1 [pw=wt], by(group) bayesboot(reps(2) saving(bsreps.dta))" == 602
erase bsreps.dta
* generate() refuses existing names; generate(, replace) overwrites those of the new replicates
* (2 of the 20 here: the others stay), and refuses the inputs of the command
rcof "noi seldecomp y x1 [pw=wt], by(group) bayesboot(reps(2) generate(bw))" == 110
gen double bw1_old = bw_1
gen double bw3_old = bw_3
seldecomp y x1 [pw=wt], by(group) bayesboot(reps(2) rseed(5) nodots generate(bw, replace))
assert "`e(bsgenerate)'" == "bw"
count if reldif(bw_1, bw1_old) > 1e-6
assert r(N) > 0
assert bw_3 == bw3_old
drop bw1_old bw3_old
gen double wtbs_1 = wt
rcof "noi seldecomp y x1 [pw=wtbs_1], by(group) bayesboot(reps(2) generate(wtbs, replace))" == 198
gen byte xr_2 = x1
rcof "noi seldecomp y xr_2 [pw=wt], by(group) bayesboot(reps(2) generate(xr, replace))" == 198
drop wtbs_1 xr_2
drop bw_*
capture confirm variable bw_1
assert _rc == 111


* ==============================================================================
* 5. ties, nodots, dots(), level(), bayesboot alone
* ==============================================================================
di as txt _n "---- 5. ties / dots / level"
seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(20) rseed(3) ties nodots)
assert e(ties) == 1
assert e(N_reps) == 20
matrix BC = e(ci_bc)
assert rowsof(BC) == 2 & colsof(BC) == 6

seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(20) rseed(3) dots(5))
assert e(N_reps) == 20

* the confidence level: bayesboot(level()) sets the bootstrap intervals; level()
* outside bayesboot() is an error (a replication prefix takes its own), and so is a level
* outside 10 to 99.99
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) level(90) bayesboot(reps(20) rseed(3))" == 198
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(20) rseed(3) level(9))" == 198
seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(20) rseed(3) level(90))
assert e(level) == 90
matrix b  = e(b)
matrix SE = e(se)
matrix CN = e(ci_normal)
local z = invnormal(0.95)
forvalues j = 1/6 {
    assert reldif(el(CN, 1, `j'), el(b, 1, `j') - `z' * el(SE, 1, `j')) < `itol'
    assert reldif(el(CN, 2, `j'), el(b, 1, `j') + `z' * el(SE, 1, `j')) < `itol'
}
seldecomp, level(95)
* abbreviated
seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(20) rseed(3) l(80))
assert e(level) == 80

* bayesboot alone: 100 replications
seldecomp y x1 [pw=wt], by(group) bayesboot
assert e(N_reps) == 100
assert e(N_misreps) == 0
assert e(seed) < .
assert wordcount("`e(bsnames)'") == 6


* ==============================================================================
* 6. sequentially: every stage bootstrapped; the blocks not yet added are free by default
*    (the reference's case F) and held at the Focal Group's margins with fixothermargins (case G)
* ==============================================================================
di as txt _n "---- 6. sequentially"
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially bayesboot(reps(20) rseed(3))
assert "`e(stagetype)'" == "sequentially" & "`e(fixothermargins)'" == ""
assert e(N_stages) == 3
assert e(N_reps) == 20
matrix b = e(b)
assert colsof(b) == 18
matrix V = e(V)
assert rowsof(V) == 18 & colsof(V) == 18
forvalues j = 1/18 {
    assert el(V, `j', `j') > 0 & el(V, `j', `j') < .
}
assert wordcount("`e(bsnames)'") == 18
foreach m in ci_normal ci_percentile ci_bc {
    matrix M = e(`m')
    assert rowsof(M) == 2 & colsof(M) == 18
}
matrix SE = e(se)
assert colsof(SE) == 18
matrix R = e(reps)
assert colsof(R) == 18
* point estimates are those of the reference's case F: sequentially alone leaves the blocks
* not yet added free, so stage 2 is F's, not G's (where they are held)
assert_close e(explained) exp_F3_explained `tol' "sequentially explained with bayesboot"
matrix S = e(bystage)
assert reldif(el(S, 2, 2), exp_F2_explained) < `tol'
assert reldif(exp_F2_explained, exp_G2_explained) > 1e-3
* the overall explained part and the last cumulative explained part are the same statistic
local bnames : colfullnames b
local c1 : list posof "overall:explained" in bnames
local c2 : list posof "cum_expl:x4" in bnames
assert `c1' == 5 & `c2' == 12
assert reldif(el(V, `c1', `c1'), el(V, `c2', `c2')) < `itol'
* increments sum to the overall explained part in every replicate:
* Var(overall:explained) = Var(sum of increments)
local i1 : list posof "marg_expl:x1" in bnames
local i2 : list posof "marg_expl:x2" in bnames
local i3 : list posof "marg_expl:x4" in bnames
assert `i1' == 7 & `i2' == 8 & `i3' == 9
local vsum = el(V,`i1',`i1') + el(V,`i2',`i2') + el(V,`i3',`i3') ///
    + 2 * (el(V,`i1',`i2') + el(V,`i1',`i3') + el(V,`i2',`i3'))
assert reldif(el(V, `c1', `c1'), `vsum') < `itol'
seldecomp

* fixothermargins with every stage bootstrapped: the blocks not yet added are held at Group
* 1's margins (the reference's case G: stage 2 is G's)
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially fixothermargins bayesboot(reps(10) rseed(3) nodots)
assert "`e(stagetype)'" == "sequentially" & "`e(fixothermargins)'" == "fixothermargins"
assert e(N_reps) == 10
matrix V = e(V)
forvalues j = 1/15 {
    assert el(V, `j', `j') > 0 & el(V, `j', `j') < .
}
assert_close e(explained) exp_G3_explained `tol' "sequentially fixothermargins explained with bayesboot"
matrix S = e(bystage)
assert reldif(el(S, 2, 2), exp_G2_explained) < `tol'

* the words varbyvar and sequential inside bayesboot() are refused (198)
tempfile plog
log using `"`plog'"', text replace name(bsparse)
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) sequentially bayesboot(varbyvar reps(2) nodots)" == 198
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) sequentially bayesboot(sequential reps(2) nodots)" == 198
log close bsparse
log_has `"`plog'"' "bayesboot(): option varbyvar not allowed"
assert r(count) == 1
log_has `"`plog'"' "bayesboot(): option sequential not allowed"
assert r(count) == 1


* ==============================================================================
* 7. bca on a small subsample (delete-one jackknife acceleration)
* ==============================================================================
di as txt _n "---- 7. bca"
preserve
keep if _n <= 150
capture erase bsbca.dta
seldecomp y x1 [pw=wt], by(group) bayesboot(reps(20) rseed(5) bca saving(bsbca.dta, replace))
assert e(bca) == 1
assert e(N_reps) == 20
matrix A = e(accel)
assert rowsof(A) == 1 & colsof(A) == 6
matrix BCA = e(ci_bca)
assert rowsof(BCA) == 2 & colsof(BCA) == 6
forvalues j = 1/6 {
    if !missing(el(BCA, 1, `j')) & !missing(el(BCA, 2, `j')) {
        assert el(BCA, 1, `j') < el(BCA, 2, `j')
    }
}
seldecomp
* outside check: e(accel) equals the acceleration of a delete-one jackknife made by hand, each
* observation dropped in turn with -if- and the rest decomposed by a plain run (a deletion
* that fails or does not converge is left out, as in seldecomp), to 1e-8
matrix JK = J(150, 7, .)
forvalues i = 1/150 {
    capture seldecomp y x1 [pw=wt] if _n != `i', by(group)
    if _rc == 0 & e(converged) == 1 {
        matrix JK[`i', 1] = (e(mean_focal), e(mean_raked_focal), e(mean_ref), e(difference), ///
            e(explained), e(unexplained), e(pct_explained))
    }
}
jk_accel JK AH
forvalues j = 1/6 {
    assert reldif(el(A, 1, `j'), el(AH, 1, `j')) < 1e-8
}
restore
* the saving() file carries the same accelerations, the explained share's included
preserve
use bsbca.dta, clear
local j 0
foreach v in focal_group raked_focal ref_group difference explained unexplained pct_explained {
    local ++j
    local acc : char `v'[acceleration]
    assert reldif(`acc', el(AH, 1, `j')) < 1e-8
}
restore
erase bsbca.dta

* the same with clusters and a staged run: one deletion per cluster, stage columns included
preserve
keep if cluster <= 20
seldecomp y x1 x2 [pw=wt], by(group) sequentially fixothermargins bayesboot(reps(10) rseed(5) bca cluster(cluster) nodots)
assert e(N_clust) == 20
matrix A = e(accel)
assert colsof(A) == 14
local an : colfullnames A
assert "`: word 7 of `an''" == "marg_expl:x1" & "`: word 12 of `an''" == "cum_unexpl:x2"
assert "`: word 14 of `an''" == "raked_focal:x2"
matrix JK = J(20, 14, .)
forvalues c = 1/20 {
    capture seldecomp y x1 x2 [pw=wt] if cluster != `c', by(group) sequentially fixothermargins
    if _rc continue
    matrix S = e(bystage)
    local cv = colnumb(S, "converged")
    if S[1, `cv'] == 0 | S[2, `cv'] == 0 continue
    matrix JK[`c', 1] = (e(mean_focal), e(mean_raked_focal), e(mean_ref), e(difference), e(explained), e(unexplained))
    forvalues s = 1/2 {
        matrix JK[`c', 6 + `s']  = S[`s', colnumb(S, "increment")]
        matrix JK[`c', 8 + `s']  = S[`s', colnumb(S, "cum_expl")]
        matrix JK[`c', 10 + `s'] = S[`s', colnumb(S, "cum_unexpl")]
        matrix JK[`c', 12 + `s'] = S[`s', colnumb(S, "raked_focal")]
    }
}
jk_accel JK AH
forvalues j = 1/14 {
    assert reldif(el(A, 1, `j'), el(AH, 1, `j')) < 1e-8
}
restore

* bca with clusters: one jackknife deletion per cluster
seldecomp y x1 [pw=wt], by(group) bayesboot(reps(10) rseed(5) bca cluster(cluster) nodots)
assert e(bca) == 1
assert e(N_clust) == 60
matrix A = e(accel)
assert colsof(A) == 6


* ==============================================================================
* 8. bootstrap options combined with generate()
* ==============================================================================
di as txt _n "---- 8. combined"
seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(10) rseed(11) nodots) generate(rwb)
confirm variable rwb
sum wt if group == 0, meanonly
local W1 = r(sum)
sum rwb if group == 0, meanonly
local srwb = r(sum)
assert_close `srwb' `W1' `tol' "raked weights with bootstrap"
count if group == 1 & rwb != wt                  // Reference Group keeps its base weight
assert r(N) == 0
drop rwb

* nothing left on disk
capture erase bsreps.dta

* ---- calibrated draws: the SE of a plain mean matches its unbiased estimate
*      (Gamma((m-2)/m,1) draws, deterministic with rseed()) ----
seldecomp y x1 in 1/300, by(group) bayesboot(reps(200) rseed(11) nodots)
matrix SE = e(se)
qui sum y if group == 0 in 1/300
assert reldif(SE[1, 1], r(sd) / sqrt(r(N))) < 0.15


* ==============================================================================
* 9. outside checks: saved replicates against plain runs; failed replicates
* ==============================================================================
di as txt _n "---- 9a. saving() rows = plain runs on the generate() replicate weights"
sort id             // the draws follow the data order: make this section independent of the ones above
* every column of every row of the saving() file equals a plain seldecomp run with that
* replicate's weights (bayesboot(generate())) as pweights, to 1e-10; categorical raking
* variables only (the cut points of c. variables are held fixed across replicates)
capture erase f4reps.dta
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially fixothermargins ///
    bayesboot(reps(6) rseed(41) nodots saving(f4reps.dta, replace) generate(fw))
assert e(N_stages) == 3 & e(N_reps) == 6 & e(N_misreps) == 0
preserve
use f4reps.dta, clear
assert _N == 6 & c(k) == 19
unab svars : _all
mkmat `svars', matrix(SV)
restore
forvalues r = 1/6 {
    qui seldecomp y x1 x2 x4 [pw=fw_`r'], by(group) sequentially fixothermargins
    row_vs_e SV `r' 1e-10
}
drop fw_*
* clusters, one-stage replicates
seldecomp y x1 x2 [pw=wt], by(group) ///
    bayesboot(reps(5) rseed(43) cluster(cluster) nodots saving(f4reps.dta, replace) generate(fc))
preserve
use f4reps.dta, clear
assert _N == 5 & c(k) == 7
unab svars : _all
mkmat `svars', matrix(SV)
restore
forvalues r = 1/5 {
    qui seldecomp y x1 x2 [pw=fc_`r'], by(group)
    row_vs_e SV `r' 1e-10
}
drop fc_*

di as txt _n "---- 9b. failed replicates"
* a replicate fails when its raking misses a target (as the point estimate would stop with
* 430). Focal Group lacks the cell a = 1, b = 1, so its margins can be met only while Reference Group's
* share of (1, 1) stays below its share of (0, 0): 24 against 30 of 200 observations here, so
* the point estimate is feasible, and a replicate whose draws reverse the order has margins
* that cannot be met. The number of failures depends on the draws, so only ranges are
* asserted; at least 2 of 20 (10%), the plural of the note
local R 20
tempfile nlog
local ls = c(linesize)
tempfile f4main
qui save "`f4main'"
clear
set obs 400
gen byte group = _n > 200
gen int k = cond(group == 0, _n, _n - 200)
gen byte a = cond(group == 0, k > 130, cond(k <= 30, 0, cond(k <= 54, 1, !mod(k, 2))))
gen byte b = cond(group == 0, inrange(k, 61, 130), cond(k <= 30, 0, cond(k <= 54, 1, mod(k, 2))))
gen double y = a + 2 * b + mod(_n * 7, 11) / 10
set linesize 120
log using "`nlog'", text replace name(f4note)
noisily seldecomp y a b, by(group) bayesboot(reps(`R') rseed(7) nodots saving(f4reps.dta, replace) generate(ff))
log close f4note
set linesize `ls'
local nm = e(N_misreps)
local nr = e(N_reps)
assert e(converged) == 1 & e(badcontrols) == 0
assert `nr' + `nm' == `R'
assert inrange(`nm', 2, `R' - 2)
assert e(N_iterlimit) == 0
assert "`e(missing)'" == "missing"
matrix RR = e(reps)
forvalues j = 1/6 {
    assert el(RR, 1, `j') == `nr'
}
matrix SE = e(se)
matrix BB = e(b_bs)
matrix DR = e(bsdraws)
matrix DRP = e(bsdraws_pct)
* one note: the count of failed replicates of R (plural: nm >= 2) and what failed them (the
* margins missed; in some the raking stopped at the iterate() limit of 10000 first), with no
* remedy; at replay, the one-line footnote
log_has "`nlog'" "note: `nm' of `R' replications failed and were dropped: in them, raking did not reach its targets"
assert r(count) == 1
log_has "`nlog'" "sparse levels"
assert r(found) == 0
log_has "`nlog'" "raise iterate()"
assert r(found) == 0
log_has "`nlog'" "more than 5%"
assert r(found) == 0
log_has "`nlog'" "replications failed and were dropped)"
assert r(found) == 0
log_has "`nlog'" "(raking failed or missed its targets)"
assert r(found) == 0
set linesize 120
log using "`nlog'", text replace name(f4note)
seldecomp
log close f4note
set linesize `ls'
log_has "`nlog'" "(`nm' of `R' replications failed and were dropped; bootstrap results may be incorrect)"
assert r(count) == 1
log_has "`nlog'" "note: `nm' of"
assert r(found) == 0
* the saving() file keeps a row of missing values for each failed replicate, and its
* characteristics give the number of completed ones
preserve
use f4reps.dta, clear
assert _N == `R'
unab svars : _all
local nv : word count `svars'
egen int nmiss = rowmiss(`svars')
assert nmiss == 0 | nmiss == `nv'
count if nmiss == `nv'
assert r(N) == `nm'
assert "`: char _dta[N_reps]'" == "`nr'"
drop nmiss
mkmat `svars', matrix(SV)
restore
* e(bsdraws) and e(bsdraws_pct) hold the same rows as the saving() file, missing ones included
assert "`: word 7 of `svars''" == "pct_explained" & `nv' == 7
matrix SV6 = SV[1..., 1..6]
matrix SV7 = SV[1..., 7]
mat_same DR SV6 0
mat_same DRP SV7 0
* a replicate fails exactly when a plain run on its weights stops with 430 (a missed target);
* the completed rows equal the plain runs
local nfail 0
forvalues r = 1/`R' {
    capture seldecomp y a b [pw=ff_`r'], by(group)
    if _rc {
        assert _rc == 430
        local ++nfail
        assert missing(SV[`r', 1]) & missing(SV[`r', 5])
    }
    else row_vs_e SV `r' 1e-10
}
assert `nfail' == `nm'
drop ff_*
* -bstat using- reads the same counts and, from the completed rows, the same results
bstat using f4reps.dta, notable noheader
assert e(N_reps) == `nr'
assert e(N_misreps) == `nm'
matrix bse = e(se)
matrix bbb = e(b_bs)
local c = colnumb(bse, "overall:explained")
assert reldif(el(bse, 1, `c'), el(SE, 1, 5)) < 1e-8
assert reldif(el(bbb, 1, `c'), el(BB, 1, 5)) < 1e-8
erase f4reps.dta
* the BCa jackknife deletions follow the same rule: cluster 11 holds 10 of the Reference Group's 30 (0, 0)
* observations and nothing else, so deleting it leaves 20 against 24 (1, 1), margins that
* cannot be met: that deletion fails, and only that one (the plain runs without each cluster);
* the acceleration uses the others, without a note
gen int cl = cond(group == 1 & k <= 10, 11, 1 + mod(k, 10))
set linesize 120
log using "`nlog'", text replace name(f4note)
noisily seldecomp y a b, by(group) bayesboot(reps(`R') rseed(1) nodots bca cluster(cl))
log close f4note
set linesize `ls'
assert e(N_reps) + e(N_misreps) == `R' & e(N_iterlimit) == 0
log_has "`nlog'" "jackknife deletion"
assert r(found) == 0
forvalues j = 1/11 {
    capture seldecomp y a b if cl != `j', by(group)
    assert _rc == cond(`j' == 11, 430, 0)
}
* fewer failures (1 of 20 with these draws, once 4 of the Reference Group's (1, 1) observations move to
* the other cells): the same note, with the singular; checked in whichever branch the draws
* give
replace a = !mod(k, 2) if group == 1 & inrange(k, 51, 54)
replace b = mod(k, 2) if group == 1 & inrange(k, 51, 54)
set linesize 120
log using "`nlog'", text replace name(f4note)
noisily seldecomp y a b, by(group) bayesboot(reps(`R') rseed(7) nodots)
log close f4note
set linesize `ls'
local nm = e(N_misreps)
assert e(N_reps) + `nm' == `R' & `nm' < `R' - 1
local ws = cond(`nm' == 1, "was", "were")
local wt = cond(`nm' == 1, "it", "them")
log_has "`nlog'" "note: `nm' of `R' replications failed and `ws' dropped: in `wt', raking "
assert r(found) == (`nm' > 0)
log_has "`nlog'" "sparse levels"
assert r(found) == 0
use "`f4main'", clear

di as txt _n "---- 9b2. replicates kept at the iterate() limit"
* a replicate whose raking stops at the iterate() limit before the weights converge but
* meets every target is kept, as a point estimate in that state is (a note, e(converged) = 0,
* rc 0), and counted in e(N_iterlimit). With iterate(11) tolerance(1e-10) the
* point estimate converges at the 11th cycle and some replicates need more (8 of 20 with
* these draws; a range is asserted); their margins are met to about 1e-10
local ipf iterate(11) tolerance(1e-10)
set linesize 120
log using "`nlog'", text replace name(f4note)
noisily seldecomp y x1 x2 x4 [pw=wt], by(group) `ipf' ///
    bayesboot(reps(`R') rseed(7) nodots saving(f4reps.dta, replace) generate(ff))
log close f4note
set linesize `ls'
local ni = e(N_iterlimit)
assert e(converged) == 1
assert e(N_reps) == `R' & e(N_misreps) == 0 & "`e(missing)'" == ""
assert inrange(`ni', 2, `R' - 2)
log_has "`nlog'" "note: `ni' of `R' replications stopped at the iterate() limit with the raking targets met and were kept; raise iterate()"
assert r(count) == 1
log_has "`nlog'" "replications stopped at iterate()"
assert r(found) == 0
log_has "`nlog'" "failed"
assert r(found) == 0
* at replay, the one-line footnote in place of the note
set linesize 120
log using "`nlog'", text replace name(f4note)
seldecomp
log close f4note
set linesize `ls'
log_has "`nlog'" "(`ni' of `R' replications stopped before the weights converged)"
assert r(count) == 1
log_has "`nlog'" "note: `ni' of"
assert r(found) == 0
preserve
use f4reps.dta, clear
assert _N == `R'
unab svars : _all
mkmat `svars', matrix(SV)
restore
erase f4reps.dta
* they are the replicates whose plain run stops at the limit with its targets met (rc 0,
* e(converged) = 0, e(badcontrols) = 0); every row equals its plain run
local nit 0
forvalues r = 1/`R' {
    qui seldecomp y x1 x2 x4 [pw=ff_`r'], by(group) `ipf'
    if e(converged) == 0 {
        assert e(badcontrols) == 0
        local ++nit
    }
    row_vs_e SV `r' 1e-10
}
assert `nit' == `ni'
drop ff_*
* the BCa jackknife deletions follow the same rule. With iterate(9) tolerance(1e-14) every
* raking stops at the limit with its targets met (to about 1e-12): the point estimate (its
* note, e(converged) = 0, shown again at replay), all 10 replicates and all 150 deletions are
* kept; the point estimate's note counts the replications (one message on the limit, at
* estimation and at replay), nothing on the deletions; the same with sequentially
* fixothermargins and reverse
set linesize 120
log using "`nlog'", text replace name(f4note)
noisily seldecomp y x1 x2 [pw=wt] in 1/150, by(group) iterate(9) tolerance(1e-14) ///
    bayesboot(reps(10) rseed(3) bca nodots)
seldecomp
log close f4note
set linesize `ls'
assert e(converged) == 0 & e(badcontrols) == 0
assert e(N_reps) == 10 & e(N_misreps) == 0 & e(N_iterlimit) == 10
log_has "`nlog'" "jackknife deletion"
assert r(found) == 0
log_has "`nlog'" "note: 10 of 10 replications stopped"
assert r(found) == 0
log_has "`nlog'" "note: raking stopped at the iterate() limit, also in 10 of 10 replications; raise iterate()"
assert r(count) == 2
log_has "`nlog'" "before the weights converged"
assert r(found) == 0
log_has "`nlog'" "failed"
assert r(found) == 0
log_has "`nlog'" "(10 of 10 replications"
assert r(found) == 0
seldecomp y (x1) (x2) [pw=wt] in 1/150, by(group) sequentially fixothermargins reverse ///
    iterate(9) tolerance(1e-14) bayesboot(reps(10) rseed(3) bca nodots)
assert e(N_reps) == 10 & e(N_misreps) == 0 & e(N_iterlimit) == 10

di as txt _n "---- 9c. failed jackknife deletions"
* bca with cluster(): a level of x found in the Focal Group only in cluster 1 (and in the Reference Group in
* clusters 5 and 6) is absent from the Focal Group when cluster 1 is deleted, so that deletion
* fails (Reference Group has a level Focal Group lacks) while no replicate does; a second such level,
* in cluster 2, makes two. The acceleration still uses the other deletions, with a caution counting failures
preserve
clear
set obs 400
gen int cl = ceil(_n / 40)
gen byte g = mod(_n, 2)
gen byte x = 1 + mod(floor(_n / 2), 2)
replace x = 3 if (g == 0 & cl == 1) | (g == 1 & inlist(cl, 5, 6))
gen double yy = x + mod(_n * 7, 11) / 10
foreach k in 1 2 {
    if `k' == 2 replace x = 4 if (g == 0 & cl == 2) | (g == 1 & inlist(cl, 7, 8))
    log using "`nlog'", text replace name(f4note)
    noisily seldecomp yy x, by(g) bayesboot(reps(10) rseed(1) bca cluster(cl) nodots)
    log close f4note
    assert e(N_reps) == 10 & e(N_misreps) == 0
    assert matmissing(e(accel)) == 0
    assert e(N_jkfail) == `k'
    log_has "`nlog'" "Caution: `k' jackknife replicates failed; BCa intervals may be incorrect."
    assert r(count) == 1
}
restore


* ==============================================================================
* 10. e(bsdraws): a replay at another level recomputes every interval exactly as the
*     estimation at that level does; e() is left unchanged, r() holds them
* ==============================================================================
di as txt _n "---- 10. replay at another level"
local ivs ci_normal ci_percentile ci_bc ci_bca pct_explained_ci
* one stage, bca and ties, on a 100-observation subsample (one jackknife deletion each)
seldecomp y x1 x2 [pw=wt] in 1/100, by(group) bayesboot(reps(20) rseed(21) bca ties nodots)
matrix D = e(bsdraws)
matrix DP = e(bsdraws_pct)
assert rowsof(D) == 20 & colsof(D) == 6 & rowsof(DP) == 20 & colsof(DP) == 1
local dn : colfullnames D
assert "`dn'" == "`e(bsnames)'"
local dn : colnames DP
assert "`dn'" == "pct_explained"
* the columns are the replicates behind e(b_bs) and e(se)
mata: st_matrix("MD", mean(st_matrix("D")))
mata: st_matrix("SD", sqrt(diagonal(variance(st_matrix("D"))))')
mat_same MD e(b_bs) 1e-12
mat_same SD e(se) 1e-10
* a replay at the estimation level returns e()'s intervals
seldecomp
assert r(level) == 95
foreach m of local ivs {
    mat_same r(`m') e(`m') 0
}
seldecomp, level(90) citable
assert r(level) == 90 & e(level) == 95
foreach m of local ivs {
    matrix R90_`m' = r(`m')
}
estimates store bsd
* the fresh estimation at level(90) with the same seed
seldecomp y x1 x2 [pw=wt] in 1/100, by(group) bayesboot(reps(20) rseed(21) bca ties nodots level(90))
foreach m of local ivs {
    mat_same R90_`m' e(`m') 1e-12
}
local n1 : colfullnames R90_ci_bca
local n2 : colfullnames e(ci_bca)
local r1 : rownames R90_ci_bca
assert "`n1'" == "`n2'" & "`r1'" == "ll ul"
* estimates store / restore keep the replicates: the same intervals at replay
seldecomp y x1 [pw=wt], by(group)
estimates restore bsd
mat_same e(bsdraws) D 0
mat_same e(bsdraws_pct) DP 0
seldecomp, level(90)
foreach m of local ivs {
    mat_same r(`m') R90_`m' 0
}
estimates drop bsd
* every stage bootstrapped, clusters and bca (3 clusters: 3 jackknife deletions): the replay
* at level(80) against a fresh run
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially fixothermargins ///
    bayesboot(reps(20) rseed(22) cluster(strata) bca nodots)
assert colsof(e(bsdraws)) == 18 & rowsof(e(bsdraws)) == 20
seldecomp, level(80) citable
foreach m of local ivs {
    matrix R80_`m' = r(`m')
}
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially fixothermargins ///
    bayesboot(reps(20) rseed(22) cluster(strata) bca nodots level(80))
foreach m of local ivs {
    mat_same R80_`m' e(`m') 1e-12
}


* ==============================================================================
* bootstrap edge cases: a difference of 0, failed replicates in the saving() file
* ==============================================================================
di as txt _n "---- bootstrap edge cases"
preserve

* an observed difference of exactly 0: the explained share is undefined, so it has no
* bootstrap SE or CI, and the saving() file (without pct_explained) is readable by bstat
clear
set obs 40
gen byte g = _n > 20
gen byte x = mod(_n, 2)
gen byte yb = mod(_n, 4) < 2
capture erase r2zero.dta
seldecomp yb x, by(g) bayesboot(reps(20) rseed(1) nodots saving(r2zero, replace))
assert e(difference) == 0 & missing(e(pct_explained))
assert missing(e(pct_explained_se))
matrix PC = e(pct_explained_ci)
assert missing(PC[1, 1]) & missing(PC[1, 2])
* nor replicates of it in e(); the other replicates are kept, and a replay at another level
* recomputes their intervals only
confirm matrix e(bsdraws)
capture confirm matrix e(bsdraws_pct)
assert _rc != 0
seldecomp, level(90)
assert r(level) == 90
capture confirm matrix r(pct_explained_ci)
assert _rc != 0
bstat using r2zero, notable noheader
assert e(N_reps) == 20 & colsof(e(b)) == 6
erase r2zero.dta

* the saving() file has one row per requested replicate, missing for a failed one (the
* margins of 9b that some replicates cannot meet: 9 of 30 with these draws)
clear
set obs 400
gen byte group = _n > 200
gen int k = cond(group == 0, _n, _n - 200)
gen byte a = cond(group == 0, k > 130, cond(k <= 30, 0, cond(k <= 54, 1, !mod(k, 2))))
gen byte b = cond(group == 0, inrange(k, 61, 130), cond(k <= 30, 0, cond(k <= 54, 1, mod(k, 2))))
gen double y = a + 2 * b + mod(_n * 7, 11) / 10
capture erase r2fail.dta
seldecomp y a b, by(group) bayesboot(reps(30) rseed(6) nodots saving(r2fail, replace))
local nm = e(N_misreps)
assert `nm' > 0 & `nm' < 30 & e(N_reps) + `nm' == 30
use r2fail, clear
assert _N == 30
count if missing(explained)
assert r(N) == `nm'
erase r2fail.dta
restore

* ==============================================================================
* the explained share with heavy tails; fewer than 10 clusters
* ==============================================================================
preserve
tempfile s1log
capture program drop s1_note
program define s1_note, rclass
    // r(found): the number of the e(note#) that contains the text, 0 if none
    args text
    local f 0
    forvalues m = 1/`e(N_notes)' {
        if strpos(`"`e(note`m')'"', `"`text'"') local f `m'
    }
    return scalar found = `f'
end

* at level(90) the CI of the difference excludes 0 but the share's
* std. err. exceeds 100: the std. err. and the percentile interval, which the heavy tails
* hardly move, without comment; the same reporting at 95% when the difference's CI includes
* zero. Neither replay changes the stored estimates or uncertainty calculations.
sysuse nlsw88, clear
seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ///
    bayesboot(reps(100) rseed(1) nodots)
local pse = e(pct_explained_se)
matrix PCI0 = e(pct_explained_ci)
assert `pse' > 100
local sse = strtrim(string(`pse', "%9.1f"))
local sh = strtrim(string(e(pct_explained), "%9.1f"))
log using `"`s1log'"', text replace name(s1u10)
seldecomp, level(90)
matrix P90 = r(pct_explained_ci)
log close s1u10
local pl = strtrim(string(P90[1, 1], "%9.1f"))
local pu = strtrim(string(P90[1, 2], "%9.1f"))
log_has `"`s1log'"' "Explained share of the difference: `sh'%, bootstrap std. err. `sse' 90% percentile CI: `pl'% to `pu'%"
assert r(found) == 1
log_has `"`s1log'"' "unstable"
assert r(found) == 0
log_has `"`s1log'"' "heavy-tailed"
assert r(found) == 0
assert e(pct_explained_se) == `pse' & mreldif(e(pct_explained_ci), PCI0) == 0
log using `"`s1log'"', text replace name(s1u10)
seldecomp
log close s1u10
local pl = strtrim(string(PCI0[1, 1], "%9.1f"))
local pu = strtrim(string(PCI0[1, 2], "%9.1f"))
log_has `"`s1log'"' "Explained share of the difference: `sh'%, bootstrap std. err. `sse' 95% percentile CI: `pl'% to `pu'%"
assert r(found) == 1
log_has `"`s1log'"' "unstable"
assert r(found) == 0
log_has `"`s1log'"' "heavy-tailed"
assert r(found) == 0

* the data of the checks with few clusters
clear
set obs 300
gen byte group = mod(_n, 2)
gen byte x = 1 + mod(floor(_n / 3), 2)
gen double y = x + group + mod(_n * 7, 11) / 10

* fewer than 10 clusters: the header gives their number, and no note cautions at estimation
* or replay (the help's Remarks give the caution)
gen int c8 = ceil(_n / 37.5)
gen int c10 = ceil(_n / 30)
log using `"`s1log'"', text replace name(s1u28)
seldecomp y x, by(group) bayesboot(reps(3) nodots rseed(1) cluster(c8))
seldecomp
log close s1u28
assert e(N_clust) == 8
s1_note "clusters"
assert r(found) == 0
log_has `"`s1log'"' "under-cover"
assert r(found) == 0
log_has `"`s1log'"' "(Replications based on 8 clusters in c8)"
assert r(count) == 2
seldecomp y x, by(group) bayesboot(reps(3) nodots rseed(1) cluster(c10))
assert e(N_clust) == 10
s1_note "clusters: intervals under-cover"
assert r(found) == 0
restore

* ==============================================================================
* estat bootstrap after bayesboot(): the intervals and digits of citable, Stata's layout, the
* options, the errors, the level of the last display, estimates store/restore, and the prefixes
* ==============================================================================
di as txt _n "---- estat bootstrap"
preserve
tempfile e8c e8e
capture mata: mata drop e8_toks() e8_width() e8_ends()
mata:
string scalar e8_toks(string scalar fn)
{
    // the rows of the bootstrap-interval tables in the log (ending in (N), (P), (BC) or
    // (BCa)), each as the numbers and key after its "|", joined by " / "
    string colvector l
    string scalar out
    real scalar i
    l = cat(fn)
    l = select(l, ustrregexm(l, "[(](N|P|BC|BCa)[)]$"))
    out = ""
    for (i = 1; i <= rows(l); i++) {
        out = out + " / " + strtrim(stritrim(substr(l[i], strpos(l[i], "|") + 1, .)))
    }
    return(out)
}
string scalar e8_ends(string scalar fn, string scalar txt)
{
    // the first line of the log that contains txt: the columns where its words (runs of
    // nonblanks) end, and its length, as "e1 e2 ... / length"
    string colvector l
    string scalar s, out
    real scalar i
    l = cat(fn)
    l = select(l, strpos(l, txt) :> 0)
    if (!rows(l)) return("")
    s = l[1]
    out = ""
    for (i = 1; i <= strlen(s); i++) {
        if (substr(s, i, 1) != " " & (i == strlen(s) | substr(s, i + 1, 1) == " ")) {
            out = out + (out == "" ? "" : " ") + strofreal(i)
        }
    }
    return(out + " / " + strofreal(udstrlen(s)))
}
real scalar e8_width(string scalar fn)
{
    // the widest line of the log, command echoes and the log's own header lines aside
    string colvector l
    l = cat(fn)
    l = select(l, substr(l, 1, 2) :!= ". " :& substr(l, 1, 2) :!= "> ")
    l = select(l, !ustrregexm(l, "^ *(name|log|log type|opened on|closed on):") :&
        !ustrregexm(l, "^-+$"))
    return(rows(l) ? max(udstrlen(l)) : 0)
}
end
capture program drop e8_same
program define e8_same
    // e8_same "file1" "file2"  --  the interval rows of the two logs are identical, and there
    // is at least one
    args f1 f2
    mata: st_local("t1", e8_toks(st_local("f1")))
    mata: st_local("t2", e8_toks(st_local("f2")))
    if `"`t1'"' == "" | `"`t1'"' != `"`t2'"' {
        di as err "e8_same: the interval rows differ"
        di as txt `"  citable: `t1'"'
        di as txt `"  estat:   `t2'"'
        exit 9
    }
end
capture program drop e8_log
program define e8_log
    // e8_log "file" command  --  runs the command with its output logged to the file, then
    // again quietly, which leaves its r() (log close replaces r())
    gettoken fn 0 : 0
    log using `"`fn'"', text replace name(e8)
    `0'
    log close e8
    quietly `0'
end
capture program drop e8_count
program define e8_count, rclass
    // e8_count "file" "text"  --  r(n) = the number of lines of the log that contain the text
    // as it stands (blanks included)
    args fn txt
    local n 0
    mata: st_local("n", strofreal(sum(strpos(cat(st_local("fn")), st_local("txt")) :> 0)))
    return scalar n = `n'
end

* (a) bayesboot() with bca and clusters: estat bootstrap, all equals citable to the digit
* (the rows after "|", both at %9.0g) and e() exactly; normal, so that citable shows the normal
* intervals as estat does (with clusters the default is t). Explicit types are in
* estat's canonical order for this layout comparison; typed-order behavior is tested below.
seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(30) rseed(12345) normal percentile bc bca cluster(cluster) nodots)
assert "`e(estat_cmd)'" == "seldecomp_estat"
local kb : word count `e(bsnames)'
e8_log `"`e8c'"' seldecomp, citable
tempname Rn Rp Rb Ra Rs E R
matrix `Rn' = r(ci_normal)
matrix `Ra' = r(ci_bca)
e8_log `"`e8e'"' estat bootstrap, all
assert r(level) == e(level) & r(level) == 95
foreach m in ci_normal ci_percentile ci_bc ci_bca {
    matrix `E' = e(`m')
    matrix `R' = r(`m')
    mat_same `E' `R' 0
}
matrix `R' = r(ci_normal)
mat_same `Rn' `R' 0
matrix `R' = r(ci_bca)
mat_same `Ra' `R' 0
e8_same `"`e8c'"' `"`e8e'"'
* Stata's layout (_bs_display on Stata 19): the header, the
* legend (the definitions), the column headings, a row, the key; for six lines the columns
* where their words end, and their length
foreach t in "Raking decomposition (Pylypchuk-Selden)" "Number of obs     =      2,000"   ///
    "Replications      =         30" "Difference = focal - reference. Explained = focal - raked focal."                    ///
    "Estimate"                                                                               ///
    "(Replications based on 60 clusters in cluster)" "Key:   N: Normal"                    ///
    "       P: Percentile" "      BC: Bias-corrected" "     BCa: Bias-corrected and accelerated" {
    e8_count `"`e8e'"' "`t'"
    assert r(n) == 1
}
local lay1 `""Number of obs" "6 20 39 54 57 61 67 78 / 78""'
local lay2 `""Replications      =" "60 67 78 / 78""'
local lay3 `""(Replications based on" "45 51 54 57 66 69 78 / 78""'
local lay4 `""Bootstrap" "18 50 / 50""'
local lay5 `""Estimate" "16 18 29 40 45 50 56 62 72 / 72""'
local lay6 `""      Focal mean |" "11 16 18 29 40 50 61 72 78 / 78""'
forvalues i = 1/6 {
    gettoken t e : lay`i'
    gettoken e : e
    mata: st_local("got", e8_ends(st_local("e8e"), st_local("t")))
    if `"`got'"' != `"`e'"' {
        di as err `"layout of the line with `t': `got' (expected `e')"'
        exit 9
    }
}
foreach k in "   (N)" "   (P)" "  (BC)" " (BCa)" {
    e8_count `"`e8e'"' "`k'"
    assert r(n) == `kb'
}
* one equation: no equation row; no command line without verbose; 78 columns
foreach t in "overall" "Command:" {
    e8_count `"`e8e'"' "`t'"
    assert r(n) == 0
}
mata: st_local("w", strofreal(e8_width(st_local("e8e"))))
assert `w' == 78

* (b) the options: bc by default; normal, percentile, bc, bca alone; nolegend, verbose
e8_log `"`e8e'"' estat bootstrap
e8_count `"`e8e'"' "  (BC)"
assert r(n) == `kb'
foreach k in "(N)" "(P)" "(BCa)" {
    e8_count `"`e8e'"' "`k'"
    assert r(n) == 0
}
e8_count `"`e8e'"' "Key: BC: Bias-corrected"
assert r(n) == 1
foreach o in normal percentile bc bca {
    log using `"`e8e'"', text replace name(e8)
    estat bootstrap, `o'
    log close e8
    foreach k2 in normal percentile bc bca {
        local k = cond("`k2'" == "normal", "(N)", cond("`k2'" == "percentile", "(P)", ///
            cond("`k2'" == "bc", "(BC)", "(BCa)")))
        e8_count `"`e8e'"' "`k'"
        assert r(n) == cond("`k2'" == "`o'", `kb', 0)
    }
}
e8_log `"`e8e'"' estat bootstrap, nolegend
e8_count `"`e8e'"' "Difference = focal - reference."
assert r(n) == 0
e8_count `"`e8e'"' "Number of obs"
assert r(n) == 1
e8_log `"`e8e'"' estat bootstrap, verbose
e8_count `"`e8e'"' "      Command: seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(30)"
assert r(n) == 1
e8_count `"`e8e'"' "Difference = focal - reference. Explained = focal - raked focal."
assert r(n) == 1
* the subcommand abbreviated to 4 letters, not 3
estat boot
assert r(level) == 95
rcof "noi estat boo" == 321

* (c) the errors: level(), an option estat bootstrap does not take (noheader among them), bca
* without bca; the other subcommands (estat vce is Stata's own)
log using `"`e8e'"', text replace name(e8)
rcof "noi estat bootstrap, level(90)" == 198
rcof "noi estat bootstrap, title(x)" == 198
rcof "noi estat bootstrap, noheader" == 198
rcof "noi estat summarize" == 321
rcof "noi estat ic" == 321
rcof "noi estat nonsense" == 321
log close e8
foreach t in "option level() not allowed: estat bootstrap uses the level of the last" ///
    "display; type seldecomp, level(#) first" "estat summarize not valid"                 ///
    "estat ic not valid" "option noheader not allowed" {
    e8_count `"`e8e'"' "`t'"
    assert r(n) == 1
}
estat vce
seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(10) rseed(1) nodots)
e8_log `"`e8e'"' rcof "noi estat bootstrap, bca" == 198
e8_count `"`e8e'"' "option bca requires results from bayesboot() with bca"
assert r(n) == 1
estat bootstrap, all
assert !missing(el(r(ci_bc), 1, 1))
capture confirm matrix r(ci_bca)
assert _rc
seldecomp y x1 x2 [pw=wt], by(group)
assert "`e(estat_cmd)'" == "seldecomp_estat"
e8_log `"`e8e'"' rcof "noi estat bootstrap" == 321
e8_count `"`e8e'"' "estat bootstrap requires results from bayesboot() or the bootstrap prefix"
assert r(n) == 1

* (d) the level of the last display: a replay at level(90) gives estat its intervals, e()
* keeps the estimation's; a replay at the estimation level brings estat back to it
seldecomp y (a: x1) (b: x2 x4) [pw=wt], by(group) sequentially fixothermargins outcomes(ybin) ///
    bayesboot(reps(20) rseed(3) ties nodots)
local kb : word count `e(bsnames)'
tempname B0
matrix `B0' = e(ci_bc)
assert e(level_replay) >= .
e8_log `"`e8c'"' seldecomp, level(90) citable
matrix `Rn' = r(ci_normal)
matrix `Rp' = r(ci_percentile)
matrix `Rb' = r(ci_bc)
matrix `Rs' = r(pct_explained_ci)
assert e(level) == 95 & e(level_replay) == 90
matrix `E' = e(ci_bc)
mat_same `B0' `E' 0
e8_log `"`e8e'"' estat bootstrap, all
assert r(level) == 90
matrix `R' = r(ci_normal)
mat_same `Rn' `R' 0
matrix `R' = r(ci_percentile)
mat_same `Rp' `R' 0
matrix `R' = r(ci_bc)
mat_same `Rb' `R' 0
matrix `R' = r(pct_explained_ci)
mat_same `Rs' `R' 0
e8_same `"`e8c'"' `"`e8e'"'
* several equations: a row for each; several outcomes: no depvar in the corner; ties
foreach t in "Outcome          | ybin" "[90% conf. interval]"                       ///
    "     BC: Bias-corrected (adjusted for ties)" "Difference = focal - reference. Explained = focal - raked focal." {
    e8_count `"`e8e'"' "`t'"
    assert r(n) == 1
}
local t "Estimate"
local e "18 29 40 45 50 56 62 72 / 72"
mata: st_local("got", e8_ends(st_local("e8e"), st_local("t")))
if `"`got'"' != `"`e'"' {
    di as err `"layout of the line with `t': `got' (expected `e')"'
    exit 9
}
foreach k in "   (N)" "   (P)" "  (BC)" {
    e8_count `"`e8e'"' "`k'"
    assert r(n) == `kb'
}
mata: st_local("w", strofreal(e8_width(st_local("e8e"))))
assert `w' <= 78
* a replay at another level again, then at the estimation level
quietly seldecomp, level(80)
assert e(level_replay) == 80
estat bootstrap
assert r(level) == 80
* the graph at replay is no display of the intervals: the record stays
quietly seldecomp, makeplot
assert e(level_replay) == 80
graph drop seldecomp
quietly seldecomp
assert e(level_replay) == 95
estat bootstrap, all
assert r(level) == 95
matrix `R' = r(ci_bc)
mat_same `B0' `R' 0
* a level with decimals: _bs_display's heading
quietly seldecomp, level(97.5)
e8_log `"`e8e'"' estat bootstrap, nolegend
e8_count `"`e8e'"' "std. err.   [ 97.5% conf. int.]"
assert r(n) == 1

* (e) after estimates store and restore: the results as stored (stored before a replay: the
* estimation level; stored after one: its level), and estat equals citable after a replay
seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(20) rseed(7) nodots)
estimates store e8s
tempname S0
matrix `S0' = e(ci_percentile)
quietly seldecomp, level(90)
estimates store e8t
seldecomp y x1 [pw=wt], by(group)
estimates restore e8s
estat bootstrap, all
assert r(level) == 95
matrix `R' = r(ci_percentile)
mat_same `S0' `R' 0
e8_log `"`e8c'"' seldecomp, citable
e8_log `"`e8e'"' estat bootstrap, all
e8_same `"`e8c'"' `"`e8e'"'
estimates restore e8t
estat bootstrap, all
assert r(level) == 90
matrix `Rp' = r(ci_percentile)
e8_log `"`e8c'"' seldecomp, level(90) citable
matrix `R' = r(ci_percentile)
mat_same `Rp' `R' 0
e8_log `"`e8e'"' estat bootstrap, all
e8_same `"`e8c'"' `"`e8e'"'
* a replay at the estimation level brings the restored results back to it
quietly seldecomp
assert e(level_replay) == 95
estat bootstrap, percentile
assert r(level) == 95
matrix `R' = r(ci_percentile)
mat_same `S0' `R' 0
estimates drop e8s e8t
* estimates replay: results without a stored name are
* replayed in place and keep the record; with a stored name (estimates replay name, or
* estimates replay alone after estimates store) the stored set is replayed and gets the
* record, which estimates restore from other results brings back, while the active
* results keep their own (also through an estimates restore of their own name, a no-op);
* a record changes e(), and Stata then drops the stored name
* the named form: the stored set gets the record, the active results keep theirs
quietly seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(20) rseed(7) nodots)
estimates store e8s
quietly estimates replay e8s, level(80)
assert e(level_replay) >= .
quietly estimates restore e8s
assert e(level_replay) >= .
quietly estat bootstrap
assert r(level) == 95
quietly seldecomp y x1 [pw=wt], by(group)
quietly estimates restore e8s
assert e(level_replay) == 80
quietly estat bootstrap
assert r(level) == 80
estimates drop e8s
* the bare form of results never stored: in place; of results just stored: the stored set
quietly seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(20) rseed(7) nodots)
quietly estimates replay, level(90)
assert e(level_replay) == 90
quietly estat bootstrap
assert r(level) == 90
estimates store e8s
quietly estimates replay, level(85)
assert e(level_replay) == 90
quietly estat bootstrap
assert r(level) == 90
quietly seldecomp y x1 [pw=wt], by(group)
quietly estimates restore e8s
assert e(level_replay) == 85
quietly estat bootstrap
assert r(level) == 85
* restored from a set a replay changed: no stored name, so in place again
assert `"`e(_estimates_name)'"' == ""
quietly estimates replay, level(75)
assert e(level_replay) == 75
quietly estat bootstrap
assert r(level) == 75
quietly seldecomp y x1 [pw=wt], by(group)
quietly estimates restore e8s
assert e(level_replay) == 85
estimates drop e8s

* (f) the prefixes: after bootstrap: Stata's own estat bootstrap; after
* jackknife: the message; svy replaces e(estat_cmd) (test_svy.do: estat effects)
bootstrap, reps(5) seed(1) nodots: seldecomp y x1 x2, by(group)
assert "`e(estat_cmd)'" == "seldecomp_estat" & "`e(prefix)'" == "bootstrap"
e8_log `"`e8e'"' estat bootstrap, all
e8_count `"`e8e'"' "  (BC)"
assert r(n) == 6
if c(stata_version) >= 19 {
    bayesboot, reps(5) rseed(1) nodots: seldecomp y x1 x2, by(group)
    assert "`e(estat_cmd)'" == "seldecomp_estat" & "`e(prefix)'" == "bayesboot"
    e8_log `"`e8e'"' estat bootstrap, all
    e8_count `"`e8e'"' "  (BC)"
    assert r(n) == 6
}
jackknife, cluster(cluster) nodots: seldecomp y x1 x2, by(group)
rcof "noi estat bootstrap" == 321
estat vce

* (g) with chain and with detail: every
* element of e(bsnames), the intervals of e(), and the digits of the citable table. Under
* chain a value is shown as 0 below its pair's rounding threshold (the total's: the largest
* of the pairs'), as in the citable table: here the first period's mean is about 1e6 and the
* others' about 1e-6, so pair 1's threshold would show pair 3's std. errors as 0
clear
set seed 20260927
set obs 800
gen int period = 1995 + 2 * ceil(_n / 200)
gen double wt = 50 + 100 * runiform()
gen byte a = 1 + (runiform() < .3 + .1 * (period - 1997) / 2)
gen byte b = 1 + (runiform() < .4) + (runiform() < .2 + .05 * (period - 1997) / 2)
gen double y0 = 1 + .5 * a - .3 * b + .2 * rnormal()
gen double y = y0 * cond(period == 1997, 1e6, 1e-6)
seldecomp y (ea: a) (eb: b) [pw=wt], by(period) chain sequentially fixothermargins ///
    bayesboot(reps(20) rseed(3) nodots)
assert _se[p3:difference] > 0 & _se[p3:difference] < 1e-13 * abs(_b[p1:focal_group])
local kb : word count `e(bsnames)'
e8_log `"`e8c'"' seldecomp, citable
e8_log `"`e8e'"' estat bootstrap, all
assert r(level) == 95
foreach m in ci_normal ci_percentile ci_bc {
    matrix `E' = e(`m')
    matrix `R' = r(`m')
    mat_same `E' `R' 0
}
e8_same `"`e8c'"' `"`e8e'"'
e8_count `"`e8e'"' "  (BC)"
assert r(n) == `kb'
* All pair and total panels remain, with descriptive stage equation headings.
foreach q in p1 p2 p3 total {
    local t = substr("`q'" + 17 * " ", 1, 17) + "|"
    e8_count `"`e8e'"' "`t'"
    assert r(n) == 1
}
e8_count `"`e8e'"' "By stage         | Cumulative unexplained"
assert r(n) == 3
e8_count `"`e8e'"' "By stage         | Marginal explained"
assert r(n) == 4
* a replay at another level: estat follows it, with the citable table's digits
quietly seldecomp, level(90)
assert e(level_replay) == 90
e8_log `"`e8c'"' seldecomp, level(90) citable
e8_log `"`e8e'"' estat bootstrap, all
assert r(level) == 90
e8_same `"`e8c'"' `"`e8e'"'
* detail: the terms' contributions (equation detail) are rows of the table too
seldecomp y0 (ab: a b) [pw=wt] if period <= 1999, by(period) sequentially fixothermargins detail ///
    bayesboot(reps(20) rseed(3) nodots)
local kb : word count `e(bsnames)'
assert strpos("`e(bsnames)'", "detail:a detail:b")
e8_log `"`e8c'"' seldecomp, citable
e8_log `"`e8e'"' estat bootstrap, all
foreach m in ci_normal ci_percentile ci_bc {
    matrix `E' = e(`m')
    matrix `R' = r(`m')
    mat_same `E' `R' 0
}
e8_same `"`e8c'"' `"`e8e'"'
e8_count `"`e8e'"' "  (BC)"
assert r(n) == `kb'
e8_count `"`e8e'"' "By term          |"
assert r(n) == 1
mata: mata drop e8_toks() e8_width() e8_ends()
restore


* ==============================================================================
* 11. the interval types: bayesboot(normal tinterval percentile bc bca citable)
* ==============================================================================
* one type: the coefficient table's interval column, headed by it, and no interval table; two
* or more: the first in the table, all in the interval table in explicitly requested order;
* citable: percentile and bc with the default of normal and t (t with cluster(), its degrees
* of freedom the clusters less 1; else normal, the t intervals' those of the observations
* less 1) unless one of them is requested; e(citypes), e(df_bs), e(ci_t), and e(df_r) when t
* is the family of the tests
di as txt _n "---- 11. interval types"
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
recast double y x3 wt
tempfile ilg
local ilf `"`ilg'.log"'
capture log close ilog
capture program drop ci_case
program define ci_case, rclass
    // ci_case logfile "bayesboot() suboptions" ["replay options"]: estimate (and replay) with
    // the output in the log; r(nmain) the coefficient tables' interval headings found, r(ntab)
    // the interval tables
    args lf bopts ropts
    quietly log using `"`lf'"', text replace name(ilog)
    seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(20) rseed(4) nodots `bopts')
    if `"`ropts'"' != "" seldecomp, `ropts'
    quietly log close ilog
end
* the default: normal, no interval table; t intervals with N - 1 degrees of freedom stored
ci_case `"`ilf'"' ""
assert "`e(citypes)'" == "normal" & e(df_bs) == e(N) - 1 & e(df_r) >= .
log_has `"`ilf'"' "Normal-based"
assert r(found) == 1
log_has `"`ilf'"' "Bayesian bootstrap confidence intervals"
assert r(found) == 0
matrix B = e(b)
matrix SE = e(se)
matrix CT = e(ci_t)
local q = invttail(e(df_bs), .025)
forvalues j = 1/6 {
    assert reldif(CT[1, `j'], B[1, `j'] - `q' * SE[1, `j']) < 1e-12
    assert reldif(CT[2, `j'], B[1, `j'] + `q' * SE[1, `j']) < 1e-12
}
* percentile alone: its intervals in the coefficient table
ci_case `"`ilf'"' "percentile"
assert "`e(citypes)'" == "percentile"
log_has `"`ilf'"' "Percentile"
assert r(found) == 1
log_has `"`ilf'"' "Normal-based"
assert r(found) == 0
log_has `"`ilf'"' "Bayesian bootstrap confidence intervals"
assert r(found) == 0
matrix P = e(ci_percentile)
local ll : di %9.0g P[1, 4]
log_has `"`ilf'"' "`=strtrim("`ll'")'"
assert r(found) == 1
* percentile and bc: percentile in the coefficient table, both in the interval table
ci_case `"`ilf'"' "percentile bc"
assert "`e(citypes)'" == "percentile bc"
log_has `"`ilf'"' "Key: P: percentile, BC: bias-corrected"
assert r(found) == 1
foreach t in "(N)" "(T)" "(BCa)" {
    log_has `"`ilf'"' "`t'"
    assert r(found) == 0
}
log_has `"`ilf'"' "(P)"
assert r(count) == 6
* citable = normal percentile bc without clusters (the same output after the command line)
ci_case `"`ilf'"' "citable"
assert "`e(citypes)'" == "normal percentile bc"
local o1 : copy local ilf
tempfile ilg2
local ilf2 `"`ilg2'.log"'
ci_case `"`ilf2'"' "normal percentile bc"
assert "`e(citypes)'" == "normal percentile bc"
mata: _a = cat(st_local("o1")); _b = cat(st_local("ilf2"))
mata: _a = select(_a, !strpos(_a, "seldecomp y x1") :& !strpos(_a, "> ") :& !strpos(_a, " on: "))
mata: _b = select(_b, !strpos(_b, "seldecomp y x1") :& !strpos(_b, "> ") :& !strpos(_b, " on: "))
mata: assert(_a == _b)
mata: mata drop _a _b
* normal percentile: the interval table without the bias-corrected rows
ci_case `"`ilf'"' "normal percentile"
log_has `"`ilf'"' "Key: N: normal, P: percentile"
assert r(found) == 1
log_has `"`ilf'"' "(BC)"
assert r(found) == 0
* bca alone: its intervals in the coefficient table; with citable a row of the table
ci_case `"`ilf'"' "bca"
assert "`e(citypes)'" == "bca"
log_has `"`ilf'"' "BCa"
assert r(found) == 1
log_has `"`ilf'"' "Bayesian bootstrap confidence intervals"
assert r(found) == 0
ci_case `"`ilf'"' "citable bca"
assert "`e(citypes)'" == "bca normal percentile bc"
log_has `"`ilf'"' "(BCa)"
assert r(count) == 6
* tinterval: t intervals and tests, N - 1 degrees of freedom, said under the table
ci_case `"`ilf'"' "tinterval"
assert "`e(citypes)'" == "t" & e(df_r) == e(N) - 1
log_has `"`ilf'"' "t-based"
assert r(found) == 1
log_has `"`ilf'"' "P>|t|"
assert r(found) == 1
log_has `"`ilf'"' "(t: `=e(N) - 1' degrees of freedom, the number of observations less 1)"
assert r(found) == 1
* normal and tinterval: both, normal first (the tests z)
ci_case `"`ilf'"' "normal tinterval"
assert "`e(citypes)'" == "normal t" & e(df_r) >= .
log_has `"`ilf'"' "Key: N: normal, T: t with `=e(N) - 1' degrees of freedom"
assert r(found) == 1
* clusters: t by default (the clusters less 1); citable adds t, not normal; normal citable
* keeps normal and adds no t
quietly log using `"`ilf'"', text replace name(ilog)
seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(20) rseed(4) nodots cluster(cluster))
quietly log close ilog
assert "`e(citypes)'" == "t" & e(df_bs) == e(N_clust) - 1 & e(df_r) == e(df_bs)
log_has `"`ilf'"' "(t: `=e(N_clust) - 1' degrees of freedom, the number of clusters less 1)"
assert r(found) == 1
matrix B = e(b)
matrix SE = e(se)
matrix CT = e(ci_t)
local q = invttail(e(N_clust) - 1, .025)
assert reldif(CT[1, 4], B[1, 4] - `q' * SE[1, 4]) < 1e-12
seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(20) rseed(4) nodots cluster(cluster) citable)
assert "`e(citypes)'" == "t percentile bc"
seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(20) rseed(4) nodots cluster(cluster) normal citable)
assert "`e(citypes)'" == "normal percentile bc" & e(df_r) >= .
* a replay with citable adds the interval table (percentile alone at estimation: the default
* normal and bc join it), and returns the t intervals at another level
ci_case `"`ilf'"' "percentile" "citable"
log_has `"`ilf'"' "Key: P: percentile, BC: bias-corrected, N: normal"
assert r(found) == 1
assert "`e(citypes)'" == "percentile"
quietly seldecomp, level(90)
matrix R90 = r(ci_t)
local q = invttail(e(df_bs), .05)
assert reldif(R90[1, 1], el(e(b), 1, 1) - `q' * el(e(se), 1, 1)) < 1e-12
* the suboptions outside bayesboot() say where they go
foreach o in normal tinterval percentile bc {
    rcof "noi seldecomp y x1 [pw=wt], by(group) bayesboot(reps(2)) `o'" == 198
}
capture erase `"`ilf'"'
capture erase `"`ilf2'"'

di as res _n "PASSED: test_bootstrap.do"

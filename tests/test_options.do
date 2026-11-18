* test_options.do -- option handling, error messages, sample handling, generate(),
* the raking options, covariatedetail, display options, string and labelled by() variables,
* blocks; the staged modes (sequentially leaves the other blocks free, fixothermargins
* holds them at the Focal Group's margins), their parsing rules, and the refusal of the option
* words varbyvar, cumulative, isolated,
* nohold, holdothers, refgroup(), refminusfocal, noheader, notable, nodefinitions,
* raker() and ipfraking(), at estimation and at replay.
*
* Return codes are checked with the built-in -rcof-:  rcof "noi cmd" == #
*
* Run from the tests directory:   . do test_options.do

version 16
clear
set more off
set varabbrev off
discard             // drop any seldecomp already in memory; the repository copy is what gets tested

* ---- tolerances ---------------------------------------------------------------
local tol      1e-6    // reldif() tolerance for means and totals
local sharetol 1e-5    // raked shares vs Reference Group shares (the raking stops when the
                       // relative weight change is below tolerance(), seldecomp's
                       // default 1e-8; margins are checked with ctrltolerance(1e-6))
local itol     1e-10   // identities / identical reruns

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

capture program drop assert_mat
program define assert_mat
    // assert_mat matname row col expected [tol]   (row/col: stripe names or numbers)
    args mat row col expect tol
    if "`tol'" == "" local tol 1e-6
    capture confirm integer number `row'
    if _rc {
        local rnames : rowfullnames `mat'
        local r : list posof "`row'" in rnames
    }
    else local r `row'
    capture confirm integer number `col'
    if _rc {
        local cnames : colfullnames `mat'
        local c : list posof "`col'" in cnames
    }
    else local c `col'
    if `r' < 1 | `r' > rowsof(`mat') | `c' < 1 | `c' > colsof(`mat') {
        di as err "assert_mat: `mat'[`row', `col'] not found (row `r', column `c')"
        exit 9
    }
    tempname v
    scalar `v' = el(`mat', `r', `c')
    if missing(`v') {
        di as err "`mat'[`row', `col'] is missing; expected " `expect'
        exit 9
    }
    if reldif(`v', `expect') > `tol' {
        di as err "`mat'[`row', `col']: " %18.0g `v' " vs " %18.0g `expect' ///
            "  (reldif = " %9.3g reldif(`v', `expect') " > `tol')"
        exit 9
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

count if group == 0
local N0 = r(N)
count if group == 1
local N1 = r(N)
sum wt if group == 0, meanonly
local W1 = r(sum)
sum wt if group == 1, meanonly
local W2 = r(sum)


* ==============================================================================
* 1. error handling
* ==============================================================================
di as txt _n "---- 1. errors"

* fweights are not allowed (101)
gen int wti = ceil(wt)
assert wti >= 1
rcof "noi seldecomp y x1 [fw=wti], by(group)" == 101

* by() must take exactly two values (459)
gen byte g3 = group + (id <= 100)                 // 0, 1, 2
assert inlist(g3, 0, 1, 2)
rcof "noi seldecomp y x1, by(g3)" == 459
rcof "noi seldecomp y x1, by(x4)" == 459          // x4 has four levels
rcof "noi seldecomp y x1, by(x4) focalgroup(north)" == 459

* focalgroup() must be a value of by()
rcof "noi seldecomp y x1, by(group) focalgroup(5)" == 459
rcof "noi seldecomp y x1, by(group) focalgroup(abc)" == 198

* continuous variables
rcof "noi seldecomp y c.x3 x2, by(group)" == 198                     // c. without contvars()
rcof "noi seldecomp y c.x3 x2, by(group) contvars(x3 1)" == 198      // fewer than 2 groups
rcof "noi seldecomp y x3 x2, by(group) contvars(x3 4)" == 198        // not marked c.
rcof "noi seldecomp y c.x3 x2, by(group) contvars(x3 4, pooled focalgroup)" == 198
rcof "noi seldecomp y c.x3 x2, by(group) contvars(x3 four)" == 198
rcof "noi seldecomp y c.x4, by(group) contvars(x4 4)" == 108        // string in contvars()
rcof "noi seldecomp y c.x3 x2, by(group) contvars(x3 3 x3 5)" == 198  // listed twice
rcof "noi seldecomp y c.x3 x2, by(group) contvars(x3 4 c.x3 4)" == 198

* options that require sequentially (or individually); the two modes exclude each other
rcof "noi seldecomp y x1 x2, by(group) fixothermargins" == 198
rcof "noi seldecomp y x1 x2, by(group) sequentially individually" == 198
rcof "noi seldecomp y x1 x2, by(group) makeplot" == 198

* replace is a suboption of generate(), and level() and citable of bayesboot()
rcof "noi seldecomp y x1 x2, by(group) replace" == 198
rcof "noi seldecomp y x1 x2, by(group) generate(rw) replace" == 198
rcof "noi seldecomp y x1 x2, by(group) level(90)" == 198
rcof "noi seldecomp y x1 x2, by(group) bayesboot(reps(5)) citable" == 198

* mincell() must be a nonnegative integer
rcof "noi seldecomp y x1, by(group) mincell(-1)" == 198
* a noninteger is rejected by -syntax- itself (mincell(integer 10)); its return
* code is syntax's, so accept the codes syntax uses for a malformed number
capture noisily seldecomp y x1, by(group) mincell(2.5)
assert inlist(_rc, 7, 126, 198)

* the raking options: invalid values; an unknown option
rcof "noi seldecomp y x1 x2, by(group) foo" == 198
rcof "noi seldecomp y x1 x2, by(group) iterate(0)" == 198
rcof "noi seldecomp y x1 x2, by(group) tolerance(0)" == 198
rcof "noi seldecomp y x1 x2, by(group) trimfrequency(often)" == 198   // needs a trim option
rcof "noi seldecomp y x1 x2, by(group) trimhirel(5) trimfrequency(never)" == 198

* bayesboot() suboptions
rcof "noi seldecomp y x1, by(group) bayesboot(reps(1))" == 198
rcof "noi seldecomp y x1, by(group) bayesboot(foo)" == 198
rcof "noi seldecomp y x1, by(group) bayesboot(rseed(1.5))" == 198

* raking variable list
rcof "noi seldecomp y ib2.x1, by(group)" == 198          // only i. and c. operators
rcof "noi seldecomp y, by(group)" == 102                 // no raking variables
rcof "noi seldecomp y nosuchvar, by(group)" == 111
rcof "noi seldecomp x4 x1, by(group)" == 108             // string depvar
rcof "noi seldecomp y (a: x1) (a: x2), by(group) sequentially" == 198        // duplicate block name
rcof "noi seldecomp y (overall: x1) x2, by(group) sequentially" == 198       // reserved block name
rcof "noi seldecomp y () x2, by(group) sequentially" == 198                  // empty block

* negative weights (402)
gen double wneg = wt
replace wneg = -1 in 1
rcof "noi seldecomp y x1 [pw=wneg], by(group)" == 402
drop wneg


* ==============================================================================
* 2. support failure and zero-weight levels
* ==============================================================================
di as txt _n "---- 2. support"

* a level present in the Reference Group (group==1) but absent from the focal Focal Group (group==0)
* cannot be reached by raking: hard error 459 naming the level
gen x5 = x1
replace x5 = 9 if group == 1 & _n < 20
count if group == 1 & x5 == 9
assert r(N) > 0
count if group == 0 & x5 == 9
assert r(N) == 0
rcof "noi seldecomp y x5, by(group)" == 459
rcof "noi seldecomp y x5 [pw=wt], by(group)" == 459
rcof "noi seldecomp y x1 x5 [pw=wt], by(group) sequentially" == 459
* every Focal Group observation carries a level absent from the Reference Group on one of two variables while
* every Reference Group level is still present in the Focal Group: no support failure, but nothing is left to
* rake once the zero-target observations are set aside (Engine rc 2000).
gen byte x6 = cond(group == 0 & mod(_n, 2) == 0, 9, x1)
gen byte x7 = cond(group == 0 & mod(_n, 2) == 1, 9, x2)
rcof "noi seldecomp y x6 x7, by(group)" == 2000
drop x6 x7
* with the focal group swapped the same level is in the Focal Group only: allowed (see below)
seldecomp y x5, by(group) focalgroup(1)
assert e(N_focal_zero) > 0

* mirror case: a level present only in the focal group gets zero raked weight,
* the run succeeds with a note, and generate() stores 0 for those observations
replace x5 = x1
replace x5 = 9 if group == 0 & _n < 20
count if group == 0 & x5 == 9
local nz = r(N)
assert `nz' > 0
seldecomp y x5 [pw=wt], by(group) generate(rw5)
assert e(N_focal_zero) == `nz'
assert e(N_raked)  == e(N_focal) - `nz'
assert e(N_focal)      == `N0'
assert e(N_notes)  > 0
count if group == 0 & x5 == 9 & rw5 == 0
assert r(N) == `nz'
count if group == 0 & x5 != 9 & !(rw5 > 0 & rw5 < .)
assert r(N) == 0
* Reference Group observations receive their base weight
count if group == 1 & rw5 != wt
assert r(N) == 0
* out-of-support accounting: explained = change of support + reweighting within support
assert e(pct_W_focal_zero) > 0 & e(pct_W_focal_zero) < 100
assert !missing(e(mean_focal_zero)) & !missing(e(mean_focal_cs))
assert reldif(e(explained_support) + e(explained_cs), e(explained)) < `itol'
assert reldif(e(mean_focal) - e(mean_focal_cs), e(explained_support)) < `itol'
sum y [aw=wt] if group == 0 & x5 == 9, meanonly
local m1z = r(mean)
assert_close e(mean_focal_zero) `m1z' `tol' "mean of y among zero-weight Focal Group obs"
sum y [aw=wt] if group == 0 & x5 != 9, meanonly
local m1cs = r(mean)
assert_close e(mean_focal_cs) `m1cs' `tol' "Focal Group mean on the raking sample (base weights)"
* the raked weights still sum to the focal group's total weight
sum rw5 if group == 0, meanonly
local srw5 = r(sum)
assert_close `srw5' `W1' `tol' "sum of raked weights with zero-weight level"
* the raked shares of the remaining levels match the Reference Group
foreach l in 1 2 3 {
    sum rw5 if group == 0 & x5 == `l', meanonly
    local s1 = r(sum) / `W1'
    sum wt if group == 1 & x5 == `l', meanonly
    local s2 = r(sum) / `W2'
    assert reldif(`s1', `s2') < `sharetol'
}
drop rw5 x5


* ==============================================================================
* 3. estimation sample: if, missing values, zero weights
* ==============================================================================
di as txt _n "---- 3. sample"

seldecomp y x1 if x2 == 1, by(group)
count if x2 == 1
assert e(N) == r(N)
count if e(sample)
assert r(N) == e(N)
count if e(sample) & x2 == 0
assert r(N) == 0

seldecomp y x1 in 1/1000, by(group)
assert e(N) == 1000

preserve
replace y = . in 1/5
seldecomp y x1, by(group)
assert e(N) == 1995
replace x1 = . in 6/10
seldecomp y x1, by(group)
assert e(N) == 1990
replace group = . in 11/15
seldecomp y x1, by(group)
assert e(N) == 1985
replace x4 = "" in 16/20
seldecomp y x1, by(group)
assert e(N) == 1985                      // x4 is not in the model
seldecomp y x4, by(group)
assert e(N) == 1985                      // 5 (y) + 5 (group) + 5 (empty x4); x1 is not in the model
restore

gen double wz = wt
replace wz = 0 in 1/5
seldecomp y x1 [pw=wz], by(group)
assert e(N) == 1995                      // zero-weight observations are excluded
count if e(sample)
assert r(N) == 1995
drop wz


* ==============================================================================
* 4. generate() and replace
* ==============================================================================
di as txt _n "---- 4. generate()"

seldecomp y x1 x2 [pw=wt], by(group) generate(rw)
assert "`e(generate)'" == "rw"
confirm variable rw
assert "`: type rw'" == "double"
* base weights for the Reference Group, positive raked weights for the Focal Group
count if group == 1 & rw != wt
assert r(N) == 0
count if group == 0 & !(rw > 0 & rw < .)
assert r(N) == 0
* total weight preserved
sum rw if group == 0, meanonly
local srw = r(sum)
assert_close `srw' `W1' `tol' "sum of raked weights = W_focal"
assert_close e(W_focal) `W1' `tol' "e(W_focal)"
assert_close e(W_ref) `W2' `tol' "e(W_ref)"
* raked shares of the Focal Group equal Reference Group's shares, variable by variable
foreach v in x1 x2 {
    levelsof `v', local(lv)
    foreach l of local lv {
        sum rw if group == 0 & `v' == `l', meanonly
        local s1 = r(sum) / `srw'
        sum wt if group == 1 & `v' == `l', meanonly
        local s2 = r(sum) / `W2'
        if reldif(`s1', `s2') > `sharetol' {
            di as err "raked share of `v' == `l': `s1' vs Reference Group share `s2'"
            exit 9
        }
    }
}
* the raked mean computed from rw equals e(mean_raked_focal)
sum y [aw=rw] if group == 0, meanonly
local mrw = r(mean)
assert_close `mrw' e(mean_raked_focal) `itol' "mean of y under rw"
* e(raked_*) diagnostics agree with the generated variable
sum rw if group == 0
local rw_min  = r(min)
local rw_max  = r(max)
local rw_mean = r(mean)
local rw_sd   = r(sd)
assert_close `rw_min'  e(raked_min)  `itol' "raked_min"
assert_close `rw_max'  e(raked_max)  `itol' "raked_max"
assert_close `rw_mean' e(raked_mean) `itol' "raked_mean"
assert_close `rw_sd'   e(raked_sd)   `itol' "raked_sd"
* Kish effective sample sizes on the raking sample: base weights, raked weights, DEFF
gen double wt2 = wt^2
sum wt2 if group == 0, meanonly
local neffb = `W1'^2 / r(sum)
gen double rw2 = rw^2
sum rw2 if group == 0, meanonly
local neffr = `srw'^2 / r(sum)
assert_close e(N_eff_base)  `neffb' `itol' "N_eff_base"
assert_close e(N_eff_raked) `neffr' `itol' "N_eff_raked"
assert_close e(ess_ratio) `=`neffb' / `neffr'' `itol' "ess_ratio"
drop wt2 rw2
* provenance characteristics attached to the generated variable
assert "`: char rw[seldecomp_stage]'"     == "1"
assert "`: char rw[seldecomp_converged]'" == "1"
assert "`: char rw[seldecomp_rakevars]'"  == "x1 x2"
assert strpos(`"`: char rw[seldecomp_cmd]'"', "seldecomp y x1 x2 [pw=wt], by(group) generate(rw)") == 1
assert "`: variable label rw'" == "Raked weights (focal group), base weights (reference group)"

* rerun without replace fails (110), with replace succeeds
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) generate(rw)" == 110
seldecomp y x1 x2 [pw=wt], by(group) generate(rw, replace)
confirm variable rw
drop rw

* sequentially with three blocks: rw_s1, rw_s2 and rw (stage 3 = full specification)
seldecomp y x1 [pw=wt], by(group) generate(rw1)
seldecomp y x1 x2 x4 [pw=wt], by(group) generate(rwfull)
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially generate(rw)
confirm variable rw_s1 rw_s2 rw
capture confirm variable rw_s3
assert _rc == 111
assert reldif(rw_s1, rw1)  < `tol' if group == 0
assert reldif(rw, rwfull)  < `tol' if group == 0
count if group == 1 & (rw_s1 != wt | rw_s2 != wt | rw != wt)
assert r(N) == 0
assert "`: char rw_s1[seldecomp_stage]'" == "1" & "`: char rw_s2[seldecomp_stage]'" == "2" & "`: char rw[seldecomp_stage]'" == "3"
assert "`: char rw_s1[seldecomp_rakevars]'" == "x1"
assert "`: char rw_s2[seldecomp_rakevars]'" == "x1 x2"
assert "`: char rw[seldecomp_rakevars]'"    == "x1 x2 x4"
rcof "noi seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially generate(rw)" == 110
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially generate(rw, replace)
confirm variable rw_s1 rw_s2 rw
* generate() is checked before any computation: an existing stage name blocks the run
drop rw rw_s2
rcof "noi seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially generate(rw)" == 110
drop rw_s1 rw1 rwfull
capture noisily seldecomp y x1 x2, by(group) generate(1bad)
assert inlist(_rc, 7, 198)               // invalid name, rejected by -syntax- or the ado


* ==============================================================================
* 5. the raking options
* ==============================================================================
di as txt _n "---- 5. raking options"

* one iteration is not enough for two margins: the targets are missed, which is an error
* unless force is given; with force the results are shown and e(converged) == 0
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) iterate(1)" == 430
seldecomp y x1 x2 [pw=wt], by(group) iterate(1) force
assert e(converged) == 0
assert e(N_notes) > 0
assert strpos("`e(rakeopts)'", "iterate(1)") > 0
* replay still works after non-convergence
seldecomp

* tight tolerance
seldecomp y x1 x2 [pw=wt], by(group) tolerance(1e-10) ctrltolerance(1e-8)
assert e(converged) == 1
assert e(badcontrols) == 0
assert strpos("`e(rakeopts)'", "tolerance(1e-10)") > 0
assert_close e(explained) exp_C_explained `tol' "explained with tolerance(1e-10)"

* the trimming options are accepted
seldecomp y x1 x2 [pw=wt], by(group) trimhirel(10) trimlorel(0.1) trimfrequency(sometimes)
assert strpos("`e(rakeopts)'", "trimhirel(10)") > 0
seldecomp y x1 x2 [pw=wt], by(group) trimhiabs(20) trimloabs(0.05)
* every bound binds: a weight is at most the smaller of trimhirel() times its base
* weight and trimhiabs(), and at least the larger of trimlorel() times it and trimloabs(); the
* base weights run from .5 to 3, so each of the four is the binding one for some observations
seldecomp y x1 x2 [pw=wt], by(group) trimhirel(1.25) trimhiabs(2.6) trimlorel(.8) trimloabs(.6) ///
    generate(rwb) force
gen double hi = min(1.25 * wt, 2.6)
gen double lo = max(.8 * wt, .6)
count if group == 0 & e(sample) & (rwb > hi * (1 + 1e-9) | rwb < lo * (1 - 1e-9))
assert r(N) == 0
count if group == 0 & e(sample) & abs(rwb / (1.25 * wt) - 1) < 1e-9 & 1.25 * wt < 2.6
assert r(N) > 0
count if group == 0 & e(sample) & abs(rwb / 2.6 - 1) < 1e-9 & 2.6 < 1.25 * wt
assert r(N) > 0
count if group == 0 & e(sample) & abs(rwb / (.8 * wt) - 1) < 1e-9 & .8 * wt > .6
assert r(N) > 0
count if group == 0 & e(sample) & abs(rwb / .6 - 1) < 1e-9 & .6 > .8 * wt
assert r(N) > 0
drop rwb hi lo
* bounds that contradict for some observation (a lower bound above an upper one) are refused:
* trimloabs() above trimhirel() times the smallest base weights, trimlorel() times the largest
* above trimhiabs(); in a bayesboot() replicate they fail it (below, with the messages)
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) trimhirel(1.25) trimloabs(.9)" == 198
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) trimlorel(.8) trimhiabs(2)" == 198
* a bound that binds on feasible margins and keeps them from their targets gives rc 430, as
* any other missed margin; force keeps the results, the miss reported by
* e(targets_met) = 0, e(badcontrols) and e(converged) = 0
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) trimhirel(1.05)" == 430
seldecomp y x1 x2 [pw=wt], by(group) trimhirel(1.05) force
assert e(badcontrols) > 0 & e(N_at_bound) > 0
assert e(converged) == 0 & e(targets_met) == 0
* margins that cannot be met (Focal Group has a == 1 <=> b == 0):
* the raking settles after one cycle with the adjustment factors .81 and 1.19, a miss
preserve
clear
set obs 400
gen byte grp = _n > 200
gen byte a = cond(grp, mod(_n, 10) < 3, mod(_n, 2))
gen byte b = cond(grp, mod(floor(_n / 3), 10) < 6, 1 - mod(_n, 2))
gen double yy = a + b + mod(_n, 7) / 7
* bounds that never bind give 430 as without them, whatever trimfrequency(); also at stage 2
* of sequentially fixothermargins
rcof "noi seldecomp yy i.a i.b, by(grp) trimhirel(1e6)" == 430
rcof "noi seldecomp yy i.a i.b, by(grp) trimhirel(1.5) trimhiabs(1.3)" == 430
rcof "noi seldecomp yy i.a i.b, by(grp) trimhiabs(1e6) trimloabs(1e-6) trimfrequency(once)" == 430
rcof "noi seldecomp yy i.a i.b, by(grp) trimhirel(1e6) trimfrequency(often)" == 430
rcof "noi seldecomp yy (a) (b), by(grp) sequentially fixothermargins trimhirel(1e6)" == 430
* a bound that the factor .81 sits on (e(N_at_bound) counts the weights at a bound within a
* relative 1e-6) but that changes no weight: 430 too
rcof "noi seldecomp yy i.a i.b, by(grp) trimlorel(.81)" == 430
seldecomp yy i.a i.b, by(grp) trimlorel(.81) force
assert e(N_at_bound) > 0 & e(badcontrols) > 0 & e(converged) == 0 & e(targets_met) == 0
* with force the results are those without the trimming option
seldecomp yy i.a i.b, by(grp) trimhirel(1e6) force generate(rwt)
assert e(badcontrols) > 0 & e(N_at_bound) == 0
assert e(converged) == 0
assert el(e(bystage), 1, colnumb(e(bystage), "converged")) == 0
assert "`: char rwt[seldecomp_converged]'" == "0"
matrix B1 = e(b)
seldecomp yy i.a i.b, by(grp) force generate(rwn)
assert mreldif(e(b), B1) < 1e-12
assert reldif(rwt, rwn) < 1e-12
drop rwt rwn
* bounds that bind (factors 1.19 -> 1.1 and .81 -> .9): rc 430 too; with
* force one note (the margins trimming left unmet, the worst, the weights at a bound, kept
* by force), the imbalance left to unexplained
rcof "noi seldecomp yy i.a i.b, by(grp) trimhirel(1.1) trimlorel(.9) generate(rwt)" == 430
capture confirm variable rwt
assert _rc
seldecomp yy i.a i.b, by(grp) trimhirel(1.1) trimlorel(.9) generate(rwt) force
assert e(badcontrols) > 0 & e(N_at_bound) > 0 & e(converged) == 0 & e(targets_met) == 0
assert el(e(bystage), 1, colnumb(e(bystage), "converged")) == 0
assert "`: char rwt[seldecomp_converged]'" == "0"
assert e(N_notes) == 1
local nb = e(badcontrols)
local nt = strtrim(string(e(N_at_bound), "%12.0fc"))
assert strpos(`"`e(note1)'"', "trimming left `nb' raking margin") == 1
local ntail "; `nt' weights at a bound); results kept because force was specified"
assert substr(`"`e(note1)'"', -strlen("`ntail'"), .) == "`ntail'"
assert strpos(`"`e(note1)'"', "the imbalance is in unexplained") == 0 & strpos(`"`e(note1)'"', "widen the bounds") == 0
sum rwt if grp == 0
assert r(min) >= .9 * (1 - 1e-6) & r(max) <= 1.1 * (1 + 1e-6)
drop rwt
* ... per stage: stage 1 of sequentially fixothermargins (a at the Reference Group's shares, b held at
* Focal Group's) ends each cycle at the base weights, so the same bounds change nothing there: 430
* at stage 1
rcof "noi seldecomp yy (a) (b), by(grp) sequentially fixothermargins trimhirel(1.1) trimlorel(.9)" == 430
* the same rule in bayesboot() replicates: with a bound that never binds and with a binding
* one every replicate misses its targets, and force keeps them (e(N_forced), and a note)
tempfile lgforce
log using `"`lgforce'"', text replace name(sdforce)
seldecomp yy i.a i.b, by(grp) trimhirel(1e6) force bayesboot(reps(5) rseed(1) nodots)
log close sdforce
assert e(N_reps) == 5 & e(N_misreps) == 0 & e(N_forced) == 5
mata: assert(sum(strpos(cat("`lgforce'"), "note: force kept 5 of 5 replications that missed raking targets") :> 0) == 1)
rcof "noi seldecomp yy i.a i.b, by(grp) trimhirel(1.1) trimlorel(.9) bayesboot(reps(5) rseed(1) nodots)" == 430
seldecomp yy i.a i.b, by(grp) trimhirel(1.1) trimlorel(.9) force bayesboot(reps(5) rseed(1) nodots)
assert e(N_reps) == 5 & e(N_misreps) == 0 & e(N_forced) == 5
restore
* a bound that the full sample can just satisfy (trimhirel(2.4): weights at the bound, every
* margin met, so no force is needed; 2.3 leaves margins unmet) and some replicates cannot: the
* replicates whose margins the bound leaves unmet fail and are dropped; force keeps them,
* and without force e(N_forced) is not posted
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) trimhirel(2.3)" == 430
seldecomp y x1 x2 [pw=wt], by(group) trimhirel(2.4) bayesboot(reps(40) rseed(11) nodots)
assert e(converged) == 1 & e(targets_met) == 1 & e(N_at_bound) > 0
assert e(N_forced) == .
local nmis = e(N_misreps)
assert `nmis' > 0 & `nmis' < 40 & e(N_reps) == 40 - `nmis'
seldecomp y x1 x2 [pw=wt], by(group) trimhirel(2.4) force bayesboot(reps(40) rseed(11) nodots)
assert e(converged) == 1 & e(targets_met) == 1
assert e(N_reps) == 40 & e(N_misreps) == 0 & e(N_forced) == `nmis'

* default options recorded in e(rakeopts): iterate(10000) tolerance(1e-8) ctrltolerance(1e-6)
seldecomp y x1 x2 [pw=wt], by(group)
assert "`e(rakeopts)'" == "iterate(10000) tolerance(1e-8) ctrltolerance(1e-6)"

* noisily shows the raking's details
seldecomp y x1 x2 [pw=wt], by(group) noisily
assert_close e(explained) exp_C_explained `tol' "explained with noisily"


* ==============================================================================
* 6. covariatedetail
* ==============================================================================
di as txt _n "---- 6. covariatedetail"

seldecomp y x1 x2 [pw=wt], by(group) covariatedetail
matrix C = e(covariates)
assert rowsof(C) == 5                     // 3 levels of x1 + 2 levels of x2
assert colsof(C) == 5                     // focal_group ref_group raked_focal N_focal N_ref
local cn : colnames C
assert "`cn'" == "focal_group ref_group raked_focal N_focal N_ref"
forvalues r = 1/5 {
    assert el(C, `r', 1) > 0 & el(C, `r', 1) < 1
    assert el(C, `r', 2) > 0 & el(C, `r', 2) < 1
    if reldif(el(C, `r', 2), el(C, `r', 3)) > `sharetol' {
        di as err "e(covariates) row `r': raked share " el(C, `r', 3) " vs Reference Group share " el(C, `r', 2)
        exit 9
    }
}
* shares sum to one within each variable
assert reldif(el(C,1,1) + el(C,2,1) + el(C,3,1), 1) < `itol'
assert reldif(el(C,1,2) + el(C,2,2) + el(C,3,2), 1) < `itol'
assert reldif(el(C,1,3) + el(C,2,3) + el(C,3,3), 1) < `itol'
assert reldif(el(C,4,1) + el(C,5,1), 1) < `itol'
* column 1 is the weighted Focal Group share, column 2 the weighted Reference Group share
sum wt if group == 0 & x1 == 1, meanonly
local sh11 = r(sum) / `W1'
assert_close el(C,1,1) `sh11' `tol' "e(covariates)[1,1]"
sum wt if group == 1 & x2 == 1, meanonly
local sh52 = r(sum) / `W2'
assert_close el(C,5,2) `sh52' `tol' "e(covariates)[5,2]"
* columns 4 and 5 are the unweighted counts
count if group == 0 & x1 == 1
assert el(C,1,4) == r(N)
count if group == 1 & x1 == 1
assert el(C,1,5) == r(N)
count if group == 0 & x2 == 1
assert el(C,5,4) == r(N)
assert el(C,1,4) + el(C,2,4) + el(C,3,4) == `N0'
assert el(C,1,5) + el(C,2,5) + el(C,3,5) == `N1'
* replay with covariatedetail
seldecomp, covariatedetail
* sequentially: rows = 3 + 2 + 4 levels
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially covariatedetail
matrix C = e(covariates)
assert rowsof(C) == 9 & colsof(C) == 5
forvalues r = 1/9 {
    if reldif(el(C, `r', 2), el(C, `r', 3)) > `sharetol' {
        di as err "e(covariates) row `r' (sequentially): raked share " el(C, `r', 3) " vs Reference Group share " el(C, `r', 2)
        exit 9
    }
}
* interaction term: 6 joint cells
seldecomp y x1##x2 [pw=wt], by(group) covariatedetail
matrix C = e(covariates)
assert rowsof(C) == 6


* ==============================================================================
* 7. display options and replay
* ==============================================================================
di as txt _n "---- 7. display"

* runs without error, at estimation and at replay
capture noisily seldecomp y x1 x2 [pw=wt], by(group)
local rc = _rc
if `rc' {
    di as err "seldecomp ... failed with rc `rc'"
    exit 9
}
assert_close e(explained) exp_C_explained `tol' "explained, default display"
capture noisily seldecomp
local rc = _rc
if `rc' {
    di as err "replay: seldecomp, failed with rc `rc'"
    exit 9
}
* these words are not display options: refused (198), at estimation and at replay, and the
* results in memory stay as they are
foreach o in noheader notable nodefinitions nodiagnostics fulldiagnostics {
    rcof "noi seldecomp y x1 x2 x4 [pw=wt], by(group) `o'" == 198
    assert "`e(rakevars)'" == "x1 x2"
    rcof "noi seldecomp, `o'" == 198
}
seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(5) rseed(1) nodots level(90))
assert e(level) == 90
seldecomp, level(90)
seldecomp, level(99)
seldecomp
seldecomp
assert e(level) == 90                      // replay does not change the stored level
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially
seldecomp
* The compact definition note is shown at estimation and replay, and searched for
* in a text log (the same technique as
* test_reference.do uses for the string-conversion message)
capture program drop log_has
program define log_has, rclass
    // log_has "filename" "text"  --  r(found) = 1 if some line of the file contains text
    args fn txt
    tempname fh
    local found 0
    file open `fh' using `"`fn'"', read text
    file read `fh' line
    while r(eof) == 0 {
        if strpos(`"`macval(line)'"', `"`txt'"') local found 1
        file read `fh' line
    }
    file close `fh'
    return scalar found = `found'
end
tempfile lg
local lgf `"`lg'.log"'
capture log close sdlog
log using `"`lgf'"', text replace name(sdlog)
seldecomp y x1 x2 [pw=wt], by(group)
log close sdlog
log_has `"`lgf'"' "Difference = focal - reference. Explained = focal - raked focal."
assert r(found) == 1
log using `"`lgf'"', text replace name(sdlog)
seldecomp                   // the replay shows them too
log close sdlog
log_has `"`lgf'"' "Difference = focal - reference. Explained = focal - raked focal."
assert r(found) == 1
log using `"`lgf'"', text replace name(sdlog)
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially
log close sdlog
log_has `"`lgf'"' "Difference = focal - reference. Explained = focal - raked focal."
assert r(found) == 1
capture erase `"`lgf'"'
* replay is refused when the last estimates are not seldecomp's
regress y x1
rcof "noi seldecomp" == 301
rcof "noi seldecomp, level(90)" == 301
* cmdline and estimates store/replay
seldecomp y x1 x2 [pw=wt], by(group)
assert strpos(`"`e(cmdline)'"', "seldecomp y x1 x2 [pw=wt], by(group)") == 1
estimates store sd1
regress y x1
estimates restore sd1
assert "`e(cmd)'" == "seldecomp"
assert_close e(explained) exp_C_explained `tol' "explained after estimates restore"
* replay of the restored results
capture noisily seldecomp
local rc = _rc
if `rc' {
    di as err "replay after estimates restore failed with rc `rc'"
    exit 9
}
assert "`e(cmd)'" == "seldecomp"
estimates replay sd1
estimates drop sd1
* the same with the conventional name m1 and a staged model in between
seldecomp y x1 x2 [pw=wt], by(group)
estimates store m1
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially
regress y x1
estimates restore m1
assert "`e(cmd)'" == "seldecomp"
assert e(N_stages) == 1
assert_close e(explained) exp_C_explained `tol' "explained after estimates restore m1"
capture noisily seldecomp
local rc = _rc
if `rc' {
    di as err "replay after estimates restore m1 failed with rc `rc'"
    exit 9
}
estimates drop m1


* ==============================================================================
* 8. blocks, named blocks, parentheses without sequentially
* ==============================================================================
di as txt _n "---- 8. blocks"

seldecomp y (educ: x1 x2) x4 [pw=wt], by(group) sequentially
assert e(N_stages) == 2 & e(N_blocks) == 2 & e(N_terms) == 3
assert "`e(blocks)'" == "educ x4"
assert "`e(stagenames)'" == "educ x4"
assert "`e(stagetype)'" == "sequentially" & "`e(fixothermargins)'" == ""
assert "`e(block1)'" == "x1 x2"
assert "`e(block2)'" == "x4"
assert "`e(stage1)'" == "x1 x2"
assert "`e(stage2)'" == "x1 x2 x4"
matrix S = e(bystage)
assert rowsof(S) == 2
assert_mat S 1 cum_expl exp_C_explained  `tol'      // stage 1 = x1 x2 jointly = case C
assert_mat S 2 cum_expl exp_F3_explained `tol'      // stage 2 = full specification
matrix b = e(b)
assert colsof(b) == 14
assert_mat b 1 marg_expl:educ exp_C_explained  `tol'
assert_mat b 1 marg_expl:x4   exp_F3_increment `tol'
assert_mat b 1 cum_expl:x4    exp_F3_explained `tol'
assert_close e(explained) exp_F3_explained `tol' "explained with blocks"

* unnamed block takes the name of its first term
seldecomp y (x1 x2) x4 [pw=wt], by(group) sequentially
assert "`e(blocks)'" == "x1 x4"

* an interaction is one block; its default block name replaces # by _
seldecomp y x1##x2 x4 [pw=wt], by(group) sequentially
assert e(N_stages) == 2
assert "`e(blocks)'" == "x1_x2 x4"
matrix S = e(bystage)
assert_mat S 1 cum_expl exp_D_explained `tol'

* parentheses are ignored without sequentially or individually
seldecomp y (x1 x2) x4 [pw=wt], by(group)
assert e(N_stages) == 1
assert "`e(blocks)'" == "x1 x4"
assert "`e(stagenames)'" == "overall"
assert "`e(stagetype)'" == ""
assert_close e(explained) exp_F3_explained `tol' "explained, parentheses without sequentially"

* sequentially with a single block equals the unstaged run
seldecomp y x1 [pw=wt], by(group) sequentially
assert e(N_stages) == 1
assert_close e(explained) exp_B_explained `tol' "sequentially with one block"
matrix b = e(b)
assert colsof(b) == 10
assert_mat b 1 marg_expl:x1 exp_B_explained `tol'


* ==============================================================================
* 9. string by() variable and focalgroup(string)
* ==============================================================================
di as txt _n "---- 9. string by()"

gen g2 = cond(group == 1, "treated", "control")
confirm string variable g2
seldecomp y x1, by(g2)
assert "`e(focal_group)'" == "control"          // alphabetically first = default focal group
assert "`e(ref_group)'" == "treated"
assert "`e(focal_group_label)'" == "control"
assert e(N_focal) == `N0'
local expl_c = e(explained)
seldecomp y x1, by(group)
assert_close e(explained) `expl_c' `itol' "string by() = numeric by()"

seldecomp y x1, by(g2) focalgroup(treated)
assert "`e(focal_group)'" == "treated"
assert "`e(ref_group)'" == "control"
assert e(N_focal) == `N1'
local expl_t = e(explained)
seldecomp y x1, by(group) focalgroup(1)
assert_close e(explained) `expl_t' `itol' "focalgroup(treated) = focalgroup(1)"
rcof "noi seldecomp y x1, by(g2) focalgroup(other)" == 459
rcof "noi seldecomp y x1, by(g2) focalgroup(1)" == 459
* focalgroup() naming the default Focal Group gives the default
seldecomp y x1, by(g2) focalgroup(control)
assert "`e(focal_group)'" == "control" & "`e(ref_group)'" == "treated"
assert_close e(explained) `expl_c' `itol' "focalgroup(control) = default"
drop g2


* ==============================================================================
* 10. value-labelled by() variable
* ==============================================================================
di as txt _n "---- 10. value labels"

preserve
label define gl 0 "men" 1 "women"
label values group gl
seldecomp y x1, by(group)
assert "`e(focal_group)'" == "0"
assert "`e(ref_group)'" == "1"
assert strpos(`"`e(focal_group_label)'"', "men") > 0
assert strpos(`"`e(focal_group_label)'"', "women") == 0
assert strpos(`"`e(ref_group_label)'"', "women") > 0
assert_close e(explained) `expl_c' `itol' "labelled by() = unlabelled by()"
seldecomp y x1, by(group) focalgroup(1)
assert strpos(`"`e(focal_group_label)'"', "women") > 0
* focalgroup() also takes the value label
seldecomp y x1, by(group) focalgroup(women)
assert "`e(focal_group)'" == "1" & strpos(`"`e(focal_group_label)'"', "women") > 0
assert_close e(explained) `expl_t' `itol' "focalgroup(women) = focalgroup(1)"
* value labels on a raking variable appear in its level labels
label define x1l 1 "low" 2 "mid" 3 "high"
label values x1 x1l
seldecomp y x1, by(group)
assert strpos(`"`e(levels1)'"', "high") > 0
restore

* ==============================================================================
* 11. mincell(): thin Focal Group cells (in the first 100 rows, x1 == 3 has 8
*     Focal Group observations and the other Focal Group cells have at least 10;
*     Reference Group cells are not checked -- x1 == 1 has only 6 there)
* ==============================================================================
di as txt _n "---- 11. mincell()"

count if group == 0 & x1 == 3 in 1/100
local nthin = r(N)
assert `nthin' < 10 & `nthin' > 0
count if group == 0 & x1 == 1 in 1/100
assert r(N) >= 10
count if group == 0 & x1 == 2 in 1/100
assert r(N) >= 10
seldecomp y x1 in 1/100, by(group)
assert e(mincell) == 10
assert e(N_focal_mincell) == 1
assert e(N_notes) == 2                    // one heading line, one cell line
assert `"`e(note1)'"' == "thin target level:"            // the definition is in the diagnostics footnote
assert strpos(`"`e(note2)'"', "x1 = 3") > 0 & strpos(`"`e(note2)'"', "(`nthin' obs)") > 0
matrix S = e(bystage)
assert el(S, 1, colnumb(S, "N_focal_mincell")) == 1
* the note does not change the estimates
local expl_100 = e(explained)
seldecomp y x1 in 1/100, by(group) mincell(0)
assert e(mincell) == 0
assert e(N_focal_mincell) == 0
assert e(N_notes) == 0
assert_close e(explained) `expl_100' `itol' "mincell(0) leaves the estimates unchanged"
seldecomp y x1 in 1/100, by(group) mincell(`nthin')
assert e(N_focal_mincell) == 0                // "fewer than #" is strict
seldecomp y x1 in 1/100, by(group) mincell(`=`nthin' + 1')
assert e(N_focal_mincell) == 1
* the full sample has no thin target level (see the cell counts in section 6)
seldecomp y x1 x2, by(group)
assert e(N_focal_mincell) == 0
assert e(N_notes) == 0


* ==============================================================================
* 12. by prefix (byable(recall)): one decomposition per by-group; generate(),
*     makeplot and bayesboot(saving()|generate()) are refused (190).
*     This section sorts the data (bysort); sort id restores the order after it.
* ==============================================================================
di as txt _n "---- 12. by prefix"

rcof "noi bysort x2: seldecomp y x1, by(group) generate(rwby)" == 190
rcof "noi bysort x2: seldecomp y x1, by(group) sequentially makeplot" == 190
rcof "noi bysort x2: seldecomp y x1, by(group) bayesboot(reps(2) saving(byreps))" == 190
rcof "noi bysort x2: seldecomp y x1, by(group) bayesboot(reps(2) generate(bwby))" == 190
capture confirm variable rwby
assert _rc == 111
* the last by-group's results are left behind and equal the corresponding if-run
seldecomp y x1 [pw=wt] if x2 == 1, by(group)
local expl_x2 = e(explained)
local N_x2 = e(N)
bysort x2: seldecomp y x1 [pw=wt], by(group)
assert "`e(cmd)'" == "seldecomp"
assert e(N) == `N_x2'
count if e(sample)
assert r(N) == `N_x2'
assert_close e(explained) `expl_x2' `itol' "by x2: last group = if x2 == 1"
* the internal by-group options handed to Estimate never reach e(cmdline)
assert strpos(`"`e(cmdline)'"', "seldby") == 0
* bayesboot() itself is allowed under by; the by-group restriction must survive Bootstrap's own
* sorts and preserve/restore (draws differ from an if-run because the sort order differs)
bysort x2: seldecomp y x1 [pw=wt], by(group) bayesboot(reps(3) rseed(1) nodots)
assert e(N) == `N_x2'
assert e(N_reps) + e(N_misreps) == 3
* replay under by is refused (190)
rcof "noi bysort x2: seldecomp" == 190
sort id

* ==============================================================================
* 13. miscellaneous checks
* ==============================================================================
di as txt _n "---- 13. miscellaneous checks"

* bayesboot() and makeplot() with empty parentheses are not ignored
seldecomp y x1 [pw=wt] in 1/400, by(group) bayesboot()
assert "`e(vce)'" == "bayesboot" & e(N_reps) == 100
seldecomp y x1 x2 in 1/400, by(group) sequentially makeplot()
graph describe seldecomp
graph drop _all
* non-integer double by() values are matched exactly (levelsof text would not round-trip)
preserve
gen double gd = cond(group, 0.1 + 0.2, 1/7)
seldecomp y x1 [pw=wt], by(gd)
assert e(N_focal) == 888 & e(N_ref) == 1112
seldecomp y x1 [pw=wt], by(gd) focalgroup(.3)
assert e(N_focal) == 1112
restore
* a generated block name longer than 32 characters falls back to block#
preserve
gen byte abcdefghijklmnopqrstuvwxyz_abcd = x1
gen byte abcdefghijklmnopqrstuvwxyz_abce = x2
seldecomp y abcdefghijklmnopqrstuvwxyz_abcd#abcdefghijklmnopqrstuvwxyz_abce x4, by(group) sequentially
assert "`: word 1 of `e(blocks)''" == "block1"
restore
* predict is refused (301); rseed() above 2^31-1 is refused up front (198)
seldecomp y x1, by(group)
rcof "noi predict double pr" == 301
rcof "noi seldecomp y x1, by(group) bayesboot(reps(2) rseed(3000000000))" == 198
* sequentially alone leaves the blocks not yet added free (stage 2 = x1 and x2 raked, x4 not
* constrained: case F); fixothermargins holds them at the Focal Group's own margins (case G,
* the papers' procedure); the minimal abbreviations seq and fix; fixothermargins needs a mode,
* and the two modes exclude each other
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially
assert "`e(stagetype)'" == "sequentially" & "`e(fixothermargins)'" == ""
assert_mat e(bystage) 2 cum_expl exp_F2_explained `tol'
matrix Sf = e(bystage)
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially fixothermargins
assert "`e(stagetype)'" == "sequentially" & "`e(fixothermargins)'" == "fixothermargins"
assert_mat e(bystage) 2 cum_expl exp_G2_explained `tol'
assert mreldif(e(bystage), Sf) > 1e-6
matrix Sv = e(bystage)
seldecomp y x1 x2 x4 [pw=wt], by(group) seq fix
assert mreldif(e(bystage), Sv) < 1e-12
seldecomp y x1 x2 x4 [pw=wt], by(group) seq
assert mreldif(e(bystage), Sf) < 1e-12
rcof "noi seldecomp y x1 x2 x4 [pw=wt], by(group) fixothermargins" == 198
rcof "noi seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially individually" == 198
rcof "noi seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially individually fixothermargins" == 198
* fixothermargins with a variable in two blocks is refused (198); without it the other blocks
* are free, and a variable in two blocks is allowed
rcof "noi seldecomp y (x1) (x1#x2), by(group) sequentially fixothermargins" == 198
seldecomp y (x1) (x1#x2), by(group) sequentially
assert e(N_stages) == 2 & "`e(fixothermargins)'" == ""
* unmet raking targets (infeasible margins) are an error unless force is given
preserve
clear
set obs 400
set seed 4
gen byte grp = _n > 200
gen byte a = cond(grp, mod(_n, 10) < 3, mod(_n, 2))          // focal: a == 1 <=> b == 0
gen byte b = cond(grp, mod(floor(_n / 3), 10) < 6, 1 - mod(_n, 2))
gen double yy = a + b + rnormal()
rcof "noi seldecomp yy i.a i.b, by(grp)" == 430
seldecomp yy i.a i.b, by(grp) force
assert e(converged) == 0 & e(badcontrols) > 0
restore
* the difference is always Focal Group - Reference Group (e(diffsign)), whichever group focalgroup() names
seldecomp y x1 x2 [pw=wt], by(group) sequentially
assert "`e(diffsign)'" == "focal_group - ref_group"
assert reldif(e(difference), e(mean_focal) - e(mean_ref)) < 1e-12
seldecomp y x1 x2 [pw=wt], by(group) sequentially focalgroup(1)
assert "`e(diffsign)'" == "focal_group - ref_group"
assert reldif(e(difference), e(mean_focal) - e(mean_ref)) < 1e-12
* the support-failure listing is capped and the bayesboot(saving()) run prints no stray lines
tempfile lgc
log using `lgc', text replace name(c2)
capture noisily seldecomp y x3, by(group)
seldecomp y x1 [pw=wt] in 1/300, by(group) bayesboot(reps(5) rseed(1) nodots saving(bsx, replace))
log close c2
mata: assert(sum(strpos(cat("`lgc'"), "x3 = ") :> 0) <= 12)
mata: assert(!any(strpos(cat("`lgc'"), "invalid name")))
capture erase bsx.dta
* makeplot at replay, also after a replication prefix (only e(b) survives it)
seldecomp y x1 x2 in 1/300, by(group) sequentially bayesboot(reps(10) rseed(1) nodots)
seldecomp, makeplot(, twoway(name(mpc, replace)))
graph describe mpc
bootstrap _b, reps(3) seed(1) nodots: seldecomp y x1 x2 in 1/300, by(group) sequentially
seldecomp, makeplot(, twoway(name(mpb, replace)))
graph describe mpb
graph drop _all
* contvars() naming a variable that is not a c. raking variable is refused (198)
rcof "noi seldecomp y x1, by(group) contvars(x3 3)" == 198
* a bare replay uses the estimation level
seldecomp y x1 [pw=wt] in 1/300, by(group) bayesboot(reps(5) rseed(1) nodots level(90))
seldecomp
assert e(level) == 90
* by() inside makeplot(, twoway()) is refused before estimating (198); a zero difference displays n/a
rcof "noi seldecomp y x1 x2, by(group) sequentially makeplot(, twoway(by(x2)))" == 198
* sdcheckonly is not a twoway option (the plotspec check before estimating is internal):
* twoway refuses it (198) before the estimation, as the twoway options are tried first, so
* the results in memory are kept; and at replay; no graph
capture graph drop _all
qui seldecomp y x1 x2 in 1/300, by(group) sequentially fixothermargins
matrix b_sdc = e(b)
rcof "noi seldecomp y x1 in 1/300, by(group) sequentially fixothermargins makeplot(, twoway(sdcheckonly))" == 198
assert "`e(cmd)'" == "seldecomp" & "`e(blocks)'" == "x1 x2" & mreldif(e(b), b_sdc) == 0
rcof "noi seldecomp, makeplot(, twoway(sdcheckonly))" == 198
rcof "noi seldecomp, makeplot(waterfalllevels, twoway(sdcheckonly))" == 198
capture graph describe seldecomp
assert _rc
capture graph drop _all
preserve
gen double yc = 1
seldecomp yc x1 x2, by(group) sequentially
assert e(pct_explained) == .
restore
* fewer than 2 clusters is refused (459)
gen byte onecl = 1
rcof "noi seldecomp y x1, by(group) bayesboot(reps(2) rseed(1) nodots cluster(onecl))" == 459
drop onecl
* no worst-fitting control is reported when every margin was met
seldecomp y x1 x2 [pw=wt], by(group)
assert e(badcontrols) == 0 & "`e(worstvar)'" == ""
* value labels and string values with quote, pipe, SMCL, dollar and backquote characters
* ($ and ` are written with char() so that this do-file's own macro expansion cannot touch them)
preserve
keep in 1/400
label define xl 1 `"a"b"' 2 "c|d" 3 "{bf:B} {hline}"
label values x1 xl
* (a do-file line re-expands `=char(36)' into $US before -label define- runs, so use Mata)
mata: st_vlmodify("gl", (0 \ 1), ("cost " + char(36) + "US" \ "it" + char(96) + "s"))
label values group gl
seldecomp y x1, by(group) covariatedetail
mata: assert(strpos(st_global("e(focal_group_label)"), char(36) + "US") > 0)
mata: assert(strpos(st_global("e(ref_group_label)"), char(96)) > 0)
mata: assert(strpos(st_global("e(levels1)"), "c/d") > 0 & strpos(st_global("e(levels1)"), "{hline}") > 0)
seldecomp
gen str12 gs = "b"
replace gs = "a " + char(36) + "US" if group == 0
seldecomp y x1, by(gs)
mata: assert(st_global("e(focal_group)") == "a " + char(36) + "US")
restore
* level(99.9) survives the replay
seldecomp y x1 [pw=wt] in 1/300, by(group) bayesboot(reps(5) rseed(1) nodots level(99.9))
assert e(level) == 99.9
seldecomp
* a failing run leaves the previous estimation results in place
gen byte x9 = x1
replace x9 = 4 if group == 1 & mod(_n, 50) == 0      // level 4 in the Reference Group only: support failure
regress y x3
rcof "noi seldecomp y (s1: x1) (s2: x9), by(group) sequentially" == 459
assert "`e(cmd)'" == "regress"
drop x9


* ==============================================================================
* messages for slips on the command line, group-value errors and the trimming
* note. The message texts are read from a text log written at linesize 255, so none is
* wrapped.
* ==============================================================================
di as txt _n "---- messages: command-line slips, group values, trimming"

capture program drop u2_msg
program define u2_msg
    // u2_msg "command" rc "text" ["text2" ["text3"]]: the command fails with rc and prints
    // text (and text2, text3); a text that begins with ! must not be printed
    args cmd rc t1 t2 t3
    tempfile lg
    local ls = c(linesize)
    set linesize 255
    log using `"`lg'"', text replace name(u2msg)
    capture noisily `cmd'
    local r = _rc
    log close u2msg
    set linesize `ls'
    if `r' != `rc' {
        di as err `"`cmd': rc `r', expected `rc'"'
        exit 9
    }
    foreach t in t1 t2 t3 {
        if `"``t''"' == "" continue
        local neg = substr(`"``t''"', 1, 1) == "!"
        mata: st_local("n", strofreal(sum(strpos(cat(st_local("lg")), substr(st_local("`t'"), 1 + `neg', .)) :> 0)))
        if (`n' == 0) != `neg' {
            di as err `"`cmd': message text "' cond(`neg', "found", "not found") `": ``t''"'
            exit 9
        }
    }
end

* the dependent variable cannot be a raking variable, also as c., inside an
* interaction or a block, and before the c./contvars() rules; the by() variable is checked
* here for its message (below for its return codes)
u2_msg "seldecomp y x1 c.y, by(group) contvars(y 4)" 198 "the dependent variable y cannot be a raking variable"
rcof "noi seldecomp y c.y, by(group)" == 198
rcof "noi seldecomp y x1#y, by(group)" == 198
rcof "noi seldecomp y (a: x2 y) x1, by(group) sequentially fixothermargins" == 198
u2_msg "seldecomp y group x1, by(group)" 198 "the by() variable group cannot be a raking variable"
* every Focal Group observation has a level absent from the Reference Group
preserve
keep if mod(_n, 10) == 0
* (the message names the variable and gives no hint about c. or contvars())
u2_msg "seldecomp y id, by(group)" 2000 "every focal group observation has a level of id absent from the reference group" "!contvars(id #)"
restore

* errors inside bayesboot() and contvars() name the option list
u2_msg "seldecomp y x1, by(group) bayesboot(cluser(cluster))" 198 "bayesboot(): option cluser() not allowed"
u2_msg "seldecomp y x1, by(group) bayesboot(seed(1))" 198 "bayesboot(): option seed() not allowed; use rseed()"
u2_msg "seldecomp y x1, by(group) bayesboot(5)" 198 "bayesboot(5) not allowed; type bayesboot(reps(5))"
u2_msg "seldecomp y x1, by(group) bayesboot(50 rseed(1))" 198 "bayesboot(50) not allowed; type bayesboot(reps(50))"
* a known suboption with a bad argument keeps -syntax-'s message and code, and the list is named
u2_msg "seldecomp y x1, by(group) bayesboot(reps(abc))" 198 "option reps() incorrectly specified" "(error in option bayesboot())"
u2_msg "seldecomp y x1, by(group) bayesboot(cluster(nosuchvar))" 111 "(error in option bayesboot())"
* the raking options are options of their own
u2_msg "seldecomp y x1 x2, by(group) maxiter(5)" 198 "option maxiter() not allowed"
u2_msg "seldecomp y x1 x2, by(group) iterate(abc)" 198 "option iterate() incorrectly specified"
u2_msg "seldecomp y c.x3 x1, by(group) contvars(x3 4 pooled)" 198 "contvars(): put a comma before pooled, as in contvars(x3 4, pooled)"
rcof "noi seldecomp y c.x3 x1, by(group) contvars(x3 4 refgroup)" == 198
u2_msg "seldecomp y c.x3 x1, by(group) contvars(x3 4, poled)" 198 "contvars(): option poled not allowed"

* group-value errors list the values of by() in the estimation sample; more than two values
* do not point to if inlist()
u2_msg "seldecomp y x1, by(group) focalgroup(5)" 459 "focalgroup(5) is not a value of group; its values are 0 and 1"
u2_msg "seldecomp y x1, by(group) focalgroup(abc)" 198 "focalgroup(abc) must be a number; the values of group are 0 and 1"
preserve
* a numeric by() with a value label takes the value or the label text: a
* value first, else a label (case-sensitive); unknown or ambiguous text lists the
* values and their labels
label define u2gl 0 "Zero" 1 "One"
label values group u2gl
seldecomp y x1, by(group) focalgroup(1)
local e1 = e(explained)
foreach o in "focalgroup(One)" `"focalgroup("One")"' {
    seldecomp y x1, by(group) `o'
    assert "`e(focal_group)'" == "1" & "`e(focal_group_label)'" == "1 (One)" & e(explained) == `e1'
}
u2_msg "seldecomp y x1, by(group) focalgroup(one)" 459 "focalgroup(one) is neither a value nor a value label of group; its values are 0 (Zero) and 1 (One)"
u2_msg "seldecomp y x1, by(group) focalgroup(5)" 459 "focalgroup(5) is neither a value nor a value label of group; its values are 0 (Zero) and 1 (One)"
* the label of a value outside the estimation sample names no group
label define u2gl 2 "Two", add
u2_msg "seldecomp y x1, by(group) focalgroup(Two)" 459 "focalgroup(Two) is neither a value nor a value label of group"
* a label with a blank
label define u2gl 0 "Not one", modify
seldecomp y x1, by(group) focalgroup(Not one)
assert "`e(focal_group)'" == "0" & "`e(focal_group_label)'" == "0 (Not one)"
* ambiguous: a value and the other value's label, or the label of both values
label define u2amb 0 "1" 1 "0"
label values group u2amb
u2_msg "seldecomp y x1, by(group) focalgroup(1)" 198 "focalgroup(1) is ambiguous: a value of group and the value label of the other value; its values are 0 (1) and 1 (0)"
* 1.0 is the value 1 and no label
seldecomp y x1, by(group) focalgroup(1.0)
assert "`e(focal_group)'" == "1" & e(explained) == `e1'
label define u2same 0 "Same" 1 "Same"
label values group u2same
u2_msg "seldecomp y x1, by(group) focalgroup(Same)" 198 "focalgroup(Same) is ambiguous: the value label of both values; its values are 0 (Same) and 1 (Same)"
* a label that reads as a number but is not a value: taken as the label, with a note
label define u2num 0 "7" 1 "8"
label values group u2num
u2_msg "seldecomp y x1, by(group) focalgroup(8)" 0 "note: focalgroup(8): 8 is not a value of group;" "label of 1"
assert "`e(focal_group)'" == "1" & e(explained) == `e1'
label values group
* a string by() takes the value that sorts first as Focal Group (uppercase before lowercase)
gen str1 sg = cond(group, "a", "B")
seldecomp y x1, by(sg)
assert "`e(focal_group)'" == "B" & e(N_focal) == `N0'
u2_msg "seldecomp y x1, by(sg) focalgroup(b)" 459 `"focalgroup(b) is not a value of sg; its values are "B" and "a""'
gen byte g3u = group + (id <= 100)
u2_msg "seldecomp y x1, by(g3u)" 459 "by(g3u) takes 3 values (0, 1, and 2); it must take two" "!inlist("
u2_msg "seldecomp y x1, by(x4)" 459 `"by(x4) takes 4 values ("east", "north", "south", and "west"); it must take two"' "!inlist("
gen byte g9 = mod(id, 9)
u2_msg "seldecomp y x1, by(g9)" 459 "by(g9) takes 9 values (0, 1, 2, 3, 4, 5, ...); it must take two"
restore

* a binding trimhiabs(): the raked weights stay at or below the bound (in the units of
* the user's weights) and differ from the untrimmed ones. A trimming note only when
* trimming left a margin unmet (kept by force), none for a bound that binds with every
* margin met or that does not bind
seldecomp y x1 x2 [pw=wt], by(group) generate(u2rw0)
forvalues m = 1/`e(N_notes)' {
    assert strpos(`"`e(note`m')'"', "trimming") == 0
}
sum u2rw0 if group == 0, meanonly
local bound = 0.8 * r(max)
seldecomp y x1 x2 [pw=wt], by(group) trimhiabs(`bound') generate(u2rw1) force
sum u2rw1 if group == 0, meanonly
assert r(max) <= `bound' * (1 + 1e-6)
count if group == 0 & reldif(u2rw1, u2rw0) > 1e-6
assert r(N) > 0
local hit 0
forvalues m = 1/`e(N_notes)' {
    if strpos(`"`e(note`m')'"', "trimming left") local hit = `hit' + 1
}
assert `hit' == (e(badcontrols) > 0)
seldecomp y x1 x2 [pw=wt], by(group) trimhirel(1e6)
assert e(N_at_bound) == 0 & e(badcontrols) == 0 & e(N_notes) == 0
seldecomp y x1 x2 [pw=wt], by(group) trimhirel(1e6) trimfrequency(once)
assert e(N_notes) == 0
drop u2rw0 u2rw1
* the 430 of a trimming option that changed no weight: the message does not mention the
* bounds; under noisily Engine reports the comparison with the untrimmed raking (the margins
* of section 5 that cannot be met)
preserve
clear
set obs 400
gen byte grp = _n > 200
gen byte a = cond(grp, mod(_n, 10) < 3, mod(_n, 2))
gen byte b = cond(grp, mod(floor(_n / 3), 10) < 6, 1 - mod(_n, 2))
gen double yy = a + b + mod(_n, 7) / 7
u2_msg "seldecomp yy i.a i.b, by(grp) trimhirel(1e6)" 430 "specify force to keep the results" "!the trimming bounds"
u2_msg "seldecomp yy i.a i.b, by(grp) trimhirel(1e6) noisily" 430 "trimming changed the raked weights: no"
* bounds that changed the weights: the 430 names them
u2_msg "seldecomp yy i.a i.b, by(grp) trimhirel(1.1) trimlorel(.9) noisily" 430 "trimming changed the raked weights: yes" "; the trimming bounds are binding: relax them, or specify force to keep the results" "!raking stopped at iterate()"
u2_msg "seldecomp yy i.a i.b, by(grp) trimhirel(1.1) trimlorel(.9) force" 0 "!the trimming bounds are binding"
* trimming bounds that contradict: the stage's error, naming the bounds; in the replicates of
* bayesboot(), whose base weights are redrawn, they fail every one here, force or not (force
* keeps only a replicate that missed raking targets; the full sample's
* smallest base weight is .5, so trimloabs(.6) is just within trimhirel(1.25) times it)
restore
preserve
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
u2_msg "seldecomp y x1 x2 [pw=wt], by(group) trimhirel(1.25) trimloabs(.9)" 198 "stage 1: the trimming bounds contradict each other: trimloabs() is above trimhirel() times the base weight of" "focal group observations" "!loosen one of them"
u2_msg "seldecomp y x1 x2 [pw=wt], by(group) trimlorel(.8) trimhiabs(2)" 198 "trimlorel() times the base weight of" "is above trimhiabs()"
u2_msg "seldecomp y x1 x2 [pw=wt], by(group) trimhirel(1.25) trimloabs(.6) force bayesboot(reps(4) rseed(1) nodots)" 430 "bayesboot(): 0 of 4 replications succeeded; 2 or more are required. Failures: the trimming bounds contradicted each other"
restore


* ==============================================================================
* parsing and group logic
* ==============================================================================
di as txt _n "---- parsing and group logic"

* the command line is not macro-expanded a second time: a literal $name in an if
* expression, a group value or e(cmdline) survives ($ is written \$ in a do-file)
preserve
gen str10 cur = cond(mod(_n, 2), "a" + char(36) + "b", "zz")
count if cur != "a\$b"
local nok = r(N)
assert `nok' < _N
seldecomp y x1 if cur != "a\$b", by(group) generate(rwd)
assert e(N) == `nok'
mata: assert(strpos(st_global("e(cmdline)"), "a" + char(36) + "b") > 0)
mata: assert(strpos(st_global("rwd[seldecomp_cmd]"), "a" + char(36) + "b") > 0)
gen str10 sd = cond(group, "b", "a" + char(36) + "US")
seldecomp y x1, by(sd) focalgroup(a\$US)
mata: assert(st_global("e(focal_group)") == "a" + char(36) + "US")
assert e(N_focal) == `N0'
restore

* a variable listed twice in contvars() is an error
rcof "noi seldecomp y c.x3, by(group) contvars(x3 4 x3 5)" == 198
rcof "noi seldecomp y c.x3, by(group) contvars(x3 4 c.x3 4)" == 198

* malformed interactions and repeated terms are errors
foreach t in "x1#" "#x1" "x1###x2" "x1#x1" "i.x1##x1" "x1#x2 x2#x1" {
    rcof "noi seldecomp y `t', by(group)" == 198
}
rcof "noi seldecomp y c.x3#x1 x1#c.x3, by(group) contvars(x3 4)" == 198
rcof "noi seldecomp y (a: x1#x2) (b: x2#x1), by(group) sequentially" == 198

* parentheses or commas inside a term are errors
rcof "noi seldecomp y i.(x1 x2), by(group)" == 198
rcof "noi seldecomp y (a: x1, x2), by(group) sequentially fixothermargins" == 198
rcof "noi seldecomp y ((x1) x2), by(group) sequentially fixothermargins" == 198

* the by() variable cannot be a raking variable
rcof "noi seldecomp y group x1, by(group)" == 198
rcof "noi seldecomp y x1#group, by(group)" == 198

* a negative weight on an observation outside the estimation sample (here: a missing
* raking variable) is not an error, as in -regress-
preserve
gen double wn = wt
replace wn = -1 in 1
replace x1 = . in 1
seldecomp y x1 [pw=wn], by(group)
assert e(N) == _N - 1
replace x1 = 1 in 1
rcof "noi seldecomp y x1 [pw=wn], by(group)" == 402
restore
* the other ways out of the sample ignore the negative weight
* too: a missing outcome, by() value or bayesboot(cluster()) value, and another by-group under
* the by prefix; inside the sample it is 402 in every path
preserve
gen long obs = _n
sum obs if x2 == 1, meanonly
local i = r(min)
gen double wn = wt
replace wn = -1 in `i'
gen double yn = y
replace yn = . in `i'
seldecomp yn x1 [pw=wn], by(group)
assert e(N) == _N - 1
gen double gn = group
replace gn = . in `i'
seldecomp y x1 [pw=wn], by(gn)
assert e(N) == _N - 1
gen double cn = cluster
replace cn = . in `i'
seldecomp y x1 [pw=wn], by(group) bayesboot(reps(2) rseed(1) nodots cluster(cn))
assert e(N) == _N - 1
count if x2 == 0
local n0 = r(N)
ereturn clear
rcof "noi by x2, sort: seldecomp y x1 [pw=wn], by(group)" == 402
assert e(N) == `n0'                 // by-group x2 == 0 ran; by-group x2 == 1 holds the weight
sort obs
rcof "noi seldecomp y (x1) (x2) [pw=wn], by(group) sequentially fixothermargins" == 402
rcof "noi seldecomp y x1 [pw=wn], by(group) bayesboot(reps(2) rseed(1) nodots)" == 402
restore

* operators on the dependent variable give 101
rcof "noi seldecomp L.y x1, by(group)" == 101
rcof "noi seldecomp i.y x1, by(group)" == 101
rcof "noi seldecomp y#x1 x2, by(group)" == 101
rcof "noi seldecomp y##x1 x2, by(group)" == 101

* generate() and bayesboot(generate()) naming the same variable are refused before
* estimating; no rw_2 is left in the data
rcof "noi seldecomp y x1 in 1/300, by(group) generate(rw_1) bayesboot(reps(2) rseed(1) nodots generate(rw))" == 198
capture confirm variable rw_2
assert _rc == 111

* a trimming bound of -1 (the internal "absent" marker) is refused like any other
* nonpositive bound
foreach o in trimhirel trimhiabs trimlorel trimloabs {
    rcof "noi seldecomp y x1 x2 [pw=wt], by(group) `o'(-1)" == 198
}


* ==============================================================================
* notes on excluded observations: the note on the observations excluded for a zero base
* weight and e(N_zerowt); the error when a group's observations all have zero weight; one
* message per tolerance
* ==============================================================================
di as txt _n "---- notes on excluded observations"

capture program drop u5_count
program define u5_count, rclass
    // u5_count "command" rc "text": the command fails with rc (0: runs); r(n) is the number
    // of lines of its output (at linesize 255) that contain text, a note or warning wrapped at
    // 78 columns counted as one line (its indented continuation lines joined to it)
    args cmd rc t
    tempfile lg
    local ls = c(linesize)
    set linesize 255
    log using `"`lg'"', text replace name(u5log)
    capture noisily `cmd'
    local r = _rc
    log close u5log
    set linesize `ls'
    if `r' != `rc' {
        di as err `"`cmd': rc `r', expected `rc'"'
        exit 9
    }
    mata: _u5_l = cat(st_local("lg"))
    mata: _u5_r = J(0, 1, "")
    mata: _u5_p = 0
    mata: for (_u5_i = 1; _u5_i <= rows(_u5_l); _u5_i++) { ///
        if (_u5_p & substr(_u5_l[_u5_i], 1, 4) == "    " & strtrim(_u5_l[_u5_i]) != "") ///
            _u5_r[rows(_u5_r)] = _u5_r[rows(_u5_r)] + " " + strtrim(_u5_l[_u5_i]); ///
        else { ///
            _u5_r = _u5_r \ _u5_l[_u5_i]; ///
            _u5_p = (substr(_u5_l[_u5_i], 1, 5) == "note:" | substr(_u5_l[_u5_i], 1, 8) == "warning:"); ///
        } ///
    }
    mata: st_local("n", strofreal(sum(strpos(_u5_r, st_local("t")) :> 0)))
    mata: mata drop _u5_l _u5_r _u5_p _u5_i
    return scalar n = `n'
end

import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
recast double y x3 wt

* the zero-weight note: its text, the counts by group, e(N_zerowt), shown at estimation and replay
gen double wz = wt
replace wz = 0 in 1/9
count if wz == 0 & group == 0
local z1 = r(N)
count if wz == 0 & group == 1
local z2 = r(N)
assert `z1' > 0 & `z2' > 0 & `z1' != `z2'      // unequal: a swap of the groups shows
local ztxt "9 observations with zero weight excluded (`z1' in the focal group, `z2' in the reference group)"
u5_count "seldecomp y x1 [pw=wz], by(group)" 0 "note: `ztxt'"
assert r(n) == 1
assert e(N_zerowt) == 9 & e(N) == 1991 & e(N_focal) == `N0' - `z1' & e(N_ref) == `N1' - `z2'
assert e(N_notes) == 1 & `"`e(note1)'"' == "`ztxt'"
u5_count "seldecomp" 0 "note: `ztxt'"
assert r(n) == 1
* the other group raked: the counts follow the groups
seldecomp y x1 [pw=wz], by(group) focalgroup(1)
assert `"`e(note1)'"' == "9 observations with zero weight excluded (`z2' in the focal group, `z1' in the reference group)"
* the same sample with the zero-weight observations dropped by if: equal results, no note
seldecomp y x1 [pw=wz], by(group)
matrix u5b = e(b)
seldecomp y x1 [pw=wt] if wz > 0, by(group)
assert mreldif(e(b), u5b) < `itol'
assert e(N_zerowt) == 0 & e(N_notes) == 0
matrix drop u5b
* one observation
gen double wz1 = wt
replace wz1 = 0 in 1
seldecomp y x1 [pw=wz1], by(group)
assert e(N_zerowt) == 1
assert `"`e(note1)'"' == "1 observation with zero weight excluded (1 in the focal group, 0 in the reference group)"

* no note and e(N_zerowt) = 0 without zero weights, weighted or not
seldecomp y x1 [pw=wt], by(group)
assert e(N_zerowt) == 0 & e(N_notes) == 0
seldecomp y x1, by(group)
assert e(N_zerowt) == 0 & e(N_notes) == 0

* only the observations that the weight alone excludes are counted: a missing outcome,
* outcomes() or ratio() variable, raking variable or group takes an observation out
* first; aweights and iweights are counted as pweights
gen double yz = y
replace yz = . in 1/2
gen double ybz = ybin
replace ybz = . in 3
gen double den = 1 + x2
replace den = . in 4
gen double x1z = x1
replace x1z = . in 5
gen double gz = group
replace gz = . in 6
count if wz == 0 & group == 0 & !inrange(_n, 1, 6)
local zz1 = r(N)
seldecomp yz x1z [pw=wz], by(gz) outcomes(ybz) ratio(den)
assert e(N_zerowt) == 3
* (the two observations missing ybz or den are noted first)
count if group == 0 & inlist(_n, 3, 4)
assert `"`e(note1)'"' == "2 observations with a missing value of an outcomes() or the ratio() variable excluded (`r(N)' in the focal group, `=2 - r(N)' in the reference group); every outcome uses the same observations"
assert `"`e(note2)'"' == "3 observations with zero weight excluded (`zz1' in the focal group, `=3 - `zz1'' in the reference group)"
seldecomp y x1 [aw=wz], by(group)
assert e(N_zerowt) == 9 & `"`e(note1)'"' == "`ztxt'"
seldecomp y x1 [iw=wz], by(group)
assert e(N_zerowt) == 9 & `"`e(note1)'"' == "`ztxt'"
drop yz ybz den x1z gz

* insupport: the zero-weight observations leave first, the screen runs on the rest; the
* zero-weight note comes before the insupport note; e(N_focal_all) excludes them
gen byte x1s = x1
replace x1s = 9 if group == 0 & inrange(_n, 11, 60)
count if x1s == 9
local n9 = r(N)
seldecomp y x1s [pw=wz], by(group) insupport
assert e(N_zerowt) == 9 & `"`e(note1)'"' == "`ztxt'"
assert `"`e(note2)'"' == "insupport excluded x1s = 9 (`n9' obs)"
assert e(N_focal_all) == `N0' - `z1' & e(N_focal) == `N0' - `z1' - `n9'
drop x1s

* string by() values, compared as data ($ in a value), and a third by() value whose
* observations all have zero weight (in neither group)
gen str10 sg = cond(group, "b", "a" + char(36) + "US")
seldecomp y x1 [pw=wz], by(sg)
assert e(N_zerowt) == 9 & `"`e(note1)'"' == "`ztxt'"
gen byte g3 = group
replace g3 = 2 in 1/4
count if wz == 0 & g3 == 0
local y1 = r(N)
count if wz == 0 & g3 == 1
local y2 = r(N)
seldecomp y x1 [pw=wz], by(g3)
assert e(N_zerowt) == 9
assert `"`e(note1)'"' == "9 observations with zero weight excluded (`y1' in the focal group, `y2' in the reference group, 4 with other values of g3)"

* the by prefix: a note and a count for each by-group
count if wz == 0 & x2 == 0
local zx0 = r(N)
count if wz == 0 & x2 == 1
local zx1 = r(N)
assert `zx0' > 0 & `zx1' > 0
preserve
u5_count "by x2, sort: seldecomp y x1 [pw=wz], by(group)" 0 "with zero weight excluded"
assert r(n) == 2
assert e(N_zerowt) == `zx1'                 // the last by-group, x2 == 1
restore

* bayesboot(): the note belongs to the point estimate (once, as without bayesboot())
u5_count "seldecomp y x1 [pw=wz], by(group) bayesboot(reps(5) rseed(1) nodots)" 0 "with zero weight excluded"
assert r(n) == 1
assert e(N_zerowt) == 9 & e(N_notes) == 1 & `"`e(note1)'"' == "`ztxt'"
assert e(N) == 1991 & e(N_reps) == 5

* the replication prefixes: the replicate weights of svy brr are 0 for the half-sample
* left out, so every replicate call excludes half of the observations; no replicate
* call makes the note (they are detected by _loop in c(prefix)), and the full-sample
* call has no zero weights. Called without the prefix, the same weights give the note
preserve
clear
set seed 20260926
set obs 384
gen int stratum = ceil(_n / 24)
gen byte psu = 1 + mod(ceil(_n / 12) - 1, 2)
gen int psuid = stratum * 10 + psu
gen double wt = 50 + 100 * runiform()
gen byte group = runiform() < .5
gen byte x1 = 1 + (runiform() < .4) + (runiform() < .3 + .3 * group)
gen double y = 2 + x1 + .8 * group + rnormal()
mata: Hd = (1); for (i = 1; i <= 4; i++) Hd = (Hd, Hd \ Hd, -Hd); st_matrix("U5HAD", Hd)
mata: mata drop Hd i
forvalues r = 1/16 {
    gen double brr`r' = 2 * wt * ((U5HAD[`r', stratum] == 1) == (psu == 1))
}
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr16)
u5_count "svy brr, noisily: seldecomp y i.x1, by(group)" 0 "with zero weight excluded"
assert r(n) == 0
assert e(N_zerowt) == 0 & e(N_notes) == 0 & e(N_reps) == 16
u5_count "svy brr, noisily: seldecomp y i.x1, by(group)" 0 "Number of obs     =        192"
assert r(n) == 16                            // the 16 replicate calls, each on a half-sample
seldecomp y i.x1 [pw=brr1], by(group)
assert e(N_zerowt) == 192 & e(N) == 192
assert strpos(`"`e(note1)'"', "192 observations with zero weight excluded (") == 1
* zero base weights (7): the full-sample call of svy brr, svy sdr and jackknife: receives
* them and makes the note; svy jackknife leaves them out before it calls seldecomp, so there
* is no note and e(N_zerowt) is 0, and svy's e(N), which counts them, exceeds N_focal + N_ref
* (help, weights); the same estimates
replace wt = 0 in 1/7
forvalues r = 1/16 {
    replace brr`r' = 0 in 1/7
}
count if wt == 0 & group == 0
local s1 = r(N)
local stxt "7 observations with zero weight excluded (`s1' in the focal group, `=7 - `s1'' in the reference group)"
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr16) sdrweight(brr1-brr16)
svy brr: seldecomp y i.x1, by(group)
assert e(N_zerowt) == 7 & e(N_notes) == 1 & `"`e(note1)'"' == "`stxt'"
assert e(N) == 384 & e(N_focal) + e(N_ref) == 377
matrix u5b = e(b)
svy sdr: seldecomp y i.x1, by(group)
assert e(N_zerowt) == 7 & `"`e(note1)'"' == "`stxt'"
jackknife, cluster(psuid): seldecomp y i.x1 [pw=wt], by(group)
assert e(N_zerowt) == 7 & `"`e(note1)'"' == "`stxt'"
svy jackknife: seldecomp y i.x1, by(group)
assert e(N_zerowt) == 0 & e(N_notes) == 0
assert e(N) == 384 & e(N_focal) + e(N_ref) == 377
assert mreldif(e(b), u5b) < `itol'
u5_count "seldecomp" 0 "with zero weight excluded"
assert r(n) == 0
matrix drop U5HAD u5b
restore

* a group whose observations all have zero weight: the error says so (459), naming the
* group as focalgroup() (a value or a value label) makes it, else "the other group";
* with every weight zero, 2000
gen double w0 = wt * (group == 1)
gen double w1 = wt * (group == 0)
u2_msg "seldecomp y x1 [pw=w0], by(group)" 459 "by(group) takes one value (1): all `N0' focal group observations (group = 0) with nonmissing data have zero weight"
u2_msg "seldecomp y x1 [pw=w0], by(group) focalgroup(1)" 459 "all `N0' reference group observations (group = 0) with nonmissing data have zero weight"
u2_msg "seldecomp y x1 [pw=w0], by(group) focalgroup(7)" 459 "all `N0' observations of the other group (group = 0) with nonmissing data have zero weight"
u2_msg "seldecomp y x1 [pw=w1], by(group)" 459 "by(group) takes one value (0): all `N1' reference group observations (group = 1) with nonmissing data have zero weight"
u2_msg "seldecomp y x1 [pw=w1], by(group) focalgroup(0)" 459 "all `N1' reference group observations (group = 1) with nonmissing data have zero weight"
label define u5gl 0 "Zero" 1 "One"
label values group u5gl
u2_msg "seldecomp y x1 [pw=w0], by(group) focalgroup(One)" 459 "all `N0' reference group observations (group = 0 (Zero)) with nonmissing data have zero weight"
u2_msg "seldecomp y x1 [pw=w0], by(group) focalgroup(Zero)" 459 "all `N0' focal group observations (group = 0 (Zero)) with nonmissing data have zero weight"
label values group
gen str1 sb = cond(group, "b", "a")
u2_msg "seldecomp y x1 [pw=w0], by(sb)" 459 `"by(sb) takes one value ("b"): all `N0' focal group observations (sb = "a") with nonmissing data have zero weight"'
u2_msg "seldecomp y x1 [pw=w1], by(sb) focalgroup(b)" 459 `"all `N1' focal group observations (sb = "b") with nonmissing data have zero weight"'
u2_msg "seldecomp y x1 [pw=w1], by(sb) focalgroup(a)" 459 `"all `N1' reference group observations (sb = "b") with nonmissing data have zero weight"'
* one observation, and several other values
gen long obs = _n
sum obs if group == 0, meanonly
local o1 = r(min)
u2_msg "seldecomp y x1 [pw=w0] if group == 1 | _n == `o1', by(group)" 459 "the one focal group observation (group = 0) with nonmissing data has zero weight"
replace g3 = 2 in 1/4
gen double w3 = wt * (g3 == 1)
u2_msg "seldecomp y x1 [pw=w3], by(g3)" 459 "the observations with the other values of g3 (0 and 2) and nonmissing data all have zero weight"
* the general message when the other group is missing for another reason
u2_msg "seldecomp y x1 [pw=wz] if group == 1, by(group)" 459 "the other group has no observations with nonmissing data and a positive weight"
* every weight zero
gen double wnil = 0
u2_msg "seldecomp y x1 [pw=wnil], by(group)" 2000 "no observations: all 2000 observations with nonmissing data have zero weight"
* the count is of the observations with nonmissing data: Focal Group observations with a
* positive weight and a missing outcome are not among them (the text says "with nonmissing
* data")
gen double ym = y
replace ym = . if group == 0 & _n <= 100
count if group == 0 & !missing(ym)
local nm0 = r(N)
gen double wm = wt
replace wm = 0 if group == 0 & !missing(ym)
u5_count "seldecomp ym x1 [pw=wm], by(group)" 459 "all `nm0' focal group observations (group = 0) with nonmissing data have zero weight"
assert r(n) == 1
u5_count "seldecomp ym x1 [pw=wm], by(group)" 459 "observations (group = 0) have zero weight"
assert r(n) == 0
drop ym wm
* by() with three values: focalgroup() names the value whose observations
* all have zero weight, and two values are left; the 459 says so: a value, a value label,
* a string value; one observation; with no zero weights, or naming another value, the
* message says only that the value is not in the estimation sample
replace g3 = group
replace g3 = 2 in 1/300
gen double w32 = wt * (g3 != 2)
local z32 "all 300 observations with g3 = 2 and nonmissing data have zero weight"
u2_msg "seldecomp y x1 [pw=w32], by(g3) focalgroup(2)" 459 "focalgroup(2) is not a value of g3: `z32'; its values are 0 and 1"
u2_msg "seldecomp y x1 [pw=w32], by(g3) focalgroup(2.00000000001)" 459 "`z32'"
label define u5g3 0 "Zero" 1 "One" 2 "Two"
label values g3 u5g3
u2_msg "seldecomp y x1 [pw=w32], by(g3) focalgroup(Two)" 459 "focalgroup(Two) is neither a value nor a value label of g3: all 300 observations with g3 = 2 (Two) and nonmissing data have zero weight; its values are 0 (Zero) and 1 (One)"
u2_msg "seldecomp y x1 [pw=w32], by(g3) focalgroup(2)" 459 "focalgroup(2) is neither a value nor a value label of g3: all 300 observations with g3 = 2 (Two) and nonmissing data have zero weight"
label values g3
gen str1 s3 = cond(g3 == 2, "c", cond(group, "b", "a"))
u2_msg "seldecomp y x1 [pw=w32], by(s3) focalgroup(c)" 459 `"focalgroup(c) is not a value of s3: all 300 observations with s3 = "c" and nonmissing data have zero weight; its values are "a" and "b""'
u2_msg "seldecomp y x1 [pw=w32] if g3 < 2 | _n == 1, by(g3) focalgroup(2)" 459 "focalgroup(2) is not a value of g3: the one observation with g3 = 2 and nonmissing data has zero weight; its values are 0 and 1"
u5_count "seldecomp y x1 [pw=w32], by(g3) focalgroup(3)" 459 "focalgroup(3) is not a value of g3; its values are 0 and 1"
assert r(n) == 1
u5_count "seldecomp y x1 [pw=wt] if g3 < 2, by(g3) focalgroup(2)" 459 "focalgroup(2) is not a value of g3; its values are 0 and 1"
assert r(n) == 1
u5_count "seldecomp y x1 if g3 < 2, by(s3) focalgroup(c)" 459 `"focalgroup(c) is not a value of s3; its values are "a" and "b""'
assert r(n) == 1
drop w32 s3
label drop u5g3
drop w0 w1 w3 wnil obs sg sb g3 wz wz1

* one message per tolerance, naming the option that is wrong (rc 198)
u2_msg "seldecomp y x1 x2, by(group) tolerance(0)" 198 "tolerance() must be a positive number"
u2_msg "seldecomp y x1 x2, by(group) ctrltolerance(-1e-6)" 198 "ctrltolerance() must be a positive number"
u5_count "seldecomp y x1 x2, by(group) ctrltolerance(0)" 198 " tolerance()"
assert r(n) == 0
u5_count "seldecomp y x1 x2, by(group) tolerance(-1)" 198 "ctrltolerance"
assert r(n) == 0
u2_msg "seldecomp y x1 x2, by(group) iterate(0)" 198 "iterate() must be a positive integer"


* ==============================================================================
* messages for common slips at the command line, as asserts
* ==============================================================================
di as txt _n "---- messages for common command-line slips"

capture program drop s2_msg
program define s2_msg, rclass
    // s2_msg "command" rc ["text" ...]: the command fails with rc (0: runs), and its output,
    // its lines stripped and joined by one blank (a wrapped paragraph matches), contains each
    // text; a text that begins with ! must not occur. r(first): the first nonblank line;
    // r(all): the output so joined
    gettoken cmd 0 : 0
    gettoken rc 0 : 0
    tempfile lg
    log using `"`lg'"', text replace name(s2msg)
    capture noisily `cmd'
    local r = _rc
    log close s2msg
    if `r' != `rc' {
        di as err `"`cmd': rc `r', expected `rc'"'
        exit 9
    }
    // the output: the nonblank lines after the log header (which ends with "opened on:")
    mata: _s2_l = strtrim(cat(st_local("lg"))); _s2_n = select(1::rows(_s2_l), strpos(_s2_l, "opened on:") :== 1)
    mata: _s2_n = (rows(_s2_n) ? _s2_n[1] : 0)
    mata: _s2_l = (_s2_n < rows(_s2_l) ? _s2_l[(_s2_n + 1)::rows(_s2_l)] : J(0, 1, ""))
    // (and before the log's closing lines, which begin with "name:")
    mata: _s2_n = select(1::rows(_s2_l), strpos(_s2_l, "name:") :== 1); _s2_n = (rows(_s2_n) ? _s2_n[rows(_s2_n)] : 0)
    mata: _s2_l = (_s2_n > 1 ? _s2_l[1::(_s2_n - 1)] : (_s2_n == 1 ? J(0, 1, "") : _s2_l))
    mata: _s2_l = select(_s2_l, _s2_l :!= "")
    mata: st_local("first", rows(_s2_l) ? _s2_l[1] : ""); _s2_a = stritrim(invtokens(_s2_l')); _s2_t = ""
    while `"`0'"' != "" {
        gettoken t 0 : 0
        mata: _s2_t = st_local("t"); _s2_n = (substr(_s2_t, 1, 1) == "!")
        mata: st_local("bad", strofreal((strpos(_s2_a, substr(_s2_t, 1 + _s2_n, .)) > 0) == _s2_n))
        if `bad' {
            di as err `"`cmd': message text "' cond(substr(`"`t'"', 1, 1) == "!", "found", "not found") `": `t'"'
            exit 9
        }
    }
    mata: st_local("all", _s2_a)
    mata: mata drop _s2_l _s2_a _s2_t _s2_n
    return local first `"`first'"'
    return local all `"`all'"'
end

preserve
sysuse nlsw88, clear
local base "wage i.collgrad i.south, by(married)"

* fweights are refused (101), with the remedy and without the reason; the expression typed
* is echoed in expand
gen int fw = 1 + mod(_n, 3)
s2_msg "seldecomp wage i.collgrad i.south [fw=fw], by(married)" 101 "fweights not allowed" ///
    "expand the data first: expand fw" "!overstates the standard errors"
s2_msg "seldecomp wage i.collgrad [fweight = 2*fw], by(married)" 101 "expand the data first: expand 2*fw"

* bootstrap and raking suboptions, and oaxaca's options, typed at the top level: the
* seldecomp equivalent
s2_msg "seldecomp `base' bayesboot(reps(20)) rseed(1)" 198 "rseed() is a suboption of bayesboot(): type bayesboot(reps(20) rseed(1))"
s2_msg "seldecomp `base' bayesboot reps(20)" 198 "reps() is a suboption of bayesboot(): type bayesboot(reps(20))"
s2_msg "seldecomp `base' bayesboot(reps(20)) cluster(industry)" 198 "cluster() is a suboption of bayesboot(): type bayesboot(reps(20) cluster(industry))"
s2_msg "seldecomp `base' bayesboot(reps(20)) nodots" 198 "nodots is a suboption of bayesboot(): type bayesboot(reps(20) nodots)"
* level() and citable are suboptions of bayesboot() (a replication prefix takes its own level()),
* replace of generate()
s2_msg "seldecomp `base' bayesboot(reps(20)) level(90)" 198 "level() is a suboption of bayesboot(): type bayesboot(reps(20) level(90)); with a replication prefix, give level() to the prefix"
s2_msg "seldecomp `base' bayesboot(reps(20)) citable" 198 "citable is a suboption of bayesboot(): type bayesboot(reps(20) citable)"
s2_msg "seldecomp `base' generate(rw) replace" 198 "replace is a suboption of generate(): type generate(rw, replace)"
s2_msg "seldecomp `base' replace" 198 "replace is a suboption of generate(): type generate(newvar, replace)"
s2_msg "seldecomp `base' generate(rw, replace nodots)" 198 "generate(): type generate(newvar) or generate(newvar, replace)"
s2_msg "seldecomp `base' bayesboot(reps(2) generate(bw replace))" 198 "bayesboot(generate()): type generate(stub) or generate(stub, replace)"
s2_msg "seldecomp `base' bayesboot(reps(20)) seed(3) nodots" 198 "rseed() and nodots are suboptions of bayesboot(): type bayesboot(reps(20) rseed(3) nodots); use rseed(), not seed()"
s2_msg "seldecomp `base' rep(20) dots(5) cluster(race) saving(x) bca" 198 "reps(), dots(), cluster(), saving() and bca are suboptions of bayesboot(): type bayesboot(rep(20) dots(5) cluster(race) saving(x) bca)"
* the raking options, abbreviated
s2_msg "seldecomp `base' iter(50) tol(1e-9) ctrltol(1e-7)" 0
assert "`e(rakeopts)'" == "iterate(50) tolerance(1e-9) ctrltolerance(1e-7)"
* (without sequentially, detail says what it needs)
s2_msg "seldecomp `base' detail" 198 "detail requires sequentially" "!each term's" "!seldecomp's detailed decomposition"
s2_msg "seldecomp `base' swap" 198 "swap not allowed; to rake the other group, use focalgroup(1)"
s2_msg "seldecomp `base' weight(0)" 198 "weight(0) not allowed; to rake the other group, use focalgroup(1)"
s2_msg "seldecomp `base' weight(1)" 198 "weight(1) not allowed; it is the default"
s2_msg "seldecomp `base' weight(0.5)" 198 "weight() not allowed; the default corresponds to weight(1), and focalgroup(1) to weight(0)"
s2_msg "seldecomp `base' pooled" 198 "pooled not allowed; pooled is a suboption of contvars()"
s2_msg "seldecomp `base' vce(bootstrap, reps(50))" 198 "vce(bootstrap) not allowed; use bayesboot() or the bootstrap: prefix"
s2_msg "seldecomp `base' vce(cluster industry)" 198 "vce(cluster) not allowed; use bayesboot(cluster(industry))"
s2_msg "seldecomp `base' vce(jackknife)" 198 "vce(jackknife) not allowed; use the jackknife: prefix"
s2_msg "seldecomp `base' vce(robust)" 198 "vce() not allowed; use bayesboot() or a replication prefix"
gen str3 ms = cond(married, "yes", "no")
s2_msg "seldecomp wage i.collgrad, by(ms) swap" 198 "focalgroup(<the other group>)"
* the command shown is accepted: a suboption that takes an argument typed
* without one is shown with it (dots(1), the default; a placeholder otherwise), and one typed both
* inside and outside the option, or twice outside it, is shown once
s2_msg "seldecomp `base' dots" 198 "dots() is a suboption of bayesboot(): type bayesboot(dots(1))"
s2_msg "seldecomp `base' bayesboot(dots(1) reps(2) rseed(1))" 0
s2_msg "seldecomp `base' reps" 198 "reps() is a suboption of bayesboot(): type bayesboot(reps(#))"
s2_msg "seldecomp `base' seed" 198 "rseed() is a suboption of bayesboot(): type bayesboot(rseed(#)); use rseed(), not seed()"
s2_msg "seldecomp `base' bayesboot(reps(2) rseed(1)) rseed(2)" 198 ///
    "rseed() is a suboption of bayesboot(): type bayesboot(reps(2) rseed(1)); rseed(2) is left out, as bayesboot() has rseed(1) already"
s2_msg "seldecomp `base' bayesboot(reps(2) rseed(1))" 0
s2_msg "seldecomp `base' bayesboot(r(2)) reps(30) rseed(4)" 198 ///
    "reps() and rseed() are suboptions of bayesboot(): type bayesboot(r(2) rseed(4)); reps(30) is left out, as bayesboot() has r(2) already"
s2_msg "seldecomp `base' bayesboot(r(2) rseed(4))" 0
s2_msg "seldecomp `base' reps(20) reps(30)" 198 "type bayesboot(reps(20)); reps(30) is left out, as bayesboot() has reps(20) already"
s2_msg "seldecomp `base' bayesboot(reps(2) rseed(1)) reps(3) rseed(2)" 198 ///
    "type bayesboot(reps(2) rseed(1)); reps(3) and rseed(2) are left out, as bayesboot() has reps(2) and rseed(1) already"
s2_msg "seldecomp `base' bayesboot(reps(2) rseed(1)) seed(2) dots" 198 ///
    "rseed() and dots() are suboptions of bayesboot(): type bayesboot(reps(2) rseed(1) dots(1)); seed(2) is left out, as bayesboot() has rseed(1) already"
* swap and weight(0) with focalgroup() typed: focalgroup() with the other
* value, in place of the one typed
s2_msg "seldecomp `base' focalgroup(1) swap" 198 "swap not allowed; to rake the other group, use focalgroup(0) in place of focalgroup(1)"
s2_msg "seldecomp `base' focalgroup(1) weight(0)" 198 "weight(0) not allowed; to rake the other group, use focalgroup(0) in place of focalgroup(1)"
s2_msg "seldecomp `base' focalgroup(1) weight(0.5)" 198 "and focalgroup(0) in place of focalgroup(1) to weight(0)"
s2_msg "seldecomp wage i.collgrad, by(ms) focalgroup(yes) swap" 198 "focalgroup(no) in place of focalgroup(yes)"
* the swap makes the other group Focal Group: focalgroup(1) swapped is focalgroup(0), Focal Group = 0
s2_msg "seldecomp `base' focalgroup(0)" 0
assert e(focal_group) == "0"
* anything else keeps -syntax-'s own message: an unknown word, a known option with a bad
* argument, a missing by()
s2_msg "seldecomp `base' noisy" 198 "option noisy not allowed" "!suboption"
s2_msg "seldecomp `base' mincell(x)" 198 "option mincell() incorrectly specified"
s2_msg "seldecomp `base' bayesboot(reps(2) level(123))" 198 "level() must be between 10 and 99.99 inclusive"
s2_msg "seldecomp wage i.collgrad i.south, over(married)" 198 "option by() required"

* an option typed before the comma: "option by() must follow the comma:" and the command,
* also after an if or in qualifier and under by:
s2_msg "seldecomp wage i.collgrad i.south by(married)" 198 "option by() must follow the comma: seldecomp wage i.collgrad i.south, by(married)"
s2_msg "seldecomp wage i.collgrad i.south by(married) focalgroup(1), sequentially fixothermargins" 198 "seldecomp wage i.collgrad i.south, by(married) focalgroup(1) sequentially fixothermargins"
s2_msg "seldecomp wage (edu: i.collgrad) i.south focalgroup(1), by(married)" 198 "option focalgroup() must follow the comma: seldecomp wage (edu: i.collgrad) i.south, focalgroup(1) by(married)"
s2_msg "seldecomp wage i.collgrad [pw=hours] if inlist(race, 1, 2) by(married), sequentially fixothermargins" 198 "option by() must follow the comma: seldecomp wage i.collgrad [pw=hours] if inlist(race, 1, 2), by(married) sequentially fixothermargins"
s2_msg "seldecomp wage i.collgrad in 1/1000 by(married)" 198 "option by() must follow the comma: seldecomp wage i.collgrad in 1/1000, by(married)"
s2_msg "by south, sort: seldecomp wage i.collgrad by(married)" 198 "option by() must follow the comma: seldecomp wage i.collgrad, by(married)" "!seldby"
s2_msg "by south, sort: seldecomp wage i.collgrad" 198 "option by() required" "!comma"
* the words after the option typed before the comma, and seldecomp's options
* without argument (sequentially, fixothermargins, force), go after the comma too, and the
* command shown runs; a variable stays a raking variable wherever it is typed
s2_msg "seldecomp wage i.collgrad by(married) sequentially fixothermargins" 198 "option by() must follow the comma: seldecomp wage i.collgrad, by(married) sequentially fixothermargins"
s2_msg "seldecomp wage i.collgrad, by(married) sequentially fixothermargins" 0
s2_msg "seldecomp wage i.collgrad sequentially fixothermargins by(married), force" 198 "option by() must follow the comma: seldecomp wage i.collgrad, sequentially fixothermargins by(married) force"
s2_msg "seldecomp wage i.collgrad, sequentially fixothermargins by(married) force" 0
s2_msg "seldecomp wage i.collgrad by(married) force" 198 "option by() must follow the comma: seldecomp wage i.collgrad, by(married) force"
s2_msg "seldecomp wage i.collgrad, by(married) force" 0
s2_msg "seldecomp wage i.collgrad i.south by(married) focalgroup(1) sequentially fixothermargins" 198 "option by() must follow the comma: seldecomp wage i.collgrad i.south, by(married) focalgroup(1) sequentially fixothermargins"
s2_msg "seldecomp wage i.collgrad i.south, by(married) focalgroup(1) sequentially fixothermargins" 0
s2_msg "by south, sort: seldecomp wage i.collgrad by(married) sequentially fixothermargins" 198 "option by() must follow the comma: seldecomp wage i.collgrad, by(married) sequentially fixothermargins" "!seldby"
s2_msg "by south, sort: seldecomp wage i.collgrad, by(married) sequentially fixothermargins" 0
s2_msg "seldecomp wage i.collgrad by(married) south" 198 "option by() must follow the comma: seldecomp wage i.collgrad south, by(married)"
s2_msg "seldecomp wage by(married) i.collgrad" 198 "option by() must follow the comma: seldecomp wage i.collgrad, by(married)"
s2_msg "seldecomp wage i.collgrad by(married) detail" 198 "option by() must follow the comma: seldecomp wage i.collgrad, by(married) detail"
s2_msg "seldecomp wage i.collgrad, by(married) detail" 198 "detail requires sequentially"
* no option typed before the comma: seldecomp's options without argument
* that are no variables, typed before it, get "... must follow the comma:" with
* the command, which runs; abbreviations as typed.
* by() after an if or in qualifier takes the option words before the qualifier along.
* The outcome is never moved; a word that is no option, a variable named like an
* option, and a command without a comma do not get that message
s2_msg "seldecomp wage i.collgrad sequentially, by(married)" 198 "sequentially must follow the comma: seldecomp wage i.collgrad, by(married) sequentially" "!not found"
s2_msg "seldecomp wage i.collgrad, by(married) sequentially" 0
s2_msg "seldecomp wage i.collgrad sequentially fixothermargins, by(married)" 198 "sequentially and fixothermargins must follow the comma: seldecomp wage i.collgrad, by(married) sequentially fixothermargins" "!not found"
s2_msg "seldecomp wage i.collgrad, by(married) sequentially fixothermargins" 0
s2_msg "seldecomp wage i.collgrad seq fix i.south force, by(married) focalgroup(1)" 198 "seq, fix and force must follow the comma: seldecomp wage i.collgrad i.south, by(married) focalgroup(1) seq fix force"
s2_msg "seldecomp wage i.collgrad i.south, by(married) focalgroup(1) seq fix force" 0
s2_msg "by south, sort: seldecomp wage i.collgrad sequentially fixothermargins, by(married)" 198 "sequentially and fixothermargins must follow the comma: seldecomp wage i.collgrad, by(married) sequentially fixothermargins" "!seldby"
s2_msg "seldecomp wage i.collgrad sequentially fixothermargins if south by(married)" 198 "option by() must follow the comma (and sequentially and fixothermargins): seldecomp wage i.collgrad if south, by(married) sequentially fixothermargins"
s2_msg "seldecomp wage i.collgrad if south, by(married) sequentially fixothermargins" 0
s2_msg "seldecomp wage i.collgrad force in 1/1000 by(married), sequentially fixothermargins" 198 "option by() must follow the comma (and force): seldecomp wage i.collgrad in 1/1000, by(married) sequentially fixothermargins force"
s2_msg "seldecomp wage i.collgrad in 1/1000, by(married) sequentially fixothermargins force" 0
s2_msg "seldecomp wage i.collgrad sequentially fixothermargins [pw=hours] if inlist(race, 1, 2) by(married)" 198 "option by() must follow the comma (and sequentially and fixothermargins): seldecomp wage i.collgrad if inlist(race, 1, 2) [pweight= hours], by(married) sequentially fixothermargins"
s2_msg "seldecomp wage i.collgrad if inlist(race, 1, 2) [pweight= hours], by(married) sequentially fixothermargins" 0
s2_msg "seldecomp sequentially i.collgrad, by(married)" 111 "dependent variable sequentially not found" "!must follow the comma"
s2_msg "seldecomp wage i.collgrad sequentialy, by(married)" 111 "variable sequentialy not found" "!must follow the comma"
s2_msg "seldecomp wage i.collgrad sequentially fixothermargins" 198 "option by() required" "!must follow the comma"
* sequential is an abbreviation of sequentially (the other blocks free); typed before the
* comma it is moved after it like any option
s2_msg "seldecomp wage i.collgrad, by(married) sequential" 0
assert "`e(stagetype)'" == "sequentially" & "`e(fixothermargins)'" == ""
s2_msg "seldecomp wage i.collgrad sequential, by(married)" 198 "sequential must follow the comma" "!not found"
gen byte force = mod(_n, 2)
s2_msg "seldecomp wage i.collgrad force, by(married)" 0
assert "`e(rakevars)'" == "collgrad force"
drop force

* by() that is a variable of the by prefix
s2_msg "bysort married: seldecomp wage i.collgrad, by(married)" 190 "married may not be both in by() and in the by prefix" "!nonmissing"
s2_msg "bysort race married: seldecomp wage i.collgrad, by(married)" 190 "married may not be both in by() and in the by prefix"
s2_msg "bysort south: seldecomp wage i.collgrad, by(married)" 0

* bayesboot() refuses the words varbyvar and sequential (198)
s2_msg "seldecomp `base' sequentially bayesboot(varbyvar reps(20))" 198 "bayesboot(): option varbyvar not allowed"
s2_msg "seldecomp `base' sequentially bayesboot(sequential reps(20))" 198 "bayesboot(): option sequential not allowed"

* the staged modes: fixothermargins needs sequentially or individually, the two
* modes may not be combined, and the option words in the lists below are refused
* (198, -syntax-'s own message)
s2_msg "seldecomp `base' fixothermargins" 198 "fixothermargins requires sequentially or individually"
s2_msg "seldecomp `base' sequentially individually" 198 "sequentially and individually may not be combined"
s2_msg "seldecomp `base' individually fixothermargins" 0
assert "`e(stagetype)'" == "individually" & "`e(fixothermargins)'" == "fixothermargins"
s2_msg "seldecomp `base' individually" 0
assert "`e(stagetype)'" == "individually" & "`e(fixothermargins)'" == ""
foreach o in varbyvar cumulative isolated nohold holdothers refminusfocal noheader notable nodefinitions {
    s2_msg "seldecomp `base' sequentially `o'" 198 "option `o' not allowed" "!must follow the comma"
}
foreach o in "refgroup(1)" "raker(mata)" "raker(ipfraking)" "ipfraking(iterate(50))" {
    mata: st_local("on", ustrregexra(st_local("o"), "[(].*$", "()"))
    s2_msg "seldecomp `base' `o'" 198 "option `on' not allowed" "!must follow the comma"
}
* the ipfraking options nograph and double are not raking options
s2_msg "seldecomp `base' nograph" 198 "option nograph not allowed"
s2_msg "seldecomp `base' double" 198 "option double not allowed"

* a continuous variable raked as categorical: the support failure, with no hint about c. or
* contvars(); the lists cut at 10 levels; no focalgroup() advice
s2_msg "seldecomp wage collgrad union ttl_exp, by(married)" 459 ///
    "collapse or drop these levels" "!if it is continuous" "!contvars(ttl_exp" "!focalgroup("
assert strpos(`"`r(first)'"', "raking support failure: these levels of the reference group") == 1
u5_count "seldecomp wage collgrad union ttl_exp, by(married)" 459 "ttl_exp = "
assert r(n) == 20
* a c. variable cut too finely: the support failure, with no hint about fewer quantile groups
s2_msg "seldecomp wage c.ttl_exp i.south, by(married) contvars(ttl_exp 200)" 459 ///
    "collapse or drop these levels" "!quantile groups"

* focalgroup() suggested only when raking the other group has support, nothing said when it
* has none; the set-aside count, "1 ... observation"
local spec "wage collgrad i.occupation i.industry union c.ttl_exp, by(married) contvars(ttl_exp 4)"
s2_msg "seldecomp `spec'" 459 "1 focal group observation set aside (levels absent from the reference group):" ///
    "collapse or drop these levels" "!or rake the other group" "!would fail as well"
s2_msg "seldecomp `spec' focalgroup(1)" 459 "collapse or drop these levels" "!or rake the other group" "!would fail as well"
* one direction fails, the other works: suggested, and it does work
gen byte occ9 = cond(occupation == 12, 13, occupation)
s2_msg "seldecomp wage i.occ9, by(married)" 459 "collapse or drop these levels, or rake the other group with focalgroup(1)"
s2_msg "seldecomp wage i.occ9, by(married) focalgroup(1)" 0
* the other group as the user chose the groups: focalgroup() typed, the
* advice names the other value "in place of" the one typed
s2_msg "seldecomp wage i.occ9, by(married) focalgroup(0)" 459 "or rake the other group with focalgroup(1) in place of focalgroup(0)" "!refgroup("
* set aside for a level of another variable
gen byte z = cond(married == 0 & _n <= 40, 7, south)
count if z == 7 & !missing(wage, occ9)
local nz = r(N)
s2_msg "seldecomp wage i.occ9 i.z, by(married)" 459 "`nz' focal group observations set aside (levels absent from the reference group):" "z = 7 (`nz' obs)" "!or rake the other group"

* "every focal group observation has a level ... absent from the reference group" names
* the variable, without asking whether it duplicates by(); several variables: each with its
* count
gen byte sep = married
s2_msg "seldecomp wage i.sep i.south, by(married)" 2000 ///
    "stage 1: every focal group observation has a level of sep absent from the reference group" ///
    "!duplicate"
s2_msg "seldecomp wage i.south i.sep, by(married) sequentially fixothermargins" 2000 "stage 2: every focal group observation has a level of sep"
s2_msg "seldecomp wage i.sep i.south, by(married) insupport" 2000 ///
    "insupport excludes every focal group observation (levels of sep absent from" ///
    "!duplicate"
gen byte x6 = cond(married == 0 & mod(_n, 2) == 0, 9, south)
gen byte x7 = cond(married == 0 & mod(_n, 2) == 1, 9, collgrad)
count if married == 0 & x6 == 9 & !missing(wage)
local n6 = strtrim(string(r(N), "%12.0fc"))
count if married == 0 & x7 == 9 & !missing(wage)
local n7 = strtrim(string(r(N), "%12.0fc"))
s2_msg "seldecomp wage x6 x7, by(married)" 2000 "every focal group observation has a level absent from the reference group (x6: `n6' obs; x7: `n7' obs)" "!duplicate"
s2_msg "seldecomp wage x6 x7, by(married) insupport" 2000 "insupport excludes every focal group observation (levels absent from the reference group: x6: `n6' obs; x7: `n7' obs)"
* 430 names the worst margin, with its raked share and its target, and
* says that the raking stopped at iterate() only when it did
s2_msg "seldecomp wage i.collgrad i.south i.race, by(married) iterate(1)" 430 ///
    "2 margins missed by more than ctrltolerance() (raking stopped at iterate()); worst: south = 1 (South)"
restore
preserve
clear
set obs 400
gen byte g = _n > 200
gen byte a = mod(_n, 2)
gen byte b = a
replace b = 1 - a if g == 1 & _n > 350
replace a = 1 if g == 1 & _n <= 300
gen double yy = mod(_n, 7) + a + b
* the shares printed are those of the raked weights kept with force
quietly seldecomp yy i.a i.b, by(g) force generate(s2rw)
sum s2rw if g == 0, meanonly
local tot = r(sum)
sum s2rw if g == 0 & a == 0, meanonly
local rsh = strtrim(string(100 * r(sum) / `tot', "%9.2f"))
count if g == 1
local n2 = r(N)
count if g == 1 & a == 0
local tsh = strtrim(string(100 * r(N) / `n2', "%9.2f"))
assert "`rsh'" == "50.00" & "`tsh'" == "25.00"
s2_msg "seldecomp yy i.a i.b, by(g)" 430 ///
    "raking did not reach its targets: 1 margin missed by more than ctrltolerance(); worst: a = 0 (raked share `rsh'%, target `tsh'%); specify force to keep the results" ///
    "!collapse levels" "!iterate()"
s2_msg "seldecomp yy (A: i.a) (B: i.b), by(g) sequentially fixothermargins" 430 "raking did not reach its targets at stage 2: 1 margin missed" "worst: a = 0"
restore

* with cut points from one group's sample: raking the other group cuts
* at the other group's quantiles, and the advice follows that run;
* when that run fails nothing is said about it
preserve
clear
set obs 700
gen byte g = _n > 400
gen double x = 4 * (_n - 0.5) / 400 if g == 0
replace x = (_n - 400.5) / 100 if inrange(_n, 401, 500)
replace x = 3 + (_n - 500.5) / 100 if inrange(_n, 501, 600)
replace x = 10 + (_n - 600.5) / 100 if inrange(_n, 601, 700)
gen double yy = x + mod(_n, 3)
s2_msg "seldecomp yy c.x, by(g) contvars(x 4, refgroup)" 459 "x = > 10.25" ///
    "collapse or drop these levels" "!or rake the other group"
s2_msg "seldecomp yy c.x, by(g) contvars(x 4, refgroup) focalgroup(1)" 459 "x = (1, 2]" ///
    "collapse or drop these levels" "!or rake the other group"
* the other group too small for the quantile groups: that run fails (2001)
keep if g == 1 | _n <= 3
s2_msg "seldecomp yy c.x, by(g) contvars(x 4, refgroup)" 459 "collapse or drop these levels" "!or rake the other group"
s2_msg "seldecomp yy c.x, by(g) contvars(x 4, refgroup) focalgroup(1)" 2001
restore
preserve
clear
set obs 800
gen byte g = _n > 400
gen double x = 4 * (_n - 0.5) / 400 if g == 0
replace x = 4 * (_n - 400.5) / 100 if inrange(_n, 401, 500)
replace x = 6 + 4 * (_n - 500.5) / 300 if _n > 500
gen double yy = x + mod(_n, 3)
* suggested (the other group's quartiles have both groups), and it runs
s2_msg "seldecomp yy c.x, by(g) contvars(x 4, refgroup)" 459 "or rake the other group with focalgroup(1)"
s2_msg "seldecomp yy c.x, by(g) contvars(x 4, refgroup) focalgroup(1)" 0
* cut points from the Focal Group's sample: the other run cuts at the Reference Group's
* quantiles, and that run works
clear
set obs 90
set seed 43
gen byte g = runiform() < 0.5
gen double x = cond(g, 2 * runiform()^2, runiform() * 1.5 + (runiform() < 0.2) * 2)
gen byte x1 = 1 + floor(runiform()^2 * 5)
gen double yy = x + x1
local cmd "seldecomp yy (A: i.x1) (B: c.x), by(g) contvars(x 5, focalgroup) mincell(1) sequentially fixothermargins"
s2_msg "`cmd'" 459 "or rake the other group with focalgroup(1)" "!would fail"
s2_msg "`cmd' focalgroup(1)" 0
* a replicate of a prefix uses the cut points of the full-sample call, and those of the other
* run are not known there: nothing is said about the other group (the jackknife replicate
* without the one Focal Group observation of the top bin fails); a replicate call prints nothing;
* the prefix's display gives one warning; the level is a thin target level, and the warning
* names it, while the thin-cell note is listed under the raking diagnostics, after their
* footnote
clear
set obs 81
gen byte g = _n > 41
gen double x = 4 * (_n - 0.5) / 40 if _n <= 40
replace x = 10 in 41
replace x = 4 * (_n - 41.5) / 30 if inrange(_n, 42, 71)
replace x = 9 + 2 * (_n - 71.5) / 10 if _n > 71
gen double yy = x + mod(_n, 3)
s2_msg "jackknife, nodots: seldecomp yy c.x, by(g) contvars(x 4, refgroup) diagnostics" 0 ///
    "1 of 81 jackknife replicates failed and was dropped: in it, x = > 6.516667 had no focal group weight. Collapse that level." ///
    "Thin target levels (<10 focal obs): 1." ///
    "!in the note above" "!note: thin target level" ///
    "!replicate failed" "!raking support failure" ///
    "!rake the other group" "!would fail as well"
restore

* the 430 of the reverse order names the worst margin and its shares,
* as the forward order does: the run with the blocks typed in reverse order gives the same
* message
preserve
clear
set obs 2000
gen byte g = _n > 1000
gen byte a = cond(_n <= 500 | inrange(_n, 1001, 1200), 0, 1)
gen byte c = cond(_n <= 800 | inrange(_n, 1001, 1100) | inrange(_n, 1201, 1400), 0, 1)
gen double yy = mod(_n, 7) + a + 2 * c
s2_msg "seldecomp yy (A: i.a) (C: i.c), by(g) sequentially fixothermargins" 0
s2_msg "seldecomp yy (C: i.c) (A: i.a), by(g) sequentially fixothermargins" 430 ///
    "raking did not reach its targets at stage 1: 1 margin missed by more than ctrltolerance(); worst: c = 0 (raked share 50.00%, target 30.00%)"
local fwd `"`r(all)'"'
s2_msg "seldecomp yy (A: i.a) (C: i.c), by(g) sequentially fixothermargins reverse" 430 "!margin(s)"
assert `"`r(all)'"' == `"reverse order: `fwd'"'
s2_msg "seldecomp yy (A: i.a) (C: i.c), by(g) sequentially fixothermargins reverse force" 0 ///
    "Raking failure counts" "*not*" "!warning"
restore

* an estimation option at replay: the options replay takes; a bad argument keeps its message
preserve
sysuse nlsw88, clear
quietly seldecomp wage i.collgrad i.south, by(married)
local rl "replay takes level(), diagnostics, covariatedetail, citable, varlabels, cformat(), pformat(), sformat(), and makeplot()"
* (the words noheader, notable, nodefinitions, varbyvar, refminusfocal and refgroup() are
* refused at replay too)
foreach o in "by(married)" sequentially fixothermargins individually "bayesboot(reps(20))" "focalgroup(1)" "contvars(ttl_exp 4)" "mincell(5)" foo noheader notable nodefinitions varbyvar refminusfocal "refgroup(1)" {
    mata: st_local("on", ustrregexra(st_local("o"), "[(].*$", "()"))
    s2_msg "seldecomp, `o'" 198 "option `on' not allowed at replay; `rl'"
}
s2_msg "seldecomp, level(abc)" 198 "!not allowed at replay"
s2_msg "seldecomp, cformat(%9.3f)" 0 "!not allowed"
restore

* ==============================================================================
* generate() never overwrites an input of the command; the observations left out for a
* missing value are noted with their counts; the full diagnostics when trimming leaves a
* margin unmet
* ==============================================================================
di as txt _n "---- generate() and the inputs, exclusion notes"

* an input named in generate() is refused (198) with or without replace, before any work;
* the data are untouched
gen double den = 1 + x1
gen byte z_s1 = x2
qui datasignature
local sig0 `"`r(datasignature)'"'
u2_msg "seldecomp y i.x1 [pw=wt], by(group) generate(y, replace)" 198 "generate(): y is the dependent variable; choose another name"
u2_msg "seldecomp y i.x1 [pw=wt], by(group) generate(wt, replace)" 198 "generate(): wt is the weight variable; choose another name"
u2_msg "seldecomp y i.x1 [pw=2*wt+0], by(group) generate(wt, replace)" 198 "generate(): wt is the weight variable"
u2_msg "seldecomp y i.x1 [pw=wt], by(group) generate(wt)" 198 "generate(): wt is the weight variable"
u2_msg "seldecomp y i.x1 [pw=wt], by(group) generate(group, replace)" 198 "generate(): group is the by() variable"
u2_msg "seldecomp y i.x1 [pw=wt], by(group) generate(x1, replace)" 198 "generate(): x1 is a raking variable"
u2_msg "seldecomp y i.x1#i.x2 [pw=wt], by(group) generate(x2, replace)" 198 "generate(): x2 is a raking variable"
u2_msg "seldecomp y i.x1 c.x3 [pw=wt], by(group) contvars(x3 4) generate(x3, replace)" 198 "generate(): x3 is a contvars() variable"
u2_msg "seldecomp y i.x1 [pw=wt], by(group) outcomes(ybin) generate(ybin, replace)" 198 "generate(): ybin is an outcomes() variable"
u2_msg "seldecomp y i.x1 [pw=wt], by(group) ratio(den) generate(den, replace)" 198 "generate(): den is the ratio() variable"
u2_msg "seldecomp y i.x1 [pw=wt], by(group) bayesboot(reps(2) cluster(cluster)) generate(cluster, replace)" 198 "generate(): cluster is the bayesboot(cluster()) variable"
* the stage variables of sequentially (newvar_s#) too
u2_msg "seldecomp y (a: x1) (b: z_s1) [pw=wt], by(group) sequentially fixothermargins generate(z, replace)" 198 "generate(): z_s1 is a raking variable"
qui datasignature
assert `"`r(datasignature)'"' == `"`sig0'"'
* other names are replaced with replace
gen double rwold = 0
seldecomp y i.x1 [pw=wt], by(group) generate(rwold, replace)
assert rwold > 0 if e(sample)
drop rwold z_s1
* with varabbrev on, replace overwrites only a variable of exactly that name: generate(keep,
* replace) leaves keepme, and generate(wtc, replace) the weight variable wtcopy, although
* drop takes an abbreviation
gen double keepme = 1
gen double wtcopy = wt
local vab = c(varabbrev)
set varabbrev on
seldecomp y i.x1 [pw=wtcopy], by(group) generate(keep, replace)
seldecomp y i.x1 [pw=wtcopy], by(group) generate(wtc, replace)
set varabbrev `vab'
confirm variable keepme wtcopy keep wtc, exact
assert keepme == 1 & wtcopy == wt
assert keep > 0 & wtc > 0 if e(sample)
drop keepme wtcopy keep wtc

capture program drop s1_count
program define s1_count, rclass
    // s1_count "command" rc "text": the command fails with rc (0: runs); r(n) is the number of
    // times its output contains text once the lines are joined and runs of blanks squeezed
    // (the notes are {p} paragraphs, wrapped at 78 columns whatever the line size)
    args cmd rc t
    tempfile lg
    log using `"`lg'"', text replace name(s1log)
    capture noisily `cmd'
    local r = _rc
    log close s1log
    if `r' != `rc' {
        di as err `"`cmd': rc `r', expected `rc'"'
        exit 9
    }
    mata: _t1s = stritrim(invtokens(cat(st_local("lg"))')); _t1t = st_local("t")
    mata: st_local("n", strofreal((strlen(_t1s) - strlen(subinstr(_t1s, _t1t, ""))) / strlen(_t1t)))
    capture mata: mata drop _t1s _t1t
    return scalar n = `n'
end
* notes for a missing weight, a
* missing value of an outcomes() or ratio() variable, a missing bayesboot() cluster id; with
* the counts by group, as the zero-weight note; in that order when several apply
gen double wm = wt
replace wm = . in 1/5
count if group == 0 in 1/5
local a1 = r(N)
s1_count "seldecomp y i.x1 [pw=wm], by(group)" 0 "note: 5 observations with a missing weight excluded (`a1' in the focal group, `=5 - `a1'' in the reference group)"
assert r(n) == 1
assert e(N) == 1995 & e(N_focal) + e(N_ref) == 1995 & e(N_zerowt) == 0
* the counts of the notes have thousands separators, as the tables' numbers do
* (the data doubled, so that both notes count more than 1,000)
preserve
expand 2
gen double wk = wt
replace wk = . in 1/1200
replace wk = 0 in 1201/2300
count if group == 0 in 1/1200
local k1 = r(N)
count if group == 0 in 1201/2300
local k2 = r(N)
seldecomp y i.x1 [pw=wk], by(group)
assert `"`e(note1)'"' == "1,200 observations with a missing weight excluded (`=strtrim(string(`k1', "%12.0fc"))' in the focal group, `=strtrim(string(1200 - `k1', "%12.0fc"))' in the reference group)"
assert `"`e(note2)'"' == "1,100 observations with zero weight excluded (`=strtrim(string(`k2', "%12.0fc"))' in the focal group, `=strtrim(string(1100 - `k2', "%12.0fc"))' in the reference group)"
assert e(N_zerowt) == 1100
restore
gen str1 sb2 = cond(group, "b", "a")
s1_count "seldecomp y i.x1 [pw=wm], by(sb2)" 0 "note: 5 observations with a missing weight excluded (`a1' in the focal group, `=5 - `a1'' in the reference group)"
assert r(n) == 1
drop sb2
gen double ym = y
replace ym = . in 6/16
count if group == 0 in 6/16
local a2 = r(N)
s1_count "seldecomp y i.x1 [pw=wm], by(group) outcomes(ym)" 0 "note: 11 observations with a missing value of an outcomes() variable excluded (`a2' in the focal group, `=11 - `a2'' in the reference group); every outcome uses the same observations"
assert r(n) == 1
assert e(N) == 1984 & e(N_notes) >= 2
assert strpos(`"`e(note1)'"', "outcomes()") & strpos(`"`e(note2)'"', "missing weight")
replace den = . in 20/22
count if group == 0 in 20/22
local a3 = r(N)
s1_count "seldecomp y i.x1, by(group) ratio(den)" 0 "note: 3 observations with a missing value of den (the ratio() variable) excluded (`a3' in the focal group, `=3 - `a3'' in the reference group)"
assert r(n) == 1
s1_count "seldecomp y i.x1, by(group) ratio(den) outcomes(ym)" 0 "note: 14 observations with a missing value of an outcomes() or the ratio() variable excluded"
assert r(n) == 1
gen int cm = cluster
replace cm = . in 30/40
count if group == 0 in 30/40
local a4 = r(N)
s1_count "seldecomp y i.x1, by(group) bayesboot(reps(2) cluster(cm) rseed(1) nodots)" 0 "note: 11 observations with a missing bayesboot() cluster id excluded (`a4' in the focal group, `=11 - `a4'' in the reference group)"
assert r(n) == 1
* no note in a replicate call of a prefix (c(prefix) set as the prefix sets it)
capture program drop s1_fake
program define s1_fake, eclass
    gettoken words 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    set prefix `words'
    `0'
end
s1_fake bootstrap _loop _loop_bs : seldecomp y i.x1 [pw=wm], by(group) outcomes(ym)
forvalues m = 1/`e(N_notes)' {
    assert !strpos(`"`e(note`m')'"', "missing")
}
assert e(N) == 1984
drop wm ym den cm

* a margin that trimming left unmet: rc 430, the message naming the bounds; with force
* e(converged) 0, e(targets_met) 0 and one note, the trimming note; a miss kept by force after
* the raking stopped at the iterate() limit: one note, which says both
s1_count "seldecomp y x1 x2 [pw=wt], by(group) trimhirel(1.05)" 430 "the trimming bounds are binding: relax them, or specify force to keep the results"
assert r(n) == 1
s1_count "seldecomp y x1 x2 [pw=wt], by(group) trimhirel(1.05) force" 0 "note: trimming left"
assert r(n) == 0
assert e(converged) == 0 & e(targets_met) == 0 & e(badcontrols) > 0 & e(N_notes) == 1
assert strpos(`"`e(note1)'"', "trimming left ") == 1
assert strpos(`"`e(note1)'"', "); results kept because force was specified") > 0
s1_count "seldecomp" 0 "Converged with Bad"
assert r(n) == 0
u5_count "seldecomp y x1 x2 [pw=wt], by(group) iterate(1) force" 0 "note:"
assert r(n) == 0
assert e(converged) == 0 & e(targets_met) == 0 & e(N_notes) == 1
assert strpos(`"`e(note1)'"', "raking stopped at the iterate() limit with `=e(badcontrols)' raking margin") == 1
assert substr(`"`e(note1)'"', -strlen("), kept by force"), .) == "), kept by force"
assert strpos(`"`e(note1)'"', "raise iterate()") == 0
seldecomp y x1 x2 [pw=wt], by(group)
assert e(converged) == 1 & e(targets_met) == 1

* e(converged) is the minimum over the stages, and with reverse over the
* reverse order's too: 1 only when every raking converged; column
* converged of e(bystage) and e(bystage_rev) keeps each stage's value. Margins that trimming
* left unmet at both stages of sequentially fixothermargins (rc 430 without force), kept by
* force: 0 at both
rcof "noi seldecomp y (A: x1) (B: x2) [pw=wt], by(group) sequentially fixothermargins trimhirel(1.05)" == 430
seldecomp y (A: x1) (B: x2) [pw=wt], by(group) sequentially fixothermargins trimhirel(1.05) force
assert el(e(bystage), 1, 6) == 0 & el(e(bystage), 2, 6) == 0
assert el(e(bystage), 1, 7) > 0 & el(e(bystage), 2, 7) > 0
assert e(converged) == 0 & e(targets_met) == 0
* Focal Group has a == 1 <=> b == 0 (shares .5), Reference Group a .3 and b .7: stage 1 of sequentially
* fixothermargins (a raked to .3, b held at the Focal Group's .5) cannot be met, the last stage can
preserve
clear
set obs 400
gen byte grp = _n > 200
gen byte a = cond(grp, mod(_n, 10) < 3, mod(_n, 2))
gen byte b = 1 - a
gen double yy = a + b + mod(_n, 7) / 7
rcof "noi seldecomp yy (A: a) (B: b), by(grp) sequentially fixothermargins" == 430
seldecomp yy (A: a) (B: b), by(grp) sequentially fixothermargins force generate(rc7)
assert el(e(bystage), 1, 6) == 0 & el(e(bystage), 2, 6) == 1
assert e(converged) == 0 & e(targets_met) == 0 & e(badcontrols) == 0
assert "`: char rc7_s1[seldecomp_converged]'" == "0" & "`: char rc7[seldecomp_converged]'" == "1"
drop rc7*
* reverse: Focal Group lacks the cell a == 1 & b == 1 (shares a .4, b .3), Reference Group has a .3 and b
* .65. Both forward stages can be met (.3 + .3 and .3 + .65 are at most 1); stage 1 of the
* reverse order (b raked to .65, a held at .4) cannot
clear
input byte(grp a b) int n
0 1 0 80
0 0 1 60
0 0 0 60
1 1 0 20
1 1 1 40
1 0 1 90
1 0 0 50
end
expand n
gen double yy = a + 2 * b + mod(_n, 7) / 7
seldecomp yy (A: a) (B: b), by(grp) sequentially fixothermargins
assert e(converged) == 1 & e(targets_met) == 1
rcof "noi seldecomp yy (A: a) (B: b), by(grp) sequentially fixothermargins reverse" == 430
seldecomp yy (A: a) (B: b), by(grp) sequentially fixothermargins reverse force
assert el(e(bystage), 1, 6) == 1 & el(e(bystage), 2, 6) == 1
assert el(e(bystage_rev), 1, 6) == 0 & el(e(bystage_rev), 2, 6) == 1
assert e(converged) == 0
restore


* ==============================================================================
* the lower trimming bounds bind on the Focal Group observations with a level
* absent from the Reference Group: they stay in the raking, their level's target 0, and the
* bound holds them. The margin so missed has its own messages: a warning (not with
* force), then error 430; force keeps the results as not converged. trimsupportoverride
* gives them zero raked weight whatever the bounds, and requires a lower bound
* ==============================================================================
di as txt _n "---- lower bounds and out-of-support observations"
* the error messages unwrapped in the logs (the warning wraps at 78 columns whatever the linesize)
local ls0 = c(linesize)
set linesize 255
preserve
clear
set seed 3
set obs 600
gen byte g = _n > 300
gen byte x = ceil(3 * runiform())
replace x = 4 in 1/10
gen byte z = ceil(2 * runiform())
gen double y = x + z + g + rnormal()
gen double w = 1 + runiform()
* no lower bound: zero raked weight, no warning
log using `"`lgf'"', text replace name(sdlog)
seldecomp y i.x i.z [pw=w], by(g) generate(t0)
log close sdlog
assert e(N_focal_zero) == 10 & e(N_raked) == 290 & e(targets_met) == 1 & e(converged) == 1
count if t0 == 0 in 1/10
assert r(N) == 10
log_has `"`lgf'"' "warning:"
assert r(found) == 0
* trimloabs(), trimlorel(), both, every trimfrequency(), staged: the warning, then error 430
* naming the cause and the remedies
log using `"`lgf'"', text replace name(sdlog)
capture noisily seldecomp y i.x i.z [pw=w], by(g) trimloabs(0.5)
local rc = _rc
log close sdlog
assert `rc' == 430
foreach t in "warning: 10 focal group observations have a level absent from the reference" ///
    "(x = 4: 10 obs). trimloabs() keeps their raked weights" "trimsupportoverride sets them to 0." ///
    "the lower trimming bounds keep focal group observations" "trimsupportoverride, insupport, or force" {
    log_has `"`lgf'"' "`t'"
    assert r(found) == 1
}
log_has `"`lgf'"' "the trimming bounds did not change the weights"
assert r(found) == 0
log using `"`lgf'"', text replace name(sdlog)
capture noisily seldecomp y i.x i.z [pw=w], by(g) trimloabs(0.5) trimlorel(0.2)
local rc = _rc
log close sdlog
assert `rc' == 430
log_has `"`lgf'"' "trimloabs() and trimlorel() keep their raked"
assert r(found) == 1
rcof "noi seldecomp y i.x i.z [pw=w], by(g) trimlorel(0.3)" == 430
foreach f in often once {
    rcof "noi seldecomp y i.x i.z [pw=w], by(g) trimloabs(0.5) trimfrequency(`f')" == 430
}
rcof "noi seldecomp y (a: i.x) (b: i.z) [pw=w], by(g) sequentially trimloabs(0.5)" == 430
rcof "noi seldecomp y (a: i.x) (b: i.z) [pw=w], by(g) individually trimloabs(0.5)" == 430
* force: no warning; the results kept, the observations held at the bound (in the raking
* sample, not counted as zero weight), not converged, targets not met, the force note
log using `"`lgf'"', text replace name(sdlog)
seldecomp y i.x i.z [pw=w], by(g) generate(t1) trimloabs(0.5) force
log close sdlog
assert e(converged) == 0 & e(targets_met) == 0 & e(N_focal_zero) == 0 & e(N_raked) == 300
assert el(e(bystage), 1, colnumb(e(bystage), "converged")) == 0
assert missing(e(explained_support))
count if reldif(t1, 0.5) < 1e-8 in 1/10
assert r(N) == 10
count if t1 < 0.5 - 1e-8 & g == 0
assert r(N) == 0
log_has `"`lgf'"' "warning:"
assert r(found) == 0
log_has `"`lgf'"' "results kept because force"
assert r(found) == 0
log_has `"`lgf'"' "Raking failure counts"
assert r(found) == 1
log_has `"`lgf'"' "absent from the reference group"
assert r(found) == 0
seldecomp y i.x i.z [pw=w], by(g) generate(t2) trimlorel(0.3) force
assert e(converged) == 0
count if reldif(t2, 0.3 * w) < 1e-8 in 1/10
assert r(N) == 10
seldecomp y (a: i.x) (b: i.z) [pw=w], by(g) sequentially fixothermargins reverse detail trimloabs(0.5) force
assert e(converged) == 0
* bayesboot(): every replicate misses the target of x = 4 too, and force keeps them all
* (e(N_forced), and a note)
log using `"`lgf'"', text replace name(sdlog)
seldecomp y i.x i.z [pw=w], by(g) trimloabs(0.5) force bayesboot(reps(5) rseed(1) nodots)
log close sdlog
assert e(N_reps) == 5 & e(N_misreps) == 0 & e(N_forced) == 5
log_has `"`lgf'"' "note: force kept 5 of 5 replications that missed raking targets"
assert r(found) == 1
* trimsupportoverride: zero raked weight whatever the bounds, the other weights bounded, no
* warning, converged; recorded in e(rakeopts)
log using `"`lgf'"', text replace name(sdlog)
seldecomp y i.x i.z [pw=w], by(g) generate(t3) trimloabs(0.5) trimsupportoverride
log close sdlog
assert e(N_focal_zero) == 10 & e(N_raked) == 290 & e(converged) == 1 & e(targets_met) == 1
assert strpos("`e(rakeopts)'", "trimsupportoverride") > 0
count if t3 == 0 in 1/10
assert r(N) == 10
count if t3 < 0.5 - 1e-8 & g == 0 & _n > 10
assert r(N) == 0
log_has `"`lgf'"' "warning:"
assert r(found) == 0
seldecomp y i.x i.z [pw=w], by(g) trimlorel(0.3) trimsupportoverride bayesboot(reps(5) rseed(1) nodots)
assert e(N_focal_zero) == 10 & e(N_reps) == 5
* it requires trimloabs() or trimlorel()
rcof "noi seldecomp y i.x i.z [pw=w], by(g) trimsupportoverride" == 198
rcof "noi seldecomp y i.x i.z [pw=w], by(g) trimhirel(3) trimhiabs(5) trimsupportoverride" == 198
* a lower bound without out-of-support observations: no warning
log using `"`lgf'"', text replace name(sdlog)
seldecomp y i.z [pw=w], by(g) trimloabs(0.5)
log close sdlog
log_has `"`lgf'"' "warning:"
assert r(found) == 0
* chain: x = 4 only in the first value of g3; the pair's warning and error, force, the override
gen byte g3 = mod(_n, 3)
replace x = 1 if x == 4
replace x = 4 if g3 == 0 & _n <= 60
log using `"`lgf'"', text replace name(sdlog)
capture noisily seldecomp y i.x i.z [pw=w], by(g3) chain trimloabs(0.5)
local rc = _rc
log close sdlog
assert `rc' == 430
log_has `"`lgf'"' "warning: pair 1: 20 focal group observations"
assert r(found) == 1
seldecomp y i.x i.z [pw=w], by(g3) chain trimloabs(0.5) force
assert e(converged) == 0
seldecomp y i.x i.z [pw=w], by(g3) chain trimloabs(0.5) trimsupportoverride
assert e(converged) == 1
restore
set linesize `ls0'

di as res _n "PASSED: test_options.do"

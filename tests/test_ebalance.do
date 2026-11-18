* test_ebalance.do -- cross-check of the raking against entropy balancing: ebalance
* (Hainmueller and Xu, SSC), someone else's code for the same optimisation problem.
*
* The equivalence.  Raking Focal Group to the Reference Group's margins (iterative proportional fitting)
* converges to the weights w_i = d_i * prod_k a_k(x_ik) that minimise the Kullback-Leibler
* divergence sum_i w_i log(w_i / d_i) from the base weights d_i subject to the margins
* (Ireland and Kullback 1968; the raking case of Deville and Sarndal 1992).  Entropy
* balancing (Hainmueller 2012) minimises the same divergence from the base weights subject
* to the weighted means of the covariates and the weight total, and its solution is
* w_i = d_i exp(c_i'lambda) / normalisation.  When the covariates c_i are the indicators
* of the levels of each raking variable (one level of each left out, as it is implied by
* the total), exp(c_i'lambda) is a product of one factor per level, the means are the
* level shares and the constraints are the margins: the two problems are the same, the
* solution is unique, and the weights agree up to a constant factor.  So
*     raked Focal Group mean (seldecomp)  =  Focal Group mean under the ebalance weights
* and the raked weights are proportional to the ebalance weights.
*
* When it holds: when the margins can be met with positive weights (an interior solution).
* It is tested here with one categorical raking variable (unweighted and weighted), with
* two and with three (weighted), with a full interaction (the joint cells), and for stage 1
* of sequentially fixothermargins (the added block at the Reference Group's shares and the others held
* at the Focal Group's own, entered as ebalance's manualtargets(); unweighted, because ebalance
* 1.5.5 sets every base weight to 1 in its manualtargets() mode) and of sequentially alone
* (the added block at the Reference Group's shares, the others free: the targets of the added block
* only).  Not covered: a Focal Group level absent from
* Reference Group (raking gives those observations weight 0, a boundary point that entropy
* balancing, whose weights are all positive, cannot reach), and trimming (a different
* problem).
*
* Tolerances: with tolerance(1e-12) and ebalance tolerance(1e-10) the means
* agree to 1e-9 (reldif) and the normalised weights to 1e-7; with seldecomp's default
* tolerance(1e-8) the means agree to 1e-6.
*
* ebalance is optional: run_tests.do reports this file as SKIPPED when it is not
* installed (ssc install ebalance).
*
* Run from the tests directory:   . do test_ebalance.do

version 16
clear
set more off
set varabbrev off
discard             // drop any seldecomp already in memory; the repository copy is what gets tested

* ---- tolerances (also written out in eb_check) ----------------------------------
local mtol  1e-9    // raked mean vs ebalance mean, tight tolerances on both sides
local wtol  1e-7    // normalised raked weights vs normalised ebalance weights
local dtol  1e-6    // raked mean at seldecomp's default tolerance vs ebalance mean

* ---- locate seldecomp.ado and ebalance ----------------------------------------
if `"$SELDECOMP_ROOT"' != "" adopath ++ `"$SELDECOMP_ROOT"'
else adopath ++ `"`c(pwd)'/.."'
capture which seldecomp
if _rc {
    di as err "seldecomp.ado not found: run this do-file from the tests directory " ///
        "or set the global SELDECOMP_ROOT to the repository root"
    exit 601
}
which seldecomp        // record in the log which copy is being tested
capture which ebalance
if _rc {
    di as err "ebalance is required by this test; install it with:  ssc install ebalance"
    exit 111
}
which ebalance

* ---- helpers ------------------------------------------------------------------
capture program drop eb_check
program define eb_check
    // eb_check "covariates" "seldecomp options" [weight] -- runs ebalance (Reference Group =
    // treated: its weighted means are the targets; Focal Group = controls, reweighted) and
    // compares with seldecomp y <raking terms>: the raked mean (tight and default
    // tolerances) and the weights
    args covs terms wexp
    local wopt
    local ebw
    if "`wexp'" != "" {
        local wopt [pw=`wexp']
        local ebw basewt(`wexp') wttreat
    }
    else {
        capture drop one_
        gen double one_ = 1
        local wexp one_
    }
    quietly seldecomp y `terms' `wopt', by(group) tolerance(1e-12) generate(rk_)
    local m_tight = e(mean_raked_focal)
    quietly seldecomp y `terms' `wopt', by(group)
    local m_def = e(mean_raked_focal)
    capture drop eb_
    quietly ebalance D `covs', `ebw' tolerance(1e-10) maxiter(200) generate(eb_)
    assert e(convg) == 1
    quietly sum y [aw=eb_] if D == 0
    local m_eb = r(mean)
    di as txt %-24s "`terms'" "  raked " %18.0g `m_tight' "  ebalance " %18.0g `m_eb' ///
        "  reldif " %9.2e reldif(`m_tight', `m_eb')
    assert reldif(`m_tight', `m_eb') < 1e-9
    assert reldif(`m_def', `m_eb') < 1e-6
    * the weights: rk_ / sum(rk_) = eb_ / sum(eb_) in the Focal Group
    quietly sum rk_ if D == 0
    local srk = r(sum)
    quietly sum eb_ if D == 0
    local seb = r(sum)
    tempvar rd
    quietly gen double `rd' = reldif(rk_ / `srk', eb_ / `seb') if D == 0
    quietly sum `rd'
    di as txt _col(26) "max reldif of the normalised weights " %9.2e r(max)
    assert r(max) < 1e-7 & r(N) == 888
    drop rk_ eb_
    capture drop one_
end

* ---- data ------------------------------------------------------------------------
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
recast double y x3 wt
assert _N == 2000
gen byte D = group == 1             // ebalance's treated group: seldecomp's Reference Group
tabulate x1, generate(d1_)          // level indicators
tabulate x2, generate(d2_)
tabulate x4, generate(d4_)
egen c12 = group(x1 x2)
tabulate c12, generate(c_)
* doubles: in its manualtargets() mode ebalance writes the targets into the covariates
* themselves (in an added observation), and replace would store them in a byte indicator
* promoted to float, rounding the targets to 7 digits
recast double d1_* d2_* d4_* c_*


* ==============================================================================
* 1. overall results: one, two and three raking variables, and a full interaction
* ==============================================================================
di as txt _n "---- 1. raked Focal Group mean = mean under the ebalance weights"
eb_check "d1_2 d1_3" "x1"
eb_check "d1_2 d1_3" "x1" wt
eb_check "d1_2 d1_3 d2_2" "x1 x2" wt
eb_check "d1_2 d1_3 d2_2 d4_2 d4_3 d4_4" "x1 x2 x4" wt
eb_check "c_2 c_3 c_4 c_5 c_6" "x1##x2" wt
* all levels given to ebalance (one per variable is collinear with the total and dropped)
eb_check "d1_1 d1_2 d1_3 d2_1 d2_2" "x1 x2" wt


* ==============================================================================
* 2. stage 1 of sequentially fixothermargins: block 1 at the Reference Group's shares, block 2 held at
*    Focal Group's own
* ==============================================================================
di as txt _n "---- 2. stage 1 of sequentially fixothermargins (manualtargets)"
seldecomp y (a: x1) (b: x2 x4), by(group) sequentially fixothermargins tolerance(1e-12)
matrix G = e(bystage)
local m_s1 = G[1, colnumb(G, "raked_focal")]
* targets: the x1 shares of the Reference Group, the x2 and x4 shares of the Focal Group (unweighted)
local tg
foreach v in d1_2 d1_3 {
    quietly sum `v' if D == 1
    local tg `tg' `: di %21.17g r(mean)'
}
foreach v in d2_2 d4_2 d4_3 d4_4 {
    quietly sum `v' if D == 0
    local tg `tg' `: di %21.17g r(mean)'
}
quietly ebalance d1_2 d1_3 d2_2 d4_2 d4_3 d4_4 if D == 0, manualtargets(`tg') ///
    tolerance(1e-10) maxiter(200) generate(eb_s1)
assert e(convg) == 1
quietly sum y [aw=eb_s1] if D == 0
di as txt "stage 1: raked " %18.0g `m_s1' "  ebalance " %18.0g r(mean) ///
    "  reldif " %9.2e reldif(`m_s1', r(mean))
assert reldif(`m_s1', r(mean)) < `mtol'
* and the last stage is the overall raking
eb_check "d1_2 d1_3 d2_2 d4_2 d4_3 d4_4" "x1 x2 x4"


* ==============================================================================
* 3. stage 1 of sequentially alone: block 1 at the Reference Group's shares, block 2 free (the
*    targets of block 1 only)
* ==============================================================================
di as txt _n "---- 3. stage 1 of sequentially (manualtargets, block 1 only)"
seldecomp y (a: x1) (b: x2 x4), by(group) sequentially tolerance(1e-12)
matrix G = e(bystage)
local m_f1 = G[1, colnumb(G, "raked_focal")]
local tg
foreach v in d1_2 d1_3 {
    quietly sum `v' if D == 1
    local tg `tg' `: di %21.17g r(mean)'
}
quietly ebalance d1_2 d1_3 if D == 0, manualtargets(`tg') ///
    tolerance(1e-10) maxiter(200) generate(eb_f1)
assert e(convg) == 1
quietly sum y [aw=eb_f1] if D == 0
di as txt "stage 1: raked " %18.0g `m_f1' "  ebalance " %18.0g r(mean) ///
    "  reldif " %9.2e reldif(`m_f1', r(mean))
assert reldif(`m_f1', r(mean)) < `mtol'
* the two modes differ at stage 1: x2 and x4 are free here, held in section 2
assert reldif(`m_f1', `m_s1') > 1e-4

di as res _n "PASSED: test_ebalance.do"

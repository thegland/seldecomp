* test_oaxaca.do -- equivalence of seldecomp with oaxaca (Jann) where raking is
* post-stratification.
*
* seldecomp's Focal Group (raked) must be oaxaca's group 1 (the lower
* by() value, or the higher one with -swap-), and oaxaca must use group 1's
* coefficients as the reference coefficients (-weight(1)-).  Then
*     seldecomp difference  = oaxaca overall:difference
*     seldecomp explained   = oaxaca overall:explained   = (X1 - X2) b1
*     seldecomp unexplained = oaxaca overall:unexplained = X2 (b1 - b2)
* exactly when the raking targets saturate the cells (one categorical raking
* variable, or a full interaction matched by a fully interacted dummy set) and
* only approximately when two or more categorical variables are raked on their
* marginals.  oaxaca does not accept factor-variable
* notation, so dummies are created with -tabulate, generate()-.
*
* seldecomp's difference is always Focal Group - Reference Group, the Focal Group being the
* group raked. oaxaca's weight(0) (group 2's coefficients as the reference) therefore
* corresponds to seldecomp with focalgroup() naming oaxaca's group 2, and oaxaca's swap weight(0) to
* seldecomp's default: the groups trade places and difference, explained and unexplained
* are the negatives of oaxaca's (section vi).
*
* Run from the tests directory:   . do test_oaxaca.do

version 16
clear
set more off
set varabbrev off
discard             // drop any seldecomp already in memory; the repository copy is what gets tested

* ---- tolerances ---------------------------------------------------------------
local tol  1e-6     // exact cases: reldif() of every overall entry
local atol 0.15     // additive case: |explained_raking - explained_OLS| must be below this
local itol 1e-10

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

capture program drop assert_same
program define assert_same
    // assert_same matA matB colname [tol]  --  element (1, colname) equal in A (oaxaca's e(b))
    // and B (seldecomp's), where oaxaca's group_1 and group_2 are focal_group and ref_group
    args A B name tol
    if "`tol'" == "" local tol 1e-6
    local nameb : subinstr local name "group_1" "focal_group"
    local nameb : subinstr local nameb "group_2" "ref_group"
    local anames : colfullnames `A'
    local bnames : colfullnames `B'
    local ca : list posof "`name'" in anames
    local cb : list posof "`nameb'" in bnames
    if `ca' < 1 | `cb' < 1 {
        di as err "assert_same: column `name' not found in `A' (col `ca') or `B' (col `cb')"
        exit 9
    }
    tempname va vb
    scalar `va' = el(`A', 1, `ca')
    scalar `vb' = el(`B', 1, `cb')
    if missing(`va') | missing(`vb') {
        di as err "`name': missing value (`A' = " `va' ", `B' = " `vb' ")"
        exit 9
    }
    if reldif(`va', `vb') > `tol' {
        di as err "`name': `A' = " %18.0g `va' "  `B' = " %18.0g `vb' ///
            "  (reldif = " %9.3g reldif(`va', `vb') " > `tol')"
        exit 9
    }
    di as txt "  ok  " %-22s "`name'" %14.0g `va' "  vs " %14.0g `vb'
end

capture program drop assert_neg
program define assert_neg
    // assert_neg oaxaca_b seldecomp_b tol -- seldecomp's focal_group/ref_group are oaxaca's
    // group_2/group_1, and its difference, explained and unexplained oaxaca's negatives
    args O S tol
    local onames : colfullnames `O'
    local snames : colfullnames `S'
    foreach p in "focal_group group_2" "ref_group group_1" "difference difference" ///
        "explained explained" "unexplained unexplained" {
        gettoken sn on : p
        local on = strtrim("`on'")
        local cs : list posof "overall:`sn'" in snames
        local co : list posof "overall:`on'" in onames
        if `cs' < 1 | `co' < 1 {
            di as err "assert_neg: overall:`sn' (seldecomp) or overall:`on' (oaxaca) not found"
            exit 9
        }
        tempname sv ov
        scalar `sv' = el(`S', 1, `cs')
        scalar `ov' = cond("`sn'" == "`on'", -1, 1) * el(`O', 1, `co')
        if missing(`sv') | missing(`ov') | reldif(`sv', `ov') > `tol' {
            di as err "seldecomp overall:`sn' " %18.0g `sv' " vs " ///
                cond("`sn'" == "`on'", "minus ", "") "oaxaca overall:`on' " %18.0g `ov'
            exit 9
        }
        di as txt "  ok  " %-22s "overall:`sn'" %14.0g `sv' "  vs " %14.0g `ov'
    }
end

* ---- locate seldecomp.ado and oaxaca ------------------------------------------
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
capture which oaxaca
if _rc {
    di as err "oaxaca is required by this test; install it with:  ssc install oaxaca"
    exit 111
}

* ---- data and expected values -------------------------------------------------
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
recast double y x3 wt
confirm string variable x4
assert _N == 2000
do expected_values.do

* dummy sets for oaxaca
tabulate x1, generate(d1_)                  // d1_1 d1_2 d1_3
assert r(r) == 3
egen long cell = group(x1 x2)               // joint cells of x1 x2
tabulate cell, generate(dc_)                // dc_1 ... dc_6
assert r(r) == 6
tabulate x4, generate(d4_)                  // d4_1 ... d4_4 (east north south west)
assert r(r) == 4

local overall group_1 group_2 difference explained unexplained


* ==============================================================================
* (i) single categorical variable, unweighted
* ==============================================================================
di as txt _n "---- (i)  oaxaca y d1_2 d1_3, by(group) weight(1)   vs   seldecomp y i.x1, by(group)"
oaxaca y d1_2 d1_3, by(group) weight(1)
assert real("`e(group_1)'") == 0            // oaxaca stores the by() value as text
local oN  = e(N)
local oN1 = e(N_1)
local oN2 = e(N_2)
matrix ob = e(b)
* oaxaca reproduces the reference OLS decomposition
assert_mat ob 1 overall:explained exp_A_ob_explained `tol'

seldecomp y i.x1, by(group)
assert "`e(focal_group)'" == "0"
assert e(N)   == `oN'
assert e(N_focal) == `oN1'
assert e(N_ref) == `oN2'
matrix sb = e(b)
foreach s of local overall {
    assert_same ob sb overall:`s' `tol'
}


* ==============================================================================
* (ii) single categorical variable, pweights
* ==============================================================================
di as txt _n "---- (ii) oaxaca y d1_2 d1_3 [pw=wt], by(group) weight(1)   vs   seldecomp y i.x1 [pw=wt], by(group)"
oaxaca y d1_2 d1_3 [pw=wt], by(group) weight(1)
local oN  = e(N)
local oN1 = e(N_1)
local oN2 = e(N_2)
matrix ob = e(b)
assert_mat ob 1 overall:explained exp_B_ob_explained `tol'

seldecomp y i.x1 [pw=wt], by(group)
assert e(N)   == `oN'
assert e(N_focal) == `oN1'
assert e(N_ref) == `oN2'
matrix sb = e(b)
foreach s of local overall {
    assert_same ob sb overall:`s' `tol'
}

* the same with the string variable x4 (encoded by seldecomp) as the single raking variable
di as txt _n "---- (ii b) oaxaca y d4_2 d4_3 d4_4 [pw=wt], by(group) weight(1)   vs   seldecomp y x4 [pw=wt], by(group)"
oaxaca y d4_2 d4_3 d4_4 [pw=wt], by(group) weight(1)
matrix ob = e(b)
seldecomp y x4 [pw=wt], by(group)
matrix sb = e(b)
foreach s of local overall {
    assert_same ob sb overall:`s' `tol'
}

* and with a binary outcome (linear probability model in oaxaca)
di as txt _n "---- (ii c) oaxaca ybin d1_2 d1_3, by(group) weight(1)   vs   seldecomp ybin i.x1, by(group)"
oaxaca ybin d1_2 d1_3, by(group) weight(1)
matrix ob = e(b)
seldecomp ybin i.x1, by(group)
matrix sb = e(b)
foreach s of local overall {
    assert_same ob sb overall:`s' `tol'
}


* ==============================================================================
* (iii) full interaction x1##x2 vs one dummy per joint cell
* ==============================================================================
di as txt _n "---- (iii) oaxaca y dc_2-dc_6 [pw=wt], by(group) weight(1)   vs   seldecomp y x1##x2 [pw=wt], by(group)"
oaxaca y dc_2-dc_6 [pw=wt], by(group) weight(1)
local oN  = e(N)
local oN1 = e(N_1)
local oN2 = e(N_2)
matrix ob = e(b)
assert_mat ob 1 overall:explained exp_D_ob_explained `tol'

seldecomp y x1##x2 [pw=wt], by(group)
assert e(N)   == `oN'
assert e(N_focal) == `oN1'
assert e(N_ref) == `oN2'
matrix sb = e(b)
foreach s of local overall {
    assert_same ob sb overall:`s' `tol'
}

* unweighted, with the # form
di as txt _n "---- (iii b) unweighted: oaxaca y dc_2-dc_6, by(group) weight(1)   vs   seldecomp y x1#x2, by(group)"
oaxaca y dc_2-dc_6, by(group) weight(1)
matrix ob = e(b)
seldecomp y x1#x2, by(group)
matrix sb = e(b)
foreach s of local overall {
    assert_same ob sb overall:`s' `tol'
}


* ==============================================================================
* (iv) focal group = 1: oaxaca -swap- vs seldecomp focalgroup(1)
* ==============================================================================
di as txt _n "---- (iv) oaxaca y d1_2 d1_3 [pw=wt], by(group) weight(1) swap   vs   seldecomp y i.x1 [pw=wt], by(group) focalgroup(1)"
oaxaca y d1_2 d1_3 [pw=wt], by(group) weight(1) swap
assert real("`e(group_1)'") == 1
local oN  = e(N)
local oN1 = e(N_1)
local oN2 = e(N_2)
matrix ob = e(b)

seldecomp y i.x1 [pw=wt], by(group) focalgroup(1)
assert "`e(focal_group)'" == "1"
assert "`e(ref_group)'" == "0"
assert e(N)   == `oN'
assert e(N_focal) == `oN1'
assert e(N_ref) == `oN2'
matrix sb = e(b)
foreach s of local overall {
    assert_same ob sb overall:`s' `tol'
}

* full interaction with the focal group swapped
di as txt _n "---- (iv b) oaxaca y dc_2-dc_6 [pw=wt], by(group) weight(1) swap   vs   seldecomp y x1##x2 [pw=wt], by(group) focalgroup(1)"
oaxaca y dc_2-dc_6 [pw=wt], by(group) weight(1) swap
matrix ob = e(b)
seldecomp y x1##x2 [pw=wt], by(group) focalgroup(1)
matrix sb = e(b)
foreach s of local overall {
    assert_same ob sb overall:`s' `tol'
}


* ==============================================================================
* (v) additive case: two categorical variables raked on their marginals vs OLS
*     on additive dummies -- structurally different, only close
* ==============================================================================
di as txt _n "---- (v) oaxaca y d1_2 d1_3 x2 [pw=wt], by(group) weight(1)   vs   seldecomp y x1 x2 [pw=wt], by(group)"
oaxaca y d1_2 d1_3 x2 [pw=wt], by(group) weight(1)
local oN  = e(N)
local oN1 = e(N_1)
local oN2 = e(N_2)
matrix ob = e(b)
local obnames : colfullnames ob
local c_expl   : list posof "overall:explained" in obnames
local c_unexpl : list posof "overall:unexplained" in obnames
assert `c_expl' > 0 & `c_unexpl' > 0
local ob_expl   = el(ob, 1, `c_expl')
local ob_unexpl = el(ob, 1, `c_unexpl')
* oaxaca reproduces the reference (weighted OLS) decomposition
assert_mat ob 1 overall:explained exp_C_ob_explained `tol'

seldecomp y x1 x2 [pw=wt], by(group)
assert e(N)   == `oN'
assert e(N_focal) == `oN1'
assert e(N_ref) == `oN2'
matrix sb = e(b)
* observation counts, group means and the raw gap match exactly ...
assert_same ob sb overall:group_1    `tol'
assert_same ob sb overall:group_2    `tol'
assert_same ob sb overall:difference `tol'
* ... and seldecomp reproduces the reference raking decomposition
assert_close e(explained) exp_C_explained `tol' "(v) seldecomp explained vs reference"

local sd_expl   = e(explained)
local sd_unexpl = e(unexplained)
local gap = abs(`sd_expl' - `ob_expl')
di as txt _n "Additive case: x1 and x2 raked on their marginals vs OLS on additive dummies"
di as txt "  seldecomp explained   = " as res %12.6f `sd_expl'   as txt "   unexplained = " as res %12.6f `sd_unexpl'
di as txt "  oaxaca    explained   = " as res %12.6f `ob_expl'   as txt "   unexplained = " as res %12.6f `ob_unexpl'
di as txt "  absolute difference   = " as res %12.6f `gap'       as txt "   (tolerance `atol')"
di as txt "  This discrepancy is expected and is not a convergence problem: raking weights"
di as txt "  are products of one factor per raking variable (log-linear), OLS residuals are"
di as txt "  orthogonal only to the additive dummies, so the two counterfactual means differ"
di as txt "  by a second-order term.  The Python reference found 0.064."
if `gap' >= `atol' {
    di as err "(v) additive case: |seldecomp explained - oaxaca explained| = `gap' >= `atol'"
    exit 9
}
* the unexplained parts differ by the same amount with the opposite sign (same raw gap)
assert reldif(`sd_unexpl' - `ob_unexpl', `ob_expl' - `sd_expl') < 1e-8


* ==============================================================================
* (vi) oaxaca weight(0): group 2's coefficients.  seldecomp rakes oaxaca's group 2 to
*      group 1's distribution with focalgroup() naming it (and oaxaca's swap weight(0)
*      corresponds to seldecomp's default); its difference is its own Focal Group - Reference Group,
*      so the groups trade places and the three parts are the negatives of oaxaca's
* ==============================================================================
di as txt _n "---- (vi) oaxaca y d1_2 d1_3 [pw=wt], by(group) weight(0)   vs   seldecomp y i.x1 [pw=wt], by(group) focalgroup(1)"
oaxaca y d1_2 d1_3 [pw=wt], by(group) weight(0)
assert real("`e(group_1)'") == 0
matrix ob = e(b)
seldecomp y i.x1 [pw=wt], by(group) focalgroup(1)
assert "`e(focal_group)'" == "1" & "`e(ref_group)'" == "0"
assert "`e(diffsign)'" == "focal_group - ref_group"
matrix sb = e(b)
assert_neg ob sb `tol'

* full interaction
di as txt _n "---- (vi b) oaxaca y dc_2-dc_6 [pw=wt], by(group) weight(0)   vs   seldecomp y x1##x2 [pw=wt], by(group) focalgroup(1)"
oaxaca y dc_2-dc_6 [pw=wt], by(group) weight(0)
matrix ob = e(b)
seldecomp y x1##x2 [pw=wt], by(group) focalgroup(1)
matrix sb = e(b)
assert_neg ob sb `tol'

* oaxaca swap weight(0): oaxaca's group 1 is group 1 there; the Focal Group is group 0 here, raked to group 1
di as txt _n "---- (vi c) oaxaca y d1_2 d1_3 [pw=wt], by(group) weight(0) swap   vs   seldecomp y i.x1 [pw=wt], by(group)"
oaxaca y d1_2 d1_3 [pw=wt], by(group) weight(0) swap
assert real("`e(group_1)'") == 1
matrix ob = e(b)
seldecomp y i.x1 [pw=wt], by(group)
assert "`e(focal_group)'" == "0" & "`e(ref_group)'" == "1"
matrix sb = e(b)
assert_neg ob sb `tol'

di as res _n "PASSED: test_oaxaca.do"

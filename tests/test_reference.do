* test_reference.do -- seldecomp point estimates against the Python reference
*
* Every scalar exp_<case>_<quantity> in expected_values.do was produced by
* tests/reference/seldecomp_ref.py from tests/data/seldecomp_test.csv (see the
* main() function of that script for the exact definition of every case).  This
* file runs the matching seldecomp command for each case and asserts equality:
* one raking (cases A-E, H, I), the stages of sequentially, whose blocks not yet
* added are free by default (cases F1-F3), and of sequentially fixothermargins,
* which holds them at the Focal Group's own shares (cases G1-G3), and a level found only
* in the Focal Group under sequentially fixothermargins (cases K and L, expected_support.do).
* It also checks the stored results of each mode (e(stagetype), e(fixothermargins),
* e(rakeopts)) and that fixothermargins alone, and sequentially with individually,
* are refused with their messages.
*
* Run from the tests directory:   . do test_reference.do

version 16
clear
set more off
set varabbrev off
discard             // drop any seldecomp already in memory; the repository copy is what gets tested

* ---- tolerances ---------------------------------------------------------------
local tol  1e-6     // reldif() tolerance for means, differences, shares, percentages
local wtol 1e-5     // reldif() tolerance for individual raked-weight extremes
                    // (a single weight is more sensitive to the raking's tolerance()
                    // than a weighted mean is)
local itol 1e-10    // identities that must hold to floating-point accuracy

* ---- helpers (redefined in every test file so each file runs standalone) ------
capture program drop assert_close
program define assert_close
    // assert_close exp1 exp2 [tol] ["label"]  --  asserts reldif(exp1, exp2) <= tol
    // exp1 and exp2 are expressions without spaces, e.g. e(explained) exp_A_explained
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
    // assert_mat matname row col expected [tol]
    //   row and col are stripe names or numbers, e.g.  assert_mat b 1 overall:explained exp_A_explained
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

* ---- locate seldecomp.ado (see run_tests.do) ----------------------------------
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
confirm numeric variable id group y ybin x1 x2 x3 wt cluster strata
assert _N == 2000
do expected_values.do

count if group == 0
local N0 = r(N)
count if group == 1
local N1 = r(N)
assert `N0' + `N1' == 2000


* ==============================================================================
* Case A: one categorical raking variable, unweighted, focal group = 0 (default)
* ==============================================================================
di as txt _n "---- case A: seldecomp y i.x1, by(group)"
seldecomp y i.x1, by(group)

assert "`e(cmd)'"        == "seldecomp"
assert "`e(depvar)'"     == "y"
assert "`e(by)'"         == "group"
assert "`e(focal_group)'"    == "0"
assert "`e(ref_group)'"    == "1"
assert "`e(rakevars)'"   == "x1"
assert "`e(blocks)'"     == "x1"
assert "`e(stagenames)'" == "overall"
assert "`e(wtype)'"      == ""
assert "`e(fixothermargins)'" == ""
assert "`e(diffsign)'"   == "focal_group - ref_group"
assert "`e(rakeopts)'"   == "iterate(10000) tolerance(1e-8) ctrltolerance(1e-6)"
assert "`e(properties)'" == "b"
assert "`e(title)'"      != ""
assert e(N)        == 2000
assert e(N_focal)      == `N0'
assert e(N_ref)      == `N1'
assert e(N_stages) == 1
assert e(N_terms)  == 1
assert e(N_blocks) == 1
assert e(level)    == 95
assert e(converged)   == exp_A_converged
assert e(badcontrols) == 0
assert e(N_focal_zero)    == 0
assert e(N_raked)     == `N0'
assert e(N_notes)     == 0
assert "`e(stagetype)'" == ""            // only stored with sequentially or individually
assert e(mincell)      == 10
assert e(N_focal_mincell)  == 0
assert e(N_at_bound)    == 0
assert e(pct_W_focal_zero) == 0
assert missing(e(explained_support)) & missing(e(explained_cs))   // only when e(N_focal_zero) > 0
assert "`e(worstvar)'" == ""             // set only when a margin was missed
assert "`e(worstcat)'" == ""
count if e(sample)
assert r(N) == 2000

assert_close e(mean_focal)        exp_A_mean_1        `tol'  "A mean_focal"
assert_close e(mean_ref)        exp_A_mean_2        `tol'  "A mean_ref"
assert_close e(mean_raked_focal)  exp_A_raked_1       `tol'  "A mean_raked_focal"
assert_close e(difference)    exp_A_difference    `tol'  "A difference"
assert_close e(explained)     exp_A_explained     `tol'  "A explained"
assert_close e(unexplained)   exp_A_unexplained   `tol'  "A unexplained"
assert_close e(pct_explained) exp_A_pct_explained `tol'  "A pct_explained"
assert_close e(W_focal)           exp_A_W1            `tol'  "A W_focal"
assert_close e(W_ref)           `N1'                `tol'  "A W_ref"
assert_close e(raked_min)     exp_A_raked_min     `wtol' "A raked_min"
assert_close e(raked_max)     exp_A_raked_max     `wtol' "A raked_max"
assert_close e(raked_mean)    1                   `tol'  "A raked_mean (raked weights sum to W_focal = N_focal)"
* Kish effective N: unweighted base weights give N_focal; raked weights cannot exceed it
assert_close e(N_eff_base)    `N0'                `tol'  "A N_eff_base (unweighted = N_focal)"
assert e(N_eff_raked) > 0 & e(N_eff_raked) <= `N0'
assert reldif(e(ess_ratio), e(N_eff_base) / e(N_eff_raked)) < `itol'
assert e(ess_ratio) >= 1
assert reldif(e(mean_focal_cs), e(mean_focal)) < `itol'               // no zero-weight observations
assert reldif(e(pct_increment), e(pct_explained)) < `itol'    // one stage: increment = explained

* identities
assert reldif(e(explained) + e(unexplained), e(difference)) < `itol'
assert reldif(e(mean_focal) - e(mean_ref), e(difference)) < `itol'
assert reldif(e(mean_focal) - e(mean_raked_focal), e(explained)) < `itol'
assert reldif(100 * e(explained) / e(difference), e(pct_explained)) < `itol'

* e(b): the overall equation
matrix b = e(b)
assert colsof(b) == 6
assert_mat b 1 overall:focal_group     exp_A_mean_1      `tol'
assert_mat b 1 overall:raked_focal     exp_A_raked_1     `tol'
assert_mat b 1 overall:ref_group     exp_A_mean_2      `tol'
assert_mat b 1 overall:difference  exp_A_difference  `tol'
assert_mat b 1 overall:explained   exp_A_explained   `tol'
assert_mat b 1 overall:unexplained exp_A_unexplained `tol'
* the same entries through _b[]
assert reldif(_b[overall:explained], e(explained)) < `itol'
assert reldif(_b[overall:difference], e(difference)) < `itol'

* e(bystage): one row named overall, columns as laid out in the ado (30 columns)
matrix S = e(bystage)
assert rowsof(S) == 1
assert colsof(S) == 30
assert rownumb(S, "overall") == 1
local exp_cols raked_focal cum_expl cum_unexpl pct_explained increment converged badcontrols       ///
    maxreldif maxctrl raked_mean raked_min raked_max raked_sd raked_cv factor_mean factor_min ///
    factor_max factor_sd factor_cv N_raked N_focal_zero N_at_bound N_eff_base N_eff_raked           ///
    ess_ratio N_focal_mincell pct_increment pct_W_focal_zero mean_focal_zero mean_focal_cs
local exp_cols : list retokenize exp_cols    // the /// continuation lines leave extra blanks
local cn : colnames S
assert "`cn'" == "`exp_cols'"
assert_mat S 1 raked_focal       exp_A_raked_1       `tol'
assert_mat S 1 cum_expl     exp_A_explained     `tol'
assert_mat S 1 cum_unexpl   exp_A_unexplained   `tol'
assert_mat S 1 pct_explained exp_A_pct_explained `tol'
assert_mat S 1 increment     exp_A_increment     `tol'
assert el(S, 1, colnumb(S, "converged")) == 1
assert el(S, 1, colnumb(S, "N_raked"))   == `N0'
assert el(S, 1, colnumb(S, "N_focal_zero"))  == 0
assert el(S, 1, colnumb(S, "N_focal_mincell")) == 0
assert el(S, 1, colnumb(S, "pct_W_focal_zero")) == 0
assert el(S, 1, colnumb(S, "mean_focal_zero")) == .          // no zero-weight observations
assert_mat S 1 mean_focal_cs   exp_A_mean_1        `tol'
assert_mat S 1 N_eff_base  `N0'                `tol'
assert_mat S 1 pct_increment exp_A_pct_explained `tol'

* replay and ereturn list must work
seldecomp
seldecomp, level(90)
seldecomp
ereturn list


* ==============================================================================
* Case B: one categorical raking variable, pweights
* ==============================================================================
di as txt _n "---- case B: seldecomp y x1 [pw=wt], by(group)"
seldecomp y x1 [pw=wt], by(group)

assert "`e(wtype)'" == "pweight"
assert strpos("`e(wexp)'", "wt") > 0
assert "`e(rakevars)'" == "x1"
assert e(N) == 2000 & e(N_focal) == `N0' & e(N_ref) == `N1'
assert e(converged) == exp_B_converged
assert_close e(mean_focal)        exp_B_mean_1        `tol'  "B mean_focal"
assert_close e(mean_ref)        exp_B_mean_2        `tol'  "B mean_ref"
assert_close e(mean_raked_focal)  exp_B_raked_1       `tol'  "B mean_raked_focal"
assert_close e(difference)    exp_B_difference    `tol'  "B difference"
assert_close e(explained)     exp_B_explained     `tol'  "B explained"
assert_close e(unexplained)   exp_B_unexplained   `tol'  "B unexplained"
assert_close e(pct_explained) exp_B_pct_explained `tol'  "B pct_explained"
assert_close e(W_focal)           exp_B_W1            `tol'  "B W_focal"
assert_close e(raked_min)     exp_B_raked_min     `wtol' "B raked_min"
assert_close e(raked_max)     exp_B_raked_max     `wtol' "B raked_max"
matrix b = e(b)
assert_mat b 1 overall:focal_group     exp_B_mean_1      `tol'
assert_mat b 1 overall:raked_focal     exp_B_raked_1     `tol'
assert_mat b 1 overall:ref_group     exp_B_mean_2      `tol'
assert_mat b 1 overall:difference  exp_B_difference  `tol'
assert_mat b 1 overall:explained   exp_B_explained   `tol'
assert_mat b 1 overall:unexplained exp_B_unexplained `tol'
matrix bB = e(b)

* aweights and iweights give the same point estimates as pweights
seldecomp y x1 [aw=wt], by(group)
assert "`e(wtype)'" == "aweight"
matrix bB2 = e(b)
assert mreldif(bB, bB2) < `itol'
seldecomp y x1 [iw=wt], by(group)
assert "`e(wtype)'" == "iweight"
matrix bB2 = e(b)
assert mreldif(bB, bB2) < `itol'


* ==============================================================================
* Case C: two categorical raking variables (marginals), pweights
* ==============================================================================
di as txt _n "---- case C: seldecomp y x1 x2 [pw=wt], by(group)"
seldecomp y x1 x2 [pw=wt], by(group)

assert "`e(rakevars)'" == "x1 x2"
assert e(N_terms) == 2 & e(N_blocks) == 2 & e(N_stages) == 1
assert e(converged) == exp_C_converged
assert e(badcontrols) == 0
assert "`e(worstvar)'" == ""             // all margins met: no worst-fitting control is reported
assert "`e(worstcat)'" == ""
assert_close e(mean_focal)        exp_C_mean_1        `tol'  "C mean_focal"
assert_close e(mean_ref)        exp_C_mean_2        `tol'  "C mean_ref"
assert_close e(mean_raked_focal)  exp_C_raked_1       `tol'  "C mean_raked_focal"
assert_close e(difference)    exp_C_difference    `tol'  "C difference"
assert_close e(explained)     exp_C_explained     `tol'  "C explained"
assert_close e(unexplained)   exp_C_unexplained   `tol'  "C unexplained"
assert_close e(pct_explained) exp_C_pct_explained `tol'  "C pct_explained"
assert_close e(W_focal)           exp_C_W1            `tol'  "C W_focal"
assert_close e(raked_min)     exp_C_raked_min     `wtol' "C raked_min"
assert_close e(raked_max)     exp_C_raked_max     `wtol' "C raked_max"
matrix b = e(b)
assert_mat b 1 overall:explained   exp_C_explained   `tol'
assert_mat b 1 overall:unexplained exp_C_unexplained `tol'
* i.x1 i.x2 is the same specification
matrix bC = e(b)
seldecomp y i.x1 i.x2 [pw=wt], by(group)
matrix bC2 = e(b)
assert mreldif(bC, bC2) < `itol'


* ==============================================================================
* Case D: full interaction x1##x2 (one joint raking variable), pweights
* ==============================================================================
di as txt _n "---- case D: seldecomp y x1##x2 [pw=wt], by(group)"
seldecomp y x1##x2 [pw=wt], by(group)

assert "`e(rakevars)'" == "x1#x2"
assert e(N_terms) == 1
assert e(converged) == exp_D_converged
assert_close e(mean_focal)        exp_D_mean_1        `tol'  "D mean_focal"
assert_close e(mean_ref)        exp_D_mean_2        `tol'  "D mean_ref"
assert_close e(mean_raked_focal)  exp_D_raked_1       `tol'  "D mean_raked_focal"
assert_close e(difference)    exp_D_difference    `tol'  "D difference"
assert_close e(explained)     exp_D_explained     `tol'  "D explained"
assert_close e(unexplained)   exp_D_unexplained   `tol'  "D unexplained"
assert_close e(pct_explained) exp_D_pct_explained `tol'  "D pct_explained"
assert_close e(raked_min)     exp_D_raked_min     `wtol' "D raked_min"
assert_close e(raked_max)     exp_D_raked_max     `wtol' "D raked_max"
* the joint variable has 3 x 2 = 6 levels (levels are separated by | in e(levels1))
local L `"`e(levels1)'"'
assert strlen(`"`L'"') - strlen(subinstr(`"`L'"', "|", "", .)) == 5
matrix bD1 = e(b)

* x1#x2 must be identical to x1##x2 (raking to the joint distribution matches
* all lower-order marginals automatically)
seldecomp y x1#x2 [pw=wt], by(group)
assert "`e(rakevars)'" == "x1#x2"
matrix bD2 = e(b)
assert mreldif(bD1, bD2) < `itol'
* and so is i.x1##i.x2
seldecomp y i.x1##i.x2 [pw=wt], by(group)
assert "`e(rakevars)'" == "x1#x2"
matrix bD3 = e(b)
assert mreldif(bD1, bD3) < `itol'


* ==============================================================================
* Case E: focal group = 1, string raking variable x4 plus x1, pweights
* ==============================================================================
di as txt _n "---- case E: seldecomp y x4 x1 [pw=wt], by(group) focalgroup(1)"
* capture the output in a text log to check the string-conversion message
tempfile lg
local lgf `"`lg'.log"'
capture log close sdlog
log using `"`lgf'"', text replace name(sdlog)
seldecomp y x4 x1 [pw=wt], by(group) focalgroup(1)
log close sdlog

assert "`e(focal_group)'" == "1"
assert "`e(ref_group)'" == "0"
assert e(N_focal) == `N1'
assert e(N_ref) == `N0'
assert "`e(rakevars)'" == "x4 x1"
assert e(converged) == exp_E_converged
assert_close e(mean_focal)        exp_E_mean_1        `tol'  "E mean_focal"
assert_close e(mean_ref)        exp_E_mean_2        `tol'  "E mean_ref"
assert_close e(mean_raked_focal)  exp_E_raked_1       `tol'  "E mean_raked_focal"
assert_close e(difference)    exp_E_difference    `tol'  "E difference"
assert_close e(explained)     exp_E_explained     `tol'  "E explained"
assert_close e(unexplained)   exp_E_unexplained   `tol'  "E unexplained"
assert_close e(pct_explained) exp_E_pct_explained `tol'  "E pct_explained"
assert_close e(W_focal)           exp_E_W1            `tol'  "E W_focal"
assert_close e(raked_min)     exp_E_raked_min     `wtol' "E raked_min"
assert_close e(raked_max)     exp_E_raked_max     `wtol' "E raked_max"
matrix b = e(b)
assert_mat b 1 overall:focal_group   exp_E_mean_1    `tol'
assert_mat b 1 overall:explained exp_E_explained `tol'
* x4 levels are the four region names
assert strpos(`"`e(levels1)'"', "north") > 0 & strpos(`"`e(levels1)'"', "west") > 0
* the difference is the Focal Group - Reference Group whichever group is the Focal Group
assert "`e(diffsign)'" == "focal_group - ref_group"

* no conversion message (a string raking variable can only be categorical)
tempname fh
local found 0
file open `fh' using `"`lgf'"', read text
file read `fh' line
while r(eof) == 0 {
    if strpos(`"`macval(line)'"', "Converting string variable") local found 1
    file read `fh' line
}
file close `fh'
capture erase `"`lgf'"'
if `found' {
    di as err "case E: the message 'Converting string variable ...' is printed"
    exit 9
}


* ==============================================================================
* Cases F1-F3: sequentially x1, x2, x4 (the default: stage k rakes on blocks 1..k only,
* and the blocks not yet added are free)
* ==============================================================================
di as txt _n "---- cases F1-F3: seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially"
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially

assert "`e(stagetype)'" == "sequentially"
assert "`e(fixothermargins)'" == ""
assert e(N_stages) == 3 & e(N_blocks) == 3 & e(N_terms) == 3
assert "`e(blocks)'" == "x1 x2 x4"
assert "`e(stage1)'" == "x1"
assert "`e(stage2)'" == "x1 x2"
assert "`e(stage3)'" == "x1 x2 x4"
assert "`e(block1)'" == "x1" & "`e(block2)'" == "x2" & "`e(block3)'" == "x4"
assert "`e(stagenames)'" == "x1 x2 x4"
assert e(k_eq) == 5

matrix S = e(bystage)
assert rowsof(S) == 3 & colsof(S) == 30
assert rownumb(S, "x1") == 1 & rownumb(S, "x2") == 2 & rownumb(S, "x4") == 3
forvalues s = 1/3 {
    assert_mat S `s' raked_focal       exp_F`s'_raked_1       `tol'
    assert_mat S `s' cum_expl     exp_F`s'_explained     `tol'
    assert_mat S `s' cum_unexpl   exp_F`s'_unexplained   `tol'
    assert_mat S `s' pct_explained exp_F`s'_pct_explained `tol'
    assert_mat S `s' increment     exp_F`s'_increment     `tol'
    assert_mat S `s' raked_min     exp_F`s'_raked_min     `wtol'
    assert_mat S `s' raked_max     exp_F`s'_raked_max     `wtol'
    assert el(S, `s', colnumb(S, "converged")) == exp_F`s'_converged
    assert el(S, `s', colnumb(S, "N_raked"))   == `N0'
    assert el(S, `s', colnumb(S, "N_focal_zero"))  == 0
    assert_mat S `s' pct_increment `=100 * exp_F`s'_increment / exp_F`s'_difference' `tol'
}
* e(pct_increment) is the last stage's increment as a share of the difference
assert_close e(pct_increment) `=100 * exp_F3_increment / exp_F3_difference' `tol' "F3 pct_increment"
* stage 1 = case B and stage 2 = case C (sequentially alone does not hold the
* not-yet-added blocks fixed, so a stage equals the non-sequential run on its blocks)
assert_mat S 1 cum_expl exp_B_explained `tol'
assert_mat S 2 cum_expl exp_C_explained `tol'

* overall scalars are the last stage
assert_close e(mean_raked_focal)  exp_F3_raked_1       `tol' "F3 mean_raked_focal (overall)"
assert_close e(explained)     exp_F3_explained     `tol' "F3 explained (overall)"
assert_close e(unexplained)   exp_F3_unexplained   `tol' "F3 unexplained (overall)"
assert_close e(pct_explained) exp_F3_pct_explained `tol' "F3 pct_explained (overall)"
assert_close e(raked_max)     exp_F3_raked_max     `wtol' "F3 raked_max (overall)"
assert e(converged) == exp_F3_converged

* e(b): overall + marginal explained, cumulative explained/unexplained and raked mean per stage
matrix b = e(b)
assert colsof(b) == 18
assert_mat b 1 overall:explained   exp_F3_explained   `tol'
assert_mat b 1 overall:unexplained exp_F3_unexplained `tol'
assert_mat b 1 overall:raked_focal     exp_F3_raked_1     `tol'
assert_mat b 1 marg_expl:x1        exp_F1_increment   `tol'
assert_mat b 1 marg_expl:x2        exp_F2_increment   `tol'
assert_mat b 1 marg_expl:x4        exp_F3_increment   `tol'
assert_mat b 1 cum_expl:x1         exp_F1_explained   `tol'
assert_mat b 1 cum_expl:x2         exp_F2_explained   `tol'
assert_mat b 1 cum_expl:x4         exp_F3_explained   `tol'
assert_mat b 1 cum_unexpl:x1       exp_F1_unexplained `tol'
assert_mat b 1 cum_unexpl:x2       exp_F2_unexplained `tol'
assert_mat b 1 cum_unexpl:x4       exp_F3_unexplained `tol'
assert_mat b 1 raked_focal:x1      exp_F1_raked_1     `tol'
assert_mat b 1 raked_focal:x2      exp_F2_raked_1     `tol'
assert_mat b 1 raked_focal:x4      exp_F3_raked_1     `tol'
* increments sum to the overall explained part
assert reldif(_b[marg_expl:x1] + _b[marg_expl:x2] + _b[marg_expl:x4], e(explained)) < `itol'
* replay with stages
seldecomp


* ==============================================================================
* Cases G1-G3: sequentially fixothermargins (the not-yet-added blocks held at the Focal Group's
* own shares at every stage: the papers' procedure)
* ==============================================================================
di as txt _n "---- cases G1-G3: seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially fixothermargins"
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially fixothermargins

assert "`e(fixothermargins)'" == "fixothermargins"
assert "`e(stagetype)'" == "sequentially"
assert e(N_stages) == 3
matrix S = e(bystage)
assert rowsof(S) == 3
forvalues s = 1/3 {
    assert_mat S `s' raked_focal       exp_G`s'_raked_1       `tol'
    assert_mat S `s' cum_expl     exp_G`s'_explained     `tol'
    assert_mat S `s' cum_unexpl   exp_G`s'_unexplained   `tol'
    assert_mat S `s' pct_explained exp_G`s'_pct_explained `tol'
    assert_mat S `s' increment     exp_G`s'_increment     `tol'
    assert_mat S `s' raked_min     exp_G`s'_raked_min     `wtol'
    assert_mat S `s' raked_max     exp_G`s'_raked_max     `wtol'
    assert el(S, `s', colnumb(S, "converged")) == exp_G`s'_converged
}
* the final stage coincides with sequentially alone (no block is left to hold)
assert_mat S 3 cum_expl exp_F3_explained `tol'
assert_close e(explained) exp_G3_explained `tol' "G3 explained (overall)"
matrix b = e(b)
assert_mat b 1 marg_expl:x1  exp_G1_increment `tol'
assert_mat b 1 marg_expl:x2  exp_G2_increment `tol'
assert_mat b 1 marg_expl:x4  exp_G3_increment `tol'
assert_mat b 1 cum_expl:x4   exp_G3_explained `tol'
assert_mat b 1 raked_focal:x1 exp_G1_raked_1  `tol'
assert_mat b 1 raked_focal:x2 exp_G2_raked_1  `tol'
assert reldif(_b[marg_expl:x1] + _b[marg_expl:x2] + _b[marg_expl:x4], e(explained)) < `itol'
* fixothermargins holds the blocks not yet added: stage 1 differs from sequentially alone
* (case F1, x2 and x4 free), and its x2 share in the Focal Group is the Focal Group's own
assert reldif(el(S, 1, colnumb(S, "cum_expl")), exp_F1_explained) > 1e-4
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially fixothermargins generate(rwg)
quietly summarize x2 [aw=wt] if group == 0, meanonly
local x2own = r(mean)
quietly summarize x2 [aw=rwg_s1] if group == 0, meanonly
assert reldif(r(mean), `x2own') < 1e-5
* and sequentially alone lets it drift at stage 1
seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially generate(rwf)
quietly summarize x2 [aw=rwf_s1] if group == 0, meanonly
assert reldif(r(mean), `x2own') > 1e-4
drop rwg* rwf*

* fixothermargins needs a staged mode, and the two modes exclude each other (error 198,
* with the ado's messages)
tempfile lg2
local lgf2 `"`lg2'.log"'
foreach c in "fixothermargins" "sequentially individually" {
    capture log close sdlog
    log using `"`lgf2'"', text replace name(sdlog)
    capture noisily seldecomp y x1 x2 x4 [pw=wt], by(group) `c'
    local rc = _rc
    log close sdlog
    assert `rc' == 198
    if "`c'" == "fixothermargins" local msg "fixothermargins requires sequentially or individually"
    else local msg "sequentially and individually may not be combined"
    tempname fh2
    local found 0
    file open `fh2' using `"`lgf2'"', read text
    file read `fh2' line
    while r(eof) == 0 {
        if strpos(`"`macval(line)'"', "`msg'") local found 1
        file read `fh2' line
    }
    file close `fh2'
    capture erase `"`lgf2'"'
    if !`found' {
        di as err "`c': the message '`msg'' is not printed"
        exit 9
    }
}


* ==============================================================================
* Cases H, H2, H3: continuous x3 in quartiles (pooled / focal / reference cut points)
* ==============================================================================
di as txt _n "---- case H: seldecomp y c.x3 x2 [pw=wt], by(group) contvars(x3 4)"
seldecomp y c.x3 x2 [pw=wt], by(group) contvars(x3 4)

assert "`e(contvars)'" == "x3 4"
assert "`e(rakevars)'" == "x3 x2"
assert e(converged) == exp_H_converged
assert_close e(mean_focal)        exp_H_mean_1        `tol'  "H mean_focal"
assert_close e(mean_ref)        exp_H_mean_2        `tol'  "H mean_ref"
assert_close e(mean_raked_focal)  exp_H_raked_1       `tol'  "H mean_raked_focal"
assert_close e(difference)    exp_H_difference    `tol'  "H difference"
assert_close e(explained)     exp_H_explained     `tol'  "H explained"
assert_close e(unexplained)   exp_H_unexplained   `tol'  "H unexplained"
assert_close e(pct_explained) exp_H_pct_explained `tol'  "H pct_explained"
assert_close e(raked_min)     exp_H_raked_min     `wtol' "H raked_min"
assert_close e(raked_max)     exp_H_raked_max     `wtol' "H raked_max"
* the quartile cut points appear (formatted %9.0g) in the level labels of term 1
* (4 quantile groups -> 3 cut points -> 4 levels separated by |)
local L `"`e(levels1)'"'
assert strlen(`"`L'"') - strlen(subinstr(`"`L'"', "|", "", .)) == 3
forvalues j = 1/3 {
    local cs = strtrim("`: di %9.0g exp_H_cut`j''")
    if strpos(`"`L'"', "`cs'") == 0 {
        di as err `"case H: cut point `j' (`cs') not found in e(levels1) = `L'"'
        exit 9
    }
}
* the cut points are also stored numerically, with the binning sample and the number of groups
assert "`e(cutvars)'" == "x3"
assert "`e(cutq)'" == "4"
assert "`e(contsample)'" == "pooled"
matrix CUT = e(cuts1)
assert rowsof(CUT) == 1 & colsof(CUT) == 3
local cutn : colnames CUT
assert "`cutn'" == "cut1 cut2 cut3"
forvalues j = 1/3 {
    assert_mat CUT 1 `j' exp_H_cut`j' `tol'
}
matrix bH = e(b)
* pooled is the default binning sample
seldecomp y c.x3 x2 [pw=wt], by(group) contvars(x3 4, pooled)
matrix bHp = e(b)
assert mreldif(bH, bHp) < `itol'
* the c. prefix is accepted inside contvars()
seldecomp y c.x3 x2 [pw=wt], by(group) contvars(c.x3 4)
matrix bHc = e(b)
assert mreldif(bH, bHc) < `itol'

di as txt _n "---- case H2: contvars(x3 4, focalgroup)"
seldecomp y c.x3 x2 [pw=wt], by(group) contvars(x3 4, focalgroup)
assert e(converged) == exp_H2_converged
assert "`e(contsample)'" == "focalgroup"
matrix CUT = e(cuts1)
forvalues j = 1/3 {
    assert_mat CUT 1 `j' exp_H2_cut`j' `tol'
}
assert_close e(mean_raked_focal)  exp_H2_raked_1       `tol'  "H2 mean_raked_focal"
assert_close e(difference)    exp_H2_difference    `tol'  "H2 difference"
assert_close e(explained)     exp_H2_explained     `tol'  "H2 explained"
assert_close e(unexplained)   exp_H2_unexplained   `tol'  "H2 unexplained"
assert_close e(pct_explained) exp_H2_pct_explained `tol'  "H2 pct_explained"
assert_close e(raked_min)     exp_H2_raked_min     `wtol' "H2 raked_min"
assert_close e(raked_max)     exp_H2_raked_max     `wtol' "H2 raked_max"
local L `"`e(levels1)'"'
forvalues j = 1/3 {
    local cs = strtrim("`: di %9.0g exp_H2_cut`j''")
    if strpos(`"`L'"', "`cs'") == 0 {
        di as err `"case H2: cut point `j' (`cs') not found in e(levels1) = `L'"'
        exit 9
    }
}

di as txt _n "---- case H3: contvars(x3 4, refgroup)"
seldecomp y c.x3 x2 [pw=wt], by(group) contvars(x3 4, refgroup)
assert e(converged) == exp_H3_converged
assert "`e(contsample)'" == "refgroup"
matrix CUT = e(cuts1)
forvalues j = 1/3 {
    assert_mat CUT 1 `j' exp_H3_cut`j' `tol'
}
assert_close e(mean_raked_focal)  exp_H3_raked_1       `tol'  "H3 mean_raked_focal"
assert_close e(difference)    exp_H3_difference    `tol'  "H3 difference"
assert_close e(explained)     exp_H3_explained     `tol'  "H3 explained"
assert_close e(unexplained)   exp_H3_unexplained   `tol'  "H3 unexplained"
assert_close e(pct_explained) exp_H3_pct_explained `tol'  "H3 pct_explained"
assert_close e(raked_min)     exp_H3_raked_min     `wtol' "H3 raked_min"
assert_close e(raked_max)     exp_H3_raked_max     `wtol' "H3 raked_max"
local L `"`e(levels1)'"'
forvalues j = 1/3 {
    local cs = strtrim("`: di %9.0g exp_H3_cut`j''")
    if strpos(`"`L'"', "`cs'") == 0 {
        di as err `"case H3: cut point `j' (`cs') not found in e(levels1) = `L'"'
        exit 9
    }
}


* ==============================================================================
* Case I: binary outcome, unweighted, x1 x2 x4
* ==============================================================================
di as txt _n "---- case I: seldecomp ybin x1 x2 x4, by(group)"
seldecomp ybin x1 x2 x4, by(group)

assert "`e(depvar)'" == "ybin"
assert "`e(rakevars)'" == "x1 x2 x4"
assert e(N) == 2000 & e(N_focal) == `N0' & e(N_ref) == `N1'
assert e(converged) == exp_I_converged
assert_close e(mean_focal)        exp_I_mean_1        `tol'  "I mean_focal"
assert_close e(mean_ref)        exp_I_mean_2        `tol'  "I mean_ref"
assert_close e(mean_raked_focal)  exp_I_raked_1       `tol'  "I mean_raked_focal"
assert_close e(difference)    exp_I_difference    `tol'  "I difference"
assert_close e(explained)     exp_I_explained     `tol'  "I explained"
assert_close e(unexplained)   exp_I_unexplained   `tol'  "I unexplained"
assert_close e(pct_explained) exp_I_pct_explained `tol'  "I pct_explained"
assert_close e(W_focal)           exp_I_W1            `tol'  "I W_focal"
assert_close e(raked_min)     exp_I_raked_min     `wtol' "I raked_min"
assert_close e(raked_max)     exp_I_raked_max     `wtol' "I raked_max"
matrix b = e(b)
assert_mat b 1 overall:focal_group     exp_I_mean_1      `tol'
assert_mat b 1 overall:explained   exp_I_explained   `tol'
assert_mat b 1 overall:unexplained exp_I_unexplained `tol'
* the raked mean of a 0/1 outcome is a probability
assert e(mean_raked_focal) > 0 & e(mean_raked_focal) < 1
* nothing about binning is stored without c. variables
assert "`e(cutvars)'" == "" & "`e(cutq)'" == "" & "`e(contsample)'" == ""
capture confirm matrix e(cuts1)
assert _rc != 0

* ==============================================================================
* Cases K and L: a level found only in the Focal Group under sequentially fixothermargins
* (expected_support.do). The
* level is zeroed at the stage that rakes its variable; the held blocks' shares are those of
* the observations that keep a positive weight at that stage (K: x1z first, so stages 1 and
* 2 hold x2 and x4 without level 9's observations; L: x1z second, so stage 1 holds x1z at the
* whole Focal Group's shares and stage 2 holds x4 without them)
* ==============================================================================
di as txt _n "---- cases K and L: a level found only in the Focal Group, sequentially fixothermargins"
do expected_support.do
gen byte x1z = x1
replace x1z = 9 if group == 0 & id <= 60
foreach case in K L {
    if "`case'" == "K" seldecomp y x1z x2 x4 [pw=wt], by(group) sequentially fixothermargins
    else seldecomp y x2 x1z x4 [pw=wt], by(group) sequentially fixothermargins
    matrix S = e(bystage)
    forvalues s = 1/3 {
        assert_close S[`s',1] exp_`case'`s'_raked_1     `tol' "`case'`s' raked_focal"
        assert_close S[`s',2] exp_`case'`s'_explained   `tol' "`case'`s' explained"
        assert_close S[`s',3] exp_`case'`s'_unexplained `tol' "`case'`s' unexplained"
        assert_close S[`s',5] exp_`case'`s'_increment   `tol' "`case'`s' increment"
        assert S[`s',6] == exp_`case'`s'_converged
    }
    assert e(N_focal_zero) == 30
    assert S[1,21] == cond("`case'" == "K", 30, 0)
}
drop x1z

di as res _n "PASSED: test_reference.do"

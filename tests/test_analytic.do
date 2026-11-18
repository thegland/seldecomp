* test_analytic.do -- the worked example of Pylypchuk and Selden (2008), footnote 18
* (p. 1114), against its closed-form solution.  With two binary raking variables the
* raked joint distribution is the unique 2x2 table with the target margins and the focal
* group's odds ratio (iterative proportional fitting preserves the odds ratio), so every
* stage of both modes of sequentially can be computed by hand, without the Python
* reference: sequentially fixothermargins (the papers' procedure: the block not yet added
* is held at the Focal Group's own margin) and sequentially alone (the default: it is left free
* and drifts, as footnote 18 warns).  The difference is the Focal Group - Reference Group.
*
* Run from the tests directory:   . do test_analytic.do

version 16
clear
set more off
set varabbrev off
discard             // drop any seldecomp already in memory; the repository copy is what gets tested

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

* ---- the footnote-18 populations ----------------------------------------------
* minority (Focal Group, raked): POOR mean .4, IMMIG mean .5; joint cells 300/100/200/400 of 1000
* white (Reference Group): POOR mean .2, IMMIG mean .1; joint cells 20/180/80/720 of 1000
* outcome: a coverage indicator whose cell means are exact by construction
input byte white byte poor byte immig int n double ybar
0 1 1 300 .40
0 1 0 100 .55
0 0 1 200 .60
0 0 0 400 .80
1 1 1  20 .50
1 1 0 180 .65
1 0 1  80 .70
1 0 0 720 .90
end
expand n
bysort white poor immig: gen int k = _n
gen byte covered = k <= round(ybar * n)
label define wl 0 "Minority" 1 "White"
label values white wl

* ---- closed form -----------------------------------------------------------------
* p11 = a solves a * p00 = OR * p10 * p01 with p10 = mp - a, p01 = mi - a,
* p00 = 1 - mp - mi + a, and OR = (.3 * .4) / (.1 * .2) = 6 in the minority sample
mata:
real scalar orsolve(real scalar mp, real scalar mi, real scalar OR)
{
    real scalar A, B, C, r1, r2
    A = 1 - OR
    B = (1 - mp - mi) + OR * (mp + mi)
    C = -OR * mp * mi
    r1 = (-B + sqrt(B^2 - 4*A*C)) / (2*A)
    r2 = (-B - sqrt(B^2 - 4*A*C)) / (2*A)
    return(r1 > 0 & r1 < min((mp, mi)) ? r1 : r2)
}
real scalar ymean(real scalar mp, real scalar mi, real scalar a)
{
    // minority cell means: (1,1) .40, (1,0) .55, (0,1) .60, (0,0) .80
    return(a*.40 + (mp - a)*.55 + (mi - a)*.60 + (1 - mp - mi + a)*.80)
}
m1 = ymean(.4, .5, .3)
// fixothermargins stage 1: POOR -> .2 with IMMIG held at .5
a1 = orsolve(.2, .5, 6); h1 = ymean(.2, .5, a1)
// both variants end at the white margins (.2, .1)
a2 = orsolve(.2, .1, 6); h2 = ymean(.2, .1, a2)
// sequentially alone, stage 1: rake on POOR only, cells scaled by .5 and 4/3
c1 = .3*.5*.40 + .1*.5*.55 + .2*(4/3)*.60 + .4*(4/3)*.80
st_numscalar("cf_m1", m1); st_numscalar("cf_h1", h1); st_numscalar("cf_h2", h2); st_numscalar("cf_c1", c1)
end
scalar cf_m2 = (20*.50 + 180*.65 + 80*.70 + 720*.90) / 1000

* ---- sequentially fixothermargins: the papers' procedure -------------------------
seldecomp covered (poor) (immig), by(white) sequentially fixothermargins generate(rwh) tolerance(1e-12)
assert "`e(stagetype)'" == "sequentially" & "`e(fixothermargins)'" == "fixothermargins"
assert "`e(rakeopts)'" == "iterate(10000) tolerance(1e-12) ctrltolerance(1e-6)"
matrix Sh = e(bystage)
assert reldif(e(mean_focal), cf_m1) < 1e-12
assert reldif(e(mean_ref), cf_m2) < 1e-12
* the sign is always Focal Group (Minority) - Reference Group (White); the papers report the
* opposite sign (White - Minority): negate the three parts for their tables
assert "`e(diffsign)'" == "focal_group - ref_group"
assert reldif(e(difference), cf_m1 - cf_m2) < 1e-12
assert reldif(Sh[1, 1], cf_h1) < 1e-9
assert reldif(Sh[2, 1], cf_h2) < 1e-9
assert reldif(e(explained), cf_m1 - cf_h2) < 1e-9
assert reldif(e(unexplained), cf_h2 - cf_m2) < 1e-9
* the stage-1 weights reproduce the footnote's intermediate margins (.2, .5)
sum poor [aw=rwh_s1] if white == 0, meanonly
assert reldif(r(mean), .2) < 1e-9
sum immig [aw=rwh_s1] if white == 0, meanonly
assert reldif(r(mean), .5) < 1e-9

* ---- sequentially alone (the default: the other block free): IMMIG drifts at stage 1
* (4/3 * .2 + .5 * .3 = .41667), as fn 18 warns
seldecomp covered (poor) (immig), by(white) sequentially generate(rwc) tolerance(1e-12)
assert "`e(stagetype)'" == "sequentially" & "`e(fixothermargins)'" == ""
matrix Sc = e(bystage)
assert reldif(Sc[1, 1], cf_c1) < 1e-9
assert reldif(Sc[2, 1], cf_h2) < 1e-9
sum immig [aw=rwc_s1] if white == 0, meanonly
assert reldif(r(mean), 4/3*.2 + .5*.3) < 1e-9

di as res _n "PASSED: test_analytic.do"

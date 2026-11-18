* test_svy.do -- the svy replication prefixes on a synthetic stratified design, without
* network data: svy brr against a hand-computed BRR variance; the papers' route (svy brr with
* fay() on sequentially fixothermargins: e(stagetype), e(fixothermargins), and the difference
* always Focal Group - Reference Group); svy jackknife with
* subpop(), replay after each, and the accounting of a replicate in which raking fails; the
* prefix safety (bayesboot() and makeplot refused under a prefix,
* generate() under a prefix, the pre-flight check and the messages on failed replicates,
* margins kept by force, also those that trimming left unmet). The last block (the pre-flight
* check on nhanes2brr) uses the network when it is available and is skipped otherwise.
*
* Run from the tests directory:   . do test_svy.do

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

* ---- design: 16 strata x 2 PSUs x 12 observations, Hadamard BRR half-samples ----
set seed 20260925
set obs 384
gen int stratum = ceil(_n / 24)
gen byte psu = 1 + mod(ceil(_n / 12) - 1, 2)
gen int psuid = stratum * 10 + psu
gen double wt = 50 + 100 * runiform()
gen byte group = runiform() < .5
gen byte x1 = 1 + (runiform() < .4) + (runiform() < .3 + .3 * group)
gen byte x2 = runiform() < .6 - .2 * group
gen double y = 2 + x1 + x2 + .8 * group + rnormal()
* Hadamard matrix of order 16: replicate r keeps PSU 1 of stratum h when H[r, h] == 1
mata: Hd = (1); for (i = 1; i <= 4; i++) Hd = (Hd, Hd \ Hd, -Hd); st_matrix("HAD", Hd)
forvalues r = 1/16 {
    gen double brr`r' = 2 * wt * ((HAD[`r', stratum] == 1) == (psu == 1))
}
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr16)

* ---- point estimate is unaffected by the prefix ----------------------------------
seldecomp y i.x1 i.x2 [pw=wt], by(group)
matrix b_plain = e(b)
svy brr, mse: seldecomp y i.x1 i.x2, by(group)
assert "`e(vce)'" == "brr"
assert e(N_reps) == 16 & e(N_misreps) == 0
assert mreldif(e(b), b_plain) < 1e-10
matrix V_brr = e(V)
seldecomp                                   // replay after svy brr
assert "`e(cmd)'" == "seldecomp"

* ---- svy brr, mse against a hand-computed BRR variance: V = sum_r (b_r - b)^2 / R -----
* (each replicate is seldecomp with the BRR half-sample weights as pweights)
local k = colsof(b_plain)
matrix acc = J(1, `k', 0)
forvalues r = 1/16 {
    qui seldecomp y i.x1 i.x2 [pw=brr`r'], by(group)
    matrix d = e(b) - b_plain
    forvalues j = 1/`k' {
        matrix acc[1, `j'] = acc[1, `j'] + d[1, `j']^2
    }
}
forvalues j = 1/`k' {
    assert reldif(sqrt(V_brr[`j', `j']), sqrt(acc[1, `j'] / 16)) < 1e-6
}

* ---- the papers' inference route: BRR with Fay adjustment, sequentially fixothermargins --
svy brr, hadamard(HAD) fay(.5): seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins
assert e(N_reps) == 16 & colsof(e(b)) == 14
assert "`e(stagetype)'" == "sequentially" & "`e(fixothermargins)'" == "fixothermargins"
assert "`e(diffsign)'" == "focal_group - ref_group"
seldecomp
seldecomp, makeplot(, twoway(name(svyplot, replace)))
graph describe svyplot
graph drop _all

* ---- svy jackknife with subpop() ----------------------------------------------------
svy jackknife, subpop(if x2 == 1): seldecomp y i.x1, by(group)
assert "`e(vce)'" == "jackknife"
assert e(N_sub) > 0 & e(N_misreps) == 0
seldecomp

* ---- a replicate that empties a Focal Group level fails and is counted -----------------
* level 4 of x1 occurs in the Focal Group only in PSU 1 of stratum 2 (column 2 of the Hadamard matrix
* alternates, column 1 is all +1): the 8 half-samples without that PSU cannot be raked
* (support failure), so svy drops them
gen byte x3 = x1
replace x3 = 4 if group == 1 & mod(_n, 7) == 0          // Reference Group: level 4 in many PSUs
replace x3 = 4 if group == 0 & stratum == 2 & psu == 1 & mod(_n, 3) == 0   // Focal Group: one PSU
svy brr, mse: seldecomp y i.x3, by(group)
assert e(N_misreps) == 8 & e(N_reps) == 8

* ==============================================================================
* Prefix safety, on this design
* ==============================================================================
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
capture program drop log_count
program define log_count, rclass
    // log_count "file" "text"  --  r(n), the number of lines of the text log with the text
    args file text
    mata: st_local("n", strofreal(sum(strpos(cat(st_local("file")), st_local("text")) :> 0)))
    return scalar n = `n'
end

* ---- the preamble before svy's table ----------------------------------------------------
* svy shows its own coefficient table, with e(b)'s names. Before it seldecomp prints, once
* each: the two groups, the staging, the order of the difference, the raking status and the
* pointer to its own tables; not the weights and inference lines of its header, which svy's
* header gives; and nothing under quietly
tempfile pre
quietly log using `"`pre'"', text replace name(prelog)
svy brr, hadamard(HAD) fay(.5): seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins
quietly log close prelog
foreach t in "Focal Group (raked): group = 0" "Reference Group: group = 1"            ///
    "Blocks added sequentially (others held at focal group margins):"                  ///
    "Difference = focal - reference. Explained = focal - raked focal."                 ///
    "Raking (all stages): weight iterations converged; all targets met."               ///
    "(type seldecomp to redisplay these results in seldecomp's own tables)" {
    log_has `"`pre'"' "`t'"
    assert r(count) == 1
}
foreach t in "Weights: " "Inference: " {
    log_has `"`pre'"' "`t'"
    assert r(found) == 0
}
quietly log using `"`pre'"', text replace name(prelog)
quietly svy brr, hadamard(HAD) fay(.5): seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins
quietly log close prelog
foreach t in "Focal Group (raked)" "Raking (all stages)" "type seldecomp to redisplay" {
    log_has `"`pre'"' "`t'"
    assert r(found) == 0
}

capture program drop s1_refused
program define s1_refused
    // s1_refused "<prefix>" "<its name>" "<log>": bayesboot() and makeplot under the prefix
    // are refused (198) in its full-sample call, before any work: the messages, no graph, no
    // Bayesian bootstrap run
    args pfx pname lg
    capture graph drop _all
    quietly log using "`lg'", text replace name(s1log)
    capture noisily `pfx': seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins bayesboot(reps(3))
    local rc1 = _rc
    capture noisily `pfx': seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins makeplot(marginal)
    local rc2 = _rc
    capture noisily `pfx': seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins makeplot
    local rc3 = _rc
    quietly log close s1log
    assert `rc1' == 198 & `rc2' == 198 & `rc3' == 198
    log_has "`lg'" "bayesboot() may not be combined with the `pname' prefix"
    assert r(found) == 1
    log_has "`lg'" "makeplot may not be combined with the `pname' prefix; afterward, type seldecomp, makeplot(marginal)"
    assert r(found) == 1
    log_has "`lg'" "makeplot may not be combined with the `pname' prefix; afterward, type seldecomp, makeplot"
    assert r(count) == 2
    log_has "`lg'" "Bayesian bootstrap replications"
    assert r(found) == 0
    quietly graph dir
    assert strtrim(`"`r(list)'"') == ""
end
tempfile s1lg
local lg `"`s1lg'"'

* every replication prefix
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr16)
s1_refused "svy brr" "svy brr" `"`lg'"'
s1_refused "svy brr, nodots" "svy brr" `"`lg'"'
svyset psuid [pw=wt], strata(stratum)
s1_refused "svy jackknife" "svy jackknife" `"`lg'"'
svyset psuid [pw=wt], strata(stratum) sdrweight(brr1-brr16)
s1_refused "svy sdr" "svy sdr" `"`lg'"'
svyset psuid [pw=wt], strata(stratum) bsrweight(brr1-brr16)
s1_refused "svy bootstrap" "svy bootstrap" `"`lg'"'
svyset, clear
s1_refused "bootstrap, reps(3)" "bootstrap" `"`lg'"'
s1_refused "bootstrap, reps(3) bca" "bootstrap" `"`lg'"'
s1_refused "jackknife" "jackknife" `"`lg'"'
if c(stata_version) >= 19 {
    s1_refused "bayesboot, reps(3)" "bayesboot" `"`lg'"'
}
* under a prefix the refusal comes before the plotspec is checked
* and its twoway options are tried (Plot_spec, Plot_try), and it shows the plotspec as typed:
* braces are not read as SMCL
capture graph drop _all
quietly log using "`lg'", text replace name(s1log)
capture noisily bootstrap, reps(3): seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins makeplot(contibutions)
local rc1 = _rc
capture noisily jackknife: seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins makeplot(waterfall, blocklabels(a "{bf:Demo}" b Geo))
local rc2 = _rc
quietly log close s1log
assert `rc1' == 198 & `rc2' == 198
log_has "`lg'" "makeplot may not be combined with the bootstrap prefix; afterward, type seldecomp, makeplot(contibutions)"
assert r(found) == 1
log_has "`lg'" "twoway returned error"
assert r(found) == 0
log_has "`lg'" `"afterward, type seldecomp, makeplot(waterfall, blocklabels(a "{bf:Demo}" b Geo))"'
assert r(found) == 1
quietly graph dir
assert strtrim(`"`r(list)'"') == ""
* without a prefix both run, and the replay after a prefix draws the graph
seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins makeplot(, twoway(name(s1g, replace)))
graph describe s1g
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr16)
svy brr, nodots: seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins
seldecomp, makeplot(marginal, twoway(name(s1g, replace)))
graph describe s1g
capture graph drop _all

* generate() under a prefix is made by the full-sample call only: no "already defined",
* no r(2000), and the variables hold the full-sample weights; a variable that exists before
* the run still needs replace
seldecomp y (a: x1) (b: x2) [pw=wt], by(group) sequentially fixothermargins generate(gw)
seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins generate(gu)
quietly log using `"`lg'"', text replace name(s1log)
svy brr, nodots: seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins generate(g1)
assert e(N_reps) == 16 & e(N_misreps) == 0
bootstrap, reps(4) seed(1) nodots: seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins generate(g2)
jackknife, cluster(psuid) nodots: seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins generate(g3)
assert e(N_misreps) == 0
* bca: the jackknife's full-sample call of the same run replaces the variables
bootstrap, reps(4) seed(1) bca cluster(psuid) nodots: seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins generate(g4)
if c(stata_version) >= 19 {
    // the bayesboot prefix makes a dry full-sample call first, then bootstrap's
    bayesboot, reps(3) nodots: seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins generate(g5)
}
else {
    clonevar g5 = gu
    clonevar g5_s1 = gu_s1
}
quietly log close s1log
log_has `"`lg'"' "already defined"
assert r(found) == 0
assert reldif(g1, gw) < 1e-12 & reldif(g1_s1, gw_s1) < 1e-12
foreach g in g2 g3 g4 g5 {
    assert reldif(`g', gu) < 1e-12 & reldif(`g'_s1, gu_s1) < 1e-12
}
assert "`: char g1[seldecomp_stage]'" == "2" & "`: char g1_s1[seldecomp_stage]'" == "1"
capture svy brr, nodots: seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins generate(g1)
assert _rc == 110
svy brr, nodots: seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins generate(g1, replace)
assert reldif(g1, gw) < 1e-12
drop gw* gu* g1* g2* g3* g4* g5*

* under svy the svyset weight, design and replicate-weight variables are inputs too (the
* weight the full-sample call gets is a temporary copy)
foreach v in wt brr3 psuid stratum {
    capture noisily svy brr, nodots: seldecomp y i.x1 i.x2, by(group) generate(`v', replace)
    assert _rc == 198
}
quietly log using `"`lg'"', text replace name(s1log)
capture noisily svy brr, nodots: seldecomp y i.x1 i.x2, by(group) generate(brr3, replace)
quietly log close s1log
log_has `"`lg'"' "generate(): brr3 is a variable of the svyset design; choose another name"
assert r(found) == 1
assert brr3 == 2 * wt * ((HAD[3, stratum] == 1) == (psu == 1))

* the prefix's own variables: generate(clb, replace) would write the raked weights into the
* cluster() variable of bootstrap or jackknife before the replicates read it, and into svy's
* subpop() variable. The command cannot see the prefix's options, so under a prefix replace
* overwrites only a variable that generate() made (char seldecomp_cmd)
gen long clb = psuid
gen byte sv = x2 == 1
gen double keepold = 1
gen double fpcw = 4
gen byte psb = 1 + (x2 == 1)
gen double pstotb = .
forvalues p = 1/2 {
    qui sum wt if psb == `p', meanonly
    qui replace pstotb = r(sum) if psb == `p'
}
qui datasignature
local sig0 `"`r(datasignature)'"'
quietly log using `"`lg'"', text replace name(s1log)
capture noisily bootstrap, reps(3) seed(5) cluster(clb) nodots: seldecomp y i.x1 i.x2, by(group) generate(clb, replace)
assert _rc == 198
capture noisily jackknife, cluster(clb) nodots: seldecomp y i.x1 i.x2, by(group) generate(clb, replace)
assert _rc == 198
capture noisily svy brr, nodots subpop(sv): seldecomp y i.x1, by(group) generate(sv, replace)
assert _rc == 198
capture noisily bootstrap, reps(3) seed(5) nodots: seldecomp y i.x1 i.x2, by(group) generate(keepold, replace)
assert _rc == 198
quietly log close s1log
log_has `"`lg'"' "generate(..., replace): under the bootstrap prefix, replace overwrites only variables made by generate(); clb was not"
assert r(found) == 1
log_has `"`lg'"' "generate(..., replace): under the jackknife prefix, replace overwrites only variables made by generate(); clb was not"
assert r(found) == 1
log_has `"`lg'"' "generate(..., replace): under the svy brr prefix, replace overwrites only variables made by generate(); sv was not"
assert r(found) == 1
* under svy, fpc(), poststrata()/postweight() and calibration variables are design
* variables too
svyset psuid [pw=wt], strata(stratum) fpc(fpcw)
capture noisily svy jackknife, nodots: seldecomp y i.x1 i.x2, by(group) generate(fpcw, replace)
assert _rc == 198
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr16) poststrata(psb) postweight(pstotb)
foreach v in psb pstotb {
    capture noisily svy brr, nodots: seldecomp y i.x1, by(group) generate(`v', replace)
    assert _rc == 198
}
qui sum wt, meanonly
local t0 = r(sum)
qui sum wt if x2 == 1, meanonly
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr16) rake(i.x2, totals(_cons=`t0' 1.x2=`r(sum)'))
quietly log using `"`lg'"', text replace name(s1log)
capture noisily svy brr, nodots: seldecomp y i.x1, by(group) generate(x2, replace)
assert _rc == 198
quietly log close s1log
log_has `"`lg'"' "generate(): x2 is a variable of the svyset design; choose another name"
assert r(found) == 1
qui datasignature
assert `"`r(datasignature)'"' == `"`sig0'"'
* without a prefix replace still overwrites keepold, and then, under a prefix, as a variable
* that generate() made (as g1 above)
seldecomp y i.x1 i.x2, by(group) generate(keepold, replace)
assert "`: char keepold[seldecomp_cmd]'" != ""
bootstrap, reps(3) seed(5) nodots: seldecomp y i.x1 i.x2, by(group) generate(keepold, replace)
drop clb sv keepold fpcw psb pstotb
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr16)

* the sample: variables made by an earlier full-sample call of the same prefix run
* are replaced only if its specification and estimation sample are the same. A program that
* calls seldecomp twice with one specification on two samples fails on the second call (110),
* as without a prefix; under bootstrap too, where the second sample's weights would otherwise
* replace the first's silently
capture program drop s1_two
program define s1_two, eclass
    args a b
    seldecomp y i.x1 i.x2 if `a', by(group) generate(gt)
    seldecomp y i.x1 i.x2 if `b', by(group) generate(gt)
end
capture noisily s1_two "x2 == 0 | x2 == 1" "stratum <= 8"
assert _rc == 110
drop gt
capture noisily bootstrap _b, reps(3) seed(1) nodots: s1_two "x2 == 0 | x2 == 1" "stratum <= 8"
assert _rc == 110
assert !missing(gt)
drop gt
bootstrap _b, reps(3) seed(1) nodots: s1_two "stratum <= 8" "stratum <= 8"
drop gt
* the sample is compared, not its text: bootstrap bca's jackknife and Stata 19's bayesboot
* dry run pass if/in and the weight as temporary variables; with
* missing values and zero weights the variables still equal a plain run's
gen double x1m = x1 if mod(_n, 17) != 0
gen double wz = wt * (mod(_n, 23) != 0)
seldecomp y (a: x1m) (b: x2) [pw=wz] if stratum <= 12, by(group) sequentially fixothermargins generate(hw)
seldecomp y (a: x1m) (b: x2) if stratum <= 12, by(group) sequentially fixothermargins generate(hu)
bootstrap, reps(3) seed(1) bca cluster(psuid) force nodots: seldecomp y (a: x1m) (b: x2) [pw=wz] if stratum <= 12, by(group) sequentially fixothermargins generate(h1)
if c(stata_version) >= 19 {
    bayesboot, reps(3) nodots: seldecomp y (a: x1m) (b: x2) if stratum <= 12, by(group) sequentially fixothermargins generate(h2)
}
else {
    clonevar h2 = hu
    clonevar h2_s1 = hu_s1
}
assert (missing(h1) & missing(hw)) | reldif(h1, hw) < 1e-12
assert (missing(h1_s1) & missing(hw_s1)) | reldif(h1_s1, hw_s1) < 1e-12
assert (missing(h2) & missing(hu)) | reldif(h2, hu) < 1e-12
assert (missing(h2_s1) & missing(hu_s1)) | reldif(h2_s1, hu_s1) < 1e-12
count if missing(h1) & stratum <= 12
assert r(N) > 0
drop x1m wz hw* hu* h1* h2*

* failed replicates (the x3 design above: 8 of 16 half-samples fail). One message per cause:
* the full-sample call's pre-flight warning, which names the level and the count; no line per
* replicate, with nodots too; nothing after the loop; at replay a one-line footnote
local w8 "warning: x3 = 4 has no focal group weight in 8 of 16 BRR replicates, which will fail. Collapse that level, or use Fay's BRR (fay())."
quietly log using `"`lg'"', text replace name(s1log)
svy brr, nodots: seldecomp y i.x3, by(group)
quietly log close s1log
assert e(N_misreps) == 8 & e(N_reps) == 8 & e(N_prefail) == 8
log_count `"`lg'"' "replicate failed"
assert r(n) == 0
log_has `"`lg'"' "focalgroup()"
assert r(found) == 0
log_has `"`lg'"' "raking support failure"
assert r(found) == 0
log_has `"`lg'"' "`w8'"
assert r(count) == 1
log_count `"`lg'"' "warning:"
assert r(n) == 1
log_has `"`lg'"' "too small"
assert r(found) == 0
quietly log using `"`lg'"', text replace name(s1log)
seldecomp
quietly log close s1log
log_has `"`lg'"' "(8 of 16 BRR replicates failed)"
assert r(count) == 1
log_has `"`lg'"' "warning"
assert r(found) == 0
* with dots, the failures are the prefix's x's; the warning is the same
quietly log using `"`lg'"', text replace name(s1log)
svy brr: seldecomp y i.x3, by(group)
quietly log close s1log
log_count `"`lg'"' "replicate failed"
assert r(n) == 0
log_has `"`lg'"' "`w8'"
assert r(count) == 1
log_count `"`lg'"' "warning:"
assert r(n) == 1
* svy jackknife (PSU 1 of stratum 2 deleted: the pre-flight names the PSU; svy's own note
* says that the std. errs. are missing, and seldecomp says nothing about them: under svy its
* messages say that a replicate fails, not that it is dropped) and
* jackknife: (no pre-flight: one warning after the loop, which names the level and does not
* refer to the thin-cell note, which is under the raking diagnostics);
* without Fay; at replay the footnote, with "replicates" plural, without a reference
svyset psuid [pw=wt], strata(stratum)
quietly log using `"`lg'"', text replace name(s1log)
svy jackknife, nodots: seldecomp y i.x3, by(group)
assert e(N_misreps) == 1 & e(N_prefail) == 1
jackknife, cluster(psuid) nodots: seldecomp y i.x3, by(group)
assert e(N_misreps) == 1 & missing(e(N_prefail))
seldecomp
quietly log close s1log
log_has `"`lg'"' "warning: x3 = 4 has no focal group weight when psuid 21 of stratum 2 is deleted: that jackknife replicate will fail. Collapse that level."
assert r(count) == 1
log_has `"`lg'"' "warning: 1 of 32 jackknife replicates failed and was dropped: in it, x3 = 4 had no focal group weight. Collapse that level."
assert r(count) == 1
log_has `"`lg'"' "(1 of 32 jackknife replicates failed and was dropped)"
assert r(count) == 1
log_has `"`lg'"' "the thin"
assert r(found) == 0
log_count `"`lg'"' "warning:"
assert r(n) == 2
log_has `"`lg'"' "Fay"
assert r(found) == 0
log_has `"`lg'"' "too small"
assert r(found) == 0
log_count `"`lg'"' "replicate failed"
assert r(n) == 0
* svy bootstrap with bsrweight(): the pre-flight warning, then the footnote at replay
svyset psuid [pw=wt], strata(stratum) bsrweight(brr1-brr16)
quietly log using `"`lg'"', text replace name(s1log)
svy bootstrap, nodots: seldecomp y i.x3, by(group)
assert e(N_misreps) == 8 & e(N_prefail) == 8
seldecomp
quietly log close s1log
log_has `"`lg'"' "warning: x3 = 4 has no focal group weight in 8 of 16 bootstrap replicates, which will fail. Collapse that level."
assert r(count) == 1
log_has `"`lg'"' "(8 of 16 bootstrap replicates failed)"
assert r(count) == 1
log_count `"`lg'"' "warning:"
assert r(n) == 1
* bootstrap: (no pre-flight): one warning after the loop, which names the level once, and
* the footnote at replay; neither refers to the thin-cell note
svyset, clear
quietly log using `"`lg'"', text replace name(s1log)
bootstrap, reps(40) cluster(psuid) strata(stratum) seed(11) nodots: seldecomp y i.x3, by(group)
local nm = e(N_misreps)
assert `nm' > 0 & `nm' < . & missing(e(N_prefail))
seldecomp
quietly log close s1log
log_has `"`lg'"' "warning: `nm' of 40 bootstrap replicates failed and were dropped: in them, x3 = 4 had no focal group weight. Collapse that level."
assert r(count) == 1
log_has `"`lg'"' "(`nm' of 40 bootstrap replicates failed and were dropped)"
assert r(count) == 1
log_has `"`lg'"' "had no focal group weight"
assert r(count) == 1
log_has `"`lg'"' "the thin"
assert r(found) == 0
log_count `"`lg'"' "warning:"
assert r(n) == 1
log_count `"`lg'"' "replicate failed"
assert r(n) == 0

* the pre-flight check: in the full-sample call of svy brr, svy
* jackknife, svy sdr and svy bootstrap, the replicates that will fail are found from the
* svyset replicate weights (svy jackknife without jkrweight(): from the PSUs) and the levels
* named, once, before the replicates; e(N_prefail) is their number, e(N_misreps) the number
* dropped. The x3 design: PSU 1 of stratum 2 holds all of the Focal Group's x3 = 4
capture program drop c7_order
program define c7_order
    // c7_order "file" "first" "second": the first line with "first" precedes the first with
    // "second"
    args file a b
    mata: _c7l = cat(st_local("file"))
    mata: st_local("pa", strofreal(min(select(1::rows(_c7l), strpos(_c7l, st_local("a")) :> 0))))
    mata: st_local("pb", strofreal(min(select(1::rows(_c7l), strpos(_c7l, st_local("b")) :> 0))))
    mata: mata drop _c7l
    assert `pa' < `pb' & `pb' < .
end
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr16)
quietly log using `"`lg'"', text replace name(s1log)
svy brr: seldecomp y i.x3, by(group)
seldecomp
quietly log close s1log
assert e(N_prefail) == 8 & e(N_misreps) == 8
log_has `"`lg'"' "`w8'"
assert r(count) == 1
c7_order `"`lg'"' "has no focal group weight in 8 of 16" "BRR replications (16)"
* the same weights as SDR, bootstrap and jackknife replicate weights
foreach m in sdr bsr jkr {
    local mn = cond("`m'" == "sdr", "sdr", cond("`m'" == "bsr", "bootstrap", "jackknife"))
    local mw = cond("`m'" == "sdr", "SDR", "`mn'")
    svyset psuid [pw=wt], strata(stratum) `m'weight(brr1-brr16)
    quietly log using `"`lg'"', text replace name(s1log)
    svy `mn', nodots: seldecomp y i.x3, by(group)
    quietly log close s1log
    assert e(N_prefail) == 8 & e(N_misreps) == 8
    log_has `"`lg'"' "warning: x3 = 4 has no focal group weight in 8 of 16 `mw' replicates, which will fail. Collapse that level."
    assert r(count) == 1
    log_has `"`lg'"' "Fay"
    assert r(found) == 0
}
* svy jackknife without jkrweight(): the delete-one-PSU replicates, the PSU named
svyset psuid [pw=wt], strata(stratum)
quietly log using `"`lg'"', text replace name(s1log)
svy jackknife: seldecomp y i.x3, by(group)
quietly log close s1log
assert e(N_prefail) == 1 & e(N_misreps) == 1
log_has `"`lg'"' "warning: x3 = 4 has no focal group weight when psuid 21 of stratum 2 is deleted: that jackknife replicate will fail. Collapse that level."
assert r(count) == 1
c7_order `"`lg'"' "has no focal group weight when" "Jackknife replications (32)"
* what the failed replicate does to the standard errors is svy's to say, in the
* notes under its table (here they are missing: its stratum is left with one replicate);
* seldecomp's messages under svy do not say "dropped", at estimation or at replay
log_has `"`lg'"' "dropped"
assert r(found) == 0
quietly log using `"`lg'"', text replace name(s1log)
seldecomp
quietly log close s1log
log_has `"`lg'"' "(1 of 32 jackknife replicates failed)"
assert r(count) == 1
log_has `"`lg'"' "dropped"
assert r(found) == 0
* with subpop(): the subpopulation is the sample (the PSUs of stratum 5 are left out)
svy jackknife, nodots subpop(if stratum != 5): seldecomp y i.x3, by(group)
assert e(N_prefail) == 1 & e(N_misreps) == 1
* Fay's BRR keeps every PSU in every replicate: nothing fails, nothing is said
forvalues r = 1/16 {
    gen double fay`r' = .5 * wt + .5 * brr`r'
}
svyset psuid [pw=wt], strata(stratum) brrweight(fay1-fay16) fay(.5)
quietly log using `"`lg'"', text replace name(s1log)
svy brr, nodots: seldecomp y i.x3, by(group)
quietly log close s1log
assert e(N_prefail) == 0 & e(N_misreps) == 0
log_has `"`lg'"' "warning"
assert r(found) == 0
drop fay*
* a level of the Reference Group held by one PSU (x4 = 5; Focal Group has it everywhere): in the 8 half-samples
* without that PSU its Focal Group observations get zero raked weight; nothing fails, and one
* warning names the level
gen byte x4 = x1
replace x4 = 5 if group == 0 & mod(_n, 5) == 0
replace x4 = 5 if group == 1 & stratum == 3 & psu == 2 & mod(_n, 2) == 0
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr16)
quietly log using `"`lg'"', text replace name(s1log)
svy brr, nodots: seldecomp y i.x4, by(group)
quietly log close s1log
assert e(N_prefail) == 0 & e(N_misreps) == 0
log_has `"`lg'"' "warning: 8 of 16 BRR replicates give focal group obs zero weight; collapse this level, or use Fay's BRR (fay()): x4 = 5: no reference group weight in 8 replicates"
assert r(count) == 1
log_count `"`lg'"' "warning:"
assert r(n) == 1
log_has `"`lg'"' "will fail"
assert r(found) == 0
* a level (x5 = 5, tied to x6 = 9 in PSU 1 of
* stratum 4) that has no focal group weight in 4 half-samples and no reference group weight in 4 others:
* one warning, the level listed once with both counts
gen byte x6 = x2
replace x6 = 9 if stratum == 4 & psu == 1
gen byte x5 = x1
replace x5 = 5 if group == 0 & x6 == 9 & mod(_n, 2) == 0
gen long o = sum(group == 1 & _n > 150)
replace x5 = 5 if group == 1 & _n > 150 & o <= 3
quietly log using `"`lg'"', text replace name(s1log)
svy brr, nodots: seldecomp y i.x5 i.x6, by(group)
quietly log close s1log
assert e(N_prefail) > 0 & e(N_misreps) == e(N_prefail)
log_count `"`lg'"' "warning:"
assert r(n) == 1
log_has `"`lg'"' "x5 = 5: no focal group weight in `=e(N_prefail)' replicates, no reference group weight in"
assert r(count) == 1
drop x5 x6 o
qui seldecomp y i.x4 [pw=brr1], by(group)
assert e(N_focal_zero) > 0 & missing(e(N_prefail))
* not checked: svy brr with hadamard() (the matrix is invisible to the command), a plain run,
* the bootstrap: prefix
svyset psuid [pw=wt], strata(stratum)
svy brr, nodots hadamard(HAD): seldecomp y i.x3, by(group)
assert e(N_misreps) == 8 & missing(e(N_prefail))
seldecomp y i.x3 [pw=wt], by(group)
assert missing(e(N_prefail))
bootstrap, reps(3) seed(1) nodots: seldecomp y i.x1, by(group)
assert missing(e(N_prefail))
* hadamard() typed while brrweight() variables are svyset: svy uses the matrix
* and not the weights, and the check must not read them. HAD2 puts stratum 2 on the all-ones
* column, so PSU 21 is in every half-sample and nothing fails, while 8 fail under the svyset
* weights: no warning, e(N_prefail) not posted (the route is recognised from the sort by the
* svyset strata and PSU and the prefix's working variables)
mata: H2 = st_matrix("HAD"); H2[., (1, 2)] = H2[., (2, 1)]; st_matrix("HAD2", H2)
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr16)
quietly log using `"`lg'"', text replace name(s1log)
svy brr, nodots hadamard(HAD2): seldecomp y i.x3, by(group)
quietly log close s1log
assert missing(e(N_prefail)) & e(N_misreps) == 0 & e(N_reps) == 16
log_has `"`lg'"' "note: using hadamard() option instead of replicate weights"
assert r(found) == 1
log_has `"`lg'"' "will fail"
assert r(found) == 0
svy brr, nodots hadamard(HAD): seldecomp y i.x3, by(group)
assert missing(e(N_prefail)) & e(N_misreps) == 8
* the same with the data already sorted as the prefix sorts it; and the weights route with
* that sort order and no hadamard() is still checked (the working variables tell them apart)
sort stratum psuid, stable
svy brr, nodots hadamard(HAD2): seldecomp y i.x3, by(group)
assert missing(e(N_prefail)) & e(N_misreps) == 0
svy brr, nodots: seldecomp y i.x3, by(group)
assert e(N_prefail) == 8 & e(N_misreps) == 8
* a level label with SMCL braces, quotes and a dollar sign: the pre-flight
* warning and the warning after the loop (hadamard(): no pre-flight) show it as text, and the
* output after them renders
label define x3l 1 "one" 2 "two" 3 "three"
mata: st_vlmodify("x3l", 4, "the " + char(34) + "four" + char(34) + " " + char(36) + "x {p}")
label values x3 x3l
quietly log using `"`lg'"', text replace name(s1log)
svy brr, nodots: seldecomp y i.x3, by(group)
assert e(N_prefail) == 8 & e(N_misreps) == 8
svy brr, nodots hadamard(HAD): seldecomp y i.x3, by(group)
assert missing(e(N_prefail)) & e(N_misreps) == 8
quietly log close s1log
log_has `"`lg'"' `"x {p}) has no focal group weight in 8 of 16 BRR replicates, which will fail."'
assert r(count) == 1
log_has `"`lg'"' `"warning: 8 of 16 BRR replicates failed: in them, x3 = 4 (the "four" "'
assert r(count) == 1
log_has `"`lg'"' `"x {p}) had no focal group weight. Collapse that level, or use Fay's BRR (fay())."'
assert r(count) == 1
log_has `"`lg'"' "replicate failed"
assert r(found) == 0
log_has `"`lg'"' "{p_end}"
assert r(found) == 0
log_has `"`lg'"' "{col "
assert r(found) == 0
label values x3
label drop x3l
drop x4

* after svy, the design effects are missing (no linearized variance), not 0; citable
* with nothing to show says so
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr16)
svy brr, nodots: seldecomp y i.x1 i.x2, by(group)
estat effects
matrix DE = r(deff)
matrix DT = r(deft)
mata: assert(allof(st_matrix("DE"), .) & allof(st_matrix("DT"), .))
assert colsof(DE) == colsof(e(b))
quietly log using `"`lg'"', text replace name(s1log)
seldecomp, citable
quietly log close s1log
log_has `"`lg'"' "(citable: no bayesboot() replicates; the CIs above are BRR)"
assert r(found) == 1
* svy with subpop(): the subpopulation count below the header N
svy brr, nodots subpop(if x2 == 1): seldecomp y i.x1, by(group)
assert e(N_sub) == e(N_focal) + e(N_ref) & e(N) == _N
estat effects
matrix DE = r(deffsub)
mata: assert(allof(st_matrix("DE"), .))
quietly log using `"`lg'"', text replace name(s1log)
seldecomp
quietly log close s1log
log_has `"`lg'"' "Number of obs = `=strtrim(string(e(N), "%12.0fc"))' Subpop. no. obs = `=strtrim(string(e(N_sub), "%12.0fc"))' Focal Group (raked):"
assert r(found) == 1
* without subpop() no such line
svy brr, nodots: seldecomp y i.x1, by(group)
quietly log using `"`lg'"', text replace name(s1log)
seldecomp
quietly log close s1log
log_has `"`lg'"' "Subpop."
assert r(found) == 0

* trimming bounds that the full sample can satisfy and some replicates cannot (trimhirel(1.98);
* below 1.93 the full sample misses its targets too): those replicates fail, with the one
* warning of failed replicates, which does not blame the iterate() limit; force keeps them
quietly log using `"`lg'"', text replace name(s1log)
svy brr, nodots: seldecomp y i.x1 i.x2, by(group) trimhirel(1.98)
quietly log close s1log
assert e(converged) == 1 & e(targets_met) == 1
local nmis = e(N_misreps)
assert `nmis' > 0 & `nmis' < 16 & e(N_reps) == 16 - `nmis'
log_has `"`lg'"' "warning: `nmis' of 16 BRR replicates failed: in them, raking did not reach its targets."
assert r(count) == 1
log_has `"`lg'"' "iterate()"
assert r(found) == 0
rcof "noi svy brr, nodots: seldecomp y i.x1 i.x2, by(group) trimhirel(1.5)" == 430
svy brr, nodots: seldecomp y i.x1 i.x2, by(group) trimhirel(1.98) force
assert e(N_misreps) == 0 & e(N_reps) == 16 & e(converged) == 1 & e(targets_met) == 1

* svy prints "Convergence not achieved." when a miss is kept by force, also a miss that
* the trimming bounds cause (without force that miss is rc 430, from the full-sample call).
* Margins that cannot all be met (Focal Group has a == 1 <=> b == 0), 16 strata of 2 PSUs,
* jackknife by PSU
preserve
clear
set obs 400
gen byte grp = _n > 200
gen byte a = cond(grp, mod(_n, 10) < 3, mod(_n, 2))
gen byte b = cond(grp, mod(floor(_n / 3), 10) < 6, 1 - mod(_n, 2))
gen double yy = a + b + mod(_n, 7) / 7
gen int st = ceil(_n / 25)
gen byte pu = 1 + (mod(_n, 4) < 2)
gen byte one = 1
svyset pu [pw=one], strata(st)
quietly log using `"`lg'"', text replace name(s1log)
capture noisily svy jackknife, nodots: seldecomp yy i.a i.b, by(grp) trimhirel(1.1) trimlorel(.9)
local rc = _rc
quietly log close s1log
assert `rc' == 430
log_has `"`lg'"' "the trimming bounds are binding: relax them, or specify force to keep the results"
assert r(count) == 1
quietly log using `"`lg'"', text replace name(s1log)
svy jackknife, nodots: seldecomp yy i.a i.b, by(grp) trimhirel(1.1) trimlorel(.9) force
quietly log close s1log
assert e(converged) == 0 & e(targets_met) == 0 & e(badcontrols) > 0 & e(N_misreps) == 0
log_has `"`lg'"' "Convergence not achieved"
assert r(found) == 1
* the trimming note, which svy does not show: once, as a warning of the full-sample call
log_has `"`lg'"' "warning: trimming left"
assert r(count) == 1
log_has `"`lg'"' "results kept because force was specified"
assert r(count) == 1
log_count `"`lg'"' "warning:"
assert r(n) == 1
quietly log using `"`lg'"', text replace name(s1log)
svy jackknife, nodots: seldecomp yy i.a i.b, by(grp) trimlorel(.81) force
quietly log close s1log
assert e(converged) == 0 & e(targets_met) == 0 & e(badcontrols) > 0
log_has `"`lg'"' "Convergence not achieved"
assert r(found) == 1
* the note of the margins kept by force, which svy does not show: once, as a warning of the
* full-sample call
log_has `"`lg'"' "results kept because force was specified"
assert r(count) == 1
log_count `"`lg'"' "warning:"
assert r(n) == 1
* e(converged) is the minimum over the stages: under sequentially
* fixothermargins a stage 1 miss kept by force makes it 0 although the last stage converged,
* and svy warns. b2 = 1 - a: stage 1 (a raked to the Reference Group's .3, b2 held at the
* Focal Group's .5 by fixothermargins) cannot be met, stage 2 can. Margins that trimming
* leaves unmet at every stage, kept by force, make it 0 too, and svy warns (block A holds
* the margins that cannot all be met; d is held at stage 1)
gen byte b2 = 1 - a
gen byte d = mod(floor(_n / 5), 2)
quietly log using `"`lg'"', text replace name(s1log)
svy jackknife, nodots: seldecomp yy (A: a) (B: b2), by(grp) sequentially fixothermargins force
quietly log close s1log
assert e(converged) == 0 & e(targets_met) == 0 & e(badcontrols) == 0
* (every replicate posts e(converged) = 0 too, and svy keeps them all)
assert e(N_misreps) == 0 & e(N_reps) == 32
assert el(e(bystage), 1, 6) == 0 & el(e(bystage), 2, 6) == 1
log_has `"`lg'"' "Convergence not achieved"
assert r(found) == 1
quietly log using `"`lg'"', text replace name(s1log)
svy jackknife, nodots: seldecomp yy (A: a b) (D: d), by(grp) sequentially fixothermargins trimhirel(1.1) trimlorel(.9) force
quietly log close s1log
assert e(converged) == 0 & e(targets_met) == 0
assert el(e(bystage), 1, 6) == 0 & el(e(bystage), 2, 6) == 0
assert el(e(bystage), 1, 7) > 0 & el(e(bystage), 2, 7) > 0
log_has `"`lg'"' "Convergence not achieved"
assert r(found) == 1
drop b2 d
quietly log using `"`lg'"', text replace name(s1log)
svy jackknife, nodots: seldecomp yy i.a i.b, by(grp) trimlorel(.81) force
quietly log close s1log
* a replicate that keeps its missed margins by force prints nothing, and
* so does one that misses them without force (a replicate call made by hand: c(prefix) as the
* prefix sets it); the full-sample call's warning is the one line
log_has `"`lg'"' "because force was specified"
assert r(count) == 1
capture program drop s1_fake
program define s1_fake, eclass
    gettoken words 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    set prefix `words'
    `0'
end
quietly log using `"`lg'"', text replace name(s1log)
capture noisily s1_fake svy jackknife _loop _loop_rw : seldecomp yy i.a i.b [pw=one], by(grp)
local rc = _rc
quietly log close s1log
assert `rc' == 430
log_has `"`lg'"' "replicate failed"
assert r(found) == 0
log_has `"`lg'"' "raking did not reach its targets"
assert r(found) == 0
log_has `"`lg'"' "force"
assert r(found) == 0
restore
svyset, clear

* the pre-flight check on nhanes2brr (the Black cell of West x not in an SMSA has 11
* observations; needs the network, skipped without it)
preserve
capture webuse nhanes2brr, clear
if _rc == 0 {
    quietly log using `"`lg'"', text replace name(s1log)
    svy brr, nodots: seldecomp highbp (ses: i.hsizgp) (geo: i.region#i.smsa) (age: i.agegrp) ///
        (sex: i.sex) if race < 3, by(race) focalgroup(Black) sequentially fixothermargins
    quietly log close s1log
    di as txt "SDW_C7_NHANES N_prefail=" e(N_prefail) " N_misreps=" e(N_misreps) " N_reps=" e(N_reps)
    assert e(N_prefail) == e(N_misreps) & e(N_prefail) > 0
    log_has `"`lg'"' "region#smsa = "
    assert r(found) == 1
}
else di as txt "SDW_C7_NHANES skipped (nhanes2brr not available)"
restore

* a replicate that misses a target, which the pre-flight check does not
* foresee (it looks for support only): the margins of test_bootstrap.do 9b (Focal Group lacks the
* cell a = 1, b = 1; Reference Group's share of it is a little below its share of (0, 0)); cluster 11
* holds 10 of the Reference Group's 30 (0, 0) observations, so the jackknife replicate without it cannot
* meet the margins and fails. The check runs and finds nothing (e(N_prefail) 0), so the one
* warning after the loop names the cause
preserve
clear
set obs 400
gen byte group = _n > 200
gen int k = cond(group == 0, _n, _n - 200)
gen byte a = cond(group == 0, k > 130, cond(k <= 30, 0, cond(k <= 54, 1, !mod(k, 2))))
gen byte b = cond(group == 0, inrange(k, 61, 130), cond(k <= 30, 0, cond(k <= 54, 1, mod(k, 2))))
gen double y = a + 2 * b + mod(_n * 7, 11) / 10
gen int cl = cond(group == 1 & k <= 10, 11, 1 + mod(k, 10))
svyset cl
quietly log using `"`lg'"', text replace name(s1log)
svy jackknife, nodots: seldecomp y a b, by(group)
quietly log close s1log
assert e(N_prefail) == 0 & e(N_misreps) == 1
log_has `"`lg'"' "warning: 1 of 11 jackknife replicates failed: in it, raking did not reach its targets."
assert r(count) == 1
* (no remedy: nothing certain to suggest for a missed target)
log_has `"`lg'"' "sparse levels"
assert r(found) == 0
log_count `"`lg'"' "warning:"
assert r(n) == 1
log_has `"`lg'"' "further"
assert r(found) == 0
restore

di as res _n "PASSED: test_svy.do"

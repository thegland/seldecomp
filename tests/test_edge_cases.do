* test_edge_cases.do -- focalgroup() with two values within a relative 1e-7 (large codes:
* yyyymmdd dates, tract ids); the svy pre-flight warning's number of replicates under
* subpop(); labels with SMCL in the error messages; the labels of non-integer levels of a
* float variable; factor-variable operators with parentheses; the stage table's percent
* columns; bayesboot(saving(, replace)) checked before the replications; a Break in the
* reverse order; the remedy's echo of a focalgroup() text with $ or a backquote; a Break
* inside the Mata raker, the BCa Break with bayesboot(generate()) included (section 10);
* estat bootstrap after a replay at a level such as 99.9, a constant raking variable with
* insupport, a number typed with a plus sign, the one observation with zero weight
* (section 11); a constant raking variable with a support failure, a number of quantile
* groups typed with a plus sign (section 12).
* The staged runs use sequentially fixothermargins.
*
* Run from the tests directory:   . do test_edge_cases.do

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

* ---- helpers --------------------------------------------------------------------
* rv_log "command": runs the command (capture noisily) with a text log; r(rc) its return
* code; the Mata string scalar _rv_txt holds its output, lines joined by one blank, runs of
* blanks squeezed (a wrapped paragraph matches); _rv_lines the lines themselves
capture program drop rv_log
program define rv_log, rclass
    args cmd
    tempfile lg
    quietly log using `"`lg'"', text replace name(rvlog) nomsg
    capture noisily `cmd'
    local rc = _rc
    quietly log close rvlog
    mata: _rv_lines = cat(st_local("lg")); _rv_txt = stritrim(invtokens(strtrim(_rv_lines)'))
    return scalar rc = `rc'
end
* rv_has "text" [0]: the last rv_log output contains the text (0: does not)
capture program drop rv_has
program define rv_has
    args text want
    if "`want'" == "" local want 1
    mata: st_local("f", strofreal(strpos(_rv_txt, st_local("text")) > 0))
    if `f' != `want' {
        di as err `"output "' cond(`want', "lacks", "contains") `": `text'"'
        exit 9
    }
end
* rv_width: r(w), the widest line of the last rv_log output that is not a command echo
capture program drop rv_width
program define rv_width, rclass
    mata: _rv_l = select(_rv_lines, substr(_rv_lines, 1, 2) :!= ". " :& substr(_rv_lines, 1, 2) :!= "> ")
    mata: st_local("w", strofreal(max(udstrlen(_rv_l))))
    mata: mata drop _rv_l
    return scalar w = `w'
end


* ================================================================================
* 1. focalgroup(): a value typed exactly is that value, also when the two values are within
*    a relative 1e-7 of each other; a number between them is no value
* ================================================================================
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
qui seldecomp y x1 x2 [pw=wt], by(group) focalgroup(1)
matrix B1 = e(b)
qui seldecomp y x1 x2 [pw=wt], by(group) focalgroup(0)
matrix B0 = e(b)
assert mreldif(B1, B0) > .1
* yyyymmdd dates (long) and 11-digit tract codes (double); group == 1 is the higher value
gen long day = cond(group == 1, 20010102, 20010101)
gen double tract = cond(group == 1, 36061000200, 36061000100)
local hi_day 20010102
local lo_day 20010101
local mid_day 20010103
local hi_tract 36061000200
local lo_tract 36061000100
local mid_tract 36061000150
foreach v in day tract {
    assert reldif(`hi_`v'', `lo_`v'') < 1e-7
    qui seldecomp y x1 x2 [pw=wt], by(`v') focalgroup(`hi_`v'')
    assert real("`e(focal_group)'") == `hi_`v'' & mreldif(e(b), B1) == 0
    qui seldecomp y x1 x2 [pw=wt], by(`v') focalgroup(`lo_`v'')
    assert real("`e(focal_group)'") == `lo_`v'' & mreldif(e(b), B0) == 0
    // a number within 1e-7 of both values is neither
    rcof "noi seldecomp y x1 x2 [pw=wt], by(`v') focalgroup(`mid_`v'')" == 459
}
* a float by() typed as a decimal still matches (the tolerance's purpose), as does its
* full double expansion
gen float gf = cond(group == 1, .2, .1)
qui seldecomp y x1 x2 [pw=wt], by(gf) focalgroup(0.2)
assert mreldif(e(b), B1) == 0
qui seldecomp y x1 x2 [pw=wt], by(gf) focalgroup(.1)
assert mreldif(e(b), B0) == 0
qui seldecomp y x1 x2 [pw=wt], by(gf) focalgroup(.2000000029802322)
assert mreldif(e(b), B1) == 0
drop day tract gf
di as txt "section 1: focalgroup() with close values passed"


* ================================================================================
* 2. the svy pre-flight warning: the replicates svy runs. With replicate weights,
*    svy jackknife skips a variable equal to the full-sample weight on the subpopulation;
*    with the PSU jackknife the number is not known in the call (under subpop(), svy has
*    the PSUs of the whole strata), and the warning gives the count alone
* ================================================================================
* 8 strata x 3 PSUs x 10 observations; two levels of x3 that Focal Group (group == 0) has in
* one PSU each; subpop(sub) leaves out stratum 5 and PSUs 2 and 3 of stratum 6
clear
set seed 4711
set obs 240
gen int stratum = ceil(_n / 30)
gen byte psu = 1 + mod(ceil(_n / 10) - 1, 3)
gen double wt = 20 + 60 * runiform()
gen byte group = runiform() < .5
gen byte x1 = 1 + (runiform() < .4) + (runiform() < .3 + .3 * group)
gen double y = 2 + x1 + .8 * group + rnormal()
gen byte x3 = x1
replace x3 = 4 if group == 1 & mod(_n, 7) == 0
replace x3 = 5 if group == 1 & mod(_n, 7) == 3
replace x3 = 4 if group == 0 & stratum == 2 & psu == 1 & mod(_n, 3) == 0
replace x3 = 5 if group == 0 & stratum == 3 & psu == 2 & mod(_n, 3) == 0
gen byte sub = stratum != 5 & !(stratum == 6 & psu > 1)
local r 0
forvalues h = 1/8 {
    forvalues i = 1/3 {
        local ++r
        gen double jk`r' = cond(stratum == `h', cond(psu == `i', 0, wt * 3 / 2), wt)
    }
}
* the PSU jackknife, with and without subpop(): the count alone
svyset psu [pw=wt], strata(stratum) vce(jackknife)
foreach sp in "subpop(sub)" "" {
    rv_log "svy jackknife, nodots `sp': seldecomp y i.x3, by(group)"
    assert r(rc) == 0
    assert e(N_prefail) == 2 & e(N_misreps) == 2
    rv_has "warning: 2 jackknife replicates will fail; collapse these levels:"
    rv_has "of 19 jackknife" 0
    rv_has "of 24 jackknife" 0
}
* replicate weights: 21 of the 24 variables differ from wt on the subpopulation, which is
* what svy runs; without subpop(), 24
svyset psu [pw=wt], strata(stratum) jkrweight(jk1-jk24, multiplier(`=2/3')) vce(jackknife)
rv_log "svy jackknife, nodots subpop(sub): seldecomp y i.x3, by(group)"
assert r(rc) == 0
assert e(N_prefail) == 2 & e(N_misreps) == 2 & e(N_reps) + e(N_misreps) == 21
rv_has "warning: 2 of 21 jackknife replicates will fail; collapse these levels:"
rv_log "svy jackknife, nodots: seldecomp y i.x3, by(group)"
assert r(rc) == 0 & e(N_reps) + e(N_misreps) == 24
rv_has "warning: 2 of 24 jackknife replicates will fail; collapse these levels:"
* SDR runs every variable (the rule is svy jackknife's)
svyset psu [pw=wt], strata(stratum) sdrweight(jk1-jk24) vce(sdr)
rv_log "svy sdr, nodots subpop(sub): seldecomp y i.x3, by(group)"
assert r(rc) == 0 & e(N_reps) + e(N_misreps) == 24
rv_has "warning: 2 of 24 SDR replicates will fail; collapse these levels:"
* chain: x6 = 6 in period 1 (Focal Group of pair 1) only in PSU 1 of stratum 2
gen byte period = 1 + mod(_n, 3)
gen byte x6 = x1
replace x6 = 6 if period > 1 & mod(_n, 5) == 0
replace x6 = 6 if period == 1 & stratum == 2 & psu == 1 & mod(_n, 2) == 0
svyset psu [pw=wt], strata(stratum) vce(jackknife)
rv_log "svy jackknife, nodots subpop(sub): seldecomp y i.x6, by(period) chain"
assert r(rc) == 0 & e(N_prefail) == 1 & e(N_misreps) == 1
rv_has "warning: p1: x6 = 6 has no focal group weight in 1 jackknife replicate, which will fail. Collapse that level."
svyset psu [pw=wt], strata(stratum) jkrweight(jk1-jk24, multiplier(`=2/3')) vce(jackknife)
rv_log "svy jackknife, nodots subpop(sub): seldecomp y i.x6, by(period) chain"
assert r(rc) == 0 & e(N_prefail) == 1 & e(N_misreps) == 1 & e(N_reps) + e(N_misreps) == 21
rv_has "warning: p1: x6 = 6 has no focal group weight in 1 of 21 jackknife replicates, which will fail. Collapse that level."
svyset, clear
di as txt "section 2: the pre-flight warning's number of replicates passed"


* ================================================================================
* 3. labels with SMCL in the error messages are shown as typed, as in the header,
*    the notes and covariatedetail
* ================================================================================
capture program drop rv_lbldata
program define rv_lbldata
    clear
    set seed 4242
    set obs 480
    gen g = _n > 240
    gen x = mod(_n, 8) + 1
    gen z = mod(floor(_n / 8), 2)
    gen y = x / 3 + z + rnormal() + g
    mata: st_vlmodify("gl", (0 \ 1), ("single {it:x}" \ "married {c -(}lit{c )-}"))
    mata: st_vlmodify("xl", 9, "nine {p_end}{dup 70:=}")
    label values g gl
    label values x xl
end
* the support failure: x = 9 in the Reference Group only
rv_lbldata
replace x = 9 in 471/480
rv_log "seldecomp y i.x z, by(g)"
assert r(rc) == 459
rv_has "these levels of the reference group (g = 1 (married {c -(}lit{c )-})) have no focal group observations (g = 0 (single {it:x})):"
rv_has "x = 9 (nine {p_end}{dup 70:=})"
rv_has "===========" 0
rv_width
assert r(w) <= 78
* three values of by(), and a focalgroup() that is none of them
rv_lbldata
replace g = 2 in 1/5
mata: st_vlmodify("gl", 2, "other {it:y}")
rv_log "seldecomp y i.x z, by(g)"
assert r(rc) == 459
rv_has "by(g) takes 3 values (0 (single {it:x}), 1 (married {c -(}lit{c )-}), and 2 (other {it:y})); it must take two"
rv_lbldata
rv_log "seldecomp y i.x z, by(g) focalgroup(nosuch)"
assert r(rc) == 459
rv_has "its values are 0 (single {it:x}) and 1 (married {c -(}lit{c )-})"
* insupport: once x = 9 (Focal Group only) is excluded, z = 1 occurs in the Reference Group alone
rv_lbldata
replace x = 9 if g == 0 & _n <= 20
replace z = (x == 9) if g == 0
mata: st_vlmodify("zl", 1, "one {it:z}")
label values z zl
rv_log "seldecomp y i.x i.z, by(g) insupport"
assert r(rc) == 459
rv_has "insupport: after the exclusion, these levels of the reference group (g = 1 (married {c -(}lit{c )-})) have no focal group observations:"
rv_has "z = 1 (one {it:z})"
* the margins that cannot be met (430): the worst level
clear
set seed 8
set obs 400
gen g = _n > 200
gen a = cond(g == 0, _n <= 100, runiform() < .3)
gen b = cond(g == 0, 1 - a, runiform() < .3)
gen y = a + b + rnormal() + g
mata: st_vlmodify("al", (0 \ 1), ("a0 {it:zero}" \ "a1 {p_end}{dup 70:-}"))
label values a al
rv_log "seldecomp y a b, by(g)"
assert r(rc) == 430
rv_has "worst: a = 1 (a1 {p_end}{dup 70:-})"
rv_has "-----------" 0
rv_width
assert r(w) <= 78
di as txt "section 3: labels in the error messages passed"


* ================================================================================
* 4. non-integer levels of a float variable are labelled as Stata shows them:
*    the shortest text of 7 to 9 significant digits that is that float
* ================================================================================
clear
set seed 5
set obs 600
gen byte g = _n > 300
gen byte x = 1 + mod(_n, 3)
gen byte z = runiform() < .5
gen double y = x + z + g + rnormal()
gen float xf = x / 10
gen double xd = x / 10
qui seldecomp y i.xf z, by(g) covariatedetail
assert `"`e(levels1)'"' == ".1 | .2 | .3"
qui seldecomp y i.xd z, by(g)
assert `"`e(levels1)'"' == ".1 | .2 | .3"
* a level found in the Focal Group only: its note
replace xf = .4 in 1/5
qui seldecomp y i.xf z, by(g)
local nn = e(N_notes)
local found 0
forvalues m = 1/`nn' {
    if strpos(`"`e(note`m')'"', "xf = .4 (5 obs)") local found 1
}
assert `found'
* two floats that 7 digits do not tell apart keep their own texts
replace xf = cond(x == 1, 1234567.5, cond(x == 2, 1234568, 2.25))
qui seldecomp y i.xf z, by(g)
assert `"`e(levels1)'"' == "2.25 | 1234567.5 | 1234568"
di as txt "section 4: float levels passed"


* ================================================================================
* 5. factor-variable operators with parentheses are refused as operators, not taken for
*    an option typed before the comma
* ================================================================================
clear
set seed 3
set obs 200
gen byte g = _n > 100
gen byte x = 1 + mod(_n, 3)
gen byte z = mod(_n, 2)
gen double y = x + z + g + rnormal()
foreach spec in "ib(first).x z" "z ib(last).x" "ib(freq).x z" "ib(#2).x z" "i(1 2).x z" ///
    "z ib(first).x#z" "(d: ib(first).x z)" {
    rv_log "seldecomp y `spec', by(g)"
    assert r(rc) == 198
    mata: st_local("op", ustrregexm(st_local("spec"), "[a-z]+[(][^()]*[)][.]x") ? ustrregexs(0) : "")
    rv_has "`op': only the i. and c. operators are allowed in the raking variable list"
    rv_has "must follow the comma" 0
}
* an option typed before the comma still gets its message
rv_log "seldecomp y x z by(g)"
assert r(rc) == 198
rv_has "option by() must follow the comma"
di as txt "section 5: factor-variable operators passed"


* ================================================================================
* 6. the stage table's percent columns stay within their 6 columns: shares of
*    1e5% and more (a difference near 0 against sizeable block contributions)
* ================================================================================
clear
set seed 7
set obs 4000
gen byte grp = _n > 2000
gen byte x = runiform() < cond(grp, .7, .3)
gen byte z = runiform() < cond(grp, .4, .6)
gen double y = 5 + x - z + rnormal()
qui summarize y if grp == 0
local m0 = r(mean)
qui summarize y if grp == 1
local m1 = r(mean)
replace y = y - (`m1' - `m0') - 3e-4 if grp == 1
local ls0 = c(linesize)
set linesize 120
foreach spec in "(a: x) (b: z)" "(a: x z)" {
    local o = cond("`spec'" == "(a: x z)", "detail", "")
    rv_log "seldecomp y `spec', by(grp) sequentially fixothermargins `o'"
    assert r(rc) == 0
    assert abs(el(e(bystage), 1, 4)) >= 1e5
    rv_width
    assert r(w) <= 78
    rv_has "e+0" 0
    rv_has "-1.4e5"
}
set linesize `ls0'
di as txt "section 6: stage table widths passed"


* ================================================================================
* 7. bayesboot(saving(, replace)): a file that cannot be written is refused before any
*    replication; an existing file is not touched by the check
* ================================================================================
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
qui seldecomp y x1 x2 [pw=wt], by(group)
local cmd0 `"`e(cmdline)'"'
local k0 = c(k)
rv_log "seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(5) rseed(1) generate(bw) saving(no_such_dir_rvt/sub/bsfile, replace))"
assert r(rc) == 603
rv_has "bayesboot(saving()): file no_such_dir_rvt/sub/bsfile.dta cannot be written"
rv_has "Bayesian bootstrap" 0
assert c(k) == `k0' & `"`e(cmdline)'"' == `"`cmd0'"'
capture describe bw_*
assert _rc != 0
* an existing file stays as it was when the command fails after the check
capture erase rvt_keep.dta
capture erase rvt_new.dta
preserve
clear
set obs 3
gen k = _n
qui save rvt_keep
restore
gen byte xbad = cond(group == 1, 1 + mod(_n, 3), 1 + mod(_n, 2))
rcof "noi seldecomp y xbad [pw=wt], by(group) bayesboot(reps(5) rseed(1) saving(rvt_keep, replace))" == 459
rcof "noi seldecomp y xbad [pw=wt], by(group) bayesboot(reps(5) rseed(1) saving(rvt_new, replace))" == 459
preserve
use rvt_keep, clear
assert _N == 3 & k[3] == 3
restore
capture confirm file rvt_new.dta
assert _rc != 0
* and is replaced by a run that works
qui seldecomp y x1 x2 [pw=wt], by(group) bayesboot(reps(5) rseed(1) nodots saving(rvt_keep, replace))
preserve
use rvt_keep, clear
assert _N == 5
restore
erase rvt_keep.dta
drop xbad
di as txt "section 7: bayesboot(saving(, replace)) passed"


qui seldecomp y x1 x2 [pw=wt], by(group)
local cmd0 `"`e(cmdline)'"'

* ================================================================================
* 8. a Break in the reverse order stops without a message, as in the forward order
*    and the detail chain. The Break is simulated by the Mata raker's hook: with the
*    external seldecomp_brk = n the (n+1)th raking stops with return code 1
* ================================================================================
foreach n in 1 2 3 {
    mata: seldecomp_brk = `n'
    rv_log "seldecomp y (a: x1) (b: x2 x4) [pw=wt], by(group) sequentially fixothermargins reverse"
    assert r(rc) == 1
    rv_has "raking failed" 0
    rv_has "reverse order:" 0
    mata: st_local("left", strofreal(findexternal("seldecomp_brk") != NULL))
    assert `left' == 0
}
assert `"`e(cmdline)'"' == `"`cmd0'"'
di as txt "section 8: Break in the reverse order passed"


* ================================================================================
* 9. the support failure's remedy shows the focalgroup() text as typed: a $ or a
*    backquote in a string value of by()
* ================================================================================
clear
set obs 60
gen double y = mod(_n, 7) + .5 * (_n > 30)
gen str10 sg = cond(_n <= 30, "A" + char(36) + "S_OS", "B")
gen byte x = cond(_n <= 30, 1 + mod(_n, 2), 1 + mod(_n, 3))
mata: st_local("fg", "A" + char(36) + "S_OS")
tempfile lg
quietly log using `"`lg'"', text replace name(rvlog) nomsg
capture noisily seldecomp y x, by(sg) focalgroup(`macval(fg)')
local rc = _rc
quietly log close rvlog
assert `rc' == 459
mata: _rv_txt = stritrim(invtokens(strtrim(cat(st_local("lg")))'))
mata: st_local("ok", strofreal(strpos(_rv_txt, "in place of focalgroup(A" + char(36) + "S_OS)") > 0))
assert `ok'
replace sg = "A" + char(96) + "q" + char(39) if _n <= 30
mata: st_local("fg", "A" + char(96) + "q" + char(39))
quietly log using `"`lg'"', text replace name(rvlog) nomsg
capture noisily seldecomp y x, by(sg) focalgroup(`macval(fg)')
local rc = _rc
quietly log close rvlog
assert `rc' == 459
mata: _rv_txt = stritrim(invtokens(strtrim(cat(st_local("lg")))'))
mata: st_local("ok", strofreal(strpos(_rv_txt, "in place of focalgroup(A" + char(96) + "q" + char(39) + ")") > 0))
assert `ok'
di as txt "section 9: the remedy's focalgroup() text passed"


* ================================================================================
* 10. a Break inside the Mata raker leaves nothing behind. The Break is simulated
*     by the hook of seldecomp_rake: with the Mata external seldecomp_brk = n the (n+1)th
*     raking stops after its first cycle with return code 1, exactly where a Break would stop
*     it, and the hook removes the external. The data in a random order, a previous
*     estimation result, and the checks of tests/README.md ("Manual check: Break"): the
*     staged run, bayesboot(generate() saving()), and the BCa jackknife with
*     bayesboot(generate())
* ================================================================================
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
set seed 11
gen double u = runiform()
sort u
drop u
gen long ord = _n
qui seldecomp y x1 x2 [pw=wt], by(group)
local cmd0 `"`e(cmdline)'"'
qui datasignature
local sig `r(datasignature)'
local k0 = c(k)
local n0 = c(N)
capture erase rvtbreak.dta
local calls0 0 1 5 11
local cmd1 seldecomp y (a: x1) (b: x2 x4) [pw=wt], by(group) sequentially fixothermargins generate(rk)
local cmd2 seldecomp y (a: x1) (b: x2 x4) [pw=wt], by(group) sequentially fixothermargins generate(rk)
local cmd3 seldecomp y x1 x2 x4 [pw=wt], by(group) generate(rk) bayesboot(reps(20) rseed(1) nodots generate(bw) saving(rvtbreak, replace))
local cmd4 seldecomp y x1 x2 [pw=wt] in 1/200, by(group) bayesboot(reps(10) rseed(1) bca nodots generate(bw))
forvalues j = 1/4 {
    local n : word `j' of `calls0'
    mata: seldecomp_brk = `n'
    capture noisily `cmd`j''
    assert _rc == 1
    mata: st_local("left", strofreal(findexternal("seldecomp_brk") != NULL))
    assert `left' == 0
    qui datasignature
    assert "`r(datasignature)'" == "`sig'"
    assert c(k) == `k0' & c(N) == `n0'
    assert ord == _n
    assert "`: sortedby'" == ""
    capture confirm variable rk
    assert _rc != 0
    capture describe bw_*
    assert _rc != 0
    capture confirm file rvtbreak.dta
    assert _rc != 0
    mata: st_local("left", strofreal(findexternal("seldecomp_JK") != NULL))
    assert `left' == 0
    assert `"`e(cmdline)'"' == `"`cmd0'"'
    di as txt "Break simulated in raking call `=`n' + 1': nothing left behind"
}
drop ord
di as txt "section 10: Break inside the Mata raker passed"

* Numeric focal-group selection must use the selected index, not rounded label text.
preserve
clear
set obs 30
gen byte x = mod(_n, 3)
gen double y = _n
gen double g = cond(_n <= 15, 1, 1 + 2e-16)
quietly seldecomp y x, by(g) focalgroup(1.0000000000000002) mincell(0)
assert e(mean_focal) == 23 & e(mean_ref) == 8 & e(difference) == 15
quietly seldecomp y x, by(g) mincell(0)
assert e(mean_focal) == 8 & e(mean_ref) == 23 & e(difference) == -15
replace g = cond(_n <= 15, 10000000000000000, 10000000000000002)
quietly seldecomp y x, by(g) focalgroup(10000000000000002) mincell(0)
assert e(mean_focal) == 23 & e(mean_ref) == 8 & e(difference) == 15
restore

* 11. estat bootstrap after a replay at a level that has no exact binary form; a raking
*     variable that is constant in the sample when insupport excludes observations; a
*     number typed with a plus sign; the message when the one observation left has zero
*     weight
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
quietly seldecomp y x1 x2 [pw=wt] if id <= 300, by(group) bayesboot(reps(20) rseed(7) nodots)
foreach L in 99.9 80.1 66.6 99.99 97.3 {
    quietly seldecomp, level(`L')
    rv_log "estat bootstrap"
    assert r(rc) == 0
    rv_has "level() can have" 0
    rv_has "`L'% conf. int"
}
gen byte one = 1
gen byte xs = x1
replace xs = 9 if group == 0 & mod(id, 40) == 0
rv_log "seldecomp y one xs [pw=wt], by(group) insupport"
assert r(rc) == 0
rv_has "note: one is constant in the estimation sample"
tempname b1
matrix `b1' = e(b)
assert e(N_terms) == 2
quietly seldecomp y xs [pw=wt], by(group) insupport
assert mreldif(e(b), `b1') < 1e-12
quietly seldecomp y (one) (xs x2) [pw=wt], by(group) insupport sequentially fixothermargins
assert e(N_stages) == 2
quietly seldecomp y x1 x2 [pw=wt], by(group) trimhirel(50)
matrix `b1' = e(b)
quietly seldecomp y x1 x2 [pw=wt], by(group) tolerance(+1e-8) ctrltolerance(+1e-6) trimhirel(+50)
assert mreldif(e(b), `b1') == 0
assert "`e(rakeopts)'" == "iterate(10000) tolerance(1e-8) ctrltolerance(1e-6) trimhirel(50)"
gen double wz = wt
replace wz = 0 if id == 5
rv_log "seldecomp y x1 x2 [pw=wz] if id == 5, by(group)"
assert r(rc) == 2000
rv_has "no observations: the one observation with nonmissing data has zero weight"
replace wz = 0 if id == 6
rv_log "seldecomp y x1 x2 [pw=wz] if id == 5 | id == 6, by(group)"
assert r(rc) == 2000
rv_has "no observations: all 2 observations with nonmissing data have zero weight"

* 12. a raking variable that is constant in the sample when the raking fails for support:
*     the message and return code 459, with and without stages; a constant raking variable
*     alone under svy jackknife; the number of quantile groups typed with a plus sign
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
gen byte one = 1
gen byte xr = x1
replace xr = 8 if group == 1 & mod(id, 50) == 0
rv_log "seldecomp y one xr [pw=wt], by(group)"
assert r(rc) == 459
rv_has "raking support failure: these levels of the reference group"
rv_has "xr = 8"
rv_has "conformability" 0
rv_log "seldecomp y one xr x2 [pw=wt], by(group) sequentially"
assert r(rc) == 459
rv_has "raking support failure at stage 2"
rv_has "conformability" 0
rv_log "seldecomp y one xr [pw=wt], by(group) focalgroup(1)"
assert r(rc) == 0
preserve
svyset cluster [pw=wt]
rv_log "svy jackknife: seldecomp y one, by(group)"
assert r(rc) == 0
rv_has "conformability" 0
restore
quietly seldecomp y x1 c.x3 [pw=wt], by(group) contvars(x3 4)
tempname b2
matrix `b2' = e(b)
rv_log "seldecomp y x1 c.x3 [pw=wt], by(group) contvars(x3 +4)"
assert r(rc) == 0
assert mreldif(e(b), `b2') == 0
assert "`e(cutq)'" == "4"
rv_has "x3: 4 quantile bins"
rv_has "+4" 0

capture mata: mata drop _rv_txt _rv_lines
di as res _n "PASSED: test_edge_cases.do"

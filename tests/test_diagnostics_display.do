* Diagnostics visibility, status/ordering, saved-result replay and prefixes.
* Numerical target counting is independently tested in test_diagnostic_counts.do.
version 16
clear
set more off
set varabbrev off
set linesize 78
discard
if `"$SELDECOMP_ROOT"' != "" adopath ++ `"$SELDECOMP_ROOT"'
else adopath ++ `"`c(pwd)'/.."'

capture program drop dr_log
program define dr_log
    gettoken fn 0 : 0
    log using `"`fn'"', text replace name(dr_test) nomsg
    capture noisily `0'
    local rc = _rc
    log close dr_test
    if `rc' exit `rc'
end
capture program drop dr_find
program define dr_find, rclass
    args fn txt
    mata: _dr_l = cat(st_local("fn"))
    mata: _dr_i = select(1::rows(_dr_l), strpos(_dr_l, st_local("txt")) :> 0)
    mata: st_local("n", strofreal(rows(_dr_i)))
    mata: st_local("i", strofreal(rows(_dr_i) ? _dr_i[1] : 0))
    return scalar n = `n'
    return scalar i = `i'
    if `n' {
        mata: st_local("line", _dr_l[_dr_i[1]])
        return local line `"`macval(line)'"'
    }
    mata: mata drop _dr_l _dr_i
end
capture program drop dr_visibility
program define dr_visibility
    args fn diag fail factor
    dr_find `"`fn'"' "Raking diagnostics"
    assert r(n) == `diag'
    dr_find `"`fn'"' "Raking failure counts"
    assert r(n) == `fail'
    dr_find `"`fn'"' "Raking factors"
    assert r(n) == `factor'
    mata: _dr_s = cat(st_local("fn")); _dr_k = select(1::rows(_dr_s), regexm(_dr_s, "^Raking(:| [(])"))
    mata: st_local("status_n", strofreal(rows(_dr_k)))
    mata: st_local("status_i", strofreal(rows(_dr_k) ? _dr_k[1] : 0))
    mata: mata drop _dr_s _dr_k
    assert `status_n' == 1
    dr_find `"`fn'"' "Focal mean |"
    if `diag' & !`fail' assert `status_i' < r(i)
    else assert `status_i' > r(i)
    if `factor' {
        dr_find `"`fn'"' "Raking factors"
        assert r(i) == `status_i' + 1
    }
    foreach old in "difference: focal_group" "explained: focal_group" "unexplained: raked_focal" {
        dr_find `"`fn'"' "`old'"
        assert r(n) == 0
    }
end

tempfile out saved
sysuse nlsw88, clear
quietly seldecomp wage race collgrad, by(married)
matrix dr_b = e(b)
quietly datasignature
local sig `"`r(datasignature)'"'
* nodiagnostics, fulldiagnostics and their abbreviations are refused (198), data and e() unchanged.
foreach o in nodiagnostics nodiag nodiagn fulldiagnostics fulldiag {
    capture seldecomp wage race, by(married) `o' generate(should_not_exist)
    assert _rc == 198
    assert mreldif(e(b), dr_b) == 0
    capture confirm variable should_not_exist, exact
    assert _rc == 111
    capture seldecomp, `o'
    assert _rc == 198
    assert mreldif(e(b), dr_b) == 0
    quietly datasignature
    assert `"`r(datasignature)'"' == `"`sig'"'
}
* Clean default/explicit displays and replay retain the same estimates/checks.
dr_log `"`out'"' seldecomp wage race collgrad, by(married)
dr_visibility `"`out'"' 0 0 1
matrix dr_c = e(raking_checks)
assert dr_c[1,1] == 1 & dr_c[1,2] == 0 & dr_c[1,3] == 0 & dr_c[1,4] == 5
foreach label in "Focal mean |" "Raked focal mean |" "Reference mean |" "Overall" {
    dr_find `"`out'"' "`label'"
    assert r(n) == 1
}
dr_log `"`out'"' seldecomp, diagnostics
dr_visibility `"`out'"' 1 0 0
assert mreldif(e(b), dr_b) == 0 & mreldif(e(raking_checks), dr_c) == 0
foreach label in "Zero wt: zero raked weights." "At bound: distinct weights at a trimming bound." {
    dr_find `"`out'"' "`label'"
    assert r(n) == 1
}
dr_find `"`out'"' "Raking diagnostics"
local before = r(i)
dr_find `"`out'"' "Decomposition"
assert `before' > 0 & `before' < r(i)

* Success and target failures, without and with stages. Once-only trimming leaves
* iteration convergence true while clipping target totals for BOTH variables (rc 430
* without force; force keeps the results).
foreach stage in "" sequentially {
    dr_log `"`out'"' seldecomp wage race collgrad, by(married) `stage'
    dr_visibility `"`out'"' 0 0 1
    matrix dr_s = e(bystage)
    mata: st_local("lo", strtrim(strofreal(min(st_matrix("dr_s")[,16]), "%9.4g")))
    mata: st_local("hi", strtrim(strofreal(max(st_matrix("dr_s")[,17]), "%9.4g")))
    local thin = dr_s[rowsof(dr_s),26]
    dr_find `"`out'"' "`lo'-`hi'."
    assert r(n) == 1
    dr_find `"`out'"' "Thin target levels (<`e(mincell)' focal obs): `thin'."
    assert r(n) == 1
    rcof "noi seldecomp wage race collgrad, by(married) `stage' trimhirel(.8) trimfrequency(once)" == 430
    dr_log `"`out'"' seldecomp wage race collgrad, by(married) `stage' ///
        trimhirel(.8) trimfrequency(once) force
    dr_visibility `"`out'"' 1 1 0
    matrix dr_c = e(raking_checks)
    local ns = rowsof(dr_c)
    assert dr_c[`ns',1] == 1 & dr_c[`ns',2] == 2 & dr_c[`ns',3] == 4
    assert e(converged) == 0 & e(targets_met) == 0
    dr_find `"`out'"' "*not*"
    assert r(n) >= 1
    dr_find `"`out'"' "note: trimming left"
    assert r(n) == 0
    matrix dr_m = e(raking_misses)
    matrix dr_b = e(b)
    estimates save `"`saved'"', replace
    preserve
    clear
    estimates use `"`saved'"'
    dr_log `"`out'"' seldecomp
    dr_visibility `"`out'"' 1 1 0
    assert mreldif(e(b), dr_b) == 0 & mreldif(e(raking_checks), dr_c) == 0
    assert mreldif(e(raking_misses), dr_m) == 0
    restore
}
* A real iteration failure: stage 1 meets its targets but stage 2 does not.
dr_log `"`out'"' seldecomp wage race collgrad, by(married) sequentially iterate(1) force
dr_visibility `"`out'"' 1 1 0
matrix dr_c = e(raking_checks)
assert dr_c[1,1] == 0 & dr_c[1,3] == 0
assert dr_c[2,1] == 0 & dr_c[2,2] == 1 & dr_c[2,3] == 3
assert el(e(bystage), 1, 8) > 1e-8
assert el(e(bystage), 2, 8) > 1e-8

* Inference preserves all stage equations. Reverse and term detail have their own
* diagnostics, restored alongside the primary matrices.
dr_log `"`out'"' seldecomp wage (demo: race collgrad) (geo: south smsa), ///
    by(married) sequentially reverse detail diagnostics bayesboot(reps(5) rseed(571) nodots)
matrix dr_b = e(b)
matrix dr_rev = e(raking_checks_rev)
matrix dr_det = e(raking_checks_detail)
dr_find `"`out'"' "By stage"
assert r(n) >= 8
dr_find `"`out'"' "By term"
assert r(n) >= 1
dr_find `"`out'"' "reverse order"
assert r(n) >= 1
dr_find `"`out'"' "P>|z|"
assert r(n) == 1
estimates store dr_saved
quietly regress wage grade
estimates restore dr_saved
dr_log `"`out'"' seldecomp, diagnostics level(90)
assert mreldif(e(b), dr_b) == 0
assert mreldif(e(raking_checks_rev), dr_rev) == 0
assert mreldif(e(raking_checks_detail), dr_det) == 0
dr_find `"`out'"' "[90% conf. interval]"
assert r(n) == 1

* Prefix redisplay propagates diagnostics, and ordinary replay returns to default.
dr_log `"`out'"' bootstrap, reps(5) seed(101) nodots: ///
    seldecomp wage race collgrad, by(married) diagnostics
assert "`e(prefix_dispopts)'" == "diagnostics"
dr_visibility `"`out'"' 1 0 0
matrix dr_c = e(raking_checks)
assert dr_c[1,4] == 5
dr_log `"`out'"' seldecomp
dr_visibility `"`out'"' 0 0 1
assert mreldif(e(raking_checks), dr_c) == 0

* Chain replay keeps both pairs' checks and summarizes each failure panel.
clear
set obs 300
gen byte g = floor((_n - 1)/100)
gen byte x = mod(_n,3)
gen byte z = mod(_n,2)
gen double y = 1 + x + z + g
dr_log `"`out'"' seldecomp y x z, by(g) chain sequentially
assert e(N_pairs) == 2
matrix dr_c = e(raking_checks)
matrix dr_c2 = e(raking_checks2)
dr_find `"`out'"' "Raking diagnostics"
assert r(n) == 0
dr_log `"`out'"' seldecomp, diagnostics
assert mreldif(e(raking_checks), dr_c) == 0 & mreldif(e(raking_checks2), dr_c2) == 0
dr_find `"`out'"' "Pair p1"
assert r(n) == 1
dr_find `"`out'"' "Pair p2"
assert r(n) == 1
dr_log `"`out'"' seldecomp y x z, by(g) chain sequentially trimhirel(.8) trimfrequency(once) force
dr_find `"`out'"' "Raking failure counts"
assert r(n) == 2
dr_find `"`out'"' "Raking factors"
assert r(n) == 0
dr_find `"`out'"' "*not*"
assert r(n) >= 1
capture estimates drop dr_saved
di as result "PASSED: test_diagnostics_display.do"

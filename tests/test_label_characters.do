* Literal labels in tables/replay/graphs, and unconditional reporting of available
* explained-share uncertainty. Run from tests/, or set SELDECOMP_ROOT.
version 16
clear
set more off
discard
if `"$SELDECOMP_ROOT"' != "" adopath ++ `"$SELDECOMP_ROOT"'
else adopath ++ `"`c(pwd)'/.."'

capture program drop df_log
program define df_log
    gettoken fn 0 : 0
    log using `"`fn'"', text replace name(df_test) nomsg
    capture noisily `0'
    local rc = _rc
    log close df_test
    if `rc' exit `rc'
end
capture program drop df_contains
program define df_contains, rclass
    args fn txt
    mata: st_local("found", strofreal(any(strpos(cat(st_local("fn")), st_local("txt")) :> 0)))
    return scalar found = `found'
end
tempfile output

* Quoted outcome/term/group labels are data, including macro-like text and SMCL braces.
set obs 60
gen byte g = _n > 30
gen byte x = mod(_n, 3)
gen byte z = mod(_n, 2)
gen double y = g + x + z / 10
gen double y2 = 2 * y
gen double den = 1 + z / 5
mata: st_varlabel("y", "Income " + char(34) + "gross" + char(34))
mata: st_varlabel("y2", "Other " + char(34) + "mean" + char(34))
mata: st_varlabel("x", char(34) + "X" + char(34) + " {new}")
mata: st_varlabel("z", "Z " + char(96) + "q'")
mata: st_varlabel("g", "Group " + char(34) + char(36) + "DFBAD" + char(34))
global DFBAD SHOULD_NOT_APPEAR
quietly seldecomp y x z, by(g) sequentially
matrix df_b = e(b)
df_log `"`output'"' seldecomp, varlabels
assert mreldif(e(b), df_b) == 0
df_contains `"`output'"' `"Income "gross""'
assert r(found)
df_contains `"`output'"' "SHOULD_NOT_APPEAR"
assert !r(found)
mata: st_local("labelwant", "Group " + char(34) + char(36) + "DFBAD" + char(34))
df_contains `"`output'"' `"`macval(labelwant)'"'
assert r(found)
mata: st_local("labelwant", "Z " + char(96) + "q'")
df_contains `"`output'"' `"`macval(labelwant)'"'
assert r(found)
df_log `"`output'"' seldecomp y x z, by(g) sequentially varlabels mincell(60) diagnostics
df_contains `"`output'"' `""X" {new}"'
assert r(found)
df_contains `"`output'"' "{c -(}"
assert !r(found)
df_log `"`output'"' seldecomp y x z, by(g) sequentially varlabels ///
    outcomes(y2) ratio(den) bayesboot(reps(10) rseed(37) nodots citable)
df_log `"`output'"' seldecomp, varlabels citable
df_log `"`output'"' estat bootstrap, all

* A saved result can be replayed after loading data with a quoted label on the same name.
estimates store df_saved
df_log `"`output'"' bootstrap, reps(5) seed(41) nodots: ///
    seldecomp y x z, by(g) sequentially varlabels
df_contains `"`output'"' `"Income "gross""'
assert r(found)
estimates restore df_saved
preserve
clear
set obs 1
gen double y = 0
mata: st_varlabel("y", "Income " + char(34) + "gross" + char(34))
estimates restore df_saved
df_log `"`output'"' seldecomp, varlabels
restore

* All graph types keep quoted labels as one text line; no subtitle options become text.
mata: st_vlmodify("df_groups", (0 \ 1), (char(34) + "Single" + char(34) + " women" \ "Married"))
label values g df_groups
quietly seldecomp y x z, by(g) sequentially
foreach view in stages marginal waterfall {
    seldecomp, makeplot(`view', twoway(name(df_graph, replace))) varlabels
    assert `.df_graph.subtitle.text.arrnels' == 1
    assert `"`.df_graph.subtitle.text[1]'"' == "{c 34}Single{c 34} women - Married"
}
seldecomp, makeplot(waterfalllevels, twoway(name(df_levels, replace))) varlabels
assert `"`.df_levels.xaxis1.major.dlg_ticklabel[1]'"' == "{c 34}Single{c 34} women"
graph drop df_graph df_levels

* A difference whose normal CI includes zero still reports the share's percentile CI.
clear
set obs 20
gen byte g = _n > 10
gen byte x = cond(g, _n > 15, _n == 10)
gen double y = cond(g, 0, cond(x, 100, 1))
df_log `"`output'"' seldecomp y x, by(g) ///
    bayesboot(reps(100) rseed(318) percentile nodots)
assert el(e(ci_normal), 1, colnumb(e(ci_normal), "overall:difference")) < 0
assert el(e(ci_percentile), 1, colnumb(e(ci_percentile), "overall:difference")) > 0
matrix df_ci = e(pct_explained_ci)
scalar df_se = e(pct_explained_se)
foreach txt in "bootstrap std. err." "95% percentile CI:" {
    df_contains `"`output'"' "`txt'"
    assert r(found)
}
foreach txt in "unstable" "the CI of the difference includes 0" {
    df_contains `"`output'"' "`txt'"
    assert !r(found)
}
df_log `"`output'"' seldecomp
assert mreldif(e(pct_explained_ci), df_ci) == 0 & e(pct_explained_se) == df_se
df_contains `"`output'"' "95% percentile CI:"
assert r(found)

* Prefix and chain displays use the same reporting policy, including chain bootstrap CIs.
quietly bootstrap, reps(8) seed(51) nodots: seldecomp y x, by(g)
df_log `"`output'"' seldecomp
df_contains `"`output'"' "unstable"
assert !r(found)
clear
set obs 90
gen byte g = floor((_n - 1) / 30)
gen byte x = mod(_n, 3)
gen double y = 1 + x + .01*g
mata: st_varlabel("y", "Income " + char(34) + "gross" + char(34))
df_log `"`output'"' seldecomp y x, by(g) chain varlabels ///
    bayesboot(reps(20) rseed(62) nodots citable)
df_contains `"`output'"' "95% percentile CI"
assert r(found)
df_contains `"`output'"' "unstable"
assert !r(found)
df_log `"`output'"' seldecomp, varlabels citable
df_contains `"`output'"' "unstable"
assert !r(found)
macro drop DFBAD
capture estimates drop df_saved
di as result "PASSED: test_label_characters.do"

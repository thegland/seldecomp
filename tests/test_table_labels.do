* Descriptive table labels, preserved equation distinctions, and 78-column formatting.
version 16
clear
set more off
set linesize 200
discard
if `"$SELDECOMP_ROOT"' != "" adopath ++ `"$SELDECOMP_ROOT"'
else adopath ++ `"`c(pwd)'/.."'

capture program drop tl_capture
program define tl_capture
    gettoken fn 0 : 0
    log using `"`fn'"', text replace name(tl_test) nomsg
    di "TABLE_LABELS_BEGIN"
    capture noisily `0'
    local rc = _rc
    di "TABLE_LABELS_END"
    log close tl_test
    if `rc' exit `rc'
end
capture program drop tl_assert
program define tl_assert
    args fn text
    mata: st_local("found", strofreal(any(strpos(cat(st_local("fn")), st_local("text")) :> 0)))
    assert `found'
end
capture program drop tl_width
program define tl_width
    args fn
    mata: _tl_l = cat(st_local("fn"))
    mata: _tl_b = select(1::rows(_tl_l), _tl_l :== "TABLE_LABELS_BEGIN")
    mata: _tl_e = select(1::rows(_tl_l), _tl_l :== "TABLE_LABELS_END")
    mata: _tl_l = _tl_l[(_tl_b[1] + 1)::(_tl_e[1] - 1)]
    mata: st_local("width", strofreal(max(udstrlen(_tl_l))))
    assert `width' <= 78
    mata: mata drop _tl_l _tl_b _tl_e
end

set obs 240
gen byte g = _n > 120
gen byte x = mod(_n, 3)
gen byte z = mod(floor((_n - 1)/3), 2)
gen byte q = mod(floor((_n - 1)/6), 2)
gen double y = 1 + x + .5*z + .2*q + .3*g + sin(_n)
gen double y2 = 2*y + cos(_n)
gen double den = 1 + .2*z
mata: st_varlabel("y", "Income " + char(34) + "gross" + char(34))
tempfile out

quietly seldecomp y x z, by(g) sequentially bayesboot(reps(8) rseed(31) nodots)
matrix tl_b = e(b)
matrix tl_V = e(V)
local stored : colfullnames e(b)
tl_capture `"`out'"' seldecomp, varlabels citable
foreach name in "Focal mean" "Raked focal mean" "Reference mean" "Difference" "Explained" "Unexplained" "Overall" "By stage" "Raked mean" {
    tl_assert `"`out'"' "`name'"
}
tl_assert `"`out'"' `"Income "gross""'
mata: _tl_l = cat(st_local("out"))
mata: _tl_l = select(_tl_l, (strpos(_tl_l, "Income") :> 0) :& (strpos(_tl_l, "|") :> 0))
mata: st_local("aligned", strofreal(rows(_tl_l) == 2 & all(strpos(_tl_l, "|") :== 18)))
assert `aligned'
mata: mata drop _tl_l
tl_assert `"`out'"' "Difference = focal - reference. Explained = focal - raked focal."
tl_assert `"`out'"' "Raked focal group mean"
tl_width `"`out'"'
assert mreldif(e(b), tl_b) == 0 & mreldif(e(V), tl_V) == 0
local after : colfullnames e(b)
assert "`stored'" == "`after'"
assert colnumb(e(b), "overall:focal_group") < .
assert colnumb(e(b), "marg_expl:x") < .

tl_capture `"`out'"' seldecomp, citable cformat(%9.3f) pformat(%5.2f) sformat(%8.3f)
tl_width `"`out'"'
tl_capture `"`out'"' estat bootstrap, all
tl_assert `"`out'"' "Raked focal mean"
tl_assert `"`out'"' "By stage"
tl_width `"`out'"'

* Every equation survives outcomes(), detail, reverse order, and ratios.
quietly seldecomp y (demo: x z) q, by(g) sequentially detail reverse ///
    outcomes(y2) ratio(den) bayesboot(reps(8) rseed(32) nodots)
tl_capture `"`out'"' seldecomp, citable
foreach name in "Focal ratio" "Raked ratio" "Reference ratio" "By term" "reverse order" "Raked focal group ratio" {
    tl_assert `"`out'"' "`name'"
}
tl_width `"`out'"'
tl_capture `"`out'"' estat bootstrap, all
tl_assert `"`out'"' "Outcome"
tl_assert `"`out'"' "By term"
tl_assert `"`out'"' "reverse order"
tl_assert `"`out'"' "Raked ratio"
tl_width `"`out'"'

quietly seldecomp y x z, by(g) individually bayesboot(reps(8) rseed(33) nodots)
tl_capture `"`out'"' seldecomp, citable
tl_assert `"`out'"' "By block"
tl_assert `"`out'"' "Contribution of each block raked alone, and the remainder"
tl_width `"`out'"'

* Chain pair names remain visible alongside descriptive equation/row labels.
replace g = floor((_n - 1)/80)
quietly seldecomp y x z, by(g) chain sequentially bayesboot(reps(8) rseed(34) nodots)
tl_capture `"`out'"' seldecomp, citable
foreach name in "p1" "p2" "total" "Overall" "By stage" "Raked focal mean" {
    tl_assert `"`out'"' "`name'"
}
tl_width `"`out'"'
tl_capture `"`out'"' estat bootstrap, all
foreach name in "p1" "p2" "total" "Overall" "By stage" "Raked focal mean" {
    tl_assert `"`out'"' "`name'"
}
tl_width `"`out'"'
set linesize 78
di as result "PASSED: test_table_labels.do"

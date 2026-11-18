* test_bootstrap_cleanup.do -- interrupted bootstrap cleanup; bca alone and for a chain.
* Run from the tests directory, or set SELDECOMP_ROOT to the repository root.

version 16
clear
set more off
set varabbrev on
discard
if `"$SELDECOMP_ROOT"' != "" adopath ++ `"$SELDECOMP_ROOT"'
else adopath ++ `"`c(pwd)'/.."'
which seldecomp

set obs 30
gen byte group = _n > 15
gen byte period = ceil(_n / 10)
gen byte x = mod(_n, 3)
gen double y = _n
gen long ord = _n
gen double bb_2_extra = 42
gen double bb_4_extra = 44
quietly seldecomp y x, by(group) mincell(0)
matrix original_b = e(b)
local original_cmd `"`e(cmdline)'"'
quietly datasignature
local original_data `r(datasignature)'

* The raker's hook stops the call after n further raking calls.
* 1: first replicate; 3: third replicate after partial completion; 6: first BCa deletion.
foreach n in 1 3 6 {
    local bca
    if `n' == 6 local bca bca
    mata: seldecomp_brk = `n'
    capture noisily seldecomp y x, by(group) mincell(0) ///
        bayesboot(reps(5) rseed(123) nodots generate(bb) `bca')
    assert _rc == 1
    mata: st_local("left", strofreal(findexternal("seldecomp_brk") != NULL))
    assert `left' == 0
    confirm variable bb_2_extra bb_4_extra, exact
    assert bb_2_extra == 42 & bb_4_extra == 44
    forvalues j = 1/5 {
        capture confirm variable bb_`j', exact
        assert _rc == 111
    }
    quietly datasignature
    assert "`r(datasignature)'" == "`original_data'"
    assert ord == _n
    assert mreldif(e(b), original_b) == 0
    assert `"`e(cmdline)'"' == `"`original_cmd'"'
}

* With replace, a later existing replicate variable is untouched if the loop never reaches it.
gen double bb_5 = 105
quietly datasignature
local original_data `r(datasignature)'
mata: seldecomp_brk = 1
capture noisily seldecomp y x, by(group) mincell(0) ///
    bayesboot(reps(5) rseed(123) nodots generate(bb, replace))
assert _rc == 1
confirm variable bb_5 bb_2_extra bb_4_extra, exact
assert bb_5 == 105
quietly datasignature
assert "`r(datasignature)'" == "`original_data'"
assert mreldif(e(b), original_b) == 0
assert `"`e(cmdline)'"' == `"`original_cmd'"'

* A refused bayesboot() leaves seeds, generated variables, files, and estimation results
* unchanged, and the chain path refuses it in its initial probe too. The vehicle is the
* suboption strata(), refused (198) when bayesboot() is parsed.
tempfile artifact
local original_rng `"`c(rngstate)'"'
forvalues j = 1/2 {
    local byvar group
    local chain
    if `j' == 2 {
        local byvar period
        local chain chain
    }
    capture noisily seldecomp y x, by(`byvar') `chain' mincell(0) ///
        bayesboot(reps(5) rseed(987) nodots bca strata(x) ///
        generate(bb, replace) saving(`"`artifact'"', replace))
    assert _rc == 198
    quietly datasignature
    assert "`r(datasignature)'" == "`original_data'"
    assert `"`c(rngstate)'"' == `"`original_rng'"'
    assert mreldif(e(b), original_b) == 0
    assert `"`e(cmdline)'"' == `"`original_cmd'"'
    capture confirm file `"`artifact'"'
    assert _rc == 601
}

* bca runs, including for a chain.
quietly seldecomp y x, by(group) mincell(0) ///
    bayesboot(reps(5) rseed(123) nodots bca)
assert e(bca) == 1 & e(N_reps) == 5
confirm matrix e(accel)
quietly seldecomp y x, by(period) chain mincell(0) ///
    bayesboot(reps(5) rseed(123) nodots bca)
assert e(bca) == 1 & e(N_pairs) == 2 & e(N_reps) == 5
confirm matrix e(accel)

set varabbrev off
di as res _n "PASSED: test_bootstrap_cleanup.do"

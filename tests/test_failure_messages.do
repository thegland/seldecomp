* test_failure_messages.do -- concise cautions, with unchanged inference.
version 16
clear
set more off
set varabbrev off
discard
if `"$SELDECOMP_ROOT"' != "" adopath ++ `"$SELDECOMP_ROOT"'
else adopath ++ `"`c(pwd)'/.."'
which seldecomp

capture program drop fm_count
program define fm_count, rclass
    args path needle
    tempname fh
    file open `fh' using `"`path'"', read text
    local full
    file read `fh' line
    while r(eof) == 0 {
        local full `"`macval(full)' `macval(line)'"'
        file read `fh' line
    }
    file close `fh'
    // Paragraph wrapping varies with Results-window width; compare normalized output.
    mata: st_local("full", ustrregexra(st_local("full"), "[[:space:]]+", " "))
    mata: st_local("n", strofreal((strlen(st_local("full")) - strlen(subinstr(st_local("full"), st_local("needle"), ""))) / strlen(st_local("needle"))))
    return scalar n = `n'
end

* One unsupported deletion, although all positive-weight bootstrap draws succeed.
set obs 20
gen long id = _n
gen byte group = id > 10
gen byte x = id == 1 | inrange(id, 11, 14)
gen double y = id^2/100 + 7*x + mod(id, 3)
tempfile msg saved
log using `"`msg'"', text replace name(fm)
noisily seldecomp y x, by(group) mincell(0) ///
    bayesboot(reps(40) rseed(921) bca nodots)
log close fm
assert e(N_jkfail) == 1 & e(N_misreps) == 0 & e(N_reps) == 40
* A caution only: the acceleration is the partial-jackknife one.
assert abs(el(e(accel), 1, 1) - (-.05339335)) < 1e-8
fm_count `"`msg'"' "Caution: 1 jackknife replicates failed; BCa intervals may be incorrect."
assert r(n) == 1
fm_count `"`msg'"' "bootstrap results may be incorrect"
assert r(n) == 0
matrix fm_b = e(b)
matrix fm_V = e(V)
matrix fm_a = e(accel)
matrix fm_ci = e(ci_bca)
estimates save `"`saved'"', replace
estimates clear
estimates use `"`saved'"'
foreach cmd in "seldecomp" "estat bootstrap, bca" {
    log using `"`msg'"', text replace name(fm)
    noisily `cmd'
    log close fm
    fm_count `"`msg'"' "Caution: 1 jackknife replicates failed; BCa intervals may be incorrect."
    assert r(n) == 1
    assert e(N_jkfail) == 1
    assert mreldif(e(b), fm_b) == 0 & mreldif(e(V), fm_V) == 0
    assert mreldif(e(accel), fm_a) == 0 & mreldif(e(ci_bca), fm_ci) == 0
}

* No duplicate caution for multiple outcomes, reverse order, or detailed stages.
gen byte z = 1
gen double y2 = 2*y + 1
log using `"`msg'"', text replace name(fm)
noisily seldecomp y (main: x z), by(group) mincell(0) outcomes(y2) ///
    sequentially reverse detail bayesboot(reps(10) rseed(921) bca nodots)
log close fm
assert e(N_jkfail) == 1
fm_count `"`msg'"' "Caution: 1 jackknife replicates failed; BCa intervals may be incorrect."
assert r(n) == 1

* Chains retain the failed-deletion count and repeat the caution on replay.
clear
set obs 30
gen long id = _n
gen byte period = ceil(id/10)
gen byte x = id == 1 | inrange(id, 11, 14) | inrange(id, 21, 24)
gen double y = id^2/100 + 7*x + mod(id, 3)
log using `"`msg'"', text replace name(fm)
noisily seldecomp y x, by(period) chain mincell(0) ///
    bayesboot(reps(20) rseed(921) bca nodots)
log close fm
assert e(N_jkfail) == 1
fm_count `"`msg'"' "Caution: 1 jackknife replicates failed; BCa intervals may be incorrect."
assert r(n) == 1
log using `"`msg'"', text replace name(fm)
noisily seldecomp
log close fm
fm_count `"`msg'"' "Caution: 1 jackknife replicates failed; BCa intervals may be incorrect."
assert r(n) == 1

* Complete jackknives are quiet, regardless of the progress-dot setting.
replace x = mod(id, 2)
foreach dots in "nodots" "dots(1)" {
    log using `"`msg'"', text replace name(fm)
    noisily seldecomp y x if period <= 2, by(period) mincell(0) ///
        bayesboot(reps(10) rseed(921) bca `dots')
    log close fm
    assert e(N_jkfail) == 0 & e(N_misreps) == 0
    fm_count `"`msg'"' "Caution:"
    assert r(n) == 0
    fm_count `"`msg'"' "bootstrap results may be incorrect"
    assert r(n) == 0
}

* A clustered jackknife reports the actual number of failed cluster deletions.
clear
set obs 400
gen int cl = ceil(_n/40)
gen byte group = mod(_n,2)
gen byte x = 1 + mod(floor(_n/2),2)
replace x = 3 if (group==0 & cl==1) | (group==1 & inlist(cl,5,6))
replace x = 4 if (group==0 & cl==2) | (group==1 & inlist(cl,7,8))
gen double y = x + mod(_n*7,11)/10
log using `"`msg'"', text replace name(fm)
noisily seldecomp y x, by(group) mincell(0) ///
    bayesboot(reps(10) rseed(1) bca cluster(cl) nodots)
log close fm
assert e(N_jkfail) == 2 & e(N_misreps) == 0
fm_count `"`msg'"' "Caution: 2 jackknife replicates failed; BCa intervals may be incorrect."
assert r(n) == 1

* Selected successful draws: the caution is appended to the count/cause note.
clear
set obs 200
gen long id = _n
gen byte group = id > 100
gen int within = cond(group, id-100, id)
gen byte x1 = cond(group, within>29 & within<=80, within>49)
gen byte x2 = cond(group, within>51, within>51)
gen double y = 10*(x1-x2) + sin(id)/10
log using `"`msg'"', text replace name(fm)
noisily seldecomp y x1 x2, by(group) iterate(200) ctrltolerance(.0001) ///
    bayesboot(reps(200) rseed(3401) nodots)
log close fm
assert e(N_misreps) == 100 & e(N_reps) == 100
fm_count `"`msg'"' "bootstrap results may be incorrect"
assert r(n) == 1
fm_count `"`msg'"' "of them missed their targets"
assert r(n) == 1
matrix fm_b = e(b)
matrix fm_V = e(V)
foreach cmd in "seldecomp" "estat bootstrap" {
    log using `"`msg'"', text replace name(fm)
    noisily `cmd'
    log close fm
    fm_count `"`msg'"' "bootstrap results may be incorrect"
    assert r(n) == 1
    assert mreldif(e(b), fm_b) == 0 & mreldif(e(V), fm_V) == 0
}

* Chain failure notes and replay carry the same qualification, once each.
set obs 300
replace id = _n
gen byte period = ceil(id/100)
replace within = mod(id-1,100)+1
replace x1 = within>29 & within<=80 if period==3
replace x2 = within>51 if period==3
replace y = 10*(x1-x2) + sin(id)/10
log using `"`msg'"', text replace name(fm)
noisily seldecomp y x1 x2, by(period) chain iterate(200) ctrltolerance(.0001) ///
    bayesboot(reps(40) rseed(3401) nodots)
log close fm
assert e(N_misreps)>0 & e(N_misreps)<40 & e(N_reps)>=2
fm_count `"`msg'"' "bootstrap results may be incorrect"
assert r(n) == 1
log using `"`msg'"', text replace name(fm)
noisily seldecomp
log close fm
fm_count `"`msg'"' "bootstrap results may be incorrect"
assert r(n) == 1

* A point-estimate iteration note absorbs the replicate count and caution.
drop if period==3
replace x1 = within>29 & within<=82 if group==1
replace x2 = within>53 if group==1
log using `"`msg'"', text replace name(fm)
noisily seldecomp y x1 x2, by(group) iterate(1) ctrltolerance(.08) ///
    bayesboot(reps(40) rseed(3401) nodots)
log close fm
assert e(converged)==0 & e(badcontrols)==0
assert e(N_misreps)==23 & e(N_iterlimit)==17
fm_count `"`msg'"' "bootstrap results may be incorrect"
assert r(n)==1
fm_count `"`msg'"' "also in 40 of 40 replications"
assert r(n)==1
log using `"`msg'"', text replace name(fm)
noisily seldecomp
log close fm
fm_count `"`msg'"' "bootstrap results may be incorrect"
assert r(n)==1
fm_count `"`msg'"' "also in 40 of 40 replications"
assert r(n)==1

* Under a prefix: the one warning on failed replicates, without an additional caution.
clear
set obs 20
gen long id = _n
gen byte group = id > 10
gen byte x = id == 1 | inrange(id, 11, 14)
gen double y = id^2/100 + 7*x + mod(id, 3)
foreach prefix in "bootstrap, reps(40) seed(921) nodots" "jackknife, nodots" {
    log using `"`msg'"', text replace name(fm)
    noisily `prefix': seldecomp y x, by(group) mincell(0)
    log close fm
    assert e(N_misreps) > 0 & e(N_misreps) < .
    local what = cond(strpos("`prefix'", "bootstrap") == 1, "bootstrap", "jackknife")
    fm_count `"`msg'"' "`what' replicates failed and"
    assert r(n) == 1
    fm_count `"`msg'"' "bootstrap results may be incorrect"
    assert r(n) == 0
    fm_count `"`msg'"' "Caution:"
    assert r(n) == 0
}

di as result "PASSED: test_failure_messages.do"

* test_single_level.do -- trimming when every active control is constant.
* Run from tests with the repository on adopath, as for the other test files.
version 16
clear
set more off
set varabbrev off
discard
if `"$SELDECOMP_ROOT"' != "" adopath ++ `"$SELDECOMP_ROOT"'
else adopath ++ `"`c(pwd)'/.."'

* A single retained level still requires the focal group's original weight total.
input byte(g x) double(y w)
0 7 0 1
0 7 10 10
1 7 0 1
1 7 10 1
end
quietly seldecomp y x [pw=w], by(g) generate(rw) mincell(1)
assert rw == w
assert reldif(_b[overall:raked_focal], 100/11) < 1e-12
matrix untrimmed = e(b)
quietly seldecomp y x [pw=w], by(g) trimhiabs(20) generate(rw, replace) mincell(1)
assert rw == w
assert e(N_at_bound) == 0 & e(targets_met) == 1 & e(converged) == 1
assert mreldif(e(b), untrimmed) < 1e-12

* Clipping once has a direct oracle; the other frequencies repeatedly fit the total.
* A bound that leaves the total unmet is rc 430 without force; force keeps the results,
* with e(converged) = 0.
foreach freq in once sometimes often {
    rcof "noi seldecomp y x [pw=w], by(g) trimhiabs(2) trimfrequency(`freq') mincell(1)" == 430
    quietly seldecomp y x [pw=w], by(g) trimhiabs(2) trimfrequency(`freq') generate(rw, replace) mincell(1) force
    assert rw <= 2 + 1e-12 if g == 0
    assert e(targets_met) == 0 & e(badcontrols) == 1 & e(converged) == 0
    assert reldif(rw, 2) < 1e-12 if g == 0 & w == 10
    if "`freq'" == "once" {
        assert rw == 1 if g == 0 & w == 1
        assert reldif(_b[overall:raked_focal], 20/3) < 1e-12
        assert missing(e(maxreldif)) & e(N_at_bound) == 1
    }
    else {
        assert reldif(rw, 2) < 1e-12 if g == 0
        assert reldif(_b[overall:raked_focal], 5) < 1e-12
        assert e(N_at_bound) == 2
    }
    * (trimmed once, the lower bound leaves the total unmet; trimmed in the raking, it is met,
    * and no force is needed)
    capture seldecomp y x [pw=w], by(g) trimloabs(2) trimfrequency(`freq') mincell(1)
    assert _rc == cond("`freq'" == "once", 430, 0)
    quietly seldecomp y x [pw=w], by(g) trimloabs(2) trimfrequency(`freq') generate(rw, replace) mincell(1) force
    assert rw >= 2 - 1e-12 if g == 0
    assert e(N_at_bound) == 1
    if "`freq'" == "once" {
        assert rw == 10 if g == 0 & w == 10
        assert e(targets_met) == 0 & e(badcontrols) == 1 & e(converged) == 0
    }
    else {
        assert reldif(rw, 9) < 1e-7 if g == 0 & w == 10
        assert e(targets_met) == 1 & e(badcontrols) == 0 & e(converged) == 1
    }
    rcof "noi seldecomp y x [pw=w], by(g) trimhirel(.5) trimfrequency(`freq') mincell(1)" == 430
    quietly seldecomp y x [pw=w], by(g) trimhirel(.5) trimfrequency(`freq') generate(rw, replace) mincell(1) force
    assert reldif(rw, .5*w) < 1e-12 if g == 0
    assert e(targets_met) == 0 & e(N_at_bound) == 2 & e(converged) == 0
    rcof "noi seldecomp y x [pw=w], by(g) trimlorel(1.5) trimfrequency(`freq') mincell(1)" == 430
    quietly seldecomp y x [pw=w], by(g) trimlorel(1.5) trimfrequency(`freq') generate(rw, replace) mincell(1) force
    assert reldif(rw, 1.5*w) < 1e-12 if g == 0
    assert e(targets_met) == 0 & e(N_at_bound) == 2 & e(converged) == 0
}
capture noisily seldecomp y x [pw=w], by(g) trimhiabs(2) trimlorel(.5) mincell(1)
assert _rc == 198
capture noisily seldecomp y x [pw=w], by(g) trimloabs(2) trimhirel(1.5) mincell(1)
assert _rc == 198

* Constant-covariate bootstrap draws must be clipped too, checked against direct means.
tempfile draws
quietly seldecomp y x [pw=w], by(g) trimhiabs(2) trimfrequency(once) mincell(1) force ///
    bayesboot(reps(4) rseed(91) generate(bw) saving(`"`draws'"', replace) nodots)
assert e(N_reps) == 4 & e(N_misreps) == 0
preserve
use `"`draws'"', clear
mkmat raked_focal, matrix(raked_draws)
restore
forvalues r = 1/4 {
    generate double clipped = min(bw_`r', 2)
    quietly summarize y [iw=clipped] if g == 0, meanonly
    assert reldif(raked_draws[`r', 1], r(mean)) < 1e-10
    drop clipped
}

* Support screening leaves level 2, not level 1. Targets retain the original W1=12.
clear
input byte(g x) double(y w)
0 2 0 1
0 2 10 10
0 1 20 1
1 2 0 1
1 2 10 1
end
quietly seldecomp y x [pw=w], by(g) generate(rw) mincell(1)
assert reldif(rw, w*12/11) < 1e-12 if g == 0 & x == 2
assert rw == 0 if g == 0 & x == 1
quietly seldecomp y x [pw=w], by(g) trimhiabs(2) trimfrequency(once) generate(rw, replace) mincell(1) force
assert reldif(rw, min(w*12/11, 2)) < 1e-12 if g == 0 & x == 2
assert rw == 0 if g == 0 & x == 1
assert e(targets_met) == 0 & e(badcontrols) == 1 & e(N_at_bound) == 1
quietly seldecomp y x [pw=w], by(g) trimhirel(1) generate(rw, replace) mincell(1) force
assert reldif(rw, w) < 1e-12 if g == 0 & x == 2
assert rw == 0 if g == 0 & x == 1
assert e(targets_met) == 0 & e(N_at_bound) == 2
quietly seldecomp y x [pw=w], by(g) trimloabs(2) trimsupportoverride generate(rw, replace) mincell(1)
assert reldif(rw, 2) < 1e-12 if g == 0 & w == 1 & x == 2
assert reldif(rw, 10) < 1e-7 if g == 0 & w == 10
assert rw == 0 if g == 0 & x == 1
assert e(targets_met) == 1 & e(N_at_bound) == 1
capture noisily seldecomp y x [pw=w], by(g) trimhiabs(2) trimlorel(.5) trimsupportoverride mincell(1)
assert _rc == 198

* The first of two stages is constant; the second has a nonconstant control.
drop if x == 1
replace x = (_n == 2 | _n == 4)
generate byte constant = 7
rcof "noi seldecomp y constant x [pw=w], by(g) sequentially trimhiabs(2) trimfrequency(once) mincell(1)" == 430
quietly seldecomp y constant x [pw=w], by(g) sequentially trimhiabs(2) trimfrequency(once) ///
    generate(stage) mincell(1) force
assert reldif(stage_s1, min(w, 2)) < 1e-12 if g == 0
assert reldif(stage, 2) < 1e-12 if g == 0
matrix stages = e(bystage)
assert stages[1, 7] == 1 & stages[2, 7] == 1
assert stages[1, 22] == 1 & stages[2, 22] == 2

di as result "test_single_level.do PASSED"

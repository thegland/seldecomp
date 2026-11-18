* Diagnostics: unique bound counts and explicit ESS ratio.
version 16
clear all
set more off
set varabbrev off
if `"$SELDECOMP_ROOT"' != "" adopath ++ `"$SELDECOMP_ROOT"'
else adopath ++ `"`c(pwd)'/.."'

* Coincident absolute/relative bounds count each focal observation once.
set obs 20
gen byte g = _n > 10
gen byte x = 1
gen double y = _n
quietly seldecomp y x, by(g) mincell(0)
matrix unbounded = e(b)
foreach freq in once sometimes often {
    foreach bounds in "trimhirel(1) trimhiabs(1)" "trimlorel(1) trimloabs(1)" {
        quietly seldecomp y x, by(g) `bounds' trimfrequency(`freq') diagnostics generate(rw, replace)
        assert e(N_at_bound) == 10
        assert rw == 1
        assert mreldif(e(b), unbounded) < 1e-12
        assert e(ess_ratio) == 1
        local posted : e(scalars)
        assert !strpos(" `posted' ", " N_trimmed ")
        assert !strpos(" `posted' ", " deff_raking ")
        matrix S = e(bystage)
        assert colnumb(S, "N_at_bound") == 22
        assert colnumb(S, "ess_ratio") == 25
        assert missing(colnumb(S, "N_trimmed"))
        assert missing(colnumb(S, "deff_raking"))
    }
}

* Independent union-of-bounds count for nonconstant controls and unequal weights.
clear
set obs 96
gen byte g = _n > 48
gen byte x = mod(_n, 3)
gen byte z = mod(floor((_n-1)/3), 2)
gen double w = 1 + mod(_n, 5)
gen double y = x + z + sin(_n)
foreach freq in once sometimes often {
    quietly seldecomp y x z [pw=w], by(g) sequentially trimhiabs(4) trimhirel(1.2) ///
        trimloabs(.5) trimlorel(.6) trimfrequency(`freq') generate(rw, replace) force
    matrix S = e(bystage)
    forvalues s=1/2 {
        local rv = cond(`s' == 1, "rw_s1", "rw")
        quietly count if g == 0 & (`rv'/4 > 1-1e-6 | `rv'/.5 < 1+1e-6 | ///
            `rv'/w > (1-1e-6)*1.2 | `rv'/w < (1+1e-6)*.6)
        assert S[`s',22] == r(N)
        assert S[`s',22] <= S[`s',20]
    }
    assert e(N_at_bound) == S[2,22]
    quietly summarize w if g == 0, meanonly
    scalar sw = r(sum)
    gen double sq = w^2
    quietly summarize sq if g == 0, meanonly
    scalar eb = sw^2/r(sum)
    quietly summarize rw if g == 0, meanonly
    scalar sr = r(sum)
    replace sq = rw^2
    quietly summarize sq if g == 0, meanonly
    scalar er = sr^2/r(sum)
    assert reldif(e(ess_ratio),eb/er) < 1e-12
    drop sq
}

* Zeroed out-of-support weights are excluded from the final raking sample/count.
clear
input byte(g x) double y
0 1 1
0 1 2
0 2 3
0 2 4
1 1 5
1 1 6
end
quietly seldecomp y x, by(g) trimloabs(2) trimlorel(2) trimsupportoverride generate(rw)
assert e(N_focal_zero) == 2 & e(N_raked) == 2
assert e(N_at_bound) == 2
assert rw == 0 if g == 0 & x == 2
assert rw == 2 if g == 0 & x == 1
display as result "PASSED: test_diagnostic_bounds.do"

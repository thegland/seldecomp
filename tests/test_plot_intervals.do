* Plot interval families and replay levels: compare actual graph coordinates to
* coefficient intervals, including asymmetric/reflected waterfall steps.
version 16
clear
set more off
set varabbrev off
discard
if `"$SELDECOMP_ROOT"' != "" adopath ++ `"$SELDECOMP_ROOT"'
else adopath ++ `"`c(pwd)'/.."'
which seldecomp

capture program drop pi_check
program define pi_check
    syntax name(name=G), TYPE(string) CI(name) [ EQ(string) PFX(string) REVerse ]
    if "`eq'" == "" local eq overall
    tempname B E
    matrix `B' = e(b)
    local blocks `e(blocks)'
    local rp
    if "`reverse'" != "" {
        local original `blocks'
        local blocks
        forvalues s = `: word count `original''(-1)1 {
            local blocks `blocks' `: word `s' of `original''
        }
        local rp rev_
    }
    if "`e(stagetype)'" == "individually" local blocks `blocks' _remainder
    local nb : word count `blocks'
    if "`type'" == "stages" {
        matrix `E' = J(2*`nb', 4, 0)
        forvalues s = 1/`nb' {
            local bn : word `s' of `blocks'
            foreach v in expl unexpl {
                local r = 2*`s' - ("`v'" == "expl")
                local c = colnumb(`ci', "`pfx'`rp'cum_`v':`bn'")
                matrix `E'[`r',1] = `ci'[1,`c']
                matrix `E'[`r',2] = `ci'[2,`c']
                matrix `E'[`r',3] = `s' + cond("`v'" == "expl", -.2, .2)
            }
        }
    }
    else {
        local names
        if "`type'" == "marginal" local names `eq':difference `eq':explained
        if "`type'" == "waterfalllevels" local names `eq':focal_group
        foreach bn of local blocks {
            local names `names' `pfx'`rp'marg_expl:`bn'
        }
        if "`type'" == "marginal" local names `names' `eq':unexplained
        if "`type'" == "waterfall" local names `names' `eq':explained `eq':unexplained `eq':difference
        if "`type'" == "waterfalllevels" local names `names' `eq':raked_focal `eq':unexplained `eq':ref_group
        local nr : word count `names'
        matrix `E' = J(`nr',4,0)
        local previous 0
        forvalues r = 1/`nr' {
            local nm : word `r' of `names'
            local c = colnumb(`ci', "`nm'")
            local b = colnumb(`B', "`nm'")
            local lo = `ci'[1,`c']
            local hi = `ci'[2,`c']
            if "`type'" == "waterfall" {
                local offset 0
                if `r' <= `nb' {
                    local offset `previous'
                    local previous = `previous' + `B'[1,`b']
                }
                if `r' == `nb'+2 local offset = `B'[1,colnumb(`B',"`eq':explained")]
                local lo = `lo' + `offset'
                local hi = `hi' + `offset'
            }
            if "`type'" == "waterfalllevels" {
                if inrange(`r',2,`nb'+1) | `r' == `nb'+3 {
                    local hi = `previous' - `ci'[1,`c']
                    local lo = `previous' - `ci'[2,`c']
                    local previous = `previous' - `B'[1,`b']
                }
                else local previous = `B'[1,`b']
            }
            matrix `E'[`r',1] = `lo'
            matrix `E'[`r',2] = `hi'
            matrix `E'[`r',3] = cond("`type'" == "marginal", ///
                cond(`r'==1,1,`r'+cond(`r'==`nr',.4,.2)), `r')
        }
    }
    local found 0
    preserve
    forvalues k = 1/30 {
        local P .`G'.plotregion1.plot`k'
        if "``P'.classname'" == "" continue, break
        if "``P'.type.snm'" != "rcap" continue
        local id ``P'.serset.id'
        qui serset set `id'
        qui serset use, clear
        mata: for (i=1; i<=st_nvar(); i++) st_varrename(i,"pi_v"+strofreal(i))
        local lo pi_v``P'.yvar'
        local hi pi_v``P'.y2var'
        local pos pi_v``P'.xvar'
        forvalues j = 1/`=_N' {
            if missing(`lo'[`j'],`hi'[`j']) continue
            local match 0
            forvalues r = 1/`=rowsof(`E')' {
                if reldif(`pos'[`j'],`E'[`r',3]) < 1e-12 local match `r'
            }
            assert `match' > 0
            assert reldif(`lo'[`j'],`E'[`match',1]) < 1e-9
            assert reldif(`hi'[`j'],`E'[`match',2]) < 1e-9
            matrix `E'[`match',4] = `E'[`match',4] + 1
            local ++found
        }
    }
    restore
    local expected 0
    forvalues r = 1/`=rowsof(`E')' {
        if `E'[`r',1] < . & `E'[`r',2] < . {
            assert `E'[`r',4] == 1
            local ++expected
        }
    }
    assert `found' == `expected' & `found' > 0
end

* Compare label anchors with the actual bar endpoints and CI spikes. A percentile
* or BC interval need not contain the estimate: labels must clear both objects.
capture program drop pi_check_anchors
program define pi_check_anchors
    syntax name(name=G), TYPE(string) [ SWAPaxis ]
    local horizontal = (("`type'" == "marginal") != ("`swapaxis'" != ""))
    tempname E
    matrix `E' = J(100,4,.)
    local nr 0
    local checked 0
    local inward 0
    local positive 0
    local negative 0
    preserve
    forvalues k = 1/30 {
        local P .`G'.plotregion1.plot`k'
        if "``P'.classname'" == "" continue, break
        local kind ``P'.type.snm'
        if !inlist("`kind'","rbar","rcap","scatter") continue
        local id ``P'.serset.id'
        qui serset set `id'
        qui serset use, clear
        mata: for (i=1; i<=st_nvar(); i++) st_varrename(i,"pa_v"+strofreal(i))
        local y pa_v``P'.yvar'
        local x pa_v``P'.xvar'
        if inlist("`kind'","rbar","rcap") {
            local y2 pa_v``P'.y2var'
            forvalues j = 1/`=_N' {
                if missing(`y'[`j'],`y2'[`j'],`x'[`j']) continue
                if "`kind'" == "rbar" {
                    local ++nr
                    matrix `E'[`nr',1] = `x'[`j']
                    matrix `E'[`nr',2] = `y2'[`j']
                }
                else {
                    local match 0
                    forvalues r = 1/`nr' {
                        if reldif(`x'[`j'],`E'[`r',1]) < 1e-12 local match `r'
                    }
                    assert `match' > 0
                    matrix `E'[`match',3] = `y'[`j']
                    matrix `E'[`match',4] = `y2'[`j']
                }
            }
        }
        else {
            local p pa_v``P'.pvar'
            local pos = cond(`horizontal',"`y'","`x'")
            local anchor = cond(`horizontal',"`x'","`y'")
            forvalues j = 1/`=_N' {
                if missing(`pos'[`j'],`anchor'[`j']) continue
                local match 0
                forvalues r = 1/`nr' {
                    if reldif(`pos'[`j'],`E'[`r',1]) < 1e-12 local match `r'
                }
                assert `match' > 0
                local end = `E'[`match',2]
                local up = inlist(`p'[`j'],3,12)
                local bound = cond(`up',max(`end',`E'[`match',4]),min(`end',`E'[`match',3]))
                assert cond(`up',`anchor'[`j'] >= `bound'-1e-9,`anchor'[`j'] <= `bound'+1e-9)
                * Stages may move labels farther out to clear the total-difference line.
                if "`type'" != "stages" assert reldif(`anchor'[`j'],`bound') < 1e-9
                local inward = `inward' + cond(`up',`E'[`match',4] < `end',`E'[`match',3] > `end')
                local positive = `positive' + `up'
                local negative = `negative' + !`up'
                local ++checked
            }
        }
    }
    restore
    assert `checked' == `nr' & `inward' > 0
    assert `positive' > 0 & `negative' > 0
end

capture program drop pi_inward_intervals
program define pi_inward_intervals, eclass
    tempname B C
    matrix `B' = e(b)
    matrix `C' = e(ci_percentile)
    local names : colfullnames `C'
    forvalues j=1/`=colsof(`C')' {
        local c = colnumb(`B',"`: word `j' of `names''")
        local b = `B'[1,`c']
        matrix `C'[1,`j'] = min(.1*`b',.2*`b')
        matrix `C'[2,`j'] = max(.1*`b',.2*`b')
    }
    ereturn matrix ci_percentile = `C'
end

set seed 4721
set obs 120
gen long id = _n
gen byte g = _n > 60
gen byte x = mod(_n,3)
gen byte z = mod(floor((_n-1)/3),2)
gen double y = exp(rnormal()/2) + x*(1-g/2) + z + g
gen double y2 = -y + rnormal()/4
gen byte period = floor((_n-1)/40)
gen long cluster = ceil(_n/4)
tempfile original saved
save `original'

* Store independently estimated 80% and 95% limits on identical draws. The table
* and graph default must retain the first explicitly typed inference interval.
foreach lv in 80 95 {
    quietly seldecomp y x z, by(g) sequentially reverse ///
        bayesboot(reps(99) rseed(19) percentile bca normal level(`lv') nodots)
    assert "`e(citypes)'" == "percentile bca normal"
    foreach ty in normal t percentile bc bca {
        matrix PI`lv'_`ty' = e(ci_`ty')
    }
}
matrix PIb=e(b)
matrix PIV=e(V)

* Check percentile endpoints directly against empirical quantiles of saved draws.
matrix PIdraws=e(bsdraws)
preserve
clear
svmat double PIdraws, names(draw)
forvalues j=1/`=colsof(PIdraws)' {
    quietly _pctile draw`j', p(10 90)
    assert reldif(r(r1),PI80_percentile[1,`j']) < 1e-12
    assert reldif(r(r2),PI80_percentile[2,`j']) < 1e-12
}
restore

foreach ty in normal t percentile bc bca {
    local option = cond("`ty'"=="t","tinterval","`ty'")
    foreach shape in stages marginal waterfall waterfalllevels {
        quietly seldecomp, makeplot(`shape', ci `option' twoway(name(pifigure, replace)))
        pi_check pifigure, type(`shape') ci(PI95_`ty')
        quietly seldecomp, level(80) makeplot(`shape', ci `option' swapaxis reverse twoway(name(pifigure, replace)))
        pi_check pifigure, type(`shape') ci(PI80_`ty') reverse
    }
}
quietly seldecomp, level(80) makeplot(marginal, ci twoway(name(pifigure,replace)))
pi_check pifigure, type(marginal) ci(PI80_percentile)
assert `"`.pifigure.note.text[1]'"' == "80% percentile intervals"
assert mreldif(PIb,e(b))==0 & mreldif(PIV,e(V))==0
cf _all using `original'

* Saved-result replay needs no original data, and does not change e(level).
estimates save `saved'
preserve
clear
estimates use `saved'
quietly seldecomp, level(80) makeplot(waterfall, ci bca twoway(name(pifigure,replace)))
pi_check pifigure, type(waterfall) ci(PI80_bca)
assert e(level)==95
restore

* Invalid requests are errors; no unavailable intervals silently disappear.
capture seldecomp, makeplot(marginal, ci normal bc)
assert _rc==198
capture seldecomp, makeplot(marginal, percentile)
assert _rc==198
capture seldecomp, level(0) makeplot(marginal, ci)
assert _rc==198
quietly seldecomp y x z, by(g) sequentially bayesboot(reps(39) percentile rseed(2) nodots)
capture seldecomp, makeplot(marginal, ci bca)
assert _rc==198
matrix PIsaved=e(b)
capture seldecomp y x z, by(g) sequentially bayesboot(reps(39) nodots) makeplot(marginal, ci bca)
assert _rc==198
assert mreldif(PIsaved,e(b))==0
quietly seldecomp y x z, by(g) sequentially
capture seldecomp, makeplot(marginal, ci)
assert _rc==198
capture seldecomp y x z, by(g) sequentially makeplot(marginal, ci)
assert _rc==198

* Explicit order with citable is retained; cluster default and overrides agree
* with their inference matrices, including when the table selects percentile.
quietly seldecomp y x z, by(g) sequentially bayesboot(reps(39) bc percentile citable rseed(2) nodots)
assert "`e(citypes)'"=="bc percentile normal"
quietly seldecomp y x z, by(g) sequentially bayesboot(reps(39) cluster(cluster) rseed(2) nodots)
assert "`e(citypes)'"=="t"
matrix PIcluster=e(ci_t)
quietly seldecomp, makeplot(marginal, ci twoway(name(pifigure,replace)))
pi_check pifigure, type(marginal) ci(PIcluster)
quietly seldecomp y x z, by(g) individually bayesboot(reps(39) percentile rseed(2) nodots)
matrix PIindividual=e(ci_percentile)
quietly seldecomp, level(95) makeplot(marginal, ci twoway(name(pifigure,replace)))
pi_check pifigure, type(marginal) ci(PIindividual)
quietly seldecomp y x z, by(g) sequentially outcomes(y2) bayesboot(reps(39) percentile rseed(2) nodots)
matrix PIout=e(ci_percentile)
quietly seldecomp, makeplot(waterfalllevels, ci outcome(y2) twoway(name(pifigure,replace)))
pi_check pifigure, type(waterfalllevels) ci(PIout) eq(y2) pfx(y2_)
quietly seldecomp y x z, by(period) chain sequentially bayesboot(reps(39) percentile rseed(2) nodots)
matrix PIchain=e(ci_percentile)
quietly seldecomp, makeplot(waterfall, ci twoway(name(pifigure,replace)))
pi_check pifigure, type(waterfall) ci(PIchain) eq(total) pfx(total_)
quietly seldecomp, makeplot(stages, ci chainstep(2) twoway(name(pifigure,replace)))
pi_check pifigure, type(stages) ci(PIchain) eq(p2) pfx(p2_)

* Native prefix intervals are usable when stored; a new level cannot be
* reconstructed from endpoints alone. Normal and t still use covariance.
quietly bootstrap, reps(39) seed(7) nodots: seldecomp y x z, by(g) sequentially
matrix PIprefix=e(ci_percentile)
quietly seldecomp, makeplot(marginal, ci percentile twoway(name(pifigure,replace)))
pi_check pifigure, type(marginal) ci(PIprefix)
capture seldecomp, level(80) makeplot(marginal, ci percentile)
assert _rc==198
quietly seldecomp, level(80) makeplot(marginal, ci normal)
quietly seldecomp, level(80) makeplot(marginal, ci tinterval)
capture seldecomp, makeplot(marginal, ci bca)
assert _rc==198

* Deterministic rendering fixture: contributions have both signs, and stored
* percentile limits lie entirely between zero and each estimate. These manually
* replaced limits exercise legal non-Wald interval geometry, not inference.
use `original', clear
replace x = mod(_n,10) < cond(g,7,3)
replace z = mod(floor((_n-1)/10),6) < cond(g,2,4)
replace y = 10 + 3*x + 5*z + 2*g
quietly seldecomp y x z, by(g) sequentially reverse bayesboot(reps(19) percentile rseed(3) nodots)
pi_inward_intervals
matrix PIinward=e(ci_percentile)
foreach shape in stages marginal waterfall waterfalllevels {
    foreach swap in none swapaxis {
        local axisopt = cond("`swap'"=="none","","swapaxis")
        quietly seldecomp, makeplot(`shape', ci percentile `axisopt' twoway(name(pifigure,replace)))
        pi_check pifigure, type(`shape') ci(PIinward)
        pi_check_anchors pifigure, type(`shape') `axisopt'
        quietly seldecomp, makeplot(`shape', ci percentile `axisopt' reverse twoway(name(pifigure,replace)))
        pi_check pifigure, type(`shape') ci(PIinward) reverse
        pi_check_anchors pifigure, type(`shape') `axisopt'
    }
}
capture graph drop pifigure
display as result "PASSED: test_plot_intervals.do"

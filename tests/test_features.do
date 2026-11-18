* test_features.do -- ratio(denvar), outcomes(varlist) and reverse:
* each against the plain runs it must reproduce (the ratio of the means
* computed from the generate() weights; each outcome against its own single-outcome run;
* the reverse order against a run with the blocks typed in reverse order, with the blocks
* not yet added held (sequentially fixothermargins) and free (sequentially)), with and
* without bayesboot() (which bootstraps every stage), under the prefixes (by,
* bootstrap, svy jackknife), at replay and after estimates restore, the display, makeplot,
* and the error paths, those of the mode options included (reverse requires sequentially and
* refuses individually; fixothermargins requires a mode; sequentially and individually exclude
* each other; varbyvar, cumulative, holdothers, isolated, nohold and refminusfocal are refused).
*
* Run from the tests directory:   . do test_features.do

version 16
clear
set more off
set varabbrev off
discard             // drop any seldecomp already in memory; the repository copy is what gets tested

* ---- tolerances ---------------------------------------------------------------
local tol   1e-6     // means from other commands (summarize) vs seldecomp
local same  1e-12    // the same computation reached by two routes (the features' contract)
local rtol  1e-7     // final stages of the two orders: the same targets reached along two
                     // raking paths, so they agree to the raking's tolerance(1e-8)

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
capture program drop assert_close
program define assert_close
    // assert_close exp1 exp2 [tol] ["label"]  --  asserts reldif(exp1, exp2) <= tol
    args a b tol name
    if "`tol'" == "" local tol 1e-6
    if "`name'" == "" local name "`a' vs `b'"
    if missing(`a') | missing(`b') {
        di as err "`name': missing value (`a' = " `a' ", `b' = " `b' ")"
        exit 9
    }
    if reldif(`a', `b') > `tol' {
        di as err "`name': " %18.0g `a' " vs " %18.0g `b' ///
            "  (reldif = " %9.3g reldif(`a', `b') " > `tol')"
        exit 9
    }
end

capture program drop mat_close
program define mat_close
    // mat_close A B [tol]  --  same dimensions and every element equal within tol (an
    // element missing in both is equal)
    args A B tol
    if "`tol'" == "" local tol 1e-12
    if rowsof(`A') != rowsof(`B') | colsof(`A') != colsof(`B') {
        di as err "mat_close: `A' is " rowsof(`A') " x " colsof(`A') ", `B' is " rowsof(`B') " x " colsof(`B')
        exit 9
    }
    forvalues i = 1/`=rowsof(`A')' {
        forvalues j = 1/`=colsof(`A')' {
            if missing(`A'[`i', `j']) & missing(`B'[`i', `j']) continue
            if reldif(`A'[`i', `j'], `B'[`i', `j']) > `tol' | missing(reldif(`A'[`i', `j'], `B'[`i', `j'])) {
                di as err "mat_close: [`i', `j'] " %18.0g `A'[`i', `j'] " vs " %18.0g `B'[`i', `j']
                exit 9
            }
        }
    }
end

capture program drop match_b
program define match_b, rclass
    // match_b A B, map(from=to ...) [tol(#)]  --  A and B are row vectors with equation
    // stripes; every column eq:name of B whose equation is a "from" in map() must equal the
    // column to:name of A; the other columns of B are skipped. r(n): columns compared
    syntax namelist(min=2 max=2), map(string) [ tol(real 1e-12) ]
    gettoken A B : namelist
    local B = strtrim("`B'")
    local bn : colfullnames `B'
    local n 0
    forvalues j = 1/`=colsof(`B')' {
        local fn : word `j' of `bn'
        gettoken eq nm : fn, parse(":")
        gettoken colon nm : nm, parse(":")
        local to
        foreach pr of local map {
            gettoken f t : pr, parse("=")
            gettoken eqsign t : t, parse("=")
            if "`f'" == "`eq'" local to `t'
        }
        if "`to'" == "" continue
        local ja = colnumb(`A', "`to':`nm'")
        if missing(`ja') {
            di as err "match_b: `to':`nm' (for `fn') is not a column of `A'"
            exit 9
        }
        if reldif(`A'[1, `ja'], `B'[1, `j']) > `tol' | missing(reldif(`A'[1, `ja'], `B'[1, `j'])) {
            di as err "match_b: `to':`nm' " %18.0g `A'[1, `ja'] " vs `fn' " %18.0g `B'[1, `j']
            exit 9
        }
        local ++n
    }
    if `n' == 0 {
        di as err "match_b: no column compared"
        exit 9
    }
    return scalar n = `n'
end

capture program drop log_count
program define log_count, rclass
    // log_count "file" "text"  --  r(n) = number of lines of the file that contain text
    args fn txt exact
    local n 0
    if "`exact'" == "exact" {
        mata: st_local("n", strofreal(sum(strrtrim(cat(st_local("fn"))) :== st_local("txt"))))
    }
    else {
        mata: st_local("n", strofreal(sum(strpos(cat(st_local("fn")), st_local("txt")) :> 0)))
    }
    return scalar n = `n'
end

capture program drop log_width
program define log_width, rclass
    // log_width "file"  --  r(max) = widest line between the markers SDW_W_BEGIN and
    // SDW_W_END, not counting command echoes (". " and "> ") and the marker commands
    args fn
    local w 0
    mata: _sd_l = cat(st_local("fn"))
    mata: _sd_b = select(1::rows(_sd_l), _sd_l :== "SDW_W_BEGIN")
    mata: _sd_e = select(1::rows(_sd_l), _sd_l :== "SDW_W_END")
    mata: _sd_l = _sd_l[(_sd_b[1] + 1)::(_sd_e[1] - 1)]
    mata: _sd_l = select(_sd_l, substr(_sd_l, 1, 2) :!= ". " :& substr(_sd_l, 1, 2) :!= "> ")
    mata: st_local("w", strofreal(max(udstrlen(_sd_l))))
    mata: mata drop _sd_l _sd_b _sd_e
    return scalar max = `w'
end

capture program drop stage_rows
program define stage_rows, rclass
    // stage_rows "file"  --  the rows of the stage tables (lines "    # block ..."), joined
    // with "|" in r(rows); r(n) their number
    args fn
    mata: _sd_l = cat(st_local("fn"))
    mata: _sd_l = select(_sd_l, regexm(_sd_l, "^    [0-9] [a-z] "))
    mata: st_local("n", strofreal(rows(_sd_l)))
    mata: st_local("rows", invtokens(_sd_l', "|"))
    mata: mata drop _sd_l
    return scalar n = `n'
    return local rows `"`rows'"'
end

capture program drop drop_stagemats
program define drop_stagemats, eclass
    // leave out the stage matrices, as a prefix that does not carry them would
    foreach m in bystage bystage2 bystage_rev bystage_rev2 {
        mata: st_matrix("e(`m')", J(0, 0, .))
    }
end

capture program drop drop_bsdraws
program define drop_bsdraws, eclass
    // leave out the replicates, as a run with more of them than c(max_matdim) does
    foreach m in bsdraws bsdraws_pct {
        mata: st_matrix("e(`m')", J(0, 0, .))
    }
end

capture program drop rbar_has
program define rbar_has, rclass
    // rbar_has graphname matname  --  r(found) = 1 if an rbar plot of the graph runs from 0 to
    // the values of the column vector matname, stage by stage (its serset, sorted by x)
    args G M
    local found 0
    preserve
    forvalues k = 1/30 {
        if "`.`G'.plotregion1.plot`k'.classname'" == "" continue, break
        if "`.`G'.plotregion1.plot`k'.type.snm'" != "rbar" continue
        local P .`G'.plotregion1.plot`k'
        qui serset set ``P'.serset.id'
        qui serset use, clear
        mata: for (i = 1; i <= st_nvar(); i++) st_varrename(i, "v" + strofreal(i))
        local vy2 v``P'.y2var'
        local vx v``P'.xvar'
        if _N != rowsof(`M') continue
        sort `vx'
        local ok 1
        forvalues s = 1/`=rowsof(`M')' {
            if reldif(`vy2'[`s'], `M'[`s', 1]) > 1e-10 local ok 0
        }
        if `ok' local found 1
    }
    restore
    return scalar found = `found'
end

tempfile lg bsf1 bsf2
local lgf `"`lg'.log"'
capture log close flog

* ---- data ---------------------------------------------------------------------
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
recast double y x3 wt
assert _N == 2000
gen double d  = 1 + exp(x3 / 3) + x2          // a positive denominator
gen double y3 = y^2 / 10                       // a third outcome
save `"`bsf1'.dta"', replace
local dta `"`bsf1'.dta"'
tempname B B1 B2 S S1 S2 V V1 V2 SE SE1 SE2


* ==============================================================================
* 1. ratio(denvar)
* ==============================================================================
di as txt _n "---- 1. ratio()"

* the three means are ratios of weighted means, the raked one under the generate() weights
seldecomp y x1 x2 [pw=wt], by(group) ratio(d) generate(rw)
assert "`e(ratio)'" == "d"
foreach g in 0 1 {
    qui sum y [aw=wt] if group == `g'
    local ny`g' = r(mean)
    qui sum d [aw=wt] if group == `g'
    local nd`g' = r(mean)
}
qui sum y [aw=rw] if group == 0
local nyr = r(mean)
qui sum d [aw=rw] if group == 0
local ndr = r(mean)
assert_close e(mean_focal) `ny0'/`nd0' `tol' "ratio: focal_group"
assert_close e(mean_ref) `ny1'/`nd1' `tol' "ratio: ref_group"
assert_close e(mean_raked_focal) `nyr'/`ndr' `tol' "ratio: raked_focal"
assert_close e(difference) e(mean_focal)-e(mean_ref) `same'
assert_close e(explained) e(mean_focal)-e(mean_raked_focal) `same'
assert_close _b[overall:explained] e(explained) `same'
assert_close _b[overall:raked_focal] e(mean_raked_focal) `same'
* e(ratio_parts): numerator and denominator means of each row
matrix `S' = e(ratio_parts)
assert rowsof(`S') == 3 & colsof(`S') == 2
assert "`: rowfullnames `S''" == "y:focal_group y:raked_focal y:ref_group"
assert "`: colnames `S''" == "numerator denominator"
assert_close `S'[1,1] `ny0' `tol'
assert_close `S'[1,2] `nd0' `tol'
assert_close `S'[2,1] `nyr' `tol'
assert_close `S'[2,2] `ndr' `tol'
assert_close `S'[3,1] `ny1' `tol'
assert_close `S'[3,2] `nd1' `tol'
forvalues r = 1/3 {
    assert_close `S'[`r',1]/`S'[`r',2] el(e(b),1,`r') `same' "ratio_parts row `r'"
}
* the raked weights do not depend on the outcome: the same as without ratio()
seldecomp y x1 x2 [pw=wt], by(group) generate(rw0)
assert rw == rw0
assert !missing(e(mean_focal)) & "`e(ratio)'" == ""
capture confirm matrix e(ratio_parts)
assert _rc
drop rw rw0

* a constant denominator 2 halves every element of e(b) (sequentially fixothermargins)
gen double two = 2
seldecomp y (a: x1) (b: x2 x4) [pw=wt], by(group) sequentially fixothermargins
matrix `B1' = e(b)
matrix `S1' = e(bystage)
seldecomp y (a: x1) (b: x2 x4) [pw=wt], by(group) sequentially fixothermargins ratio(two)
matrix `B2' = e(b)
matrix `S2' = e(bystage)
matrix `B1' = `B1' / 2
mat_close `B1' `B2' `same'
assert_close e(pct_explained) `S1'[2,4] `same'
forvalues s = 1/2 {
    assert_close `S2'[`s',4] `S1'[`s',4] `same' "ratio(two): pct_explained, stage `s'"
}

* the stages: the raked ratio of every stage from its generate() weights
seldecomp y (a: x1) (b: x2 x4) [pw=wt], by(group) sequentially fixothermargins ratio(d) generate(rs)
matrix `B1' = e(b)
matrix `S' = e(bystage)
foreach s in 1 2 {
    local v = cond(`s' == 1, "rs_s1", "rs")
    qui sum y [aw=`v'] if group == 0
    local a = r(mean)
    qui sum d [aw=`v'] if group == 0
    assert_close `S'[`s',1] `a'/r(mean) `tol' "ratio: stage `s' raked"
}
assert_close `S'[1,5]+`S'[2,5] e(explained) `same'
drop rs rs_s1

* the change-of-support split holds for ratios: a level found only in the Focal Group
gen byte x1z = x1
replace x1z = 9 if group == 0 & id <= 60
seldecomp y x1z x2 [pw=wt], by(group) ratio(d)
assert e(N_focal_zero) > 0
assert_close e(explained_support)+e(explained_cs) e(explained) `same'
qui sum y [aw=wt] if group == 0 & x1z == 9
local a = r(mean)
qui sum d [aw=wt] if group == 0 & x1z == 9
assert_close e(mean_focal_zero) `a'/r(mean) `tol' "ratio: mean_focal_zero"
qui sum y [aw=wt] if group == 0 & x1z != 9
local a = r(mean)
qui sum d [aw=wt] if group == 0 & x1z != 9
assert_close e(mean_focal_cs) `a'/r(mean) `tol' "ratio: mean_focal_cs"

* contvars() cut points are unaffected (they do not depend on the outcome)
seldecomp y c.x3 x2, by(group) contvars(x3 4)
matrix `S1' = e(cuts1)
seldecomp y c.x3 x2, by(group) contvars(x3 4) ratio(d)
matrix `S2' = e(cuts1)
mat_close `S1' `S2' 0

* bayesboot(): every replicate recomputes the ratio (saving() rows = plain runs on the
* generate() weights)
seldecomp y (a: x1) (b: x2) [pw=wt], by(group) sequentially fixothermargins ratio(d) ///
    bayesboot(reps(3) rseed(5) generate(bw) saving(`"`bsf2'"', replace) nodots)
local bsn `e(bsnames)'
assert wordcount("`bsn'") == 6 + 4 * 2
preserve
use `"`bsf2'"', clear
mkmat _all, matrix(`V')
restore
forvalues r = 1/3 {
    qui seldecomp y (a: x1) (b: x2) [iw=bw_`r'], by(group) sequentially fixothermargins ratio(d)
    matrix `B' = e(b)
    local j 0
    foreach fn of local bsn {
        local ++j
        local jb = colnumb(`B', "`fn'")
        assert_close `V'[`r',`j'] `B'[1,`jb'] 1e-10 "bayesboot ratio: replicate `r', `fn'"
    }
    assert_close `V'[`r',`j'+1] e(pct_explained) 1e-10
}
drop bw_*

* the bootstrap prefix recomputes the ratio in every replicate
bootstrap, reps(4) seed(12) nodots: seldecomp y x1 x2, by(group) ratio(d)
matrix `B2' = e(b)
matrix `V' = e(V)
assert `V'[5, 5] > 0
qui seldecomp y x1 x2, by(group) ratio(d)
matrix `B1' = e(b)
mat_close `B1' `B2' `same'

* errors
rcof "noi seldecomp y x1 x2, by(group) ratio(x4)" == 108              // string
gen double dmiss = .
rcof "noi seldecomp y x1 x2, by(group) ratio(dmiss)" == 2000          // all missing
gen double dzero = cond(group == 0, 0, 1)
rcof "noi seldecomp y x1 x2, by(group) ratio(dzero)" == 459           // zero mean in the Focal Group
rcof "noi seldecomp y x1 x2, by(group) ratio(nosuchvar)" == 111
* a denominator missing for some observations marks them out of the sample
gen double dm = d
replace dm = . in 1/10
seldecomp y x1 x2, by(group) ratio(dm)
assert e(N) == 1990
drop dm dmiss dzero


* ==============================================================================
* 2. outcomes(varlist)
* ==============================================================================
di as txt _n "---- 2. outcomes()"

* every outcome's rows equal its single-outcome run (sequentially fixothermargins, weights,
* a zeroed level)
local rspec (a: x1z) (b: x2 x4)
seldecomp y `rspec' [pw=wt], by(group) sequentially fixothermargins outcomes(ybin y3) generate(rwm) covariatedetail
matrix `B' = e(b)
assert "`e(depvar)'" == "y ybin y3"
assert "`e(outcomes)'" == "y ybin y3"
assert e(N_outcomes) == 3
assert e(k_eq) == 15
assert colsof(`B') == 3 * (6 + 4 * 2)
local cn : colfullnames `B'
assert "`: word 1 of `cn''" == "y:focal_group"
assert colnumb(`B', "ybin_marg_expl:a") < . & colnumb(`B', "y3_cum_unexpl:b") < .
assert colnumb(`B', "overall:focal_group") == .
matrix `S' = e(by_outcome)
matrix `V' = e(covariates)
matrix `SE' = e(bystage)
matrix `S2' = e(bystage3)
local j 0
foreach yv in y ybin y3 {
    local ++j
    qui seldecomp `yv' `rspec' [pw=wt], by(group) sequentially fixothermargins generate(rws) covariatedetail
    matrix `B1' = e(b)
    match_b `B' `B1', map(overall=`yv' marg_expl=`yv'_marg_expl cum_expl=`yv'_cum_expl ///
        cum_unexpl=`yv'_cum_unexpl raked_focal=`yv'_raked_focal) tol(`same')
    assert r(n) == 14
    local c 0
    foreach q in mean_focal mean_raked_focal mean_ref difference explained unexplained pct_explained {
        assert_close `S'[`j',`++c'] e(`q') `same' "by_outcome[`yv', `q']"
    }
    assert missing(`S'[`j', 8])
    matrix `S1' = e(bystage)
    if `j' == 1 mat_close `SE' `S1' `same'
    if `j' == 3 mat_close `S2' `S1' `same'
    * the raked weights and the covariate shares do not depend on the outcome
    assert reldif(rwm, rws) == 0 if !missing(rwm)
    matrix `S1' = e(covariates)
    mat_close `V' `S1' 0
    drop rws rws_s1
}
drop rwm rwm_s1
* the scalars are the first outcome's
qui seldecomp y `rspec' [pw=wt], by(group) sequentially fixothermargins outcomes(ybin y3)
assert_close e(explained_support)+e(explained_cs) e(explained) `same'
local es1 = e(explained_support)
qui seldecomp y `rspec' [pw=wt], by(group) sequentially fixothermargins
assert_close e(explained_support) `es1' `same'

* ratio() with outcomes(): every outcome over the same denominator
seldecomp y x1 x2 [pw=wt], by(group) outcomes(ybin) ratio(d)
matrix `B' = e(b)
matrix `S' = e(ratio_parts)
assert rowsof(`S') == 6
assert "`: rowfullnames `S''" == "y:focal_group y:raked_focal y:ref_group ybin:focal_group ybin:raked_focal ybin:ref_group"
matrix `S1' = `S'[1..3, 2]
matrix `S2' = `S'[4..6, 2]
mat_close `S1' `S2' `same'
foreach yv in y ybin {
    qui seldecomp `yv' x1 x2 [pw=wt], by(group) ratio(d)
    matrix `B1' = e(b)
    match_b `B' `B1', map(overall=`yv') tol(`same')
}

* joint inference: one set of draws per replicate; with the same seed each outcome's
* standard errors equal those of its single-outcome run
seldecomp y (a: x1) (b: x2) [pw=wt], by(group) sequentially fixothermargins outcomes(ybin) ///
    bayesboot(reps(10) rseed(7) nodots)
matrix `SE' = e(se)
matrix `V' = e(V)
matrix `S' = e(by_outcome)
assert wordcount("`e(bsnames)'") == 2 * 14
assert colsof(`V') == 28
* the covariance between the two outcomes' explained parts is estimated (not zero)
local p1 = colnumb(`V', "y:explained")
local p2 = colnumb(`V', "ybin:explained")
assert `V'[`p1', `p2'] != 0
local j 0
foreach yv in y ybin {
    local ++j
    qui seldecomp `yv' (a: x1) (b: x2) [pw=wt], by(group) sequentially fixothermargins ///
        bayesboot(reps(10) rseed(7) nodots)
    matrix `SE1' = e(se)
    match_b `SE' `SE1', map(overall=`yv' marg_expl=`yv'_marg_expl cum_expl=`yv'_cum_expl ///
        cum_unexpl=`yv'_cum_unexpl raked_focal=`yv'_raked_focal) tol(`same')
    assert r(n) == 14
    assert_close `S'[`j',8] e(pct_explained_se) `same' "pct_explained_se of `yv'"
    matrix `V1' = e(pct_explained_ci)
    assert_close `S'[`j',9] `V1'[1,1] `same'
    assert_close `S'[`j',10] `V1'[1,2] `same'
}

* bca and saving(): the jackknife acceleration of each outcome equals its single run's;
* bstat reads the file
preserve
keep if cluster <= 20
seldecomp y x1 x2 [pw=wt], by(group) outcomes(ybin) ///
    bayesboot(reps(20) rseed(3) bca cluster(cluster) saving(`"`bsf2'"', replace) nodots)
matrix `SE' = e(se)
matrix `S' = e(accel)
assert colsof(`S') == 12
qui seldecomp ybin x1 x2 [pw=wt], by(group) bayesboot(reps(20) rseed(3) bca cluster(cluster) nodots)
matrix `S1' = e(accel)
match_b `S' `S1', map(overall=ybin) tol(`same')
bstat using `"`bsf2'"'
matrix `SE1' = e(se)
assert colsof(`SE1') == 14                 // 12 statistics and the two explained shares
match_b `SE1' `SE', map(y=y ybin=ybin) tol(1e-10)
restore

* a prefix replicates the whole e(b): svy jackknife with outcomes(), reverse and ratio();
* replay and makeplot afterwards
preserve
keep if cluster <= 15
svyset cluster [pw=wt]
svy jackknife, nodots: seldecomp y (a: x1) (b: x2) , by(group) sequentially fixothermargins outcomes(ybin) ///
    reverse ratio(d)
matrix `B2' = e(b)
matrix `V' = e(V)
assert colsof(`B2') == 2 * (6 + 8 * 2)
forvalues j = 1/`=colsof(`B2')' {
    assert `V'[`j', `j'] > 0
}
qui seldecomp y (a: x1) (b: x2) [pw=wt], by(group) sequentially fixothermargins outcomes(ybin) reverse ratio(d)
matrix `B1' = e(b)
mat_close `B1' `B2' `same'
svy jackknife, nodots: seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins outcomes(ybin) reverse ratio(d)
log using `"`lgf'"', text replace name(flog)
seldecomp
seldecomp, makeplot(stages, outcome(ybin) twoway(name(fsvy, replace)))
log close flog
log_count `"`lgf'"' "Outcomes (one raking): y ybin"
assert r(n) == 1
* the header's line on the reverse order, and each outcome's reverse-order table (its title,
* longer than 78 columns with the outcome's name, wraps: its first part is looked for)
log_count `"`lgf'"' "and in reverse order: + b + a"
assert r(n) == 1
log_count `"`lgf'"' "Decomposition by raking stage (sequentially, fixothermargins), reverse"
assert r(n) == 2
matrix `S' = J(2, 1, .)
matrix `S'[1, 1] = _b[ybin_cum_expl:a]
matrix `S'[2, 1] = _b[ybin_cum_expl:b]
rbar_has fsvy `S'
assert r(found) == 1
graph drop fsvy
svyset, clear
restore

* the bootstrap prefix
bootstrap, reps(3) seed(2) nodots: seldecomp y x1 x2, by(group) outcomes(ybin)
matrix `B2' = e(b)
assert "`: colfullnames `B2''" == "y:focal_group y:raked_focal y:ref_group y:difference y:explained y:unexplained ybin:focal_group ybin:raked_focal ybin:ref_group ybin:difference ybin:explained ybin:unexplained"

* the by prefix: each group's results are those of the if-restricted run
bysort x2: seldecomp y x1 [pw=wt], by(group) outcomes(ybin y3) ratio(d)
matrix `B2' = e(b)
qui seldecomp y x1 [pw=wt] if x2 == 1, by(group) outcomes(ybin y3) ratio(d)
matrix `B1' = e(b)
mat_close `B1' `B2' `same'

* estimates store and restore, then replay
seldecomp y (a: x1) (b: x2) [pw=wt], by(group) sequentially fixothermargins outcomes(ybin)
matrix `B1' = e(b)
estimates store feat
seldecomp y x1, by(group)
estimates restore feat
matrix `B2' = e(b)
mat_close `B1' `B2' 0
log using `"`lgf'"', text replace name(flog)
seldecomp
log close flog
log_count `"`lgf'"' "Decomposition by raking stage (sequentially, fixothermargins): "
assert r(n) == 2
estimates drop feat

* makeplot draws the first outcome and says so; makeplot(, outcome()) picks another
seldecomp y (a: x1) (b: x2) [pw=wt], by(group) sequentially fixothermargins outcomes(ybin)
log using `"`lgf'"', text replace name(flog)
seldecomp, makeplot(stages, twoway(name(fo1, replace)))
log close flog
log_count `"`lgf'"' "(makeplot: outcome y, the first of the outcomes"
assert r(n) == 1
matrix `S' = (_b[y_cum_expl:a] \ _b[y_cum_expl:b])
rbar_has fo1 `S'
assert r(found) == 1
seldecomp, makeplot(stages, outcome(ybin) twoway(name(fo2, replace)))
matrix `S' = (_b[ybin_cum_unexpl:a] \ _b[ybin_cum_unexpl:b])
rbar_has fo2 `S'
assert r(found) == 1
rcof "noi seldecomp, makeplot(, outcome(y3))" == 198
graph drop fo1 fo2

* errors
rcof "noi seldecomp y x1 x2, by(group) outcomes(y)" == 198             // the depvar again
rcof "noi seldecomp y x1 x2, by(group) outcomes(ybin ybin)" == 198     // twice
rcof "noi seldecomp y x1 x2, by(group) outcomes(x4)" == 108            // string
rcof "noi seldecomp y x1 x2, by(group) outcomes(x2)" == 198            // a raking variable
rcof "noi seldecomp y x1 x2, by(group) outcomes(L.ybin)" == 101
gen double y_marg_expl = ybin
rcof "noi seldecomp y (a: x1) (b: x2), by(group) sequentially fixothermargins outcomes(y_marg_expl)" == 198
seldecomp y x1 x2, by(group) outcomes(y_marg_expl)                       // no clash without sequentially
drop y_marg_expl
* an outcome missing for some observations marks them out of the sample
gen double ym = ybin
replace ym = . in 1/7
seldecomp y x1 x2, by(group) outcomes(ym)
assert e(N) == 1993
drop ym


* ==============================================================================
* 3. reverse
* ==============================================================================
di as txt _n "---- 3. reverse"

* the mode options: reverse requires sequentially and is refused with individually;
* fixothermargins requires sequentially or individually; the two modes exclude each other;
* the words varbyvar, cumulative, holdothers, isolated, nohold and refminusfocal are refused
log using `"`lgf'"', text replace name(flog)
rcof "noi seldecomp y x1 x2, by(group) reverse" == 198
rcof "noi seldecomp y (a: x1) (b: x2), by(group) individually reverse" == 198
rcof "noi seldecomp y (a: x1) (b: x2), by(group) individually fixothermargins reverse" == 198
rcof "noi seldecomp y (a: x1) (b: x2), by(group) fixothermargins" == 198
rcof "noi seldecomp y (a: x1) (b: x2), by(group) fixothermargins reverse" == 198
rcof "noi seldecomp y (a: x1) (b: x2), by(group) sequentially individually" == 198
rcof "noi seldecomp y (a: x1) (b: x2), by(group) sequentially individually reverse" == 198
local oldw varbyvar cumulative holdothers isolated nohold refminusfocal
foreach w of local oldw {
    rcof "noi seldecomp y (a: x1) (b: x2), by(group) sequentially `w'" == 198
}
log close flog
log_count `"`lgf'"' "reverse requires sequentially"
assert r(n) == 1
log_count `"`lgf'"' "reverse may not be combined with individually"
assert r(n) == 2
log_count `"`lgf'"' "fixothermargins requires sequentially or individually"
assert r(n) == 2
log_count `"`lgf'"' "sequentially and individually may not be combined"
assert r(n) == 2
foreach w of local oldw {
    log_count `"`lgf'"' "option `w' not allowed"
    assert r(n) == 1
}

* the reverse order equals a run with the blocks typed in reverse order, to the last digit;
* with the blocks not yet added held at the Focal Group's margins (fixothermargins) and free
foreach var in fixothermargins "" {
    seldecomp y (a: x1) (b: x2) (c: x4) [pw=wt], by(group) sequentially `var' reverse
    assert "`e(reverse)'" == "reverse"
    assert e(k_eq) == 9
    matrix `B' = e(b)
    matrix `S' = e(bystage)
    matrix `S1' = e(bystage_rev)
    assert "`: rownames `S1''" == "c b a"
    assert colsof(`B') == 6 + 8 * 3
    * the forward part is the run without reverse
    qui seldecomp y (a: x1) (b: x2) (c: x4) [pw=wt], by(group) sequentially `var'
    matrix `B1' = e(b)
    match_b `B' `B1', map(overall=overall marg_expl=marg_expl cum_expl=cum_expl ///
        cum_unexpl=cum_unexpl raked_focal=raked_focal) tol(`same')
    assert r(n) == 18
    matrix `S2' = e(bystage)
    mat_close `S' `S2' `same'
    * the reverse part is the retyped run
    qui seldecomp y (c: x4) (b: x2) (a: x1) [pw=wt], by(group) sequentially `var'
    matrix `B2' = e(b)
    matrix `S2' = e(bystage)
    mat_close `S1' `S2' `same'
    match_b `B' `B2', map(marg_expl=rev_marg_expl cum_expl=rev_cum_expl ///
        cum_unexpl=rev_cum_unexpl raked_focal=rev_raked_focal) tol(`same')
    assert r(n) == 12
    * the final stages rake to the same targets; the increments sum to the same explained part
    assert_close `S'[3,1] `S1'[3,1] `rtol' "final stages of the two orders (`var')"
    assert_close `S'[1,5]+`S'[2,5]+`S'[3,5] `S1'[1,5]+`S1'[2,5]+`S1'[3,5] `rtol'
    assert_close `S1'[1,5]+`S1'[2,5]+`S1'[3,5] _b[overall:explained] `rtol'
}

* bayesboot(): the reverse rows are bootstrapped on the same draws, so with the same
* seed they equal the retyped run's; the forward rows equal the run without reverse
seldecomp y (a: x1) (b: x2) (c: x4) [pw=wt], by(group) sequentially fixothermargins reverse ///
    bayesboot(reps(8) rseed(11) nodots)
matrix `SE' = e(se)
matrix `V' = e(V)
assert wordcount("`e(bsnames)'") == 30
assert colsof(`SE') == 30
qui seldecomp y (c: x4) (b: x2) (a: x1) [pw=wt], by(group) sequentially fixothermargins ///
    bayesboot(reps(8) rseed(11) nodots)
matrix `SE2' = e(se)
match_b `SE' `SE2', map(marg_expl=rev_marg_expl cum_expl=rev_cum_expl cum_unexpl=rev_cum_unexpl ///
    raked_focal=rev_raked_focal) tol(`same')
qui seldecomp y (a: x1) (b: x2) (c: x4) [pw=wt], by(group) sequentially fixothermargins ///
    bayesboot(reps(8) rseed(11) nodots)
matrix `SE1' = e(se)
match_b `SE' `SE1', map(overall=overall marg_expl=marg_expl cum_expl=cum_expl cum_unexpl=cum_unexpl ///
    raked_focal=raked_focal) tol(`same')

* with outcomes(): each outcome's reverse rows equal its retyped single run
seldecomp y (a: x1) (b: x2 x4) [pw=wt], by(group) sequentially fixothermargins outcomes(ybin) reverse
matrix `B' = e(b)
assert e(k_eq) == 18
matrix `S1' = e(bystage_rev2)
qui seldecomp ybin (b: x2 x4) (a: x1) [pw=wt], by(group) sequentially fixothermargins
matrix `B2' = e(b)
matrix `S2' = e(bystage)
mat_close `S1' `S2' `same'
match_b `B' `B2', map(marg_expl=ybin_rev_marg_expl cum_expl=ybin_rev_cum_expl ///
    cum_unexpl=ybin_rev_cum_unexpl raked_focal=ybin_rev_raked_focal) tol(`same')

* the stage tables are rebuilt from e(b) when the stage matrices are missing (as a prefix
* could leave them out): the same rows as from the matrices
seldecomp y (a: x1) (b: x2) [pw=wt], by(group) sequentially fixothermargins reverse outcomes(ybin)
log using `"`lgf'"', text replace name(flog)
seldecomp
log close flog
stage_rows `"`lgf'"'
local rows0 `"`r(rows)'"'
assert r(n) == 8
drop_stagemats
log using `"`lgf'"', text replace name(flog)
seldecomp
log close flog
stage_rows `"`lgf'"'
assert r(n) == 8
assert `"`r(rows)'"' == `"`rows0'"'

* makeplot draws the forward order; makeplot(, reverse) the reverse order's stages (c, b, a)
seldecomp y (a: x1) (b: x2) (c: x4) [pw=wt], by(group) sequentially fixothermargins reverse makeplot(stages, twoway(name(frev, replace)))
matrix `S' = (_b[cum_expl:a] \ _b[cum_expl:b] \ _b[cum_expl:c])
rbar_has frev `S'
assert r(found) == 1
seldecomp, makeplot(stages, reverse twoway(name(frev2, replace)))
matrix `S' = (_b[rev_cum_expl:c] \ _b[rev_cum_expl:b] \ _b[rev_cum_expl:a])
rbar_has frev2 `S'
assert r(found) == 1
* (the marginal plot: the explained part and the stages are one plot, from the top)
seldecomp, makeplot(marginal, reverse twoway(name(frev3, replace)))
matrix `S' = (_b[overall:explained] \ _b[rev_marg_expl:c] \ _b[rev_marg_expl:b] \ _b[rev_marg_expl:a])
rbar_has frev3 `S'
assert r(found) == 1
graph drop frev frev2 frev3


* ==============================================================================
* 4. display
* ==============================================================================
di as txt _n "---- 4. display"
* the output of the three options together, at linesize 160: one coefficient table, stage
* table and reverse stage table per outcome, the header lines, at most 78 columns; the
* bootstrap table only with citable (no pointer line to it); with
* inference the stage rows of every stage equation, forward and reverse
local ls0 = c(linesize)
local fcmd seldecomp y (a: x1) (b: x2) [pw=wt], by(group) sequentially fixothermargins reverse outcomes(ybin) ratio(d)
local fbs bayesboot(reps(10) rseed(2) nodots)
set linesize 160
log using `"`lgf'"', text replace name(flog)
di "SDW_W_BEGIN"
`fcmd' `fbs'
di "SDW_W_END"
log close flog
set linesize `ls0'
log_width `"`lgf'"'
di as txt "widest display line: " r(max)
assert r(max) <= 78
log_count `"`lgf'"' "Outcomes (one raking): y ybin"
assert r(n) == 1
log_count `"`lgf'"' "Ratio of means: each outcome / d"
assert r(n) == 1
log_count `"`lgf'"' "and in reverse order: + b + a"
assert r(n) == 1
* (a line naming ybin also contains the text naming y; a reverse-order title with an
* outcome is longer than 78 columns and wraps: its first part is looked for)
log_count `"`lgf'"' "Decomposition by raking stage (sequentially, fixothermargins): y"
assert r(n) == 2
log_count `"`lgf'"' "Decomposition by raking stage (sequentially, fixothermargins): ybin"
assert r(n) == 1
log_count `"`lgf'"' "Decomposition by raking stage (sequentially, fixothermargins), reverse"
assert r(n) == 2
* one coefficient table per outcome, headed by the outcome; no bootstrap table
log_count `"`lgf'"' "           y |"
assert r(n) == 1
log_count `"`lgf'"' "        ybin |"
assert r(n) == 1
log_count `"`lgf'"' "Explained share of the difference:"
assert r(n) == 2
log_count `"`lgf'"' "Bayesian bootstrap confidence intervals"
assert r(n) == 0
* no pointer line to citable
log_count `"`lgf'"' "intervals: seldecomp, citable"
assert r(n) == 0
local strows `""By stage         | Marginal explained" "By stage         | Cumulative explained""'
local strows `"`strows' "By stage         | Cumulative unexplained" "By stage         | Raked focal group ratio""'
local strows `"`strows' "By stage         | Marginal explained, reverse order" "By stage         | Cumulative explained, reverse order""'
local strows `"`strows' "By stage         | Cumulative unexplained, reverse order" "By stage         | Raked focal group ratio, reverse order""'
foreach t of local strows {
    log_count `"`lgf'"' "`t'" exact
    assert r(n) == 2
}
* citable, at estimation and at replay: the bootstrap table of each outcome, the stage rows
* of both orders in both tables, no pointer line
foreach how in estimation replay {
    set linesize 160
    log using `"`lgf'"', text replace name(flog)
    di "SDW_W_BEGIN"
    if "`how'" == "estimation" `fcmd' bayesboot(reps(10) rseed(2) nodots citable)
    else seldecomp, citable
    di "SDW_W_END"
    log close flog
    set linesize `ls0'
    log_width `"`lgf'"'
    assert r(max) <= 78
    log_count `"`lgf'"' "           y |"
    assert r(n) == 2
    log_count `"`lgf'"' "        ybin |"
    assert r(n) == 2
    log_count `"`lgf'"' "Explained share of the difference:"
    assert r(n) == 2
    log_count `"`lgf'"' "Bayesian bootstrap confidence intervals"
    assert r(n) == 2
    log_count `"`lgf'"' "Percentile and BC intervals: seldecomp, citable"
    assert r(n) == 0
    foreach t of local strows {
        log_count `"`lgf'"' "`t'" exact
        assert r(n) == 4
    }
}
* one outcome without inference: the stage tables hold the stage values, and the coefficient
* table has the overall equation alone
qui seldecomp y (a: x1) (b: x2) [pw=wt], by(group) sequentially fixothermargins reverse
log using `"`lgf'"', text replace name(flog)
seldecomp
log close flog
log_count `"`lgf'"' "Overall          |"
assert r(n) == 1
foreach t of local strows {
    log_count `"`lgf'"' "`t'" exact
    assert r(n) == 0
}

* a replay at another level with outcomes() recomputes every interval, the explained
* shares' too (one row per outcome), exactly as a fresh estimation at that level with the
* same seed; e() keeps the estimation level's
preserve
keep if cluster <= 20
local fbs bayesboot(reps(20) rseed(4) bca cluster(cluster) nodots)
qui `fcmd' `fbs'
matrix `S' = e(ci_bca)
matrix `S1' = e(by_outcome)
log using `"`lgf'"', text replace name(flog)
seldecomp, level(90) citable
assert r(level) == 90 & e(level) == 95
foreach m in ci_normal ci_percentile ci_bc ci_bca pct_explained_ci {
    tempname R_`m'
    matrix `R_`m'' = r(`m')
}
log close flog
assert "`: rownames `R_pct_explained_ci''" == "y ybin"
assert "`: colnames `R_pct_explained_ci''" == "ll ul"
matrix `S2' = e(ci_bca)
mat_close `S' `S2' 0
matrix `S2' = e(by_outcome)
mat_close `S1' `S2' 0
* the displayed share intervals are those of r(), each on its outcome's table
forvalues k = 1/2 {
    local pl = strtrim(string(`R_pct_explained_ci'[`k', 1], "%9.1f"))
    local pu = strtrim(string(`R_pct_explained_ci'[`k', 2], "%9.1f"))
    log_count `"`lgf'"' "    90% percentile CI: `pl'% to `pu'%"
    assert r(n) == 1
}
log_count `"`lgf'"' "[90% conf. interval]"
assert r(n) == 4
qui `fcmd' bayesboot(reps(20) rseed(4) bca cluster(cluster) nodots level(90))
foreach m in ci_normal ci_percentile ci_bc ci_bca {
    matrix `S2' = e(`m')
    mat_close `R_`m'' `S2' `same'
}
matrix `S2' = e(by_outcome)
matrix `S2' = `S2'[1..., 9..10]
mat_close `R_pct_explained_ci' `S2' `same'
matrix `S2' = e(pct_explained_ci)
matrix `S1' = `R_pct_explained_ci'[1, 1..2]
mat_close `S1' `S2' `same'
* without the replicates in e() (as with more of them than c(max_matdim)), a replay at
* another level shows the estimation level's intervals in the bootstrap tables and says so
* once, under the first outcome's table
drop_bsdraws
log using `"`lgf'"', text replace name(flog)
seldecomp, level(80) citable
log close flog
log_count `"`lgf'"' "Bayesian bootstrap confidence intervals"
assert r(n) == 2
* BCa is the first requested type, so both main and citable intervals retain 90%.
log_count `"`lgf'"' "[90% conf. interval]"
assert r(n) == 4
log_count `"`lgf'"' "(intervals at the estimation level, 90%: level() cannot change them"
assert r(n) == 1
restore

* makeplot(marginal): the first outcome with its message, outcome() another; the value
* axis of a ratio says so
seldecomp y (a: x1) (b: x2) [pw=wt], by(group) sequentially fixothermargins outcomes(ybin) ratio(d) ///
    bayesboot(reps(10) rseed(2) nodots)
log using `"`lgf'"', text replace name(flog)
seldecomp, makeplot(marginal, twoway(name(fc1, replace)))
log close flog
log_count `"`lgf'"' "(makeplot: outcome y, the first of the outcomes"
assert r(n) == 1
* (the explained part and the stages are one plot, from the top)
matrix `S' = (_b[y:explained] \ _b[y_marg_expl:a] \ _b[y_marg_expl:b])
rbar_has fc1 `S'
assert r(found) == 1
seldecomp, makeplot(marginal, outcome(ybin) twoway(name(fc2, replace)))
matrix `S' = (_b[ybin:explained] \ _b[ybin_marg_expl:a] \ _b[ybin_marg_expl:b])
rbar_has fc2 `S'
assert r(found) == 1
foreach q in difference unexplained {
    matrix `S' = (_b[ybin:`q'])
    rbar_has fc2 `S'
    assert r(found) == 1
}
* (the value axis reads "<label>: difference in ratios")
assert `"`.fc2.xaxis1.title.text[1]'"' == "ybin / d: difference in ratios"
seldecomp, makeplot(, outcome(ybin) twoway(name(fc3, replace)))
assert `"`.fc3.yaxis1.title.text[1]'"' == "ybin / d: difference in ratios"
rcof "noi seldecomp, makeplot(marginal, outcome(y3))" == 198
graph drop fc1 fc2 fc3

* a single-outcome run without the options shows none of it
log using `"`lgf'"', text replace name(flog)
seldecomp y (a: x1) (b: x2) [pw=wt], by(group) sequentially fixothermargins
log close flog
foreach t in "Outcomes (one raking)" "Ratio of means" "reverse order" "rev_" {
    log_count `"`lgf'"' "`t'"
    assert r(n) == 0
}

capture erase `"`lgf'"'
capture erase `"`bsf2'.dta"'
capture erase `"`bsf2'"'
di as res _n "PASSED: test_features.do"

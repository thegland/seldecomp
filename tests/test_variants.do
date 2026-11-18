* test_variants.do -- the staged mode individually, with and without fixothermargins, and
* the insupport screen: each against the runs it must
* reproduce. individually fixothermargins (the others held at Focal Group's margins):
* raking k equals stage 1 of a sequentially fixothermargins run with block k typed first
* and the other blocks after it (1e-12); individually alone (the other blocks free, the
* default): raking k equals the plain run on block k (1e-12); both: the full raking equals
* the one-stage run, the remainder identity, the order (exact without fixothermargins, to
* the raking tolerance with it), outcomes(), ratio(), generate(), the bootstrap
* (bayesboot(): SEs equal those of the corresponding runs with the same seed, saving()
* rows equal plain runs on the replicate weights), the prefixes, the display (the mode line
* names individually and whether the other blocks are fixed; no refused option word in the
* output), makeplot, the parsing rules (sequentially and individually exclude each other;
* reverse, detail and the block name _remainder are refused with individually; a variable in
* two blocks only with fixothermargins; the words isolated and nohold are refused) and the
* errors. insupport: the same as the run restricted
* with -if- to the screened sample, the same as the run without it when no level is missing
* from the Reference Group (every e(), elist-style), the held shares, bayesboot(), by:, svy jackknife,
* generate(), covariatedetail, outcomes(), ratio(), the cut points, the display and the
* errors.
*
* Run from the tests directory:   . do test_variants.do

version 16
clear
set more off
set varabbrev off
discard             // drop any seldecomp already in memory; the repository copy is what gets tested

* ---- tolerances ---------------------------------------------------------------
local tol   1e-6     // shares of the raked weights against their targets (ctrltolerance 1e-6)
local mtol  1e-10    // means from summarize under the generate() weights vs seldecomp
local same  1e-12    // the same computation reached by two routes
local rtol  1e-7     // the same targets reached along two raking paths (tolerance(1e-8))

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
    // assert_close exp1 exp2 [tol] ["label"]  --  asserts reldif(exp1,exp2) <= tol
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

capture program drop bcol
program define bcol
    // bcol lname matname "eq:name"  --  the element of row 1 of matname named eq:name (asserted
    // to exist) in the local lname
    args lname M nm
    local c = colnumb(`M', "`nm'")
    if missing(`c') {
        di as err "bcol: `nm' is not a column of `M'"
        exit 9
    }
    c_local `lname' = `M'[1, `c']
end

capture program drop log_count
program define log_count, rclass
    // log_count "file" "text"  --  r(n) = number of lines of the file that contain text
    args fn txt
    local n 0
    mata: st_local("n", strofreal(sum(strpos(cat(st_local("fn")), st_local("txt")) :> 0)))
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

capture program drop contrib_lines
program define contrib_lines, rclass
    // contrib_lines "file"  --  the lines of the contributions table (from its title to the
    // footnote), joined with "|" in r(lines)
    args fn
    mata: _sd_l = cat(st_local("fn"))
    mata: _sd_b = select(1::rows(_sd_l), strpos(_sd_l, "Contribution of each block, raked alone") :== 1)
    mata: _sd_e = select(1::rows(_sd_l), strpos(_sd_l, "%: percent of the difference. Other blocks") :== 1)
    mata: st_local("lines", invtokens(_sd_l[_sd_b[1]::_sd_e[1]]', "|"))
    mata: mata drop _sd_l _sd_b _sd_e
    return local lines `"`lines'"'
end

capture program drop drop_stagemats
program define drop_stagemats, eclass
    // leave out the stage matrices, as a prefix that does not carry them would
    foreach m in bystage bystage2 {
        mata: st_matrix("e(`m')", J(0, 0, .))
    }
end

capture program drop e_dump
program define e_dump
    // e_dump "file" [, skip(names)]  --  every e() scalar, macro and matrix, one line each,
    // leaving out e(cmdline) and the names in skip()
    syntax anything(name=fn) [, skip(string)]
    local skip `skip' cmdline
    tempname fh M
    file open `fh' using `fn', write text replace
    foreach t in scalars macros matrices {
        local names : e(`t')
        local names : list sort names
        local names : list names - skip
        foreach n of local names {
            if "`t'" == "scalars" file write `fh' "s `n' " %21.15g (e(`n')) _n
            else if "`t'" == "macros" {
                mata: st_local("v", st_global("e(`n')"))
                file write `fh' `"m `n' `macval(v)'"' _n
            }
            else {
                matrix `M' = e(`n')
                file write `fh' "M `n' `: rowfullnames `M'' | `: colfullnames `M''" _n
                forvalues i = 1/`=rowsof(`M')' {
                    forvalues j = 1/`=colsof(`M')' {
                        file write `fh' %21.15g (`M'[`i', `j']) " "
                    }
                    file write `fh' _n
                }
            }
        }
    }
    file close `fh'
end

capture program drop same_file
program define same_file
    // same_file "a" "b"  --  the two text files have the same lines
    args a b
    mata: _sd_a = cat(st_local("a")); _sd_b = cat(st_local("b"))
    mata: st_local("n", strofreal(rows(_sd_a) == rows(_sd_b) ? sum(_sd_a :!= _sd_b) : -1))
    if `n' != 0 {
        mata: st_local("n1", strofreal(rows(_sd_a))); st_local("n2", strofreal(rows(_sd_b)))
        di as err "same_file: the files differ (`n1' and `n2' lines)"
        mata: _sd_k = min((rows(_sd_a), rows(_sd_b)))
        mata: _sd_i = select(1::_sd_k, _sd_a[1::_sd_k] :!= _sd_b[1::_sd_k])
        mata: if (rows(_sd_i)) printf("first difference, line %g:\n  %s\n  %s\n", _sd_i[1], _sd_a[_sd_i[1]], _sd_b[_sd_i[1]])
        mata: mata drop _sd_a _sd_b _sd_k _sd_i
        exit 9
    }
    mata: mata drop _sd_a _sd_b
end

capture program drop plot_data
program define plot_data
    // plot_data graphname k  --  the serset of plot k in memory, variables v1, v2, ...; locals
    // v_y v_y2 v_x v_type (as tests/test_plot.do reads them)
    args G k
    local P .`G'.plotregion1.plot`k'
    local id ``P'.serset.id'
    confirm integer number `id'
    qui serset set `id'
    qui serset use, clear
    mata: for (i = 1; i <= st_nvar(); i++) st_varrename(i, "v" + strofreal(i))
    foreach m in y y2 x {
        local i ``P'.`m'var'
        if "`i'" == "" | "`i'" == "." c_local v_`m'
        else c_local v_`m' v`i'
    }
    c_local v_type ``P'.type.snm'
end

capture program drop check_contrib
program define check_contrib
    // check_contrib graphname [, spikes vertical]  --  the marginal plot of individually: rows
    // (y-axis labels from the top; x-axis labels from the left with vertical) Total diff, Expl.,
    // the blocks, Remainder, Unexpl.; one bar per row from 0 to its e(b) value; spikes on every
    // row with spikes (makeplot(, ci))
    syntax name(name=G) [, spikes vertical ]
    local ax = cond("`vertical'" != "", "xaxis1", "yaxis1")
    tempname b R
    matrix `b' = e(b)
    local names overall:difference overall:explained
    local rls `""Total diff" "Expl.""'
    foreach bn in `e(blocks)' {
        local names `names' marg_expl:`bn'
        local rls `"`rls' "`bn'""'
    }
    local names `names' marg_expl:_remainder overall:unexplained
    local rls `"`rls' "Remainder" "Unexpl.""'
    local nr : word count `names'
    assert `nr' == e(N_stages) + 3
    matrix `R' = J(`nr', 3, 0)
    forvalues r = 1/`nr' {
        local nm : word `r' of `names'
        matrix `R'[`r', 1] = `b'[1, colnumb(`b', "`nm'")]
        local t `.`G'.`ax'.major.dlg_ticklabel[`r']'
        local rl : word `r' of `rls'
        assert `"`t'"' == `"`rl'"'
        matrix `R'[`r', 2] = `.`G'.`ax'.major.dlg_tickpos[`r']'
    }
    local np 0
    forvalues k = 1/30 {
        if "`.`G'.plotregion1.plot`k'.classname'" == "" continue, break
        local np `k'
    }
    preserve
    local ncap 0
    forvalues k = 1/`np' {
        plot_data `G' `k'
        if !inlist("`v_type'", "rbar", "rcap") continue
        forvalues i = 1/`=_N' {
            local row 0
            forvalues r = 1/`nr' {
                if reldif(`v_x'[`i'], `R'[`r', 2]) < 1e-12 local row `r'
            }
            assert `row' > 0
            if "`v_type'" == "rbar" {
                assert `v_y'[`i'] == 0
                assert reldif(`v_y2'[`i'], `R'[`row', 1]) < 1e-10
                matrix `R'[`row', 3] = `R'[`row', 3] + 1
            }
            else local ++ncap
        }
    }
    restore
    * every row has exactly one bar
    forvalues r = 1/`nr' {
        assert `R'[`r', 3] == 1
    }
    if "`spikes'" != "" assert `ncap' == `nr'
    else assert `ncap' == 0
end

tempfile lg bsf e1 e2
local lgf `"`lg'.log"'
local e1f `"`e1'.txt"'
local e2f `"`e2'.txt"'
capture log close flog

* ---- data ---------------------------------------------------------------------
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
recast double y x3 wt
assert _N == 2000
gen double d = 1 + exp(x3 / 3) + x2           // a positive denominator
gen byte x1z = x1
replace x1z = 9 if group == 0 & id <= 60      // a level found only in the Focal Group (30 observations)
count if x1z == 9
assert r(N) == 30
tempname B B1 B2 S S1 S2 V SE SE1 C C1 C2
local spec (a: x1) (b: x2 c.x3) (c: x4)
local cv contvars(x3 3)


* ==============================================================================
* 1. parsing and errors
* ==============================================================================
di as txt _n "---- 1. parsing and errors"
* sequentially and individually exclude each other; reverse, detail and the block name
* _remainder are refused with individually, with or without fixothermargins
rcof "noi seldecomp y (a: x1) (b: x2), by(group) individually sequentially" == 198
rcof "noi seldecomp y (a: x1) (b: x2), by(group) sequentially individually fixothermargins" == 198
foreach v in individually "individually fixothermargins" {
    rcof "noi seldecomp y (a: x1) (b: x2), by(group) `v' reverse" == 198
    rcof "noi seldecomp y (a: x1) (b: x2), by(group) `v' detail" == 198
    rcof "noi seldecomp y (a: x1) (_remainder: x2), by(group) `v'" == 198
}
* the words isolated and nohold are refused, and so are cumulative, holdothers and hold next
* to individually
foreach v in isolated nohold "individually cumulative" "individually holdothers" "individually hold" {
    rcof "noi seldecomp y (a: x1) (b: x2), by(group) `v'" == 198
}
* a variable in two blocks: refused with fixothermargins (it would be raked and held at
* once), allowed without it, as with sequentially alone
rcof "noi seldecomp y (a: x1) (b: x1#x2), by(group) individually fixothermargins" == 198
qui seldecomp y (a: x1) (b: x1#x2), by(group) individually
assert "`e(stagetype)'" == "individually" & "`e(fixothermargins)'" == ""
* _remainder is reserved under individually only
qui seldecomp y (a: x1) (_remainder: x2), by(group) sequentially fixothermargins
assert "`e(blocks)'" == "a _remainder"
* the minimal abbreviations ind and fix; e(stagetype) and e(fixothermargins), and no
* e(varbyvar), e(sequential) or e(holdothers)
qui seldecomp y (a: x1) (b: x2), by(group) ind fix
assert "`e(stagetype)'" == "individually" & "`e(fixothermargins)'" == "fixothermargins"
assert "`e(varbyvar)'`e(sequential)'`e(holdothers)'" == ""
qui seldecomp y (a: x1) (b: x2), by(group) ind
assert "`e(stagetype)'" == "individually" & "`e(fixothermargins)'" == ""
* the messages (at linesize 160, so that each is on one line of the log)
local ls0 = c(linesize)
set linesize 160
log using `"`lgf'"', text replace name(flog)
capture noisily seldecomp y (a: x1) (b: x2), by(group) individually reverse
capture noisily seldecomp y (a: x1) (b: x1#x2), by(group) individually fixothermargins
capture noisily seldecomp y (a: x1) (b: x2), by(group) individually sequentially
log close flog
set linesize `ls0'
log_count `"`lgf'"' "reverse may not be combined with individually"
assert r(n) == 1
* the variable and its two blocks; no remedy
log_count `"`lgf'"' "with fixothermargins, a variable may be in one block only: x1 is in blocks a and b"
assert r(n) == 1
log_count `"`lgf'"' "drop fixothermargins"
assert r(n) == 0
log_count `"`lgf'"' "sequentially and individually may not be combined"
assert r(n) == 1


* ==============================================================================
* 2. individually fixothermargins (the others held at Focal Group's margins)
* ==============================================================================
di as txt _n "---- 2. individually fixothermargins"
seldecomp y `spec' [pw=wt], by(group) `cv' individually fixothermargins generate(ri)
matrix `B' = e(b)
matrix `S' = e(bystage)
assert "`e(stagetype)'" == "individually" & "`e(fixothermargins)'" == "fixothermargins"
assert e(N_stages) == 4 & e(N_blocks) == 3 & e(k_eq) == 3
assert "`e(stagenames)'" == "a b c overall"
assert "`e(stage1)'" == "x1" & "`e(stage2)'" == "x2 x3" & "`e(stage4)'" == "x1 x2 x3 x4"
assert "`: colfullnames `B''" == "overall:focal_group overall:raked_focal overall:ref_group overall:difference overall:explained overall:unexplained marg_expl:a marg_expl:b marg_expl:c marg_expl:_remainder raked_focal:a raked_focal:b raked_focal:c"
assert "`: rownames `S''" == "a b c overall"
local cn : colnames `S'
assert "`: word 2 of `cn''" == "explained" & "`: word 3 of `cn''" == "unexplained"
assert "`: word 5 of `cn''" == "contribution" & "`: word 27 of `cn''" == "pct_contribution"
assert colsof(`S') == 30 & rowsof(`S') == 4
* the stage matrix: column 5 is the block's own explained part, the remainder last; column 5
* sums to the explained part; the overall row is the full raking
local sum 0
forvalues k = 1/3 {
    assert_close `S'[`k',5] `S'[`k',2] 0 "individually fixothermargins: contribution = explained, block `k'"
    local bn : word `k' of a b c
    bcol v `B' marg_expl:`bn'
    assert_close `v' `S'[`k',5] 0 "individually fixothermargins: e(b) block `k'"
    assert_close `S'[`k',4] 100*`S'[`k',2]/e(difference) `same'
    assert_close `S'[`k',27] 100*`S'[`k',5]/e(difference) `same'
    assert_close `S'[`k',3] `S'[`k',1]-e(mean_ref) `same'
    local sum = `sum' + `S'[`k', 5]
}
assert_close `S'[4,5] `S'[4,2]-`sum' `same' "individually fixothermargins: remainder"
bcol v `B' marg_expl:_remainder
assert_close `v' `S'[4,5] 0
assert_close `S'[4,2] e(explained) 0
assert_close `S'[4,1] e(mean_raked_focal) 0
assert_close `sum'+`S'[4,5] e(explained) `same' "individually fixothermargins: the contributions and the remainder sum to explained"
assert_close e(pct_increment) `S'[4,27] 0
* raking k is stage 1 of sequentially fixothermargins with block k first and the other blocks
* after it in their order, to the last digit
local r1 (a: x1) (b: x2 c.x3) (c: x4)
local r2 (b: x2 c.x3) (a: x1) (c: x4)
local r3 (c: x4) (a: x1) (b: x2 c.x3)
forvalues k = 1/3 {
    local bn : word `k' of a b c
    qui seldecomp y `r`k'' [pw=wt], by(group) `cv' sequentially fixothermargins
    bcol v `B' marg_expl:`bn'
    assert_close _b[marg_expl:`bn'] `v' `same' "individually fixothermargins block `bn' vs sequentially fixothermargins stage 1"
    assert_close el(e(bystage),1,1) `S'[`k',1] `same'
}
* the full raking is the one-stage run, to the last digit
seldecomp y x1 x2 c.x3 x4 [pw=wt], by(group) `cv' generate(rp)
matrix `B1' = e(b)
matrix `B2' = `B'[1, 1..6]
mat_close `B1' `B2' `same'
assert reldif(ri, rp) < `same' if !missing(rp)
assert missing(ri) == missing(rp)
* generate(): ri_s# holds raking #: its raked mean is that row's raked_focal; block k at the Reference Group's
* shares, the other blocks at the Focal Group's own (base-weight) shares
forvalues k = 1/3 {
    qui sum y [aw=ri_s`k'] if group == 0
    assert_close r(mean) `S'[`k',1] `mtol' "generate(): raking `k'"
    assert "`: char ri_s`k'[seldecomp_stage]'" == "`k'"
}
foreach v in x1 x2 x4 {
    qui levelsof `v', local(levs)
    foreach l of local levs {
        if "`v'" == "x4" local cond x4 == "`l'"
        else local cond `v' == `l'
        qui sum wt if group == 1 & `cond', meanonly
        local s2 = r(sum)
        qui sum wt if group == 1, meanonly
        local s2 = `s2' / r(sum)
        qui sum wt if group == 0 & `cond', meanonly
        local s1 = r(sum)
        qui sum wt if group == 0, meanonly
        local s1 = `s1' / r(sum)
        qui sum ri_s1 if group == 0 & `cond', meanonly
        local r1s = r(sum)
        qui sum ri_s1 if group == 0, meanonly
        local r1s = `r1s' / r(sum)
        * raking 1 (block a, x1): x1 at the Reference Group's shares, x2 and x4 held at the Focal Group's
        if "`v'" == "x1" assert_close `r1s' `s2' `tol' "individually fixothermargins raking 1: `v' = `l' at the Reference Group's share"
        else assert_close `r1s' `s1' `tol' "individually fixothermargins raking 1: `v' = `l' held at the Focal Group's share"
    }
}
* the order of the blocks: the same contributions and remainder, to the raking tolerance (the
* held margins reach the raking in another order)
seldecomp y (c: x4) (b: x2 c.x3) (a: x1) [pw=wt], by(group) `cv' individually fixothermargins
foreach bn in a b c _remainder {
    bcol v `B' marg_expl:`bn'
    assert_close _b[marg_expl:`bn'] `v' `rtol' "individually fixothermargins, blocks reversed: `bn'"
}
drop ri* rp


* ==============================================================================
* 3. individually (the default: the other blocks free)
* ==============================================================================
di as txt _n "---- 3. individually"
seldecomp y `spec' [pw=wt], by(group) `cv' individually generate(rn)
matrix `B2' = e(b)
matrix `S2' = e(bystage)
assert "`e(stagetype)'" == "individually" & "`e(fixothermargins)'" == ""
assert e(N_stages) == 4 & e(k_eq) == 3
assert "`: colfullnames `B2''" == "`: colfullnames `B''"
* raking k is the plain run on block k, to the last digit
* (the plain run on block k: its terms, the weights and its options)
local p1 x1 [pw=wt],
local p2 x2 c.x3 [pw=wt], `cv'
local p3 x4 [pw=wt],
forvalues k = 1/3 {
    local bn : word `k' of a b c
    qui seldecomp y `p`k'' by(group)
    bcol v `B2' marg_expl:`bn'
    assert_close _b[overall:explained] `v' `same' "individually block `bn' vs the plain run"
    assert_close e(mean_raked_focal) `S2'[`k',1] `same'
}
* the full raking is that of individually fixothermargins, to the last digit
matrix `B1' = `B2'[1, 1..6]
matrix `C' = `B'[1, 1..6]
mat_close `B1' `C' `same'
* the remainder identity
local sum 0
foreach bn in a b c {
    bcol v `B2' marg_expl:`bn'
    local sum = `sum' + `v'
}
bcol v `B2' overall:explained
local ex `v'
bcol v `B2' marg_expl:_remainder
assert_close `v' `ex'-`sum' `same' "individually: remainder"
* the other blocks drift: x2's Focal Group share under raking 1 is not its base share
qui sum wt if group == 0 & x2 == 1, meanonly
local s1 = r(sum)
qui sum wt if group == 0, meanonly
local s1 = `s1' / r(sum)
qui sum rn_s1 if group == 0 & x2 == 1, meanonly
local r1s = r(sum)
qui sum rn_s1 if group == 0, meanonly
assert reldif(`r1s' / r(sum), `s1') > 1e-3
* the order: exactly the same contributions (each raking involves its block alone); the
* remainder to the raking tolerance (the full raking visits the margins in another order)
seldecomp y (c: x4) (b: x2 c.x3) (a: x1) [pw=wt], by(group) `cv' individually
foreach bn in a b c {
    bcol v `B2' marg_expl:`bn'
    assert_close _b[marg_expl:`bn'] `v' `same' "individually, blocks reversed: `bn'"
}
bcol v `B2' marg_expl:_remainder
assert_close _b[marg_expl:_remainder] `v' `rtol'
* a variable in two blocks: each block is its plain run
seldecomp y (a: x1) (b: x1#x2) [pw=wt], by(group) individually
matrix `B1' = e(b)
qui seldecomp y x1#x2 [pw=wt], by(group)
bcol v `B1' marg_expl:b
assert_close _b[overall:explained] `v' `same'
drop rn*


* ==============================================================================
* 4. outcomes(), ratio()
* ==============================================================================
di as txt _n "---- 4. outcomes(), ratio()"
local sp3 (a: x1) (b: x2) (c: x4)
foreach v in "individually fixothermargins" individually {
    qui seldecomp y `sp3' [pw=wt], by(group) `v'
    matrix `B1' = e(b)
    * outcomes(): each outcome's equations are its own run's
    qui seldecomp y `sp3' [pw=wt], by(group) `v' outcomes(ybin)
    matrix `B2' = e(b)
    assert e(k_eq) == 6
    assert "`: coleq `B2''" == "y y y y y y y_marg_expl y_marg_expl y_marg_expl y_marg_expl y_raked_focal y_raked_focal y_raked_focal ybin ybin ybin ybin ybin ybin ybin_marg_expl ybin_marg_expl ybin_marg_expl ybin_marg_expl ybin_raked_focal ybin_raked_focal ybin_raked_focal"
    matrix `C' = `B2'[1, 1..13]
    mat_close `C' `B1' `same'
    matrix `S1' = e(bystage2)
    assert "`: rownames `S1''" == "a b c overall"
    assert "`: word 5 of `: colnames `S1'''" == "contribution"
    qui seldecomp ybin `sp3' [pw=wt], by(group) `v'
    matrix `B1' = e(b)
    matrix `S' = e(bystage)
    matrix `C' = `B2'[1, 14..26]
    mat_close `C' `B1' `same'
    matrix `C' = `S1'[1..., 1..5]
    matrix `C1' = `S'[1..., 1..5]
    mat_close `C' `C1' `same'
    * ratio(): raking k's raked ratio from its generate() weights
    qui seldecomp y `sp3' [pw=wt], by(group) `v' ratio(d) generate(rr)
    matrix `S' = e(bystage)
    forvalues k = 1/4 {
        local wv = cond(`k' < 4, "rr_s`k'", "rr")
        qui sum y [aw=`wv'] if group == 0
        local a = r(mean)
        qui sum d [aw=`wv'] if group == 0
        assert_close `a'/r(mean) `S'[`k',1] `mtol' "`v' ratio(): raking `k'"
    }
    assert_close `S'[4,5] `S'[4,2]-`S'[1,5]-`S'[2,5]-`S'[3,5] `same'
    drop rr*
}
* individually fixothermargins, ratio(): raking 2 is stage 1 of the ratio run with block b
* first (sequentially fixothermargins)
qui seldecomp y `sp3' [pw=wt], by(group) individually fixothermargins ratio(d)
matrix `B1' = e(b)
qui seldecomp y (b: x2) (a: x1) (c: x4) [pw=wt], by(group) sequentially fixothermargins ratio(d)
bcol v `B1' marg_expl:b
assert_close _b[marg_expl:b] `v' `same'


* ==============================================================================
* 5. bayesboot()
* ==============================================================================
di as txt _n "---- 5. bayesboot()"
* individually: the SE of block k's contribution is that of the explained part of the plain
* run on block k with the same seed (the same draws); the overall SEs are the one-stage run's
seldecomp y `spec' [pw=wt], by(group) `cv' individually bayesboot(reps(20) rseed(11) nodots)
matrix `SE' = e(se)
local bsn `e(bsnames)'
assert wordcount("`bsn'") == 6 + 4 + 3
assert "`: word 10 of `bsn''" == "marg_expl:_remainder"
assert "`: word 7 of `bsn''" == "marg_expl:a"
assert "`: word 11 of `bsn''" == "raked_focal:a" & "`: word 13 of `bsn''" == "raked_focal:c"
assert e(N_misreps) == 0
matrix `V' = e(V)
forvalues j = 1/13 {
    assert `V'[`j', `j'] > 0
}
forvalues k = 1/3 {
    local bn : word `k' of a b c
    qui seldecomp y `p`k'' by(group) bayesboot(reps(20) rseed(11) nodots)
    bcol v `SE' marg_expl:`bn'
    assert_close el(e(se),1,5) `v' `same' "individually bootstrap SE, block `bn'"
}
qui seldecomp y x1 x2 c.x3 x4 [pw=wt], by(group) `cv' bayesboot(reps(20) rseed(11) nodots)
matrix `SE1' = e(se)
matrix `C' = `SE'[1, 1..6]
mat_close `C' `SE1' `same'
* individually fixothermargins: the SE of block k's contribution is that of
* marg_expl:<block k> in the block-k-first sequentially fixothermargins run with
* bayesboot() and the same seed
seldecomp y `spec' [pw=wt], by(group) `cv' individually fixothermargins bayesboot(reps(20) rseed(11) nodots)
matrix `SE' = e(se)
matrix `C' = `SE'[1, 1..6]
mat_close `C' `SE1' `same'
forvalues k = 2/3 {
    local bn : word `k' of a b c
    qui seldecomp y `r`k'' [pw=wt], by(group) `cv' sequentially fixothermargins bayesboot(reps(20) rseed(11) nodots)
    matrix `C2' = e(se)
    bcol w `C2' marg_expl:`bn'
    bcol v `SE' marg_expl:`bn'
    assert_close `w' `v' `same' "individually fixothermargins bootstrap SE, block `bn'"
}
* every saving() row is a plain individually fixothermargins run on that replicate's
* generate() weights
seldecomp y `sp3' [pw=wt], by(group) individually fixothermargins ///
    bayesboot(reps(3) rseed(5) generate(bw) saving(`"`bsf'"', replace) nodots)
local bsn `e(bsnames)'
preserve
use `"`bsf'"', clear
assert _N == 3
mkmat _all, matrix(`C')
restore
forvalues r = 1/3 {
    qui seldecomp y `sp3' [iw=bw_`r'], by(group) individually fixothermargins
    matrix `B1' = e(b)
    local j 0
    foreach fn of local bsn {
        local ++j
        bcol v `B1' `fn'
        assert_close `C'[`r',`j'] `v' 1e-10 "saving() row `r', `fn'"
    }
}
drop bw_*
* bca with clusters: an acceleration for every replicated element
qui seldecomp y `sp3' [pw=wt], by(group) individually bayesboot(reps(10) rseed(3) cluster(cluster) bca nodots)
matrix `C' = e(accel)
assert colsof(`C') == 13
forvalues j = 1/13 {
    assert !missing(`C'[1, `j'])
}


* ==============================================================================
* 6. prefixes
* ==============================================================================
di as txt _n "---- 6. prefixes"
preserve
keep if cluster <= 20
svyset cluster [pw=wt]
svy jackknife, nodots: seldecomp y `sp3', by(group) individually fixothermargins
matrix `B2' = e(b)
matrix `V' = e(V)
assert colsof(`B2') == 13
forvalues j = 1/13 {
    assert `V'[`j', `j'] > 0
}
qui seldecomp y `sp3' [pw=wt], by(group) individually fixothermargins
matrix `B1' = e(b)
mat_close `B1' `B2' `same'
svy jackknife, nodots: seldecomp y `sp3', by(group) individually
log using `"`lgf'"', text replace name(flog)
seldecomp
seldecomp, makeplot(, ci twoway(name(vsvy, replace)))
log close flog
log_count `"`lgf'"' "Contribution of each block, raked alone (individually)"
assert r(n) == 1
check_contrib vsvy, spikes
graph drop vsvy
svyset, clear
restore
* bootstrap:, then the plot at replay (every row with a spike)
bootstrap, reps(3) seed(2) nodots: seldecomp y `sp3', by(group) individually fixothermargins
matrix `V' = e(V)
forvalues j = 1/10 {
    assert `V'[`j', `j'] > 0
}
seldecomp, makeplot(, ci twoway(name(vbs, replace)))
check_contrib vbs, spikes
graph drop vbs
* by: each group's results are those of the if-restricted run
bysort x2: seldecomp y (a: x1) (c: x4) [pw=wt], by(group) individually
matrix `B2' = e(b)
qui seldecomp y (a: x1) (c: x4) [pw=wt] if x2 == 1, by(group) individually
matrix `B1' = e(b)
mat_close `B1' `B2' `same'


* ==============================================================================
* 7. display and makeplot
* ==============================================================================
di as txt _n "---- 7. display and makeplot"
local ls0 = c(linesize)
set linesize 160
log using `"`lgf'"', text replace name(flog)
di "SDW_W_BEGIN"
seldecomp y `spec' [pw=wt], by(group) `cv' individually fixothermargins bayesboot(reps(10) rseed(2) nodots citable)
di "SDW_W_END"
log close flog
set linesize `ls0'
log_width `"`lgf'"'
di as txt "widest display line: " r(max)
assert r(max) <= 78
foreach t in "Blocks raked individually (others held at focal group margins):" ///
    "block 1: a (x1)" "block 2: b (x2 x3)" "Inference: Bayesian bootstrap" ///
    "Contribution of each block, raked alone (individually, fixothermargins)" "Sum of the blocks" ///
    "All blocks (overall)" "Remainder" "  _remainder |" ///
    "%: percent of the difference. Other blocks are held at their focal group" {
    log_count `"`lgf'"' "`t'"
    if r(n) < 1 {
        di as err `"display: "`t'" not found"'
        exit 9
    }
}
* (and none of the refused option words; no section heading "Decomposition", which is printed
* only when no stage table follows, and no other wording of the header line)
foreach t in "cum_expl" "Decomposition" "Blocks added" ///
    "isolated" "nohold" "varbyvar" "holdothers" "Block contributions" {
    log_count `"`lgf'"' "`t'"
    assert r(n) == 0
}
* the table lines: the rows of e(bystage); the same when rebuilt from e(b) (a prefix that does
* not carry the stage matrix)
log using `"`lgf'"', text replace name(flog)
seldecomp
log close flog
contrib_lines `"`lgf'"'
local tl1 `"`r(lines)'"'
drop_stagemats
log using `"`lgf'"', text replace name(flog)
seldecomp
log close flog
contrib_lines `"`lgf'"'
assert `"`r(lines)'"' == `"`tl1'"'
* with bayesboot(): the header names the inference
log using `"`lgf'"', text replace name(flog)
seldecomp y `sp3' [pw=wt], by(group) individually bayesboot(reps(10) rseed(2) nodots)
log close flog
log_count `"`lgf'"' "Inference: Bayesian bootstrap"
assert r(n) == 1
log_count `"`lgf'"' "Blocks raked individually (others unconstrained):"
assert r(n) == 1
log_count `"`lgf'"' "%: percent of the difference. Other blocks are unconstrained. Remainder:"
assert r(n) == 1
log_count `"`lgf'"' "Contribution of each block, raked alone (individually)"
assert r(n) == 1
* makeplot: the marginal plot with a remainder row, the default and the only plot type, at
* estimation and at replay, swapped to vertical bars with swapaxis; no note; spikes on every
* row with bayesboot() and ci, none without ci
seldecomp y `sp3' [pw=wt], by(group) individually fixothermargins makeplot(, twoway(name(v1, replace)))
check_contrib v1
assert `"`.v1.note.text[1]'"' == ""
graph drop v1
seldecomp, makeplot(marginal, twoway(name(v2, replace)))
check_contrib v2
graph drop v2
seldecomp y `sp3' [pw=wt], by(group) individually makeplot(marginal, swapaxis twoway(name(v3, replace)))
check_contrib v3, vertical
graph drop v3
foreach t in stages waterfall waterfalllevels {
    rcof "noi seldecomp, makeplot(`t')" == 198
    rcof "noi seldecomp y `sp3' [pw=wt], by(group) individually makeplot(`t')" == 198
}
seldecomp y `sp3' [pw=wt], by(group) individually bayesboot(reps(10) rseed(2) nodots)
seldecomp, makeplot(, ci twoway(name(v4, replace)))
check_contrib v4, spikes
seldecomp, makeplot(, twoway(name(v5, replace)))
check_contrib v5
graph drop v4 v5


* ==============================================================================
* 8. insupport
* ==============================================================================
di as txt _n "---- 8. insupport"
gen byte excl = group == 0 & x1z == 9
* the screened sample: the same as the run restricted with -if-, to the last digit, with the
* unrestricted Focal Group kept in e(N_focal_all), e(W_focal_all), e(mean_focal_all)
seldecomp y x1z x2 x4 [pw=wt], by(group) insupport generate(rg) covariatedetail
matrix `B1' = e(b)
matrix `S1' = e(bystage)
matrix `C1' = e(covariates)
local lv1 `"`e(levels1)'"'
assert "`e(insupport)'" == "insupport"
assert e(N_focal) == 888 - 30 & e(N_focal_all) == 888 & e(N) == 2000 - 30
assert e(N_focal_zero) == 0 & `S1'[1, 21] == 0
capture confirm scalar e(explained_support)
assert _rc
gen byte s1 = e(sample)
qui sum y [aw=wt] if group == 0
assert_close e(mean_focal_all) r(mean) `mtol'
assert_close e(W_focal_all) r(sum_w) `mtol'
qui sum y [aw=wt] if group == 0 & !excl
assert_close e(mean_focal) r(mean) `mtol' "insupport: focal_group is the mean on the screened sample"
assert "`e(note1)'" == "insupport excluded x1z = 9 (30 obs)"
assert "`e(note2)'" != "    x1z = 9 (30 obs)"
seldecomp y x1z x2 x4 [pw=wt] if !excl, by(group) generate(rh) covariatedetail
matrix `B2' = e(b)
matrix `S2' = e(bystage)
matrix `C2' = e(covariates)
mat_close `B1' `B2' 0
mat_close `S1' `S2' 0
mat_close `C1' `C2' 0
assert `"`lv1'"' == `"`e(levels1)'"' & `"`lv1'"' == "1 | 2 | 3"
assert s1 == e(sample)
assert rg == rh if !missing(rh)
assert missing(rg) == missing(rh)
assert missing(rg) if excl
drop s1 rg rh
* without insupport: the change of support (zero raked weight) instead
qui seldecomp y x1z x2 x4 [pw=wt], by(group)
assert e(N_focal_zero) == 30 & e(N_focal) == 888
* no level absent from the Reference Group: every e() the same with and without insupport, apart from
* the four insupport results (elist-style), with sequentially fixothermargins,
* bayesboot(), generate() and covariatedetail
local nsc seldecomp y (a: x1) (b: x2 c.x3) (c: x4) [pw=wt], by(group) `cv' sequentially fixothermargins ///
    bayesboot(reps(10) rseed(4) nodots) covariatedetail generate(rq)
qui `nsc'
gen byte s1 = e(sample)
rename (rq rq_s1 rq_s2) (rq0 rq0_s1 rq0_s2)
e_dump `"`e1f'"'
qui `nsc' insupport
assert e(N_focal_all) == e(N_focal) & e(W_focal_all) == e(W_focal) & e(mean_focal_all) == e(mean_focal)
e_dump `"`e2f'"', skip(insupport N_focal_all W_focal_all mean_focal_all)
same_file `"`e1f'"' `"`e2f'"'
assert s1 == e(sample)
foreach v in "" _s1 _s2 {
    assert rq`v' == rq0`v' | (missing(rq`v') & missing(rq0`v'))
}
drop s1 rq*
qui seldecomp y `sp3' [pw=wt], by(group) individually fixothermargins
e_dump `"`e1f'"'
qui seldecomp y `sp3' [pw=wt], by(group) individually fixothermargins insupport
e_dump `"`e2f'"', skip(insupport N_focal_all W_focal_all mean_focal_all)
same_file `"`e1f'"' `"`e2f'"'
* the held shares: with the Focal Group level in block 2, stage 1 of sequentially fixothermargins
* holds block 2 at the screened sample's shares with insupport; without it, at the whole Group
* 1's (nothing is zeroed at stage 1: the rule restricts only at the stage that rakes the
* variable), and at stage 2 the held block 3 at the shares of the observations kept
seldecomp y (a: x2) (b: x1z) (c: x4) [pw=wt], by(group) sequentially fixothermargins insupport generate(ra)
seldecomp y (a: x2) (b: x1z) (c: x4) [pw=wt], by(group) sequentially fixothermargins generate(rb)
assert el(e(bystage), 1, 21) == 0 & el(e(bystage), 2, 21) == 30
foreach v in x1z x4 {
    qui levelsof `v' if group == 0, local(levs)
    foreach l of local levs {
        if "`v'" == "x4" local cond x4 == "`l'"
        else local cond `v' == `l'
        * the screened Focal Group's base shares, and the whole Focal Group's
        qui sum wt if group == 0 & !excl & `cond', meanonly
        local sk = r(sum)
        qui sum wt if group == 0 & !excl, meanonly
        local sk = `sk' / r(sum)
        qui sum wt if group == 0 & `cond', meanonly
        local sa = r(sum)
        qui sum wt if group == 0, meanonly
        local sa = `sa' / r(sum)
        * insupport, stage 1: held at the screened shares
        qui sum ra_s1 if group == 0 & `cond', meanonly
        local q = r(sum)
        qui sum ra_s1 if group == 0, meanonly
        if "`v'`l'" == "x1z9" assert `q' == 0
        else assert_close `q'/r(sum) `sk' `tol' "insupport, stage 1: `v' = `l'"
        * default, stage 1: held at the whole Focal Group's shares (level 9 included)
        qui sum rb_s1 if group == 0 & `cond', meanonly
        local q = r(sum)
        qui sum rb_s1 if group == 0, meanonly
        assert_close `q'/r(sum) `sa' `tol' "default, stage 1: `v' = `l'"
        * default, stage 2: x1z raked (level 9 zeroed), x4 held at the kept observations' shares
        if "`v'" == "x4" {
            qui sum rb_s2 if group == 0 & `cond', meanonly
            local q = r(sum)
            qui sum rb_s2 if group == 0, meanonly
            assert_close `q'/r(sum) `sk' `tol' "default, stage 2: `v' = `l'"
        }
    }
}
qui count if group == 0 & x1z == 9 & rb_s1 > 0
assert r(N) == 30
drop ra* rb*
* individually fixothermargins: no zeroing at any raking
qui seldecomp y (a: x1z) (b: x2) (c: x4) [pw=wt], by(group) individually fixothermargins insupport
matrix `S' = e(bystage)
forvalues k = 1/4 {
    assert `S'[`k', 21] == 0
}
matrix `B1' = e(b)
qui seldecomp y (a: x1z) (b: x2) (c: x4) [pw=wt] if !excl, by(group) individually fixothermargins
matrix `B2' = e(b)
mat_close `B1' `B2' 0
* bayesboot(): the screen is done once, on the estimation sample; the draws are those of the
* if-restricted run with the same seed
qui seldecomp y (a: x1z) (b: x2) [pw=wt], by(group) sequentially fixothermargins insupport ///
    bayesboot(reps(10) rseed(8) nodots)
matrix `V' = e(V)
qui seldecomp y (a: x1z) (b: x2) [pw=wt] if !excl, by(group) sequentially fixothermargins ///
    bayesboot(reps(10) rseed(8) nodots)
matrix `C' = e(V)
mat_close `V' `C' `same'
* by: each by-group is screened against its own Reference Group
bysort strata: seldecomp y x1z x2 [pw=wt], by(group) insupport
matrix `B2' = e(b)
local n1a = e(N_focal_all)
qui seldecomp y x1z x2 [pw=wt] if strata == 3, by(group) insupport
matrix `B1' = e(b)
mat_close `B1' `B2' 0
assert e(N_focal_all) == `n1a'
* svy jackknife: every replicate screens again; the point estimates are the plain run's
preserve
keep if cluster <= 20
svyset cluster [pw=wt]
svy jackknife, nodots: seldecomp y x1z x2, by(group) insupport
matrix `B2' = e(b)
matrix `V' = e(V)
forvalues j = 1/6 {
    assert `V'[`j', `j'] > 0
}
qui seldecomp y x1z x2 [pw=wt], by(group) insupport
matrix `B1' = e(b)
mat_close `B1' `B2' `same'
assert "`e(insupport)'" == "insupport"
svyset, clear
restore
* outcomes() and ratio(): the unrestricted Focal Group mean of each outcome, a ratio with ratio()
qui seldecomp y x1z x2 [pw=wt], by(group) insupport outcomes(ybin)
matrix `C' = e(by_outcome)
assert colsof(`C') == 11 & "`: word 11 of `: colnames `C'''" == "mean_focal_all"
qui sum ybin [aw=wt] if group == 0
assert_close `C'[2,11] r(mean) `mtol'
assert_close `C'[1,11] e(mean_focal_all) 0
qui seldecomp y x1z x2 [pw=wt], by(group) outcomes(ybin)
assert colsof(e(by_outcome)) == 10
qui seldecomp y x1z x2 [pw=wt], by(group) insupport ratio(d)
qui sum y [aw=wt] if group == 0
local a = r(mean)
qui sum d [aw=wt] if group == 0
assert_close e(mean_focal_all) `a'/r(mean) `mtol'
* the cut points are those of the unscreened sample (the screen needs the bins)
qui seldecomp y x1z c.x3 [pw=wt], by(group) contvars(x3 4)
matrix `C1' = e(cuts1)
qui seldecomp y x1z c.x3 [pw=wt], by(group) contvars(x3 4) insupport
matrix `C2' = e(cuts1)
mat_close `C1' `C2' 0
assert e(N_focal) == 858
* the display: the header line and the note; nothing excluded
log using `"`lgf'"', text replace name(flog)
seldecomp y x1z x2 x4 [pw=wt], by(group) insupport
seldecomp y x1 x2 [pw=wt], by(group) insupport
log close flog
log_count `"`lgf'"' "Support (insupport): 30 of 888 focal group obs excluded (levels absent from"
assert r(n) == 1
log_count `"`lgf'"' "x1z = 9 (30 obs)"
assert r(n) == 1
log_count `"`lgf'"' "Support (insupport): no observation excluded"
assert r(n) == 1
* errors: the exclusion empties Focal Group (2000), or a level found in the Reference Group (459; the default
* rule fails at the full raking too)
gen byte z = cond(group == 0, 1, 2)
rcof "noi seldecomp y x1 z, by(group) insupport" == 2000
gen byte u = 1
replace u = 2 if (group == 0 & x1z == 9) | (group == 1 & id <= 100)
rcof "noi seldecomp y x1z u, by(group) insupport" == 459
rcof "noi seldecomp y x1z u, by(group)" == 459
* (at linesize 160, so that each message is on one line of the log)
local ls0 = c(linesize)
set linesize 160
log using `"`lgf'"', text replace name(flog)
capture noisily seldecomp y x1z u, by(group) insupport
capture noisily seldecomp y x1 z, by(group) insupport
log close flog
set linesize `ls0'
log_count `"`lgf'"' "insupport: after the exclusion, these levels of the reference group (group = 1) have no focal group observations:"
assert r(n) == 1
log_count `"`lgf'"' "    u = 2"
assert r(n) == 1
* no remedy after the list
log_count `"`lgf'"' "without insupport the full raking fails"
assert r(n) == 0
* the variable is named; no guess at its cause (whether z, which takes different values in
* the two groups, duplicates by())
log_count `"`lgf'"' "insupport excludes every focal group observation (levels of z absent from the reference group)"
assert r(n) == 1
log_count `"`lgf'"' "does it duplicate"
assert r(n) == 0
drop z u excl

di as res _n "PASSED: test_variants.do"

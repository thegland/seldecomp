* test_detail.do -- the option detail: the terms of
* each block added one at a time, under sequentially, with fixothermargins (the later terms and
* blocks held at the Focal Group's margins) and without it (left free, the default). The term increments
* sum to the block increments (1e-12); each term's increment equals that of the run with the
* block split into single-term blocks up to that term, or entirely (1e-12); every e() of the run
* without detail is unchanged by it, apart from the added detail-specific diagnostics, and
* the layout without detail has no detail entries; the
* bootstrap (bayesboot()) SEs of the detail rows equal the split run's with the same seed,
* those of the block rows the run's without detail; saving() rows; outcomes(), ratio(), reverse,
* insupport, recut, generate(); the prefixes; the display (the term rows, rebuilt from e(b)
* alone); the parsing rules (detail requires sequentially and is refused with individually; a
* variable in two terms of a block only without fixothermargins); the error paths (a term stage
* that cannot meet its margins under fixothermargins; force, also under bayesboot()).
*
* Run from the tests directory:   . do test_detail.do

version 16
clear
set more off
set varabbrev off
discard             // drop any seldecomp already in memory; the repository copy is what gets tested

* ---- tolerances ---------------------------------------------------------------
local same  1e-12    // the same computation reached by two routes
local mtol  1e-10    // means from summarize under the generate() weights vs seldecomp

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

capture program drop vcell
program define vcell
    // vcell lname matname "eq:name1" "eq:name2"  --  the element of the square matrix matname in
    // row eq:name1, column eq:name2
    args lname M n1 n2
    local r = rownumb(`M', "`n1'")
    local c = colnumb(`M', "`n2'")
    if missing(`r') | missing(`c') {
        di as err "vcell: `n1' or `n2' is not in `M'"
        exit 9
    }
    c_local `lname' = `M'[`r', `c']
end

capture program drop log_count
program define log_count, rclass
    // log_count "file" "text" [unwrap] -- count matching lines, or occurrences after
    // joining wrapped display paragraphs when unwrap is specified.
    args fn txt unwrap
    local n 0
    if "`unwrap'" == "unwrap" {
        mata: _sd_l = stritrim(invtokens(strtrim(cat(st_local("fn")))'))
        mata: st_local("n", strofreal((strlen(_sd_l) - strlen(subinstr(_sd_l, st_local("txt"), ""))) / strlen(st_local("txt"))))
        mata: mata drop _sd_l
    }
    else mata: st_local("n", strofreal(sum(strpos(cat(st_local("fn")), st_local("txt")) :> 0)))
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

capture program drop stage_lines
program define stage_lines, rclass
    // stage_lines "file" "title"  --  the lines of the first stage table with that title (from
    // the title to the rule that closes it), joined with "|" in r(lines), and r(n) their number
    args fn ttl
    mata: _sd_l = cat(st_local("fn"))
    mata: _sd_b = select(1::rows(_sd_l), _sd_l :== st_local("ttl"))
    mata: _sd_h = select(1::rows(_sd_l), _sd_l :== 78 * "-")
    mata: _sd_h = select(_sd_h, _sd_h :> _sd_b[1])
    mata: st_local("lines", invtokens(_sd_l[_sd_b[1]::_sd_h[3]]', "|"))
    mata: st_local("n", strofreal(_sd_h[3] - _sd_b[1] + 1))
    mata: mata drop _sd_l _sd_b _sd_h
    return local lines `"`lines'"'
    return scalar n = `n'
end

capture program drop drop_mats
program define drop_mats, eclass
    // leave out the stage matrices, as a prefix that does not carry them would
    foreach m in bystage bystage2 bystage_detail bystage_detail2 {
        mata: st_matrix("e(`m')", J(0, 0, .))
    }
end

capture program drop e_dump
program define e_dump
    // e_dump "file" [, skip(names) skipdetail] -- every e() scalar, macro and matrix.
    // skipdetail omits only detail-specific diagnostic sidecars and their level labels;
    // shared forward/reverse diagnostics remain part of the exact comparison.
    syntax anything(name=fn) [, skip(string) SKIPDETAIL]
    local skip `skip' cmdline
    if "`skipdetail'" != "" local skip `skip' raking_checks_detail raking_misses_detail
    tempname fh M
    file open `fh' using `fn', write text replace
    foreach t in scalars macros matrices {
        local names : e(`t')
        local names : list sort names
        local names : list names - skip
        foreach n of local names {
            if "`skipdetail'" != "" & regexm("`n'", "^raking_miss_detail[0-9]+_(var|level)$") continue
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
        mata: st_local("ndiff", strofreal(rows(_sd_i)))
        if `ndiff' > 0 mata: printf("first difference, line %g:\n  %s\n  %s\n", _sd_i[1], _sd_a[_sd_i[1]], _sd_b[_sd_i[1]])
        mata: mata drop _sd_a _sd_b _sd_k _sd_i
        exit 9
    }
    mata: mata drop _sd_a _sd_b
end

capture program drop sum_check
program define sum_check
    // sum_check [prefix]  --  in e(b) (equations prefixed by prefix with outcomes()), for every
    // block: its terms' detail elements sum to its explained element (1e-12, and 1e-12 of the
    // difference for a block increment near 0), and a single-term block's element is its block's
    // exactly; every term has a detail element, in the order typed; e(bystage_detail) holds them
    // in column 5
    args pfx
    tempname b
    matrix `b' = e(b)
    local dnames : rownames e(bystage_detail)
    local eq0 = cond("`pfx'" == "", "overall", substr("`pfx'", 1, strlen("`pfx'") - 1))
    local tt 0
    forvalues k = 1/`=e(N_blocks)' {
        local bn : word `k' of `e(blocks)'
        bcol be `b' `pfx'marg_expl:`bn'
        bcol dd `b' `eq0':difference
        local nbt : word count `e(block`k')'
        local s 0
        forvalues u = 1/`nbt' {
            local ++tt
            local dn : word `tt' of `dnames'
            bcol de `b' `pfx'detail:`dn'
            local s = `s' + `de'
        }
        if `nbt' == 1 assert `s' == `be'
        else if abs(`be') > 1e-6 * abs(`dd') assert_close `s' `be' 1e-12 "block `bn': sum of its terms"
        else assert abs(`s' - `be') < 1e-12 * abs(`dd')
    }
    assert `tt' == e(N_terms)
end

tempfile lg bsf bsf2 e1 e2
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
* a pair tied within each group (b = 1 - a): feasible as one block, not term by term under
* fixothermargins (a raked to the Reference Group with b held at the Focal Group's share); tc a copy of ta
gen byte ta = x2
gen byte tb = 1 - x2
gen byte tc = x2
tempname B B1 B2 S S1 S2 D D1 V V1 V2 SE SE1 C C1 C2 A A1
local cv contvars(x3 3)
local spec (a: x1) (b: x2 c.x3) (c: x4)
local split (a: x1) x2 c.x3 (c: x4)
local spec3 (a: x1) (b: x2 c.x3 x4)


* ==============================================================================
* 1. parsing and errors
* ==============================================================================
di as txt _n "---- 1. parsing and errors"
rcof "noi seldecomp y `spec', by(group) `cv' detail" == 198
rcof "noi seldecomp y `spec', by(group) `cv' fixothermargins detail" == 198
rcof "noi seldecomp y `spec', by(group) `cv' individually fixothermargins detail" == 198
rcof "noi seldecomp y `spec', by(group) `cv' individually detail" == 198
rcof "noi seldecomp y `spec', by(group) `cv' sequentially individually det" == 198
* a variable in two terms of a block: refused with fixothermargins, allowed without it (and
* without detail the block is one raking)
rcof "noi seldecomp y (a: x4) (b: x2 x1#x2), by(group) sequentially fixothermargins detail" == 198
qui seldecomp y (a: x4) (b: x2 x1#x2), by(group) sequentially fixothermargins
qui seldecomp y (a: x4) (b: x2 x1#x2), by(group) sequentially detail
assert "`e(detail)'" == "detail" & "`e(stagetype)'" == "sequentially" & "`e(fixothermargins)'" == ""
matrix `B' = e(b)
bcol v1 `B' detail:x1_x2
sum_check
* the abbreviation det; replay does not take it
qui seldecomp y `spec', by(group) `cv' sequentially fixothermargins det
assert "`e(detail)'" == "detail"
rcof "noi seldecomp, detail" == 198
* outcomes(): the equation y_detail would be taken twice
gen double y_detail = y + 1
rcof "noi seldecomp y `spec', by(group) `cv' sequentially fixothermargins detail outcomes(y_detail)" == 198
qui seldecomp y `spec', by(group) `cv' sequentially fixothermargins outcomes(y_detail)
drop y_detail
* the messages (at linesize 160, so that each is on one line of the log)
local ls0 = c(linesize)
set linesize 160
log using `"`lgf'"', text replace name(flog)
capture noisily seldecomp y `spec', by(group) `cv' detail
capture noisily seldecomp y `spec', by(group) `cv' individually detail
capture noisily seldecomp y (a: x4) (b: x2 x1#x2), by(group) sequentially fixothermargins detail
log close flog
set linesize `ls0'
foreach t in "detail requires sequentially" ///
    "detail may not be combined with individually" ///
    "with detail and fixothermargins, a variable may be in one term only: x2 is in x2 and x1#x2" {
    log_count `"`lgf'"' "`t'"
    if r(n) != 1 {
        di as err `"message not found: `t'"'
        exit 9
    }
}
* no remedy
log_count `"`lgf'"' "drop fixothermargins"
assert r(n) == 0
* detail is an option of seldecomp: it gets no message as an oaxaca habit
log_count `"`lgf'"' "seldecomp's detailed decomposition is"
assert r(n) == 0
* term names: # becomes _, a name that is not a valid name of at most 32 characters or repeats
* an earlier one is term#, with a note
gen byte averyveryverylongname_1 = x1
gen byte averyveryverylongname_2 = x2
gen byte x1_x2 = x4 == "east"
log using `"`lgf'"', text replace name(flog)
seldecomp y (a: x4 averyveryverylongname_1#averyveryverylongname_2) (b: x1_x2 c.x3#x1), ///
    by(group) sequentially detail contvars(x3 2)
seldecomp y (a: x1_x2) (b: x1#x2), by(group) sequentially fixothermargins detail
log close flog
matrix `B' = e(b)
bcol v `B' detail:x1_x2
bcol v `B' detail:term2
log_count `"`lgf'"' "note: term 2 (averyveryverylongname_1#averyveryverylongname_2) is shown as"
assert r(n) == 1
log_count `"`lgf'"' "note: term 2 (x1#x2) is shown as term2 in detail"
assert r(n) == 1
drop averyveryverylongname_1 averyveryverylongname_2 x1_x2


* ==============================================================================
* 2. the layout: without detail nothing changes, with it one equation more
* ==============================================================================
di as txt _n "---- 2. layout"
* (fixothermargins, and the default of the other blocks free, alone and with the options)
foreach opt in "fixothermargins" "" "fixothermargins reverse" "fixothermargins outcomes(ybin)" "fixothermargins ratio(d)" {
    qui seldecomp y `spec' [pw=wt], by(group) `cv' sequentially `opt'
    matrix `B1' = e(b)
    local k1 = e(k_eq)
    assert "`e(detail)'" == ""
    capture confirm matrix e(bystage_detail)
    assert _rc
    assert !strpos("`: coleq e(b)'", "detail")
    e_dump `"`e1f'"', skip(b k_eq)
    qui seldecomp y `spec' [pw=wt], by(group) `cv' sequentially `opt' detail
    matrix `B2' = e(b)
    local nout = cond(strpos("`opt'", "outcomes(ybin)"), 2, 1)
    assert e(k_eq) == `k1' + `nout'
    assert colsof(`B2') == colsof(`B1') + 4 * `nout'
    assert rowsof(e(raking_checks_detail)) == 4 & colsof(e(raking_checks_detail)) == 5
    e_dump `"`e2f'"', skip(b k_eq detail bystage_detail bystage_detail2) skipdetail
    same_file `"`e1f'"' `"`e2f'"'
    * every element of the run without detail, with the same value (exactly)
    local n1 : colfullnames `B1'
    foreach fn of local n1 {
        bcol v1 `B1' `fn'
        bcol v2 `B2' `fn'
        assert `v1' == `v2'
    }
    * the detail equation follows the forward equations (and comes before the reverse ones);
    * without reverse, every existing column keeps its position
    local eqs : coleq `B2'
    local eqs : list uniq eqs
    if strpos("`opt'", "reverse") assert "`eqs'" == "overall marg_expl cum_expl cum_unexpl raked_focal detail rev_marg_expl rev_cum_expl rev_cum_unexpl rev_raked_focal"
    else if strpos("`opt'", "outcomes(ybin)") assert "`eqs'" == "y y_marg_expl y_cum_expl y_cum_unexpl y_raked_focal y_detail ybin ybin_marg_expl ybin_cum_expl ybin_cum_unexpl ybin_raked_focal ybin_detail"
    else {
        assert "`eqs'" == "overall marg_expl cum_expl cum_unexpl raked_focal detail"
        forvalues j = 1/`=colsof(`B1')' {
            assert `B1'[1, `j'] == `B2'[1, `j']
        }
    }
    if strpos("`opt'", "outcomes(ybin)") {
        sum_check y_
        sum_check ybin_
        confirm matrix e(bystage_detail2)
    }
    else sum_check
    assert "`: colnames e(bystage_detail)'" == "`: colnames e(bystage)'"
    assert "`: rownames e(bystage_detail)'" == "x1 x2 x3 x4"
}


* ==============================================================================
* 3. each term against the run with its block split
* ==============================================================================
di as txt _n "---- 3. against the split runs"
foreach var in fixothermargins "" {
    foreach w in "" "[pw=wt]" {
        qui seldecomp y `spec' `w', by(group) `cv' sequentially `var' detail
        matrix `B' = e(b)
        matrix `D' = e(bystage_detail)
        matrix `S' = e(bystage)
        sum_check
        * the whole split: every term a block; its stage matrix is the term chain
        qui seldecomp y `split' `w', by(group) `cv' sequentially `var'
        matrix `B1' = e(b)
        matrix `S1' = e(bystage)
        mat_close `D' `S1' 0
        foreach t in x1 x2 x3 x4 {
            local bt = cond("`t'" == "x1", "a", cond("`t'" == "x4", "c", "`t'"))
            bcol v `B' detail:`t'
            bcol v1 `B1' marg_expl:`bt'
            assert_close `v' `v1' `same' "detail:`t' vs split run"
        }
        * the block split up to x2 only: x2's increment is that run's
        qui seldecomp y (a: x1) x2 (b2: c.x3) (c: x4) `w', by(group) `cv' sequentially `var'
        matrix `B2' = e(b)
        bcol v1 `B2' marg_expl:x2
        bcol v `B' detail:x2
        assert_close `v' `v1' `same' "detail:x2 vs the run split up to x2"
        * the row of a block's last term is the block's stage (the same raking): every column but
        * the increment and its percent, exactly
        local r 0
        foreach e in 1 3 4 {
            local ++r
            forvalues c = 1/30 {
                if inlist(`c', 5, 27) continue
                if missing(`D'[`e', `c']) & missing(`S'[`r', `c']) continue
                assert `D'[`e', `c'] == `S'[`r', `c']
            }
        }
    }
}
* a block of three terms; a single-term block's element equals its block's increment
qui seldecomp y `spec3' [pw=wt], by(group) `cv' sequentially fixothermargins detail
sum_check
matrix `B' = e(b)
qui seldecomp y (a: x1) x2 c.x3 x4 [pw=wt], by(group) `cv' sequentially fixothermargins
matrix `B1' = e(b)
foreach t in x2 x3 x4 {
    bcol v `B' detail:`t'
    bcol v1 `B1' marg_expl:`t'
    assert_close `v' `v1' `same' "detail:`t' (block of three) vs split run"
}
* every block one term: the chain is the block stages
qui seldecomp y x1 x2 c.x3 (c: x4) [pw=wt], by(group) `cv' sequentially fixothermargins detail
matrix `S' = e(bystage)
matrix `D' = e(bystage_detail)
mat_close `S' `D' 0
assert "`: rownames e(bystage_detail)'" == "x1 x2 x3 x4"
sum_check


* ==============================================================================
* 4. the other options, per outcome
* ==============================================================================
di as txt _n "---- 4. options"
* outcomes(): each outcome as its own run
qui seldecomp y `spec' [pw=wt], by(group) `cv' sequentially fixothermargins detail outcomes(ybin)
matrix `B' = e(b)
matrix `D' = e(bystage_detail2)
qui seldecomp ybin `spec' [pw=wt], by(group) `cv' sequentially fixothermargins detail
matrix `B1' = e(b)
matrix `D1' = e(bystage_detail)
mat_close `D' `D1' `same'
foreach t in x1 x2 x3 x4 {
    bcol v `B' ybin_detail:`t'
    bcol v1 `B1' detail:`t'
    assert_close `v' `v1' `same' "outcomes(): ybin detail:`t'"
}
* ratio(): the term stages' raked ratios from the weights of the split run
qui seldecomp y `spec' [pw=wt], by(group) `cv' sequentially fixothermargins detail ratio(d)
matrix `D' = e(bystage_detail)
sum_check
qui seldecomp y `split' [pw=wt], by(group) `cv' sequentially fixothermargins ratio(d) generate(rr)
matrix `S1' = e(bystage)
mat_close `D' `S1' 0
qui sum y [iw=rr_s2] if group == 0
local n1 = r(mean)
qui sum d [iw=rr_s2] if group == 0
local r1 = `D'[2, 1]
local r2 = `n1' / r(mean)
assert_close `r1' `r2' `mtol' "ratio(): raked ratio at term x2"
drop rr rr_s*
* generate(): the block stages' weights, as without detail
qui seldecomp y `spec' [pw=wt], by(group) `cv' sequentially fixothermargins generate(g0)
qui seldecomp y `spec' [pw=wt], by(group) `cv' sequentially fixothermargins detail generate(g1)
foreach s in "" _s1 _s2 {
    assert g0`s' == g1`s'
}
capture confirm variable g1_s3
assert _rc
drop g0* g1*
* insupport, and the change of support: the notes are the block stages' (one note, not two)
qui seldecomp y (a: x1z x2) (b: c.x3 x4) [pw=wt], by(group) `cv' sequentially fixothermargins
e_dump `"`e1f'"', skip(b k_eq)
local nn = e(N_notes)
qui seldecomp y (a: x1z x2) (b: c.x3 x4) [pw=wt], by(group) `cv' sequentially fixothermargins detail
e_dump `"`e2f'"', skip(b k_eq detail bystage_detail) skipdetail
same_file `"`e1f'"' `"`e2f'"'
assert e(N_notes) == `nn'
matrix `D' = e(bystage_detail)
assert `D'[1, 21] > 0 & `D'[1, 21] == `D'[2, 21]
sum_check
qui seldecomp y (a: x1z x2) (b: c.x3 x4) [pw=wt], by(group) `cv' sequentially fixothermargins detail insupport
assert e(N_focal_zero) == 0
matrix `D' = e(bystage_detail)
assert `D'[1, 21] == 0
sum_check
* reverse: the detail is the forward order's; the reverse elements are those of the run without it
qui seldecomp y `spec' [pw=wt], by(group) `cv' sequentially fixothermargins reverse
matrix `B1' = e(b)
qui seldecomp y `spec' [pw=wt], by(group) `cv' sequentially fixothermargins reverse detail
matrix `B2' = e(b)
foreach bn in a b c {
    bcol v1 `B1' rev_marg_expl:`bn'
    bcol v2 `B2' rev_marg_expl:`bn'
    assert `v1' == `v2'
}
sum_check
* makeplot: the blocks, as without detail
qui seldecomp y `spec' [pw=wt], by(group) `cv' sequentially fixothermargins detail makeplot(marginal, twoway(name(dp, replace)))
assert `"`.dp.yaxis1.major.dlg_ticklabel[3]'"' == "a"
assert `"`.dp.yaxis1.major.dlg_ticklabel[5]'"' == "c"
graph drop dp


* ==============================================================================
* 5. bayesboot()
* ==============================================================================
di as txt _n "---- 5. bayesboot()"
* bayesboot(): the detail SEs equal those of marg_expl:<term> in the split run, the
* other SEs (and covariances) those of the run without detail, with the same seed
qui seldecomp y `spec' [pw=wt], by(group) `cv' sequentially fixothermargins detail ///
    bayesboot(reps(20) rseed(11) saving(`"`bsf'"', replace) nodots)
matrix `B' = e(b)
matrix `V' = e(V)
local bsn `e(bsnames)'
assert wordcount("`bsn'") == 22
assert e(N_misreps) == 0
qui seldecomp y `split' [pw=wt], by(group) `cv' sequentially fixothermargins ///
    bayesboot(reps(20) rseed(11) saving(`"`bsf2'"', replace) nodots)
matrix `V1' = e(V)
local bsn2 `e(bsnames)'
foreach t in x1 x2 x3 x4 {
    local bt = cond("`t'" == "x1", "a", cond("`t'" == "x4", "c", "`t'"))
    foreach u in x1 x2 x3 x4 {
        local bu = cond("`u'" == "x1", "a", cond("`u'" == "x4", "c", "`u'"))
        vcell v `V' detail:`t' detail:`u'
        vcell v1 `V1' marg_expl:`bt' marg_expl:`bu'
        assert_close `v' `v1' `same' "cov(detail:`t', detail:`u') vs split run"
    }
}
* the saving() files: the detail columns are the split run's explained columns, row by row
preserve
use `"`bsf'"', clear
mkmat _all, matrix(`C')
use `"`bsf2'"', clear
mkmat _all, matrix(`C1')
restore
local j 0
foreach fn of local bsn {
    local ++j
    if substr("`fn'", 1, 7) != "detail:" continue
    local t = substr("`fn'", 8, .)
    local bt = cond("`t'" == "x1", "a", cond("`t'" == "x4", "c", "`t'"))
    local j1 : list posof "marg_expl:`bt'" in bsn2
    assert `j1' > 0
    forvalues r = 1/20 {
        assert `C'[`r', `j'] == `C1'[`r', `j1']
    }
}
* the rest of e(V): as without detail (the block rows from the rows of the blocks' last terms)
qui seldecomp y `spec' [pw=wt], by(group) `cv' sequentially fixothermargins bayesboot(reps(20) rseed(11) nodots)
matrix `V2' = e(V)
local n2 : colfullnames `V2'
foreach f1 of local n2 {
    foreach f2 of local n2 {
        vcell v `V' `f1' `f2'
        vcell v2 `V2' `f1' `f2'
        if `v2' == 0 assert `v' == 0
        else assert_close `v' `v2' `same' "cov(`f1', `f2') with and without detail"
    }
}
* every saving() row is a plain detail run on that replicate's generate() weights (categorical
* terms: the plain run would cut a c. variable at its own quantiles)
qui seldecomp y (a: x1) (b: x2 x4) [pw=wt], by(group) sequentially fixothermargins detail ///
    bayesboot(reps(3) rseed(5) generate(bw) saving(`"`bsf'"', replace) nodots)
local bsn `e(bsnames)'
assert wordcount("`bsn'") == 17
preserve
use `"`bsf'"', clear
mkmat _all, matrix(`C')
restore
forvalues r = 1/3 {
    qui seldecomp y (a: x1) (b: x2 x4) [iw=bw_`r'], by(group) sequentially fixothermargins detail
    matrix `B1' = e(b)
    local j 0
    foreach fn of local bsn {
        local ++j
        bcol v `B1' `fn'
        assert_close `C'[`r',`j'] `v' 1e-10 "saving() row `r', `fn'"
    }
}
drop bw_*
* bca with clusters: the acceleration of a detail element is the split run's
qui seldecomp y `spec' [pw=wt], by(group) `cv' sequentially fixothermargins detail ///
    bayesboot(reps(10) rseed(3) cluster(cluster) bca nodots)
matrix `A' = e(accel)
bcol a2 `A' detail:x2
bcol a3 `A' detail:x3
qui seldecomp y `split' [pw=wt], by(group) `cv' sequentially fixothermargins ///
    bayesboot(reps(10) rseed(3) cluster(cluster) bca nodots)
matrix `A1' = e(accel)
bcol b2 `A1' marg_expl:x2
bcol b3 `A1' marg_expl:x3
assert_close `a2' `b2' `same' "accel detail:x2"
assert_close `a3' `b3' `same' "accel detail:x3"
* reverse and outcomes() with bayesboot(): the reverse rows as without detail
qui seldecomp y `spec' [pw=wt], by(group) `cv' sequentially fixothermargins reverse outcomes(ybin) ///
    bayesboot(reps(10) rseed(4) nodots)
matrix `SE1' = e(se)
qui seldecomp y `spec' [pw=wt], by(group) `cv' sequentially fixothermargins reverse outcomes(ybin) detail ///
    bayesboot(reps(10) rseed(4) nodots)
matrix `SE' = e(se)
foreach fn in `: colfullnames `SE1'' {
    bcol v1 `SE1' `fn'
    bcol v `SE' `fn'
    assert_close `v' `v1' `same' "se(`fn') with reverse outcomes()"
}
bcol v `SE' ybin_detail:x3
assert `v' > 0
* contvars(..., recut): the detail SEs are the split run's
qui seldecomp y `spec' [pw=wt], by(group) contvars(x3 3, recut) sequentially fixothermargins detail ///
    bayesboot(reps(10) rseed(6) nodots)
matrix `SE' = e(se)
qui seldecomp y `split' [pw=wt], by(group) contvars(x3 3, recut) sequentially fixothermargins ///
    bayesboot(reps(10) rseed(6) nodots)
matrix `SE1' = e(se)
bcol v `SE' detail:x3
bcol v1 `SE1' marg_expl:x3
assert_close `v' `v1' `same' "recut: se(detail:x3)"


* ==============================================================================
* 6. prefixes
* ==============================================================================
di as txt _n "---- 6. prefixes"
preserve
keep if cluster <= 20
svyset cluster [pw=wt]
svy jackknife, nodots: seldecomp y `spec', by(group) `cv' sequentially fixothermargins detail
matrix `B2' = e(b)
matrix `V' = e(V)
foreach t in x1 x2 x3 x4 {
    vcell v `V' detail:`t' detail:`t'
    assert `v' > 0
}
qui seldecomp y `spec' [pw=wt], by(group) `cv' sequentially fixothermargins detail
matrix `B1' = e(b)
mat_close `B1' `B2' `same'
* the jackknife SEs of the detail elements equal those of the split run's explained elements
svy jackknife, nodots: seldecomp y `split', by(group) `cv' sequentially fixothermargins
matrix `V1' = e(V)
vcell v `V' detail:x2 detail:x2
vcell v1 `V1' marg_expl:x2 marg_expl:x2
assert_close `v' `v1' 1e-10 "svy jackknife: var(detail:x2)"
svyset, clear
restore
bootstrap, reps(3) seed(2) nodots: seldecomp y `spec', by(group) `cv' sequentially fixothermargins detail
matrix `V' = e(V)
vcell v `V' detail:x3 detail:x3
assert `v' > 0
* by: each group's results are those of the if-restricted run
bysort x2: seldecomp y (a: x1) (b: c.x3 x4) [pw=wt], by(group) `cv' sequentially fixothermargins detail
matrix `B2' = e(b)
qui seldecomp y (a: x1) (b: c.x3 x4) [pw=wt] if x2 == 1, by(group) `cv' sequentially fixothermargins detail
matrix `B1' = e(b)
mat_close `B1' `B2' `same'


* ==============================================================================
* 7. display
* ==============================================================================
di as txt _n "---- 7. display"
local ls0 = c(linesize)
set linesize 160
log using `"`lgf'"', text replace name(flog)
di "SDW_W_BEGIN"
seldecomp y `spec3' [pw=wt], by(group) `cv' sequentially fixothermargins detail reverse bayesboot(reps(10) rseed(2) nodots citable)
di "SDW_W_END"
log close flog
set linesize `ls0'
log_width `"`lgf'"'
di as txt "widest display line: " r(max)
assert r(max) <= 78
foreach t in "  and term by term within each block (detail)" ///
    "Indented: the block's terms, added in turn (Cumulative: within the block)." {
    log_count `"`lgf'"' "`t'"
    if r(n) != 1 {
        di as err `"display: "`t'" not found once"'
        exit 9
    }
}
* the detail equation has a descriptive title in both inference tables
log_count `"`lgf'"' "Marginal explained, by term"
assert r(n) == 2
* the forward stage table: the block of three terms followed by its terms, indented, with
* raked_focal, the explained part within the block, the increment and its percent; the reverse table
* has no term rows
stage_lines `"`lgf'"' "Decomposition by raking stage (sequentially, fixothermargins)"
local fl `"`r(lines)'"'
assert r(n) == 11
matrix `D' = e(bystage_detail)
matrix `S' = e(bystage)
local k 0
foreach t in x2 x3 x4 {
    local ++k
    local r = `k' + 1
    mata: st_local("ln", ustrsplit(st_local("fl"), "[|]")[7 + `k'])
    assert substr(`"`ln'"', 1, 8) == "        " & substr(`"`ln'"', 9, strlen("`t'")) == "`t'"
    local v1 = strtrim(string(`D'[`r', 1], "%9.0g"))
    local v2 = strtrim(string(`D'[`r', 2] - `S'[1, 2], "%9.0g"))
    local v5 = strtrim(string(`D'[`r', 5], "%9.0g"))
    local v27 = strtrim(string(`D'[`r', 27], "%6.1f"))
    local want `"`v1' `v2' `v5' `v27'"'
    local got = stritrim(strtrim(substr(`"`ln'"', 9 + strlen("`t'"), .)))
    if `"`got'"' != `"`want'"' {
        di as err `"term row `t': `got' (expected `want')"'
        exit 9
    }
}
* the last term's explained part within the block is the block's increment
local w1 = `D'[4, 2] - `S'[1, 2]
local w2 = `S'[2, 5]
assert `w1' == `w2'
stage_lines `"`lgf'"' "Decomposition by raking stage (sequentially, fixothermargins), reverse order"
assert r(n) == 8
* the rows rebuilt from e(b) when the stage matrices are not posted: the same lines
log using `"`lgf'"', text replace name(flog)
seldecomp
log close flog
stage_lines `"`lgf'"' "Decomposition by raking stage (sequentially, fixothermargins)"
local tl1 `"`r(lines)'"'
drop_mats
log using `"`lgf'"', text replace name(flog)
seldecomp
log close flog
stage_lines `"`lgf'"' "Decomposition by raking stage (sequentially, fixothermargins)"
assert `"`r(lines)'"' == `"`tl1'"'
* single-term blocks have no term rows, and then no footnote line; outcomes(): each outcome's
* table has its term rows, the footnote once
log using `"`lgf'"', text replace name(flog)
seldecomp y x1 x2 (c: x4) [pw=wt], by(group) sequentially fixothermargins detail
seldecomp y `spec' [pw=wt], by(group) `cv' sequentially fixothermargins detail outcomes(ybin)
log close flog
log_count `"`lgf'"' "Indented: the block's terms"
assert r(n) == 1
log_count `"`lgf'"' "        x3   "
assert r(n) == 2
log_count `"`lgf'"' "Decomposition by raking stage (sequentially, fixothermargins): ybin"
assert r(n) == 1
* one block of two terms: its stage table (one stage) is shown with the term rows; without
* detail, or with one term, there is no stage table for one stage
log using `"`lgf'"', text replace name(flog)
seldecomp y (ab: x1 x2) [pw=wt], by(group) sequentially fixothermargins detail
log close flog
stage_lines `"`lgf'"' "Decomposition by raking stage (sequentially, fixothermargins)"
assert r(n) == 9
log_count `"`lgf'"' "        x2   "
assert r(n) == 1
log using `"`lgf'"', text replace name(flog)
seldecomp y (ab: x1 x2) [pw=wt], by(group) sequentially fixothermargins
seldecomp y x1 [pw=wt], by(group) sequentially fixothermargins detail
log close flog
log_count `"`lgf'"' "Decomposition by raking stage"
assert r(n) == 0


* ==============================================================================
* 8. a term stage that cannot meet its margins
* ==============================================================================
di as txt _n "---- 8. term stages that fail"
* (ab: ta tb) with tb = 1 - ta: one raking as a block; term by term under fixothermargins, ta
* raked to the Reference Group's share with tb held at the Focal Group's cannot be met; without fixothermargins tb
* is left free, and a term stage rakes a subset of its block stage's margins
qui seldecomp y (a: x1) (ab: ta tb) [pw=wt], by(group) sequentially fixothermargins
rcof "noi seldecomp y (a: x1) (ab: ta tb) [pw=wt], by(group) sequentially fixothermargins detail" == 430
qui seldecomp y (a: x1) (ab: ta tb) [pw=wt], by(group) sequentially detail
assert e(targets_met) == 1
sum_check
* (the error at linesize 160, so that it is on one line of the log)
local ls0 = c(linesize)
log using `"`lgf'"', text replace name(flog)
set linesize 160
capture noisily seldecomp y (a: x1) (ab: ta tb) [pw=wt], by(group) sequentially fixothermargins detail
set linesize `ls0'
seldecomp y (a: x1) (ab: ta tb) [pw=wt], by(group) sequentially fixothermargins detail force
log close flog
* one message: the term, the block, the worst margin, and force
log_count `"`lgf'"' "detail: adding ta (block ab) misses 1 raking margin (worst: ta = 0); specify force to keep the results"
assert r(n) == 1
* force is the only remedy named: neither dropping fixothermargins nor reordering the block's
* terms (reordering cannot help when a later term is a function of the term: tb = 1 - ta fails
* in either order; it can when the tie is partial: tz = 1 only where ta = 1, one order runs,
* the other is 430)
log_count `"`lgf'"' "drop fixothermargins"
assert r(n) == 0
log_count `"`lgf'"' "reorder the block's"
assert r(n) == 0
log_count `"`lgf'"' "raking did not reach its targets"
assert r(n) == 0
* the kept result reports the failed detail stage in the diagnostic block
log_count `"`lgf'"' "Raking failure counts"
assert r(n) == 1
log_count `"`lgf'"' "Worst missed target for each affected variable"
assert r(n) == 1
log_count `"`lgf'"' "detail: stage 2" unwrap
assert r(n) == 1
assert e(converged) == 0 & e(targets_met) == 1
matrix `D' = e(bystage_detail)
assert `D'[2, 6] == 0 & `D'[2, 7] > 0 & `D'[3, 6] == 1
* the tied pair fails in either order; a partial tie can run in one
rcof "noi seldecomp y (a: x1) (ab: tb ta) [pw=wt], by(group) sequentially fixothermargins detail" == 430
gen byte tz = ta == 1 & mod(id, 3) != 0
qui seldecomp y (a: x1) (ab: ta tz) [pw=wt], by(group) sequentially fixothermargins detail
sum_check
rcof "noi seldecomp y (a: x1) (ab: tz ta) [pw=wt], by(group) sequentially fixothermargins detail" == 430
drop tz
* with a forward stage that misses its margins too (tc = ta in its own block, held at the Focal Group's
* share at stage 1): both the block and detail diagnostics show their failures
log using `"`lgf'"', text replace name(flog)
seldecomp y (cc: tc) (ab: ta tb) [pw=wt], by(group) sequentially fixothermargins detail force
log close flog
log_count `"`lgf'"' "targets *not* met"
assert r(n) == 1
log_count `"`lgf'"' "Raking failure counts"
assert r(n) == 2
log_count `"`lgf'"' "detail: stages 1, 2" unwrap
assert r(n) == 1
* a miss at a block's stage only is correctly identified at its matching detail stage
log using `"`lgf'"', text replace name(flog)
seldecomp y (cc: tc) (ab: ta x1) [pw=wt], by(group) sequentially fixothermargins detail force
log close flog
log_count `"`lgf'"' "targets *not* met"
assert r(n) == 1
log_count `"`lgf'"' "detail: stage 1" unwrap
assert r(n) == 1
matrix `D' = e(bystage_detail)
assert `D'[1, 7] > 0 & `D'[2, 7] == 0 & `D'[3, 7] == 0
* trimming that leaves margins unmet at a term stage (tb3 = 1 - ta but for 1 observation in 23):
* rc 430, the message naming the bounds, and kept by force with a note, as at a block's stage.
* The block's stage meets its margins, so the run without detail needs no force and has no
* note; the split run's stage 2 is the same raking as the term stage, with the block stages'
* note
gen byte tb3 = 1 - x2
replace tb3 = x2 if mod(id, 23) == 0
gen byte te3 = 1 - x2
replace te3 = x2 if mod(id, 19) == 0
local trn "detail, adding ta (block ab): trimming left 3 raking margins unmet (worst: ta = 0; 17 weights at a bound); results kept because force was specified"
qui seldecomp y (a: x1) (ab: ta tb3) [pw=wt], by(group) sequentially fixothermargins trimhirel(3)
assert e(N_notes) == 0 & e(targets_met) == 1
rcof "noi seldecomp y (a: x1) ta tb3 [pw=wt], by(group) sequentially fixothermargins trimhirel(3)" == 430
qui seldecomp y (a: x1) ta tb3 [pw=wt], by(group) sequentially fixothermargins trimhirel(3) force
assert e(N_notes) == 1 & e(targets_met) == 0
assert strpos(`"`e(note1)'"', "stage 2: trimming left 3 raking margins unmet (worst: ta = 0; 17 weights") == 1
matrix `B1' = e(b)
bcol v1 `B1' marg_expl:ta
local ls0 = c(linesize)
* (the error at linesize 160, so that it is on one line of the log): one message, the term,
* the block, the worst margin, the bounds and force
log using `"`lgf'"', text replace name(flog)
set linesize 160
capture noisily seldecomp y (a: x1) (ab: ta tb3) [pw=wt], by(group) sequentially fixothermargins detail trimhirel(3)
local rc = _rc
set linesize `ls0'
log close flog
assert `rc' == 430
log_count `"`lgf'"' "detail: adding ta (block ab) misses 3 raking margins (worst: ta = 0); the trimming bounds are binding: relax them, or specify force to keep the results"
assert r(n) == 1
set linesize 160
log using `"`lgf'"', text replace name(flog)
di "SDW_W_BEGIN"
seldecomp y (a: x1) (ab: ta tb3) [pw=wt], by(group) sequentially fixothermargins detail trimhirel(3) force
seldecomp
di "SDW_W_END"
log close flog
set linesize `ls0'
assert e(N_notes) == 1
assert `"`e(note1)'"' == "`trn'"
assert e(converged) == 0 & e(targets_met) == 1 & e(badcontrols) == 0
matrix `D' = e(bystage_detail)
assert `D'[2, 7] == 3 & `D'[2, 6] == 0 & `D'[2, 22] == 17 & `D'[3, 7] == 0 & `D'[1, 7] == 0
matrix `B' = e(b)
bcol v `B' detail:ta
assert_close `v' `v1' `same' "detail:ta vs the split run's marg_expl:ta, trimmed"
sum_check
* the failed detail stage triggers the same diagnostics at estimation and replay
log_count `"`lgf'"' "Raking failure counts"
assert r(n) == 2
log_count `"`lgf'"' "Worst missed target for each affected variable"
assert r(n) == 2
log_count `"`lgf'"' "detail: stage 2" unwrap
assert r(n) == 2
* no remedy after the note
log_count `"`lgf'"' "widen the bounds"
assert r(n) == 0
log_width `"`lgf'"'
assert r(max) <= 78
* two term stages of a block: one note, naming the last (as the error does)
rcof "noi seldecomp y (a: x1) (ab: ta tb3 te3) [pw=wt], by(group) sequentially fixothermargins detail trimhirel(3)" == 430
qui seldecomp y (a: x1) (ab: ta tb3 te3) [pw=wt], by(group) sequentially fixothermargins detail trimhirel(3) force
assert e(N_notes) == 1 & e(targets_met) == 1
matrix `D' = e(bystage_detail)
assert `D'[2, 7] > 0 & `D'[3, 7] > 0 & `D'[4, 7] == 0
assert `D'[2, 6] == 0 & `D'[3, 6] == 0 & e(converged) == 0
assert strpos(`"`e(note1)'"', "detail, adding tb3 (block ab): trimming left " + string(`D'[3, 7]) + " raking " + plural(`D'[3, 7], "margin") + " unmet (worst: ") == 1
local tl = "; " + string(`D'[3, 22]) + " " + plural(`D'[3, 22], "weight") + " at a bound); results kept because force was specified"
assert substr(`"`e(note1)'"', -strlen("`tl'"), .) == "`tl'"
sum_check
* the block's own stage has the trimming note, kept by force: it says "also in the detail",
* and there is no note more for the block's terms
rcof "noi seldecomp y (a: x1) (ab: ta tb3 te3) [pw=wt], by(group) sequentially fixothermargins detail trimhirel(2)" == 430
qui seldecomp y (a: x1) (ab: ta tb3 te3) [pw=wt], by(group) sequentially fixothermargins detail trimhirel(2) force
matrix `S' = e(bystage)
matrix `D' = e(bystage_detail)
assert `S'[2, 7] > 0 & `D'[2, 7] > 0 & `D'[3, 7] > 0 & e(targets_met) == 0
local nn = e(N_notes)
local nd 0
forvalues m = 1/`nn' {
    if strpos(`"`e(note`m')'"', "detail") {
        local ++nd
        assert strpos(`"`e(note`m')'"', ": trimming left ") > 0
        assert strpos(`"`e(note`m')'"', "; also in the detail, adding tb3); results kept because force was specified") > 0
    }
}
assert `nd' == 1
qui seldecomp y (a: x1) (ab: ta tb3 te3) [pw=wt], by(group) sequentially fixothermargins trimhirel(2) force
assert e(N_notes) == `nn'
* trimming that changes no weight at the term stage: the 430 does not name the bounds
log using `"`lgf'"', text replace name(flog)
capture noisily seldecomp y (a: x1) (ab: ta tb) [pw=wt], by(group) sequentially fixothermargins detail trimhirel(3)
local rc = _rc
log close flog
assert `rc' == 430
log_count `"`lgf'"' "specify force to keep the results" unwrap
assert r(n) == 1
log_count `"`lgf'"' "trimming bounds are binding" unwrap
assert r(n) == 0
drop tb3 te3
* bayesboot(): force keeps the replicates that miss a margin (every one here: the tie holds
* on any weights), with a note, and counts them in e(N_forced)
log using `"`lgf'"', text replace name(flog)
seldecomp y (a: x1) (ab: ta tb) [pw=wt], by(group) sequentially fixothermargins detail force bayesboot(reps(3) rseed(1) nodots)
log close flog
assert e(N_reps) == 3 & e(N_misreps) == 0 & e(N_forced) == 3
log_count `"`lgf'"' "force kept 3 of 3 replications that missed raking targets"
assert r(n) == 1

di as res _n "test_detail.do: all tests passed"

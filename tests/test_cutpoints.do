* test_cutpoints.do -- the contvars() cut points in replicates:
* by default every replicate, of bayesboot() and of every replication prefix, uses the cut
* points of the full-sample run; contvars(..., recut) recomputes them in each replicate by the
* same rule.
*
* Under a prefix the replicate calls are recognised by the word _loop in c(prefix), which
* Stata's replicate loops set under every prefix; the
* full-sample call records its cut points for them. Checked here: every prefix (svy brr with
* replicate weights and with a Hadamard matrix, svy jackknife by PSU and with replicate
* weights, svy sdr, svy bootstrap, bootstrap, bootstrap with bca, jackknife, and on Stata 19
* the bayesboot prefix) gives exactly the replicates, e(b) and e(V) of the same run on the
* full-sample bins typed as a categorical variable (1e-12); recut reproduces the replicates
* of plain runs on the replicate weights or samples, which recompute the cut points; a call
* without a prefix is never affected; the record's rules (set prefix in a program stands in
* for a prefix); bayesboot(): the saved replicates equal plain runs on the replicate weights,
* with the full-sample bins by default and with c. and contvars() under recut (also
* interactions, the contvars() suboptions focalgroup and refgroup, reverse), and the BCa
* acceleration equals a jackknife made by hand. The staged runs use sequentially
* fixothermargins, with bayesboot() replicating every stage.
*
* Run from the tests directory:   . do test_cutpoints.do

version 16
clear
set more off
set varabbrev off
discard             // drop any seldecomp already in memory; the repository copy is what gets tested

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

* ---- helpers ------------------------------------------------------------------
capture program drop file_same
program define file_same
    // file_same f1 f2 tol  --  two saving() files of a prefix have the same variables, rows
    // and values (missing in the same cells), to tol
    args f1 f2 tol
    preserve
    qui use `"`f1'"', clear
    qui ds, has(type numeric)
    local v1 `r(varlist)'
    mkmat `v1', matrix(F1)
    qui use `"`f2'"', clear
    qui ds, has(type numeric)
    local v2 `r(varlist)'
    mkmat `v2', matrix(F2)
    restore
    assert "`v1'" == "`v2'"
    mata: st_local("ok", strofreal(tc_same(st_matrix("F1"), st_matrix("F2"), `tol')))
    if !`ok' {
        di as err "file_same: `f1' and `f2' differ"
        exit 9
    }
end
capture program drop row_is_b
program define row_is_b
    // row_is_b M r tol  --  row r of M (the columns of e(b), in order) equals the last run's e(b)
    args M r tol
    tempname b R
    matrix `b' = e(b)
    matrix `R' = `M'[`r', 1...]
    assert colsof(`R') == colsof(`b')
    mata: st_local("ok", strofreal(tc_same(st_matrix("`R'"), st_matrix("`b'"), `tol')))
    if !`ok' {
        di as err "row_is_b: row `r' of `M' differs from e(b)"
        matrix list `R', format(%20.0g)
        matrix list `b', format(%20.0g)
        exit 9
    }
end
capture program drop draws_are_b
program define draws_are_b
    // draws_are_b D r tol  --  row r of e(bsdraws) saved in D (columns named as e(b)) equals
    // the matching columns of the last run's e(b)
    args D r tol
    tempname b
    matrix `b' = e(b)
    local names : colfullnames `D'
    local j 0
    foreach nm of local names {
        local ++j
        local k = colnumb(`b', "`nm'")
        assert `k' < .
        if reldif(`D'[`r', `j'], `b'[1, `k']) > `tol' {
            di as err "replicate `r', `nm': " %20.0g `D'[`r', `j'] " vs a plain run " %20.0g `b'[1, `k']
            exit 9
        }
    }
end
capture program drop fake
program define fake, eclass
    // fake <words> : <command>  --  runs <command> with <words> added to c(prefix), as a
    // prefix command does (a full-sample call: bootstrap; a replicate: bootstrap _loop _loop_bs)
    gettoken words 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    set prefix `words'
    `0'
end
capture program drop cutcache_clear
program define cutcache_clear
    mata: rmexternal("_seldecomp_cutcache")
end
capture mata: mata drop tc_same()
mata:
real scalar tc_same(real matrix A, real matrix B, real scalar tol)
{
    real colvector a, b
    if (rows(A) != rows(B) | cols(A) != cols(B)) return(0)
    a = vec(A)
    b = vec(B)
    if (sum((a :>= .) :!= (b :>= .)) > 0) return(0)
    a = select(a, a :< .)
    b = select(b, b :< .)
    if (rows(a) == 0) return(1)
    return(max(abs(a - b) :/ (abs(b) :+ 1)) <= tol)
}
// the jackknife acceleration of each column of JK (rows = deletions), as in test_bootstrap.do
real rowvector tc_accel(real matrix JK)
{
    real rowvector a
    real colvector v, d
    real scalar j
    a = J(1, cols(JK), .)
    for (j = 1; j <= cols(JK); j++) {
        v = select(JK[., j], JK[., j] :< .)
        d = mean(v) :- v
        a[j] = sum(d :^ 3) / (6 * sum(d :^ 2)^1.5)
    }
    return(a)
}
end

* ---- data: 16 strata x 2 PSUs x 12 observations, a continuous raking variable xc -------
set seed 20260926
set obs 384
gen int id = _n
gen int stratum = ceil(_n / 24)
gen byte psu = 1 + mod(ceil(_n / 12) - 1, 2)
gen int psuid = stratum * 10 + psu
gen double wt = 50 + 100 * runiform()
gen byte group = runiform() < .5
gen double xc = rnormal() + .6 * group
gen double x2c = runiform() + .3 * group
gen byte x1 = 1 + (runiform() < .4) + (runiform() < .3 + .3 * group)
gen double y = 2 + xc + x1 + .8 * group + rnormal()
mata: Hd = (1); for (i = 1; i <= 4; i++) Hd = (Hd, Hd \ Hd, -Hd); st_matrix("HAD", Hd)
forvalues r = 1/16 {
    gen double brr`r' = 2 * wt * ((HAD[`r', stratum] == 1) == (psu == 1))
}
local rk  (age: c.xc) (base: i.x1)
local opt by(group) sequentially fixothermargins

* the full-sample bins, typed as categorical variables: weighted (the svy prefixes' full-sample
* call has the svyset pweight) and unweighted (bootstrap, jackknife); the pooled quantile rule
* of contvars(), x <= cut j in bin j
seldecomp y `rk' [pw=wt], `opt' contvars(xc 4)
matrix CW = e(cuts1)
matrix BW = e(b)
seldecomp y `rk', `opt' contvars(xc 4)
matrix CU = e(cuts1)
matrix BU = e(b)
gen byte binw = 4
gen byte binu = 4
forvalues j = 3(-1)1 {
    replace binw = `j' if xc <= CW[1, `j']
    replace binu = `j' if xc <= CU[1, `j']
}
seldecomp y (age: i.binw) (base: i.x1) [pw=wt], `opt'
assert mreldif(e(b), BW) < 1e-12
seldecomp y (age: i.binu) (base: i.x1), `opt'
assert mreldif(e(b), BU) < 1e-12
assert mreldif(CW, CU) > 1e-3          // the two sets of bins differ


* ==============================================================================
* 1. every prefix: the replicates use the full-sample cut points
* ==============================================================================
di as txt _n "---- 1. every prefix: the full-sample cut points in every replicate"
capture program drop pfx_same
program define pfx_same
    // pfx_same "prefix, with SAVE where its saving() goes" bins [seed]  --  the run with c.xc
    // and contvars() and the run on the full-sample bins give the same e(b), e(V) and
    // replicates (seed: set before each run, for a prefix that ignores its seed() option)
    args pfx bin seed
    local p1 : subinstr local pfx "SAVE" "saving(pc_c, replace double)"
    local p2 : subinstr local pfx "SAVE" "saving(pc_b, replace double)"
    if "`seed'" != "" set seed `seed'
    `p1': seldecomp y (age: c.xc) (base: i.x1), by(group) sequentially fixothermargins contvars(xc 4)
    tempname B V
    matrix `B' = e(b)
    matrix `V' = e(V)
    if "`seed'" != "" set seed `seed'
    `p2': seldecomp y (age: i.`bin') (base: i.x1), by(group) sequentially fixothermargins
    assert mreldif(e(b), `B') < 1e-12
    assert mreldif(e(V), `V') < 1e-12
    assert el(e(V), 1, 1) > 0
    file_same pc_c.dta pc_b.dta 1e-12
    erase pc_c.dta
    erase pc_b.dta
end
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr16)
pfx_same "svy brr, SAVE nodots" binw
svyset psuid [pw=wt], strata(stratum)
pfx_same "svy brr, hadamard(HAD) SAVE nodots" binw
pfx_same "svy jackknife, SAVE nodots" binw
svyset psuid [pw=wt], strata(stratum) jkrweight(brr1-brr16)
pfx_same "svy jackknife, SAVE nodots" binw
svyset psuid [pw=wt], sdrweight(brr1-brr16)
pfx_same "svy sdr, SAVE nodots" binw
svyset psuid [pw=wt], bsrweight(brr1-brr16)
pfx_same "svy bootstrap, SAVE nodots" binw
svyset, clear
pfx_same "bootstrap, reps(8) seed(11) SAVE nodots" binu
* bca: the bootstrap's full-sample call, then a jackknife's full-sample call and deletions,
* then the bootstrap replicates
pfx_same "bootstrap, reps(8) seed(12) cluster(psuid) bca SAVE nodots" binu
pfx_same "jackknife, cluster(psuid) SAVE nodots" binu
if c(stata_version) >= 19 {
    // (the bayesboot prefix ignores seed(): "replication weights were specified")
    pfx_same "bayesboot, reps(8) SAVE nodots" binu 13
}
* the replicates' own cut points differ: recut changes e(V) (the tests above have teeth)
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr16)
svy brr, nodots: seldecomp y `rk', `opt' contvars(xc 4)
matrix V0 = e(V)
svy brr, nodots: seldecomp y `rk', `opt' contvars(xc 4, recut)
assert mreldif(e(V), V0) > 1e-4


* ==============================================================================
* 2. recut: each replicate recomputes the cut points
* ==============================================================================
di as txt _n "---- 2. recut under the prefixes"
* svy brr: the replicate calls are seldecomp [iw=brr#]; plain runs with those weights
* recompute the cut points from them
svy brr, nodots saving(pc_r, replace double): seldecomp y `rk', `opt' contvars(xc 4, recut)
assert "`e(recut)'" == "recut"             // svy keeps the command's e() macros
preserve
use pc_r, clear
qui ds, has(type numeric)
mkmat `r(varlist)', matrix(SR)
restore
assert rowsof(SR) == 16
forvalues r = 1/16 {
    qui seldecomp y `rk' [iw=brr`r'], `opt' contvars(xc 4)
    row_is_b SR `r' 1e-12
}
* jackknife by cluster: the deletions in the order of psuid
jackknife, cluster(psuid) nodots saving(pc_r, replace double): seldecomp y `rk', `opt' contvars(xc 4, recut)
preserve
use pc_r, clear
qui ds, has(type numeric)
mkmat `r(varlist)', matrix(SR)
restore
qui levelsof psuid, local(pl)
local r 0
foreach c of local pl {
    local ++r
    qui seldecomp y `rk' if psuid != `c', `opt' contvars(xc 4)
    row_is_b SR `r' 1e-12
}
assert `r' == 32
* bootstrap: inside each replicate, the recut call's cut points are the replicate sample's
* own quantiles, and its results those of the replicate on these bins (a failure makes the
* replicate fail)
capture program drop recutchk
program define recutchk, eclass
    syntax [if] [in] [iw]
    local wt
    if "`weight'" != "" local wt [iw`exp']
    seldecomp y (age: c.xc) (base: i.x1) `wt' `if' `in', by(group) sequentially fixothermargins contvars(xc 4, recut)
    tempname B C
    tempvar s bin
    matrix `B' = e(b)
    matrix `C' = e(cuts1)
    qui gen byte `s' = e(sample)
    local aw
    if "`weight'" != "" local aw [aw`exp']
    qui _pctile xc `aw' if `s', nq(4)
    forvalues j = 1/3 {
        if `C'[1, `j'] != r(r`j') exit 9
    }
    qui gen byte `bin' = 4 if `s'
    forvalues j = 3(-1)1 {
        qui replace `bin' = `j' if `s' & xc <= `C'[1, `j']
    }
    seldecomp y (age: i.`bin') (base: i.x1) `wt' `if' `in', by(group) sequentially fixothermargins
    if mreldif(e(b), `B') > 1e-12 exit 9
end
svyset, clear
bootstrap _b, reps(8) seed(14) nodots: recutchk
assert e(N_misreps) == 0 & e(N_reps) == 8
if c(stata_version) >= 19 {
    bayesboot _b, reps(8) nodots: recutchk
    assert e(N_misreps) == 0 & e(N_reps) == 8
}


* ==============================================================================
* 3. a call without a prefix is never affected; the record's rules
* ==============================================================================
di as txt _n "---- 3. calls without a prefix; the record's rules"
* right after a prefix run, on changed data: the cut points are those of the data in memory
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr16)
svy brr, nodots: seldecomp y `rk', `opt' contvars(xc 4)
preserve
replace xc = xc^3 + 2 * (xc > 0)
drop if mod(_n, 5) == 0
seldecomp y `rk' [pw=wt], `opt' contvars(xc 4)
matrix C1 = e(cuts1)
matrix B1 = e(b)
_pctile xc [aw=wt] if e(sample), nq(4)
forvalues j = 1/3 {
    assert C1[1, `j'] == r(r`j')
}
cutcache_clear
seldecomp y `rk' [pw=wt], `opt' contvars(xc 4)
assert mreldif(e(b), B1) == 0 & mreldif(e(cuts1), C1) == 0
restore
svyset, clear
* a call without a prefix drops the record: a replicate after it has no cut points (459)
local fs seldecomp y `rk', `opt' contvars(xc 4)
fake bootstrap : `fs'
seldecomp y `rk', `opt' contvars(xc 4)
rcof `"noi fake bootstrap _loop _loop_bs : `fs'"' == 459
* so does a prefix that replicates nothing (statsby, permute, simulate, frame, collect)
fake bootstrap : `fs'
fake statsby : `fs'
rcof `"noi fake bootstrap _loop _loop_bs : `fs'"' == 459
* a replicate uses the cut points of the full-sample call, whatever its own sample
gen byte s1 = mod(_n, 3) != 0
gen byte s2 = mod(_n, 4) != 0
seldecomp y `rk' if s1, `opt' contvars(xc 4)
matrix C1 = e(cuts1)
seldecomp y `rk' if s2, `opt' contvars(xc 4)
matrix C2 = e(cuts1)
assert mreldif(C1, C2) > 1e-3
fake bootstrap : seldecomp y `rk' if s1, `opt' contvars(xc 4)
fake bootstrap _loop _loop_bs : seldecomp y `rk' if s1, `opt' contvars(xc 4)
assert mreldif(e(cuts1), C1) == 0
fake bootstrap _loop _loop_bs : seldecomp y `rk' if s2, `opt' contvars(xc 4)
assert mreldif(e(cuts1), C1) == 0
* recut: its own
fake bootstrap _loop _loop_bs : seldecomp y `rk' if s2, `opt' contvars(xc 4, recut)
assert mreldif(e(cuts1), C2) == 0
* the first full-sample call after a replicate starts a new record
fake bootstrap : seldecomp y `rk' if s2, `opt' contvars(xc 4)
fake bootstrap _loop _loop_bs : seldecomp y `rk' if s1, `opt' contvars(xc 4)
assert mreldif(e(cuts1), C2) == 0
* two specifications before the replicates (a program calling seldecomp twice): each its own
fake bootstrap : seldecomp y `rk' if s1, `opt' contvars(xc 4)
fake bootstrap : seldecomp y `rk' if s2, `opt' contvars(xc 5)
matrix C5 = e(cuts1)
fake bootstrap _loop _loop_bs : seldecomp y `rk' if s2, `opt' contvars(xc 4)
assert mreldif(e(cuts1), C1) == 0
fake bootstrap _loop _loop_bs : seldecomp y `rk' if s1, `opt' contvars(xc 5)
assert mreldif(e(cuts1), C5) == 0
* the same specification twice on different samples: the replicates cannot tell which applies
fake bootstrap : seldecomp y `rk' if s1, `opt' contvars(xc 4)
fake bootstrap : seldecomp y `rk' if s1, `opt' contvars(xc 4)
fake bootstrap _loop _loop_bs : seldecomp y `rk', `opt' contvars(xc 4)
fake bootstrap : seldecomp y `rk' if s1, `opt' contvars(xc 4)
fake bootstrap : seldecomp y `rk' if s2, `opt' contvars(xc 4)
rcof `"noi fake bootstrap _loop _loop_bs : seldecomp y `rk', `opt' contvars(xc 4)"' == 459
fake bootstrap _loop _loop_bs : seldecomp y `rk', `opt' contvars(xc 4, recut)
* ... under a real prefix, every replicate fails
capture program drop twocalls
program define twocalls, eclass
    syntax [, recut]
    if "`recut'" != "" local recut ", recut"
    seldecomp y (age: c.xc) (base: i.x1) if mod(_n, 3) != 0, by(group) sequentially fixothermargins contvars(xc 4`recut')
    seldecomp y (age: c.xc) (base: i.x1) if mod(_n, 4) != 0, by(group) sequentially fixothermargins contvars(xc 4`recut')
end
capture noisily bootstrap _b, reps(3) seed(1) nodots: twocalls
assert _rc != 0 | e(N_misreps) == 3
bootstrap _b, reps(3) seed(1) nodots: twocalls, recut
assert e(N_misreps) == 0 & e(N_reps) == 3
* the key: the weights and if/in are not part of it (the prefixes change them); contvars()
* and focalgroup() are
fake bootstrap : seldecomp y `rk' [pw=wt] if s1, `opt' contvars(xc 4)
fake bootstrap _loop _loop_bs : seldecomp y `rk' [iw=brr1] if s2, `opt' contvars(xc 4)
matrix C3 = e(cuts1)
seldecomp y `rk' [pw=wt] if s1, `opt' contvars(xc 4)
assert mreldif(e(cuts1), C3) == 0
fake bootstrap : seldecomp y `rk', `opt' contvars(xc 4)
rcof `"noi fake bootstrap _loop _loop_bs : seldecomp y `rk', `opt' contvars(xc 4, refgroup)"' == 459
fake bootstrap : seldecomp y `rk', `opt' contvars(xc 4)
rcof `"noi fake bootstrap _loop _loop_bs : seldecomp y `rk', `opt' focalgroup(1) contvars(xc 4)"' == 459
* the message
fake bootstrap : seldecomp y `rk', `opt' contvars(xc 4)
tempfile mlog
log using `"`mlog'"', text replace name(cc)
capture noisily fake bootstrap _loop _loop_bs : seldecomp y `rk', `opt' contvars(xc 4, focalgroup)
log close cc
mata: st_local("found", strofreal(strpos(stritrim(invtokens(cat(st_local("mlog"))')), ///
    "contvars(): no cut points are recorded for this replicate; specify contvars(..., recut)") > 0))
assert `found' == 1
* errors, and e(recut): recut outside contvars() names its place
tempfile rlog
log using `"`rlog'"', text replace name(rc)
rcof "noi seldecomp y i.x1 xc, by(group) recut" == 198
rcof "noi seldecomp y `rk', `opt' contvars(xc 4) recut" == 198
rcof "noi seldecomp y `rk', `opt' contvars(xc 4, refgroup) recut" == 198
rcof "noi seldecomp y `rk', `opt' contvars(xc 4 recut)" == 198
log close rc
foreach m in "recut is a suboption of contvars(): type contvars(varname #, recut)" ///
    "recut is a suboption of contvars(): type contvars(xc 4, recut)" ///
    "recut is a suboption of contvars(): type contvars(xc 4, refgroup recut)" ///
    "contvars(): put a comma before recut, as in contvars(xc 4, recut)" {
    mata: st_local("found", strofreal(strpos(stritrim(invtokens(cat(st_local("rlog"))')), st_local("m")) > 0))
    assert `found' == 1
}
rcof "noi seldecomp y `rk', `opt' contvars(xc 4, recut) insupport bayesboot(reps(2) nodots)" == 198
seldecomp y `rk', `opt' contvars(xc 4, recut) insupport
assert "`e(recut)'" == "recut"
seldecomp y `rk', `opt' contvars(xc 4)
assert "`e(recut)'" == ""


* ==============================================================================
* 4. bayesboot(): the replicates against plain runs on the replicate weights
* ==============================================================================
di as txt _n "---- 4. bayesboot(): fixed cut points, and recut"
* default: every replicate equals a plain run on its weights with the full-sample bins
seldecomp y `rk' [pw=wt], `opt' contvars(xc 4) bayesboot(reps(5) rseed(3) nodots generate(fw))
matrix D = e(bsdraws)
matrix VF = e(V)
matrix BF = e(b)
assert rowsof(D) == 5 & colsof(D) == colsof(e(b))
forvalues r = 1/5 {
    qui seldecomp y (age: i.binw) (base: i.x1) [pw=fw_`r'], `opt'
    draws_are_b D `r' 1e-10
}
* recut: every replicate equals a plain run on its weights with c.xc, whose cut points are
* those weights' quantiles; the point estimates do not change, the variance does
seldecomp y `rk' [pw=wt], `opt' contvars(xc 4, recut) bayesboot(reps(5) rseed(3) nodots generate(fr))
assert "`e(recut)'" == "recut"
assert mreldif(e(b), BF) == 0
assert mreldif(e(V), VF) > 1e-4
matrix D = e(bsdraws)
forvalues r = 1/5 {
    assert fw_`r' == fr_`r'                  // the same draws
    qui seldecomp y `rk' [pw=fr_`r'], `opt' contvars(xc 4)
    draws_are_b D `r' 1e-10
}
drop fw_* fr_*
* interactions, two c. variables, focalgroup quantiles, reverse
local rk2 (a: c.xc#i.x1) (b: c.x2c)
seldecomp y `rk2' [pw=wt], `opt' reverse contvars(xc 3 x2c 4, focalgroup recut) ///
    bayesboot(reps(4) rseed(4) nodots generate(fr))
matrix D = e(bsdraws)
assert colsof(D) == colsof(e(b))
forvalues r = 1/4 {
    qui seldecomp y `rk2' [pw=fr_`r'], `opt' reverse contvars(xc 3 x2c 4, focalgroup)
    draws_are_b D `r' 1e-10
}
drop fr_*
* without reverse: every stage replicated
seldecomp y `rk2' [pw=wt], `opt' contvars(xc 3 x2c 4, recut) bayesboot(reps(4) rseed(5) nodots generate(fr))
matrix D = e(bsdraws)
assert colsof(D) == colsof(e(b))
forvalues r = 1/4 {
    qui seldecomp y `rk2' [pw=fr_`r'], `opt' contvars(xc 3 x2c 4)
    draws_are_b D `r' 1e-10
}
drop fr_*
* BCa with clusters and recut: the acceleration equals a delete-one-cluster jackknife made by
* hand with plain runs, which recompute the cut points without the cluster
seldecomp y `rk' [pw=wt], `opt' contvars(xc 4, recut) bayesboot(reps(5) rseed(6) bca cluster(psuid) nodots)
matrix A = e(accel)
local k = colsof(A)
qui levelsof psuid, local(pl)
local np : word count `pl'
matrix JK = J(`np', `k', .)
local names : colfullnames A
local r 0
foreach c of local pl {
    local ++r
    qui seldecomp y `rk' [pw=wt] if psuid != `c', `opt' contvars(xc 4)
    local j 0
    foreach nm of local names {
        local ++j
        matrix JK[`r', `j'] = el(e(b), 1, colnumb(e(b), "`nm'"))
    }
}
mata: st_matrix("AH", tc_accel(st_matrix("JK")))
forvalues j = 1/`k' {
    assert reldif(el(A, 1, `j'), el(AH, 1, `j')) < 1e-8
}
* without recut the same jackknife on the full-sample bins
seldecomp y `rk' [pw=wt], `opt' contvars(xc 4) bayesboot(reps(5) rseed(6) bca cluster(psuid) nodots)
matrix A = e(accel)
local r 0
foreach c of local pl {
    local ++r
    qui seldecomp y (age: i.binw) (base: i.x1) [pw=wt] if psuid != `c', `opt'
    local j 0
    foreach nm of local names {
        local ++j
        matrix JK[`r', `j'] = el(e(b), 1, colnumb(e(b), "`nm'"))
    }
}
mata: st_matrix("AH", tc_accel(st_matrix("JK")))
forvalues j = 1/`k' {
    assert reldif(el(A, 1, `j'), el(AH, 1, `j')) < 1e-8
}

* a BCa deletion that leaves the quantile sample empty (Reference Group is one cluster, refgroup
* quantiles) fails and is left out, with recut as without it, silently (no note on failed
* deletions)
preserve
keep if stratum <= 4
replace group = (psuid == 11)
tempfile dlog
log using `"`dlog'"', text replace name(dl)
seldecomp y `rk', `opt' contvars(xc 3, refgroup recut) bayesboot(reps(5) rseed(7) bca cluster(psuid) nodots)
log close dl
assert e(N_reps) == 5 & e(N_clust) == 8
mata: st_local("found", strofreal(strpos(stritrim(invtokens(cat(st_local("dlog"))')), ///
    "jackknife deletion") > 0))
assert `found' == 0
matrix A = e(accel)
assert !missing(A[1, 1])
restore

di as txt _n "test_cutpoints.do: all checks passed"

* test_chain.do -- chain: the decompositions of the adjacent
* pairs of the values of by() and their total. Each pair equals the run on its two values
* (1e-12: default, sequentially with fixothermargins and without it (the other blocks
* free), ratio(), contvars(), insupport, no weights, two values), the total equals the sums;
* the e() layout (e(stagetype), e(fixothermargins), e(diffsign) focal_group - ref_group: the total
* difference is the first value's mean minus the last's); bayesboot() (every stage):
* e(V) equals the covariance of the pairs and totals rebuilt from runs on the replicate
* weights of bayesboot(generate()) (1e-8), the explained shares' draws, BCa's acceleration
* against a hand-built cluster jackknife; svy brr against a hand-built BRR variance, its
* pre-flight check; a replicate of a prefix without one of the values fails, with one
* message naming it; the by: prefix; replay, replay at another level, estimates
* store/restore; makeplot (the whole chain and chainstep()); the notes; the display (78 columns, a
* long by() name); insupport and the telescoping identity; the errors (focalgroup(),
* individually with or without fixothermargins, reverse, detail, outcomes(), generate(),
* covariatedetail refused with chain, with their messages).
*
* Run from the tests directory:   . do test_chain.do

version 16
clear
set more off
set varabbrev off
discard             // drop any seldecomp already in memory; the repository copy is what gets tested

* ---- tolerances ---------------------------------------------------------------
local same 1e-12     // the same computation reached by two routes
local vtol 1e-8      // variances rebuilt from the replicates

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
* chk_pairs "<command before the comma>" "<options>" [by]: runs the chain, then each pair
* alone (if inlist(by, a, b)) with the same options, and asserts that every element of the
* pair's e(b), its e(bystage) and its row of e(by_pair) equal the chain's (1e-12), and
* that the total equals the sums over the pairs
capture program drop chk_pairs
program define chk_pairs
    args cmd opts byv
    if "`byv'" == "" local byv period
    seldecomp `cmd', by(`byv') chain `opts'
    tempname CB PR BP S1 P T
    matrix `CB' = e(b)
    matrix `PR' = e(pairs)
    matrix `BP' = e(by_pair)
    local J = e(N_pairs)
    assert rowsof(`PR') == `J'
    forvalues j = 1/`J' {
        tempname CS`j'
        if `j' == 1 matrix `CS`j'' = e(bystage)
        else matrix `CS`j'' = e(bystage`j')
    }
    local tdiff = e(difference)
    local texpl = e(explained)
    local tunex = e(unexplained)
    local blocks `e(blocks)'
    local staged = ("`e(stagetype)'" != "")
    local sd 0
    local se 0
    local su 0
    foreach bl of local blocks {
        local sb_`bl' 0
    }
    forvalues j = 1/`J' {
        local a = `PR'[`j', 1]
        local c = `PR'[`j', 2]
        // the command and its if (a chain run takes the if of the command too)
        if strpos(`" `cmd' "', " if ") {
            local cmdj `"`cmd' & inlist(`byv', `a', `c')"'
        }
        else {
            local cmdj `"`cmd' if inlist(`byv', `a', `c')"'
        }
        quietly seldecomp `cmdj', by(`byv') `opts'
        matrix `P' = e(b)
        matrix `S1' = e(bystage)
        local cn : colfullnames `P'
        foreach nm of local cn {
            gettoken eq rest : nm, parse(":")
            local e2 = cond("`eq'" == "overall", "p`j'", "p`j'_`eq'")
            local v1 = `P'[1, colnumb(`P', "`nm'")]
            local c2 = colnumb(`CB', "`e2'`rest'")
            assert `c2' < .
            local v2 = `CB'[1, `c2']
            if reldif(`v1', `v2') > 1e-12 {
                di as err "pair `j', `nm': alone " %21.15g `v1' ", in the chain " %21.15g `v2'
                exit 9
            }
        }
        assert mreldif(`S1', `CS`j'') <= 1e-12
        local k 0
        foreach s in mean_focal mean_raked_focal mean_ref difference explained unexplained pct_explained ///
            N_focal N_ref W_focal W_ref {
            assert reldif(e(`s'), `BP'[`j', `++k']) <= 1e-12
        }
        local sd = `sd' + _b[overall:difference]
        local se = `se' + _b[overall:explained]
        local su = `su' + _b[overall:unexplained]
        if `staged' {
            foreach bl of local blocks {
                local sb_`bl' = `sb_`bl'' + _b[marg_expl:`bl']
            }
        }
    }
    assert reldif(`sd', `CB'[1, colnumb(`CB', "total:difference")]) <= 1e-12
    assert reldif(`se', `CB'[1, colnumb(`CB', "total:explained")]) <= 1e-12
    assert reldif(`su', `CB'[1, colnumb(`CB', "total:unexplained")]) <= 1e-12
    assert reldif(`sd', `tdiff') <= 1e-12 & reldif(`se', `texpl') <= 1e-12 & reldif(`su', `tunex') <= 1e-12
    if `staged' {
        foreach bl of local blocks {
            assert reldif(`sb_`bl'', `CB'[1, colnumb(`CB', "total_marg_expl:`bl'")]) <= 1e-12
        }
    }
end

* log_count logfile "text": r(n), the lines of the log that contain text
capture program drop log_count
program define log_count, rclass
    args f txt
    tempname fh
    local n 0
    file open `fh' using `"`f'"', read text
    file read `fh' line
    while r(eof) == 0 {
        if strpos(`"`macval(line)'"', `"`txt'"') local ++n
        file read `fh' line
    }
    file close `fh'
    return scalar n = `n'
end

* log_width logfile: r(w), the widest line between the markers SDW_W_BEGIN and SDW_W_END
* that is not a command echo; r(tilde), those lines that contain ~ (an abbreviation)
capture program drop log_width
program define log_width, rclass
    args fn
    mata: _t_l = cat(st_local("fn"))
    mata: _t_b = select(1::rows(_t_l), _t_l :== "SDW_W_BEGIN")
    mata: _t_e = select(1::rows(_t_l), _t_l :== "SDW_W_END")
    mata: _t_l = _t_l[(_t_b[1] + 1)::(_t_e[1] - 1)]
    mata: _t_l = select(_t_l, substr(_t_l, 1, 2) :!= ". " :& substr(_t_l, 1, 2) :!= "> ")
    mata: st_local("w", strofreal(max(udstrlen(_t_l))))
    mata: st_local("t", strofreal(sum(strpos(_t_l, "~") :> 0)))
    mata: mata drop _t_l _t_b _t_e
    return scalar w = `w'
    return scalar tilde = `t'
end

* ---- the data: four periods (two years each), covariates whose distributions shift -----
* a BRR design of 8 strata of 2 PSUs, and clusters within periods, for the inference below
clear
set seed 20260927
set obs 1600
gen int id = _n
gen int period = 1995 + 2 * ceil(_n / 400)
label define per 1997 "97/98" 1999 "99/00" 2001 "01/02" 2003 "03/04"
label values period per
gen double wt = 50 + 100 * runiform()
gen byte age = 1 + (runiform() < .3 + .05 * (period - 1997) / 2) + (runiform() < .3)
gen byte sex = runiform() < .5
gen byte ins = 1 + (runiform() < .4 + .1 * (period - 1997) / 2) + (runiform() < .2)
gen byte reg = 1 + floor(4 * runiform())
gen double inc = 10 + 5 * age + 3 * (period - 1997) / 2 + 8 * runiform()
gen double y = 2 + .5 * age + .3 * sex - .4 * ins + .2 * reg + .15 * (period - 1997) / 2 + rnormal()
gen double d = 1 + runiform()
gen byte stratum = 1 + mod(_n, 8)
gen byte psu = 1 + mod(floor(_n / 8), 2)
gen int psuid = 10 * stratum + psu
gen int cl = 10 * period + 1 + mod(_n, 4)
mata: _t_H = (1); for (i = 1; i <= 3; i++) _t_H = (_t_H, _t_H \ _t_H, -_t_H); st_matrix("HAD", _t_H)
forvalues r = 1/8 {
    gen double brr`r' = 2 * wt * ((HAD[`r', stratum] == 1) == (psu == 1))
}
mata: mata drop _t_H

* ================================================================================
* 1. each pair equals the run on its two values; the total equals the sums
* ================================================================================
chk_pairs "y age sex ins [pw=wt]"
chk_pairs "y (demo: age sex) (cov: ins) reg [pw=wt]" "sequentially fixothermargins"
chk_pairs "y (demo: age sex) (cov: ins) reg [pw=wt]" "sequentially"
chk_pairs "y (demo: age sex) (cov: ins) reg [pw=wt]" "sequentially fixothermargins ratio(d)"
chk_pairs "y age sex ins" "sequentially fixothermargins"
chk_pairs "y age sex ins [aw=wt] if reg != 4" "sequentially fixothermargins"
* contvars(): each pair's cut points from its own sample, as the run on the pair
chk_pairs "y c.inc sex [pw=wt]" "contvars(inc 4) sequentially fixothermargins"
seldecomp y c.inc sex [pw=wt], by(period) chain contvars(inc 4)
forvalues j = 1/3 {
    matrix C`j' = e(cuts1_p`j')
}
matrix PR = e(pairs)
forvalues j = 1/3 {
    local a = PR[`j', 1]
    local c = PR[`j', 2]
    quietly seldecomp y c.inc sex [pw=wt] if inlist(period, `a', `c'), by(period) contvars(inc 4)
    assert mreldif(C`j', e(cuts1)) == 0
}
assert mreldif(C1, C2) > 0
assert "`e(cutvars)'" == "inc"
* insupport: a level of reg that occurs in 1997 but not in 1999 (Focal Group of pair 1 only)
gen byte reg5 = reg
replace reg5 = 5 if period == 1997 & mod(_n, 25) == 0
chk_pairs "y age sex reg5 [pw=wt]" "insupport"
seldecomp y age sex reg5 [pw=wt], by(period) chain insupport
matrix BP = e(by_pair)
assert BP[1, colnumb(BP, "N_focal_all")] - BP[1, colnumb(BP, "N_focal")] == 16
assert BP[2, colnumb(BP, "N_focal_all")] == BP[2, colnumb(BP, "N_focal")]
assert "`e(insupport)'" == "insupport"
* e(sample) and e(N): the union of the pairs' samples, without the 16 observations of 1997
* that pair 1 excludes and no pair uses
assert e(N) == 1600 - 16
quietly count if e(sample)
assert r(N) == e(N)
quietly count if e(sample) & period == 1997 & reg5 == 5
assert r(N) == 0
gen byte chs = e(sample)
matrix PRS = e(pairs)
gen byte uni = 0
forvalues j = 1/3 {
    local a = PRS[`j', 1]
    local c = PRS[`j', 2]
    quietly seldecomp y age sex reg5 [pw=wt] if inlist(period, `a', `c'), by(period) insupport
    quietly replace uni = 1 if e(sample)
}
quietly count if uni != chs
assert r(N) == 0
* two values: the chain's e(sample) and e(N) are the run's alone
quietly seldecomp y age sex reg5 [pw=wt] if period <= 1999, by(period) insupport
replace chs = e(sample)
local chn = e(N)
quietly seldecomp y age sex reg5 [pw=wt] if period <= 1999, by(period) chain insupport
assert e(N) == `chn'
quietly count if e(sample) != chs
assert r(N) == 0
drop chs uni
* bayesboot(): the draws are made on e(sample), so the generate() weights are missing for
* the 16 observations that no pair uses, as without chain
quietly seldecomp y age sex reg5 [pw=wt], by(period) chain insupport ///
    bayesboot(reps(3) rseed(5) nodots generate(cbw))
quietly count if !missing(cbw_1)
assert r(N) == e(N)
quietly count if !missing(cbw_1) & !e(sample)
assert r(N) == 0
drop cbw_*
seldecomp y age sex reg5 [pw=wt], by(period) chain insupport
* e(mean_focal) is pair 1's Focal Group mean after its screen, e(by_pair)[1, mean_focal_all] the first
* value's mean; only pair 1 screens, so the total difference is still e(mean_focal) - e(mean_ref)
quietly sum y [aw=wt] if period == 1997
assert reldif(BP[1, colnumb(BP, "mean_focal_all")], r(mean)) < 1e-12
assert reldif(e(mean_focal), BP[1, colnumb(BP, "mean_focal")]) == 0
assert reldif(e(difference), e(mean_focal) - e(mean_ref)) < 1e-12
* an inner value screened (a level of reg6 in 1997 and 1999, not in 2001: pair 2 screens
* 1999, which pair 1 took whole as its Reference Group): the total difference, the sum of the pairs,
* departs from e(mean_focal) - e(mean_ref) by 1999's screened mean minus its full mean
gen byte reg6 = reg
replace reg6 = 5 if period <= 1999 & mod(_n, 25) == 0
chk_pairs "y age sex reg6 [pw=wt]" "insupport"
seldecomp y age sex reg6 [pw=wt], by(period) chain insupport
matrix BP = e(by_pair)
assert BP[1, colnumb(BP, "N_focal_all")] == BP[1, colnumb(BP, "N_focal")]
assert BP[2, colnumb(BP, "N_focal_all")] - BP[2, colnumb(BP, "N_focal")] == 16
quietly sum y [aw=wt] if period == 1997
assert reldif(e(mean_focal), r(mean)) < 1e-12
quietly sum y [aw=wt] if period == 2003
assert reldif(e(mean_ref), r(mean)) < 1e-12
assert reldif(BP[2, colnumb(BP, "mean_focal_all")], BP[1, colnumb(BP, "mean_ref")]) < 1e-12
local gap = BP[2, colnumb(BP, "mean_focal")] - BP[2, colnumb(BP, "mean_focal_all")]
assert abs(`gap') > 1e-6
assert abs(e(difference) - (e(mean_focal) - e(mean_ref)) - `gap') < 1e-12
drop reg6
* two values: one pair, whose total is itself, and the run without chain
chk_pairs "y (demo: age sex) (cov: ins) [pw=wt] if period <= 1999" "sequentially fixothermargins"
seldecomp y (demo: age sex) (cov: ins) [pw=wt] if period <= 1999, by(period) chain sequentially fixothermargins
assert e(N_pairs) == 1
matrix CB = e(b)
quietly seldecomp y (demo: age sex) (cov: ins) [pw=wt] if period <= 1999, by(period) sequentially fixothermargins
assert reldif(CB[1, colnumb(CB, "total:explained")], _b[overall:explained]) == 0

* ================================================================================
* 2. the e() layout
* ================================================================================
seldecomp y (demo: age sex) (cov: ins) reg [pw=wt], by(period) chain sequentially fixothermargins
assert "`e(cmd)'" == "seldecomp" & "`e(chain)'" == "chain"
assert e(N_pairs) == 3 & e(N) == 1600 & e(N_stages) == 3 & e(k_eq) == 3 * 5 + 2
assert "`e(by_levels)'" == "1997 1999 2001 2003"
assert `"`e(by_labels)'"' == "1997 (97/98) | 1999 (99/00) | 2001 (01/02) | 2003 (03/04)"
assert "`e(focal_group)'" == "1997" & "`e(ref_group)'" == "2003"
assert `"`e(focal_group_label)'"' == "1997 (97/98)" & `"`e(ref_group_label)'"' == "2003 (03/04)"
assert "`e(diffsign)'" == "focal_group - ref_group" & "`e(stagetype)'" == "sequentially"
assert "`e(fixothermargins)'" == "fixothermargins"
matrix PR = e(pairs)
assert PR[1, 1] == 1997 & PR[1, 2] == 1999 & PR[3, 1] == 2001 & PR[3, 2] == 2003
assert "`: rownames PR'" == "p1 p2 p3" & "`: colnames PR'" == "focal_group ref_group"
* e(b): per pair the ordinary layout prefixed, then the total
local cn : colfullnames e(b)
assert "`: word 1 of `cn''" == "p1:focal_group" & "`: word 7 of `cn''" == "p1_marg_expl:demo"
assert "`: word 10 of `cn''" == "p1_cum_expl:demo" & "`: word 19 of `cn''" == "p2:focal_group"
assert "`: word 55 of `cn''" == "total:difference" & "`: word 58 of `cn''" == "total_marg_expl:demo"
assert colsof(e(b)) == 60
* the telescoping identity: the total difference is the first mean minus the last
assert reldif(e(difference), e(mean_focal) - e(mean_ref)) < 1e-12
assert reldif(e(mean_focal), _b[p1:focal_group]) == 0 & reldif(e(mean_ref), _b[p3:ref_group]) == 0
assert reldif(e(pct_explained), 100 * e(explained) / e(difference)) < 1e-12
matrix S2 = e(bystage2)
assert rowsof(S2) == 3 & colsof(S2) == 30
assert reldif(S2[3, 2], _b[p2:explained]) < 1e-12
assert e(converged) == 1 & e(targets_met) == 1 & e(N_zerowt) == 0
* not posted under chain: the pair quantities
foreach s in N_focal N_ref W_focal W_ref mean_raked_focal {
    assert e(`s') >= .
}
assert "`e(levels1)'" == ""

* ================================================================================
* 3. bayesboot(): one replicate loop over the pairs
* ================================================================================
* (a) e(V) against the covariance rebuilt from runs of each pair on the replicate weights
local R 12
seldecomp y (demo: age sex) (cov: ins) reg [pw=wt], by(period) chain sequentially fixothermargins ///
    bayesboot(reps(`R') rseed(7) generate(bw) nodots)
matrix V = e(V)
matrix B = e(b)
matrix D = e(bsdraws)
matrix DP = e(bsdraws_pct)
matrix BPB = e(by_pair)
matrix PCI = e(pct_explained_ci)
local tpse = e(pct_explained_se)
local tpct = e(pct_explained)
local k = colsof(B)
assert "`e(bsnames)'" == "`: colfullnames e(b)'"
assert e(N_reps) == `R' & e(N_misreps) == 0 & "`e(vce)'" == "bayesboot"
matrix PR = e(pairs)
matrix RB = J(`R', `k', .)
matrix RP = J(`R', 4, .)
forvalues r = 1/`R' {
    local c 0
    local td 0
    local te 0
    local tu 0
    foreach bl in demo cov reg {
        local tb_`bl' 0
    }
    forvalues j = 1/3 {
        local a = PR[`j', 1]
        local cc = PR[`j', 2]
        quietly seldecomp y (demo: age sex) (cov: ins) reg [pw=bw_`r'] if inlist(period, `a', `cc'), ///
            by(period) sequentially fixothermargins
        matrix P = e(b)
        forvalues i = 1/`=colsof(P)' {
            matrix RB[`r', `++c'] = P[1, `i']
        }
        matrix RP[`r', `j'] = e(pct_explained)
        local td = `td' + _b[overall:difference]
        local te = `te' + _b[overall:explained]
        local tu = `tu' + _b[overall:unexplained]
        foreach bl in demo cov reg {
            local tb_`bl' = `tb_`bl'' + _b[marg_expl:`bl']
        }
    }
    matrix RB[`r', `++c'] = `td'
    matrix RB[`r', `++c'] = `te'
    matrix RB[`r', `++c'] = `tu'
    foreach bl in demo cov reg {
        matrix RB[`r', `++c'] = `tb_`bl''
    }
    matrix RP[`r', 4] = 100 * `te' / `td'
}
assert mreldif(RB, D) < 1e-10
assert mreldif(RP, DP) < 1e-10
mata: st_matrix("VH", variance(st_matrix("RB")))
assert mreldif(VH, V) < `vtol'
* the totals' variances, and their covariances with the pairs, are the rebuilt ones
local ct = colnumb(V, "total:explained")
assert V[`ct', `ct'] > 0 & reldif(V[`ct', `ct'], VH[`ct', `ct']) < `vtol'
assert reldif(V[`ct', colnumb(V, "p2:explained")], VH[`ct', colnumb(V, "p2:explained")]) < `vtol'
* the shares' std. errors: of each pair (e(by_pair)) and of the total
mata: st_matrix("SDP", sqrt(diagonal(variance(st_matrix("RP"))))')
forvalues j = 1/3 {
    assert reldif(BPB[`j', colnumb(BPB, "pct_explained_se")], SDP[1, `j']) < `vtol'
}
assert reldif(`tpse', SDP[1, 4]) < `vtol'
assert colsof(PCI) == 2 & PCI[1, 1] < `tpct' & PCI[1, 2] > `tpct'
* the replicate weights: drawn once on the whole chain sample (every period)
forvalues r = 1/`R' {
    quietly count if bw_`r' >= .
    assert r(N) == 0
}
drop bw_*
* the draws do not depend on the pairs: the same seed on the same sample gives the same
* weights with another specification
seldecomp y age sex [pw=wt], by(period) chain bayesboot(reps(3) rseed(7) generate(bw2) nodots)
seldecomp y ins reg [pw=wt], by(period) chain bayesboot(reps(3) rseed(7) generate(bw3) nodots)
forvalues r = 1/3 {
    assert bw2_`r' == bw3_`r'
}
* generate(, replace) overwrites them (another seed), and never an input of the command
rcof "noi seldecomp y ins reg [pw=wt], by(period) chain bayesboot(reps(3) rseed(7) generate(bw3) nodots)" == 110
seldecomp y ins reg [pw=wt], by(period) chain bayesboot(reps(3) rseed(8) generate(bw3, replace) nodots)
count if reldif(bw2_1, bw3_1) > 1e-6 & bw2_1 < .
assert r(N) > 0
gen double wb_1 = wt
rcof "noi seldecomp y ins reg [pw=wb_1], by(period) chain bayesboot(reps(3) generate(wb, replace) nodots)" == 198
drop bw2_* bw3_* wb_1

* (b) BCa: the acceleration against a hand-built jackknife over clusters within periods,
* whose deletions touch only the pairs that contain the period (the others keep their point
* estimates); the replicates are the same as without bca
seldecomp y age sex ins [pw=wt], by(period) chain bayesboot(reps(`R') rseed(3) cluster(cl) bca nodots)
matrix A = e(accel)
matrix B = e(b)
matrix CIB = e(ci_bca)
assert rowsof(CIB) == 2 & "`: colfullnames CIB'" == "`: colfullnames A'"
local kb = colsof(A)
assert `kb' == 3 * 6 + 3 & e(N_clust) == 16
levelsof cl, local(cls)
matrix JK = J(16, `kb', .)
local u 0
foreach cv of local cls {
    local ++u
    local c 0
    local td 0
    local te 0
    local tu 0
    forvalues j = 1/3 {
        local a = PR[`j', 1]
        local cc = PR[`j', 2]
        quietly seldecomp y age sex ins [pw=wt] if inlist(period, `a', `cc') & cl != `cv', by(period)
        forvalues i = 1/6 {
            matrix JK[`u', `++c'] = _b[overall:`: word `i' of focal_group raked_focal ref_group difference explained unexplained']
        }
        local td = `td' + _b[overall:difference]
        local te = `te' + _b[overall:explained]
        local tu = `tu' + _b[overall:unexplained]
    }
    matrix JK[`u', `++c'] = `td'
    matrix JK[`u', `++c'] = `te'
    matrix JK[`u', `++c'] = `tu'
}
mata: _t_J = st_matrix("JK"); _t_d = mean(_t_J) :- _t_J
mata: st_matrix("AH", colsum(_t_d:^3) :/ (6 :* colsum(_t_d:^2):^1.5))
mata: mata drop _t_J _t_d
assert mreldif(AH, A) < `vtol'

* (c) saving(): a bootstrap dataset that bstat reads, with the pairs' and the total's columns
tempfile bsf
seldecomp y age sex ins [pw=wt], by(period) chain bayesboot(reps(`R') rseed(3) saving("`bsf'") nodots)
matrix SE0 = e(se)
preserve
use "`bsf'", clear
confirm variable p1_focal_group p3_unexplained t_difference t_explained t_unexplained t_pct_explained
assert _N == `R'
bstat using "`bsf'"
matrix SE1 = e(se)
restore
matrix SE1 = SE1[1, 1..21]
assert mreldif(SE1, SE0) < 1e-10

* (d) failed replicates: pair 2 has the design of test_bootstrap.do 9b (Focal Group, period 2,
* lacks the cell a = 1, b = 1, so a replicate whose draws reverse the order of the Reference Group's
* shares of (1, 1) and (0, 0) has margins that cannot be met); pair 1 is benign. A replicate
* that fails in pair 2 is dropped as a whole; one note names the pair
preserve
clear
set obs 600
gen byte per = 1 + (_n > 200) + (_n > 400)
gen int k = mod(_n - 1, 200) + 1
gen byte a = cond(per == 1, mod(k, 2), cond(per == 2, k > 130, cond(k <= 30, 0, cond(k <= 54, 1, !mod(k, 2)))))
gen byte b = cond(per == 1, mod(floor(k / 2), 2), cond(per == 2, inrange(k, 61, 130), cond(k <= 30, 0, cond(k <= 54, 1, mod(k, 2)))))
gen double yy = a + 2 * b + mod(_n * 7, 11) / 10 + per / 4
tempfile flg
log using `"`flg'"', text replace name(flog)
seldecomp yy a b, by(per) chain bayesboot(reps(20) rseed(7) nodots)
log close flog
local nm = e(N_misreps)
assert inrange(`nm', 2, 18) & e(N_reps) + `nm' == 20 & "`e(missing)'" == "missing"
matrix RR = e(reps)
assert RR[1, 1] == e(N_reps) & RR[1, colsof(RR)] == e(N_reps)
matrix DD = e(bsdraws)
mata: st_local("nmr", strofreal(sum(rowmissing(st_matrix("DD")) :== cols(st_matrix("DD")))))
assert `nmr' == `nm'
log_count `"`flg'"' "note: `nm' of 20 replications failed and were dropped: in them, raking"
assert r(n) == 1
log_count `"`flg'"' "(pair p2)"
assert r(n) == 1
restore

* ================================================================================
* 4. svy brr: every call runs the whole chain; e(V) against a hand-built BRR variance
* ================================================================================
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr8)
svy brr, nodots: seldecomp y (demo: age sex) (cov: ins) reg, by(period) chain sequentially fixothermargins
matrix VS = e(V)
matrix BS = e(b)
local k = colsof(BS)
quietly seldecomp y (demo: age sex) (cov: ins) reg [pw=wt], by(period) chain sequentially fixothermargins
assert mreldif(e(b), BS) < `same'
matrix RR = J(8, `k', .)
forvalues r = 1/8 {
    quietly seldecomp y (demo: age sex) (cov: ins) reg [pw=brr`r'], by(period) chain sequentially fixothermargins
    matrix BR = e(b)
    forvalues i = 1/`k' {
        matrix RR[`r', `i'] = BR[1, `i']
    }
}
* svy brr's default: deviations from the replicates' mean, divided by R
mata: _t_R = st_matrix("RR"); _t_d = _t_R :- mean(_t_R); st_matrix("VB", _t_d' * _t_d / rows(_t_R))
mata: mata drop _t_R _t_d
assert mreldif(VB, VS) < `vtol'
* the replay after svy: the chain's own table
seldecomp
svyset, clear

* the pre-flight check over the pairs: a level of x5 that period 1999 has only in PSU 1 of
* stratum 2 fails pair 2 (1999 is its Focal Group) in the half-samples that drop that PSU; one
* warning, and e(N_prefail) counts each replicate once
gen byte x5 = 1 + (reg > 2)
replace x5 = 3 if mod(_n, 17) == 0 & period != 1999
replace x5 = 3 if period == 1999 & stratum == 2 & psu == 1 & mod(_n, 3) == 0
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr8)
tempfile lg
log using `"`lg'"', text replace name(pflog)
capture noisily svy brr, nodots: seldecomp y age x5, by(period) chain
log close pflog
local rc = _rc
assert `rc' == 0
assert e(N_prefail) == 4
assert e(N_misreps) == 4
log_count `"`lg'"' "warning:"
assert r(n) == 1
log_count `"`lg'"' "p2: x5 = 3"
assert r(n) == 1
svyset, clear

* the other prefixes: bootstrap: and jackknife: (by cluster), svy jackknife; contvars() under
* a prefix: each pair's replicates bin with that pair's full-sample cut points (a key of its
* own: with one key for all pairs the replicates would fail as ambiguous, 459)
quietly seldecomp y (demo: age sex) (cov: ins), by(period) chain sequentially fixothermargins
matrix B4 = e(b)
bootstrap, reps(4) seed(11) nodots: seldecomp y (demo: age sex) (cov: ins), by(period) chain sequentially fixothermargins
assert mreldif(e(b), B4) < `same' & e(N_reps) == 4 & "`e(chain)'" == "chain"
jackknife, cluster(psuid) nodots: seldecomp y (demo: age sex) (cov: ins), by(period) chain sequentially fixothermargins
assert mreldif(e(b), B4) < `same' & e(N_reps) == 16
quietly seldecomp y (demo: age sex) (cov: ins) [pw=wt], by(period) chain sequentially fixothermargins
matrix B4 = e(b)
svyset psuid [pw=wt], strata(stratum)
svy jackknife, nodots: seldecomp y (demo: age sex) (cov: ins), by(period) chain sequentially fixothermargins
assert mreldif(e(b), B4) < `same' & e(N_misreps) == 0
svyset psuid [pw=wt], strata(stratum) brrweight(brr1-brr8)
svy brr, nodots: seldecomp y c.inc sex, by(period) chain contvars(inc 4)
assert e(N_misreps) == 0
matrix VC = e(V)
* the replicates keep the pair's full-sample cut points: the replicate of half-sample 1
* equals the run on its weights with those cut points, which contvars(..., recut) would not use
svy brr, nodots: seldecomp y c.inc sex, by(period) chain contvars(inc 4, recut)
assert mreldif(e(V), VC) > 1e-6
svyset, clear

* ================================================================================
* 4b. a replicate of a prefix in which a value of by() has no observations or no weight
* fails (it would run a shorter chain whose pairs join values that are not adjacent), and
* one message names the value
* ================================================================================
* log_find logfile "text": r(n), the occurrences of text in the log with its lines joined
* (trimmed, one blank between them), so that a wrapped message is found whole
capture program drop log_find
program define log_find, rclass
    args fn txt
    mata: _t_l = strtrim(cat(st_local("fn")))
    mata: _t_s = invtokens(_t_l', " ")
    mata: _t_n = 0; _t_p = strpos(_t_s, st_local("txt")); while (_t_p) { _t_n++; _t_s = substr(_t_s, _t_p + 1, .); _t_p = strpos(_t_s, st_local("txt")); }
    mata: st_local("n", strofreal(_t_n))
    mata: mata drop _t_l _t_s _t_n _t_p
    return scalar n = `n'
end
* the design: 1999 lies entirely in PSU 1 of stratum 2 (psv 21); for the svy jackknife, 2003
* entirely in PSU 2 of stratum 3 (psw 32)
gen byte stv = cond(period == 1999, 2, stratum)
gen byte psv = cond(period == 1999, 1, psu)
gen int psvid = 10 * stv + psv
forvalues r = 1/8 {
    gen double brv`r' = 2 * wt * ((HAD[`r', stv] == 1) == (psv == 1))
}
gen byte stw = cond(period == 2003, 3, stratum)
gen byte psw = cond(period == 2003, 2, psu)
gen int pswid = 10 * stw + psw
* jackknife: by PSU: the deletion of PSU 21 leaves no 1999 and fails, also with a prefix
* expression (the shorter chain's total is not kept); the full-sample results are the
* chain's
quietly seldecomp y age sex [pw=wt] if psvid != 21, by(period) chain
assert e(N_pairs) == 2 & "`e(by_levels)'" == "1997 2001 2003"
local tedel = _b[total:explained]
quietly seldecomp y age sex [pw=wt], by(period) chain
local tefull = _b[total:explained]
tempfile jkf
jackknife te = _b[total:explained], cluster(psvid) nodots saving(`"`jkf'"', replace double): ///
    seldecomp y age sex [pw=wt], by(period) chain
assert e(N_misreps) == 1 & e(N_reps) == 15
assert reldif(_b[te], `tefull') < `same'
preserve
use `"`jkf'"', clear
quietly count if reldif(te, `tedel') < 1e-12
assert r(N) == 0
quietly count if te < .
assert r(N) == 15
restore
* the whole e(b): one warning after the loop names the value
tempfile lg
log using `"`lg'"', text replace name(vlog)
jackknife, cluster(psvid) nodots: seldecomp y age sex [pw=wt], by(period) chain
log close vlog
assert e(N_misreps) == 1
log_count `"`lg'"' "warning:"
assert r(n) == 1
log_find `"`lg'"' "1 of 16 jackknife replicates failed and was dropped: in it, period = 1999 had no weight."
assert r(n) == 1
* bootstrap: by PSU: the resamples without PSU 21 fail; the remedy is the strata
log using `"`lg'"', text replace name(vlog)
bootstrap, reps(20) seed(3) cluster(psvid) nodots: seldecomp y age sex, by(period) chain
log close vlog
assert e(N_misreps) > 0 & e(N_misreps) < 20
log_count `"`lg'"' "warning:"
assert r(n) == 1
log_find `"`lg'"' "period = 1999 had no weight. Use bootstrap, strata(period)."
assert r(n) == 1
* svy brr, a prefix expression: the pre-flight check's one warning names the value (not the
* levels of pair 2), the 4 half-samples without 1999 fail, and e(V) is the BRR variance
* over the other 4
svyset psvid [pw=wt], strata(stv) brrweight(brv1-brv8)
log using `"`lg'"', text replace name(vlog)
svy brr te = _b[total:explained], nodots: seldecomp y age sex, by(period) chain
log close vlog
assert e(N_misreps) == 4
matrix VS = e(V)
log_count `"`lg'"' "warning:"
assert r(n) == 1
log_find `"`lg'"' "4 of 8 BRR replicates will fail: in them, period = 1999 has no weight. Use Fay's BRR (fay())."
assert r(n) == 1
log_count `"`lg'"' "collapse"
assert r(n) == 0
matrix RT = J(4, 1, .)
local i 0
forvalues r = 1/8 {
    quietly count if period == 1999 & brv`r' > 0
    if r(N) == 0 continue
    quietly seldecomp y age sex [pw=brv`r'], by(period) chain
    assert e(N_pairs) == 3
    matrix RT[`++i', 1] = _b[total:explained]
}
assert `i' == 4
mata: _t_R = st_matrix("RT"); _t_d = _t_R :- mean(_t_R); st_numscalar("VH", _t_d' * _t_d / rows(_t_R))
mata: mata drop _t_R _t_d
assert reldif(VS[1, 1], scalar(VH)) < `vtol'
* svy jackknife: the last value in one PSU; its deletion fails, and the pre-flight check
* counts it in e(N_prefail) and says so
svyset pswid [pw=wt], strata(stw)
log using `"`lg'"', text replace name(vlog)
svy jackknife, nodots: seldecomp y age sex, by(period) chain
log close vlog
assert e(N_misreps) == 1 & e(N_prefail) == 1
log_count `"`lg'"' "warning:"
assert r(n) == 1
log_find `"`lg'"' "1 jackknife replicate will fail: in it, period = 2003 has no weight."
assert r(n) == 1
svyset, clear
* bayesboot() keeps every value (positive draws); a replicate with all values passes the check
seldecomp y age sex [pw=wt], by(period) chain bayesboot(reps(4) rseed(2) cluster(psvid) nodots)
assert e(N_misreps) == 0 & e(N_reps) == 4
* a program that runs the chain on two samples with the same specification under one prefix:
* each full-sample call records its values, and a replicate passes when its values are
* one of the recorded sets (the second chain's 1997-2001 do not excuse the first chain's
* deletion of 1999)
capture program drop two_chains
program define two_chains, eclass
    syntax [if]
    marksample touse
    quietly count if `touse'
    local N = r(N)
    quietly seldecomp y age sex [pw=wt] if `touse', by(period) chain
    local a = _b[total:explained]
    quietly seldecomp y age sex [pw=wt] if `touse' & period <= 2001, by(period) chain
    local c = _b[total:explained]
    tempname B
    matrix `B' = (`a', `c')
    matrix colnames `B' = all early
    ereturn post `B', obs(`N') esample(`touse')
end
jackknife, cluster(psvid) nodots: two_chains
assert e(N_misreps) == 1 & e(N_reps) == 15
* (after results of jackknife without e(cmd), a command line that starts with a comma replays)
ereturn clear
jackknife, cluster(psuid) nodots: two_chains
assert e(N_misreps) == 0 & e(N_reps) == 16
drop stv psv psvid brv* stw psw pswid

* ================================================================================
* 5. the by: prefix: each by-group runs its own chain
* ================================================================================
gen byte half = mod(_n, 2)
sort half, stable
by half: seldecomp y (demo: age sex) (cov: ins) [pw=wt], by(period) chain sequentially fixothermargins
matrix BB = e(b)
seldecomp y (demo: age sex) (cov: ins) [pw=wt] if half == 1, by(period) chain sequentially fixothermargins
assert mreldif(BB, e(b)) < `same'
* every by-group, not only the last: groups of unequal size, through statsby
gen byte third = mod(id, 3) == 0
forvalues g = 0/1 {
    quietly seldecomp y (demo: age sex) (cov: ins) [aw=wt] if third == `g', by(period) chain sequentially fixothermargins
    local n`g' = e(N)
    local d`g' = e(difference)
    local e`g' = _b[total_marg_expl:cov]
}
preserve
statsby n = e(N) d = e(difference) e = _b[total_marg_expl:cov], by(third) clear double nodots: ///
    seldecomp y (demo: age sex) (cov: ins) [aw=wt], by(period) chain sequentially fixothermargins
sort third
forvalues g = 0/1 {
    assert n[`g' + 1] == `n`g''
    assert reldif(d[`g' + 1], `d`g'') < `same' & reldif(e[`g' + 1], `e`g'') < `same'
}
assert `n0' != `n1'
restore
drop third

* ================================================================================
* 6. replay, replay at another level, estimates store and restore
* ================================================================================
seldecomp y (demo: age sex) (cov: ins) reg [pw=wt], by(period) chain sequentially fixothermargins ///
    bayesboot(reps(`R') rseed(7) nodots)
matrix B6 = e(b)
matrix V6 = e(V)
seldecomp
seldecomp, level(90)
matrix CN = r(ci_normal)
local c = colnumb(B6, "total:explained")
assert reldif(CN[1, `c'], B6[1, `c'] - invnormal(.95) * sqrt(V6[`c', `c'])) < 1e-12
matrix PCI = r(pct_explained_ci)
assert rowsof(PCI) == 4 & "`: rownames PCI'" == "p1 p2 p3 total"
assert r(level) == 90
estimates store ch1
seldecomp y age sex, by(period) chain
estimates restore ch1
assert mreldif(e(b), B6) == 0
seldecomp, citable diagnostics

* ================================================================================
* 7. makeplot: the whole chain without chainstep(), a step with chainstep()
* ================================================================================
* chain_bars G: the bar values of the marginal plot, from the top (the rbar plots 1-3)
capture program drop chain_bars
program define chain_bars
    args G
    tempname M
    preserve
    local rows
    forvalues k = 1/3 {
        local P .`G'.plotregion1.plot`k'
        local id ``P'.serset.id'
        qui serset set `id'
        qui serset use, clear
        mata: for (i = 1; i <= st_nvar(); i++) st_varrename(i, "v" + strofreal(i))
        local iy ``P'.y2var'
        local ix ``P'.xvar'
        mata: _t_m = st_data(., ("v`ix'", "v`iy'"))
        if `k' == 1 mata: _t_all = _t_m
        else mata: _t_all = _t_all \ _t_m
    }
    restore
    mata: _t_all = sort(_t_all, 1); st_matrix("CBARS", _t_all[., 2]')
    mata: mata drop _t_m _t_all
end
* the default plot of the whole chain is the waterfall
seldecomp, makeplot
assert `"`.seldecomp.title.text[1]'"' == "Decomposition of the difference by raking stage"
seldecomp, makeplot(marginal)
graph describe seldecomp
chain_bars seldecomp
matrix EX = (_b[total:difference], _b[total:explained], _b[total_marg_expl:demo], ///
    _b[total_marg_expl:cov], _b[total_marg_expl:reg], _b[total:unexplained])
assert mreldif(CBARS, EX) < 1e-10
assert strpos(`"`.seldecomp.subtitle.text[1]'"', "(sum of 3 pairs)")
* a step's default plot is the waterfall too, as a run's without chain; stages draws a step
seldecomp, makeplot(, chainstep(2))
assert `"`.seldecomp.title.text[1]'"' == "Decomposition of the difference by raking stage"
seldecomp, makeplot(stages, chainstep(2))
assert "`.seldecomp.plotregion1.plot3.classname'" == "function_g"
seldecomp, makeplot(marginal, chainstep(2))
chain_bars seldecomp
matrix EX = (_b[p2:difference], _b[p2:explained], _b[p2_marg_expl:demo], ///
    _b[p2_marg_expl:cov], _b[p2_marg_expl:reg], _b[p2:unexplained])
assert mreldif(CBARS, EX) < 1e-10
seldecomp, makeplot(waterfall)
seldecomp, makeplot(waterfall, chainstep(1))
seldecomp, makeplot(waterfalllevels, chainstep(3))
seldecomp, makeplot(stages, chainstep(1))
* chainstep(#): the #th pair of adjacent values (step 1: 1997 to 1999, by their labels)
assert `"`.seldecomp.subtitle.text[1]'"' == "97/98 - 99/00"
graph drop _all
* at estimation
seldecomp y (demo: age sex) (cov: ins) [pw=wt], by(period) chain sequentially fixothermargins makeplot(marginal, chainstep(1))
graph describe seldecomp
graph drop _all

* ================================================================================
* 8. the notes: the chain's exclusions (by value of by()), then each pair's, merged
* ================================================================================
gen double wt0 = wt
replace wt0 = 0 if inlist(id, 1, 2, 3, 801, 802)
gen byte rg = 1 + mod(_n, 3)
bysort period (rg): replace rg = 4 if _n <= 15
sort period, stable
seldecomp y age rg [pw=wt0], by(period) chain mincell(16)
assert e(N_zerowt) == 5
assert `"`e(note1)'"' == "5 observations with zero weight excluded (period = 1997: 3; 2001: 2)"
* the thin target level rg = 4 (15 obs, 12 with period 1997 left after the zero weights ...): the
* notes of the pairs start with the pairs, and the same note in several pairs is given once
local nh 0
forvalues m = 2/`e(N_notes)' {
    local nt `"`e(note`m')'"'
    if substr(`"`nt'"', 1, 4) != "    " {
        assert regexm(`"`nt'"', "^p[0-9]+(, p[0-9]+)*: ")
        local ++nh
    }
}
assert `nh' >= 1
tempfile tlg
log using `"`tlg'"', text replace name(tlog)
seldecomp y age rg [pw=wt], by(period) chain mincell(16) diagnostics
log close tlog
assert `"`e(note1)'"' == "p1, p2, p3: thin target level:"
assert `"`e(note2)'"' == "    rg = 4 (15 obs)"
assert e(N_notes) == 2
* Diagnostics summarize each pair's final thin-cell count; the detailed cell list
* remains in e(note#), outside the header and diagnostic display.
log_count `"`tlg'"' "note: "
assert r(n) == 0
log_count `"`tlg'"' "Thin target levels (<16 focal obs): 1."
assert r(n) == 3
log_count `"`tlg'"' "rg = 4 (15 obs)"
assert r(n) == 0

drop wt0 rg

* ================================================================================
* 9. the display: 78 columns, the pairs table and the panels
* ================================================================================
tempfile dlg
log using `"`dlg'"', text replace name(dlog)
di "SDW_W_BEGIN"
seldecomp y (demo: age sex) (cov: ins) reg [pw=wt], by(period) chain sequentially fixothermargins ///
    ratio(d) bayesboot(reps(`R') rseed(7) nodots citable)
seldecomp, diagnostics
di "SDW_W_END"
log close dlog
log_width `"`dlg'"'
assert r(w) <= 78 & r(tilde) == 0
foreach t in "Chain of adjacent pairs of period" "Pair  Focal Group (raked)" ///
    "total: the sum of the pairs" "Explained share of the difference (%)" ///
    "Decomposition by raking stage (sequentially, fixothermargins): p3," ///
    "Pair p2" "Ratio of means:" {
    log_count `"`dlg'"' "`t'"
    assert r(n) >= 1
}
* the pairs' rows of the coefficient table name their values; one per pair and the total
* (the coefficient table, the citable table, and the table of the replay)
log_count `"`dlg'"' "p2               | period 1999 (99/00) to 2001 (01/02)"
assert r(n) == 3
log_count `"`dlg'"' "total            | period 1997 (97/98) to 2003 (03/04), the sum of the pairs"
assert r(n) == 3
* a by() name of 21 characters: the header line follows the Pairs
* row, the pairs table shows the values (the name would hide them), the stage tables'
* titles leave the name out where it would pass 78 columns (and wrap where they still do,
* as with fixothermargins)
clonevar interview_year_pooled = period
tempfile dlg2
log using `"`dlg2'"', text replace name(dlog2)
di "SDW_W_BEGIN"
seldecomp y (demo: age sex) (cov: ins) [pw=wt], by(interview_year_pooled) chain sequentially fixothermargins ///
    diagnostics
di "SDW_W_END"
log close dlog2
log_width `"`dlg2'"'
assert r(w) <= 78 & r(tilde) == 0
log_count `"`dlg2'"' "pooledPairs"
assert r(n) == 0
foreach t in "Chain of adjacent pairs of interview_year_pooled" ///
    "  p1  1997 (97/98)               1999 (99/00)  " ///
    "  p3  2001 (01/02)               2003 (03/04)  " ///
    "Decomposition by raking stage (sequentially, fixothermargins): p2," ///
    "Pair p2" {
    log_count `"`dlg2'"' "`t'"
    assert r(n) == 1
}
log_count `"`dlg2'"' "                                                Pairs             =          3"
assert r(n) == 1
drop interview_year_pooled

* ================================================================================
* 10. the errors
* ================================================================================
rcof "noi seldecomp y age sex, by(period) chain focalgroup(1999)" == 198
rcof "noi seldecomp y (a: age) (b: sex), by(period) chain individually" == 198
rcof "noi seldecomp y (a: age) (b: sex), by(period) chain individually fixothermargins" == 198
rcof "noi seldecomp y (a: age) (b: sex), by(period) chain sequentially fixothermargins reverse" == 198
rcof "noi seldecomp y (a: age) (b: sex), by(period) chain sequentially reverse" == 198
rcof "noi seldecomp y age sex, by(period) chain outcomes(d)" == 198
rcof "noi seldecomp y age sex, by(period) chain generate(g)" == 198
rcof "noi seldecomp y age sex, by(period) chain covariatedetail" == 198
* detail: refused, with
* sequentially (with or without fixothermargins) or neither (chain's message, not detail's
* own), the results in memory kept; individually and focalgroup() get chain's messages too
seldecomp y age sex, by(period) chain
matrix B0 = e(b)
rcof "noi seldecomp y (a: age sex) (b: ins), by(period) chain sequentially fixothermargins detail" == 198
rcof "noi seldecomp y (a: age sex) (b: ins), by(period) chain sequentially detail" == 198
rcof "noi seldecomp y age sex, by(period) chain det" == 198
assert mreldif(e(b), B0) == 0 & "`e(chain)'" == "chain"
tempfile dtl
log using `"`dtl'"', text replace name(dlog)
capture noisily seldecomp y (a: age sex) (b: ins), by(period) chain sequentially fixothermargins detail
capture noisily seldecomp y age sex, by(period) chain det
capture noisily seldecomp y (a: age) (b: sex), by(period) chain individually
capture noisily seldecomp y (a: age) (b: sex), by(period) chain individually fixothermargins
capture noisily seldecomp y age sex, by(period) chain focalgroup(1999)
log close dlog
log_count `"`dtl'"' "detail may not be combined with chain"
assert r(n) == 2
* no advice to run the pair alone
log_count `"`dtl'"' "run the pair alone"
assert r(n) == 0
log_count `"`dtl'"' "requires sequentially"
assert r(n) == 0
log_count `"`dtl'"' "individually may not be combined with chain"
assert r(n) == 2
log_count `"`dtl'"' "focalgroup() may not be combined with chain"
assert r(n) == 1
log_count `"`dtl'"' "in each pair the lower"
assert r(n) == 0
log_count `"`dtl'"' "refgroup("
assert r(n) == 0
assert mreldif(e(b), B0) == 0 & "`e(chain)'" == "chain"
rcof "noi seldecomp y c.inc sex, by(period) chain contvars(inc 4, recut) bayesboot(reps(3))" == 198
rcof "noi seldecomp y (a: age) (b: sex), by(period) sequentially fixothermargins makeplot(, chainstep(1))" == 198
rcof "noi seldecomp y (a: age) (b: sex), by(period) chain sequentially fixothermargins makeplot(, chainstep(4))" == 198
rcof "noi seldecomp y (a: age) (b: sex), by(period) chain sequentially fixothermargins makeplot(stages)" == 198
rcof "noi seldecomp y (a: age) (b: sex), by(period) chain sequentially fixothermargins makeplot(waterfalllevels)" == 198
rcof "noi seldecomp y age sex if period == 1997, by(period) chain" == 459
decode period, gen(speriod)
rcof "noi seldecomp y age sex, by(speriod) chain" == 198
drop speriod
* at replay
seldecomp y (a: age) (b: sex), by(period) chain sequentially fixothermargins
rcof "noi seldecomp, makeplot(, chainstep(0))" == 198
rcof "noi seldecomp, makeplot(stages)" == 198
seldecomp y (a: age) (b: sex) if period <= 1999, by(period) sequentially fixothermargins
rcof "noi seldecomp, makeplot(, chainstep(1))" == 198
* a pair that fails stops the chain with its message; the results in memory are kept
seldecomp y age sex, by(period) chain
matrix B0 = e(b)
gen byte x6 = 1 + (reg > 2)
replace x6 = 3 if period == 2001 & mod(_n, 9) == 0
rcof "noi seldecomp y age x6, by(period) chain" == 459
assert mreldif(e(b), B0) == 0 & "`e(chain)'" == "chain"
* nothing is left in Mata: the pairs' notes of the pairs that ran; the same after a
* Break in pair 2 (the Mata raker's hook: the second raking stops)
mata: st_local("left", strofreal(findexternal("_seldecomp_chainnotes") != NULL))
assert `left' == 0
mata: seldecomp_brk = 1
rcof "noi seldecomp y age sex, by(period) chain" == 1
mata: st_local("left", strofreal(findexternal("_seldecomp_chainnotes") != NULL | findexternal("seldecomp_brk") != NULL))
assert `left' == 0
assert mreldif(e(b), B0) == 0
tempfile elg
log using `"`elg'"', text replace name(elog)
capture noisily seldecomp y age x6, by(period) chain
log close elog
log_count `"`elg'"' "raking support failure: these levels of the reference group (period ="
assert r(n) == 1
log_count `"`elg'"' "collapse or drop these levels"
assert r(n) == 1
log_count `"`elg'"' "rake the other group"
assert r(n) == 0
drop x6

di as res _n "PASSED: test_chain.do"

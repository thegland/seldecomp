* test_diagnostic_counts.do -- raw convergence and level-target diagnostic sidecars.
* Run from tests with the repository on adopath.
version 16
clear
set more off
set varabbrev off
discard

sysuse nlsw88, clear
quietly seldecomp wage race collgrad, by(married) sequentially
matrix C = e(raking_checks)
assert colsof(C) == 5 & rowsof(C) == 2
assert C[1,1] == 1 & C[2,1] == 1
assert C[1,2] == 0 & C[2,2] == 0 & C[1,3] == 0 & C[2,3] == 0
assert C[1,4] == 3 & C[2,4] == 5 & C[1,5] == 1 & C[2,5] == 1
capture confirm matrix e(raking_misses)
assert _rc != 0

* Two affected terms but four missed levels; only race is constrained at stage 1. The bound
* leaves targets unmet: rc 430 without force, and e(converged) = 0 with it.
rcof "noi seldecomp wage race collgrad, by(married) sequentially trimhirel(.8) trimfrequency(once)" == 430
quietly seldecomp wage race collgrad, by(married) sequentially ///
    trimhirel(.8) trimfrequency(once) generate(rw) force
matrix C = e(raking_checks)
matrix M = e(raking_misses)
assert C[1,1] == 1 & C[2,1] == 1
assert C[1,2] == 1 & C[1,3] == 2 & C[2,2] == 2 & C[2,3] == 4
assert rowsof(M) == 3 & colsof(M) == 6
assert missing(e(maxreldif)) & e(converged) == 0 & e(targets_met) == 0
local W1 = e(W_focal)
local W2 = e(W_ref)
forvalues m = 1/3 {
    local s = M[`m',1]
    local v `e(raking_miss`m'_var)'
    local l = M[`m',3] - ("`v'" == "collgrad")
    local rv = cond(`s' == 1, "rw_s1", "rw")
    quietly count if e(sample) & married == 1 & `v' == `l'
    local target = `W1' * r(N) / `W2'
    quietly summarize `rv' if e(sample) & married == 0 & `v' == `l', meanonly
    local achieved = r(sum)
    assert reldif(M[`m',4], `target') < 1e-12
    assert reldif(M[`m',5], `achieved') < 1e-12
    assert reldif(M[`m',6], reldif(`achieved', `target')) < 1e-12
}
assert strpos(`"`e(raking_miss1_level)'"', "White") > 0
assert "`e(raking_miss3_var)'" == "collgrad"
matrix unweighted_C = C
matrix unweighted_M = M

* Counts and discrepancies use the normalized control scale, while totals stay in user units.
generate double wt = 10
quietly seldecomp wage race collgrad [pw=wt], by(married) sequentially ///
    trimhirel(.8) trimfrequency(once) force
matrix C = e(raking_checks)
matrix M = e(raking_misses)
assert reldif(C[1,5], .1) < 1e-12 & reldif(C[2,5], .1) < 1e-12
forvalues m = 1/3 {
    assert reldif(M[`m',4], 10*unweighted_M[`m',4]) < 1e-12
    assert reldif(M[`m',5], 10*unweighted_M[`m',5]) < 1e-12
    assert reldif(M[`m',6], unweighted_M[`m',6]) < 1e-12
}

* A stop at the iteration limit and target failure are separate checks.
quietly seldecomp wage race collgrad, by(married) sequentially iterate(1) force
matrix C = e(raking_checks)
assert C[1,1] == 0 & C[1,2] == 0 & C[1,3] == 0
assert C[2,1] == 0 & C[2,2] == 1 & C[2,3] == 3

* Reverse terms are renumbered internally; stored term identities and labels stay original.
quietly seldecomp wage race collgrad, by(married) sequentially reverse ///
    trimhirel(.8) trimfrequency(once) force
matrix C = e(raking_checks_rev)
matrix M = e(raking_misses_rev)
assert C[1,4] == 2 & C[2,4] == 5
assert M[1,1] == 1 & M[1,2] == 2
assert "`e(raking_miss_rev1_var)'" == "collgrad"
quietly seldecomp wage (demo: race collgrad), by(married) sequentially detail ///
    trimhirel(.8) trimfrequency(once) force
matrix C = e(raking_checks_detail)
matrix M = e(raking_misses_detail)
assert rowsof(C) == 2 & C[1,4] == 3 & C[2,4] == 5
assert C[1,3] == 2 & C[2,3] == 4 & rowsof(M) == 3
assert "`e(raking_miss_detail3_var)'" == "collgrad"
* The detail shortcut copies the raw flag, not e(bystage)'s adjusted convergence.
quietly seldecomp wage race collgrad, by(married) sequentially detail ///
    trimhirel(.8) trimfrequency(once) force
assert mreldif(e(raking_checks), e(raking_checks_detail)) == 0
assert mreldif(e(raking_misses), e(raking_misses_detail)) == 0

* Held margins are active targets too, with their own focal shares.
quietly seldecomp wage race collgrad, by(married) sequentially fixothermargins ///
    trimhirel(.8) trimfrequency(once) force
matrix C = e(raking_checks)
matrix M = e(raking_misses)
assert C[1,4] == 5 & C[2,4] == 5
quietly count if e(sample) & married == 0 & collgrad == 0
local held = r(N)
assert M[2,1] == 1 & M[2,2] == 2 & M[2,3] == 1
assert reldif(M[2,4], `held') < 1e-12

* Each chain pair retains its sidecars and labels after subsequent pair estimates.
generate byte period = mod(_n, 3)
quietly seldecomp wage race collgrad, by(period) chain sequentially ///
    trimhirel(.8) trimfrequency(once) force
matrix C2 = e(raking_checks2)
matrix M2 = e(raking_misses2)
local pairvar `e(raking_miss2_1_var)'
local pairlevel `"`e(raking_miss2_1_level)'"'
quietly seldecomp wage race collgrad if period >= 1, by(period) sequentially ///
    trimhirel(.8) trimfrequency(once) force
assert mreldif(C2, e(raking_checks)) < 1e-12
assert mreldif(M2, e(raking_misses)) < 1e-12
assert "`pairvar'" == "`e(raking_miss1_var)'"
assert `"`pairlevel'"' == `"`e(raking_miss1_level)'"'

* The iterations stabilize on incompatible margins: force changes e(converged), not the raw flag.
clear
input byte(g a b n) double y
0 0 0 3 0
0 1 1 1 10
1 0 1 3 5
1 1 0 1 2
end
expand n
quietly seldecomp y a b, by(g) force
matrix C = e(raking_checks)
assert e(converged) == 0 & C[1,1] == 1
assert C[1,2] == 1 & C[1,3] == 2 & C[1,4] == 4

* Singleton control after support removal retains the original total and original level label.
clear
input byte(g x) double(y w)
0 2 0 1
0 2 10 10
0 1 20 1
1 2 0 1
1 2 10 1
end
quietly seldecomp y x [pw=w], by(g)
matrix C = e(raking_checks)
assert C[1,1] == 1 & C[1,3] == 0 & C[1,4] == 0
quietly seldecomp y x [pw=w], by(g) trimhiabs(2) trimfrequency(once) force
matrix C = e(raking_checks)
matrix M = e(raking_misses)
assert C[1,1] == 1 & C[1,2] == 1 & C[1,3] == 1 & C[1,4] == 1
assert reldif(C[1,5], 2/11) < 1e-12
assert M[1,3] == 2 & reldif(M[1,4], 12) < 1e-12
assert reldif(M[1,5], 12/11 + 2) < 1e-12
assert reldif(M[1,6], reldif(C[1,5]*M[1,5], C[1,5]*M[1,4])) < 1e-12
assert "`e(raking_miss1_level)'" == "2"

* Lower bounds retain an unsupported level: its zero target must be counted as missed.
quietly seldecomp y x [pw=w], by(g) trimloabs(1) force
matrix C = e(raking_checks)
matrix M = e(raking_misses)
assert C[1,2] == 1 & C[1,3] == 1 & C[1,4] == 2
assert M[1,3] == 1 & M[1,4] == 0 & M[1,5] == 1
assert "`e(raking_miss1_level)'" == "1"

di as result "test_diagnostic_counts.do PASSED"

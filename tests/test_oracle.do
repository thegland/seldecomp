* test_oracle.do -- seldecomp against an independent raking oracle written in Mata from the
* method description alone (no Python): 40 random configurations with two to four raking
* variables, an optional interaction, a level present only in the Focal Group, pweights, iweights or
* no weights, one raking or sequentially, with the blocks not yet added left free (the
* default) or held at the Focal Group's own shares (fixothermargins).  Explained parts must agree to
* 1e-7 and the generated raked weights to 1e-6.
*
* Run from the tests directory:   . do test_oracle.do

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

mata:
// one stage of the raking decomposition, written from the method description only:
// Focal Group (focal) weights are raked so that the variables in S match the Reference Group's weighted
// shares; with hold (fixothermargins), the other variables are held at the Focal Group's own
// shares, and without it they are free.  Focal Group obs
// whose level of a variable in S is absent from the Reference Group get zero weight (target share 0).
// returns (m1, m1r, explained, supportfail, converged) and fills rw
real rowvector orc_stage(real colvector y, real colvector g1, real colvector w,
    real matrix C, real rowvector S, real scalar hold, real colvector rw)
{
    real colvector A, B, keep, cur, tgt, f, wprev, c, wa
    real scalar W1, W2, k, v, L, l, it, m1, m1r, K, fail, sB, sA
    real matrix T
    A = selectindex(g1:==1); B = selectindex(g1:==0)
    W1 = sum(w[A]); W2 = sum(w[B])
    m1 = sum(w[A]:*y[A])/W1
    K = cols(C)
    keep = J(rows(A),1,1)
    fail = 0
    // zero-target levels
    for (k=1; k<=cols(S); k++) {
        v = S[k]; L = max(C[.,v])
        for (l=1; l<=L; l++) {
            sB = sum(w[B]:*(C[B,v]:==l))
            if (sB==0) keep = keep :* (C[A,v]:!=l)
        }
    }
    // support check on the remaining Focal Group obs
    for (k=1; k<=cols(S); k++) {
        v = S[k]; L = max(C[.,v])
        for (l=1; l<=L; l++) {
            sB = sum(w[B]:*(C[B,v]:==l))
            sA = sum(keep:*(C[A,v]:==l))
            if (sB>0 & sA==0) fail = 1
        }
    }
    rw = J(rows(y),1,.)
    if (fail) return((m1, ., ., 1, 0))
    wa = w[A]:*keep
    // targets: list of variables and target totals
    real rowvector vars
    vars = S
    if (hold) {
        for (v=1; v<=K; v++) if (!anyof(S, v)) vars = vars, v
    }
    real scalar Wk
    Wk = sum(wa)
    T = J(cols(vars), max(C), 0)
    for (k=1; k<=cols(vars); k++) {
        v = vars[k]; L = max(C[.,v])
        for (l=1; l<=L; l++) {
            if (anyof(S, v)) T[k,l] = W1*sum(w[B]:*(C[B,v]:==l))/W2
            else             T[k,l] = W1*sum(wa:*(C[A,v]:==l))/Wk
        }
    }
    for (it=1; it<=20000; it++) {
        wprev = wa
        for (k=1; k<=cols(vars); k++) {
            v = vars[k]; L = max(C[.,v]); c = C[A,v]
            for (l=1; l<=L; l++) {
                cur = sum(wa:*(c:==l))
                if (cur>0) wa = wa :* (1 :+ (c:==l)*(T[k,l]/cur - 1))
            }
        }
        if (max(abs(wa-wprev):/(wprev:+(wprev:==0))) < 1e-13) break
    }
    m1r = sum(wa:*y[A])/sum(wa)
    rw[A] = wa
    // converged: every target is met to 1e-6 (relative); with inconsistent margins IPF's
    // per-cycle change can decay below the tolerance while the margins stay unmet
    real scalar mx, ach
    mx = 0
    for (k=1; k<=cols(vars); k++) {
        v = vars[k]; L = max(C[.,v])
        for (l=1; l<=L; l++) {
            ach = sum(wa:*(C[A,v]:==l))
            if (T[k,l] > 0) mx = max((mx, abs(ach - T[k,l]) / T[k,l]))
        }
    }
    return((m1, m1r, m1-m1r, 0, mx < 1e-6))
}
end

tempname post
tempfile results
postfile `post' int(rep) byte(nvar nstage hold inter zsup ywt) long(n) double(maxdiff_expl maxdiff_rw maxdiff_fac) ///
    int(rc) byte(orcfail oconv) using `results', replace

set seed 424242
local R 40
forvalues rep = 1/`R' {
    clear
    local n = runiformint(200, 3000)
    qui set obs `n'
    local K = runiformint(2, 4)
    local inter = (runiform() < 0.25)
    local zsup  = (runiform() < 0.35)
    local seq   = (runiform() < 0.6)
    local hold  = `seq' & (runiform() < 0.5)
    local ywt   = runiformint(1, 3)        // 1 pw, 2 iw, 3 unweighted
    qui gen byte grp = runiform() < runiform(0.2, 0.8)
    forvalues k = 1/`K' {
        local L`k' = runiformint(2, 6)
        // group-dependent distribution, correlated with the previous variable
        local prev 0
        if `k' > 1 local prev 0.3*x`=`k'-1'/`L`=`k'-1''
        qui gen double u`k' = runiform() + 0.35*grp*runiform(-1,1) + `prev'
        qui gen byte x`k' = 1 + floor(`L`k'' * (u`k' - (-0.35)) / (1.3 + 0.35 + 0.35))
        qui replace x`k' = max(1, min(`L`k'', x`k'))
    }
    qui gen double y = 0
    forvalues k = 1/`K' {
        qui replace y = y + runiform(-1, 1) * x`k' + 0.5*(x`k' == 2)*(x1 == 1)
    }
    qui replace y = y + 0.8*grp + rnormal()
    if runiform() < 0.3 qui replace y = (y > 0)
    // MEPS-like weights: 449 to 51,662
    qui gen double wt = min(51662, max(449, exp(rnormal(8.5, 0.9))))
    // guarantee support: rows 1-36 are Focal Group (grp 0), rows 37-72 Reference Group, each with every level
    qui replace grp = 0 in 1/36
    qui replace grp = 1 in 37/72
    qui gen int ii = mod(_n - 1, 36)
    qui replace x1 = mod(ii, `L1') + 1 in 1/72
    qui replace x2 = mod(floor(ii / `L1'), `L2') + 1 in 1/72
    forvalues k = 3/`K' {
        qui replace x`k' = mod(ii, `L`k'') + 1 in 1/72
    }
    drop ii
    // a level present only in the Focal Group (focal = grp==0 by default), never in the support rows
    if `zsup' {
        qui replace x1 = 9 if grp == 0 & runiform() < 0.03 & _n > 72
        qui count if grp == 0 & x1 == 9
        if r(N) == 0 qui replace x1 = 9 in 73
        qui replace grp = 0 in 73
    }
    // spec
    local spec
    local codes
    if `inter' {
        local spec x1#x2
        qui egen long c1 = group(x1 x2)
        local codes c1
        local first 3
    }
    else local first 1
    forvalues k = `first'/`K' {
        local spec `spec' x`k'
        qui egen long c`k' = group(x`k')
        local codes `codes' c`k'
    }
    local nv : word count `codes'
    local wopt
    if `ywt' == 1 local wopt [pw=wt]
    if `ywt' == 2 local wopt [iw=wt]
    if `ywt' == 3 qui replace wt = 1
    local sopt
    if `seq' local sopt = cond(`hold', "sequentially fixothermargins", "sequentially")
    capture drop rw*
    capture noisily seldecomp y `spec' `wopt', by(grp) `sopt' generate(rw)
    local rc = _rc
    local nst = cond(`seq', `nv', 1)
    local md_e = .
    local md_w = .
    local md_f = .
    local of = 0
    local oc = 1
    // the oracle is run whenever seldecomp produced results or refused to (430: unmet targets),
    // so that the two verdicts on feasibility can be compared
    if inlist(`rc', 0, 430) {
        if `rc' == 0 matrix S = e(bystage)
        local md_e 0
        local md_w 0
        local md_f 0
        forvalues s = 1/`nst' {
            if `seq' {
                local Sset
                forvalues j = 1/`s' {
                    local Sset `Sset' `j'
                }
            }
            else {
                local Sset
                forvalues j = 1/`nv' {
                    local Sset `Sset' `j'
                }
            }
            local Sset : subinstr local Sset " " ",", all
            mata: st_view(Y=., ., "y"); st_view(G=., ., "grp"); st_view(W=., ., "wt"); st_view(CC=., ., "`codes'")
            mata: rwo = .; r = orc_stage(Y, 1:-G, W, CC, (`Sset'), `hold', rwo); st_numscalar("o_expl", r[3]); st_numscalar("o_fail", r[4]); st_numscalar("o_conv", r[5]); st_store(., st_addvar("double", "orw`s'"), rwo)
            if o_fail == 1 local of 1
            if o_fail == 0 & o_conv == 0 local oc 0
            if `rc' != 0 continue
            if o_fail == 0 {
                local d = abs(S[`s', 2] - o_expl) / max(1, abs(o_expl))
                local md_e = max(`md_e', `d')
                if `s' == `nst' local rwv rw
                else local rwv rw_s`s'
                // weights: the ado's generated weight vs the oracle (Focal Group only)
                qui gen double dd = reldif(`rwv', orw`s') if grp == 0
                qui sum dd, meanonly
                local md_w = max(`md_w', r(max))
                drop dd
            }
        }
    }
    post `post' (`rep') (`nv') (`nst') (`hold') (`inter') (`zsup') (`ywt') (`n') (`md_e') (`md_w') (`md_f') (`rc') (`of') (`oc')
}
postclose `post'
use `results', clear
tab rc
sum maxdiff_expl maxdiff_rw
* seldecomp fails exactly when the oracle finds a support failure (459) or margins that never
* settle (430, unmet targets), and agrees otherwise
list rep nvar nstage hold inter zsup n rc orcfail oconv if rc != 0, clean
assert (rc == 0 & orcfail == 0 & oconv == 1) | (rc == 459 & orcfail == 1) | (rc == 430 & oconv == 0)
count if rc == 0
assert r(N) >= 30
assert maxdiff_expl < 1e-7 if rc == 0
assert maxdiff_rw < 1e-6 if rc == 0

di as res _n "PASSED: test_oracle.do"

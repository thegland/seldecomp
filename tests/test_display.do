* test_display.do -- terminal formats, numeric precision, labels after prefixes,
* diagnostics visibility, stored notes, saved transcripts, stage equations,
* confidence intervals, width/wrapping, literal labels, and replication messages.
* Clean fits show a raking status and factor summary by default (and the thin-cell count
* when it is above zero).
* diagnostics requests the base table; failed raking checks show it automatically,
* with failure counts and per-variable detail. Stage and covariate summaries remain
* above the coefficient and optional bootstrap-CI tables.
* Output is captured in text logs and searched; the data are Stata's nlsw88 and small
* synthetic sets. With the global SELDECOMP_REGEN_TRANSCRIPTS set to 1, section 5 rewrites
* the reference transcripts instead of comparing with them (tests/README.md).
*
* Run from the tests directory:   . do test_display.do

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

* ---- helpers --------------------------------------------------------------------
capture program drop log_count
program define log_count, rclass
    // log_count "file" "text"  --  r(n) = number of lines of the file that contain text
    args fn txt
    local n 0
    mata: st_local("n", strofreal(sum(strpos(cat(st_local("fn")), st_local("txt")) :> 0)))
    return scalar n = `n'
end

capture program drop log_rows
program define log_rows, rclass
    // log_rows "file" "text"  --  r(n) lines contain text; r(row1), r(row2), ... are those lines
    args fn txt
    local n 0
    mata: _sd_l = cat(st_local("fn"))
    // Log metadata can contain a Windows short path such as THOMAS~1; it is not output.
    mata: _sd_l = select(_sd_l, !ustrregexm(_sd_l, "^ *(name|log|log type|opened on|closed on):") :& ///
        substr(_sd_l, 1, 2) :!= ". " :& substr(_sd_l, 1, 2) :!= "> ")
    mata: _sd_l = select(_sd_l, strpos(_sd_l, st_local("txt")) :> 0)
    mata: st_local("n", strofreal(rows(_sd_l)))
    forvalues i = 1/`n' {
        mata: st_local("r", _sd_l[`i'])
        return local row`i' `"`r'"'
    }
    mata: mata drop _sd_l
    return scalar n = `n'
end

capture program drop log_width
program define log_width, rclass
    // log_width "file"  --  r(max) = widest line between the markers SDW_W_BEGIN and
    // SDW_W_END, not counting command echoes (". " and "> ") and the marker commands;
    // r(tilde) = the number of those lines that contain "~" (an abbreviation)
    args fn
    local w 0
    local t 0
    mata: _sd_l = cat(st_local("fn"))
    mata: _sd_b = select(1::rows(_sd_l), _sd_l :== "SDW_W_BEGIN")
    mata: _sd_e = select(1::rows(_sd_l), _sd_l :== "SDW_W_END")
    mata: _sd_l = _sd_l[(_sd_b[1] + 1)::(_sd_e[1] - 1)]
    mata: _sd_l = select(_sd_l, substr(_sd_l, 1, 2) :!= ". " :& substr(_sd_l, 1, 2) :!= "> ")
    mata: st_local("w", strofreal(max(udstrlen(_sd_l))))
    mata: st_local("t", strofreal(sum(strpos(_sd_l, "~") :> 0)))
    mata: mata drop _sd_l _sd_b _sd_e
    return scalar max = `w'
    return scalar tilde = `t'
end

capture program drop file_width
program define file_width, rclass
    // file_width "file"  --  r(max) = widest line of the file
    args fn
    local w 0
    mata: st_local("w", strofreal(max(udstrlen(cat(st_local("fn"))) \ 0)))
    return scalar max = `w'
end

capture program drop log_msgs
program define log_msgs, rclass
    // log_msgs "file"  --  r(n) = the messages of seldecomp between the markers SDW_M_BEGIN
    // and SDW_M_END, counted by their first line (a note, a warning, a footnote "(# of", an
    // error), not counting command echoes (". " and "> "); r(pfx) = the prefixes' own notes
    // and warnings ("Note:", "Warning:")
    args fn
    local n 0
    local p 0
    mata: _sd_l = cat(st_local("fn"))
    mata: _sd_b = select(1::rows(_sd_l), _sd_l :== "SDW_M_BEGIN")
    mata: _sd_e = select(1::rows(_sd_l), _sd_l :== "SDW_M_END")
    mata: _sd_l = _sd_l[(_sd_b[1] + 1)::(_sd_e[1] - 1)]
    mata: _sd_l = select(_sd_l, substr(_sd_l, 1, 2) :!= ". " :& substr(_sd_l, 1, 2) :!= "> ")
    mata: st_local("n", strofreal(sum(regexm(_sd_l, ///
        "^(note:|warning:|[(][0-9]+ of |by[(]|bayesboot[(][)]:|raking |one of the groups|insupport)"))))
    mata: st_local("p", strofreal(sum(regexm(_sd_l, "^(Note:|Warning:)"))))
    mata: mata drop _sd_l _sd_b _sd_e
    return scalar n = `n'
    return scalar pfx = `p'
end

capture program drop log_after
program define log_after, rclass
    // log_after "file" "text" k  --  r(line) = the line k lines after (k < 0: before) the
    // first line of the file that contains text
    args fn txt k
    mata: _sd_l = cat(st_local("fn"))
    mata: _sd_i = select(1::rows(_sd_l), strpos(_sd_l, st_local("txt")) :> 0)
    mata: st_local("line", _sd_l[_sd_i[1] + `k'])
    mata: mata drop _sd_l _sd_i
    return local line `"`line'"'
end

capture program drop sep_above
program define sep_above
    // sep_above "file" "text"  --  the first line of the file that contains text starts a
    // section: above it one blank line, the rule of 78 "=" and one blank line (Display_sep),
    // above which the output goes on
    args fn txt
    local eq = "=" * 78
    foreach k in 1 2 3 4 {
        log_after `"`fn'"' `"`txt'"' -`k'
        local want = cond(`k' == 2, "`eq'", cond(`k' == 4, "<nonblank>", ""))
        local got = strrtrim(`"`r(line)'"')
        if ("`want'" == "<nonblank>" & `"`got'"' == "") | ("`want'" != "<nonblank>" & `"`got'"' != "`want'") {
            di as err `"`fn': line -`k' from the one with `txt' is |`r(line)'|, expected |`want'|"'
            exit 9
        }
    }
end

capture program drop log_index
program define log_index, rclass
    // log_index "file" "text"  --  r(i) = the number of the first line of the file that
    // contains text, 0 when none does
    args fn txt
    local i 0
    mata: _sd_l = cat(st_local("fn"))
    mata: _sd_i = select(1::rows(_sd_l), strpos(_sd_l, st_local("txt")) :> 0) \ 0
    mata: st_local("i", strofreal(_sd_i[1]))
    mata: mata drop _sd_l _sd_i
    return scalar i = `i'
end

capture program drop log_para
program define log_para, rclass
    // log_para "file" "text"  --  r(text) = the first line of the file that contains text,
    // joined by one blank with the continuation lines after it (nonblank lines that start
    // with a blank), each stripped of its blanks; r(lines) = the number of lines joined
    args fn txt
    local n 0
    mata: _sd_l = cat(st_local("fn"))
    mata: _sd_n = .
    mata: st_local("t", td_para(_sd_l, st_local("txt"), _sd_n))
    mata: st_local("n", strofreal(_sd_n))
    mata: mata drop _sd_l _sd_n
    return local text `"`t'"'
    return scalar lines = `n'
end

capture program drop tr_compare
program define tr_compare
    // tr_compare "new" "reference"  --  the two text files have the same lines once trailing
    // blanks are removed; otherwise the first difference is shown and the program fails
    args new ref
    local i 0
    mata: _sd_a = strrtrim(cat(st_local("new")))
    mata: _sd_b = strrtrim(cat(st_local("ref")))
    mata: st_local("i", strofreal(td_first_diff(_sd_a, _sd_b)))
    if `i' > 0 {
        mata: st_local("la", td_line(_sd_a, `i'))
        mata: st_local("lb", td_line(_sd_b, `i'))
        mata: mata drop _sd_a _sd_b
        di as err `"transcript `ref' differs at line `i':"'
        di as err `"  new:       `macval(la)'"'
        di as err `"  reference: `macval(lb)'"'
        exit 9
    }
    mata: mata drop _sd_a _sd_b
end
capture mata: mata drop td_first_diff() td_line() td_para()
mata:
real scalar td_first_diff(string colvector a, string colvector b)
{
    real scalar i, n
    n = min((rows(a), rows(b)))
    for (i = 1; i <= n; i++) {
        if (a[i] != b[i]) return(i)
    }
    if (rows(a) != rows(b)) return(n + 1)
    return(0)
}
string scalar td_line(string colvector a, real scalar i)
{
    if (i > rows(a)) return("(end of file)")
    return(a[i])
}
string scalar td_para(string colvector a, string scalar txt, real scalar n)
{
    real scalar i, k
    string scalar t
    n = 0
    k = 0
    for (i = 1; i <= rows(a); i++) {
        if (strpos(a[i], txt)) {
            k = i
            break
        }
    }
    if (k == 0) return("")
    t = strrtrim(a[k])
    n = 1
    while (k < rows(a)) {
        if (substr(a[k + 1], 1, 1) != " " | strtrim(a[k + 1]) == "") break
        k++
        n++
        t = t + " " + strtrim(a[k])
    }
    return(t)
}
end

tempfile lg
local lgf `"`lg'.log"'
capture log close dlog

* ==============================================================================
* 1. formats: set cformat, cformat(), pformat(), sformat(), invalid formats
* ==============================================================================
di as txt _n "---- 1. formats"
sysuse nlsw88, clear
quietly seldecomp wage i.collgrad i.south, by(married)
local g1_def = strtrim(string(e(mean_focal), "%9.0g"))
local g1_3f  = strtrim(string(e(mean_focal), "%9.3f"))
assert "`g1_def'" != "`g1_3f'"
* set cformat changes the printed focal_group (reset at once: the setting outlives this file)
log using `"`lgf'"', text replace name(dlog)
set cformat %9.3f
capture noisily seldecomp wage i.collgrad i.south, by(married)
local rc = _rc
set cformat
log close dlog
assert `rc' == 0
log_rows `"`lgf'"' "      Focal mean |"
assert r(n) == 1
assert strtrim(substr(`"`r(row1)'"', 19, .)) == "`g1_3f'"
* without it, Stata's default %9.0g (8.080765, as oaxaca and ereturn display print it)
log using `"`lgf'"', text replace name(dlog)
seldecomp
log close dlog
log_rows `"`lgf'"' "      Focal mean |"
assert strtrim(substr(`"`r(row1)'"', 19, .)) == "`g1_def'"
* the options, at estimation and at replay; a format too wide is replaced by the default
log using `"`lgf'"', text replace name(dlog)
seldecomp wage i.collgrad i.south, by(married) cformat(%9.3f)
seldecomp, cformat(%9.3f)
seldecomp, cformat(%12.3f)
log close dlog
log_rows `"`lgf'"' "      Focal mean |"
assert r(n) == 3
assert strtrim(substr(`"`r(row1)'"', 19, .)) == "`g1_3f'"
assert strtrim(substr(`"`r(row2)'"', 19, .)) == "`g1_3f'"
assert strtrim(substr(`"`r(row3)'"', 19, .)) == "`g1_def'"
log_count `"`lgf'"' "note: invalid cformat(), using default"
assert r(n) == 1
rcof "noi seldecomp, cformat(abc)" == 7
rcof "noi seldecomp wage i.collgrad, by(married) pformat(%q)" == 7
* pformat() and sformat() with the bootstrap columns: P>|z| and z follow them
log using `"`lgf'"', text replace name(dlog)
seldecomp wage i.collgrad i.south, by(married) bayesboot(reps(20) rseed(1) nodots) ///
 pformat(%5.2f) sformat(%8.3f)
log close dlog
log_rows `"`lgf'"' "      Difference |"
local se = sqrt(el(e(V), 4, 4))
local z = _b[overall:difference] / `se'
local p = 2 * normal(-abs(`z'))
assert strtrim(substr(`"`r(row1)'"', 42, 8)) == strtrim(string(`z', "%8.3f"))
assert strtrim(substr(`"`r(row1)'"', 51, 5)) == strtrim(string(`p', "%5.2f"))

* ==============================================================================
* 2. labels after a prefix: svy jackknife (t, design df) and the Inference line
* ==============================================================================
di as txt _n "---- 2. prefix labels"
sysuse nlsw88, clear
keep in 1/300
gen double one = 1
svyset idcode [pw=one]
quietly svy jackknife, nodots: seldecomp wage i.collgrad, by(married)
assert e(df_r) < .
log using `"`lgf'"', text replace name(dlog)
svy
seldecomp
log close dlog
* the six rows of svy's own table and of seldecomp's redisplay are identical: the same
* coefficient, std. err., t, P>|t| and confidence interval
local native focal_group raked_focal ref_group difference explained unexplained
local shown `""Focal mean" "Raked focal mean" "Reference mean" "Difference" "Explained" "Unexplained""'
forvalues j = 1/6 {
    local n : word `j' of `native'
    local t : word `j' of `shown'
    log_rows `"`lgf'"' " `n' |"
    assert r(n) == 1
    local nr = stritrim(substr(`"`r(row1)'"', strpos(`"`r(row1)'"', "|") + 1, .))
    log_rows `"`lgf'"' "`t' |"
    assert r(n) == 1
    local sr = stritrim(substr(`"`r(row1)'"', strpos(`"`r(row1)'"', "|") + 1, .))
    assert `"`nr'"' == `"`sr'"'
}
log_count `"`lgf'"' "P>|t|"
assert r(n) == 2
log_count `"`lgf'"' "Inference: Jackknife"
assert r(n) == 1
* at another level the redisplay uses t with the design df
log using `"`lgf'"', text replace name(dlog)
seldecomp, level(90)
log close dlog
log_rows `"`lgf'"' "      Focal mean |"
local crit = invttail(e(df_r), .05)
local ll = _b[overall:focal_group] - `crit' * _se[overall:focal_group]
assert reldif(real(word(`"`r(row1)'"', 8)), `ll') < 1e-6
* seldecomp's own bootstrap: the Inference line is the same at estimation and on replay
sysuse nlsw88, clear
log using `"`lgf'"', text replace name(dlog)
seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ///
    bayesboot(reps(10) rseed(1) nodots)
seldecomp
log close dlog
log_count `"`lgf'"' "Inference: Bayesian bootstrap"
assert r(n) == 2
log_count `"`lgf'"' "Estimate  std. err."
assert r(n) == 2
log_count `"`lgf'"' "Observed"
assert r(n) == 0
log_count `"`lgf'"' "b~p"
assert r(n) == 0
* Stata 19's bayesboot prefix: the same label at estimation (when e(vcetype) is still
* "Bootstrap") and on replay; the two-word variance type is split, never abbreviated
if c(stata_version) >= 19 {
    log using `"`lgf'"', text replace name(dlog)
    bayesboot, reps(5) rseed(1) nodots: seldecomp wage i.collgrad, by(married)
    seldecomp
    log close dlog
    log_count `"`lgf'"' "Inference: Bayesian bootstrap"
    assert r(n) == 2
    log_count `"`lgf'"' "Estimate  std. err."
    assert r(n) == 2
    log_count `"`lgf'"' "b~p"
    assert r(n) == 0
}

* ==============================================================================
* 3. diagnostics: optional for clean fits and automatic when a raking check fails
* ==============================================================================
di as txt _n "---- 3. diagnostics"
sysuse nlsw88, clear
log using `"`lgf'"', text replace name(dlog)
seldecomp wage i.collgrad i.south, by(married)
log close dlog
log_count `"`lgf'"' "Raking diagnostics"
assert r(n) == 0
log_count `"`lgf'"' "Raking factors"
assert r(n) == 1
foreach t in "Max rel. chg" "Effective N (Kish)" "Raked weights" {
    log_count `"`lgf'"' "`t'"
    assert r(n) == 0
}
* The same compact table is used for weighted/unweighted fits and replay.
gen double wt = 1 + grade / 20
foreach w in "" "[pw=wt]" {
    log using `"`lgf'"', text replace name(dlog)
    seldecomp wage i.collgrad i.south `w', by(married) diagnostics
    seldecomp, diagnostics
    log close dlog
    foreach t in "Raking diagnostics" "Zero wt: zero raked weights." "Thin target levels (" {
        log_count `"`lgf'"' "`t'"
        assert r(n) == 2
    }
    foreach t in "Raking failure counts" "Raking factors" "Raked weights" {
        log_count `"`lgf'"' "`t'"
        assert r(n) == 0
    }
}
* Unattainable held margins: diagnostics and failure counts appear automatically.
preserve
clear
set obs 400
set seed 4
gen byte grp = _n > 200
gen byte a = cond(grp, mod(_n, 10) < 3, mod(_n, 2))
gen byte b = cond(grp, mod(floor(_n / 3), 10) < 6, 1 - mod(_n, 2))
gen double yy = a + b + rnormal()
log using `"`lgf'"', text replace name(dlog)
seldecomp yy (a) (b), by(grp) sequentially fixothermargins force
log close dlog
restore
assert e(converged) == 0 & e(targets_met) == 0
foreach t in "Raking diagnostics" "Raking failure counts" "*not*" {
    log_count `"`lgf'"' "`t'"
    assert r(n) >= 1
}
foreach t in "Raking factors" "note: stages 1-2:" {
    log_count `"`lgf'"' "`t'"
    assert r(n) == 0
}
* Stored notes remain available even though target-failure header notes are suppressed.
assert e(N_notes) == 1
assert strpos(`"`e(note1)'"', "stages 1-2: 1 raking margin not met (worst: ") == 1
assert strpos(`"`e(note1)'"', "); results kept because force was specified") > 0

* ==============================================================================
* 4. definitions and notes
* ==============================================================================
di as txt _n "---- 4. definitions and notes"
sysuse nlsw88, clear
* the definitions are always shown, at estimation and at replay, in the one order of the
* difference, focal group - reference group (e(diffsign))
log using `"`lgf'"', text replace name(dlog)
seldecomp wage i.collgrad i.south, by(married)
seldecomp
log close dlog
assert "`e(diffsign)'" == "focal_group - ref_group"
foreach t in "Difference = focal - reference. Explained = focal - raked focal." ///
    "Focal Group (raked): married = 0 (Single)" "Reference Group:     married = 1 (Married)" {
    log_count `"`lgf'"' "`t'"
    assert r(n) == 2
}
foreach t in "difference: focal_group" "explained: focal_group" "unexplained: raked_focal" {
    log_count `"`lgf'"' "`t'"
    assert r(n) == 0
}
* help example 2: the thin target level race = 3 (8 obs) is raked at all three stages and is
* reported once (two note lines, not one note per stage; the definition of a thin target level is in
* the raking diagnostics' footnote, not in the note)
seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) ///
    sequentially fixothermargins
assert e(N_notes) == 2
assert `"`e(note1)'"' == "thin target level, stages 1-3:"
assert `"`e(note2)'"' == "    race = 3 (Other) (8 obs)"
assert e(N_focal_mincell) == 1
* several thin target levels in one note, counted; no stage text in a one-stage run
seldecomp wage i.industry, by(married) mincell(20)
assert e(N_notes) == 5
assert `"`e(note1)'"' == "4 thin target levels:"
* a level absent from the reference group from stage 2 on: one zero-weight note with its stages; a cell
* whose count changes between stages (zero-weight observations leave the raking sample) is
* listed under each of its stage ranges
gen byte race2 = race
replace race2 = 9 if married == 0 & _n <= 30
count if married == 0 & race2 == 9
local nz = r(N)
seldecomp wage (a: collgrad) (b: race2) (c: south), by(married) sequentially fixothermargins mincell(300)
assert `"`e(note1)'"' == "level absent from the reference group (zero raked weight), stages 2-3:"
assert `"`e(note2)'"' == "    race2 = 9 (`nz' obs)"
assert `"`e(note3)'"' == "thin target level, stage 1:"
assert `"`e(note5)'"' == "2 thin target levels, stages 2-3:"
assert e(N_notes) == 7

* ==============================================================================
* 5. reference transcripts: the runs of display_cases.do, line by line
* ==============================================================================
di as txt _n "---- 5. reference transcripts"
* display_cases.do writes one text log per case; each must equal its file in
* reference/transcripts up to trailing blanks. With the global SELDECOMP_REGEN_TRANSCRIPTS
* set to 1 the logs are written into reference/transcripts instead (tests/README.md).
local tref "reference/transcripts"
if "$SELDECOMP_REGEN_TRANSCRIPTS" == "1" {
    capture mkdir "reference"
    capture mkdir "`tref'"
    local tdir "`tref'"
    di as txt "(regenerating the reference transcripts in `tref')"
}
else {
    tempfile tdd
    local tdir `"`tdd'_d"'
    mkdir `"`tdir'"'
}
do display_cases.do `"`tdir'"'
local want : dir "`tref'" files "*.txt"
local got : dir `"`tdir'"' files "*.txt"
local want : list sort want
local got : list sort got
assert `: word count `want'' == 10
assert `"`got'"' == `"`want'"'
foreach f of local want {
    tr_compare `"`tdir'/`f'"' `"`tref'/`f'"'
    file_width `"`tdir'/`f'"'
    assert r(max) <= 78
    * none of these words in the output
    foreach w in varbyvar isolated nohold holdothers refminusfocal refgroup ipfraking raker {
        log_count `"`tdir'/`f'"' "`w'"
        assert r(n) == 0
    }
}

* ==============================================================================
* 6. the explained share, the stage equations (with inference only), the CI table, the stage
*    table; on the transcripts of section 5 and a few further runs
* ==============================================================================
di as txt _n "---- 6. share, stage rows, citable, stage table"
local t1  `"`tdir'/01_default.txt"'
local t2  `"`tdir'/02_weighted.txt"'
local t3  `"`tdir'/03_sequentially_fixed.txt"'
local t4  `"`tdir'/04_sequentially_fixed_bs.txt"'
local t5  `"`tdir'/05_replay.txt"'
local t6  `"`tdir'/06_bca_cluster.txt"'
local t7  `"`tdir'/07_sequentially.txt"'
local t10 `"`tdir'/10_zero_share.txt"'
local hl = "-" * 78

* The share retains its std. err. and percentile CI even when the difference's CI includes 0.
log_rows `"`t4'"' "      Difference |"
assert r(n) == 1
assert real(word(`"`r(row1)'"', 7)) < 0 & real(word(`"`r(row1)'"', 8)) > 0
log_rows `"`t4'"' "Explained share of the difference:"
assert r(n) == 1
assert regexm(strtrim(`"`r(row1)'"'), "^Explained share of the difference: -?[0-9]+[.][0-9]%, bootstrap std[.] err[.] [0-9]+[.][0-9]$")
log_count `"`t4'"' "unstable"
assert r(n) == 0
log_count `"`t4'"' "percentile CI"
assert r(n) == 1
* its replay at level(90) (transcript 5): the 90% CI of the difference excludes 0, but the
* share's std. err. exceeds 100 (heavy tails): the std. err. and the
* percentile interval, without comment
log_rows `"`t5'"' "      Difference |"
assert real(word(`"`r(row1)'"', 7)) > 0
log_rows `"`t5'"' "Explained share of the difference:"
assert r(n) == 1
assert regexm(strtrim(`"`r(row1)'"'), "^Explained share of the difference: -?[0-9]+[.][0-9]%, bootstrap std[.] err[.] [0-9]+[.][0-9]$")
log_count `"`t5'"' "    90% percentile CI: "
assert r(n) == 1
log_count `"`t5'"' "heavy-tailed"
assert r(n) == 0
log_count `"`t5'"' "unstable"
assert r(n) == 0
* a stable share: std. err. and percentile CI, one decimal and % units
sysuse nlsw88, clear
log using `"`lgf'"', text replace name(dlog)
seldecomp wage i.race i.south, by(collgrad) bayesboot(reps(30) rseed(1) nodots)
log close dlog
assert abs(_b[overall:difference] / _se[overall:difference]) > invnormal(.975)
assert e(pct_explained_se) < 100
local sh = strtrim(string(e(pct_explained), "%9.1f"))
local se = strtrim(string(e(pct_explained_se), "%9.1f"))
local pl = strtrim(string(el(e(pct_explained_ci), 1, 1), "%9.1f"))
local pu = strtrim(string(el(e(pct_explained_ci), 1, 2), "%9.1f"))
foreach x in sh pl pu {
    if real("``x''") == 0 local `x' 0.0
}
log_count `"`lgf'"' "Explained share of the difference: `sh'%, bootstrap std. err. `se'"
assert r(n) == 1
log_count `"`lgf'"' "    95% percentile CI: `pl'% to `pu'%"
assert r(n) == 1
log_count `"`lgf'"' "unstable"
assert r(n) == 0
* never -0.0: a share of -0.04% (transcript 10) is 0.0%, in the share line and the stage table
preserve
clear
set obs 2000
gen byte grp = _n > 1000
gen byte x = cond(grp, _n <= 1498, _n <= 502)
gen byte z = mod(_n, 2)
gen double y = 10 * grp + x
quietly seldecomp y (a: x) (b: z), by(grp) sequentially fixothermargins
assert e(pct_explained) < 0 & e(pct_explained) > -.05
restore
log_count `"`t10'"' "-0.0"
assert r(n) == 0
log_count `"`t10'"' "Explained share of the difference: 0.0%"
assert r(n) == 1
log_rows `"`t10'"' "    1 a    "
assert r(n) == 1
assert substr(`"`r(row1)'"', 55, 6) == "   0.0" & substr(`"`r(row1)'"', 73, 6) == "   0.0"

* the stage equations in the coefficient table only with inference: without it
* the stage table holds their values, and the table has the overall equation alone, at
* estimation and at replay (transcript 3 and lgf); with it every equation, titled, in the
* order marg_expl, cum_expl, cum_unexpl, raked_focal (transcript 4; transcript 5, the replay
* with citable: in the coefficient table and in the CI table)
log using `"`lgf'"', text replace name(dlog)
seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) ///
    sequentially fixothermargins
seldecomp
log close dlog
local eqrows `""By stage         | Marginal explained" "By stage         | Cumulative explained""'
local eqrows `"`eqrows' "By stage         | Cumulative unexplained" "By stage         | Raked focal group mean""'
foreach t of local eqrows {
    log_count `"`lgf'"' "`t'"
    assert r(n) == 0
    log_count `"`t3'"' "`t'"
    assert r(n) == 0
    log_count `"`t4'"' "`t'"
    assert r(n) == 1
    log_count `"`t5'"' "`t'"
    assert r(n) == 2
}
log_count `"`lgf'"' "Overall          |"
assert r(n) == 2
log_count `"`t3'"' "Overall          |"
assert r(n) == 1
* transcript 4 (two blocks): after the overall equation's six rows and a separator, each
* equation's title row, its two rows and a separator
local k 0
foreach t of local eqrows {
    log_after `"`t4'"' "Overall          |" `=8 + 4 * `k''
    assert `"`r(line)'"' == "`t'"
    local ++k
}

* the stage table of help example 2 (transcript 3), against e(bystage) of the same run
matrix S = e(bystage)
log_count `"`t3'"' "                                -------- Cumulative -------- ---- Marginal ---"
assert r(n) == 1
log_count `"`t3'"' "Stage Block        Raked mean   Explained Unexplained      %  Explained      %"
assert r(n) == 1
forvalues s = 1/3 {
    local nm : word `s' of `e(blocks)'
    log_rows `"`t3'"' "    `s' `nm'   "
    assert r(n) == 1
    local row `"`r(row1)'"'
    assert strtrim(substr(`"`row'"', 22, 9)) == strtrim(string(S[`s', 1], "%9.0g"))
    assert strtrim(substr(`"`row'"', 33, 9)) == strtrim(string(S[`s', 2], "%9.0g"))
    assert strtrim(substr(`"`row'"', 45, 9)) == strtrim(string(S[`s', 3], "%9.0g"))
    assert strtrim(substr(`"`row'"', 55, 6)) == strtrim(string(S[`s', 4], "%6.1f"))
    assert strtrim(substr(`"`row'"', 63, 9)) == strtrim(string(S[`s', 5], "%9.0g"))
    assert strtrim(substr(`"`row'"', 73, 6)) == strtrim(string(S[`s', 27], "%6.1f"))
}
* the footnote is below the closing rule, and a blank line under it (the main table follows):
* one line, the same under sequentially fixothermargins (transcript 3) and under sequentially
* alone (transcript 7), without the Pylypchuk-Selden references and without a sentence on
* the blocks not yet added
foreach t in 3 7 {
    log_after `"`t`t''"' "%: percent of the difference." -1
    assert `"`r(line)'"' == "`hl'"
    log_after `"`t`t''"' "%: percent of the difference." 0
    assert `"`r(line)'"' == "%: percent of the difference."
    log_after `"`t`t''"' "%: percent of the difference." 1
    assert strtrim(`"`r(line)'"') == ""
}
foreach t in 3 7 {
    log_count `"`t`t''"' "Pylypchuk and Selden 2008"
    assert r(n) == 0
    log_count `"`t`t''"' "Pylypchuk-Selden procedure"
    assert r(n) == 0
}

* the CI table only with citable; no pointer line otherwise (the help documents citable)
foreach t in 4 6 {
    log_count `"`t`t''"' "Bayesian bootstrap confidence intervals"
    assert r(n) == 0
    log_count `"`t`t''"' "intervals: seldecomp, citable"
    assert r(n) == 0
}
log_count `"`t5'"' "Bayesian bootstrap confidence intervals"
assert r(n) == 1
log_count `"`t5'"' "intervals: seldecomp, citable"
assert r(n) == 0
* the replay at level(90): both tables at 90%
log_count `"`t5'"' "[90% conf. interval]"
assert r(n) == 2
* citable at estimation (with bca); the replay without it has no CI table and no pointer
sysuse nlsw88, clear
log using `"`lgf'"', text replace name(dlog)
seldecomp wage i.collgrad i.south, by(married) bayesboot(reps(20) rseed(1) cluster(industry) ///
    bca nodots citable)
seldecomp
log close dlog
log_count `"`lgf'"' "Bayesian bootstrap confidence intervals"
assert r(n) == 1
log_count `"`lgf'"' " (BCa)"
assert r(n) == 6
log_count `"`lgf'"' "intervals: seldecomp, citable"
assert r(n) == 0

* ==============================================================================
* 7. width: the display stays within 78 columns at any linesize
* ==============================================================================
di as txt _n "---- 7. width"
sysuse nlsw88, clear
* a value label longer than a line: covariatedetail wraps it
label define lrace3 1 "White: a label that is much longer than the 78 columns of the output, to be wrapped" ///
    2 "Black" 3 "Other"
gen byte race3 = race
label values race3 lrace3
gen byte race2 = race
replace race2 = 9 if married == 0 & _n <= 30
local ls0 = c(linesize)
set linesize 200
log using `"`lgf'"', text replace name(dlog)
di "SDW_W_BEGIN"
seldecomp wage i.collgrad i.south c.ttl_exp, by(married) contvars(ttl_exp 4)
seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ///
    bayesboot(reps(10) rseed(1) nodots citable) diagnostics
seldecomp, level(90) citable
seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) ///
    sequentially fixothermargins
seldecomp wage i.industry, by(married) mincell(20) covariatedetail
seldecomp wage i.race3 i.south, by(married) covariatedetail
seldecomp wage (a: collgrad) (b: race2) (c: south), by(married) sequentially fixothermargins mincell(300)
seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) individually fixothermargins
seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially reverse
quietly jackknife, cluster(industry) nodots: seldecomp wage i.collgrad i.south if industry < ., by(married)
seldecomp
* replicates kept at the iterate() limit (the synthetic data and raking options of
* test_bootstrap.do, 9b2), and failed replicates and a failed BCa deletion (the margins
* of 9b that some replicates cannot meet)
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
recast double y x3 wt
seldecomp y x1 x2 x4 [pw=wt], by(group) iterate(11) tolerance(1e-10) ///
    bayesboot(reps(20) rseed(7) nodots)
seldecomp y x1 x2 [pw=wt] in 1/150, by(group) iterate(9) tolerance(1e-14) ///
    bayesboot(reps(10) rseed(3) bca nodots)
clear
set obs 400
gen byte group = _n > 200
gen int k = cond(group == 0, _n, _n - 200)
gen byte a = cond(group == 0, k > 130, cond(k <= 30, 0, cond(k <= 54, 1, !mod(k, 2))))
gen byte b = cond(group == 0, inrange(k, 61, 130), cond(k <= 30, 0, cond(k <= 54, 1, mod(k, 2))))
gen double y = a + 2 * b + mod(_n * 7, 11) / 10
gen int cl = cond(group == 1 & k <= 10, 11, 1 + mod(k, 10))
seldecomp y a b, by(group) bayesboot(reps(20) rseed(1) nodots bca cluster(cl))
* a margin not met, kept with force (the synthetic data of section 3)
clear
set obs 400
set seed 4
gen byte grp = _n > 200
gen byte a = cond(grp, mod(_n, 10) < 3, mod(_n, 2))
gen byte b = cond(grp, mod(floor(_n / 3), 10) < 6, 1 - mod(_n, 2))
gen double yy = a + b + rnormal()
seldecomp yy (a) (b), by(grp) sequentially fixothermargins force
di "SDW_W_END"
log close dlog
set linesize `ls0'
log_width `"`lgf'"'
di as txt "widest display line: " r(max)
assert r(max) <= 78
* the runs above reach the paragraphs in question: one note on the replications, none on the
* jackknife deletions, the force note and no warning; the point estimate's note of
* the iterate() limit counts the replications that stopped there
foreach t in "replications failed and were dropped: in them, raking" "(change of support" ///
    "White: a label that is much longer" ///
    "note: 8 of 20 replications stopped at the iterate() limit" ///
    "note: raking stopped at the iterate() limit, also in 10 of 10 replications;" {
    log_count `"`lgf'"' "`t'"
    assert r(n) == 1
}
log_count `"`lgf'"' "Raking failure counts"
assert r(n) == 2
foreach t in "more than 5%" "jackknife deletion" "warning:" "kept replications stopped" ///
    "note: 10 of 10 replications" "(raking failed or missed" {
    log_count `"`lgf'"' "`t'"
    assert r(n) == 0
}
log_count `"`lgf'"' "Bayesian bootstrap confidence intervals"
assert r(n) == 2
* the stage rows of the bootstrapped staged run: its coefficient and CI tables, at estimation
* and at the replay at level(90)
log_count `"`lgf'"' "By stage         | Cumulative explained"
assert r(n) == 4
* the binned variables are listed under the raking variables (or the stages)
log_count `"`lgf'"' "  (ttl_exp: 4 quantile bins, cut points from the pooled sample)"
assert r(n) == 1
log_count `"`lgf'"' "  (age: 5 quantile bins, cut points from the pooled sample)"
assert r(n) == 1

* ==============================================================================
* 8. long names and labels in the header: a left-column
*    entry that does not fit (47 columns) leaves its slot blank, and the line follows the
*    rows in full, as a paragraph beyond 78 columns; nothing is truncated or abbreviated
*    (the tables are shown too: the width check covers them)
* ==============================================================================
di as txt _n "---- 8. long names and labels"
local ls0 = c(linesize)
set linesize 200
* (a) a 16-character by() variable with 20-character value labels, at estimation and at
* replay: the group lines below the rows, Focal obs and Reference obs with the left column blank
sysuse nlsw88, clear
rename married married_status_x
label define mlong 0 "Never married single" 1 "Currently married xx"
label values married_status_x mlong
log using `"`lgf'"', text replace name(dlog)
di "SDW_W_BEGIN"
seldecomp wage i.collgrad i.south, by(married_status_x)
seldecomp
di "SDW_W_END"
log close dlog
log_width `"`lgf'"'
assert r(max) <= 78 & r(tilde) == 0
foreach t in "Focal Group (raked): married_status_x = 0 (Never married single)" ///
    "Reference Group:     married_status_x = 1 (Currently married xx)" {
    log_rows `"`lgf'"' "`t'"
    assert r(n) == 2
    assert strrtrim(`"`r(row1)'"') == "`t'" & strrtrim(`"`r(row2)'"') == "`t'"
}
foreach t in "Focal obs" "Reference obs" {
    log_rows `"`lgf'"' "`t'"
    assert r(n) == 2
    forvalues i = 1/2 {
        assert strtrim(substr(`"`r(row`i')'"', 1, 48)) == "" & substr(`"`r(row`i')'"', 49, strlen("`t'")) == "`t'"
    }
}
log_after `"`lgf'"' "Reference obs" 1
assert substr(`"`r(line)'"', 1, 21) == "Focal Group (raked): "
* without inference no Inference line: the raking variables follow the group lines
log_after `"`lgf'"' "Reference Group:     " 1
assert substr(`"`r(line)'"', 1, 18) == "Raking variables: "
* (b) a by() name longer than 16 characters is not abbreviated; it stays in the rows when
* its line fits
rename married_status_x married_status_xyz20
label values married_status_xyz20
log using `"`lgf'"', text replace name(dlog)
di "SDW_W_BEGIN"
seldecomp wage i.collgrad i.south, by(married_status_xyz20)
di "SDW_W_END"
log close dlog
log_width `"`lgf'"'
assert r(max) <= 78 & r(tilde) == 0
log_rows `"`lgf'"' "Focal Group (raked): married_status_xyz20 = 0"
assert r(n) == 1
assert substr(`"`r(row1)'"', 49, 9) == "Focal obs"
* (c) a string by() with long values: one line at full width, and one wrapped under the
* by() variable (hanging indent 21), braces included
sysuse nlsw88, clear
local v1 "A string group value with {braces} that is long enough to be wrapped over two lines"
local v2 "Northeast region"
gen str90 mstr_long_name16 = cond(married, "`v2'", "`v1'")
log using `"`lgf'"', text replace name(dlog)
di "SDW_W_BEGIN"
seldecomp wage i.collgrad i.south, by(mstr_long_name16)
di "SDW_W_END"
log close dlog
assert `"`e(focal_group_label)'"' == "`v1'"
log_width `"`lgf'"'
assert r(max) <= 78 & r(tilde) == 0
log_para `"`lgf'"' "Focal Group (raked):"
assert r(lines) == 2
assert `"`r(text)'"' == "Focal Group (raked): mstr_long_name16 = `v1'"
log_after `"`lgf'"' "Focal Group (raked):" 1
assert strlen(`"`r(line)'"') - strlen(strltrim(`"`r(line)'"')) == 21
log_after `"`lgf'"' "Focal Group (raked):" 2
assert strrtrim(`"`r(line)'"') == "Reference Group:     mstr_long_name16 = `v2'"
* (d) a 32-character by() variable with the labels of (a), wrapped, and a 32-character
* weight variable: Focal wgt sum and 2 keep their rows; the group lines and then the
* weights line follow them, then the Inference line with its replications
sysuse nlsw88, clear
rename married b_this_is_a_32_character_name_xx
label define mlong 0 "Never married single" 1 "Currently married xx"
label values b_this_is_a_32_character_name_xx mlong
gen double wt_a_very_long_weight_variable_x = 1 + (age - 34) / 20
log using `"`lgf'"', text replace name(dlog)
di "SDW_W_BEGIN"
seldecomp wage i.collgrad i.south [pw=wt_a_very_long_weight_variable_x], ///
    by(b_this_is_a_32_character_name_xx) bayesboot(reps(5) rseed(1) nodots cluster(industry)) 
di "SDW_W_END"
log close dlog
log_width `"`lgf'"'
assert r(max) <= 78 & r(tilde) == 0
log_para `"`lgf'"' "Focal Group (raked):"
assert r(lines) == 2
assert `"`r(text)'"' == "Focal Group (raked): b_this_is_a_32_character_name_xx = 0 (Never married single)"
log_para `"`lgf'"' "Reference Group:"
assert r(lines) == 2
assert `"`r(text)'"' == "Reference Group:     b_this_is_a_32_character_name_xx = 1 (Currently married xx)"
foreach t in "Focal obs" "Reference obs" "Focal wgt sum" "Reference wgt sum" {
    log_rows `"`lgf'"' "`t'"
    assert r(n) == 1
    assert strtrim(substr(`"`r(row1)'"', 1, 48)) == ""
}
log_after `"`lgf'"' "Reference wgt sum" 1
assert substr(`"`r(line)'"', 1, 21) == "Focal Group (raked): "
log_rows `"`lgf'"' "Weights: "
assert r(n) == 1
assert strrtrim(`"`r(row1)'"') == "Weights: pweight = wt_a_very_long_weight_variable_x"
log_after `"`lgf'"' "Weights: " 1
assert substr(`"`r(line)'"', 1, 29) == "Inference: Bayesian bootstrap" & substr(`"`r(line)'"', 49, 12) == "Replications"
log_count `"`lgf'"' "(Replications based on 12 clusters in industry)"
assert r(n) == 1
* (e) a 32-character depvar and denominator with ratio(), a 32-character raking variable,
* and the insupport paragraph: the paragraphs wrap, the names are whole
sysuse nlsw88, clear
gen double y_this_is_a_32_character_name_xx = wage
gen double d_this_is_a_32_character_name_xx = hours
gen byte r_this_is_a_32_character_name_xx = south
gen byte race2 = race
replace race2 = 9 if married == 0 & _n <= 30
log using `"`lgf'"', text replace name(dlog)
di "SDW_W_BEGIN"
seldecomp y_this_is_a_32_character_name_xx i.collgrad i.r_this_is_a_32_character_name_xx, ///
    by(married) ratio(d_this_is_a_32_character_name_xx)
seldecomp y_this_is_a_32_character_name_xx i.collgrad i.race2, by(married) insupport ///
    ratio(d_this_is_a_32_character_name_xx)
di "SDW_W_END"
log close dlog
* the header abbreviates nothing; the one "~" of each run is the coefficient table's, whose
* heading gives the outcome 12 columns (as Stata's own tables do)
log_width `"`lgf'"'
assert r(max) <= 78 & r(tilde) == 2
log_rows `"`lgf'"' "~"
assert r(n) == 2
forvalues i = 1/2 {
    assert substr(`"`r(row`i')'"', 1, 18) == abbrev("y_this_is_a_32_character_name_xx", 16) + " |"
}
log_para `"`lgf'"' "Ratio of means:"
assert r(lines) == 2
assert `"`r(text)'"' == "Ratio of means: y_this_is_a_32_character_name_xx / d_this_is_a_32_character_name_xx"
log_para `"`lgf'"' "Raking variables: collgrad r_"
assert `"`r(text)'"' == "Raking variables: collgrad r_this_is_a_32_character_name_xx"
log_para `"`lgf'"' "Support (insupport):"
assert strpos(`"`r(text)'"', "; focal group ratio of y_this_is_a_32_character_name_xx before the exclusion: ")
* the group lines of by(married) are unchanged
log_count `"`lgf'"' "Focal Group (raked): married = 0 (Single)       Focal obs"
assert r(n) == 2
* (f) an Inference line longer than the left column (no prefix sets such a variance type;
* e(vcetype) is overwritten before the replay): Replications keeps its row
quietly bootstrap, reps(3) seed(1) nodots: seldecomp wage i.collgrad, by(married)
capture program drop sdw_vcetype
program define sdw_vcetype, eclass
    ereturn local vcetype "A variance type long enough not to fit the left column"
end
sdw_vcetype
log using `"`lgf'"', text replace name(dlog)
di "SDW_W_BEGIN"
seldecomp
di "SDW_W_END"
log close dlog
program drop sdw_vcetype
log_width `"`lgf'"'
assert r(max) <= 78 & r(tilde) == 0
log_rows `"`lgf'"' "Replications"
assert r(n) == 1
assert strtrim(substr(`"`r(row1)'"', 1, 48)) == ""
log_after `"`lgf'"' "Replications" 1
assert strrtrim(`"`r(line)'"') == "Inference: A variance type long enough not to fit the left column"
set linesize `ls0'

* ==============================================================================
* 9. the header's rows, rounding error shown as 0, the covariatedetail columns, the line size
* ==============================================================================
di as txt _n "---- 9. header rows, rounding error, covariatedetail, line size"
* (a) the Focal Group line follows Number of obs (there is no Raking method row:
* seldecomp's own raker is the only one), no Inference line without inference, and the full
* diagnostics are titled "Raking diagnostics" alone
sysuse nlsw88, clear
log using `"`lgf'"', text replace name(dlog)
seldecomp wage i.collgrad i.south, by(married) diagnostics
seldecomp, diagnostics
log close dlog
foreach t in "Raking method" "Raking diagnostics (" {
    log_count `"`lgf'"' "`t'"
    assert r(n) == 0
}
log_rows `"`lgf'"' "Raking diagnostics"
assert r(n) == 2
assert strrtrim(`"`r(row1)'"') == "Raking diagnostics" & strrtrim(`"`r(row2)'"') == "Raking diagnostics"
log_count `"`lgf'"' "Inference:"
assert r(n) == 0
log_after `"`lgf'"' "Number of obs" 1
assert substr(`"`r(line)'"', 1, 21) == "Focal Group (raked): "
* (b) a value smaller in magnitude than 1e-13 times the largest of the outcome's
* focal_group, raked_focal and ref_group is displayed as 0; e() keeps it. The data:
* a constant outcome, whose block contributions are rounding error; an
* outcome with means of 1e7, whose x2 contribution (-9.76e-06, 5e-13 of the means) is real
* (it stays under tolerance(1e-14)) and is shown; an outcome with means of 1e-8 (no floor);
* and a difference that is rounding error (focal group half 1/6 and half 1/2, reference group 1/3: 1e-16,
* and a share of -8.5e+15% in e()), shown as 0, with no share and no percentages
clear
set obs 400
set seed 4
gen byte g = _n > 200
gen byte x = 1 + (runiform() < cond(g, .7, .3)) + (runiform() < cond(g, .7, .3))
gen byte x2 = runiform() < .5
gen double yc = 5
gen double yb = 1e7 * (1 + g) + 1e5 * x + rnormal()
gen double ys = 1e-8 * (1 + g + x)
gen double yd = cond(g, 1/3, cond(mod(_n, 2), 1/6, 1/2))
log using `"`lgf'"', text replace name(dlog)
seldecomp yc x x2, by(g) sequentially
log close dlog
assert _b[marg_expl:x] != 0 & abs(_b[marg_expl:x]) < 1e-14
assert abs(_b[cum_unexpl:x]) < 1e-14
* (without inference the coefficient table has no stage rows: the stage table shows them)
foreach t in "               x |" "              x2 |" {
    log_rows `"`lgf'"' "`t'"
    assert r(n) == 0
}
foreach s in 1 2 {
    log_rows `"`lgf'"' "    `s' x"
    assert r(n) == 1
    assert word(`"`r(row1)'"', 4) == "0" & word(`"`r(row1)'"', 5) == "0" & word(`"`r(row1)'"', 7) == "0"
}
log_count `"`lgf'"' "e-1"
assert r(n) == 0
* with inference the stage rows are in the coefficient table, rounding error shown as 0 (the
* rows of x in marg_expl, cum_expl and cum_unexpl; raked_focal's is the mean, 5)
log using `"`lgf'"', text replace name(dlog)
seldecomp yc x x2, by(g) sequentially bayesboot(reps(5) rseed(1) nodots)
log close dlog
log_rows `"`lgf'"' "               x |"
assert r(n) == 4
forvalues i = 1/3 {
    assert strtrim(substr(`"`r(row`i')'"', 19, .)) == "0"
}
assert strtrim(substr(`"`r(row4)'"', 19, .)) == "5"
log using `"`lgf'"', text replace name(dlog)
seldecomp yb x x2, by(g) sequentially
log close dlog
assert abs(_b[marg_expl:x2]) > 1e-13 * 2e7 & abs(_b[marg_expl:x2]) < 1e-12 * 2e7
log_rows `"`lgf'"' "    2 x2"
assert r(n) == 1
assert word(`"`r(row1)'"', 7) == strtrim(string(_b[marg_expl:x2], "%9.0g"))
log using `"`lgf'"', text replace name(dlog)
seldecomp ys x x2, by(g) sequentially
log close dlog
log_rows `"`lgf'"' "       Explained |"
assert strtrim(substr(`"`r(row1)'"', 19, .)) == strtrim(string(_b[overall:explained], "%9.0g"))
assert abs(_b[overall:explained]) < 1e-7
log using `"`lgf'"', text replace name(dlog)
seldecomp yd x x2, by(g) sequentially fixothermargins
log close dlog
assert e(difference) != 0 & abs(e(difference)) < 1e-15 & e(pct_explained) < .
log_rows `"`lgf'"' "      Difference |"
assert strtrim(substr(`"`r(row1)'"', 19, .)) == "0"
log_count `"`lgf'"' "Explained share of the difference: n/a (the difference is 0)"
assert r(n) == 1
forvalues s = 1/2 {
    log_rows `"`lgf'"' "    `s' x"
    assert substr(`"`r(row1)'"', 55, 6) == "     ." & substr(`"`r(row1)'"', 73, 6) == "     ."
}
* the bootstrap tables: rounding error in the coefficients, bias, std. errors and limits is
* 0; a std. err. that is 0 leaves the row with its coefficient only
log using `"`lgf'"', text replace name(dlog)
seldecomp yc x x2, by(g) sequentially fixothermargins bayesboot(reps(10) rseed(1) nodots citable) 
log close dlog
log_count `"`lgf'"' "e-1"
assert r(n) == 0
* (x's rows: marg_expl, cum_expl, cum_unexpl and raked_focal in each table; the first of the
* CI table is row 5)
log_rows `"`lgf'"' "               x |"
assert r(n) == 8
assert strtrim(substr(`"`r(row1)'"', 19, .)) == "0"
assert real(word(`"`r(row5)'"', 3)) == 0 & real(word(`"`r(row5)'"', 6)) == 0
* individually fixothermargins: the contributions table
log using `"`lgf'"', text replace name(dlog)
seldecomp yc (a: x) (b: x2), by(g) individually fixothermargins
log close dlog
log_count `"`lgf'"' "e-1"
assert r(n) == 0
log_rows `"`lgf'"' "Remainder "
assert strtrim(substr(`"`r(row1)'"', 10, 36)) == "0"
* The base table displays factor CVs, keeping numerical residuals in e(bystage).
log using `"`lgf'"', text replace name(dlog)
seldecomp yc x x2, by(g) sequentially fixothermargins diagnostics
log close dlog
matrix S = e(bystage)
assert S[1, 9] < 1e-13
log_count `"`lgf'"' "Zero wt: zero raked weights."
assert r(n) == 1
log_count `"`lgf'"' "Raking failure counts"
assert r(n) == 0
* (c) covariatedetail's columns: observations and weighted shares of the focal group and
* reference group, the focal group share under the final raked weights, and the gap (focal group
* minus reference group, the order of the difference, never -0.0000); the cells of an interaction by their value
* labels joined by #, and the bins of a c. variable; at estimation and at replay
sysuse nlsw88, clear
local ls0 = c(linesize)
set linesize 200
log using `"`lgf'"', text replace name(dlog)
di "SDW_W_BEGIN"
seldecomp wage i.race#i.collgrad c.age i.south, by(married) contvars(age 3) covariatedetail 
matrix C1 = e(covariates)
seldecomp, covariatedetail
matrix C2 = e(covariates)
di "SDW_W_END"
log close dlog
set linesize `ls0'
log_width `"`lgf'"'
assert r(max) <= 78 & r(tilde) == 0
log_count `"`lgf'"' "                             ----- Obs. ----  --------- Weighted share -------"
assert r(n) == 2
log_count `"`lgf'"' "Level                          Focal    Ref.    Focal    Ref.   Raked      Gap"
assert r(n) == 2
local labs `""White # Not college grad" "White # College grad" "Black # Not college grad""'
local labs `"`labs' "Black # College grad" "Other # Not college grad" "Other # College grad""'
local labs `"`labs' "<= 37" "(37, 41]" "> 41" "0 (Not south)" "1 (South)""'
assert rowsof(C1) == 11
forvalues r = 1/11 {
    local lb : word `r' of `labs'
    log_rows `"`lgf'"' "  `lb'  "
    assert r(n) == 2
    forvalues i = 1/2 {
        local row `"`r(row`i')'"'
        assert substr(`"`row'"', 1, 2 + strlen("`lb'")) == "  `lb'"
        * (string() drops the blanks that pad a number to its format's width)
        assert strtrim(substr(`"`row'"', 30, 7)) == string(C`i'[`r', 4], "%7.0fc")
        assert strtrim(substr(`"`row'"', 38, 7)) == string(C`i'[`r', 5], "%7.0fc")
        assert strtrim(substr(`"`row'"', 47, 7)) == string(C`i'[`r', 1], "%7.4f")
        assert strtrim(substr(`"`row'"', 55, 7)) == string(C`i'[`r', 2], "%7.4f")
        assert strtrim(substr(`"`row'"', 63, 7)) == string(C`i'[`r', 3], "%7.4f")
        assert substr(`"`row'"', 29, 1) == " " & substr(`"`row'"', 70, 1) == " "
        local gap = C`i'[`r', 1] - C`i'[`r', 2]
        if abs(`gap') < .00005 local gap 0
        assert strtrim(substr(`"`row'"', 71, 8)) == string(`gap', "%8.4f")
        assert strlen(strrtrim(`"`row'"')) == 78
        assert !strpos(`"`row'"', "-0.0000")
    }
}
log_para `"`lgf'"' "Focal: focal group."
assert `"`r(text)'"' == "Focal: focal group. Ref.: reference group. Raked: focal group share under the final raked weights. Gap: Focal minus Reference share (the order of the difference)."
log_count `"`lgf'"' "raked_focal:"
assert r(n) == 0
log_count `"`lgf'"' "Gap: Reference minus"
assert r(n) == 0
log_count `"`lgf'"' "race=1"
assert r(n) == 0
* (d) below a line size of 78 Stata wraps the rows; at any line size there is no note about it
foreach ls in 60 77 78 80 {
    set linesize `ls'
    log using `"`lgf'"', text replace name(dlog)
    seldecomp wage i.collgrad, by(married)
    seldecomp
    log close dlog
    log_rows `"`lgf'"' "78 columns"
    assert r(n) == 0
    log_rows `"`lgf'"' "Raking decomposition (Pylypchuk-Selden)"
    assert r(n) == 2
}
set linesize `ls0'

* ==============================================================================
* 10. one message per cause: the messages between markers are counted (log_msgs), at
*     estimation and at replay
* ==============================================================================
di as txt _n "---- 10. one message per cause"
clear
set seed 20260927
set obs 600
gen byte group = runiform() < .5
gen byte x1 = 1 + (runiform() < .4) + (runiform() < .3 + .3 * group)
gen byte x2 = runiform() < .6 - .2 * group
gen byte x5 = cond(runiform() < .85, x2, 1 - x2)
gen double y = 2 + x1 + x2 + .8 * group + rnormal()
gen double wt = 50 + 100 * runiform()
gen byte x4 = 1 + (runiform() < .5)
replace x4 = 9 if group == 0 & _n <= 40 & mod(_n, 8) == 0
replace x4 = 8 if _n > 500 & mod(_n, 11) == 0
replace x4 = 7 if group == 0 & inrange(_n, 41, 60) & mod(_n, 10) == 0
set linesize 80
* (a) bootstrap: whose replicates miss their targets at the iterate() limit, where the point
* estimate stopped too: its note counts them and gives the remedy, the one message at
* estimation and at replay; with iterate(50) no replicate fails
local ipf iterate(2) ctrltolerance(.037)
log using `"`lgf'"', text replace name(dlog)
di "SDW_M_BEGIN"
bootstrap, reps(30) seed(3) nodots: seldecomp y i.x2 i.x5 i.x1, by(group) `ipf'
di "SDW_M_END"
log close dlog
local nm = e(N_misreps)
assert inrange(`nm', 2, 28)
log_msgs `"`lgf'"'
assert r(n) == 1
log_para `"`lgf'"' "note: raking stopped at the iterate() limit"
assert `"`r(text)'"' == "note: raking stopped at the iterate() limit, also in `nm' of 30 bootstrap replicates (short of the targets, dropped); raise iterate()"
assert r(lines) == 2
foreach t in "Collapse sparse" "warning:" "failed and w" {
    log_count `"`lgf'"' "`t'"
    assert r(n) == 0
}
log using `"`lgf'"', text replace name(dlog)
di "SDW_M_BEGIN"
seldecomp
di "SDW_M_END"
log close dlog
log_msgs `"`lgf'"'
assert r(n) == 1
log_para `"`lgf'"' "note: raking stopped at the iterate() limit"
assert strpos(`"`r(text)'"', ", also in `nm' of 30 bootstrap replicates (short of the targets, dropped)")
quietly bootstrap, reps(30) seed(3) nodots: seldecomp y i.x2 i.x5 i.x1, by(group) iterate(50) ctrltolerance(.037)
assert e(N_misreps) == 0
* (b) the same under svy brr, whose output hides the notes: one warning after the loop names
* the cause, without a remedy; at replay the note counts them
mata: _sd_h = (1); for (i = 1; i <= 3; i++) _sd_h = (_sd_h, _sd_h \ _sd_h, -_sd_h); st_matrix("HAD8", _sd_h)
mata: mata drop _sd_h i
gen int stratum = 1 + mod(_n, 8)
gen byte ps2 = 1 + (mod(floor(_n / 8), 3) > 0)
forvalues r = 1/8 {
    gen double brr`r' = 2 * wt * ((HAD8[`r', stratum] == 1) == (ps2 == 1))
}
svyset ps2 [pw=wt], strata(stratum) brrweight(brr1-brr8)
log using `"`lgf'"', text replace name(dlog)
di "SDW_M_BEGIN"
svy brr, nodots: seldecomp y i.x2 i.x5 i.x1, by(group) `ipf'
di "SDW_M_END"
log close dlog
local nm = e(N_misreps)
assert inrange(`nm', 1, 7)
log_msgs `"`lgf'"'
assert r(n) == 1
log_para `"`lgf'"' "warning: "
assert strpos(`"`r(text)'"', "in " + cond(`nm' == 1, "it", "them") + ", raking stopped at the iterate() limit short of its targets.")
assert substr(`"`r(text)'"', -strlen("short of its targets."), .) == "short of its targets."
assert !strpos(`"`r(text)'"', "aise iterate") & !strpos(`"`r(text)'"', "ollapse")
log using `"`lgf'"', text replace name(dlog)
di "SDW_M_BEGIN"
seldecomp
di "SDW_M_END"
log close dlog
log_msgs `"`lgf'"'
assert r(n) == 1
log_para `"`lgf'"' "note: raking stopped at the iterate() limit"
assert strpos(`"`r(text)'"', ", also in `nm' of 8 BRR replicates (short of the targets)")
svyset, clear
drop stratum ps2 brr*
matrix drop HAD8
* (c) bayesboot() at the iterate() limit, where the point estimate stopped too: its note counts
* the replications that stopped there and those dropped, the one message at estimation and at
* replay
foreach c in "60 .037" "20 .3" {
    gettoken R tol : c
    local tol = strtrim("`tol'")
    log using `"`lgf'"', text replace name(dlog)
    di "SDW_M_BEGIN"
    seldecomp y i.x2 i.x5 i.x1 [pw=wt], by(group) iterate(2) ctrltolerance(`tol') ///
        bayesboot(reps(`R') rseed(3) nodots)
    di "SDW_M_END"
    log close dlog
    local nm = e(N_misreps)
    local ni = e(N_iterlimit)
    assert `nm' + `ni' == `R'
    if `R' == 60 assert `nm' > 0
    else assert `nm' == 0
    log_msgs `"`lgf'"'
    assert r(n) == 1
    local dtx = cond(`nm' > 0, " (`nm' short of the targets, dropped)", "")
    local caution = cond(`nm' > 0, "; bootstrap results may be incorrect.", "")
    log_para `"`lgf'"' "note: raking stopped at the iterate() limit"
    assert `"`r(text)'"' == "note: raking stopped at the iterate() limit, also in `R' of `R' replications`dtx'; raise iterate()`caution'"
    log using `"`lgf'"', text replace name(dlog)
    di "SDW_M_BEGIN"
    seldecomp
    di "SDW_M_END"
    log close dlog
    log_msgs `"`lgf'"'
    assert r(n) == 1
}
* (d) trimming that leaves margins unmet (rc 430 without force), kept by
* force: the replications at the iterate() limit are the trimming's (their margins unmet, as
* the point estimate's), not counted in e(N_iterlimit) and not said to meet the targets, but
* counted in the one note of the replications that force kept; diagnostic tables in place of
* a trimming header note
rcof "noi seldecomp y i.x1 i.x2 i.x4 [pw=wt], by(group) trimhirel(1.1) bayesboot(reps(20) rseed(2) nodots)" == 430
log using `"`lgf'"', text replace name(dlog)
di "SDW_M_BEGIN"
seldecomp y i.x1 i.x2 i.x4 [pw=wt], by(group) trimhirel(1.1) force ///
    bayesboot(reps(20) rseed(2) nodots generate(ff))
di "SDW_M_END"
log close dlog
assert e(N_iterlimit) == 0 & e(N_misreps) == 0 & e(badcontrols) > 0 & e(N_forced) == 20
log_count `"`lgf'"' "note: trimming left "
assert r(n) == 0
log_count `"`lgf'"' "note: force kept 20 of 20 replications that missed raking targets"
assert r(n) == 1
log_count `"`lgf'"' "replications"
assert r(n) == 1
log_count `"`lgf'"' "met the raking targets"
assert r(n) == 0
* each replication rerun on its weights: those at the limit have margins unmet
local nlim 0
forvalues r = 1/20 {
    quietly seldecomp y i.x1 i.x2 i.x4 [pw=ff_`r'], by(group) trimhirel(1.1) force
    if e(converged) == 0 {
        assert e(badcontrols) > 0
        local ++nlim
    }
}
assert `nlim' > 0
drop ff_*
* (e) replicates that fail before the raking (a resample of clusters without one group): one
* warning after the loop that names the cause, nothing per replicate under nodots; jackknife:
* likewise
clear
set seed 20260925
set obs 384
gen double wt = 50 + 100 * runiform()
gen byte group = runiform() < .5
gen byte x1 = 1 + (runiform() < .4) + (runiform() < .3 + .3 * group)
keep in 1/96
gen byte cl4 = 1 + mod(_n, 4)
replace group = (cl4 > 2)
gen double y = x1 + group + mod(_n, 7) / 7
log using `"`lgf'"', text replace name(dlog)
di "SDW_M_BEGIN"
bootstrap, reps(40) cluster(cl4) seed(3) nodots: seldecomp y i.x1, by(group)
di "SDW_M_END"
log close dlog
local nm = e(N_misreps)
assert inrange(`nm', 1, 20)
log_msgs `"`lgf'"'
assert r(n) == 1
log_count `"`lgf'"' "takes only one value"
assert r(n) == 0
local wr = cond(`nm' == 1, "was", "were")
log_para `"`lgf'"' "warning: "
assert `"`r(text)'"' == "warning: `nm' of 40 bootstrap replicates failed and `wr' dropped: in " + cond(`nm' == 1, "it", "them") + ", by(group) took only one value. Resample within the groups: bootstrap, strata(group)."
gen byte cl3 = cond(group == 1, 3, 1 + mod(_n, 2))
log using `"`lgf'"', text replace name(dlog)
di "SDW_M_BEGIN"
jackknife, cluster(cl3) nodots: seldecomp y i.x1, by(group)
di "SDW_M_END"
log close dlog
assert e(N_misreps) == 1
log_msgs `"`lgf'"'
assert r(n) == 1
log_para `"`lgf'"' "warning: "
assert `"`r(text)'"' == "warning: 1 of 3 jackknife replicates failed and was dropped: in it, by(group) took only one value."
* (f) bayesboot() with force and a stage margin no replication can meet: force keeps the
* replications as it keeps the point estimate: the point estimate's stored note and one
* displayed note counting kept replications, e(N_forced); replay retains only the replication
* footnote alongside the failure diagnostics
clear
set obs 400
gen byte grp = _n > 200
gen byte a = cond(grp, mod(_n, 10) < 3, mod(_n, 2))
gen byte b = 1 - a
gen double yy = a + b + mod(_n, 7) / 7
log using `"`lgf'"', text replace name(dlog)
di "SDW_M_BEGIN"
capture noisily seldecomp yy (A: a) (B: b), by(grp) sequentially fixothermargins force ///
    bayesboot(reps(20) rseed(1) nodots)
local rc = _rc
di "SDW_M_END"
log close dlog
assert `rc' == 0
assert e(N_reps) == 20 & e(N_misreps) == 0 & e(N_forced) == 20
log_msgs `"`lgf'"'
assert r(n) == 1
log_para `"`lgf'"' "note: force kept"
assert `"`r(text)'"' == "note: force kept 20 of 20 replications that missed raking targets"
log_count `"`lgf'"' "note: stage 1: 1 raking margin not met"
assert r(n) == 0
foreach t in "bayesboot():" "failed" "force does not keep" "(force kept" {
    log_count `"`lgf'"' "`t'"
    assert r(n) == 0
}
log using `"`lgf'"', text replace name(dlog)
di "SDW_M_BEGIN"
seldecomp
di "SDW_M_END"
log close dlog
log_msgs `"`lgf'"'
assert r(n) == 0
log_count `"`lgf'"' "note: stage 1: 1 raking margin not met"
assert r(n) == 0
log_count `"`lgf'"' "(force kept 20 of 20 replications that missed raking targets)"
assert r(n) == 1
log_count `"`lgf'"' "note: force kept"
assert r(n) == 0
* (g) a thin focal group cell (x3 = 4, 2 obs, in one cluster) that fails replicates of
* bootstrap: and jackknife:
* the thin target level is listed at the bottom of the raking diagnostics, not among
* the header's notes, so the warning after the loop is the one message and names the level;
* at replay the footnote has the count alone; with mincell(2) there is no thin target level,
* and the warning is the same
clear
set seed 20260925
set obs 384
gen int stratum = ceil(_n / 24)
gen byte psu = 1 + mod(ceil(_n / 12) - 1, 2)
gen int psuid = stratum * 10 + psu
gen double wt = 50 + 100 * runiform()
gen byte group = runiform() < .5
gen byte x1 = 1 + (runiform() < .4) + (runiform() < .3 + .3 * group)
gen byte x2 = runiform() < .6 - .2 * group
gen double y = 2 + x1 + x2 + .8 * group + rnormal()
gen byte x3 = x1
replace x3 = 4 if group == 1 & mod(_n, 7) == 0
replace x3 = 4 if group == 0 & stratum == 2 & psu == 1 & mod(_n, 3) == 0
foreach pfx in "bootstrap, reps(40) cluster(psuid) strata(stratum) seed(11) nodots" ///
    "jackknife, cluster(psuid) nodots" {
    log using `"`lgf'"', text replace name(dlog)
    di "SDW_M_BEGIN"
    `pfx': seldecomp y i.x3, by(group) diagnostics
    di "SDW_M_END"
    log close dlog
    local nm = e(N_misreps)
    assert `nm' > 0 & `nm' < .
    local nr = e(N_reps) + `nm'
    local cnt "`nm' of `nr' `e(vce)' replicates failed and `=cond(`nm' == 1, "was", "were")' dropped"
    * the warning, the one message start, names the level; the thin-cell note is in e(note#),
    * and the prefix display summarizes thin target levels below diagnostics, without
    * listing them in the header (no "note: thin target level")
    log_msgs `"`lgf'"'
    assert r(n) == 1
    assert e(N_notes) == 2
    assert `"`e(note1)'"' == "thin target level:" & `"`e(note2)'"' == "    x3 = 4 (2 obs)"
    log_count `"`lgf'"' "note: thin target level"
    assert r(n) == 0
    log_count `"`lgf'"' "Thin target levels (<10 focal obs): 1."
    assert r(n) == 1
    log_para `"`lgf'"' "warning: "
    assert `"`r(text)'"' == "warning: `cnt': in " + cond(`nm' == 1, "it", "them") + ", x3 = 4 had no focal group weight. Collapse that level."
    assert r(lines) == 2
    log_count `"`lgf'"' "in the note above"
    assert r(n) == 0
    log_index `"`lgf'"' "warning: "
    local i1 = r(i)
    log_index `"`lgf'"' "Raking diagnostics"
    local i2 = r(i)
    log_index `"`lgf'"' "Thin target levels ("
    assert `i1' > 0 & `i1' < `i2' & `i2' < r(i)
    * at replay with: the footnote with the count alone; no thin target level listed
    log using `"`lgf'"', text replace name(dlog)
    di "SDW_M_BEGIN"
    seldecomp
    di "SDW_M_END"
    log close dlog
    log_msgs `"`lgf'"'
    assert r(n) == 1
    log_count `"`lgf'"' "x3 = 4"
    assert r(n) == 0
    log_count `"`lgf'"' "Thin target level:"
    assert r(n) == 0
    log_count `"`lgf'"' "(`cnt')"
    assert r(n) == 1
}
* diagnostics typed in the command: the prefix's display applies it (the
* full-sample call posts e(prefix_dispopts)), so requested diagnostics and the thin
* count appear; the replication warning remains the one warning
log using `"`lgf'"', text replace name(dlog)
di "SDW_M_BEGIN"
bootstrap, reps(40) cluster(psuid) strata(stratum) seed(11) nodots: seldecomp y i.x3, by(group) diagnostics
di "SDW_M_END"
log close dlog
assert "`e(prefix_dispopts)'" == "diagnostics"
foreach t in "Raking diagnostics" "Thin target levels (<10 focal obs): 1." {
    log_count `"`lgf'"' "`t'"
    assert r(n) == 1
}
log_msgs `"`lgf'"'
assert r(n) == 1
log_count `"`lgf'"' "warning: "
assert r(n) == 1
* A replay typed afterwards returns to the default factor and thin-count summary.
log using `"`lgf'"', text replace name(dlog)
seldecomp
log close dlog
log_count `"`lgf'"' "Thin target level:"
assert r(n) == 0
log using `"`lgf'"', text replace name(dlog)
di "SDW_M_BEGIN"
bootstrap, reps(40) cluster(psuid) strata(stratum) seed(11) nodots: seldecomp y i.x3, by(group) mincell(2)
di "SDW_M_END"
log close dlog
local nm = e(N_misreps)
assert `nm' > 1 & `nm' < .
log_msgs `"`lgf'"'
assert r(n) == 1
log_para `"`lgf'"' "warning: "
assert `"`r(text)'"' == "warning: `nm' of 40 bootstrap replicates failed and were dropped: in them, x3 = 4 had no focal group weight. Collapse that level."
log using `"`lgf'"', text replace name(dlog)
di "SDW_M_BEGIN"
seldecomp
di "SDW_M_END"
log close dlog
log_msgs `"`lgf'"'
assert r(n) == 1
log_count `"`lgf'"' "(`nm' of 40 bootstrap replicates failed and were dropped)"
assert r(n) == 1
set linesize `ls0'

* ==============================================================================
* 11. display order: optional diagnostics, stage/covariate tables, coefficient
*     table and CI table, followed by a concise raking status
* ==============================================================================
di as txt _n "---- 11. display order"
log_count `"`t1'"' "Raking diagnostics"
assert r(n) == 0
log_index `"`t1'"' "Decomposition"
local i1 = r(i)
log_index `"`t1'"' "            wage | "
local i2 = r(i)
log_index `"`t1'"' "Difference = focal - reference."
local i3 = r(i)
log_index `"`t1'"' "Raking factors"
assert `i1' > 0 & `i1' < `i2' & `i2' < `i3' & `i3' < r(i)
sep_above `"`t1'"' "Decomposition"
* Explicit diagnostics precede covariate detail and decomposition.
log_index `"`t2'"' "Raking diagnostics"
local i1 = r(i)
log_index `"`t2'"' "Raking variables: weighted shares by level"
local i2 = r(i)
log_index `"`t2'"' "Decomposition"
local i3 = r(i)
log_index `"`t2'"' "            wage | "
assert `i1' > 0 & `i1' < `i2' & `i2' < `i3' & `i3' < r(i)
foreach t in "Raking diagnostics" "Raking variables: weighted shares by level" "Decomposition" {
    sep_above `"`t2'"' "`t'"
}
log_count `"`t2'"' "Raking factors"
assert r(n) == 0
log_count `"`t3'"' "Blocks added sequentially (others held at focal group margins):"
assert r(n) == 1
log_count `"`t7'"' "Blocks added sequentially (others unconstrained):"
assert r(n) == 1
foreach t in 3 7 {
    local md = cond(`t' == 3, "sequentially, fixothermargins", "sequentially")
    log_rows `"`t`t''"' "Decomposition"
    assert r(n) == 1 & strrtrim(`"`r(row1)'"') == "Decomposition by raking stage (`md')"
    log_count `"`t`t''"' "Raking diagnostics"
    assert r(n) == 0
    log_index `"`t`t''"' "Decomposition by raking stage"
    local i1 = r(i)
    log_index `"`t`t''"' "            wage | "
    local i2 = r(i)
    log_index `"`t`t''"' "Explained share of the difference:"
    local i3 = r(i)
    log_index `"`t`t''"' "Raking factors"
    assert `i1' > 0 & `i1' < `i2' & `i2' < `i3' & `i3' < r(i)
    sep_above `"`t`t''"' "Decomposition by raking stage (`md')"
}
log_index `"`t7'"' "Raking variables: weighted shares by level"
local i1 = r(i)
log_index `"`t7'"' "Decomposition by raking stage (sequentially)"
assert `i1' > 0 & `i1' < r(i)
sep_above `"`t7'"' "Raking variables: weighted shares by level"
* The replay with citable keeps every table in order.
log_index `"`t5'"' "Decomposition by raking stage"
local i1 = r(i)
log_index `"`t5'"' "            wage | "
local i2 = r(i)
log_index `"`t5'"' "Bayesian bootstrap confidence intervals"
local i3 = r(i)
log_index `"`t5'"' "Raking factors"
assert `i1' > 0 & `i1' < `i2' & `i2' < `i3' & `i3' < r(i)
* Individual and reverse modes retain their full stage summaries and explanations.
sysuse nlsw88, clear
foreach f in "" fixothermargins {
    log using `"`lgf'"', text replace name(dlog)
    seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) individually `f'
    log close dlog
    assert "`e(stagetype)'" == "individually" & "`e(fixothermargins)'" == "`f'"
    local md = "individually" + cond("`f'" != "", ", fixothermargins", "")
    local how = cond("`f'" != "", "others held at focal group margins", "others unconstrained")
    log_para `"`lgf'"' "Blocks raked individually ("
    assert `"`r(text)'"' == "Blocks raked individually (`how'): block 1: demo (race collgrad) block 2: geo (south smsa)"
    assert r(lines) == 3
    log_rows `"`lgf'"' "Contribution of each block, raked alone ("
    assert r(n) == 1 & strrtrim(`"`r(row1)'"') == "Contribution of each block, raked alone (`md')"
    local how = cond("`f'" != "", "are held at their focal group distributions", "are unconstrained")
    log_para `"`lgf'"' "%: percent of the difference. Other blocks"
    assert `"`r(text)'"' == "%: percent of the difference. Other blocks `how'. Remainder: explained part of all blocks less the sum of the blocks' parts."
    log_count `"`lgf'"' "Raking diagnostics"
    assert r(n) == 0
    sep_above `"`lgf'"' "Contribution of each block, raked alone ("
    log_index `"`lgf'"' "Contribution of each block, raked alone"
    local i1 = r(i)
    log_index `"`lgf'"' "            wage | "
    assert `i1' > 0 & `i1' < r(i)
}
log using `"`lgf'"', text replace name(dlog)
seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially reverse
log close dlog
log_count `"`lgf'"' "Decomposition"
assert r(n) == 2
log_index `"`lgf'"' "Decomposition by raking stage (sequentially)"
local i1 = r(i)
log_index `"`lgf'"' "Decomposition by raking stage (sequentially), reverse order"
local i2 = r(i)
log_index `"`lgf'"' "            wage | "
assert `i1' > 0 & `i1' < `i2' & `i2' < r(i)
log_after `"`lgf'"' "Reverse order: " 0
assert strrtrim(`"`r(line)'"') == "Reverse order: blocks added from last to first."
* Thin-cell details are optional; their stored notes survive either display.
foreach o in "" diagnostics {
    log using `"`lgf'"', text replace name(dlog)
    seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) ///
        sequentially fixothermargins `o'
    log close dlog
    assert e(N_notes) == 2 & `"`e(note1)'"' == "thin target level, stages 1-3:"
    log_count `"`lgf'"' "note: thin target level"
    assert r(n) == 0
    log_count `"`lgf'"' "race = 3 (Other) (8 obs)"
    assert r(n) == 0
    log_count `"`lgf'"' "Raking factors"
    assert r(n) == ("`o'" == "")
    sep_above `"`lgf'"' "Decomposition by raking stage"
}

capture erase `"`lgf'"'

* ==============================================================================
* 13. varlabels: the variable labels in place of the variable names, in the
*     header, the notes, the covariate detail, the stage table, the coefficient tables and
*     the graphs; a named block keeps its name; without varlabels the names
* ==============================================================================
di as txt _n "---- 13. varlabels"
sysuse nlsw88, clear
label variable race "Race"
label variable collgrad "College graduate"
label variable south "Lives in the south"
label variable union "Union worker"
label variable smsa "Lives in SMSA"
label variable married "Married"
label variable wage "Hourly wage"
log using `"`lgf'"', text replace name(dlog)
seldecomp wage (demo: race collgrad) south union#smsa, by(married) sequentially fixothermargins ///
    detail covariatedetail varlabels diagnostics bayesboot(reps(5) rseed(1) nodots)
log close dlog
* Variable labels are retained in headers, stages, covariate detail, and main rows.
foreach t in "Focal Group (raked): Married = 0 (Single)" "stage 1: + demo (Race, College graduate)" ///
    "stage 2: + Lives in the south" "stage 3: + Union worker # Lives in SMSA" ///
    "    2 Lives in th~ " "College graduate |" "     Hourly wage |" "Lives in the so~ |" ///
    "Union worker # Lives in SMSA" "    1 demo " {
    log_count `"`lgf'"' "`t'"
    assert r(n) >= 1
}
foreach t in "    2 south " "married = 0" "race = 3" {
    log_count `"`lgf'"' "`t'"
    assert r(n) == 0
}
* the replay without varlabels shows the names
log using `"`lgf'"', text replace name(dlog)
seldecomp, diagnostics
log close dlog
foreach t in "Focal Group (raked): married = 0 (Single)" "stage 2: + south" ///
    "           south |" "            wage |" {
    log_count `"`lgf'"' "`t'"
    assert r(n) >= 1
}
log_count `"`lgf'"' "Lives in th~"
assert r(n) == 0
* the graphs: a block named after its terms shows their labels; a named block its name
seldecomp, makeplot(marginal, twoway(name(gvl, replace))) varlabels
assert `"`.gvl.yaxis1.major.dlg_ticklabel[3]'"' == "demo"
assert `"`.gvl.yaxis1.major.dlg_ticklabel[4]'"' == "Lives in the south"
assert `"`.gvl.yaxis1.major.dlg_ticklabel[5]'"' == "Union worker # Lives in SMSA"
seldecomp, makeplot(marginal, twoway(name(gvl, replace)))
assert `"`.gvl.yaxis1.major.dlg_ticklabel[4]'"' == "south"
graph drop gvl

di as res _n "PASSED: test_display.do"

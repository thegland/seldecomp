* test_plot.do -- sequentially / individually + makeplot: the plotted data
*
* Every graph is checked through the data it draws, not only through the return
* code: check_plot (below) reads the sersets of the drawn graph and asserts that the
* explained and unexplained bars of makeplot(stages) start at 0 and end at the stage values in
* e(b) and e(bystage), side by side; that the dashed line is the total difference; and that the
* share labels and the confidence-interval caps are right; check_contrib does the same for
* makeplot(marginal) (one bar per stage's marginal explained part, with the difference and the
* explained and unexplained parts, horizontal or, with swapaxis, vertical; spikes with ci). The
* sign cases covered are explained and unexplained of the same sign, of opposite signs (help
* example 2 on nlsw88), explained larger than the difference, and a difference of exactly 0; and
* the graph drawn at replay after a replication prefix (point estimates, not the last
* replicate), also when only e(b) survives. Also: swapaxis in every plot type, sequentially alone
* (the other blocks free) and sequentially fixothermargins (later blocks held at the Focal Group's
* shares), individually with and without fixothermargins (the marginal plot alone, with the
* remainder; the other plot types refused), one stage, twelve long block names, the titles (the
* marginal plot's shortened to the width of the graph), the colours under three schemes, and the
* data in memory left unchanged. Also the waterfall (check_waterfall: every bar from the running
* total of the steps in e(b), the spikes, labels and connectors), the waterfall in levels,
* makeplot(waterfalllevels) (check_levels: the means and the running levels from e(b), the base
* of the means' bars, spikes, labels, connectors), the stages of the reverse order in every plot
* type (makeplot(..., reverse)), the default graph name, seldecomp, and the default plot type
* (the waterfall with sequentially, also for a chain, the whole chain and a step of it; the
* marginal plot with individually). The makeplot() syntax: makeplot([plottype] [, options]), the
* plot types stages, marginal, waterfall and waterfalllevels, the options swapaxis, ci, reverse,
* chainstep(), outcome(), format(), blocklabels() and twoway(), every twoway option inside
* twoway(); confidence intervals are drawn only with ci, on every bar. CI plots identify their
* level and interval type in a note (all interval families and changed levels are also checked by
* test_plot_intervals.do). The plotspec refused before the estimation with e() left as it was,
* and at replay (198, with the message); makeplot without sequentially or individually,
* fixothermargins alone, sequentially with individually and words that are no options refused (198,
* with the message where there is one). Also: no interval unless ci asks for one; misspelt plot
* types and options, outcome() and the other options refused before the estimation; the sign of
* the shares in levels; format(); the titles, a label off the total line, one angle for all
* names; blocklabels() in every plot type, with $, `, braces and quotes. The difference is
* always Focal Group - Reference Group (e(diffsign)).
* Graphs are drawn to memory, which also works in console (batch) Stata.
*
* Run from the tests directory:   . do test_plot.do

version 16
clear
set more off
set varabbrev off
discard             // drop any seldecomp already in memory; the repository copy is what gets tested

local tol 1e-6

* ---- locate seldecomp.ado -----------------------------------------------------
* the repository copy of seldecomp.ado always goes first on the ado-path (see the -discard-
* at the top of this file)
if `"$SELDECOMP_ROOT"' != "" adopath ++ `"$SELDECOMP_ROOT"'
else adopath ++ `"`c(pwd)'/.."'
capture which seldecomp
if _rc {
    di as err "seldecomp.ado not found: run this do-file from the tests directory " ///
        "or set the global SELDECOMP_ROOT to the repository root"
    exit 601
}
which seldecomp        // record in the log which copy is being tested

* ---- plot_data: load one plot's data ------------------------------------------
* plot_data graphname #
* Replaces the data in memory with the serset of plot # of the graph, its variables
* renamed v1, v2, ... in serset order. The plot's class object (.graphname.plotregion1.plot#)
* holds the serset id and the positions of its y, y2, x, marker-label and label-position
* variables (.yvar .y2var .xvar .lvar .pvar); they are returned as v_y v_y2 v_x v_l v_p
* (empty when the plot has none), with the plot type (.type.snm) in v_type, the bar
* width in v_barsize, the marker symbol in v_symbol and the size of the marker labels in
* v_lsize.
capture program drop plot_data
program define plot_data
    args G k
    local P .`G'.plotregion1.plot`k'
    local id ``P'.serset.id'
    confirm integer number `id'
    qui serset set `id'
    qui serset use, clear
    // -serset use- names the variables like temporary variables (dropped when this program
    // ends, and not matched by rename's *): rename them by position
    mata: for (i = 1; i <= st_nvar(); i++) st_varrename(i, "v" + strofreal(i))
    foreach m in y y2 x x2 l p {
        local i ``P'.`m'var'
        if "`i'" == "" | "`i'" == "." c_local v_`m'
        else c_local v_`m' v`i'
    }
    c_local v_type ``P'.type.snm'
    c_local v_barsize ``P'.bar_size'
    c_local v_symbol ``P'.style.marker.symbol.snm'
    c_local v_lsize ``P'.style.label.textstyle.size.snm'
end

* ---- crit: the critical value of the intervals drawn ---------------------------
* crit lname: z for e(level), or t with e(df_r) when it is set, as in the coefficient table: by
* a prefix (the residual (design) degrees of freedom) or by bayesboot() when t is the family of
* its tests (e(df_bs); the default with bayesboot(cluster()))
capture program drop crit
program define crit
    args lname
    local df = e(df_r)
    if `df' < . c_local `lname' = invttail(`df', (100 - e(level)) / 200)
    else c_local `lname' = invnormal(1 - (100 - e(level)) / 200)
end

* ---- the texts that name the groups and the outcome ----------------------------------------
* The help's rules, written out here. tp_esc(s): a text as the graph stores it to show it as
* typed: {, }, $ and ` as {c -(}, {c )-}, {c S|} and {c 96}, and every double quote
* as {c 34}. tp_fit(s, n): s in at most n characters, a cut text ending with
* "..." after its last whole word when that keeps more than half of the n - 3 characters, else
* within the word, without spaces, commas, semicolons or colons before the "...". tp_wrap(s, n,
* tail): a row of a group in the horizontal levels view: s
* (trimmed when longer than n) and tail on one line when s has at most n characters, else two
* lines split at a space with both at most n, after a comma when possible, else the split with
* the shortest longer line; when none fits, the words that fit on the first line (a first word
* longer than n cut by tp_fit) and the rest cut by tp_fit, tail after the second line; stored as
* `"a"' `"b"'; a text without spaces on one line, cut
capture mata: mata drop tp_esc() tp_fit() tp_wrap()
mata:
string scalar tp_esc(string scalar s)
{
    real scalar i
    string scalar c, out
    out = ""
    for (i = 1; i <= ustrlen(s); i++) {
        c = usubstr(s, i, 1)
        if (c == "{") out = out + "{c -(}"
        else if (c == "}") out = out + "{c )-}"
        else if (c == char(36)) out = out + "{c S|}"
        else if (c == char(96)) out = out + "{c 96}"
        else if (c == char(34)) out = out + "{c 34}"
        else out = out + c
    }
    return(out)
}
string scalar tp_fit(string scalar s, real scalar n)
{
    real scalar i, p
    if (ustrlen(s) <= n) return(s)
    p = 0
    for (i = 1; i <= n - 2; i++) {
        if (usubstr(s, i, 1) == " ") p = i
    }
    if (p == n - 2 | p <= (n - 3) / 2) return(ustrregexra(usubstr(s, 1, n - 3), "[ ,;:]+$", "") + "...")
    return(ustrregexra(usubstr(s, 1, p - 1), "[ ,;:]+$", "") + "...")
}
string scalar tp_wrap(string scalar s, real scalar n, string scalar tail)
{
    real scalar i, k, cm, best, bk, bc
    string scalar a, b
    if (ustrlen(s) <= n) return(tp_esc(s + tail))
    s = ustrtrim(s)
    if (ustrlen(s) <= n) return(tp_esc(s + tail))
    if (!ustrpos(s, " ")) return(tp_esc(tp_fit(s, n) + tail))
    best = 0
    bk = .
    bc = 0
    for (i = 2; i < ustrlen(s); i++) {
        if (usubstr(s, i, 1) != " ") continue
        a = ustrrtrim(usubstr(s, 1, i - 1))
        b = ustrltrim(usubstr(s, i + 1, .))
        k = max((ustrlen(a), ustrlen(b)))
        cm = usubstr(a, ustrlen(a), 1) == ","
        if (k > n) continue
        if (cm > bc | (cm == bc & k < bk)) {
            best = i
            bk = k
            bc = cm
        }
    }
    if (!best) {
        for (i = min((n + 1, ustrlen(s))); i >= 2; i--) {
            if (usubstr(s, i, 1) == " " & !best) best = i
        }
    }
    if (best) {
        a = ustrrtrim(usubstr(s, 1, best - 1))
        b = ustrltrim(usubstr(s, best + 1, .))
    }
    else {
        a = tp_fit(usubstr(s, 1, ustrpos(s, " ") - 1), n)
        b = ustrltrim(usubstr(s, ustrpos(s, " ") + 1, .))
    }
    a = char(96) + char(34) + tp_esc(a) + char(34) + char(39)
    b = char(96) + char(34) + tp_esc(tp_fit(b, n) + tail) + char(34) + char(39)
    return(a + " " + b)
}
end

* ---- ptexts: the texts of e() that the graphs show --------------------------------------------
* ptexts [outcome]: c_local gr1 gr2, the groups as typed (the value label of e(group_#_label),
* "value (label)", else by = value, else the string value), sr1 sr2 (an unlabelled number in
* parentheses, as the subtitle and title name it), and yt0, the labels of the value axis (the
* variable label of the outcome drawn, else its name; with ratio() " / " and that of the
* denominator)
capture program drop ptexts
program define ptexts
    args outcome
    forvalues g = 1/2 {
        mata: st_local("gv", st_global("e(" + (`g' == 1 ? "focal" : "ref") + "_group)"))
        mata: st_local("gl", st_global("e(" + (`g' == 1 ? "focal" : "ref") + "_group_label)"))
        mata: st_local("gr", st_local("gl") != st_local("gv") ? usubstr(st_local("gl"), ///
            ustrlen(st_local("gv")) + 3, ustrlen(st_local("gl")) - ustrlen(st_local("gv")) - 3) : ///
            (strtoreal(st_local("gv")) < . ? st_global("e(by)") + " = " + st_local("gv") : st_local("gv")))
        mata: st_local("sr", st_local("gl") == st_local("gv") & strtoreal(st_local("gv")) < . ? ///
            "(" + st_local("gr") + ")" : st_local("gr"))
        c_local gr`g' `"`macval(gr)'"'
        c_local sr`g' `"`macval(sr)'"'
    }
    local dv `e(depvar)'
    if "`e(outcomes)'" != "" local dv : word 1 of `e(outcomes)'
    if "`outcome'" != "" local dv `outcome'
    local yt0
    foreach v in `dv' `e(ratio)' {
        local vl `v'
        capture confirm variable `v', exact
        if !_rc {
            mata: st_local("vl", st_varlabel("`v'"))
        }
        if `"`macval(vl)'"' == "" local vl `v'
        * (brace blocks: a one-line if expands the label again)
        if `"`macval(yt0)'"' == "" {
            local yt0 `"`macval(vl)'"'
        }
        else {
            local yt0 `"`macval(yt0)' / `macval(vl)'"'
        }
    }
    c_local yt0 `"`macval(yt0)'"'
end

* ---- text_lines: a text element of the graph against its expected lines --------------------------
* text_lines "object" "line" ["line"]: the text object (.G.subtitle, .G.xaxis1.title) has exactly
* the lines given
capture program drop text_lines
program define text_lines
    args o l1 l2
    local n = cond(`"`macval(l2)'"' != "", 2, 1)
    if `.`o'.text.arrnels' != `n' {
        di as err `"check `o': `.`o'.text.arrnels' lines, expected `n'"'
        exit 9
    }
    forvalues k = 1/`n' {
        if `"`.`o'.text[`k']'"' != `"`macval(l`k')'"' {
            di as err `"check `o': line `k' is |`.`o'.text[`k']'|, expected |`macval(l`k')'|"'
            exit 9
        }
    }
end

* ---- no_text: a text element of the graph without any line ----------------------------------------
* no_text "object": the text object (.G.note, .G.xaxis1.title) has no line (twoway leaves a note
* and an axis title given as "" without lines)
capture program drop no_text
program define no_text
    args o
    local n `.`o'.text.arrnels'
    if !inlist("`n'", "", "0") {
        di as err `"check `o': `n' lines, |`.`o'.text[1]'| ..., none expected"'
        exit 9
    }
end

* ---- check_texts: the subtitle, the axis titles and the note ------------------------------------
* check_texts G, view(stages|marginal|waterfall|levels) [horizontal outcome(name) reverse ci]: the
* subtitle (stages, marginal, waterfall): the groups in the order of the difference, Focal Group -
* Reference Group, with ", reverse order" after them for the reverse order's stages: "A - B" on
* one line up to 72 characters (the suffix counted), else "A -" and "B", each tp_fit to 70 (the
* suffix after B); levels: no subtitle (the title names the groups), "Reverse order" for the
* reverse order's stages; the value axis: "<yt0>: <what>" (what "difference in means" or "ratios",
* in levels "mean" or "ratio") on one line up to 40 characters, else "<yt0>:" (tp_fit to 40) and
* "<what>"; the category axis without a title; a level/type note only when confidence intervals
* are drawn
capture program drop check_texts
program define check_texts
    syntax name(name=G), view(string) [HORizontal OUTcome(string) REVerse ci]
    assert inlist("`view'", "stages", "marginal", "waterfall", "levels")
    ptexts `outcome'
    local ym = cond("`e(ratio)'" != "", "ratio", "mean")
    local sfx = cond("`reverse'" != "", ", reverse order", "")
    if "`view'" != "levels" {
        * the difference is always Focal Group - Reference Group
        assert "`e(diffsign)'" == "focal_group - ref_group"
        mata: st_local("n", strofreal(ustrlen(st_local("sr1") + " - " + st_local("sr2") + st_local("sfx"))))
        if `n' <= 72 {
            mata: st_local("s1", tp_esc(st_local("sr1") + " - " + st_local("sr2") + st_local("sfx")))
            text_lines `G'.subtitle `"`macval(s1)'"'
        }
        else {
            mata: st_local("s1", tp_esc(tp_fit(st_local("sr1"), 70) + " -"))
            mata: st_local("s2", tp_esc(tp_fit(st_local("sr2"), 70) + st_local("sfx")))
            text_lines `G'.subtitle `"`macval(s1)'"' `"`macval(s2)'"'
        }
        local what "difference in `ym's"
    }
    else {
        if "`reverse'" != "" text_lines `G'.subtitle "Reverse order"
        else no_text `G'.subtitle
        local what `ym'
    }
    local vax = cond("`horizontal'" != "", "xaxis1", "yaxis1")
    local cax = cond("`horizontal'" != "", "yaxis1", "xaxis1")
    mata: st_local("n", strofreal(ustrlen(st_local("yt0") + ": `what'")))
    if `n' <= 40 {
        mata: st_local("v1", tp_esc(st_local("yt0") + ": `what'"))
        text_lines `G'.`vax'.title `"`macval(v1)'"'
    }
    else {
        mata: st_local("v1", tp_esc(tp_fit(st_local("yt0"), 40) + ":"))
        text_lines `G'.`vax'.title `"`macval(v1)'"' "`what'"
    }
    no_text `G'.`cax'.title
    if "`ci'" == "" no_text `G'.note
    else {
        local ct : word 1 of `e(citypes)'
        if "`ct'" == "" local ct = cond((e(df_r) > 0 & e(df_r) < .) | ///
            (e(N_clust) > 0 & e(N_clust) < .), "t", "normal")
        local ct = cond("`ct'" == "bc", "bias-corrected", cond("`ct'" == "bca", "BCa", "`ct'"))
        text_lines `G'.note "`e(level)'% `ct' intervals"
    }
end

* ---- check_plot: the stages plot against e() ----------------------------------
* check_plot graphname [, ci horizontal reverse tolerance(#) subtitle() vtitle() format() blabels()
*                         notitle]
* Legend key j of the graph (.graphname.legend.labels[j]) names plot .legend.map[j].
* Asserted, for every stage s:
*   - the plot keyed "Expl." is an rbar from 0 to cum_expl of stage s and the plot keyed
*     "Unexpl." an rbar from 0 to cum_unexpl, both in the slot of stage s and side by side
*     (their centres at least a bar width apart);
*   - the plot keyed "Total diff (#)", # the difference, is flat at e(difference) and spans all
*     bars;
*   - the share labels sit at the bar ends (at the outer end of the CI cap when there
*     is one), read 100 * part / difference (to the rounding of the label, never "-0%";
*     blank when the difference is 0) and are placed above positive and below negative
*     bars (right and left with -horizontal-);
*   - with -ci-: the two rcaps run from part - z * se to part + z * se (se from e(V), z
*     for e(level), t with e(df_r) when it is set); without -ci-: no rcap, whether or not
*     e(V) is posted (the intervals are drawn only with makeplot(..., ci));
*   - a label that the dashed line would strike is moved to the line: a label runs
*     from its anchor away from its bar, and when the line lies within 5% of the value axis
*     beyond the anchor (in the horizontal view, (.02 + .018 * the label's length) of it),
*     the anchor is the line; the axis is the span the ado computes (0, the difference and
*     the anchors, with a margin of 8%, 10% horizontal, on the sides that have bars);
*   - the rows are named by the blocks, "+ " before all but the first (by -blabels()-, the
*     texts given to blocklabels(), as the graph stores them, in the order drawn, when given),
*     angled as the ado does it (beyond four stages or a name longer than 10 characters, in the
*     vertical view; a {c ...} code counts as one character);
*   - the legend has exactly the three keys, the total's value in -format()- (%9.3g by
*     default); the title "Explained and unexplained difference by raking stage" (not with
*     -notitle-, for a title() given in twoway()); the subtitle is -subtitle()- and the title of
*     the value axis is -vtitle()- when these are given; the texts as check_texts has them (no
*     title on the category axis, no note).
* With -horizontal- the stages run down the y axis: the bars' roles are unchanged
* (y = 0, y2 = value, x = slot), while the line and the labels are drawn as x = value.
* The expected values are e(b)'s cum_expl:<block> and cum_unexpl:<block> (the point
* estimates, also after a replication prefix) and, when it is posted, e(bystage)'s. With
* -reverse- (makeplot(..., reverse), results of sequentially reverse) the stages of the reverse
* order: the blocks from the last to the first, from rev_cum_expl:<block> and
* rev_cum_unexpl:<block> (e(bystage) holds the forward order's), and ", reverse order" in the
* subtitle.
capture program drop check_plot
program define check_plot
    syntax name(name=G) [, ci HORizontal TOLerance(real 1e-10) SUBtitle(string) VTItle(string) ///
        FORmat(string) BLAbels(string asis) REVerse NOTItle]
    local tol = `tolerance'
    if "`format'" == "" local format %9.3g
    local ns = e(N_stages)
    local blocks `e(blocks)'
    local rp
    if "`reverse'" != "" {
        assert "`e(reverse)'" != ""
        local rp rev_
        local blocks
        forvalues s = `ns'(-1)1 {
            local blocks `blocks' `: word `s' of `e(blocks)''
        }
    }
    local D = e(difference)
    crit z
    tempname b S V X XE XU Dsc
    scalar `Dsc' = e(difference)
    matrix `b' = e(b)
    assert reldif(scalar(`Dsc'), `b'[1, colnumb(`b', "overall:difference")]) < `tol'
    local hasS 0
    capture confirm matrix e(bystage)
    if !_rc & "`reverse'" == "" {
        matrix `S' = e(bystage)
        local hasS 1
    }
    if "`ci'" != "" matrix `V' = e(V)
    * X: 1 explained, 2 unexplained (e(b)); 3, 4 their standard errors; 5, 6 e(bystage)
    matrix `X' = J(`ns', 6, .)
    forvalues s = 1/`ns' {
        local bn : word `s' of `blocks'
        local ce = colnumb(`b', "`rp'cum_expl:`bn'")
        local cu = colnumb(`b', "`rp'cum_unexpl:`bn'")
        assert `ce' < . & `cu' < .
        matrix `X'[`s', 1] = `b'[1, `ce']
        matrix `X'[`s', 2] = `b'[1, `cu']
        if "`ci'" != "" {
            assert `V'[`ce', `ce'] > 0 & `V'[`ce', `ce'] < . & `V'[`cu', `cu'] > 0 & `V'[`cu', `cu'] < .
            matrix `X'[`s', 3] = sqrt(`V'[`ce', `ce'])
            matrix `X'[`s', 4] = sqrt(`V'[`cu', `cu'])
        }
        if `hasS' {
            matrix `X'[`s', 5] = `S'[`s', colnumb(`S', "cum_expl")]
            matrix `X'[`s', 6] = `S'[`s', colnumb(`S', "cum_unexpl")]
        }
    }
    local lpos = cond("`horizontal'" != "", 3, 12)
    local lneg = cond("`horizontal'" != "", 9, 6)

    * titles
    if `"`subtitle'"' != "" assert `"`.`G'.subtitle.text[1]'"' == `"`subtitle'"'
    local vax = cond("`horizontal'" != "", "xaxis1", "yaxis1")
    if `"`vtitle'"' != "" assert `"`.`G'.`vax'.title.text[1]'"' == `"`vtitle'"'
    check_texts `G', view(stages) `horizontal' `reverse' `ci'
    if "`notitle'" == "" text_lines `G'.title "Explained and unexplained difference by raking stage"
    * the rows: the block names or the texts of blocklabels(), "+ " before all but the first;
    * angled (30 degrees) beyond four stages or for a text of more than 10 characters, in the
    * vertical view
    local cax = cond("`horizontal'" != "", "yaxis1", "xaxis1")
    blabel_rows `"`blocks'"' `"`blabels'"'
    local long 0
    forvalues s = 1/`ns' {
        local t `"`.`G'.`cax'.major.dlg_ticklabel[`s']'"'
        assert `"`t'"' == `"`=cond(`s' > 1, "+ ", "")'`rl_`s''"'
        if `rn_`s'' > 10 local long 1
    }
    if "`horizontal'" == "" {
        assert `.`G'.xaxis1.style.majorstyle.tickangle.val' == cond(`ns' > 4 | `long', 30, 0)
    }

    * legend keys -> plots
    local Dtxt = strtrim(string(`D', "`format'"))
    if regexm("`Dtxt'", "^-[0.,]+$") local Dtxt = substr("`Dtxt'", 2, .)
    local nk = `.`G'.legend.map.arrnels'
    assert `nk' == 3
    forvalues j = 1/`nk' {
        local t `.`G'.legend.labels[`j']'
        local p `.`G'.legend.map[`j']'
        if `"`t'"' == "Expl." local pE `p'
        else if `"`t'"' == "Unexpl." local pU `p'
        else if `"`t'"' == "Total diff (`Dtxt')" local pT `p'
        else {
            di as err `"check_plot `G': unexpected legend key `t'"'
            exit 9
        }
    }
    if "`pE'" == "" | "`pU'" == "" | "`pT'" == "" {
        di as err "check_plot `G': legend keys Expl., Unexpl., Total diff (`Dtxt') expected"
        exit 9
    }
    local np 0
    forvalues k = 1/30 {
        if "`.`G'.plotregion1.plot`k'.classname'" == "" continue, break
        local np `k'
    }

    preserve
    * the bars: from 0 to the stage value, one pair per stage slot, side by side
    matrix `XE' = J(`ns', 1, .)
    matrix `XU' = J(`ns', 1, .)
    foreach part in E U {
        local col = cond("`part'" == "E", 1, 2)
        plot_data `G' `p`part''
        assert "`v_type'" == "rbar"
        assert _N == `ns'
        sort `v_x'
        assert `v_y' == 0
        forvalues s = 1/`ns' {
            assert reldif(`v_y2'[`s'], `X'[`s', `col']) < `tol'
            if `hasS' assert reldif(`v_y2'[`s'], `X'[`s', `col' + 4]) < `tol'
            assert abs(`v_x'[`s'] - `s') < .5
            matrix `X`part''[`s', 1] = `v_x'[`s']
        }
        local w`part' = `v_barsize'
    }
    local lo = .
    local hi = -1
    forvalues s = 1/`ns' {
        assert abs(`XU'[`s', 1] - `XE'[`s', 1]) >= (`wE' + `wU') / 2
        local lo = min(`lo', `XE'[`s', 1] - `wE' / 2, `XU'[`s', 1] - `wU' / 2)
        local hi = max(`hi', `XE'[`s', 1] + `wE' / 2, `XU'[`s', 1] + `wU' / 2)
    }

    * the total difference: one flat dashed line at e(difference) across all bars
    assert "`.`G'.plotregion1.plot`pT'.style.line.pattern.snm'" == "dash"
    plot_data `G' `pT'
    assert "`v_type'" == "line"
    local lval = cond("`horizontal'" != "", "`v_x'", "`v_y'")
    local lslot = cond("`horizontal'" != "", "`v_y'", "`v_x'")
    summarize `lval', meanonly
    assert r(min) == r(max)
    assert reldif(r(min), scalar(`Dsc')) < `tol'
    summarize `lslot', meanonly
    assert r(min) <= `lo' & r(max) >= `hi'

    * the anchors of the labels: the bar end, or the outer end of the CI cap; the value axis the
    * ado spans (a label within reach of the dashed line is moved to it)
    local alo = min(0, `D')
    local ahi = max(0, `D')
    local nneg 0
    local npos 0
    forvalues s = 1/`ns' {
        forvalues col = 1/2 {
            local a`s'_`col' = `X'[`s', `col'] + cond("`ci'" != "", ///
                cond(`X'[`s', `col'] >= 0, 1, -1) * `z' * `X'[`s', `col' + 2], 0)
            local alo = min(`alo', `a`s'_`col'')
            local ahi = max(`ahi', `a`s'_`col'')
            if `X'[`s', `col'] < 0 local ++nneg
            else local ++npos
        }
    }
    local pad = (`ahi' - `alo') * cond("`horizontal'" != "", .1, .08) * (`D' != 0)
    if `nneg' local alo = `alo' - `pad'
    if `npos' local ahi = `ahi' + `pad'
    local arg = `ahi' - `alo'
    local nmoved 0

    * CI caps and share labels, matched to their bar by the slot position
    local ncap 0
    local nlab 0
    forvalues k = 1/`np' {
        plot_data `G' `k'
        local islab = "`v_type'" == "scatter" & "`v_l'" != ""
        if "`v_type'" != "rcap" & !`islab' continue
        assert _N == `ns'
        * rcaps keep the bar roles; labels are drawn as (y, x) = (slot, value) when horizontal
        local slot `v_x'
        local val `v_y'
        if `islab' & "`horizontal'" != "" {
            local slot `v_y'
            local val `v_x'
        }
        sort `slot'
        local part = cond(reldif(`slot'[1], `XE'[1, 1]) < 1e-12, "E", "U")
        local col = cond("`part'" == "E", 1, 2)
        forvalues s = 1/`ns' {
            assert reldif(`slot'[`s'], `X`part''[`s', 1]) < 1e-12
        }
        if "`v_type'" == "rcap" {
            local ++ncap
            forvalues s = 1/`ns' {
                assert reldif(`v_y'[`s'], `X'[`s', `col'] - `z' * `X'[`s', `col' + 2]) < `tol'
                assert reldif(`v_y2'[`s'], `X'[`s', `col'] + `z' * `X'[`s', `col' + 2]) < `tol'
            }
        }
        else {
            local ++nlab
            forvalues s = 1/`ns' {
                * anchored at the bar end, or at the outer end of the CI cap when there is one;
                * at the line when the line lies within the label's reach beyond that end
                local anchor = `a`s'_`col''
                local sd = cond(`X'[`s', `col'] >= 0, 1, -1)
                local ext = cond("`horizontal'" != "", .02 + .018 * strlen(`v_l'[`s']), .05)
                if `D' != 0 & `sd' * (`D' - `anchor') > 0 & `sd' * (`D' - `anchor') < `ext' * `arg' {
                    local anchor = `D'
                    local ++nmoved
                }
                assert reldif(`val'[`s'], `anchor') < `tol'
                if `D' != 0 {
                    assert substr(`v_l'[`s'], -1, 1) == "%"
                    assert `v_l'[`s'] != "-0%"
                    assert abs(real(subinstr(`v_l'[`s'], "%", "", 1)) - 100 * `X'[`s', `col'] / `D') <= .5
                }
                else assert `v_l'[`s'] == ""
                assert `v_p'[`s'] == cond(`X'[`s', `col'] >= 0, `lpos', `lneg')
            }
        }
    }
    assert `ncap' == cond("`ci'" != "", 2, 0)
    assert `nlab' == 2
    restore
    c_local check_moved `nmoved'
    di as txt "check_plot `G': `ns' stages, difference `Dtxt'" cond("`ci'" != "", ", CI caps", "") ///
        cond("`horizontal'" != "", ", horizontal", "") cond("`reverse'" != "", ", reverse order", "") ///
        cond(`nmoved', ", `nmoved' label(s) at the line", "") " -- ok"
end

* ---- blabel_rows: the expected names of the block rows -----------------------------------
* blabel_rows "blocks" "texts": c_local rl_# the name of block # as the graph stores it (the
* block name, "Remainder" for _remainder, or the #th of -texts-, a list of the texts given to
* blocklabels() in the form the graph stores them, {c S|} for $ and so on) and rn_# its length
* ({c ...} one character)
capture program drop blabel_rows
program define blabel_rows
    args blocks texts
    local nb : word count `blocks'
    forvalues k = 1/`nb' {
        local rl : word `k' of `blocks'
        if "`rl'" == "_remainder" local rl Remainder
        else if `"`macval(texts)'"' != "" {
            mata: st_local("rl", tokens(st_local("texts"))[`k'])
        }
        mata: st_local("rn", strofreal(ustrlen(ustrregexra(st_local("rl"), "\{c [^}]*\}", "x"))))
        c_local rl_`k' `"`macval(rl)'"'
        c_local rn_`k' `rn'
    }
end

* ---- fmt_lab: a value as makeplot(format()) labels it ------------------------------------
* fmt_lab lname value format [plus]: the value in the format, trimmed; a value that reads 0
* without its sign; with -plus-, "+" before a positive value that does not read 0
capture program drop fmt_lab
program define fmt_lab
    args lname v fmt plus
    local s = strtrim(string(`v', "`fmt'"))
    if regexm("`s'", "^-[0.,]+$") local s = substr("`s'", 2, .)
    if "`plus'" != "" & `v' > 0 & !regexm("`s'", "^[0.,]+$") local s +`s'
    c_local `lname' `s'
end

* ---- fit_title: the title of the marginal plot -------------------------------------------
* fit_title lname xsize ysize: the first of "Marginal explanatory contribution at each raking
* stage", "... by raking stage", "... at each stage" and "... by stage" that fits on one line,
* else the last (Plot_fit in the ado): the room is 56 characters at the default size of the s1
* and s2 schemes, 5.5 by 4 inches (the st schemes' 7.5 by 4.5), and grows as xsize / min(xsize,
* ysize); xsize and ysize those of twoway(xsize() ysize()), -1 for the scheme's
capture program drop fit_title
program define fit_title
    args lname xs ys
    local st = substr(c(scheme), 1, 2) == "st"
    if `xs' <= 0 local xs = cond(`st', 7.5, 5.5)
    if `ys' <= 0 local ys = cond(`st', 4.5, 4)
    local w = 56 * (`xs' / min(`xs', `ys')) / (5.5 / 4)
    foreach t in "at each raking stage" "by raking stage" "at each stage" "by stage" {
        local ttl "Marginal explanatory contribution `t'"
        if strlen("`ttl'") <= `w' continue, break
    }
    c_local `lname' `ttl'
end

* ---- check_contrib: the marginal plot against e() ------------------------
* check_contrib graphname [, ci vertical reverse tolerance(#) subtitle() vtitle() format()
*                            blabels() xsize(#) ysize(#) notitle]
* The rows are read from the tick labels of the category axis (.yaxis1.major.dlg_ticklabel[j] at
* .dlg_tickpos[j]; the x axis with -vertical-, makeplot(marginal, swapaxis)): "Total diff",
* "Expl.", one row per block (e(blocks), in order; "Remainder" after them with individually)
* and "Unexpl.", at the positions 1, 2.2, 3.2, ..., and the last row's number + .4 (drawn top to
* bottom, left to right with -vertical-). Asserted:
*   - every row has exactly one bar, an rbar from 0 to its value in e(b): the
*     difference, the explained part, marg_expl:<block> and the unexplained part; the
*     difference a dashed outline (plot 1), the explained part and the blocks one plot (2), the
*     unexplained part plot 3 and the remainder of individually plot 4;
*   - with -ci-, a spike (rcap) on every row, from value - z * se to value + z * se (se from
*     e(V); z as in check_plot); without -ci-, none, whether or not e(V) is posted;
*   - one label per row, beyond the bar end (or the spike end on the side of the sign),
*     right of positive and left of negative values (above and below with -vertical-), reading
*     the value (3 significant digits; with -format()- exactly the value in that format, a value
*     that reads 0 without its sign) and, except on the total row and when the difference is 0,
*     "(p%)" with p = 100 * value / difference to the rounding (never "-0"); the labels small,
*     vsmall beyond six rows with -vertical-;
*   - with -vertical-, the row names angled (30 degrees) exactly beyond six rows or for a block's
*     name (text) longer than 10 characters;
*   - the block rows named by -blabels()- (as in check_plot) when given;
*   - no legend, no other plot (no markers); the title the one of fit_title (-xsize()- and
*     -ysize()- those of twoway(); not with -notitle-, for a title() given in twoway()); the
*     texts as check_texts has them.
* With -reverse- (makeplot(marginal, reverse), results of sequentially reverse) the stages of the
* reverse order: the blocks from the last to the first, from rev_marg_expl:<block>, and ", reverse
* order" in the subtitle.
capture program drop check_contrib
program define check_contrib
    syntax name(name=G) [, ci VERTical TOLerance(real 1e-10) SUBtitle(string) VTItle(string) REVerse ///
        FORmat(string) BLAbels(string asis) XSize(real -1) YSize(real -1) NOTItle]
    local tol = `tolerance'
    local ns = e(N_stages)
    local blocks `e(blocks)'
    local rp
    if "`reverse'" != "" {
        assert "`e(reverse)'" != ""
        local rp rev_
        local blocks
        forvalues s = `ns'(-1)1 {
            local blocks `blocks' `: word `s' of `e(blocks)''
        }
    }
    * individually: a row for the remainder after the blocks (N_stages counts it)
    local variant = ("`e(stagetype)'" == "individually")
    if `variant' local blocks `blocks' _remainder
    local nr = `ns' + 3
    assert `: word count `blocks'' + 3 == `nr'
    local D = e(difference)
    local Dtxt = strtrim(string(`D', "%9.3g"))
    crit z
    tempname b V R
    matrix `b' = e(b)
    if "`ci'" != "" matrix `V' = e(V)
    local names overall:difference overall:explained
    blabel_rows `"`blocks'"' `"`blabels'"'
    local rl_0 Expl.
    local rl_m1 Total diff
    local rl_`=`ns' + 1' Unexpl.
    local len 0
    foreach bn of local blocks {
        local names `names' `rp'marg_expl:`bn'
    }
    forvalues k = 1/`ns' {
        local len = max(`len', `rn_`k'')
    }
    local names `names' overall:unexplained
    * R: 1 value, 2 se with -ci-, 3 row position, 4 bars seen, 5 the plot expected to draw the
    * row's bar
    matrix `R' = J(`nr', 5, 0)
    local cax = cond("`vertical'" != "", "xaxis1", "yaxis1")
    forvalues r = 1/`nr' {
        local nm : word `r' of `names'
        local c = colnumb(`b', "`nm'")
        assert `c' < .
        matrix `R'[`r', 1] = `b'[1, `c']
        matrix `R'[`r', 2] = .
        if "`ci'" != "" {
            local c = colnumb(`V', "`nm'")
            assert `V'[`c', `c'] > 0 & `V'[`c', `c'] < .
            matrix `R'[`r', 2] = sqrt(`V'[`c', `c'])
        }
        matrix `R'[`r', 5] = cond(`r' == 1, 1, cond(`r' == `nr', 3, cond(`variant' & `r' == `nr' - 1, 4, 2)))
        local t `"`.`G'.`cax'.major.dlg_ticklabel[`r']'"'
        local k = cond(`r' == 1, "m1", string(`r' - 2))
        assert `"`t'"' == `"`rl_`k''"'
        matrix `R'[`r', 3] = `.`G'.`cax'.major.dlg_tickpos[`r']'
        assert abs(`R'[`r', 3] - cond(`r' == 1, 1, cond(`r' < `nr', `r' + .2, `r' + .4))) < 1e-9
    }
    local nci = cond("`ci'" != "", `nr', 0)
    if "`vertical'" != "" {
        assert `.`G'.xaxis1.style.majorstyle.tickangle.val' == cond(`nr' > 6 | `len' > 10, 30, 0)
    }
    local horizontal = cond("`vertical'" != "", "", "horizontal")
    local vax = cond("`vertical'" != "", "yaxis1", "xaxis1")
    if `"`subtitle'"' != "" assert `"`.`G'.subtitle.text[1]'"' == `"`subtitle'"'
    if `"`vtitle'"' != "" assert `"`.`G'.`vax'.title.text[1]'"' == `"`vtitle'"'
    check_texts `G', view(marginal) `horizontal' `reverse' `ci'
    if "`notitle'" == "" {
        fit_title ttl `xsize' `ysize'
        text_lines `G'.title "`ttl'"
    }
    assert "`.`G'.legend.draw_view.setting'" == "no"
    local lpos = cond("`vertical'" != "", 12, 3)
    local lneg = cond("`vertical'" != "", 6, 9)
    local lsz = cond("`vertical'" != "" & `nr' > 6, "vsmall", "small")
    local np 0
    forvalues k = 1/30 {
        if "`.`G'.plotregion1.plot`k'.classname'" == "" continue, break
        local np `k'
    }

    preserve
    local nbar 0
    local ncap 0
    local nlab 0
    forvalues k = 1/`np' {
        plot_data `G' `k'
        local islab = "`v_type'" == "scatter" & "`v_l'" != ""
        * only bars, spikes and labels: no markers
        assert inlist("`v_type'", "rbar", "rcap") | `islab'
        * rbar and rcap: (y, y2, x) = (from, to, row position); labels: (y, x) = (row, anchor),
        * (anchor, row) with -vertical-
        local slot = cond(`islab' & "`vertical'" == "", "`v_y'", "`v_x'")
        local lv = cond("`vertical'" != "", "`v_y'", "`v_x'")
        if `islab' assert "`v_lsize'" == "`lsz'"
        forvalues i = 1/`=_N' {
            local row 0
            forvalues r = 1/`nr' {
                if reldif(`slot'[`i'], `R'[`r', 3]) < 1e-12 local row `r'
            }
            assert `row' > 0
            local val = `R'[`row', 1]
            local se = `R'[`row', 2]
            if "`v_type'" == "rbar" {
                assert `v_y'[`i'] == 0
                * the total difference is the dashed bar
                assert ("`.`G'.plotregion1.plot`k'.style.area.linestyle.pattern.snm'" == "dash") == (`row' == 1)
                assert reldif(`v_y2'[`i'], `val') < `tol'
                assert `k' == `R'[`row', 5]
                matrix `R'[`row', 4] = `R'[`row', 4] + 1
                local ++nbar
            }
            else if "`v_type'" == "rcap" {
                assert `se' < .
                assert reldif(`v_y'[`i'], `val' - `z' * `se') < `tol'
                assert reldif(`v_y2'[`i'], `val' + `z' * `se') < `tol'
                local ++ncap
            }
            else {
                local anchor = `val' + cond(`se' < ., cond(`val' >= 0, 1, -1) * `z' * `se', 0)
                assert reldif(`lv'[`i'], `anchor') < `tol'
                assert `v_p'[`i'] == cond(`val' >= 0, `lpos', `lneg')
                local l = `v_l'[`i']
                local num `"`l'"'
                local pct
                local k2 = strpos(`"`l'"', " (")
                if `k2' {
                    local num = substr(`"`l'"', 1, `k2' - 1)
                    local pct = substr(`"`l'"', `k2' + 2, .)
                }
                if "`format'" != "" {
                    fmt_lab want `val' `format'
                    assert `"`num'"' == "`want'"
                }
                else assert abs(real(`"`num'"') - `val') <= .005 * abs(`val') + 1e-12
                if `row' > 1 & `D' != 0 {
                    assert substr(`"`pct'"', -2, 2) == "%)"
                    local pv = substr(`"`pct'"', 1, strlen(`"`pct'"') - 2)
                    assert `"`pv'"' != "-0"
                    assert abs(real(`"`pv'"') - 100 * `val' / `D') <= .5
                }
                else assert `"`pct'"' == ""
                local ++nlab
            }
        }
    }
    assert `nbar' == `nr'
    forvalues r = 1/`nr' {
        assert `R'[`r', 4] == 1
    }
    assert `ncap' == `nci'
    assert `nlab' == `nr'
    restore
    di as txt "check_contrib `G': `ns' blocks, difference `Dtxt', spikes on `nci' rows" ///
        cond("`vertical'" != "", ", vertical", "") cond("`reverse'" != "", ", reverse order", "") " -- ok"
end

* ---- check_waterfall: the waterfall graph against e() -----------------------------
* check_waterfall graphname [, horizontal ci reverse outcome(name) tolerance(#) subtitle() vtitle()
*                              format() blabels() notitle]
* The rows are read from the tick labels of the category axis (x, or y with -horizontal-,
* where they run from the top): one per block of e(blocks), then "Expl.", "Unexpl." and "Total
* diff", at increasing positions (the waterfall is drawn for sequentially only). The steps are
* marg_expl:<block> of e(b), with outcome() the outcome's own equations. Asserted:
*   - every row has exactly one bar (rbar); a block runs from the running total of the steps
*     before it to the running total including it (the running totals computed here from e(b),
*     to 1e-10; they are also cum_expl:<block>); Expl. runs from 0 to explained, Unexpl. from
*     explained to difference, Total diff from 0 to difference, with a dashed outline;
*   - with -ci-, a spike (rcap) on every bar, around its end, from end - z se to end + z se,
*     se the standard error in e(V) of the bar's value (its step; z as in check_plot); without
*     -ci-, none, whether or not e(V) is posted;
*   - one label per row beyond the bar's end, or the spike's, on the side of the step
*     (above or right when it is positive), reading the step (3 significant digits, with "+"
*     before a positive block or unexplained step; with -format()- exactly the value in that
*     format, "+" before a positive step that does not read 0, a value that reads 0 without its
*     sign) and, except on the total row and when the difference is 0, "(p%)" with p = 100 *
*     step / difference to the rounding; the labels small, vsmall beyond six bars in the
*     vertical view;
*   - the block rows named by -blabels()- (as in check_plot) when given; in the vertical view
*     the names angled exactly beyond six bars or for a block's name (text) longer than 10
*     characters; "Total diff" on one line;
*   - a connector (pcspike, dashed) between consecutive rows at the level of the first
*     one's end, from its bar's edge to the next bar's edge;
*   - no legend; the title "Decomposition of the difference by raking stage" (not with
*     -notitle-); the texts as check_texts has them.
* With -reverse-, the stages of the reverse order: the blocks from the last to the first, the
* steps rev_marg_expl:<block> (the running totals rev_cum_expl:<block>), and ", reverse order"
* in the subtitle.
capture program drop check_waterfall
program define check_waterfall
    syntax name(name=G) [, HORizontal ci OUTcome(string) TOLerance(real 1e-10) ///
        SUBtitle(string) VTItle(string) FORmat(string) BLAbels(string asis) REVerse NOTItle]
    local tol = `tolerance'
    local eq0 overall
    local pfx
    if "`outcome'" != "" {
        local eq0 `outcome'
        local pfx `outcome'_
    }
    assert "`e(stagetype)'" == "sequentially"
    local blocks `e(blocks)'
    local nb : word count `blocks'
    if "`reverse'" != "" {
        assert "`e(reverse)'" != ""
        local pfx `pfx'rev_
        local blocks
        forvalues s = `nb'(-1)1 {
            local blocks `blocks' `: word `s' of `e(blocks)''
        }
    }
    local nr = `nb' + 3
    crit z
    tempname b V R
    matrix `b' = e(b)
    local D = `b'[1, colnumb(`b', "`eq0':difference")]
    * the rounding threshold of the labels (Display_eps): 1e-13 times the largest mean
    local eps 0
    foreach m in focal_group raked_focal ref_group {
        local eps = max(`eps', 1e-13 * abs(`b'[1, colnumb(`b', "`eq0':`m'")]))
    }
    if "`ci'" != "" matrix `V' = e(V)
    local names
    blabel_rows `"`blocks'"' `"`blabels'"'
    local len 0
    forvalues k = 1/`nb' {
        local len = max(`len', `rn_`k'')
    }
    local rl_`=`nb' + 1' Expl.
    local rl_`=`nb' + 2' Unexpl.
    local rl_`nr' Total diff
    foreach bn of local blocks {
        local names `names' `pfx'marg_expl:`bn'
    }
    local names `names' `eq0':explained `eq0':unexplained `eq0':difference
    * R: 1 value (the step), 2 start, 3 end, 4 se with -ci-, 5 position, 6 bars seen,
    * 7 connectors seen
    matrix `R' = J(`nr', 7, 0)
    local cax = cond("`horizontal'" != "", "yaxis1", "xaxis1")
    if "`horizontal'" == "" {
        assert `.`G'.xaxis1.style.majorstyle.tickangle.val' == cond(`nr' > 6 | `len' > 10, 30, 0)
    }
    forvalues r = 1/`nr' {
        local nm : word `r' of `names'
        local c = colnumb(`b', "`nm'")
        assert `c' < .
        matrix `R'[`r', 1] = `b'[1, `c']
        if `r' <= `nb' {
            * the running total of the steps (in the matrix, in double precision)
            if `r' > 1 matrix `R'[`r', 2] = `R'[`r' - 1, 3]
            matrix `R'[`r', 3] = `R'[`r', 2] + `R'[`r', 1]
            local bn : word `r' of `blocks'
            assert reldif(`R'[`r', 3], `b'[1, colnumb(`b', "`pfx'cum_expl:`bn'")]) < 1e-10
        }
        else if `r' == `nb' + 1 {
            matrix `R'[`r', 3] = `R'[`r', 1]
            * the steps add up to the explained part
            assert reldif(`R'[`nb', 3], `R'[`r', 1]) < 1e-10
        }
        else if `r' == `nb' + 2 {
            matrix `R'[`r', 2] = `R'[`r' - 1, 3]
            matrix `R'[`r', 3] = `b'[1, colnumb(`b', "`eq0':difference")]
        }
        else matrix `R'[`r', 3] = `R'[`r', 1]
        matrix `R'[`r', 4] = .
        if "`ci'" != "" {
            local c = colnumb(`V', "`nm'")
            assert `V'[`c', `c'] > 0 & `V'[`c', `c'] < .
            matrix `R'[`r', 4] = sqrt(`V'[`c', `c'])
        }
        local t `"`.`G'.`cax'.major.dlg_ticklabel[`r']'"'
        assert `"`t'"' == `"`rl_`r''"'
        matrix `R'[`r', 5] = `.`G'.`cax'.major.dlg_tickpos[`r']'
        if `r' > 1 assert `R'[`r', 5] > `R'[`r' - 1, 5]
    }
    local nci = cond("`ci'" != "", `nr', 0)
    if `"`subtitle'"' != "" assert `"`.`G'.subtitle.text[1]'"' == `"`subtitle'"'
    local vax = cond("`horizontal'" != "", "xaxis1", "yaxis1")
    if `"`vtitle'"' != "" assert `"`.`G'.`vax'.title.text[1]'"' == `"`vtitle'"'
    check_texts `G', view(waterfall) `horizontal' outcome(`outcome') `reverse' `ci'
    if "`notitle'" == "" text_lines `G'.title "Decomposition of the difference by raking stage"
    assert "`.`G'.legend.draw_view.setting'" == "no"
    local lpos = cond("`horizontal'" != "", 3, 12)
    local lneg = cond("`horizontal'" != "", 9, 6)
    local lsz = cond("`horizontal'" == "" & `nr' > 6, "vsmall", "small")
    local np 0
    forvalues k = 1/30 {
        if "`.`G'.plotregion1.plot`k'.classname'" == "" continue, break
        local np `k'
    }
    preserve
    local nbar 0
    local ncap 0
    local nlab 0
    local ncon 0
    forvalues k = 1/`np' {
        plot_data `G' `k'
        local islab = "`v_type'" == "scatter" & "`v_l'" != ""
        local iscon = "`v_type'" == "pcspike"
        assert inlist("`v_type'", "rbar", "rcap", "scatter", "pcspike")
        if `islab' assert "`v_lsize'" == "`lsz'"
        * rbar and rcap: (y, y2, x) = (from, to, row position), with -horizontal- too; labels:
        * (y, x) = (value, row), (row, value) with -horizontal-; connectors: (y, x, y2, x2) =
        * (level, from, level, to), (from, level, to, level) with -horizontal-
        local slot `v_x'
        local val `v_y'
        if `islab' & "`horizontal'" != "" {
            local slot `v_y'
            local val `v_x'
        }
        if `iscon' {
            local slot = cond("`horizontal'" != "", "`v_y'", "`v_x'")
            local val = cond("`horizontal'" != "", "`v_x'", "`v_y'")
            local slot2 = cond("`horizontal'" != "", "`v_y2'", "`v_x2'")
            local val2 = cond("`horizontal'" != "", "`v_x2'", "`v_y2'")
            assert "`.`G'.plotregion1.plot`k'.style.line.pattern.snm'" == "shortdash"
        }
        forvalues i = 1/`=_N' {
            local row 0
            forvalues r = 1/`nr' {
                if `iscon' {
                    * a connector leaves the bar of its row at its edge (+ .3)
                    if abs(`slot'[`i'] - .3 - `R'[`r', 5]) < 1e-9 local row `r'
                }
                else if reldif(`slot'[`i'], `R'[`r', 5]) < 1e-12 local row `r'
            }
            assert `row' > 0
            local step = `R'[`row', 1]
            local end = `R'[`row', 3]
            local se = `R'[`row', 4]
            if "`v_type'" == "rbar" {
                assert reldif(`v_y'[`i'], `R'[`row', 2]) < `tol'
                assert reldif(`v_y2'[`i'], `end') < `tol'
                assert ("`.`G'.plotregion1.plot`k'.style.area.linestyle.pattern.snm'" == "dash") == (`row' == `nr')
                matrix `R'[`row', 6] = `R'[`row', 6] + 1
                local ++nbar
            }
            else if "`v_type'" == "rcap" {
                assert `se' < .
                assert reldif(`v_y'[`i'], `end' - `z' * `se') < `tol'
                assert reldif(`v_y2'[`i'], `end' + `z' * `se') < `tol'
                local ++ncap
            }
            else if `iscon' {
                assert `row' < `nr'
                assert reldif(`val'[`i'], `end') < `tol' & reldif(`val2'[`i'], `end') < `tol'
                assert abs(`slot2'[`i'] - .7 - `R'[`row', 5]) < 1e-9
                matrix `R'[`row', 7] = `R'[`row', 7] + 1
                local ++ncon
            }
            else if `islab' {
                local anchor = `end' + cond(`se' < ., cond(`step' >= 0, 1, -1) * `z' * `se', 0)
                assert reldif(`val'[`i'], `anchor') < `tol'
                assert `v_p'[`i'] == cond(`step' >= 0, `lpos', `lneg')
                local l = `v_l'[`i']
                local num `"`l'"'
                local pct
                local k2 = strpos(`"`l'"', " (")
                if `k2' {
                    local num = substr(`"`l'"', 1, `k2' - 1)
                    local pct = substr(`"`l'"', `k2' + 2, .)
                }
                * "+" before a positive step of a block or the unexplained part; a value below
                * the rounding threshold reads 0, with a share of 0
                local zero = abs(`step') < `eps'
                local isstep = `row' <= `nb' + 2 & `row' != `nb' + 1
                if "`format'" != "" {
                    fmt_lab want `=cond(`zero', 0, `step')' `format' `=cond(`isstep', "plus", "")'
                    assert `"`num'"' == "`want'"
                }
                else {
                    local plus = (`isstep' & `step' > 0 & !`zero')
                    assert (substr(`"`num'"', 1, 1) == "+") == `plus'
                    if `zero' assert `"`num'"' == "0"
                    else assert abs(real(subinstr(`"`num'"', "+", "", 1)) - `step') <= .005 * abs(`step') + 1e-12
                }
                if `row' < `nr' & abs(`D') >= `eps' {
                    assert substr(`"`pct'"', -2, 2) == "%)"
                    local pv = substr(`"`pct'"', 1, strlen(`"`pct'"') - 2)
                    assert `"`pv'"' != "-0"
                    if `zero' assert `"`pv'"' == "0"
                    else assert abs(real(`"`pv'"') - 100 * `step' / `D') <= .5
                }
                else assert `"`pct'"' == ""
                local ++nlab
            }
            else {
                di as err "check_waterfall `G': unexpected plot `k' (`v_type')"
                exit 9
            }
        }
    }
    assert `nbar' == `nr' & `ncap' == `nci' & `nlab' == `nr' & `ncon' == `nr' - 1
    forvalues r = 1/`nr' {
        assert `R'[`r', 6] == 1
        assert `R'[`r', 7] == (`r' < `nr')
    }
    restore
    di as txt "check_waterfall `G': `nb' steps, difference " strtrim(string(`D', "%9.3g")) ///
        ", spikes on `nci' rows" cond("`horizontal'" != "", ", horizontal", "") ///
        cond("`reverse'" != "", ", reverse order", "") " -- ok"
end

* ---- check_levels: the waterfall in levels against e() ------------------------------
* check_levels graphname [, horizontal ci reverse outcome(name) tolerance(#) vtitle() g1() g2()
*                           save(matname) notitle format() blabels() tgroups() anyticks]
* The rows, read from the tick labels of the category axis as in check_waterfall: the
* Focal Group mean, one per block of e(blocks), the raked Focal Group mean, "Unexpl." and the
* Reference Group mean (the group rows named -g1()- and -g2()- when these are given; the raked
* row "<g1>, raked", or on two lines "<g1>" "raked" in the vertical view when the names are not
* angled; drawn for sequentially only). The steps are the changes in the mean: a step of a block
* (marg_expl:<block>) or of the unexplained part is minus its value in e(b) (e(diffsign) is
* always focal_group - ref_group). Asserted:
*   - the running levels close: focal_group plus the block steps is raked_focal, and raked_focal
*     plus the unexplained step is ref_group (to 1e-10), which pins the sign of the steps;
*   - every row has exactly one bar (rbar): a block from the running level before it
*     (focal_group plus the steps before it) to the running level including it; the
*     unexplained step from raked_focal to ref_group; the three means from the base to
*     focal_group, raked_focal and ref_group of e(b); the bar ends to -tolerance()-, 1e-10;
*   - the base is 0 when 0 lies between the lowest and highest level or spike end, and
*     otherwise 0 or beyond every level on the side of 0; it is the end of the value axis
*     nearest 0 (0 is on the axis when the base is 0); no note says so;
*   - the row names are angled (30 degrees) in the vertical view exactly beyond seven
*     bars, for a block name longer than 10 characters, or for a group name longer than 7
*     (the group rows sit beside "Unexpl."); the labels vsmall beyond seven bars in the
*     vertical view, else small;
*   - with -ci-, a spike (rcap) on every bar (the three means and every step), around the bar's
*     end, from end - z se to end + z se, se the standard error in e(V) of the row's element (z
*     as in check_plot); without -ci-, none, whether or not e(V) is posted;
*   - with -format()-, every label reads exactly the value in that format (a step with "+"
*     when positive and not reading 0, a value that reads 0 without its sign), and the rule
*     of the means' decimals below is not applied;
*   - the block rows named by -blabels()- (as in check_plot) when given;
*   - the title "From the <group 1> mean to the <group 2> mean, stage by stage" (ratio with
*     ratio(); the groups as -tgroups()- gives them, "g1" "g2", else as the rows of the means
*     are named): on one line up to 56 characters ({c ...} one), with ", stage by stage" only
*     when the whole is at most 56; beyond, the two lines "From the <group 1> mean" and "to the
*     <group 2> mean"; no subtitle, "Reverse order" with -reverse- (check_texts);
*   - one label per row beyond the bar's end (or the spike's) on the side of the step,
*     or of the mean from the base (above or right when positive): a mean reads its value
*     (at least 4 significant digits, to half a unit of its last digit, no e-notation
*     below 1e15 unless the mean would need more than 20 decimals, thousands separated
*     exactly when a mean reaches 10,000), and means whose gap is at least the rounding
*     threshold read differently (unless the smallest gap needs more than 10 significant
*     digits); the three means read like a table column, each with the most decimals any
*     of them shows but not more than its own (4 significant digits or the gap's decimals;
*     a mean that reads 0 left out), and no more than needed (beyond the gap's decimals
*     some mean with the most decimals has a last digit that is not 0); a step reads its change
*     ("+" when positive, 3 significant digits) and, unless the difference is below the
*     rounding threshold, "(p%)" with p = 100 * value in e(b) / difference (the share of the
*     stage table, so on the way down a fall has a positive share); a step below the
*     threshold (1e-13 times the largest mean, as in Display_eps) reads "0" and "(0%)";
*   - a connector (pcspike) between consecutive rows at the level of the first one's end;
*   - no legend; the title (not with -notitle-, for a title() given in twoway()); no note.
* With -reverse- (makeplot(waterfalllevels, reverse)) the stages of the reverse order: the blocks
* from the last to the first, the steps from rev_marg_expl:<block> (the raked mean is the same).
* -save()- stores the bars' starts and ends and the labels' positions in a matrix, and the
* labels in a local, to read them or to compare two graphs.
capture program drop check_levels
program define check_levels
    syntax name(name=G) [, HORizontal ci OUTcome(string) TOLerance(real 1e-10) ///
        VTItle(string) g1(string) g2(string) SAVe(name) NOTItle FORmat(string) ///
        BLAbels(string asis) TGRoups(string asis) ANYticks REVerse]
    local tol = `tolerance'
    local eq0 overall
    local pfx
    if "`outcome'" != "" {
        local eq0 `outcome'
        local pfx `outcome'_
    }
    assert "`e(stagetype)'" == "sequentially"
    local blocks `e(blocks)'
    local nb : word count `blocks'
    if "`reverse'" != "" {
        assert "`e(reverse)'" != ""
        local pfx `pfx'rev_
        local blocks
        forvalues s = `nb'(-1)1 {
            local blocks `blocks' `: word `s' of `e(blocks)''
        }
    }
    local nr = `nb' + 4
    assert "`e(diffsign)'" == "focal_group - ref_group"
    crit z
    tempname b V R
    matrix `b' = e(b)
    if "`ci'" != "" matrix `V' = e(V)
    local g1v = `b'[1, colnumb(`b', "`eq0':focal_group")]
    local r1v = `b'[1, colnumb(`b', "`eq0':raked_focal")]
    local g2v = `b'[1, colnumb(`b', "`eq0':ref_group")]
    local D = `b'[1, colnumb(`b', "`eq0':difference")]
    local mx = max(abs(`g1v'), abs(`r1v'), abs(`g2v'))
    local eps = 1e-13 * `mx'
    local names `eq0':focal_group
    foreach bn of local blocks {
        local names `names' `pfx'marg_expl:`bn'
    }
    local names `names' `eq0':raked_focal `eq0':unexplained `eq0':ref_group
    * R: 1 value in e(b), 2 start (. for a mean: the base), 3 end, 4 se with -ci-, 5 position,
    * 6 bars seen, 7 connectors seen, 8 the step (. for a mean)
    matrix `R' = J(`nr', 8, 0)
    local cax = cond("`horizontal'" != "", "yaxis1", "xaxis1")
    local ang = 0
    if "`horizontal'" == "" local ang = `.`G'.xaxis1.style.majorstyle.tickangle.val'
    blabel_rows `"`blocks'"' `"`blabels'"'
    * the rows of the groups: horizontal, tp_wrap to 20 characters up to nine rows (the raked row
    * with ", raked" after the last line), one line cut to 20 characters beyond; vertical and
    * upright, the groups as typed and the raked row "<g1>" "raked"; angled, the groups cut to 18
    * characters up to five rows, 15 up to eight, 13 up to twelve, else 12 (tp_fit), the raked row
    * "<g1 cut>, raked"
    ptexts `outcome'
    if "`horizontal'" != "" & `nr' <= 9 {
        mata: st_local("eg1", tp_wrap(st_local("gr1"), 20, ""))
        mata: st_local("egr", tp_wrap(st_local("gr1"), 20, ", raked"))
        mata: st_local("eg2", tp_wrap(st_local("gr2"), 20, ""))
    }
    else if `ang' == 0 & "`horizontal'" == "" {
        mata: st_local("eg1", tp_esc(st_local("gr1")))
        mata: st_local("eg2", tp_esc(st_local("gr2")))
        local egr `"`"`macval(eg1)'"' "raked""'
    }
    else {
        local wl = cond("`horizontal'" != "", 20, cond(`nr' <= 5, 18, cond(`nr' <= 8, 15, cond(`nr' <= 12, 13, 12))))
        mata: st_local("eg1", tp_esc(tp_fit(st_local("gr1"), `wl')))
        mata: st_local("egr", tp_esc(tp_fit(st_local("gr1"), `wl') + ", raked"))
        mata: st_local("eg2", tp_esc(tp_fit(st_local("gr2"), `wl')))
    }
    * g1(), g2(): the groups' texts, as typed with tp_esc
    forvalues g = 1/2 {
        if `"`g`g''"' != "" {
            mata: st_local("x", tp_esc(st_local("gr`g'")))
            assert `"`macval(x)'"' == `"`g`g''"'
        }
    }
    forvalues r = 1/`nr' {
        local nm : word `r' of `names'
        local c = colnumb(`b', "`nm'")
        assert `c' < .
        matrix `R'[`r', 1] = `b'[1, `c']
        matrix `R'[`r', 2] = .
        matrix `R'[`r', 8] = .
        local mean = inlist(`r', 1, `nb' + 2, `nr')
        if `mean' matrix `R'[`r', 3] = `R'[`r', 1]
        else if `r' <= `nb' + 1 {
            * the running level (in the matrix, in double precision)
            matrix `R'[`r', 8] = -`R'[`r', 1]
            matrix `R'[`r', 2] = `R'[`r' - 1, 3]
            matrix `R'[`r', 3] = `R'[`r', 2] + `R'[`r', 8]
        }
        else {
            matrix `R'[`r', 8] = -`R'[`r', 1]
            matrix `R'[`r', 2] = `r1v'
            matrix `R'[`r', 3] = `g2v'
        }
        matrix `R'[`r', 4] = .
        if "`ci'" != "" {
            * a spike on every bar: the means and the steps
            local c = colnumb(`V', "`nm'")
            assert `V'[`c', `c'] > 0 & `V'[`c', `c'] < .
            matrix `R'[`r', 4] = sqrt(`V'[`c', `c'])
        }
        * in compound quotes: a label of two lines starts with one
        local t `"`.`G'.`cax'.major.dlg_ticklabel[`r']'"'
        if `r' > 1 & `r' <= `nb' + 1 {
            assert `"`t'"' == `"`rl_`=`r' - 1''"'
        }
        else if `r' == `nb' + 3 assert `"`t'"' == "Unexpl."
        else if `r' == 1 assert `"`t'"' == `"`macval(eg1)'"'
        else if `r' == `nr' assert `"`t'"' == `"`macval(eg2)'"'
        else if `r' == `nb' + 2 assert `"`t'"' == `"`macval(egr)'"'
        matrix `R'[`r', 5] = `.`G'.`cax'.major.dlg_tickpos[`r']'
        if `r' > 1 assert `R'[`r', 5] > `R'[`r' - 1, 5]
    }
    local nci = cond("`ci'" != "", `nr', 0)
    * the angle rule of the vertical view: beyond seven bars, a block name longer than 10
    * characters, a group name longer than 7
    if "`horizontal'" == "" {
        local blen 0
        forvalues k = 1/`nb' {
            local blen = max(`blen', `rn_`k'')
        }
        mata: st_local("glen", strofreal(max((ustrlen(st_local("gr1")), ustrlen(st_local("gr2"))))))
        assert `ang' == cond(`nr' > 7 | `blen' > 10 | `glen' > 7, 30, 0)
    }
    local lsz = cond("`horizontal'" == "" & `nr' > 7, "vsmall", "small")
    * the levels close: to the raked mean after the blocks, to the Reference Group mean after the
    * unexplained step
    assert reldif(`R'[`nb' + 1, 3], `r1v') < 1e-10
    assert reldif(`r1v' + `R'[`nb' + 3, 8], `g2v') < 1e-10
    local vax = cond("`horizontal'" != "", "xaxis1", "yaxis1")
    if `"`vtitle'"' != "" assert `"`.`G'.`vax'.title.text[1]'"' == `"`vtitle'"'
    * no subtitle (the title names the groups), "Reverse order" for the reverse order; no note
    check_texts `G', view(levels) `horizontal' outcome(`outcome') `reverse' `ci'
    local ym = cond("`e(ratio)'" != "", "ratio", "mean")
    * the groups as typed (or as -tgroups()- gives them, as typed), in parentheses in the title
    * for an unlabelled numeric group; a line of the two-line title cuts a group to 42
    * characters (tp_fit)
    if `"`tgroups'"' != "" {
        mata: st_local("sr1", tokens(st_local("tgroups"))[1])
        mata: st_local("sr2", tokens(st_local("tgroups"))[2])
    }
    mata: st_local("tlen", strofreal(ustrlen("From the " + st_local("sr1") + " `ym' to the " + ///
        st_local("sr2") + " `ym'")))
    if "`notitle'" == "" & `tlen' <= 56 {
        mata: st_local("ttl", tp_esc("From the " + st_local("sr1") + " `ym' to the " + st_local("sr2") + " `ym'"))
        if `tlen' + 16 <= 56 {
            local ttl `"`macval(ttl)', stage by stage"'
        }
        text_lines `G'.title `"`macval(ttl)'"'
    }
    else if "`notitle'" == "" {
        mata: st_local("ttl1", tp_esc("From the " + tp_fit(st_local("sr1"), 42) + " `ym'"))
        mata: st_local("ttl2", tp_esc("to the " + tp_fit(st_local("sr2"), 42) + " `ym'"))
        text_lines `G'.title `"`macval(ttl1)'"' `"`macval(ttl2)'"'
    }
    assert "`.`G'.legend.draw_view.setting'" == "no"
    * the lowest and highest level or spike end
    local lo = .
    local hi = .
    forvalues r = 1/`nr' {
        foreach e in 2 3 {
            local v = `R'[`r', `e']
            if `v' < . {
                local lo = min(`lo', `v')
                local hi = max(`hi', `v')
            }
        }
        if `R'[`r', 4] < . {
            local lo = min(`lo', `R'[`r', 3] - `z' * `R'[`r', 4])
            local hi = max(`hi', `R'[`r', 3] + `z' * `R'[`r', 4])
        }
    }
    local lpos = cond("`horizontal'" != "", 3, 12)
    local lneg = cond("`horizontal'" != "", 9, 6)
    local np 0
    forvalues k = 1/30 {
        if "`.`G'.plotregion1.plot`k'.classname'" == "" continue, break
        local np `k'
    }
    preserve
    local nbar 0
    local ncap 0
    local nlab 0
    local ncon 0
    local base .
    forvalues k = 1/`np' {
        plot_data `G' `k'
        local islab = "`v_type'" == "scatter" & "`v_l'" != ""
        local iscon = "`v_type'" == "pcspike"
        assert inlist("`v_type'", "rbar", "rcap", "scatter", "pcspike")
        if `islab' assert "`v_lsize'" == "`lsz'"
        local slot `v_x'
        local val `v_y'
        if `islab' & "`horizontal'" != "" {
            local slot `v_y'
            local val `v_x'
        }
        if `iscon' {
            local slot = cond("`horizontal'" != "", "`v_y'", "`v_x'")
            local val = cond("`horizontal'" != "", "`v_x'", "`v_y'")
            local slot2 = cond("`horizontal'" != "", "`v_y2'", "`v_x2'")
            local val2 = cond("`horizontal'" != "", "`v_x2'", "`v_y2'")
        }
        forvalues i = 1/`=_N' {
            local row 0
            forvalues r = 1/`nr' {
                if `iscon' {
                    if abs(`slot'[`i'] - .3 - `R'[`r', 5]) < 1e-9 local row `r'
                }
                else if reldif(`slot'[`i'], `R'[`r', 5]) < 1e-12 local row `r'
            }
            assert `row' > 0
            local step = `R'[`row', 8]
            local end = `R'[`row', 3]
            local se = `R'[`row', 4]
            if "`v_type'" == "rbar" {
                if `step' < . assert reldif(`v_y'[`i'], `R'[`row', 2]) < `tol'
                else {
                    * the means' bars share one base
                    if `base' == . local base = `v_y'[`i']
                    assert `v_y'[`i'] == `base'
                }
                assert reldif(`v_y2'[`i'], `end') < `tol'
                matrix `R'[`row', 6] = `R'[`row', 6] + 1
                local ++nbar
            }
            else if "`v_type'" == "rcap" {
                assert `se' < .
                assert reldif(`v_y'[`i'], `end' - `z' * `se') < `tol'
                assert reldif(`v_y2'[`i'], `end' + `z' * `se') < `tol'
                local ++ncap
            }
            else if `iscon' {
                assert `row' < `nr'
                assert reldif(`val'[`i'], `end') < `tol' & reldif(`val2'[`i'], `end') < `tol'
                assert abs(`slot2'[`i'] - .7 - `R'[`row', 5]) < 1e-9
                matrix `R'[`row', 7] = `R'[`row', 7] + 1
                local ++ncon
            }
            else if `islab' {
                local lab`row' = `v_l'[`i']
                local lx`row' = `val'[`i']
                local lp`row' = `v_p'[`i']
                local ++nlab
            }
            else {
                di as err "check_levels `G': unexpected plot `k' (`v_type')"
                exit 9
            }
        }
    }
    restore
    assert `nbar' == `nr' & `ncap' == `nci' & `nlab' == `nr' & `ncon' == `nr' - 1
    forvalues r = 1/`nr' {
        assert `R'[`r', 6] == 1
        assert `R'[`r', 7] == (`r' < `nr')
    }
    * the base: 0 when 0 is among the levels and spikes, else 0 or beyond them on the side of
    * 0; the end of the value axis nearest 0, or 0 on the axis
    assert `base' < .
    if `lo' <= 0 & `hi' >= 0 assert `base' == 0
    else if `lo' > 0 assert `base' >= 0 & `base' <= `lo'
    else assert `base' <= 0 & `base' >= `hi'
    local vsc = cond("`horizontal'" != "", "xscale", "yscale")
    local amin = `.`G'.plotregion1.`vsc'.curmin'
    local amax = `.`G'.plotregion1.`vsc'.curmax'
    if `base' == 0 assert `amin' <= 0 & `amax' >= 0
    else if `base' > 0 assert reldif(`amin', `base') < 1e-9
    else assert reldif(`amax', `base') < 1e-9
    * at most eight ticks on the value axis (beyond eight, the axis extended for the labels takes
    * the step of the wider span); not with a ylabel(#10) of the user's
    if "`anyticks'" == "" assert `.`G'.`vax'.major.dlg_tickpos.arrnels' <= 8
    * the labels
    forvalues r = 1/`nr' {
        local step = `R'[`r', 8]
        local end = `R'[`r', 3]
        local se = `R'[`r', 4]
        local dir = cond(`step' < ., `step', `end' - `base')
        local anchor = `end' + cond(`se' < ., cond(`dir' >= 0, 1, -1) * `z' * `se', 0)
        assert reldif(`lx`r'', `anchor') < `tol'
        assert `lp`r'' == cond(`dir' >= 0, `lpos', `lneg')
        local l `"`lab`r''"'
        local num `"`l'"'
        local pct
        local k2 = strpos(`"`l'"', " (")
        if `k2' {
            local num = substr(`"`l'"', 1, `k2' - 1)
            local pct = substr(`"`l'"', `k2' + 2, .)
        }
        if `step' == . & "`format'" != "" {
            * format(): the mean in that format
            assert `"`pct'"' == ""
            fmt_lab want `=cond(abs(`end') < `eps', 0, `end')' `format'
            assert `"`num'"' == "`want'"
            local mlab`r' `"`num'"'
        }
        else if `step' < . & "`format'" != "" {
            * format(): the step in that format, "+" when positive and not reading 0
            fmt_lab want `=cond(abs(`step') < `eps', 0, `step')' `format' plus
            assert `"`num'"' == "`want'"
            if abs(`D') >= `eps' {
                local pv = substr(`"`pct'"', 1, strlen(`"`pct'"') - 2)
                assert substr(`"`pct'"', -2, 2) == "%)" & `"`pv'"' != "-0"
                if abs(`step') < `eps' assert `"`pv'"' == "0"
                else assert abs(real(`"`pv'"') - 100 * `R'[`r', 1] / `D') <= .5
            }
            else assert `"`pct'"' == ""
        }
        else if `step' == . {
            assert substr(`"`num'"', 1, 1) != "+" & `"`pct'"' == ""
            assert abs(real(subinstr(`"`num'"', ",", "", .)) - `end') <= .0005 * abs(`end') + 1e-12
            * thousands separated when a mean reaches 10,000; to half a unit of the last digit
            * shown (e-notation, which only a mean needing more than 20 decimals may take, is
            * checked with the table column below)
            if `mx' < 1e15 & !strpos(`"`num'"', "e") {
                if `mx' < 9999.5 assert !strpos(`"`num'"', ",")
                if abs(`end') >= 9999.5 assert strpos(`"`num'"', ",") > 0
                local k3 = strpos(`"`num'"', ".")
                local nd = cond(`k3', strlen(`"`num'"') - `k3', 0)
                assert abs(real(subinstr(`"`num'"', ",", "", .)) - `end') <= ///
                    .5 * 10^(-`nd') * (1 + 1e-9) + `eps'
            }
            local mlab`r' `"`num'"'
        }
        else if abs(`step') < `eps' {
            * below the rounding threshold: 0, with a share of 0 (as in the waterfall)
            assert `"`num'"' == "0"
            if abs(`D') >= `eps' assert `"`pct'"' == "0%)"
            else assert `"`pct'"' == ""
        }
        else {
            assert (substr(`"`num'"', 1, 1) == "+") == (`step' > 0)
            assert abs(real(subinstr(`"`num'"', "+", "", 1)) - `step') <= .005 * abs(`step') + 1e-12
            if abs(`D') >= `eps' {
                assert substr(`"`pct'"', -2, 2) == "%)"
                local pv = substr(`"`pct'"', 1, strlen(`"`pct'"') - 2)
                assert `"`pv'"' != "-0"
                assert abs(real(`"`pv'"') - 100 * `R'[`r', 1] / `D') <= .5
            }
            else assert `"`pct'"' == ""
        }
    }
    * horizontal: every label fits within the value axis by the ado's measure: a label of
    * s characters at its anchor spans (s + 1) * cw * (amax - amin) away from the bar, with cw =
    * .95 / (83 - rn), rn the longest row name as the ado counts it (11, the groups cut to 20 and
    * 7 more for the raked row, the blocks' texts); the ado extends the axis by whole ticks until
    * that holds (the step following the wider span beyond eight ticks), except past an end at 0
    * (the levels then all on one side of it)
    if "`horizontal'" != "" {
        mata: st_local("rn", strofreal(max((11, min((ustrlen(st_local("gr1")), 20)) + 7, ///
            min((ustrlen(st_local("gr2")), 20))))))
        forvalues k = 1/`nb' {
            local rn = max(`rn', `rn_`k'')
        }
        local cw = .95 / max(83 - `rn', 20)
        forvalues r = 1/`nr' {
            local need = (strlen(`"`lab`r''"') + 1) * `cw' * (`amax' - `amin')
            if `lp`r'' == `lneg' assert `lx`r'' - `need' >= `amin' - 1e-9 | (`amin' == 0 & `lo' >= 0)
            else assert `lx`r'' + `need' <= `amax' + 1e-9 | (`amax' == 0 & `hi' <= 0)
        }
    }
    * means whose gap is at least the rounding threshold read differently, unless the smallest
    * such gap needs more than 10 significant digits of the largest mean
    local m1 1
    local m2 = `nb' + 2
    local m3 = `nr'
    local gp = .
    foreach p in "1 2" "2 3" "1 3" {
        local a : word 1 of `p'
        local c : word 2 of `p'
        local q = abs(`R'[`m`c'', 1] - `R'[`m`a'', 1])
        if `q' >= `eps' local gp = min(`gp', `q')
    }
    local dg 0
    if `gp' < . & `mx' < 1e15 & "`format'" == "" {
        if floor(log10(`mx')) + 1 + max(0, -floor(log10(`gp'))) <= 10 {
            local dg = max(0, -floor(log10(`gp')))
            foreach p in "1 2" "2 3" "1 3" {
                local a : word 1 of `p'
                local c : word 2 of `p'
                if abs(`R'[`m`c'', 1] - `R'[`m`a'', 1]) >= `eps' {
                    assert `"`mlab`m`a'''"' != `"`mlab`m`c'''"'
                }
            }
        }
    }
    * a mean's own decimals: d = 4 significant digits or the gap's decimals; e-notation exactly
    * where d is more than 20 (a mean below 1e-17, or a gap below 1e-20). The means read like a
    * table column: D, the most decimals any of them shows, but no mean more than its own d, so
    * each shows min(D, d) (a mean that reads 0 left out); and no more than needed: beyond the
    * gap's decimals, some mean that shows D ends in a digit that is not 0 (not with format())
    if `mx' < 1e15 & "`format'" == "" {
        local ndc 0
        foreach m in `m1' `m2' `m3' {
            local l `"`mlab`m''"'
            if `"`l'"' == "0" continue
            local d`m' = max(0, 3 - floor(log10(abs(`R'[`m', 1]))), `dg')
            assert (strpos(`"`l'"', "e") > 0) == (`d`m'' > 20)
            if strpos(`"`l'"', "e") continue
            local k3 = strpos(`"`l'"', ".")
            local nd`m' = cond(`k3', strlen(`"`l'"') - `k3', 0)
            local ndc = max(`ndc', `nd`m'')
        }
        local last0 1
        foreach m in `m1' `m2' `m3' {
            local l `"`mlab`m''"'
            if `"`l'"' == "0" | strpos(`"`l'"', "e") continue
            assert `nd`m'' == min(`ndc', `d`m'')
            if `nd`m'' == `ndc' & substr(`"`l'"', -1, 1) != "0" local last0 0
        }
        if `ndc' > `dg' assert !`last0'
    }
    * -save()-: the bars' starts and ends in matrix <save>, the labels in local <save>_labs
    if "`save'" != "" {
        matrix `save' = `R'[1..., 2..3]
        local all
        forvalues r = 1/`nr' {
            if `R'[`r', 8] == . matrix `save'[`r', 1] = `base'
            local all `"`all'|`lab`r''"'
        }
        c_local `save'_labs `"`all'"'
    }
    di as txt "check_levels `G': `nb' steps, from " strtrim(string(`g1v', "%9.4g")) " to " ///
        strtrim(string(`g2v', "%9.4g")) ", base " strtrim(string(`base', "%12.0g")) ///
        ", spikes on `nci' rows" cond("`horizontal'" != "", ", horizontal", "") ///
        cond("`reverse'" != "", ", reverse order", "") " -- ok"
end

* ---- plot_rgb: the fill colour of a plot ----------------------------------------
capture program drop plot_rgb
program define plot_rgb
    args lname G k
    c_local `lname' `.`G'.plotregion1.plot`k'.style.area.shadestyle.color.rgb'
end

* ---- drop_bystage: leave only e(b), as a prefix that does not carry e(bystage) would ----
capture program drop drop_bystage
program define drop_bystage, eclass
    mata: st_matrix("e(bystage)", J(0, 0, .))
end

* ---- set_dfr: post residual degrees of freedom, as svy and jackknife do ------------
capture program drop set_dfr
program define set_dfr, eclass
    ereturn scalar df_r = `1'
end

* ---- same_graph: two graphs that draw the same ------------------------------------------------
* same_graph G1 G2: the graph of the default plot type (G1) against that of the type named (G2):
* the same title, subtitle, axis titles and note, the same tick positions and labels on both axes,
* and the same plots, plot by plot: the class, the plot type, the roles of the variables (plot_data)
* and the serset, observation by observation and variable by variable (cf)
capture program drop same_graph
program define same_graph
    args G1 G2
    foreach o in title subtitle xaxis1.title yaxis1.title note {
        local n1 `.`G1'.`o'.text.arrnels'
        local n2 `.`G2'.`o'.text.arrnels'
        assert "`n1'" == "`n2'"
        if !inlist("`n1'", "", "0") {
            forvalues k = 1/`n1' {
                assert `"`.`G1'.`o'.text[`k']'"' == `"`.`G2'.`o'.text[`k']'"'
            }
        }
    }
    foreach ax in xaxis1 yaxis1 {
        local n1 `.`G1'.`ax'.major.dlg_tickpos.arrnels'
        local n2 `.`G2'.`ax'.major.dlg_tickpos.arrnels'
        assert "`n1'" == "`n2'"
        if !inlist("`n1'", "", "0") {
            forvalues k = 1/`n1' {
                assert "`.`G1'.`ax'.major.dlg_tickpos[`k']'" == "`.`G2'.`ax'.major.dlg_tickpos[`k']'"
                assert `"`.`G1'.`ax'.major.dlg_ticklabel[`k']'"' == `"`.`G2'.`ax'.major.dlg_ticklabel[`k']'"'
            }
        }
    }
    local np 0
    forvalues k = 1/30 {
        local c1 `.`G1'.plotregion1.plot`k'.classname'
        local c2 `.`G2'.plotregion1.plot`k'.classname'
        assert "`c1'" == "`c2'"
        if "`c1'" == "" continue, break
        local np `k'
    }
    assert `np' > 0
    preserve
    forvalues k = 1/`np' {
        local i1 `.`G1'.plotregion1.plot`k'.serset.id'
        local i2 `.`G2'.plotregion1.plot`k'.serset.id'
        assert inlist("`i1'", "", ".") == inlist("`i2'", "", ".")
        if inlist("`i1'", "", ".") continue
        plot_data `G1' `k'
        local r1 `v_type' | `v_y' | `v_y2' | `v_x' | `v_x2' | `v_l' | `v_p' | `v_barsize' | `v_symbol' | `v_lsize'
        local nv1 = c(k)
        local no1 = _N
        tempfile f
        qui save `"`f'"'
        plot_data `G2' `k'
        local r2 `v_type' | `v_y' | `v_y2' | `v_x' | `v_x2' | `v_l' | `v_p' | `v_barsize' | `v_symbol' | `v_lsize'
        assert "`r1'" == "`r2'"
        assert c(k) == `nv1' & _N == `no1'
        cf _all using `"`f'"'
    }
    restore
    di as txt "same_graph `G1' `G2': `np' plots -- ok"
end

* ---- err_msg: a command's return code and a message it printed --------------------------
* err_msg "command" "text" [rc]: the command, run with capture noisily into a text log, returns
* rc (198 by default; 0 for a message that is not an error) and what it printed contains text,
* read as one line (sd_logtext: the lines of the log joined, a line of {p} broken at a space at
* the width of the log with a space, a longer line continued after "> " without one, runs of
* spaces as one)
capture mata: mata drop sd_logtext()
mata:
string scalar sd_logtext(string scalar f)
{
    string colvector L
    string scalar s
    real scalar i
    L = cat(f)
    s = ""
    for (i = 1; i <= rows(L); i++) {
        if (substr(L[i], 1, 2) == "> ") s = s + substr(L[i], 3, .)
        else s = s + " " + L[i]
    }
    return(stritrim(s))
}
end
capture program drop err_msg
program define err_msg
    args cmd txt rc0
    if "`rc0'" == "" local rc0 198
    tempfile f
    local lf `"`f'.log"'
    quietly log using `"`lf'"', text replace name(sderrmsg)
    capture noisily `cmd'
    local rc = _rc
    quietly log close sderrmsg
    mata: st_local("hit", strofreal(strpos(sd_logtext(st_local("lf")), stritrim(st_local("txt"))) > 0))
    capture erase `"`lf'"'
    if `rc' != `rc0' | !`hit' {
        di as err `"err_msg: `cmd' returned `rc' (`rc0' expected); "`txt'" found `hit' times"'
        exit 9
    }
end

* ---- spec_refused: a plotspec refused before the estimation and at replay ----------------
* spec_refused "plotspec" "message" ["message at replay"]: seldecomp wage (demo: race collgrad),
* by(married) sequentially fixothermargins makeplot(plotspec) returns 198 with the message before
* the estimation (e(b) still the matrix b_A, e(cmdline) still $SD_CLA), and seldecomp,
* makeplot(plotspec) returns 198 with the message (the one at replay when it differs)
capture program drop spec_refused
program define spec_refused
    args spec m1 m2
    if `"`macval(m2)'"' == "" local m2 `"`macval(m1)'"'
    * another decomposition (one block), which would replace e(b) if it ran
    err_msg `"seldecomp wage (demo: race collgrad), by(married) sequentially fixothermargins makeplot(`macval(spec)')"' ///
        `"`macval(m1)'"'
    assert colsof(e(b)) == colsof(b_A)
    assert mreldif(e(b), b_A) == 0
    assert `"`e(cmdline)'"' == `"$SD_CLA"'
    err_msg `"seldecomp, makeplot(`macval(spec)')"' `"`macval(m2)'"'
    assert mreldif(e(b), b_A) == 0
    di as txt `"spec_refused: makeplot(`macval(spec)') -- ok"'
end

* ---- data ---------------------------------------------------------------------
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
recast double y x3 wt
confirm string variable x4
assert _N == 2000
do expected_values.do

capture graph drop _all

* ---- makeplot with a graph name -----------------------------------------------
di as txt _n "---- makeplot(stages, twoway(name(sdplot, replace)))"
capture noisily seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially makeplot(stages, twoway(name(sdplot, replace)))
local rc = _rc
if `rc' {
    di as err "seldecomp ... sequentially makeplot(stages, twoway(name(sdplot, replace))) failed with rc `rc'"
    exit `rc'
}
* results are intact after the plot (the plot uses preserve/restore internally); sequentially
* alone leaves the blocks not yet added free: case F of expected_values.do
assert e(N_stages) == 3
assert "`e(stagetype)'" == "sequentially" & "`e(fixothermargins)'" == ""
assert "`e(cmd)'" == "seldecomp"
assert reldif(e(explained), exp_F3_explained) < `tol'
assert reldif(_b[cum_expl:x1], exp_F1_explained) < `tol'
assert reldif(_b[cum_expl:x2], exp_F2_explained) < `tol'
assert _N == 2000
* the graph exists under the requested name
capture noisily graph describe sdplot
local rc = _rc
if `rc' {
    di as err "graph sdplot was not created (rc `rc')"
    exit 9
}
* titles: y has no variable label, group no value label; the stages plot has no title on its
* category axis, and no note
check_plot sdplot, subtitle("(group = 0) - (group = 1)") vtitle("y: difference in means")
no_text sdplot.xaxis1.title
no_text sdplot.note
* without a plot type, sequentially alone draws the waterfall
seldecomp, makeplot(, twoway(name(sdplotw, replace)))
check_waterfall sdplotw, subtitle("(group = 0) - (group = 1)") vtitle("y: difference in means")

* ---- makeplot without a plotspec, with blocks and fixothermargins ------------------
di as txt _n "---- makeplot with blocks and fixothermargins"
capture noisily seldecomp y (edu: x1 x2) x4 [pw=wt], by(group) sequentially fixothermargins makeplot
local rc = _rc
if `rc' {
    di as err "seldecomp ... sequentially fixothermargins makeplot failed with rc `rc'"
    exit `rc'
}
assert e(N_stages) == 2
assert "`e(stagetype)'" == "sequentially" & "`e(fixothermargins)'" == "fixothermargins"
assert "`e(blocks)'" == "edu x4"
* the default graph name is seldecomp, not Graph; the default plot type of sequentially
* is the waterfall (stages names the stages plot)
check_waterfall seldecomp
capture graph describe Graph
assert _rc
seldecomp, makeplot(stages)
check_plot seldecomp
* the next graph without name() in twoway() replaces it (no error), also when twoway() holds
* other options; with name() the name given is used and seldecomp is left as it was
seldecomp, makeplot(marginal)
check_contrib seldecomp
seldecomp, makeplot(waterfall, twoway(name(sdwname, replace)))
check_waterfall sdwname
check_contrib seldecomp
graph drop sdwname
seldecomp, makeplot(waterfall, twoway(title("a title of twoway()")))
check_waterfall seldecomp, notitle
text_lines seldecomp.title "a title of twoway()"
* makeplot() and makeplot(, ...) without a type: the default type, the waterfall
seldecomp, makeplot()
check_waterfall seldecomp
seldecomp, makeplot(, swapaxis)
check_waterfall seldecomp, horizontal

* ---- twoway options are passed through twoway() --------------------------------
di as txt _n "---- makeplot(stages, twoway(title() name() scheme()))"
capture noisily seldecomp y x1 x2 [pw=wt], by(group) sequentially ///
    makeplot(stages, twoway(title("seldecomp test plot") name(sdplot2, replace) scheme(s2color)))
local rc = _rc
if `rc' {
    di as err "makeplot with twoway options failed with rc `rc'"
    exit `rc'
}
capture noisily graph describe sdplot2
local rc = _rc
if `rc' {
    di as err "graph sdplot2 was not created (rc `rc')"
    exit 9
}
assert `"`.sdplot2.title.text[1]'"' == "seldecomp test plot"
check_plot sdplot2, notitle

* ---- makeplot together with the bootstrap -------------------------------------
* bayesboot(): every stage bootstrapped, but the CI caps are drawn only with ci
di as txt _n "---- makeplot with bayesboot"
capture noisily seldecomp y x1 x2 [pw=wt], by(group) sequentially ///
    bayesboot(reps(5) rseed(1) nodots) makeplot(stages, twoway(name(sdplot3, replace)))
local rc = _rc
if `rc' {
    di as err "makeplot with bayesboot failed with rc `rc'"
    exit `rc'
}
assert e(N_reps) == 5
confirm matrix e(V)
check_plot sdplot3
seldecomp, makeplot(stages, ci twoway(name(sdplot3c, replace)))
check_plot sdplot3c, ci

* every stage bootstrapped, makeplot(stages, ci): CI caps at the estimation level, bayesboot(level());
* without cluster() the intervals are normal by default, without e(df_r): z
di as txt _n "---- makeplot(stages, ci) with bayesboot(level(90))"
qui seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially fixothermargins ///
    bayesboot(reps(20) rseed(1) nodots level(90)) makeplot(stages, ci twoway(name(sdci, replace)))
assert e(level) == 90
assert "`e(citypes)'" == "normal" & e(df_r) == . & e(df_bs) == e(N) - 1
check_plot sdci, ci

* bayesboot(cluster()): the default intervals are t, with e(df_bs) = the number of
* clusters less 1 degrees of freedom, posted as e(df_r) (t is the family of the table's tests),
* so the intervals of every plot type are t-based (crit reads e(df_r): t with 59 df here, not z);
* with normal inside bayesboot() z again, without e(df_r), from the same draws (the same e(V))
di as txt _n "---- makeplot(..., ci) with bayesboot(cluster()): t intervals"
qui seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially fixothermargins ///
    bayesboot(reps(20) rseed(1) nodots cluster(cluster)) makeplot(stages, ci twoway(name(sdcit, replace)))
qui tab cluster if e(sample)
assert e(N_clust) == r(r) & e(df_bs) == r(r) - 1 & e(df_r) == e(df_bs)
assert "`e(citypes)'" == "t"
matrix V_t = e(V)
check_plot sdcit, ci
seldecomp, makeplot(stages, swapaxis ci twoway(name(sdcith, replace)))
check_plot sdcith, ci horizontal
seldecomp, makeplot(marginal, ci twoway(name(sdcitc, replace)))
check_contrib sdcitc, ci
seldecomp, makeplot(waterfall, ci twoway(name(sdcitw, replace)))
check_waterfall sdcitw, ci
seldecomp, makeplot(waterfalllevels, ci twoway(name(sdcitl, replace)))
check_levels sdcitl, ci
* the default plot type with ci: the waterfall, t-based
seldecomp, makeplot(, ci twoway(name(sdcitd, replace)))
check_waterfall sdcitd, ci
same_graph sdcitd sdcitw
qui seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially fixothermargins ///
    bayesboot(reps(20) rseed(1) nodots cluster(cluster) normal) makeplot(stages, ci twoway(name(sdcin, replace)))
assert "`e(citypes)'" == "normal" & e(df_r) == . & e(df_bs) == e(N_clust) - 1
assert mreldif(e(V), V_t) == 0
check_plot sdcin, ci
seldecomp, makeplot(waterfall, ci twoway(name(sdcinw, replace)))
check_waterfall sdcinw, ci
matrix drop V_t

* ---- sign case 1: explained and unexplained of the same sign ------------------
di as txt _n "---- same sign"
qui seldecomp y x1 x2 x4 [pw=wt], by(group) sequentially fixothermargins makeplot(stages, twoway(name(sdsame, replace)))
* fixothermargins holds the blocks not yet added at the Focal Group's shares: case G of
* expected_values.do (the stages differ from case F, the full specification does not)
assert "`e(stagetype)'" == "sequentially" & "`e(fixothermargins)'" == "fixothermargins"
assert reldif(_b[cum_expl:x1], exp_G1_explained) < `tol'
assert reldif(_b[cum_expl:x2], exp_G2_explained) < `tol'
assert reldif(e(explained), exp_G3_explained) < `tol'
assert reldif(exp_G1_explained, exp_F1_explained) > 1e-3
matrix S = e(bystage)
forvalues s = 1/3 {
    assert sign(S[`s', colnumb(S, "cum_expl")]) == sign(e(difference))
    assert sign(S[`s', colnumb(S, "cum_unexpl")]) == sign(e(difference))
}
check_plot sdsame
* the waterfall in levels: from the Focal Group mean through the steps to the Reference Group mean
seldecomp, makeplot(waterfalllevels, twoway(name(sdsamel, replace)))
check_levels sdsamel, vtitle("y: mean") g1("group = 0") g2("group = 1")

* ---- after a replication prefix: seldecomp, makeplot draws the point estimates ----
* (inside the prefix every replicate would redraw it and the last replicate would stay)
di as txt _n "---- bootstrap: then seldecomp, makeplot"
qui bootstrap, reps(3) seed(1) nodots: seldecomp y x1 x2 x4, by(group) sequentially fixothermargins
assert "`e(prefix)'" == "bootstrap"
seldecomp, makeplot(stages, ci twoway(name(sdprefix, replace)))
check_plot sdprefix, ci tolerance(1e-12)
* without ci no caps, although e(V) holds the prefix's variances
seldecomp, makeplot(stages, twoway(name(sdprefix0, replace)))
check_plot sdprefix0, tolerance(1e-12)
* the same when only e(b) survives the prefix: the stage values are rebuilt from e(b)
drop_bystage
capture confirm matrix e(bystage)
assert _rc == 111
seldecomp, makeplot(stages, ci twoway(name(sdprefix2, replace)))
check_plot sdprefix2, ci tolerance(1e-12)
* the marginal plot is drawn from e(b) and e(V) alone: with ci, spikes on every row, in both
* orientations
seldecomp, makeplot(marginal, ci twoway(name(sdprefixc, replace)))
check_contrib sdprefixc, ci tolerance(1e-12)
seldecomp, makeplot(marginal, swapaxis ci twoway(name(sdprefixv, replace)))
check_contrib sdprefixv, ci vertical tolerance(1e-12)
* and the waterfall: the steps from e(b), with ci spikes from e(V) on every bar
seldecomp, makeplot(waterfall, ci twoway(name(sdprefixw, replace)))
check_waterfall sdprefixw, ci tolerance(1e-12)
* and the waterfall in levels: the means and steps from e(b), with ci spikes on every bar
seldecomp, makeplot(waterfalllevels, ci twoway(name(sdprefixl, replace)))
check_levels sdprefixl, ci
* with residual degrees of freedom (svy, jackknife) the intervals are t-based, as in the table
set_dfr 4
seldecomp, makeplot(stages, ci twoway(name(sddf, replace)))
check_plot sddf, ci tolerance(1e-12)
seldecomp, makeplot(marginal, ci twoway(name(sddfc, replace)))
check_contrib sddfc, ci tolerance(1e-12)
seldecomp, makeplot(marginal, ci swapaxis twoway(name(sddfv, replace)))
check_contrib sddfv, ci vertical tolerance(1e-12)
seldecomp, makeplot(waterfall, swapaxis ci twoway(name(sddfw, replace)))
check_waterfall sddfw, ci horizontal tolerance(1e-12)
seldecomp, makeplot(waterfalllevels, swapaxis ci twoway(name(sddfl, replace)))
check_levels sddfl, ci horizontal

* ---- makeplot requires sequentially or individually ------------------------------
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) makeplot" == 198
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) makeplot(, twoway(name(sdplot4, replace)))" == 198
* the messages: makeplot without either, at estimation and at replay;
* fixothermargins alone; sequentially and individually together
matrix b_u = e(b)
err_msg "seldecomp y x1 x2 [pw=wt], by(group) makeplot" "makeplot requires sequentially or individually"
err_msg "seldecomp y x1 x2 [pw=wt], by(group) makeplot(marginal, swapaxis)" "makeplot requires sequentially or individually"
err_msg "seldecomp y x1 x2 [pw=wt], by(group) fixothermargins makeplot" ///
    "fixothermargins requires sequentially or individually"
err_msg "seldecomp y x1 x2 [pw=wt], by(group) sequentially individually makeplot" ///
    "sequentially and individually may not be combined"
assert mreldif(e(b), b_u) == 0
qui seldecomp y x1 x2 [pw=wt], by(group)
err_msg "seldecomp, makeplot" "makeplot requires results estimated with sequentially or individually"
err_msg "seldecomp, makeplot(waterfall, ci)" "makeplot requires results estimated with sequentially or individually"
* these words are no options: refused (198), not taken as synonyms
foreach o in varbyvar cumulative "varbyvar cumulative" isolated nohold holdothers ///
    "sequentially holdothers" "sequentially cumulative" "individually isolated" ///
    "sequentially refminusfocal" {
    rcof "noi seldecomp y x1 x2 [pw=wt], by(group) `o' makeplot" == 198
}
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) sequentially bayesboot(reps(5) varbyvar nodots) makeplot" == 198
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) sequentially bayesboot(reps(5) sequential nodots) makeplot" == 198

* ---- sign case 2: opposite signs, help example 2 --------------------------------
* explained -.224 against a difference of .489 (and unexplained .713)
di as txt _n "---- opposite signs (help example 2)"
sysuse nlsw88, clear
qui seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) ///
    sequentially fixothermargins makeplot(stages, twoway(name(sdopp, replace)))
matrix S = e(bystage)
forvalues s = 1/3 {
    assert sign(S[`s', colnumb(S, "cum_expl")]) == -sign(e(difference))
    assert sign(S[`s', colnumb(S, "cum_unexpl")]) == sign(e(difference))
}
* titles: the variable label of wage and the value labels of married; the legend's keys
* ("Expl.", "Unexpl.", "Total diff (.489)")
check_plot sdopp, subtitle("Single - Married") vtitle("Hourly wage: difference in means")
local k1 `.sdopp.legend.labels[1]'
local k2 `.sdopp.legend.labels[2]'
local k3 `.sdopp.legend.labels[3]'
assert "`k1'" == "Expl." & "`k2'" == "Unexpl." & "`k3'" == "Total diff (.489)"

* ---- swapaxis: the stages run down the y axis ------------------------------------
di as txt _n "---- makeplot(stages, swapaxis)"
seldecomp, makeplot(stages, swapaxis twoway(name(sdoh, replace)))
check_plot sdoh, horizontal subtitle("Single - Married") vtitle("Hourly wage: difference in means")

* ---- marginal: one bar per stage's marginal explained part, and the difference, explained,
* unexplained ----
di as txt _n "---- makeplot(marginal)"
seldecomp, makeplot(marginal, twoway(name(sdc, replace)))
check_contrib sdc, subtitle("Single - Married") vtitle("Hourly wage: difference in means")
assert `"`.sdc.title.text[1]'"' == "Marginal explanatory contribution at each raking stage"
assert `"`.sdc.yaxis1.major.dlg_ticklabel[1]'"' == "Total diff"
assert `"`.sdc.yaxis1.major.dlg_ticklabel[2]'"' == "Expl."
assert `"`.sdc.yaxis1.major.dlg_ticklabel[6]'"' == "Unexpl."
* horizontal is the default of marginal: swapaxis draws it vertical, the rows from the left at
* 1, 2.2, ..., 5.2 and 6.4, the labels above positive and below negative bars, the value axis
* on y with its title, no title on the category axis, no legend; six rows of short names:
* upright, small labels
seldecomp, makeplot(marginal, swapaxis twoway(name(sdc2, replace)))
check_contrib sdc2, vertical subtitle("Single - Married") vtitle("Hourly wage: difference in means")
assert `.sdc2.xaxis1.style.majorstyle.tickangle.val' == 0
assert `"`.sdc2.xaxis1.major.dlg_ticklabel[3]'"' == "demo" & `.sdc2.xaxis1.major.dlg_tickpos[3]' == 3.2
assert `"`.sdc2.xaxis1.major.dlg_ticklabel[6]'"' == "Unexpl." & `.sdc2.xaxis1.major.dlg_tickpos[6]' == 6.4
assert `"`.sdc2.yaxis1.title.text[1]'"' == "Hourly wage: difference in means"
* the marginal plot's title fits the width of the graph (Plot_fit: 56 characters at 5.5 by 4
* inches, the room growing as xsize / min(xsize, ysize)): the first of "... at each raking
* stage" (54 characters), "... by raking stage" (49), "... at each stage" (47) and "... by
* stage" (42) that fits, else the last; xsize() and ysize() of twoway()
local j 0
foreach xy in "xsize(3)" "xsize(4) ysize(4)" "xsize(5) ysize(4)" "xsize(4.7) ysize(4)" "xsize(4.2) ysize(4)" ///
    "ysize(3)" "xsize(6) ysize(5)" {
    local ++j
    local xs -1
    local ys -1
    if regexm("`xy'", "xsize\(([0-9.]+)\)") local xs = regexs(1)
    if regexm("`xy'", "ysize\(([0-9.]+)\)") local ys = regexs(1)
    seldecomp, makeplot(marginal, twoway(`xy' name(sdct`j', replace)))
    check_contrib sdct`j', xsize(`xs') ysize(`ys')
    seldecomp, makeplot(marginal, swapaxis twoway(`xy' name(sdcv`j', replace)))
    check_contrib sdcv`j', vertical xsize(`xs') ysize(`ys')
    local ttl`j' `"`.sdct`j'.title.text[1]'"'
    assert `"`.sdcv`j'.title.text[1]'"' == `"`ttl`j''"'
}
* (the same under every scheme: the sizes are given)
assert `"`ttl1'"' == "Marginal explanatory contribution by stage"
assert `"`ttl2'"' == "Marginal explanatory contribution by stage"
assert `"`ttl3'"' == "Marginal explanatory contribution by raking stage"
assert `"`ttl4'"' == "Marginal explanatory contribution at each stage"
assert `"`ttl5'"' == "Marginal explanatory contribution by stage"
assert `"`ttl6'"' == "Marginal explanatory contribution at each raking stage"
assert `"`ttl7'"' == "Marginal explanatory contribution at each stage"

* ---- confidence intervals: with ci, on every bar ---------------------------------------
* (none are drawn without ci, although e(V) is posted)
di as txt _n "---- confidence intervals with makeplot(..., ci)"
* every stage bootstrapped, at bayesboot(level(90)): spikes on every row
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ///
    bayesboot(reps(20) rseed(1) nodots level(90)) makeplot(marginal, ci twoway(name(sdcb2, replace)))
assert e(level) == 90
check_contrib sdcb2, ci
* sequentially alone (the other blocks free): CI caps and share labels of the stages plot, and
* the marginal plot and the waterfall, with ci; without ci no interval in any plot type
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially ///
    bayesboot(reps(20) rseed(1) nodots) makeplot(stages, ci twoway(name(sdcum, replace)))
assert "`e(stagetype)'" == "sequentially" & "`e(fixothermargins)'" == ""
check_plot sdcum, ci
seldecomp, makeplot(stages, swapaxis ci twoway(name(sdcumh, replace)))
check_plot sdcumh, ci horizontal
seldecomp, makeplot(marginal, ci twoway(name(sdcumc, replace)))
check_contrib sdcumc, ci
seldecomp, makeplot(waterfall, ci twoway(name(sdcumw, replace)))
check_waterfall sdcumw, ci
seldecomp, makeplot(stages, twoway(name(sdcum0, replace)))
check_plot sdcum0
seldecomp, makeplot(marginal, twoway(name(sdcumc0, replace)))
check_contrib sdcumc0
seldecomp, makeplot(waterfall, twoway(name(sdcumw0, replace)))
check_waterfall sdcumw0
seldecomp, makeplot(waterfalllevels, twoway(name(sdcuml0, replace)))
check_levels sdcuml0

* ---- colours: explained in the scheme's first bar colour, unexplained grey ----------
* under the default scheme, s1mono and s2color; the reference colour is that of a plain
* rbar drawn under the same scheme
di as txt _n "---- colours under the default scheme, s1mono and s2color"
local sch0 `c(scheme)'
foreach sch in `sch0' s1mono s2color {
    preserve
    clear
    qui set obs 1
    gen z = 0
    gen v = 1
    gen x = 1
    twoway rbar z v x, scheme(`sch') name(sdref, replace)
    plot_rgb p1 sdref 1
    restore
    seldecomp, makeplot(stages, ci twoway(scheme(`sch') name(sdsch, replace)))
    check_plot sdsch, ci
    local k1 `.sdsch.legend.labels[1]'
    local k2 `.sdsch.legend.labels[2]'
    assert "`k1'" == "Expl." & "`k2'" == "Unexpl."
    plot_rgb pe sdsch `.sdsch.legend.map[1]'
    plot_rgb pu sdsch `.sdsch.legend.map[2]'
    assert "`pe'" == "`p1'"
    assert "`pu'" == "192 192 192"
    assert "`pe'" != "`pu'"
    * the marginal plot: the explained part and the blocks (plot 2), the unexplained part (3),
    * in both orientations
    foreach sw in "" swapaxis {
        seldecomp, makeplot(marginal, ci `sw' twoway(scheme(`sch') name(sdschc, replace)))
        check_contrib sdschc, ci `=cond("`sw'" != "", "vertical", "")'
        plot_rgb pe sdschc 2
        assert "`pe'" == "`p1'"
        plot_rgb pu sdschc 3
        assert "`pu'" == "192 192 192"
    }
    * the waterfall: the blocks and the explained part (plot 1), the unexplained part (2)
    seldecomp, makeplot(waterfall, ci twoway(scheme(`sch') name(sdschw, replace)))
    check_waterfall sdschw, ci
    plot_rgb pe sdschw 1
    assert "`pe'" == "`p1'"
    plot_rgb pu sdschw 2
    assert "`pu'" == "192 192 192"
    * the waterfall in levels: the blocks (plot 1), the unexplained step (2), the means
    * unfilled (the last bar plot)
    seldecomp, makeplot(waterfalllevels, ci twoway(scheme(`sch') name(sdschl, replace)))
    check_levels sdschl, ci
    plot_rgb pe sdschl 1
    assert "`pe'" == "`p1'"
    plot_rgb pu sdschl 2
    assert "`pu'" == "192 192 192"
    assert "`.sdschl.plotregion1.plot3.style.area.shadestyle.color.snm'" == "none"
    di as txt "scheme `sch': explained `p1', unexplained `pu' -- ok"
}

* ---- waterfall: the steps from 0 to the explained part, then to the difference ----
di as txt _n "---- makeplot(waterfall)"
qui seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) ///
    sequentially fixothermargins makeplot(waterfall, twoway(name(sdw1, replace)))
check_waterfall sdw1, subtitle("Single - Married") vtitle("Hourly wage: difference in means")
assert `"`.sdw1.title.text[1]'"' == "Decomposition of the difference by raking stage"
* six bars and short names: not angled; the last three named "Expl.", "Unexpl." and "Total diff"
* (on one line)
assert `.sdw1.xaxis1.style.majorstyle.tickangle.val' == 0
assert `"`.sdw1.xaxis1.major.dlg_ticklabel[4]'"' == "Expl."
assert `"`.sdw1.xaxis1.major.dlg_ticklabel[5]'"' == "Unexpl."
assert `"`.sdw1.xaxis1.major.dlg_ticklabel[6]'"' == "Total diff"
seldecomp, makeplot(waterfall, swapaxis twoway(name(sdw2, replace)))
check_waterfall sdw2, horizontal subtitle("Single - Married") vtitle("Hourly wage: difference in means")
* every stage bootstrapped, at bayesboot(level(90)): with ci, spikes on every bar
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ///
    bayesboot(reps(20) rseed(1) nodots level(90)) makeplot(waterfall, ci twoway(name(sdw5, replace)))
check_waterfall sdw5, ci
seldecomp, makeplot(waterfall, swapaxis ci twoway(name(sdw5h, replace)))
check_waterfall sdw5h, ci horizontal
* outcomes(): the first outcome, said so, or the one outcome() names; ratio()
gen byte hiwage = wage > 10
err_msg "seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins outcomes(hiwage) makeplot(waterfall, twoway(name(sdw8, replace)))" ///
    "(makeplot: outcome wage, the first of the outcomes; makeplot(, outcome(name)) draws another)" 0
check_waterfall sdw8, outcome(wage)
seldecomp, makeplot(waterfall, outcome(hiwage) twoway(name(sdw9, replace)))
check_waterfall sdw9, outcome(hiwage)
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ratio(hours) ///
    makeplot(waterfall, twoway(name(sdw10, replace)))
check_waterfall sdw10
* (54 characters: on two lines, beyond 40)
text_lines sdw10.yaxis1.title "Hourly wage / Usual hours worked:" "difference in ratios"
* one plot type at a time: a second word after the type is refused (the options follow a comma)
err_msg "seldecomp, makeplot(waterfall marginal)" ///
    "makeplot(): put a comma between the plot type and the options, as in makeplot(waterfall, ci)"
* waterfall takes no argument (waterfall(levels) is waterfalllevels): waterfall()
* is not a plot type
err_msg "seldecomp, makeplot(waterfall())" ///
    "makeplot(): type the plot type first (stages, marginal, waterfall, or waterfalllevels), then the options after a comma"
seldecomp, makeplot(waterfall, twoway(name(sdw11, replace)))
check_waterfall sdw11

* ---- individually: the marginal plot alone -------------------------------------------------
* each block raked alone, then the remainder (the explained part less the sum of the blocks):
* the marginal plot, the default of individually, with the remainder's row; the rakings are not
* cumulative stages, so stages, waterfall and waterfalllevels are refused, at replay and before
* the estimation (e() left as it was)
di as txt _n "---- individually: the marginal plot"
qui seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) ///
    individually fixothermargins makeplot(, twoway(name(sdw6, replace)))
assert "`e(stagetype)'" == "individually" & "`e(fixothermargins)'" == "fixothermargins"
check_contrib sdw6
assert `"`.sdw6.yaxis1.major.dlg_ticklabel[6]'"' == "Remainder"
* vertical, seven rows: angled, vsmall labels (check_contrib)
seldecomp, makeplot(marginal, swapaxis twoway(name(sdw6v, replace)))
check_contrib sdw6v, vertical
assert `.sdw6v.xaxis1.style.majorstyle.tickangle.val' == 30
foreach ty in stages waterfall waterfalllevels {
    err_msg "seldecomp, makeplot(`ty')" ///
        "makeplot(`ty') not allowed with individually; use marginal"
    err_msg "seldecomp, makeplot(`ty', swapaxis twoway(name(sdw6x, replace)))" ///
        "makeplot(`ty') not allowed with individually; use marginal"
}
capture graph describe sdw6x
assert _rc
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) individually ///
    bayesboot(reps(20) rseed(2) nodots) makeplot(marginal, ci twoway(name(sdw7, replace)))
check_contrib sdw7, ci
* individually alone leaves the other blocks free
assert "`e(stagetype)'" == "individually" & "`e(fixothermargins)'" == ""
seldecomp, makeplot(marginal, swapaxis ci twoway(name(sdw7v, replace)))
check_contrib sdw7v, ci vertical
* (six rows: upright)
assert `.sdw7v.xaxis1.style.majorstyle.tickangle.val' == 0
seldecomp, makeplot(marginal, twoway(name(sdw7n, replace)))
check_contrib sdw7n
matrix b_u = e(b)
foreach ty in stages waterfall waterfalllevels {
    err_msg "seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) individually makeplot(`ty')" ///
        "makeplot(`ty') not allowed with individually; use marginal"
    err_msg "seldecomp wage (demo: race collgrad), by(married) individually fixothermargins bayesboot(reps(5) rseed(1) nodots) makeplot(`ty', ci)" ///
        "makeplot(`ty') not allowed with individually; use marginal"
    assert mreldif(e(b), b_u) == 0
}

* ---- the default plot type: the waterfall with sequentially, marginal with individually ----
* makeplot without a plot type draws the waterfall after sequentially, at the estimation and at
* replay, also for a chain: the whole chain and a step of it, chainstep(); after individually the
* marginal plot (above). same_graph: the graph of the default type draws exactly what the type
* named draws
di as txt _n "---- the default plot type"
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially ///
    makeplot(, twoway(name(sddef1, replace)))
check_waterfall sddef1, subtitle("Single - Married") vtitle("Hourly wage: difference in means")
seldecomp, makeplot(waterfall, twoway(name(sddef1w, replace)))
same_graph sddef1 sddef1w
seldecomp, makeplot(, swapaxis twoway(name(sddef1h, replace)))
check_waterfall sddef1h, horizontal
seldecomp, makeplot(waterfall, swapaxis twoway(name(sddef1hw, replace)))
same_graph sddef1h sddef1hw
* individually, the bare makeplot: the marginal plot, as makeplot(marginal) draws it
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) individually makeplot
check_contrib seldecomp
seldecomp, makeplot(marginal, twoway(name(sddef2m, replace)))
same_graph seldecomp sddef2m
* a chain of three values of by() (two steps): the whole chain, at the estimation and at
* replay, the bare makeplot too
gen byte agegrp = 1 + (age >= 38) + (age >= 41)
qui seldecomp wage (demo: collgrad c_city) (geo: south smsa), by(agegrp) chain sequentially fixothermargins ///
    makeplot(, twoway(name(sdch0, replace)))
assert "`e(chain)'" != "" & e(N_pairs) == 2
assert "`e(stagetype)'" == "sequentially"
seldecomp, makeplot(waterfall, twoway(name(sdch0w, replace)))
same_graph sdch0 sdch0w
text_lines sdch0.title "Decomposition of the difference by raking stage"
assert strpos(`"`.sdch0.subtitle.text[1]'"', "(sum of 2 pairs)")
assert `"`.sdch0.xaxis1.major.dlg_ticklabel[1]'"' == "demo"
assert `"`.sdch0.xaxis1.major.dlg_ticklabel[2]'"' == "geo"
assert `"`.sdch0.xaxis1.major.dlg_ticklabel[3]'"' == "Expl."
assert `"`.sdch0.xaxis1.major.dlg_ticklabel[5]'"' == "Total diff"
seldecomp, makeplot
same_graph seldecomp sdch0w
seldecomp, makeplot(, swapaxis twoway(name(sdch0h, replace)))
seldecomp, makeplot(waterfall, swapaxis twoway(name(sdch0hw, replace)))
same_graph sdch0h sdch0hw
* a step, chainstep(#): at replay, and at the estimation
foreach c in 1 2 {
    seldecomp, makeplot(, chainstep(`c') twoway(name(sdch`c', replace)))
    seldecomp, makeplot(waterfall, chainstep(`c') twoway(name(sdch`c'w, replace)))
    same_graph sdch`c' sdch`c'w
    text_lines sdch`c'.title "Decomposition of the difference by raking stage"
    assert `"`.sdch`c'.xaxis1.major.dlg_ticklabel[5]'"' == "Total diff"
}
assert `"`.sdch1.subtitle.text[1]'"' == "(agegrp = 1) - (agegrp = 2)"
assert `"`.sdch2.subtitle.text[1]'"' == "(agegrp = 2) - (agegrp = 3)"
qui seldecomp wage (demo: collgrad c_city) (geo: south smsa), by(agegrp) chain sequentially fixothermargins ///
    makeplot(, chainstep(2) twoway(name(sdch2e, replace)))
same_graph sdch2e sdch2w
drop agegrp

* ---- waterfalllevels: from the Focal Group mean to the Reference Group mean -----------------
* opposite signs (help example 2): explained -.224, so the steps raise Single's mean from 8.08
* to the raked 8.30 before the unexplained step down to Married's 7.59
di as txt _n "---- makeplot(waterfalllevels)"
qui seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) ///
    sequentially fixothermargins makeplot(waterfalllevels, twoway(name(sdl1, replace)))
check_levels sdl1, vtitle("Hourly wage: mean") g1("Single") ///
    g2("Married")
assert _b[explained] < 0 & _b[difference] > 0
* seven bars and short names: not angled, the raked row on two lines; the unexplained step's row
* "Unexpl."
assert `.sdl1.xaxis1.style.majorstyle.tickangle.val' == 0
assert `"`.sdl1.xaxis1.major.dlg_ticklabel[6]'"' == "Unexpl."
seldecomp, makeplot(waterfalllevels, swapaxis twoway(name(sdl2, replace)))
check_levels sdl2, horizontal vtitle("Hourly wage: mean") ///
    g1("Single") g2("Married")
* every stage bootstrapped, at bayesboot(level(90)): with ci, spikes on every bar (the means and
* the steps)
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ///
    bayesboot(reps(20) rseed(1) nodots level(90)) makeplot(waterfalllevels, ci twoway(name(sdl5, replace)))
check_levels sdl5, ci
seldecomp, makeplot(waterfalllevels, swapaxis ci twoway(name(sdl5h, replace)))
check_levels sdl5h, ci horizontal
* sequentially alone (the other blocks free)
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially ///
    makeplot(waterfalllevels, twoway(name(sdl7c, replace)))
check_levels sdl7c
* outcomes(): the first outcome, or the one outcome() names, each from its own equations
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ///
    outcomes(hiwage) makeplot(waterfalllevels, twoway(name(sdl8, replace)))
check_levels sdl8, outcome(wage) vtitle("Hourly wage: mean")
seldecomp, makeplot(waterfalllevels, outcome(hiwage) twoway(name(sdl9, replace)))
check_levels sdl9, outcome(hiwage) vtitle("hiwage: mean")
* ratio(): the levels are ratios of means
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ratio(hours) ///
    makeplot(waterfalllevels, twoway(name(sdl10, replace)))
check_levels sdl10
assert strpos(`"`.sdl10.yaxis1.title.text[1]'"', "Hourly wage / ") == 1
assert substr(`"`.sdl10.yaxis1.title.text[1]'"', -7, 7) == ": ratio"
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ratio(hours) ///
    bayesboot(reps(20) rseed(3) nodots) makeplot(waterfalllevels, ci twoway(name(sdl11, replace)))
check_levels sdl11, ci
* results of sequentially reverse: the forward order unless makeplot(..., reverse) (below)
qui seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) ///
    sequentially fixothermargins reverse makeplot(waterfalllevels, twoway(name(sdl12, replace)))
check_levels sdl12
* after bootstrap: with only e(b) left, from e(b) and e(V)
qui bootstrap, reps(3) seed(1) nodots: seldecomp wage (demo: race collgrad) (geo: south smsa), ///
    by(married) sequentially fixothermargins
drop_bystage
seldecomp, makeplot(waterfalllevels, ci twoway(name(sdl13, replace)))
check_levels sdl13, ci
seldecomp, makeplot(waterfalllevels, twoway(name(sdl13n, replace)))
check_levels sdl13n

* ---- the plotspec is checked before the estimation (results untouched) and at replay ----
di as txt _n "---- plotspec refused before the estimation"
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins
matrix b_A = e(b)
global SD_CLA `"`e(cmdline)'"'
* the syntax makeplot([plottype] [, options]): the plot type first, one of stages, marginal,
* waterfall and waterfalllevels as typed (no abbreviation, no capitals), the options after a
* comma, twoway's own options inside twoway(). The words stage, contributions, waterfall(levels),
* horizontal, stepci, pair() and twoway options among the others are refused, not taken as
* synonyms. Every spec below is refused before the estimation (another decomposition, which
* would replace e(b) if it ran) and at replay, with its message (spec_refused). Among them:
* misspelt keywords and other words twoway refuses; format() not a numeric display format;
* blocklabels() not in pairs, naming no block or a block twice; outcome() not an outcome; by();
* outcome() naming two, blocklabels() with an empty text
local mty "type the plot type first (stages, marginal, waterfall, or waterfalllevels), then the options after a comma"
foreach w in levels stage contributions contibutions Marginal marg waterfal xyz stepci horizontal {
    spec_refused "`w'" "makeplot(): `w' is not a plot type; `mty'"
}
foreach spec in "waterfall(levels)" "waterfall()" "name(g1, replace)" "by(race)" `"title("x")"' {
    spec_refused `"`spec'"' "makeplot(): `mty'"
}
foreach spec in "marginal waterfall" "waterfall marginal, swapaxis" "waterfall levels" "stages ci" ///
    "marginal swapaxis" "waterfall stepci" "stages horizontal" "waterfalllevels name(g1, replace)" {
    gettoken ty : spec
    spec_refused "`spec'" "makeplot(): put a comma between the plot type and the options, as in makeplot(`ty', ci)"
}
local mop "not allowed; the options are swapaxis, ci, normal, tinterval, percentile, bc, bca, reverse, chainstep(), outcome(), format(), blocklabels(), and twoway()"
foreach spec in "waterfall, levels" "marginal, horizontal" "waterfall, stepci" "stages, pair(1)" ", name(g1, replace)" ///
    "marginal, title(x)" "stages, swapaxs" "waterfalllevels, scheme(s1mono)" "marginal, contributions" ", nodraw" {
    local o = substr("`spec'", strpos("`spec'", ",") + 2, .)
    local o = regexr("`o'", "\(.*$", "")
    spec_refused "`spec'" "makeplot(): option `o' `mop'"
}
* twoway(): by() is refused; an option twoway does not take is refused before the estimation by
* the try (Plot_try), at replay when twoway draws the graph
spec_refused ", twoway(by(race))" "makeplot(): by() is not allowed inside twoway()"
spec_refused "waterfall, twoway(title(x) by(race, total))" "makeplot(): by() is not allowed inside twoway()"
foreach spec in "waterfall, twoway(nosuchopt(1))" "marginal, swapaxis twoway(xyz name(g1, replace))" {
    spec_refused "`spec'" "makeplot(twoway()): twoway returned error 198" ///
        "makeplot: twoway returned error 198; the results are posted; redraw the graph with seldecomp, makeplot()"
}
capture graph describe g1
assert _rc
foreach f in "%9.3s" "%td" "%tg" "abc" "%9.3" {
    spec_refused ", format(`f')" "makeplot(format()): `f' is not a numeric display format"
}
spec_refused "waterfall, format(9.3f)" "makeplot(format()): 9.3f is not a numeric display format"
spec_refused ", blocklabels(demo)" `"makeplot(blocklabels()): type pairs of a block name and its text, as in blocklabels(demo "Demographics")"'
spec_refused `", blocklabels(demo "a" geo)"' "makeplot(blocklabels()): type pairs of a block name and its text"
spec_refused `", blocklabels(nosuch "x")"' "makeplot(blocklabels()): nosuch is not a block; the blocks are demo"
spec_refused `", blocklabels(demo "a" demo "b")"' "makeplot(blocklabels()): block demo is given twice"
spec_refused `", blocklabels(demo "")"' `"makeplot(blocklabels()): the text of demo is empty; for a blank label, type " ""'
spec_refused `"waterfall, swapaxis blocklabels(demo "")"' "makeplot(blocklabels()): the text of demo is empty"
spec_refused ", outcome(hours)" "makeplot(outcome()): hours is not the outcome (wage)"
spec_refused ", outcome(wage wage)" "makeplot(outcome()): wage wage names 2 outcomes; name one"
spec_refused "waterfall, outcome(wage hours)" "makeplot(outcome()): wage hours names 2 outcomes"
* chainstep(#): a positive integer, for the results of a chain
spec_refused ", chainstep(1)" "makeplot(chainstep()) requires chain" "makeplot(chainstep()) requires results from a chain run"
foreach c in 0 -1 1.5 x {
    spec_refused "marginal, chainstep(`c')" "makeplot(chainstep()) must be a positive integer"
}
* reverse: results of sequentially reverse; ci: results with confidence intervals (bayesboot() at
* the estimation, or a replication prefix before the replay)
foreach ty in stages marginal waterfall waterfalllevels {
    spec_refused "`ty', reverse" "makeplot(reverse) requires sequentially reverse" ///
        "makeplot(reverse) requires results estimated with sequentially reverse"
    spec_refused "`ty', swapaxis ci" "makeplot(ci) requires bayesboot()" ///
        "makeplot(ci): these results have no confidence intervals"
}
macro drop SD_CLA
* the words of a plot type or an option inside twoway() are twoway's
qui seldecomp wage (demo: race collgrad), by(married) sequentially fixothermargins ///
    makeplot(waterfalllevels, twoway(title("wage levels, marginal, swapaxis") name(sdl14, replace)))
assert `"`.sdl14.title.text[1]'"' == "wage levels, marginal, swapaxis"
check_levels sdl14, notitle

* ---- options, intervals, labels and titles --------------------------------------------------
* a misspelt plot type or option is refused before the estimation (above); the twoway options of
* twoway() are tried on a graph that is not drawn and is dropped, without saving() (which would
* write the user's file); stages names the stages plot; outcome() is checked against the outcomes
di as txt _n "---- stages, outcome(), the twoway options tried before the estimation"
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ///
    makeplot(stages, twoway(name(sdu27a, replace)))
check_plot sdu27a, subtitle("Single - Married") vtitle("Hourly wage: difference in means")
seldecomp, makeplot(stages, swapaxis twoway(name(sdu27b, replace)))
check_plot sdu27b, horizontal
* the options may be abbreviated as the help shows them
seldecomp, makeplot(marginal, swap tw(name(sdu27ab, replace)))
check_contrib sdu27ab, vertical
capture graph describe _seldecomp_try
assert _rc
tempfile sav
local savf `"`sav'.gph"'
capture erase `"`savf'"'
rcof `"noi seldecomp wage (demo: race collgrad), by(married) sequentially fixothermargins makeplot(, twoway(saving("`savf'") xyz))"' == 198
confirm new file `"`savf'"'
* saving() without replace writes the file once: the try did not write it
qui seldecomp wage (demo: race collgrad), by(married) sequentially fixothermargins ///
    makeplot(, twoway(saving("`savf'") name(sdu27s, replace)))
confirm file `"`savf'"'
capture erase `"`savf'"'
* the user's nodraw (the batch idiom with saving()) and name() are not given to the try twice
qui seldecomp wage (demo: race collgrad), by(married) sequentially fixothermargins ///
    makeplot(waterfall, twoway(nodraw saving("`savf'") name(sdu27n, replace)))
confirm file `"`savf'"'
capture erase `"`savf'"'
graph describe sdu27n
* outcome(): the outcome of a run without outcomes() may be named; another variable is refused
* before the estimation, with and without outcomes()
qui seldecomp wage (demo: race collgrad), by(married) sequentially fixothermargins ///
    makeplot(stages, outcome(wage) twoway(name(sdu27o, replace)))
check_plot sdu27o
matrix b_u = e(b)
rcof "noi seldecomp wage (demo: race collgrad), by(married) sequentially fixothermargins makeplot(, outcome(hours))" == 198
rcof "noi seldecomp wage (demo: race collgrad), by(married) sequentially fixothermargins outcomes(hiwage) makeplot(, outcome(hours))" == 198
assert mreldif(e(b), b_u) == 0
qui seldecomp wage (demo: race collgrad), by(married) sequentially fixothermargins outcomes(hiwage) ///
    makeplot(waterfall, outcome(hiwage) twoway(name(sdu27p, replace)))
check_waterfall sdu27p, outcome(hiwage)
* outcome() names one outcome; two names of outcomes (each one of the outcomes) are refused
* (198) before the estimation, with e() unchanged, and at replay
matrix b_u = e(b)
local cl_u `"`e(cmdline)'"'
rcof "noi seldecomp wage (demo: race collgrad), by(married) sequentially fixothermargins outcomes(hours tenure) makeplot(waterfall, outcome(hours tenure))" == 198
rcof "noi seldecomp wage (demo: race collgrad), by(married) sequentially fixothermargins outcomes(hours) makeplot(marginal, outcome(wage hours))" == 198
rcof "noi seldecomp wage (demo: race collgrad), by(married) sequentially fixothermargins outcomes(hours) makeplot(waterfalllevels, outcome(hours hours))" == 198
assert mreldif(e(b), b_u) == 0
assert `"`e(cmdline)'"' == `"`cl_u'"'
qui seldecomp wage (demo: race collgrad), by(married) sequentially fixothermargins outcomes(hours tenure)
matrix b_u = e(b)
foreach spec in ", outcome(hours tenure)" "waterfall, outcome(wage hours)" "stages, outcome(tenure wage)" {
    rcof `"noi seldecomp, makeplot(`spec')"' == 198
    assert mreldif(e(b), b_u) == 0
}
seldecomp, makeplot(waterfall, outcome(tenure) twoway(name(sdu27q, replace)))
check_waterfall sdu27q, outcome(tenure)

* no plot draws an interval unless ci asks for them, and then every bar has its own: the stages'
* caps, the marginal plot's spikes on every row, the waterfall's around the end of every bar (the
* steps, the explained part and the difference), the levels' on the means and the steps
di as txt _n "---- the intervals only with ci, on every bar"
sysuse nlsw88, clear
recode occupation (9 10 = 9 "Farm") (11 12 = 11 "Service/household"), gen(occ)
recode industry (1 2 3 = 2 "Ag/Mining/Constr"), gen(ind)
local u3 wage (educ: collgrad) (job: i.occ i.ind) (union: union) (exp: c.ttl_exp), by(married) contvars(ttl_exp 4)
qui seldecomp `u3' sequentially fixothermargins bayesboot(reps(20) rseed(12) nodots) ///
    makeplot(waterfall, twoway(name(sdu18a, replace)))
confirm matrix e(V)
check_waterfall sdu18a
seldecomp, makeplot(waterfall, ci twoway(name(sdu18b, replace)))
check_waterfall sdu18b, ci
seldecomp, makeplot(waterfall, swapaxis ci twoway(name(sdu18c, replace)))
check_waterfall sdu18c, ci horizontal
seldecomp, makeplot(waterfall, swapaxis twoway(name(sdu18c0, replace)))
check_waterfall sdu18c0, horizontal
seldecomp, makeplot(waterfalllevels, twoway(name(sdu18d, replace)))
check_levels sdu18d
seldecomp, makeplot(waterfalllevels, ci twoway(name(sdu18d1, replace)))
check_levels sdu18d1, ci
seldecomp, makeplot(waterfalllevels, ci swapaxis twoway(name(sdu18e, replace)))
check_levels sdu18e, ci horizontal
seldecomp, makeplot(waterfalllevels, swapaxis twoway(name(sdu18e0, replace)))
check_levels sdu18e0, horizontal
seldecomp, makeplot(marginal, twoway(name(sdu18f, replace)))
check_contrib sdu18f
seldecomp, makeplot(marginal, ci twoway(name(sdu18g, replace)))
check_contrib sdu18g, ci
seldecomp, makeplot(marginal, swapaxis twoway(name(sdu18h, replace)))
check_contrib sdu18h, vertical
seldecomp, makeplot(marginal, swapaxis ci twoway(name(sdu18i, replace)))
check_contrib sdu18i, ci vertical
seldecomp, makeplot(stages, twoway(name(sdu18j, replace)))
check_plot sdu18j
seldecomp, makeplot(stages, ci twoway(name(sdu18k, replace)))
check_plot sdu18k, ci
seldecomp, makeplot(stages, swapaxis twoway(name(sdu18l, replace)))
check_plot sdu18l, horizontal

* the share label (97%) that the dashed total line would strike through is moved beyond the
* line, in both orientations; seven bars of the waterfall: the names angled, "Total diff" too
di as txt _n "---- a label moved off the total line; one angle for all names"
qui seldecomp `u3' sequentially fixothermargins makeplot(stages, twoway(name(sdu31c, replace)))
assert abs(100 * _b[cum_unexpl:educ] / e(difference) - 97) < .5
check_plot sdu31c
assert `check_moved' == 1
seldecomp, makeplot(stages, swapaxis twoway(name(sdu31h, replace)))
check_plot sdu31h, horizontal
assert `check_moved' == 1
seldecomp, makeplot(waterfall, twoway(name(sdu31d, replace)))
check_waterfall sdu31d
assert `.sdu31d.xaxis1.style.majorstyle.tickangle.val' == 30
assert `"`.sdu31d.xaxis1.major.dlg_ticklabel[7]'"' == "Total diff"
* the title of the waterfall in levels names the groups, and there is no subtitle
seldecomp, makeplot(waterfalllevels, twoway(name(sdu31a, replace)))
check_levels sdu31a, g1("Single") g2("Married") vtitle("Hourly wage: mean")
assert `"`.sdu31a.title.text[1]'"' == "From the Single mean to the Married mean, stage by stage"
* group labels of survey length: the title is on one line up to 56 characters (", stage by
* stage" only when the whole fits), else on two lines; with the labels of this data, one line
* (above); individually draws the marginal plot
qui seldecomp `u3' individually fixothermargins makeplot(, twoway(name(sdu31k, replace)))
check_contrib sdu31k
preserve
keep if race < 3
label define rl31 1 "Non-Hispanic White" 2 "Non-Hispanic Black", replace
label values race rl31
qui seldecomp wage (edu: collgrad) (geo: south smsa) (job: union), by(race) sequentially fixothermargins ///
    makeplot(waterfalllevels, twoway(name(sdu31l, replace)))
check_levels sdu31l, g1("Non-Hispanic White") g2("Non-Hispanic Black")
assert `.sdu31l.title.text.arrnels' == 2
assert `"`.sdu31l.title.text[1]'"' == "From the Non-Hispanic White mean"
assert `"`.sdu31l.title.text[2]'"' == "to the Non-Hispanic Black mean"
seldecomp, makeplot(waterfalllevels, swapaxis twoway(name(sdu31m, replace)))
check_levels sdu31m, horizontal g1("Non-Hispanic White") g2("Non-Hispanic Black")
* 35 characters each: the title in two lines, the subtitle of the other plot types in two (73
* characters)
label define rl31 1 "White women in the 1988 NLSW sample" 2 "Black women in the 1988 NLSW sample", replace
qui seldecomp wage (edu: collgrad) (geo: south smsa) (job: union), by(race) sequentially fixothermargins ///
    makeplot(waterfalllevels, twoway(name(sdu31n, replace)))
check_levels sdu31n, g1("White women in the 1988 NLSW sample") g2("Black women in the 1988 NLSW sample")
assert `"`.sdu31n.title.text[1]'"' == "From the White women in the 1988 NLSW sample mean"
assert `"`.sdu31n.title.text[2]'"' == "to the Black women in the 1988 NLSW sample mean"
seldecomp, makeplot(waterfall, twoway(name(sdu31o, replace)))
check_waterfall sdu31o
seldecomp, makeplot(marginal, twoway(name(sdu31p, replace)))
check_contrib sdu31p
text_lines sdu31p.subtitle "White women in the 1988 NLSW sample -" "Black women in the 1988 NLSW sample"
* 43 characters: one line without ", stage by stage"
label define rl31 1 "White" 2 "Black women", replace
qui seldecomp wage (edu: collgrad) (geo: south smsa) (job: union), by(race) sequentially fixothermargins ///
    makeplot(waterfalllevels, twoway(name(sdu31q, replace)))
check_levels sdu31q, g1("White") g2("Black women")
assert `.sdu31q.title.text.arrnels' == 1
assert `"`.sdu31q.title.text[1]'"' == "From the White mean to the Black women mean"
restore

* format() for the values labelled, in every plot type; the shares stay whole percentages
di as txt _n "---- format()"
qui seldecomp `u3' sequentially fixothermargins bayesboot(reps(20) rseed(12) nodots) ///
    makeplot(stages, ci format(%9.3f) twoway(name(sdu30a, replace)))
* (the stages plot labels shares; format() sets the total's value in the legend, which
* check_plot reads)
check_plot sdu30a, ci format(%9.3f)
seldecomp, makeplot(marginal, ci format(%9.3f) twoway(name(sdu30b, replace)))
check_contrib sdu30b, ci format(%9.3f)
seldecomp, makeplot(marginal, swapaxis format(%9.2f) twoway(name(sdu30f, replace)))
check_contrib sdu30f, vertical format(%9.2f)
seldecomp, makeplot(waterfall, ci format(%9.2f) twoway(name(sdu30c, replace)))
check_waterfall sdu30c, ci format(%9.2f)
seldecomp, makeplot(waterfalllevels, ci format(%9.2f) twoway(name(sdu30d, replace)))
check_levels sdu30d, ci format(%9.2f)
seldecomp, makeplot(waterfalllevels, swapaxis ci format(%10.1e) twoway(name(sdu30e, replace)))
check_levels sdu30e, ci horizontal format(%10.1e)

* in levels a step and its share can have opposite signs (south +.0238 (-5%) on the way down
* from Single to Married): the shares are the stage table's, so on the way down a fall has a
* positive share
di as txt _n "---- the sign of the shares in levels"
qui seldecomp wage collgrad south union, by(married) sequentially fixothermargins ///
    makeplot(waterfalllevels, twoway(name(sdu29a, replace)))
check_levels sdu29a, g1("Single") g2("Married") save(L_29)
assert strpos(`"`L_29_labs'"', "|+.0238 (-5%)|")
* Married raked toward Single: the way up, a rise has a positive share
qui seldecomp wage collgrad south union, by(married) focalgroup(1) sequentially fixothermargins ///
    makeplot(waterfalllevels, twoway(name(sdu29c, replace)))
check_levels sdu29c, g1("Married") g2("Single")

* blocklabels() names the blocks' rows in every plot type, at the estimation and at replay; the
* angle rules count the text; the remainder stays "Remainder"; $, `, SMCL braces and quotes are
* shown as typed
di as txt _n "---- blocklabels()"
qui seldecomp `u3' sequentially fixothermargins makeplot(stages, blocklabels(job "Occupation & industry" exp "Experience") ///
    twoway(name(sdu32a, replace)))
check_plot sdu32a, blabels(educ "Occupation & industry" union Experience)
assert `.sdu32a.xaxis1.style.majorstyle.tickangle.val' == 30
seldecomp, makeplot(stages, swapaxis blocklabels(job "Occupation & industry") twoway(name(sdu32b, replace)))
check_plot sdu32b, horizontal blabels(educ "Occupation & industry" union exp)
seldecomp, makeplot(marginal, blocklabels(exp "Experience" educ "College graduate") twoway(name(sdu32c, replace)))
check_contrib sdu32c, blabels("College graduate" job union Experience)
seldecomp, makeplot(marginal, swapaxis blocklabels(exp "Experience" educ "College graduate") ///
    twoway(name(sdu32c2, replace)))
check_contrib sdu32c2, vertical blabels("College graduate" job union Experience)
seldecomp, makeplot(waterfall, blocklabels(union "Union") twoway(name(sdu32d, replace)))
check_waterfall sdu32d, blabels(educ job Union exp)
seldecomp, makeplot(waterfalllevels, blocklabels(job "Jobs") twoway(name(sdu32e, replace)))
check_levels sdu32e, blabels(educ Jobs union exp)
* a short text keeps a two-stage graph upright; a long one angles it (the marginal plot drawn
* vertical too)
qui seldecomp wage (education_level: collgrad) union, by(married) sequentially fixothermargins ///
    makeplot(stages, blocklabels(education_level "Education") twoway(name(sdu32f, replace)))
check_plot sdu32f, blabels(Education union)
assert `.sdu32f.xaxis1.style.majorstyle.tickangle.val' == 0
seldecomp, makeplot(marginal, swapaxis twoway(name(sdu32f2, replace)))
check_contrib sdu32f2, vertical
assert `.sdu32f2.xaxis1.style.majorstyle.tickangle.val' == 30
seldecomp, makeplot(marginal, swapaxis blocklabels(education_level "Education") twoway(name(sdu32f3, replace)))
check_contrib sdu32f3, vertical blabels(Education union)
assert `.sdu32f3.xaxis1.style.majorstyle.tickangle.val' == 0
* individually fixothermargins: the remainder's row keeps its name
qui seldecomp `u3' individually fixothermargins ///
    makeplot(marginal, swapaxis blocklabels(job "Jobs") twoway(name(sdu32g, replace)))
check_contrib sdu32g, vertical blabels(educ Jobs union exp)
seldecomp, makeplot(, blocklabels(job "Jobs") twoway(name(sdu32h, replace)))
check_contrib sdu32h, blabels(educ Jobs union exp)
rcof `"noi seldecomp, makeplot(, blocklabels(_remainder "x"))"' == 198
* $, a backquote, SMCL braces and double quotes, at the estimation and at replay (written with
* char() and passed with macval() so that this do-file's own expansion cannot touch them)
mata: st_local("t1", "Cost " + char(36) + "US {bf} " + char(96) + "q' " + char(34) + "x" + char(34))
qui seldecomp `u3' sequentially fixothermargins ///
    makeplot(waterfall, blocklabels(educ `"`macval(t1)'"') twoway(name(sdu32s, replace)))
check_waterfall sdu32s, blabels(`"Cost {c S|}US {c -(}bf{c )-} {c 96}q' {c 34}x{c 34}"' job union exp)
seldecomp, makeplot(stages, blocklabels(union `"`macval(t1)'"') twoway(name(sdu32t, replace)))
check_plot sdu32t, blabels(educ job `"Cost {c S|}US {c -(}bf{c )-} {c 96}q' {c 34}x{c 34}"' exp)
seldecomp, makeplot(waterfalllevels, swapaxis blocklabels(exp `"`macval(t1)'"') twoway(name(sdu32u, replace)))
check_levels sdu32u, horizontal blabels(educ job union `"Cost {c S|}US {c -(}bf{c )-} {c 96}q' {c 34}x{c 34}"')
seldecomp, makeplot(marginal, swapaxis blocklabels(job `"`macval(t1)'"') twoway(name(sdu32v, replace)))
check_contrib sdu32v, vertical blabels(educ `"Cost {c S|}US {c -(}bf{c )-} {c 96}q' {c 34}x{c 34}"' union exp)
* a text that starts with a double quote is shown as typed too (twoway reads "Degree" holders
* as two lines, and "a" "b" as two lines a and b): every quote is written {c 34}; the angle
* rules count the text as typed
qui seldecomp `u3' sequentially fixothermargins ///
    makeplot(waterfall, blocklabels(educ `""Degree" holders"') twoway(name(sdu32q1, replace)))
check_waterfall sdu32q1, blabels(`"{c 34}Degree{c 34} holders"' job union exp)
seldecomp, makeplot(marginal, blocklabels(union `""a" "b""') twoway(name(sdu32q2, replace)))
check_contrib sdu32q2, blabels(educ job `"{c 34}a{c 34} {c 34}b{c 34}"' exp)
seldecomp, makeplot(stages, blocklabels(educ `""a" "b""' job `"  "Jobs" (all)"') twoway(name(sdu32q3, replace)))
check_plot sdu32q3, blabels(`"{c 34}a{c 34} {c 34}b{c 34}"' `"  {c 34}Jobs{c 34} (all)"' union exp)
seldecomp, makeplot(waterfalllevels, blocklabels(exp `""x""') twoway(name(sdu32q4, replace)))
check_levels sdu32q4, blabels(educ job union `"{c 34}x{c 34}"')
* quotes later in the text are also protected from twoway's option parser
seldecomp, makeplot(waterfall, blocklabels(educ `"a "b" c"') twoway(name(sdu32q5, replace)))
check_waterfall sdu32q5, blabels(`"a {c 34}b{c 34} c"' job union exp)
* an empty text is refused (twoway would show the row's axis position, "2" or "4.2", instead),
* before the estimation and at replay; " " leaves the row without a name
matrix b_u = e(b)
rcof `"noi seldecomp `u3' sequentially fixothermargins makeplot(waterfall, blocklabels(job ""))"' == 198
assert mreldif(e(b), b_u) == 0
foreach vw in stages marginal waterfall waterfalllevels {
    rcof `"noi seldecomp, makeplot(`vw', blocklabels(educ "Education" job ""))"' == 198
}
seldecomp, makeplot(waterfall, blocklabels(job " ") twoway(name(sdu32b1, replace)))
check_waterfall sdu32b1, blabels(educ " " union exp)
seldecomp, makeplot(marginal, blocklabels(job " ") twoway(name(sdu32b2, replace)))
check_contrib sdu32b2, blabels(educ " " union exp)
seldecomp, makeplot(stages, blocklabels(job " ") twoway(name(sdu32b3, replace)))
check_plot sdu32b3, blabels(educ " " union exp)
sysuse nlsw88, clear

* ---- SMCL braces in the labels; long group labels -----------------------------------------------
* (1) The group labels, the variable labels of the outcome and of ratio(), and a string by()
* value go through the escaping of blocklabels() (seldecomp_ptxt), so that twoway reads no SMCL
* tags in them (a group label "Single {bf}x" would turn the subtitle and the axis title bold from
* x on): the stored texts in every plot type (the helpers' check_texts, and literally below) and,
* in the exported SVG, the literal text and no bold. (2) Long group labels: the subtitle "A - B"
* on two lines beyond 72 characters; the levels view's rows of the groups cut with "..." when
* angled (18, 15, 13 or 12 characters by the number of rows), on two lines beyond 20 characters
* in the horizontal view; groups cut to 42 characters in the two-line levels title, to 70 in the
* two-line subtitle; the value axis's title on two lines beyond 40 characters
capture program drop svg_check
program define svg_check
    * svg_check G ["text"]: graph G exported as SVG has an element that reads text, literally,
    * and no text in bold; without a text, the Mata string SVGTXT (a text with $ or ` is built
    * there with char(), where this do-file's own expansion cannot touch it)
    args G txt
    if `"`macval(txt)'"' == "" mata: st_local("txt", SVGTXT)
    tempfile f
    qui graph export `"`f'"', name(`G') as(svg) replace
    mata: st_local("hit", strofreal(strpos(invtokens(cat(st_local("f"))'), ">" + st_local("txt") + "<") > 0))
    mata: st_local("bold", strofreal(strpos(invtokens(cat(st_local("f"))'), "font-weight:bold") > 0))
    if !`hit' | `bold' {
        di as err `"svg_check `G': the SVG does not show `txt' as typed (found `hit', bold `bold')"'
        exit 9
    }
end
di as txt _n "---- SMCL braces in the labels; long group labels"
label define p7b 0 "Single {bf}x" 1 "Married", replace
label values married p7b
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ///
    makeplot(stages, twoway(name(sdp7a, replace)))
check_plot sdp7a
text_lines sdp7a.subtitle "Single {c -(}bf{c )-}x - Married"
svg_check sdp7a "Single {bf}x - Married"
seldecomp, makeplot(stages, swapaxis twoway(name(sdp7b, replace)))
check_plot sdp7b, horizontal
svg_check sdp7b "Single {bf}x - Married"
seldecomp, makeplot(marginal, twoway(name(sdp7c, replace)))
check_contrib sdp7c
svg_check sdp7c "Single {bf}x - Married"
seldecomp, makeplot(marginal, swapaxis twoway(name(sdp7c2, replace)))
check_contrib sdp7c2, vertical
svg_check sdp7c2 "Single {bf}x - Married"
seldecomp, makeplot(waterfall, twoway(name(sdp7d, replace)))
check_waterfall sdp7d
svg_check sdp7d "Single {bf}x - Married"
seldecomp, makeplot(waterfalllevels, twoway(name(sdp7e, replace)))
check_levels sdp7e, g1("Single {c -(}bf{c )-}x") g2("Married")
text_lines sdp7e.title "From the Single {c -(}bf{c )-}x mean to the Married mean"
assert `"`.sdp7e.xaxis1.major.dlg_ticklabel[1]'"' == "Single {c -(}bf{c )-}x"
assert `"`.sdp7e.xaxis1.major.dlg_ticklabel[4]'"' == "Single {c -(}bf{c )-}x, raked"
svg_check sdp7e "Single {bf}x, raked"
svg_check sdp7e "From the Single {bf}x mean to the Married mean"
seldecomp, makeplot(waterfalllevels, swapaxis twoway(name(sdp7f, replace)))
check_levels sdp7f, horizontal g1("Single {c -(}bf{c )-}x")
svg_check sdp7f "Single {bf}x"
* the variable labels of the outcome and of ratio(), a label that starts with a brace, SMCL
* that is not a tag, individually fixothermargins (the marginal plot)
label define p7b 0 "{it}Single" 1 "Married {c 36} {sf}", replace
label var wage "Wage {bf}per hour{sf}"
label var hours "Hours {it}worked"
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) individually fixothermargins ratio(hours) ///
    makeplot(, twoway(name(sdp7g, replace)))
check_contrib sdp7g
text_lines sdp7g.subtitle "{c -(}it{c )-}Single - Married {c -(}c 36{c )-} {c -(}sf{c )-}"
text_lines sdp7g.xaxis1.title "Wage {c -(}bf{c )-}per hour{c -(}sf{c )-} / Hours {c -(}it{c )-}worked:" ///
    "difference in ratios"
svg_check sdp7g "Wage {bf}per hour{sf} / Hours {it}worked:"
svg_check sdp7g "{it}Single - Married {c 36} {sf}"
seldecomp, makeplot(marginal, swapaxis twoway(name(sdp7g2, replace)))
check_contrib sdp7g2, vertical
text_lines sdp7g2.yaxis1.title "Wage {c -(}bf{c )-}per hour{c -(}sf{c )-} / Hours {c -(}it{c )-}worked:" ///
    "difference in ratios"
* the waterfalls, of sequentially fixothermargins
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ratio(hours) ///
    makeplot(waterfalllevels, twoway(name(sdp7h, replace)))
check_levels sdp7h
svg_check sdp7h "Wage {bf}per hour{sf} / Hours {it}worked:"
text_lines sdp7h.title "From the {c -(}it{c )-}Single ratio" "to the Married {c -(}c 36{c )-} {c -(}sf{c )-} ratio"
seldecomp, makeplot(waterfall, swapaxis twoway(name(sdp7i, replace)))
check_waterfall sdp7i, horizontal
* a string by() variable
sysuse nlsw88, clear
gen str20 grp = cond(married, "Married {ul}", "Single")
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(grp) sequentially fixothermargins ///
    makeplot(stages, twoway(name(sdp7j, replace)))
check_plot sdp7j
text_lines sdp7j.subtitle "Married {c -(}ul{c )-} - Single"
seldecomp, makeplot(waterfalllevels, twoway(name(sdp7k, replace)))
check_levels sdp7k, g1("Married {c -(}ul{c )-}") g2("Single")
svg_check sdp7k "Married {ul}, raked"
* a group label that starts with a double quote: the levels view's row of the Focal Group mean
* shows it on one line, as in the other texts
sysuse nlsw88, clear
mata: st_vlmodify("p7q", (0 \ 1), (char(34) + "Single" + char(34) + " women" \ "Married"))
label values married p7q
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins makeplot(waterfalllevels, twoway(name(sdp7l, replace)))
check_levels sdp7l
assert `"`.sdp7l.xaxis1.major.dlg_ticklabel[1]'"' == `"{c 34}Single{c 34} women"'
svg_check sdp7l `""Single" women"'
seldecomp, makeplot(waterfalllevels, swapaxis twoway(name(sdp7m, replace)))
check_levels sdp7m, horizontal
svg_check sdp7m `""Single" women, raked"'
seldecomp, makeplot(stages, twoway(name(sdp7n, replace)))
check_plot sdp7n
text_lines sdp7n.subtitle `"{c 34}Single{c 34} women - Married"'

* (2) group labels of 36 and 42 characters: the subtitle on two lines, the
* rows of the levels view cut when angled and on two lines in the horizontal view
sysuse nlsw88, clear
label define p7l 0 "Unmarried and never married, in 1988" 1 "Married, spouse present in household, 1988", replace
label values married p7l
qui seldecomp wage (demo: race collgrad) (geo: south smsa) (job: union), by(married) sequentially fixothermargins ///
    makeplot(stages, twoway(name(sdp7p, replace)))
check_plot sdp7p
text_lines sdp7p.subtitle "Unmarried and never married, in 1988 -" "Married, spouse present in household, 1988"
seldecomp, makeplot(stages, swapaxis twoway(name(sdp7q, replace)))
check_plot sdp7q, horizontal
seldecomp, makeplot(marginal, twoway(name(sdp7r, replace)))
check_contrib sdp7r
text_lines sdp7r.subtitle "Unmarried and never married, in 1988 -" "Married, spouse present in household, 1988"
seldecomp, makeplot(waterfall, twoway(name(sdp7s, replace)))
check_waterfall sdp7s
text_lines sdp7s.subtitle "Unmarried and never married, in 1988 -" "Married, spouse present in household, 1988"
* levels, angled (seven rows: 15 characters)
seldecomp, makeplot(waterfalllevels, twoway(name(sdp7t, replace)))
check_levels sdp7t
assert `.sdp7t.xaxis1.style.majorstyle.tickangle.val' == 30
assert `"`.sdp7t.xaxis1.major.dlg_ticklabel[1]'"' == "Unmarried..."
assert `"`.sdp7t.xaxis1.major.dlg_ticklabel[5]'"' == "Unmarried..., raked"
assert `"`.sdp7t.xaxis1.major.dlg_ticklabel[7]'"' == "Married..."
text_lines sdp7t.title "From the Unmarried and never married, in 1988 mean" ///
    "to the Married, spouse present in household, 1988 mean"
* horizontal: two lines of at most 20 characters, after a comma when that fits
seldecomp, makeplot(waterfalllevels, swapaxis twoway(name(sdp7u, replace)))
check_levels sdp7u, horizontal
assert `"`.sdp7u.yaxis1.major.dlg_ticklabel[1]'"' == `"`"Unmarried and never"' `"married, in 1988"'"'
assert `"`.sdp7u.yaxis1.major.dlg_ticklabel[5]'"' == `"`"Unmarried and never"' `"married, in 1988, raked"'"'
assert `"`.sdp7u.yaxis1.major.dlg_ticklabel[7]'"' == `"`"Married, spouse"' `"present in..."'"'
* the value axis is extended by a tick where a label would run past its end (here the label
* "-.482 (110%)" of the unexplained step, left of the step's end at 7.412, past the tick 7.2;
* check_levels asserts the fit of every label): the axis starts at 7, where the means' bars start
assert reldif(`.sdp7u.plotregion1.xscale.curmin', 7) < 1e-9
* beside long names the plot is narrower, and the value axis leaves the labels of the bars at
* least as much room as beside short names (the axis is a tick wider than the pad asks for on
* either side when a label needs it, with short names too)
local cmx_l = `.sdp7u.plotregion1.xscale.curmax'
local cmn_l = `.sdp7u.plotregion1.xscale.curmin'
label define p7s 0 "Single" 1 "Married", replace
label values married p7s
qui seldecomp wage (demo: race collgrad) (geo: south smsa) (job: union), by(married) sequentially fixothermargins ///
    makeplot(waterfalllevels, swapaxis twoway(name(sdp7us, replace)))
check_levels sdp7us, horizontal
di as txt "axis range: long names " `cmn_l' " to " `cmx_l' ", short " ///
    `.sdp7us.plotregion1.xscale.curmin' " to " `.sdp7us.plotregion1.xscale.curmax'
assert `cmx_l' - `cmn_l' >= `.sdp7us.plotregion1.xscale.curmax' - `.sdp7us.plotregion1.xscale.curmin'
* 18-character names on one line: the axis from 7 to 8.4
label define p7n 0 "Non-Hispanic White" 1 "Non-Hispanic Black", replace
label values married p7n
qui seldecomp wage (demo: race collgrad) (geo: south smsa) (job: union), by(married) sequentially fixothermargins ///
    makeplot(waterfalllevels, swapaxis twoway(name(sdp7ue, replace)))
check_levels sdp7ue, horizontal g1("Non-Hispanic White") g2("Non-Hispanic Black")
assert `"`.sdp7ue.yaxis1.major.dlg_ticklabel[5]'"' == "Non-Hispanic White, raked"
assert reldif(`.sdp7ue.plotregion1.xscale.curmin', 7) < 1e-9
assert reldif(`.sdp7ue.plotregion1.xscale.curmax', 8.4) < 1e-9
label values married p7l
* individually fixothermargins: the marginal plot, the subtitle on two lines
qui seldecomp wage (demo: race collgrad) (geo: south smsa) (job: union), by(married) individually fixothermargins ///
    makeplot(, twoway(name(sdp7v, replace)))
check_contrib sdp7v
text_lines sdp7v.subtitle "Unmarried and never married, in 1988 -" "Married, spouse present in household, 1988"
* 61 characters: cut to 42 in the levels title; a subtitle of 73 characters in two lines; one
* block (18 characters in the angled rows)
label define p7l 0 "Married women with the spouse present in household, 1988 wave" 1 "Unmarried", replace
qui seldecomp wage (demo: race collgrad south smsa), by(married) sequentially fixothermargins ///
    makeplot(waterfalllevels, twoway(name(sdp7x, replace)))
check_levels sdp7x
text_lines sdp7x.title "From the Married women with the spouse present... mean" "to the Unmarried mean"
assert `"`.sdp7x.xaxis1.major.dlg_ticklabel[1]'"' == "Married women..."
seldecomp, makeplot(stages, twoway(name(sdp7y, replace)))
check_plot sdp7y
text_lines sdp7y.subtitle "Married women with the spouse present in household, 1988 wave -" "Unmarried"
* twelve blocks (sixteen rows): 12 characters
label define p7l 0 "Non-Hispanic White women" 1 "Married", replace
qui seldecomp wage race collgrad south smsa union c_city c.grade c.age c.ttl_exp c.tenure i.industry ///
    c.hours, by(married) sequentially fixothermargins contvars(grade 4 age 4 ttl_exp 4 tenure 4 hours 4) force ///
    makeplot(waterfalllevels, twoway(name(sdp7z, replace)))
check_levels sdp7z
assert `"`.sdp7z.xaxis1.major.dlg_ticklabel[1]'"' == "Non-Hispa..."
* horizontal, sixteen rows: one line cut to 20
seldecomp, makeplot(waterfalllevels, swapaxis twoway(name(sdp7zh, replace)))
check_levels sdp7zh, horizontal
assert `"`.sdp7zh.yaxis1.major.dlg_ticklabel[1]'"' == "Non-Hispanic..."
assert `"`.sdp7zh.yaxis1.major.dlg_ticklabel[14]'"' == "Non-Hispanic..., raked"
* the value axis beyond 40 characters: the labels cut to 40 and the measure on a line of its own
label var wage "Hourly wage in 1988 dollars, from the NLSW 1988 extract"
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ///
    makeplot(stages, twoway(name(sdp7va, replace)))
check_plot sdp7va
text_lines sdp7va.yaxis1.title "Hourly wage in 1988 dollars, from the...:" "difference in means"
seldecomp, makeplot(waterfalllevels, swapaxis twoway(name(sdp7vb, replace)))
check_levels sdp7vb, horizontal
text_lines sdp7vb.xaxis1.title "Hourly wage in 1988 dollars, from the...:" "mean"
seldecomp, makeplot(marginal, twoway(name(sdp7vc, replace)))
check_contrib sdp7vc
text_lines sdp7vc.xaxis1.title "Hourly wage in 1988 dollars, from the...:" "difference in means"

* (3) With ratio(), a closed macro reference (`usual', $US) in the variable label of the outcome
* or of the denominator is kept in the value-axis title (the title is joined in Mata): every plot
* type through check_texts, and the SVG literal ($ and ` written with char())
sysuse nlsw88, clear
mata: st_varlabel("wage", "Wage " + char(96) + "usual" + char(39) + " in " + char(36) + "US")
mata: st_varlabel("hours", "Hours " + char(96) + "usual" + char(39))
mata: SVGTXT = "Wage " + char(96) + "usual" + char(39) + " in " + char(36) + "US / Hours " + char(96) + "usual" + char(39) + ":"
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ratio(hours) ///
    makeplot(stages, twoway(name(sdp7ra, replace)))
check_plot sdp7ra
text_lines sdp7ra.yaxis1.title "Wage {c 96}usual' in {c S|}US / Hours {c 96}usual':" "difference in ratios"
svg_check sdp7ra
seldecomp, makeplot(stages, swapaxis twoway(name(sdp7rb, replace)))
check_plot sdp7rb, horizontal
text_lines sdp7rb.xaxis1.title "Wage {c 96}usual' in {c S|}US / Hours {c 96}usual':" "difference in ratios"
seldecomp, makeplot(marginal, twoway(name(sdp7rc, replace)))
check_contrib sdp7rc
svg_check sdp7rc
seldecomp, makeplot(marginal, swapaxis twoway(name(sdp7rc2, replace)))
check_contrib sdp7rc2, vertical
svg_check sdp7rc2
seldecomp, makeplot(waterfall, twoway(name(sdp7rd, replace)))
check_waterfall sdp7rd
seldecomp, makeplot(waterfalllevels, twoway(name(sdp7re, replace)))
check_levels sdp7re
text_lines sdp7re.yaxis1.title "Wage {c 96}usual' in {c S|}US / Hours {c 96}usual':" "ratio"
svg_check sdp7re
seldecomp, makeplot(waterfalllevels, swapaxis twoway(name(sdp7rf, replace)))
check_levels sdp7rf, horizontal
* without ratio(): the outcome's label alone
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ///
    makeplot(stages, twoway(name(sdp7rg, replace)))
check_plot sdp7rg
text_lines sdp7rg.yaxis1.title "Wage {c 96}usual' in {c S|}US: difference in means"
* in the horizontal levels view a group row is on two lines up to nine rows, and one line cut to
* 20 characters beyond
sysuse nlsw88, clear
label define p7l 0 "Unmarried and never married, in 1988" 1 "Married, spouse present in household, 1988", replace
label values married p7l
* nine rows (five blocks): two lines
qui seldecomp wage race collgrad south smsa union, by(married) sequentially fixothermargins ///
    makeplot(waterfalllevels, swapaxis twoway(name(sdp7rh, replace)))
check_levels sdp7rh, horizontal
assert `"`.sdp7rh.yaxis1.major.dlg_ticklabel[7]'"' == `"`"Unmarried and never"' `"married, in 1988, raked"'"'
assert `"`.sdp7rh.yaxis1.major.dlg_ticklabel[9]'"' == `"`"Married, spouse"' `"present in..."'"'
* ten rows (six blocks): one line
qui seldecomp wage (demo: race collgrad) (geo: south smsa) (job: union) (city: c_city) (ed: c.grade) (exp: c.ttl_exp), ///
    by(married) sequentially fixothermargins contvars(grade 4 ttl_exp 4) ///
    makeplot(waterfalllevels, swapaxis twoway(name(sdp7ri, replace)))
check_levels sdp7ri, horizontal
assert `"`.sdp7ri.yaxis1.major.dlg_ticklabel[1]'"' == "Unmarried and..."
assert `"`.sdp7ri.yaxis1.major.dlg_ticklabel[8]'"' == "Unmarried and..., raked"
assert `"`.sdp7ri.yaxis1.major.dlg_ticklabel[10]'"' == "Married, spouse..."
* twelve rows (eight blocks)
qui seldecomp wage race collgrad south smsa union c_city c.grade c.age, by(married) sequentially fixothermargins ///
    contvars(grade 4 age 4) makeplot(waterfalllevels, swapaxis twoway(name(sdp7rj, replace)))
check_levels sdp7rj, horizontal
assert `"`.sdp7rj.yaxis1.major.dlg_ticklabel[1]'"' == "Unmarried and..."
assert `"`.sdp7rj.yaxis1.major.dlg_ticklabel[12]'"' == "Married, spouse..."
* Beside 50-character blocklabels() texts the value axis is extended for the labels, and beyond
* eight ticks the step follows the wider span (check_levels asserts at most eight ticks and the
* fit of every label). Medium texts of 20 to 27 characters fit at their own step: the axis from
* 7 to 8.4
label define p7s 0 "Single" 1 "Married", replace
label values married p7s
qui seldecomp wage (demo: race collgrad) (geo: south smsa) (job: union), by(married) sequentially fixothermargins
seldecomp, makeplot(waterfalllevels, swapaxis blocklabels(demo "Demographic characteristics (race, college degree)" ///
    geo "Geographic region and city size (South, SMSA)" job "Union membership status of the job") ///
    twoway(name(sdp7rk, replace)))
check_levels sdp7rk, horizontal blabels("Demographic characteristics (race, college degree)" ///
    "Geographic region and city size (South, SMSA)" "Union membership status of the job")
local ntk = `.sdp7rk.xaxis1.major.dlg_tickpos.arrnels'
local stp = `.sdp7rk.xaxis1.major.dlg_tickpos[2]' - `.sdp7rk.xaxis1.major.dlg_tickpos[1]'
di as txt "axis: " `ntk' " ticks of " `stp' " from " `.sdp7rk.plotregion1.xscale.curmin' " to " `.sdp7rk.plotregion1.xscale.curmax'
assert `ntk' <= 8 & `stp' >= .5 - 1e-9
seldecomp, makeplot(waterfalllevels, swapaxis blocklabels(demo "Demographics (race, degree)" ///
    geo "Region and city size" job "Union membership status") twoway(name(sdp7rl, replace)))
check_levels sdp7rl, horizontal blabels("Demographics (race, degree)" "Region and city size" "Union membership status")
assert reldif(`.sdp7rl.plotregion1.xscale.curmin', 7) < 1e-9 & reldif(`.sdp7rl.plotregion1.xscale.curmax', 8.4) < 1e-9
sysuse nlsw88, clear

* ---- reverse (makeplot(..., reverse)): the stages of the reverse order ---------------------
* Results of sequentially reverse hold the stages of both orders. makeplot draws the forward
* order unless reverse asks for the other: in every plot type the blocks from the last to the
* first, from the rev_ equations of e(b) (and of e(V) with ci), ", reverse order" after the groups
* in the subtitle, "Reverse order" as the subtitle of the waterfall in levels (which has none
* otherwise)
di as txt _n "---- reverse: the stages of the reverse order in every plot type"
qui seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) ///
    sequentially fixothermargins reverse makeplot(marginal, twoway(name(sdrv0, replace)))
* the forward order without reverse, and no markers (check_contrib)
check_contrib sdrv0, subtitle("Single - Married")
assert `"`.sdrv0.yaxis1.major.dlg_ticklabel[3]'"' == "demo"
* the orders differ, so the plots differ
assert reldif(_b[rev_marg_expl:demo], _b[marg_expl:demo]) > 1e-3
seldecomp, makeplot(marginal, reverse twoway(name(sdrv1, replace)))
check_contrib sdrv1, reverse subtitle("Single - Married, reverse order")
assert `"`.sdrv1.yaxis1.major.dlg_ticklabel[3]'"' == "age"
assert `"`.sdrv1.yaxis1.major.dlg_ticklabel[5]'"' == "demo"
seldecomp, makeplot(marginal, reverse swapaxis twoway(name(sdrv1v, replace)))
check_contrib sdrv1v, reverse vertical
seldecomp, makeplot(stages, reverse twoway(name(sdrv2, replace)))
check_plot sdrv2, reverse subtitle("Single - Married, reverse order")
assert `"`.sdrv2.xaxis1.major.dlg_ticklabel[1]'"' == "age"
assert `"`.sdrv2.xaxis1.major.dlg_ticklabel[3]'"' == "+ demo"
seldecomp, makeplot(stages, reverse swapaxis twoway(name(sdrv2h, replace)))
check_plot sdrv2h, reverse horizontal
seldecomp, makeplot(waterfall, reverse twoway(name(sdrv3, replace)))
check_waterfall sdrv3, reverse subtitle("Single - Married, reverse order")
assert `"`.sdrv3.xaxis1.major.dlg_ticklabel[1]'"' == "age"
seldecomp, makeplot(waterfall, reverse swapaxis twoway(name(sdrv3h, replace)))
check_waterfall sdrv3h, reverse horizontal
seldecomp, makeplot(waterfalllevels, reverse twoway(name(sdrv4, replace)))
check_levels sdrv4, reverse
text_lines sdrv4.subtitle "Reverse order"
assert `"`.sdrv4.xaxis1.major.dlg_ticklabel[2]'"' == "age"
seldecomp, makeplot(waterfalllevels, reverse swapaxis twoway(name(sdrv4h, replace)))
check_levels sdrv4h, reverse horizontal
* blocklabels() name the blocks in the reverse order's rows
seldecomp, makeplot(marginal, reverse blocklabels(demo "Demographics") twoway(name(sdrv5, replace)))
check_contrib sdrv5, reverse blabels(age geo Demographics)
* without reverse, the forward order in the stages plot and the waterfall
seldecomp, makeplot(stages, twoway(name(sdrv6, replace)))
check_plot sdrv6
seldecomp, makeplot(waterfall, twoway(name(sdrv7, replace)))
check_waterfall sdrv7
* with bayesboot() and ci: the reverse order's intervals
qui seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) ///
    sequentially fixothermargins reverse bayesboot(reps(20) rseed(1) nodots) ///
    makeplot(marginal, reverse ci twoway(name(sdrv8, replace)))
check_contrib sdrv8, reverse ci subtitle("Single - Married, reverse order")
seldecomp, makeplot(stages, reverse ci twoway(name(sdrv9, replace)))
check_plot sdrv9, reverse ci
seldecomp, makeplot(waterfall, reverse ci swapaxis twoway(name(sdrv10, replace)))
check_waterfall sdrv10, reverse ci horizontal
seldecomp, makeplot(waterfalllevels, reverse ci twoway(name(sdrv11, replace)))
check_levels sdrv11, reverse ci
seldecomp, makeplot(marginal, ci twoway(name(sdrv8f, replace)))
check_contrib sdrv8f, ci
* after a replication prefix, from e(b) alone
qui bootstrap, reps(3) seed(1) nodots: seldecomp wage (demo: race collgrad) (geo: south smsa), ///
    by(married) sequentially fixothermargins reverse
drop_bystage
seldecomp, makeplot(marginal, reverse ci twoway(name(sdrv12, replace)))
check_contrib sdrv12, reverse ci tolerance(1e-12)
seldecomp, makeplot(stages, reverse ci twoway(name(sdrv13, replace)))
check_plot sdrv13, reverse ci tolerance(1e-12)
seldecomp, makeplot(stages, reverse twoway(name(sdrv14, replace)))
check_plot sdrv14, reverse tolerance(1e-12)
* outcomes(): the outcome's own equations of the reverse order
gen byte hiwage = wage > 10
qui seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins reverse ///
    outcomes(hiwage)
seldecomp, makeplot(waterfall, reverse outcome(hiwage) twoway(name(sdrv15, replace)))
check_waterfall sdrv15, reverse outcome(hiwage)
seldecomp, makeplot(waterfalllevels, reverse outcome(hiwage) twoway(name(sdrv16, replace)))
check_levels sdrv16, reverse outcome(hiwage)
drop hiwage

* ---- sign case 3: explained larger than the difference --------------------------
* the data: difference -.776, explained -1.312, unexplained .537
di as txt _n "---- explained larger than the difference"
clear
set obs 2000
set seed 11
gen byte g  = _n > 1000
gen byte x1 = runiform() < cond(g, .8, .2)
gen byte x2 = runiform() < .5
gen double y = 2*x1 + 0.4*(g == 0) + 0.3*x2 + rnormal()
qui seldecomp y x1 x2, by(g) sequentially fixothermargins makeplot(stages, twoway(name(sdbig, replace)))
matrix S = e(bystage)
forvalues s = 1/2 {
    assert sign(S[`s', colnumb(S, "cum_expl")]) == sign(e(difference))
    assert abs(S[`s', colnumb(S, "cum_expl")]) > abs(e(difference))
}
check_plot sdbig
seldecomp, makeplot(marginal, swapaxis twoway(name(sdbigv, replace)))
check_contrib sdbigv, vertical

* ---- sign case 4: a difference of exactly 0 --------------------------------------
* y = 20 x1 + 10 x2 + 8 (Focal Group): both group means are 38, while raking moves
* Focal Group's mean to 48 (stage 1) and 46 (stage 2); the share labels are blank
di as txt _n "---- difference of 0"
clear
set obs 400
gen byte g = _n > 200
gen int i = cond(g, _n - 200, _n)
gen byte x1 = cond(g, 1 + (i > 50), 1 + (i > 150))
gen byte x2 = cond(g, cond(x1 == 1, i <= 15, i > 155), mod(i, 2))
gen double y = 20*x1 + 10*x2 + 8*(g == 0)
qui seldecomp y x1 x2, by(g) sequentially fixothermargins makeplot(stages, twoway(name(sdzero, replace)))
assert e(difference) == 0
matrix S = e(bystage)
assert reldif(S[1, colnumb(S, "cum_expl")], -10) < `tol'
assert reldif(S[2, colnumb(S, "cum_expl")], -8) < `tol'
check_plot sdzero
seldecomp, makeplot(marginal, twoway(name(sdzeroc, replace)))
check_contrib sdzeroc
seldecomp, makeplot(marginal, swapaxis twoway(name(sdzerov, replace)))
check_contrib sdzerov, vertical
seldecomp, makeplot(waterfall, twoway(name(sdzerow, replace)))
check_waterfall sdzerow
* in levels: 38, up 10 to 48, down 2 to the raked 46, down 8 to 38; no shares
seldecomp, makeplot(waterfalllevels, twoway(name(sdzerol, replace)))
check_levels sdzerol, g1("g = 0") g2("g = 1")
assert reldif(_b[raked_focal], 46) < `tol' & reldif(_b[focal_group], 38) < `tol'
* negative levels (-138 to -146): the means' bars hang from the top of the axis
gen double yn = -y - 100
qui seldecomp yn x1 x2, by(g) sequentially fixothermargins makeplot(waterfalllevels, twoway(name(sdneg, replace)))
check_levels sdneg
* levels on both sides of 0 (2 to -8): the bars start at 0
gen double yz = -y + 40
qui seldecomp yz x1 x2, by(g) sequentially fixothermargins makeplot(waterfalllevels, swapaxis twoway(name(sdcross, replace)))
check_levels sdcross, horizontal
* means in the tens of thousands (an outcome in dollars): the means read 38,000 and 46,000
gen double yk = y * 1000
qui seldecomp yk x1 x2, by(g) sequentially fixothermargins makeplot(waterfalllevels, twoway(name(sdthou, replace)))
check_levels sdthou, save(L_k)
assert strpos(`"`L_k_labs'"', "|38,000|") & strpos(`"`L_k_labs'"', "|46,000|")

* ---- waterfalllevels: the base, the means' labels, zero steps, axis options, names ----
di as txt _n "---- waterfalllevels: base, labels, axis options, names"
* positive levels (2.49 to 9.35) whose padded axis reaches 0: the ticks start at 0, so the
* means' bars start at 0, although 0 is not among the levels; the mirror
clear
set obs 2000
set seed 3
gen byte g = _n > 1000
gen byte x1 = runiform() < cond(g, .8, .2)
gen byte x2 = runiform() < .5
gen double yd = .3 + 10*x1 + 1*g + .5*x2 + rnormal()
qui seldecomp yd x1 x2, by(g) sequentially fixothermargins makeplot(waterfalllevels, twoway(name(sdl0p, replace)))
assert _b[focal_group] > 2 & _b[focal_group] < 3 & _b[ref_group] > 9
check_levels sdl0p, save(L_0)
assert L_0[1, 1] == 0 & `.sdl0p.plotregion1.yscale.curmin' == 0
gen double ydn = -yd
qui seldecomp ydn x1 x2, by(g) sequentially fixothermargins makeplot(waterfalllevels, swapaxis twoway(name(sdl0n, replace)))
check_levels sdl0n, horizontal save(L_0)
assert L_0[1, 1] == 0 & `.sdl0n.plotregion1.xscale.curmax' == 0
* a step below the rounding threshold reads 0 with a share of 0, in levels as in the waterfall
* (means of a million, steps of about 1e-13 of them)
clear
set obs 400
set seed 99
gen byte g = _n > 200
gen byte x1 = runiform() < cond(g, .8, .2)
gen byte x2 = runiform() < .5
gen double yc = 1e6 + 1e-7*x1 + 1e-7*g
qui seldecomp yc x1 x2, by(g) sequentially fixothermargins makeplot(waterfalllevels, twoway(name(sdlz, replace)))
check_levels sdlz, save(L_z)
assert strpos(`"`L_z_labs'"', "|1,000,000|0 (0%)|0 (0%)|1,000,000|")
seldecomp, makeplot(waterfall, twoway(name(sdlzw, replace)))
check_waterfall sdlzw
* a mean just below 10,000 reads 10,000 (not 1.0e+04), beside 10,006 and 10,009; the mirror
clear
set obs 400
gen byte g = _n > 200
gen int i = cond(g, _n - 200, _n)
gen byte x1 = cond(g, i <= 160, i <= 40)
gen double y = 9997.6 + 10*x1 + 3*g
qui seldecomp y x1, by(g) sequentially fixothermargins makeplot(waterfalllevels, twoway(name(sdl10k, replace)))
check_levels sdl10k, save(L_10k)
assert `"`L_10k_labs'"' == "|10,000|+6 (67%)|10,006|+3 (33%)|10,009"
gen double yn = -y
qui seldecomp yn x1, by(g) sequentially fixothermargins makeplot(waterfalllevels, twoway(name(sdl10kn, replace)))
check_levels sdl10kn, save(L_10k)
assert `"`L_10k_labs'"' == "|-10,000|-6 (67%)|-10,006|-3 (33%)|-10,009"
* the means share one number of decimals, like a table column: 3.8, 4.7 and 3.815 (a gap of
* .015) read 3.800, 4.700 and 3.815, not 3.80 and 4.70 beside 3.815; the mirror
drop y yn
gen double y = 3.5 + 1.5*x1 - .885*g
qui seldecomp y x1, by(g) sequentially fixothermargins makeplot(waterfalllevels, twoway(name(sdlcol, replace)))
check_levels sdlcol, save(L_c)
assert `"`L_c_labs'"' == "|3.800|+.9 (6000%)|4.700|-.885 (-5900%)|3.815"
gen double yn = -y
qui seldecomp yn x1, by(g) sequentially fixothermargins makeplot(waterfalllevels, twoway(name(sdlcoln, replace)))
check_levels sdlcoln, save(L_c)
assert `"`L_c_labs'"' == "|-3.800|-.9 (6000%)|-4.700|+.885 (-5900%)|-3.815"
* ... but a mean gets back only the zeros it dropped, never more decimals than its own 4
* significant digits or the gap's: beside a mean near 0 the others are not padded (250.3 and
* -180.2 beside .0213, not 250.3000 and -180.2000; 123,457 and 200,000 beside .01234, not
* 123,456.80000 and 200,000.00000), nor beyond 10 significant digits (123,457 twice beside
* .00123, not 123,456.70000 and 123,456.70001), while 3.8 beside 3.815 reads 3.800 and 1234.5
* keeps 1234.50; means below 1e-17 read in e-notation (5.00e-18), the others with at most 20
* decimals; each with its mirror
capture program drop lv_case
program define lv_case
    * lv_case name m1 mr m2 "labels" "labels of the mirror": y with the Focal Group mean m1, raked
    * Focal Group mean mr and Reference Group mean m2 (x1 has shares .2 and .8), its levels graph and the
    * mirror's
    args G m1 mr m2 labs labsn
    capture drop y yn
    local b = ((`mr') - (`m1')) / .6
    qui gen double y = (`m1') - .2 * (`b') + (`b') * x1 + ((`m2') - (`mr')) * g
    qui gen double yn = -y
    qui seldecomp y x1, by(g) sequentially fixothermargins makeplot(waterfalllevels, twoway(name(`G', replace)))
    check_levels `G', save(L_c)
    assert `"`L_c_labs'"' == `"`labs'"'
    qui seldecomp yn x1, by(g) sequentially fixothermargins makeplot(waterfalllevels, twoway(name(`G'n, replace)))
    check_levels `G'n, save(L_c)
    assert `"`L_c_labs'"' == `"`labsn'"'
end
lv_case sdlnz 250.3 .0213 -180.2 "|250.3|-250 (58%)|.0213|-180 (42%)|-180.2" ///
    "|-250.3|+250 (58%)|-.0213|+180 (42%)|180.2"
lv_case sdlbs 123456.8 .01234 200000 "|123,457|-123457 (-161%)|.01234|+200000 (261%)|200,000" ///
    "|-123,457|+123457 (-161%)|-.01234|-200000 (261%)|-200,000"
lv_case sdlcap 123456.7 123456.70001 .00123 "|123,457|+.00001 (0%)|123,457|-123457 (100%)|.00123" ///
    "|-123,457|-.00001 (0%)|-123,457|+123457 (100%)|-.00123"
lv_case sdlmix 3.8 3.815 1234.5 "|3.800|+.015 (0%)|3.815|+1231 (100%)|1234.50" ///
    "|-3.800|-.015 (0%)|-3.815|-1231 (100%)|-1234.50"
lv_case sdle 5e-18 1.5e-17 3e-17 ///
    "|5.00e-18|+1.00e-17 (40%)|.000000000000000015|+1.50e-17 (60%)|.000000000000000030" ///
    "|-5.00e-18|-1.00e-17 (40%)|-.000000000000000015|-1.50e-17 (60%)|-.000000000000000030"
* means of about 1,008, 10,000 and 100,008 with steps of tenths (help example 2 shifted): the
* means read with the decimals that tell them apart
sysuse nlsw88, clear
foreach k in 1000 9991.9 100000 {
    gen double wk = wage + `k'
    qui seldecomp wk (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins ///
        makeplot(waterfalllevels, twoway(name(sdlk, replace)))
    check_levels sdlk, save(L_m)
    local labs_`=round(`k')' `"`L_m_labs'"'
    drop wk
}
assert strpos(`"`labs_1000'"', "|1008.1|") & strpos(`"`labs_1000'"', "|1007.6")
assert strpos(`"`labs_9992'"', "|10,000.0|") & strpos(`"`labs_9992'"', "|9,999.5")
assert strpos(`"`labs_100000'"', "|100,008.1|") & strpos(`"`labs_100000'"', "|100,007.6")
* a value-axis option of twoway() that takes the axis to 0: its numbers join the axis, so the
* means' bars start at 0, and the ticks run from 0 over the whole axis; options that do not
* reach 0 keep the base (swapaxis: the x axis's options)
qui seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) sequentially fixothermargins
local j 0
foreach spec in "yscale(range(0))" "ysc(r(0))" "ylabel(0(2)10)" `"ylabel(0 "zero" 5 10, angle(0))"' ///
    "ylabel(, angle(0)) ytick(0)" "swapaxis xscale(range(0))" "swapaxis xlabel(0(2)10)" {
    local ++j
    local sw
    local tw `"`spec'"'
    if substr(`"`spec'"', 1, 9) == "swapaxis " {
        local sw swapaxis
        local tw = substr(`"`spec'"', 10, .)
    }
    seldecomp, makeplot(waterfalllevels, `sw' twoway(`tw' name(sdlax`j', replace)))
    local vx = cond("`sw'" != "", "x", "y")
    check_levels sdlax`j', `=cond("`vx'" == "x", "horizontal", "")' save(L_a)
    assert L_a[1, 1] == 0
    assert `.sdlax`j'.`vx'axis1.major.dlg_tickpos[1]' == 0
    assert `.sdlax`j'.`vx'axis1.major.dlg_tickpos.arrnels' <= 8
}
* (with swapaxis, yscale() is the category axis's)
foreach spec in "yscale(range(5 12))" "ylabel(#10)" "swapaxis yscale(range(0))" {
    local ++j
    local sw
    local tw `"`spec'"'
    if substr(`"`spec'"', 1, 9) == "swapaxis " {
        local sw swapaxis
        local tw = substr(`"`spec'"', 10, .)
    }
    seldecomp, makeplot(waterfalllevels, `sw' twoway(`tw' name(sdlax`j', replace)))
    check_levels sdlax`j', `=cond("`sw'" != "", "horizontal", "")' save(L_a) ///
        `=cond(strpos("`spec'", "#10"), "anyticks", "")'
    assert L_a[1, 1] > 0
}
* group names beside "Unexpl.": longer than 7 characters, the names are angled at seven bars
label define m10 0 "Unmarried1" 1 "Married 30", replace
label values married m10
qui seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) ///
    sequentially fixothermargins makeplot(waterfalllevels, twoway(name(sdlg, replace)))
check_levels sdlg, g1("Unmarried1") g2("Married 30")
assert `.sdlg.xaxis1.style.majorstyle.tickangle.val' == 30

* ---- twelve blocks with 32-character names --------------------------------------------
* the vertical views angle the names; small contributions of the sign opposite to the
* difference round to "0%", never "-0%" (stage 1 here: .0042 against -1.34)
di as txt _n "---- twelve long block names"
clear
set obs 600
set seed 5
gen byte g = _n > 300
gen byte x1 = 1 + (runiform() < cond(g, .8, .2)) + (runiform() < cond(g, .8, .2))
gen byte x2 = runiform() < .5
gen double ys = 2*x1 - 1*g + x2 + rnormal()
forvalues j = 1/12 {
    gen byte `=substr("block_`j'_" + "x"*32, 1, 32)' = runiform() < cond(g, .3 + `j'/30, .5)
}
qui seldecomp ys block_*, by(g) sequentially fixothermargins makeplot(stages, twoway(name(sd12, replace)))
assert e(N_stages) == 12
local nm1 : word 1 of `e(blocks)'
assert 100 * _b[marg_expl:`nm1'] / e(difference) > -.5 & _b[marg_expl:`nm1'] / e(difference) < 0
check_plot sd12
seldecomp, makeplot(stages, swapaxis twoway(name(sd12h, replace)))
check_plot sd12h, horizontal
seldecomp, makeplot(marginal, twoway(name(sd12c, replace)))
check_contrib sd12c
* vertical, fifteen rows: angled, vsmall labels
seldecomp, makeplot(marginal, swapaxis twoway(name(sd12v, replace)))
check_contrib sd12v, vertical
assert `.sd12v.xaxis1.style.majorstyle.tickangle.val' == 30
seldecomp, makeplot(waterfall, twoway(name(sd12w, replace)))
check_waterfall sd12w
assert `.sd12w.xaxis1.style.majorstyle.tickangle.val' == 30
seldecomp, makeplot(waterfall, swapaxis twoway(name(sd12wh, replace)))
check_waterfall sd12wh, horizontal
* two stages, one block name longer than 10 characters: angled as well
qui seldecomp ys (education_level: x1) x2, by(g) sequentially fixothermargins makeplot(stages, twoway(name(sdlong, replace)))
assert e(N_stages) == 2
check_plot sdlong

* ---- a single stage ---------------------------------------------------------------
di as txt _n "---- one stage"
qui seldecomp ys (all: x1 x2), by(g) sequentially fixothermargins makeplot(stages, twoway(name(sdone, replace)))
assert e(N_stages) == 1
check_plot sdone
seldecomp, makeplot(stages, swapaxis twoway(name(sdoneh, replace)))
check_plot sdoneh, horizontal
seldecomp, makeplot(marginal, twoway(name(sdonec, replace)))
check_contrib sdonec
seldecomp, makeplot(marginal, swapaxis twoway(name(sdonev, replace)))
check_contrib sdonev, vertical
seldecomp, makeplot(waterfall, twoway(name(sdonew, replace)))
check_waterfall sdonew
seldecomp, makeplot(waterfalllevels, twoway(name(sdonel, replace)))
check_levels sdonel

* ---- labels with $, a backquote and double quotes ------------------------------------
* twoway expands its text options again, so $ and ` are passed as {c S|} and {c 96}: the
* graph is drawn and shows them ($ and ` are written with char() so that this do-file's
* own macro expansion cannot touch them)
di as txt _n "---- labels with special characters"
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
keep in 1/600
mata: st_vlmodify("gl", (0 \ 1), ("cost " + char(36) + "US" \ "it" + char(96) + "s " + char(34) + "q" + char(34)))
label values group gl
mata: st_varlabel("y", "Spend " + char(36) + "X " + char(96) + "a'")
qui seldecomp y x1 x2, by(group) sequentially fixothermargins makeplot(stages, twoway(name(sdlab, replace)))
check_plot sdlab, subtitle(`"cost {c S|}US - it{c 96}s {c 34}q{c 34}"') vtitle("Spend {c S|}X {c 96}a': difference in means")
seldecomp, makeplot(marginal, twoway(name(sdlabc, replace)))
check_contrib sdlabc, subtitle(`"cost {c S|}US - it{c 96}s {c 34}q{c 34}"') vtitle("Spend {c S|}X {c 96}a': difference in means")
seldecomp, makeplot(marginal, swapaxis twoway(name(sdlabv, replace)))
check_contrib sdlabv, vertical subtitle(`"cost {c S|}US - it{c 96}s {c 34}q{c 34}"') vtitle("Spend {c S|}X {c 96}a': difference in means")
seldecomp, makeplot(waterfalllevels, twoway(name(sdlabl, replace)))
check_levels sdlabl, vtitle("Spend {c S|}X {c 96}a': mean") ///
    g1("cost {c S|}US") g2(`"it{c 96}s {c 34}q{c 34}"')

capture graph drop _all

* a twoway option that twoway rejects: the error is reported, the results stay posted
import delimited using "data/seldecomp_test.csv", clear varnames(1) case(preserve) stringcols(8) asdouble
* makeplot leaves the data in memory as they were (values, order, labels)
sort x4 id
datasignature set, reset
local srt : sortedby
seldecomp y x1 x2 [pw=wt], by(group) sequentially fixothermargins
matrix b_keep = e(b)
rcof "noi seldecomp y x1 x2 [pw=wt], by(group) sequentially fixothermargins makeplot(, twoway(nosuchopt))" == 198
assert "`e(cmd)'" == "seldecomp"
rcof "noi seldecomp, makeplot(, twoway(nosuchopt(1)))" == 198
assert mreldif(e(b), b_keep) < 1e-12
seldecomp, makeplot(stages, twoway(name(after_bad, replace)))
graph describe after_bad
seldecomp, makeplot(stages, swapaxis twoway(name(after_bad2, replace)))
seldecomp, makeplot(marginal, twoway(name(after_bad3, replace)))
seldecomp, makeplot(marginal, swapaxis twoway(name(after_bad4, replace)))
seldecomp, makeplot(waterfall, twoway(name(after_bad5, replace)))
seldecomp, makeplot(waterfalllevels, swapaxis twoway(name(after_bad6, replace)))
rcof "noi seldecomp, makeplot(marginal, twoway(nosuchopt))" == 198
rcof "noi seldecomp, makeplot(marginal, nosuchopt)" == 198
datasignature confirm
assert "`: sortedby'" == "`srt'"
capture label list _seldecomp_stage
assert _rc == 111
capture graph drop _all

di as res _n "PASSED: test_plot.do"

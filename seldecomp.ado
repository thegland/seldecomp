*! version 1.0.0  18nov2026
*! seldecomp: raking (reweighting) decomposition of differences in means between two groups
*! Method: Pylypchuk and Selden (2008), Journal of Health Economics 27: 1109-1128
*! Raking: seldecomp's own Mata raking (iterative proportional fitting), trimming included

* Layout of this file, in order:
*   seldecomp                      the command: replay and its display options, else Estimate
*                                  (then Chain when chain is specified)
*   Estimate                       parses the command, prepares the data, runs Engine and
*                                  Bootstrap, posts e()
*   Chain                          chain: Estimate once as a probe and once per pair of adjacent
*                                  values of by(), one replicate loop over the pairs, e()
*   Parse_term ... Zero_named      parsing, error messages, the two groups of by()
*   Bin_quantiles, Bin_label       contvars(): quantile bins of continuous raking variables
*   Engine, Rake_mata              the decomposition for one weight vector, all stages
*   Covariate_detail               covariatedetail: the shares of each level by group
*   Draw_weights ... Recut         bayesboot(): replicate weights, replicate loop, intervals
*   Display ... Display_covdetail  the output (header, tables, diagnostics, notes) and the
*                                  messages after a replication prefix (Prefix_epilog,
*                                  Misreps_fold, Misreps_warn); Display_estat: estat bootstrap,
*                                  called from seldecomp_estat.ado
*   Plot_spec ... Make_plot        makeplot()
*   Mata, to the end of the file   the raking and trimming (seldecomp_rake, seldecomp_trim), the
*                                  targets and the support screen (seldecomp_stage1, _stage2,
*                                  _screen), the pre-flight check of svy replicates
*                                  (seldecomp_prefail), the BCa acceleration, text helpers, and
*                                  the caches that carry information between the calls of a
*                                  replication prefix (seldecomp_failcache, _gencache, _cutcache)

program define seldecomp, eclass sortpreserve byable(recall) properties(svyb svyj)
    version 16
    Reject_nodiagnostics `macval(0)'
    if replay() {
        if _by() error 190
        if "`e(cmd)'" != "seldecomp" error 301
        foreach a in make makep makepl makeplo makeplot {
            local 0 : subinstr local 0 " `a'()" " makeplot", all
            local 0 : subinstr local 0 ",`a'()" ", makeplot", all
        }
        local rspec Level(passthru) DIAGnostics COVariatedetail MAKEplot MAKEplot2(passthru) ///
                  cformat(passthru) pformat(passthru) sformat(passthru)        ///
                  CITable VARLabels
        capture syntax [, `rspec' prefixepilog estatboot(string) ]
        local rc = _rc
        // run by a replication prefix after its replicate loop (e(prefix_epilog), posted by
        // the full-sample call): never an error, which would fail the prefix's estimation
        if !`rc' & "`prefixepilog'" != "" {
            capture noisily Prefix_epilog
            exit
        }
        // estat bootstrap after bayesboot() (seldecomp_estat.ado passes its checked options)
        if !`rc' & "`estatboot'" != "" {
            Display_estat, `estatboot'
            exit
        }
        // an option replay does not take: the message lists those it takes (Replay_error)
        if `rc' {
            Replay_error `"`rspec'"' `"`macval(0)'"'
        }
        // without level(), replay at the level of the estimation, as the bootstrap CI table does
        if "`level'" == "" local level level(`=string(e(level), "%12.0g")')
        // the plotspec is user text (blocklabels(), titles): handled without macro expansion
        if `"`macval(makeplot2)'"' != "" {
            local makeplot makeplot
            mata: st_local("makeplot2", substr(st_local("makeplot2"), strpos(st_local("makeplot2"), ///
                "(") + 1, strlen(st_local("makeplot2")) - strpos(st_local("makeplot2"), "(") - 1))
        }
        if "`makeplot'" != "" {
            // the graph is drawn from e() (rebuilt from e(b) alone should a prefix drop e(bystage))
            if "`e(stagetype)'" == "" {
                di as err "makeplot requires results estimated with sequentially or individually"
                exit 198
            }
            Make_plot `=("`varlabels'" != "")' `level' `macval(makeplot2)'
            exit
        }
        // the display that bootstrap: and jackknife: run after their loop (c(prefix) holds the
        // prefix; a replay typed by the user has none): replicates that failed at the
        // iterate() limit where the point estimate stopped too are counted in its note
        if strtrim(c(prefix)) != "" Misreps_fold
        // there, the display options typed in the command (e(prefix_dispopts), posted by the
        // full-sample call; the prefix passes only its own), unless the prefix passes the same
        if strtrim(c(prefix)) != "" & `"`e(prefix_dispopts)'"' != "" {
            local dopts diagnostics covariatedetail cformat pformat sformat citable varlabels
            foreach o of local dopts {
                local t_`o' `"``o''"'
            }
            local 0 `", `e(prefix_dispopts)'"'
            capture syntax [, DIAGnostics COVariatedetail cformat(passthru) ///
                pformat(passthru) sformat(passthru) CITable VARLabels ]
            foreach o of local dopts {
                if _rc | `"`t_`o''"' != "" local `o' `"`t_`o''"'
            }
        }
        Display, `level' `diagnostics' `covariatedetail' `cformat' ///
            `pformat' `sformat' `citable' `varlabels' replay
        // bayesboot() results: estat bootstrap shows the intervals at the level of the last
        // display; a replay at another level records it (hidden: e(level) and the intervals in
        // e() stay the estimation's), and the replays after it keep the record current
        if "`e(bsnames)'" != "" {
            if r(level) != e(level) | e(level_replay) < . {
                ereturn hidden scalar level_replay = r(level)
            }
        }
        exit
    }
    // by prefix: the by-group index is only reliably visible in this (byable) program, so
    // it is handed to Estimate explicitly, as nl.ado does for nl_7 ([P] byable)
    local byopt
    if _by() local byopt seldbyvar(`_byindex') seldbyidx(`=_byindex()') seldbyvars(`_byvars')
    // macval(): the command line was expanded once already; a literal $name in an if
    // expression or in a group value must not be expanded a second time
    Estimate `macval(0)' `byopt'
    // chain: Estimate only said so; Chain runs the pairs
    if "`chainrun'" != "" {
        Chain `macval(0)' `byopt'
    }
    ereturn local cmdline `"seldecomp `macval(0)'"'
    ereturn local estat_cmd seldecomp_estat
    Display, level(`=string(e(level), "%12.0g")') `dispopts'
    if "`makeplot'" != "" {
        Make_plot `=strpos(" `dispopts' ", " varlabels ") > 0' level(`=string(e(level), "%12.0g")') `macval(makeplot2)'
    }
end


* ============================================================================
* Estimate: parse, prepare data, compute, post results
* ============================================================================
program define Estimate, eclass
    version 16
    // bayesboot() and makeplot() typed with empty parentheses: -syntax- silently drops such an
    // option (Stata 19), so turn them into the bare words before parsing
    foreach a in bayes bayesb bayesbo bayesboo bayesboot {
        local 0 : subinstr local 0 " `a'()" " bayesboot", all
        local 0 : subinstr local 0 ",`a'()" ", bayesboot", all
    }
    foreach a in make makep makepl makeplo makeplot {
        local 0 : subinstr local 0 " `a'()" " makeplot", all
        local 0 : subinstr local 0 ",`a'()" ", makeplot", all
    }
    // the options, in a local: when -syntax- rejects the command, Top_error looks the offending
    // option up (an oaxaca or bootstrap habit gets a message naming seldecomp's equivalent)
    local spec FOCALgroup(string) CONTvars(string asis) SEQuentially INDividually FIXothermargins ///
          INSUPport                                                                   ///
          BAYESboot BAYESboot2(passthru)                                              ///
          ITERate(integer 10000) TOLerance(string) CTRLTOLerance(string)              ///
          trimhirel(string) trimhiabs(string) trimlorel(string) trimloabs(string)     ///
          TRIMFREQuency(string) trimsupportoverride                                   ///
          GENerate(string) COVariatedetail MAKEplot MAKEplot2(passthru)               ///
          MINcell(integer 10) NOIsily DIAGnostics FORCE                             ///
          seldbyvar(varname) seldbyidx(integer 0)                                     ///
          cformat(passthru) pformat(passthru) sformat(passthru) VARLabels ///
          RATio(varname) OUTcomes(varlist) REVerse seldbyvars(string) DETail           ///
          CHAIN chainprobe(name) chainpair(integer 0) chainvars(string)
    // a word followed by a parenthesis before the comma is no raking term: an option typed
    // before the comma, whose absence -syntax- would report ("option by() required")
    Comma_check `"`spec'"' `"`macval(0)'"'
    capture syntax anything(id="varlist" equalok) [if] [in] [aw pw iw fw] , BY(varname) [ `spec' ]
    if _rc {
        Top_error `"`spec'"' `"`macval(0)'"'
    }
    // by(): a variable of the by prefix takes one value in each by-group
    if "`seldbyvars'" != "" {
        if `: list by in seldbyvars' {
            di as err "`by' may not be both in by() and in the by prefix"
            exit 190
        }
    }
    // chain: Chain, called next by seldecomp, runs the pairs, each through
    // Estimate in pair mode (chainpair()) after a probe (chainprobe()); the command as typed
    // only hands over to it
    if "`chain'" != "" & "`chainprobe'" == "" & `chainpair' == 0 {
        c_local chainrun 1
        exit
    }
    // pair mode: the prefix of its messages that do not name the pair's groups
    local chp = cond(`chainpair' > 0, "pair `chainpair': ", "")
    // display formats are checked before any work is done
    Display_fmts, `cformat' `pformat' `sformat' check

    // keep the outer option values safe from the inner -syntax- calls below
    local cmdline0 : copy local 0
    local o_anything   `"`anything'"'
    local o_focalgroup `"`macval(focalgroup)'"'
    // the staged decompositions: sequentially, the blocks added one stage at a time (with
    // fixothermargins, the papers' procedure); individually, every block raked alone from
    // Focal Group's starting point, then all blocks. Either sets the staged flag (the local
    // sequential); individually also the local variant. By default the other blocks (not yet
    // added, or not the one raked alone) are unconstrained; fixothermargins holds them at
    // Focal Group's own shares (the local holdothers, Engine's option word); checked against the
    // options as typed
    if "`sequentially'" != "" & "`individually'" != "" {
        di as err "sequentially and individually may not be combined"
        exit 198
    }
    local variant `individually'
    if "`variant'" != "" {
        foreach o in reverse detail {
            if "``o''" != "" {
                di as err "`o' may not be combined with individually"
                exit 198
            }
        }
    }
    if "`sequentially'`individually'" != "" local sequential sequential
    if "`fixothermargins'" != "" & "`sequential'" == "" {
        di as err "fixothermargins requires sequentially or individually"
        exit 198
    }
    local holdothers
    if "`fixothermargins'" != "" local holdothers holdothers
    local o_sequential `sequential'
    local o_variant    `variant'
    local o_insupport  `insupport'
    local o_generate   `generate'
    local o_replace
    // the confidence level: bayesboot(level()), else the default (a replication prefix takes
    // its own level() and passes it to the replay that displays its results)
    local o_level = c(level)
    local o_citable
    local o_weight     `weight'
    local o_exp        `"`exp'"'
    local o_if         `"`macval(if)'"'
    local o_in         `"`in'"'
    local o_byvar      `seldbyvar'
    local o_byidx      `seldbyidx'
    local o_force      `force'
    local o_ratio      `ratio'
    local o_outcomes   `outcomes'
    local o_reverse    `reverse'
    local o_detail     `detail'
    // generate(newvar [, replace]): the name of the raked weights; replace overwrites them
    // (parsed without -syntax-, which would clear the qualifiers and weights read below)
    if `"`o_generate'"' != "" {
        gettoken gname grest : o_generate, parse(",")
        mata: st_local("gname", strtrim(st_local("gname")))
        mata: st_local("gopt", strtrim(substr(strtrim(st_local("grest")), 2, .)))
        mata: st_local("gcom", substr(strtrim(st_local("grest")), 1, 1))
        if `: word count `gname'' != 1 | (`"`gcom'"' != "" & `"`gcom'"' != ",") | ///
            !inlist(`"`gopt'"', "", "replace") {
            di as err "generate(): type generate(newvar) or generate(newvar, replace)"
            exit 198
        }
        local o_generate `gname'
        local o_replace `gopt'
        local generate `gname'
    }
    // chain (the probe): what a chain of pairs does not take
    if "`chainprobe'" != "" {
        if `"`macval(focalgroup)'"' != "" {
            di as err "focalgroup() may not be combined with chain"
            exit 198
        }
        // detail: its term stages are the pair's alone
        foreach o in individually reverse detail outcomes generate covariatedetail {
            if "``o''" == "" continue
            local oo = cond(inlist("`o'", "outcomes", "generate"), "`o'()", "`o'")
            di as err "`oo' may not be combined with chain"
            exit 198
        }
    }

    // ---- basic option checks ------------------------------------------------
    if "`weight'" == "fweight" {
        // the remedy: expand the data; frequency weights passed as pweights or iweights give
        // the same point estimates, but bayesboot() draws one weight per record, not per unit
        mata: st_local("fwexp", strtrim(ustrregexra(st_local("exp"), "^\s*=", "")))
        di as err "fweights not allowed"
        di as err `"expand the data first: expand `macval(fwexp)'"'
        exit 101
    }
    // bayesboot() and makeplot() may be given with empty parentheses: strip the wrapper
    // (in Mata: the plotspec is user text, not to be expanded again)
    foreach o in bayesboot makeplot {
        if `"`macval(`o'2)'"' != "" {
            local `o' `o'
            mata: st_local("`o'2", substr(st_local("`o'2"), strpos(st_local("`o'2"), "(") + 1, ///
                strlen(st_local("`o'2")) - strpos(st_local("`o'2"), "(") - 1))
        }
    }
    // chain, pair mode: the plotspec was checked by the probe, and Chain draws
    if `chainpair' > 0 {
        local makeplot
        local makeplot2
    }
    if "`makeplot'" != "" & "`o_sequential'" == "" {
        di as err "makeplot requires sequentially or individually"
        exit 198
    }
    if "`reverse'" != "" & "`o_sequential'" == "" {
        di as err "reverse requires sequentially"
        exit 198
    }
    // detail: the terms of each block added one at a time
    if "`detail'" != "" & "`o_sequential'" == "" {
        di as err "detail requires sequentially"
        exit 198
    }
    // a replication prefix: Stata's replicate loops (_loop_bs, _loop_rw, _loop_brr,
    // _loop_jknife*) put the word _loop in c(prefix), which the full-sample call and a call
    // without a prefix never have (Stata resets c(prefix) when a program ends)
    local ccpfx = strtrim(stritrim(c(prefix)))
    local ccloop : list posof "_loop" in ccpfx
    // the caller's sort order, read here, before any work (the pre-flight check below reads it)
    local ccsort : sortedby
    local ccrep 0
    foreach p in bootstrap jackknife svy brr sdr bayesboot {
        if `: list posof "`p'" in ccpfx' local ccrep 1
    }
    local ccname
    if `ccrep' {
        foreach p of local ccpfx {
            if "`ccname'" == "" & inlist("`p'", "bootstrap", "jackknife", "svy", "bayesboot") local ccname `p'
        }
        if "`ccname'" == "svy" local ccname = "svy " + word("`ccpfx'", `: list posof "svy" in ccpfx' + 1)
    }
    // failed replicates: a replicate call records what failed
    // it, silently, for the one message after the prefix's loop (Misreps_warn); every other
    // call starts a new record
    if !`ccloop' mata: seldecomp_failcache("reset", "", "")
    // bayesboot() and makeplot under a replication prefix: refused
    // before any work in the prefix's first, full-sample call, whose c(prefix) holds the
    // prefix's name without _loop (bootstrap, jackknife, svy brr, svy jackknife, svy sdr, svy
    // bootstrap, bayesboot bootstrap; seen on Stata 18 and 19);
    // the prefix shows the message and stops, so no replicate is run
    if `ccrep' & ("`bayesboot'" != "" | "`makeplot'" != "") {
        if "`bayesboot'" != "" {
            di as err "bayesboot() may not be combined with the `ccname' prefix"
        }
        if "`makeplot'" != "" {
            // the plotspec is user text (blocklabels()): shown as typed, without a second
            // macro expansion, its braces written as SMCL characters
            mata: st_local("mkp", st_local("makeplot2") == "" ? "makeplot" :                 ///
                "makeplot(" + subinstr(subinstr(subinstr(subinstr(st_local("makeplot2"), "{", ///
                char(1)), "}", char(2)), char(1), "{c -(}"), char(2), "{c )-}") + ")")
            di as err `"{p 0 4 2}makeplot may not be combined with the `ccname' prefix; afterward, type seldecomp, `macval(mkp)'{p_end}"'
        }
        exit 198
    }
    // the plotspec's own conflicts and what twoway would refuse (a misspelt option inside
    // twoway()) are refused here, before any work, rather than after the results are posted;
    // outcome() and blocklabels() once the outcomes and blocks are known
    if "`makeplot'" != "" {
        Plot_spec `macval(makeplot2)'
        if "`ps_chainstep'" != "" & "`chainprobe'" == "" {
            di as err "makeplot(chainstep()) requires chain"
            exit 198
        }
        if "`o_variant'" != "" & !inlist("`ps_type'", "", "marginal") {
            di as err "makeplot(`ps_type') not allowed with individually; use marginal"
            exit 198
        }
        if "`ps_reverse'" != "" & "`reverse'" == "" {
            di as err "makeplot(reverse) requires sequentially reverse"
            exit 198
        }
        if "`ps_ci'" != "" & "`bayesboot'" == "" {
            di as err "makeplot(ci) requires bayesboot()"
            exit 198
        }
        Plot_try `macval(ps_twoway)'
    }
    if `mincell' < 0 {
        di as err "mincell() must be a nonnegative integer"
        exit 198
    }
    if "`noisily'" != "" local noi noisily
    if "`o_byvar'" != "" {
        if "`o_generate'" != "" | "`makeplot'" != "" {
            di as err "generate() and makeplot may not be combined with the by prefix"
            exit 190
        }
    }

    // ---- the raking options -------------------------------------------------
    // the numbers are parsed as text so that e(rakeopts) records what the user typed (a real()
    // option would be re-formatted, e.g. 1e-10 -> 1.00000000000e-10)
    foreach o in tolerance ctrltolerance trimhirel trimhiabs trimlorel trimloabs {
        if "``o''" != "" {
            capture confirm number ``o''
            if _rc | missing(real("``o''")) {
                di as err "`o'() must be a number"
                exit 198
            }
            if ``o'' <= 0 {
                di as err "`o'() must be a positive number"
                exit 198
            }
            // a plus sign is dropped: the text is passed on as typed, and Mata has no unary plus
            if substr("``o''", 1, 1) == "+" local `o' = substr("``o''", 2, .)
        }
    }
    if `iterate' < 1 {
        di as err "iterate() must be a positive integer"
        exit 198
    }
    if "`tolerance'" == "" local tolerance 1e-8
    if "`ctrltolerance'" == "" local ctrltolerance 1e-6
    local ipfopts iterate(`iterate') tolerance(`tolerance') ctrltolerance(`ctrltolerance')
    // the trimming bounds (-1: absent); the absolute ones are rescaled inside Engine, the
    // relative ones pass through
    local anytrim
    foreach o in trimhirel trimhiabs trimlorel trimloabs {
        if "``o''" == "" local `o' -1
        else {
            if inlist("`o'", "trimhirel", "trimlorel") local ipfopts `ipfopts' `o'(``o'')
            local anytrim anytrim
        }
    }
    if `trimhiabs' != -1 & `trimloabs' != -1 & `trimhiabs' <= `trimloabs' {
        di as err "trimhiabs() must be greater than trimloabs()"
        exit 198
    }
    if `trimhirel' != -1 & `trimlorel' != -1 & `trimhirel' <= `trimlorel' {
        di as err "trimhirel() must be greater than trimlorel()"
        exit 198
    }
    if "`trimfrequency'" != "" {
        if !inlist("`trimfrequency'", "once", "sometimes", "often") {
            di as err "trimfrequency() must be once, sometimes, or often"
            exit 198
        }
        if "`anytrim'" == "" {
            di as err "trimfrequency() requires a trimming option"
            exit 198
        }
        local ipfopts `ipfopts' trimfrequency(`trimfrequency')
    }
    // trimsupportoverride: the Focal Group observations with a level absent from the Reference
    // Group get zero raked weight whatever the lower bounds; without it, trimloabs() and
    // trimlorel() bind on them as on any other observation (Engine: o_lobind)
    if "`trimsupportoverride'" != "" & `trimloabs' == -1 & `trimlorel' == -1 {
        di as err "trimsupportoverride requires trimloabs() or trimlorel()"
        exit 198
    }
    local o_tso `trimsupportoverride'
    local o_lobind = ("`trimsupportoverride'" == "" & (`trimloabs' != -1 | `trimlorel' != -1))
    local o_lbtxt = cond(`trimloabs' != -1 & `trimlorel' != -1, "trimloabs() and trimlorel()", ///
        cond(`trimloabs' != -1, "trimloabs()", "trimlorel()"))

    // ---- bayesboot() options ------------------------------------------------
    local bsreps 0
    if "`bayesboot'" != "" {
        // bayesboot(200), by analogy with reps(#) elsewhere: point to reps()
        if regexm(strtrim(`"`bayesboot2'"'), "^([0-9]+)([ ,]|$)") {
            local n = regexs(1)
            di as err "bayesboot(`n') not allowed; type bayesboot(reps(`n'))"
            exit 198
        }
        local 0 `", `bayesboot2'"'
        local spec Reps(integer 100) CLuster(varname) RSEED(string) bca ties ///
            SAVing(string asis) GENerate(string) Level(cilevel) CITable DOTS(integer 1) NODOTS ///
            NORMal TINTerval PERCentile BC
        capture syntax [, `spec' ]
        if _rc {
            Subopt_error bayesboot `"`spec'"' `"`macval(0)'"'
        }
        if `reps' < 2 {
            di as err "bayesboot(reps()) must be at least 2"
            exit 198
        }
        if `dots' == 0 local nodots nodots
        else if `dots' < 0 {
            di as err "bayesboot(dots()) must be a nonnegative integer"
            exit 198
        }
        if `dots' == 0 local dots 1
        if `"`rseed'"' != "" {
            capture confirm integer number `rseed'
            if _rc {
                di as err "bayesboot(rseed()) must be an integer"
                exit 198
            }
            if `rseed' < 0 | `rseed' > 2147483647 {
                di as err "bayesboot(rseed()) must be an integer between 0 and 2147483647"
                exit 198
            }
        }
        // level(): the confidence level of every interval; citable: the bootstrap interval table
        local o_level   `level'
        local o_citable `citable'
        // The interval types displayed (bstypes), initially in canonical order; explicit order
        // is restored below. Include those
        // requested; citable requests percentile and bc, and the default of normal and t (t
        // with cluster(), its degrees of freedom the clusters less 1; else normal) unless one
        // of them is requested; none requested: the default alone
        local bsdflt = cond("`cluster'" != "", "t", "normal")
        if "`normal'`tinterval'" != "" local bsdflt
        local bstypes
        foreach ty in normal t percentile bc bca {
            local tw = cond("`ty'" == "t", "tinterval", "`ty'")
            if "``tw''" != "" | ("`citable'" != "" & inlist("`ty'", "`bsdflt'", "percentile", "bc")) {
                local bstypes `bstypes' `ty'
            }
        }
        // Explicit interval types retain the user's order; the table and plot use the first.
        // citable-only additions follow them in the usual normal/t/percentile/bc/bca order.
        local bsrest : copy local bayesboot2
        local bsexplicit
        while `"`macval(bsrest)'"' != "" {
            gettoken bstoken bsrest : bsrest, bind
            mata: st_local("bty", ustrregexm(st_local("bstoken"), "^norm[a-z]*$") ? "normal" : ///
                (ustrregexm(st_local("bstoken"), "^tint[a-z]*$") ? "t" : ///
                (ustrregexm(st_local("bstoken"), "^perc[a-z]*$") ? "percentile" : ///
                (anyof(("bc", "bca"), st_local("bstoken")) ? st_local("bstoken") : ""))))
            if "`bty'" != "" & !`: list bty in bsexplicit' local bsexplicit `bsexplicit' `bty'
        }
        local bsremaining : list bstypes - bsexplicit
        local bstypes `bsexplicit' `bsremaining'
        local bstypes : list retokenize bstypes
        local bsdflt = cond("`cluster'" != "", "t", "normal")
        if "`bstypes'" == "" local bstypes `bsdflt'
        if "`makeplot'" != "" & "`ps_citype'" == "bca" & "`bca'" == "" {
            di as err "makeplot(bca): BCa intervals require bayesboot(bca)"
            exit 198
        }
        // generate(stub [, replace]): the replicate weights stub_1, stub_2, ...; replace
        // overwrites them
        local bsgenreplace
        if `"`generate'"' != "" {
            local 0 `"`generate'"'
            capture syntax namelist(name=generate max=1) [, REPLACE ]
            if _rc {
                di as err "bayesboot(generate()): type generate(stub) or generate(stub, replace)"
                exit 198
            }
            local bsgenreplace `replace'
        }
        local bsreps    `reps'
        local bsclust   `cluster'
        local bsseed    `"`rseed'"'
        local bsbca     `bca'
        local bsties    `ties'
        local bsgen     `generate'
        local bsdots    `dots'
        local bsnodots  `nodots'
        // bsseq: Bootstrap's option word for replicating every stage, as a staged run does
        local bsseq
        if "`o_sequential'" != "" local bsseq sequential
        if "`o_byvar'" != "" & (`"`saving'"' != "" | "`generate'" != "") {
            di as err "{p 0 4 2}bayesboot(saving()) and bayesboot(generate()) may not be combined with the by prefix{p_end}"
            exit 190
        }
        if `"`saving'"' != "" {
            local 0 `"`saving'"'
            syntax anything(name=bsfile id="filename") [, replace]
            gettoken bsfile bsrest : bsfile
            if `"`bsrest'"' != "" {
                di as err "{p 0 4 2}bayesboot(saving()): invalid file name; put names with spaces in double quotes{p_end}"
                exit 198
            }
            local bsfilereplace `replace'
            local bsfilechk `"`bsfile'"'
            mata: st_local("bsfilechk", pathsuffix(st_local("bsfilechk")) == "" ? st_local("bsfilechk") + ".dta" : st_local("bsfilechk"))
            if "`bsfilereplace'" == "" {
                confirm new file `"`bsfilechk'"'
            }
            else {
                // with replace, the file must be writable too, before any replication runs: an
                // existing file is opened for appending (nothing written), a new one is made
                // and erased
                capture confirm file `"`bsfilechk'"'
                local bsexist = (_rc == 0)
                tempname bsfh
                capture file open `bsfh' using `"`bsfilechk'"', write append
                if _rc {
                    di as err `"{p 0 4 2}bayesboot(saving()): file `bsfilechk' cannot be written{p_end}"'
                    exit 603
                }
                file close `bsfh'
                if !`bsexist' erase `"`bsfilechk'"'
            }
        }
    }
    local generate    `o_generate'
    local replace     `o_replace'
    local level       `o_level'
    local citable     `o_citable'
    local ratio       `o_ratio'
    local outcomes    `o_outcomes'
    local reverse     `o_reverse'
    local detail      `o_detail'

    // ---- contvars() option --------------------------------------------------
    // recut: the cut points are recomputed in every replicate, of bayesboot() and of a prefix
    local contsample pooled
    local ncont 0
    local recut
    if `"`contvars'"' != "" {
        local 0 `"`contvars'"'
        capture syntax anything(id="contvars()") [, POOLED FOCALGROUP REFGROUP RECUT ]
        if _rc {
            Subopt_error contvars "POOLED FOCALGROUP REFGROUP RECUT" `"`macval(0)'"' `"anything(id="contvars()")"'
        }
        if ("`pooled'" != "") + ("`focalgroup'" != "") + ("`refgroup'" != "") > 1 {
            di as err "contvars(): only one of pooled, focalgroup, and refgroup is allowed"
            exit 198
        }
        if "`focalgroup'" != "" local contsample focalgroup
        if "`refgroup'" != ""   local contsample refgroup
        local rest `anything'
        while `"`rest'"' != "" {
            gettoken cv rest : rest
            gettoken cq rest : rest
            if substr("`cv'", 1, 2) == "c." local cv = substr("`cv'", 3, .)
            capture unab cv : `cv'
            if _rc & inlist("`cv'", "pooled", "focalgroup", "refgroup", "recut") {
                local fix
                forvalues j = 1/`ncont' {
                    local fix `fix' `cvar`j'' `cq`j''
                }
                if `ncont' == 0 local fix age 4
                di as err "{p 0 4 2}contvars(): put a comma before `cv', as in contvars(`fix', `cv'){p_end}"
                exit 198
            }
            if _rc {
                di as err "contvars(): variable `cv' not found"
                exit 111
            }
            if `: word count `cv'' != 1 {
                di as err "contvars(): one variable name at a time (found `cv')"
                exit 198
            }
            if `: list cv in contlist' {
                di as err "contvars(): `cv' listed more than once"
                exit 198
            }
            capture confirm integer number `cq'
            if _rc {
                di as err "{p 0 4 2}contvars(): `cv' must be followed by the number of quantile groups{p_end}"
                exit 198
            }
            if `cq' < 2 {
                di as err "{p 0 4 2}contvars(): number of quantile groups for `cv' must be at least 2{p_end}"
                exit 198
            }
            // the number as it is recorded and shown: 4 for +4
            local cq = `cq'
            capture confirm numeric variable `cv'
            if _rc {
                di as err "contvars(): `cv' is not numeric"
                exit 108
            }
            local ++ncont
            local cvar`ncont' `cv'
            local cq`ncont'   `cq'
            local contlist `contlist' `cv'
        }
        if `ncont' == 0 {
            di as err "contvars(): nothing specified"
            exit 198
        }
    }
    // restore outer values clobbered by the inner -syntax- calls
    local anything   `"`o_anything'"'
    local focalgroup `"`macval(o_focalgroup)'"'
    if "`recut'" != "" & "`o_insupport'" != "" & `bsreps' > 0 {
        di as err "contvars(..., recut) may not be combined with insupport and bayesboot()"
        exit 198
    }
    if "`recut'" != "" & "`chain'" != "" & `bsreps' > 0 {
        di as err "contvars(..., recut) may not be combined with chain and bayesboot()"
        exit 198
    }
    // chain, pair mode: the replicates are Chain's (Bootstrap over the pairs), the
    // bayesboot(generate()) variables too
    if `chainpair' > 0 {
        local bsreps 0
        local bsgen
    }
    // contvars() under a replication prefix: a replicate call
    // (ccloop, above) uses the cut points of the prefix's full-sample call unless
    // contvars(..., recut) is given; the full-sample call records its cut points under a key
    // made of the specification (seldecomp_cutcache, at the end of Estimate)
    mata: st_local("cckey", invtokens((st_local("o_anything"), st_local("by"),               ///
        st_local("o_focalgroup"), st_local("contvars"),                                      ///
        st_local("o_ratio"), st_local("o_outcomes"), st_local("o_insupport"),               ///
        st_local("o_byvar"), st_local("o_byidx")), char(9)))
    // chain: each pair has its own cut points (the probe bins nothing)
    if `chainpair' > 0 mata: st_local("cckey", st_local("cckey") + char(9) + "chainpair " + st_local("chainpair"))
    local ccfixed 0
    if `ccloop' & "`chainprobe'" == "" {
        if `ncont' > 0 & "`recut'" == "" {
            local ccmats
            forvalues i = 1/`ncont' {
                tempname ccm`i'
                local ccmats `ccmats' `ccm`i''
            }
            mata: seldecomp_cutcache("get", st_local("contlist"), st_local("ccmats"))
            if "`ccerr'" != "" {
                if "`ccerr'" == "ambiguous" {
                    local cctxt "the cut points of this replicate are ambiguous"
                }
                else local cctxt "no cut points are recorded for this replicate"
                di as err `"{p 0 4 2}contvars(): `cctxt'; specify contvars(..., recut){p_end}"'
                exit 459
            }
            local ccfixed 1
        }
        else mata: seldecomp_cutcache("loop", "", "")
    }

    // ---- depvar and estimation sample ---------------------------------------
    gettoken depvar rakespec : anything
    if `"`rakespec'"' == "" {
        di as err "at least one raking variable is required"
        exit 102
    }
    if strpos("`depvar'", ".") | strpos("`depvar'", "#") {
        di as err "{p 0 4 2}dependent variable `depvar': factor-variable and time-series operators are not allowed{p_end}"
        exit 101
    }
    capture unab depvar : `depvar'
    if _rc {
        di as err "dependent variable `depvar' not found"
        exit 111
    }
    if `: word count `depvar'' != 1 {
        di as err "one dependent variable allowed"
        exit 103
    }
    capture confirm numeric variable `depvar'
    if _rc {
        di as err "dependent variable `depvar' must be numeric"
        exit 108
    }
    // outcomes(): the further outcomes, raked with the same weights; ratio(): the denominator
    local ylist `depvar'
    foreach v of local outcomes {
        capture confirm numeric variable `v'
        if _rc {
            di as err "outcomes(): `v' must be numeric"
            exit 108
        }
        if `: list v in ylist' {
            if "`v'" == "`depvar'" {
                di as err "{p 0 4 2}outcomes(): `v' is the dependent variable; list only the other outcomes{p_end}"
            }
            else di as err "outcomes(): `v' listed more than once"
            exit 198
        }
        local ylist `ylist' `v'
    }
    local nout : word count `ylist'
    // makeplot(outcome()): an outcome of these results, checked before any work
    if "`makeplot'" != "" & `"`macval(ps_outcome)'"' != "" {
        if !`: list ps_outcome in ylist' {
            if `nout' > 1 di as err "{p 0 4 2}makeplot(outcome()): `ps_outcome' is not one of the outcomes (`ylist'){p_end}"
            else {
                di as err "makeplot(outcome()): `ps_outcome' is not the outcome (`depvar')"
            }
            exit 198
        }
    }
    if "`ratio'" != "" {
        capture confirm numeric variable `ratio'
        if _rc {
            di as err "ratio(): `ratio' must be a numeric variable"
            exit 108
        }
    }

    // restore the sample-defining macros for marksample (inner -syntax- calls reset them)
    local if     `"`macval(o_if)'"'
    local in     `"`o_in'"'
    // the weights are left out of marksample: with novarlist it would check them before the
    // raking variables are marked out; missing, negative and zero weights are handled below
    local weight
    local exp
    marksample touse, novarlist
    markout `touse' `depvar'
    markout `touse' `by', strok
    // by prefix: keep only the current by-group (a no-op if marksample already did so)
    if "`o_byvar'" != "" qui replace `touse' = 0 if `o_byvar' != `o_byidx'

    // weights: internal double variable; missing, negative and zero weights are handled once
    // every other markout is done (below), so that observations outside the sample are not
    // checked, and the exclusions can be counted
    tempvar w
    if "`o_weight'" == "" qui gen double `w' = 1
    else qui gen double `w' `o_exp' if `touse'
    // data svyset with a pweight but no weights typed: say so
    local svynote 0
    if "`o_weight'" == "" & `chainpair' == 0 {
        capture qui svyset
        if _rc == 0 {
            if "`r(wvar)'" != "" {
                local svynote 1
                local svywvar `r(wvar)'
            }
        }
    }

    // ---- parse the raking specification into blocks and terms ---------------
    // Results (locals): nblocks, blk#name, blk#terms; nterms, term#vars, term#kinds, term#name
    local nblocks 0
    local nterms 0
    // a factor-variable operator with parentheses (ib(first).x, ib(#2).x, i(1 2).x) is refused
    // as the other operators are (Parse_term), before the list is split into terms
    mata: st_local("fvop", ustrregexm(st_local("rakespec"), ///
        "[A-Za-z_][A-Za-z0-9_]*[(][^()]*[)][.][A-Za-z_][A-Za-z0-9_]*") ? ustrregexs(0) : "")
    if `"`fvop'"' != "" {
        di as err `"{p 0 4 2}`fvop': only the i. and c. operators are allowed in the raking "' ///
            `"variable list{p_end}"'
        exit 198
    }
    local rest `"`rakespec'"'
    while `"`rest'"' != "" {
        gettoken tok rest : rest, match(paren)
        if "`paren'" == "(" {
            // a block: optional "name:" prefix
            local bname
            if strpos(`"`tok'"', ":") {
                gettoken bname tok : tok, parse(":")
                gettoken colon tok : tok, parse(":")
                local bname = strtrim("`bname'")
                if `: word count `bname'' != 1 {
                    di as err "invalid block name in (`bname': ...)"
                    exit 198
                }
                capture confirm name `bname'
                if _rc {
                    di as err "invalid block name `bname'; block names must be valid Stata names"
                    exit 198
                }
            }
            local tok = strtrim(`"`tok'"')
            if `"`tok'"' == "" {
                di as err "empty block () in the raking variable list"
                exit 198
            }
            local ++nblocks
            local blk`nblocks'name `bname'
            local blk`nblocks'terms
            local btoks `"`tok'"'
            while `"`btoks'"' != "" {
                gettoken t1 btoks : btoks
                Parse_term `t1'
                foreach k of numlist 1/`s_nnew' {
                    local ++nterms
                    local term`nterms'vars  `s_vars`k''
                    local term`nterms'kinds `s_kinds`k''
                    local term`nterms'name  `s_name`k''
                    local blk`nblocks'terms `blk`nblocks'terms' `nterms'
                }
            }
            if "`blk`nblocks'name'" == "" {
                local first : word 1 of `blk`nblocks'terms'
                local blk`nblocks'name : subinstr local term`first'name "#" "_", all
            }
        }
        else {
            Parse_term `tok'
            foreach k of numlist 1/`s_nnew' {
                local ++nterms
                local ++nblocks
                local term`nterms'vars  `s_vars`k''
                local term`nterms'kinds `s_kinds`k''
                local term`nterms'name  `s_name`k''
                local blk`nblocks'terms `nterms'
                local blk`nblocks'name : subinstr local s_name`k' "#" "_", all
            }
        }
    }
    // a raking variable may appear only once (an interaction in any order: x1#x2 = x2#x1)
    forvalues t = 1/`nterms' {
        local tsort : list sort term`t'vars
        forvalues t2 = 1/`=`t'-1' {
            local t2sort : list sort term`t2'vars
            if "`term`t'name'" == "`term`t2'name'" | "`tsort'" == "`t2sort'" {
                di as err "raking variable `term`t'name' specified more than once"
                exit 198
            }
        }
    }
    if `nterms' == 0 {
        di as err "at least one raking variable is required"
        exit 102
    }
    // the dependent variable cannot be raked on (checked before the c./contvars() rules, so
    // that c.depvar gets this message)
    forvalues t = 1/`nterms' {
        if `: list depvar in term`t'vars' {
            di as err "the dependent variable `depvar' cannot be a raking variable"
            exit 198
        }
        foreach v of local outcomes {
            if `: list v in term`t'vars' {
                di as err "outcomes(): `v' cannot be a raking variable"
                exit 198
            }
        }
    }
    // continuous variables must be in contvars(); contvars() variables must be c. raking variables
    local cvars
    forvalues t = 1/`nterms' {
        local j 0
        foreach v of local term`t'vars {
            local ++j
            local kind : word `j' of `term`t'kinds'
            if "`kind'" == "c" {
                local cvars `cvars' `v'
                if !`: list v in contlist' {
                    di as err "{p 0 4 2}c.`v' requires an entry in contvars(), e.g. contvars(`v' 4){p_end}"
                    exit 198
                }
            }
            else {
                if `: list v in contlist' {
                    di as err "`v' is in contvars() but is not marked c.`v'"
                    exit 198
                }
            }
            markout `touse' `v', strok
        }
    }
    foreach v of local contlist {
        if !`: list v in cvars' {
            di as err "`v' is in contvars() but is not a raking variable"
            exit 198
        }
    }
    // the group variable cannot be raked on (Focal Group has none of the Reference Group's value)
    forvalues t = 1/`nterms' {
        if `: list by in term`t'vars' {
            di as err "the by() variable `by' cannot be a raking variable"
            exit 198
        }
    }
    // the observations otherwise in the sample (if/in, depvar, by(), raking variables) that are
    // left out, in this order, for a missing value of an outcomes() or ratio() variable
    // (listwise: every outcome uses the same observations; xlw), a missing weight (xmw), a
    // missing bayesboot() cluster id (xmc); counted by group and noted below, as a
    // zero weight is
    tempvar xt xlw xmw xmc
    qui gen byte `xt' = `touse'
    markout `xt' `ylist' `ratio'
    qui gen byte `xlw' = `touse' & !`xt'
    qui replace `touse' = `xt'
    qui gen byte `xmw' = `touse' & missing(`w')
    qui replace `touse' = 0 if `xmw'
    qui replace `xt' = `touse'
    if "`bsclust'" != "" markout `xt' `bsclust', strok
    qui gen byte `xmc' = `touse' & !`xt'
    qui replace `touse' = `xt'
    // weights: positive on the estimation sample; zw marks the observations left out for a
    // zero weight alone (every other markout is done above), counted by group below
    local nzw 0
    if "`o_weight'" != "" {
        qui count if `touse' & `w' < 0
        if r(N) {
            di as err "negative weights encountered"
            exit 402
        }
        tempvar zw
        qui gen byte `zw' = `touse' & `w' == 0
        qui count if `zw'
        local nzw = r(N)
        qui replace `touse' = 0 if `zw'
    }
    // generated block names must be valid names of at most 32 characters (matrix stripes)
    forvalues b = 1/`nblocks' {
        capture confirm name `blk`b'name'
        if _rc | strlen("`blk`b'name'") > 32 {
            if "`sequential'" != "" & !`ccloop' & `chainpair' == 0 {
                di as txt "{p 0 6 2 78}note: block `b' is shown as block`b'; name it with (name: ...){p_end}"
            }
            local blk`b'name block`b'
        }
    }
    // with sequential, block names become equation names: unique and not reserved
    if "`sequential'" != "" {
        forvalues b = 1/`nblocks' {
            forvalues b2 = 1/`=`b'-1' {
                if "`blk`b'name'" == "`blk`b2'name'" {
                    di as err "duplicate block name `blk`b'name'"
                    exit 198
                }
            }
            if inlist("`blk`b'name'", "overall", "explained", "marg_expl", "cum_expl", "cum_unexpl", ///
                "raked_focal") | ///
                ("`o_variant'" != "" & "`blk`b'name'" == "_remainder") {
                di as err "`blk`b'name' cannot be used as a block name"
                exit 198
            }
        }
        // fixothermargins pins the other blocks at the Focal Group shares, which contradicts a target on
        // the same variable in another block (a block added later would then contribute exactly 0)
        if "`holdothers'" != "" {
            forvalues b = 1/`nblocks' {
                local bv`b'
                foreach t of local blk`b'terms {
                    local bv`b' `bv`b'' `term`t'vars'
                }
            }
            forvalues b = 1/`nblocks' {
                forvalues b2 = `=`b' + 1'/`nblocks' {
                    local common : list bv`b' & bv`b2'
                    if "`common'" != "" {
                        di as err "{p 0 4 2}with fixothermargins, a variable may be in one block only: `: word 1 of `common'' is in blocks `blk`b'name' and `blk`b2'name'{p_end}"
                        exit 198
                    }
                }
            }
            // detail: a block's later terms are held at the Focal Group shares while its earlier ones are
            // raked to the Reference Group, so the rule holds for the terms of a block too
            if "`detail'" != "" {
                forvalues b = 1/`nblocks' {
                    local bt `blk`b'terms'
                    local nbt : word count `bt'
                    forvalues i = 1/`nbt' {
                        local t1 : word `i' of `bt'
                        forvalues i2 = `=`i' + 1'/`nbt' {
                            local t2 : word `i2' of `bt'
                            local common : list term`t1'vars & term`t2'vars
                            if "`common'" != "" {
                                di as err "{p 0 4 2}with detail and fixothermargins, a variable may be in one term only: `: word 1 of `common'' is in `term`t1'name' and `term`t2'name'{p_end}"
                                exit 198
                            }
                        }
                    }
                }
            }
        }
    }
    // makeplot(blocklabels()): every name a block, checked before any work
    if "`makeplot'" != "" & `"`macval(ps_blocklabels)'"' != "" {
        local ps_blks
        forvalues b = 1/`nblocks' {
            local ps_blks `ps_blks' `blk`b'name'
        }
        Plot_blabels, blocks(`ps_blks') labels(`macval(ps_blocklabels)')
    }

    // ---- groups -------------------------------------------------------------
    qui count if `touse'
    if r(N) == 0 {
        // a replicate of a prefix: the cause recorded, no message (Misreps_warn)
        if `ccloop' {
            mata: seldecomp_failcache("add", "other", "no observations were left")
            exit 2000
        }
        if `nzw' > 0 {
            di as err "{p 0 4 2}no observations: "                                           ///
                cond(`nzw' == 1, "the one observation with nonmissing data has",              ///
                "all `nzw' observations with nonmissing data have") " zero weight{p_end}"
            exit 2000
        }
        error 2000
    }
    capture confirm string variable `by'
    local bystr = (_rc == 0)
    if "`chainprobe'" != "" & `bystr' {
        di as err "chain requires a numeric by() variable"
        exit 198
    }
    qui levelsof `by' if `touse', local(glev)
    local ng : word count `glev'
    // chain under a replication prefix: the full-sample call records the
    // values of by(); a replicate without one of them fails (its pairs would join values that
    // are not adjacent), each value missing recorded for the message after the loop
    if "`chainprobe'" != "" & `ccloop' {
        mata: seldecomp_chainlev("check")
        if !`clok' {
            if "`clmiss'" == "" {
                mata: seldecomp_failcache("add", "other", "the values of by(`by') differed from the full sample's")
            }
            foreach v of local clmiss {
                mata: seldecomp_failcache("add", "chain", "`by' = `v'")
            }
            exit 459
        }
    }
    if "`chainprobe'" != "" & `ccrep' & !`ccloop' & `ng' >= 2 mata: seldecomp_chainlev("put")
    // no input of the command may be overwritten by generate() or bayesboot(generate()), with or
    // without replace: among them the raking variables (gnrk) and the variables of the
    // weight expression (gnw); ginputs: the inputs by role, for Gen_role (and, below, under an
    // svy prefix the svyset design, whose weight is a temporary copy)
    local gnw
    local wtoks = ustrregexra(`"`o_exp'"', "[^A-Za-z0-9_]", " ")
    foreach tk of local wtoks {
        capture unab tv : `tk'
        if !_rc & `: word count `tv'' == 1 local gnw `gnw' `tv'
    }
    local gnrk
    forvalues t = 1/`nterms' {
        local gnrk `gnrk' `term`t'vars'
    }
    local ginputs `""the dependent variable" "`depvar'" "an outcomes() variable" "`outcomes'""'
    local ginputs `"`ginputs' "the ratio() variable" "`ratio'" "the by() variable" "`by'""'
    local ginputs `"`ginputs' "a contvars() variable" "`contlist'" "a raking variable" "`gnrk'""'
    local ginputs `"`ginputs' "the weight variable" "`gnw'" "the bayesboot(cluster()) variable" "`bsclust'""'
    // chain, the probe: two values or more; the plotspec's chainstep(), the
    // bayesboot(generate()) names, the exclusions noted by value of by() (one note per kind),
    // the estimation sample of the whole chain and the parsed options go to Chain
    if "`chainprobe'" != "" & `ng' >= 2 {
        if "`makeplot'" != "" {
            if "`ps_chainstep'" != "" {
                if real("`ps_chainstep'") > `ng' - 1 {
                    di as err "makeplot(chainstep()) must be 1 to `=`ng' - 1'"
                    exit 198
                }
            }
            else if inlist("`ps_type'", "stages", "waterfalllevels") {
                di as err "makeplot(`ps_type') requires chainstep(#) with chain"
                exit 198
            }
        }
        if "`bsgen'" != "" {
            forvalues r = 1/`bsreps' {
                if "`bsgenreplace'" == "" confirm new variable `bsgen'_`r'
                else Gen_role `bsgen'_`r' `"`ginputs'"' `"bayesboot(generate(`bsgen', replace))"'
            }
        }
        local nnote 0
        foreach k in xlw xmw xmc zw {
            if "``k''" == "" continue
            qui count if ``k''
            local nx = r(N)
            if `nx' == 0 continue
            if "`k'" == "xlw" local what "a missing value of `ratio' (the ratio() variable)"
            else if "`k'" == "xmw" local what "a missing weight"
            else if "`k'" == "zw" local what "zero weight"
            else {
                local what "a missing bayesboot() cluster id"
            }
            mata: st_local("xby", seldecomp_chaincnt("`by'", "``k''"))
            local ++nnote
            c_local ch_note`nnote' "`=strtrim(string(`nx', "%12.0fc"))' `=plural(`nx', "observation")' with `what' excluded (`xby')"
        }
        c_local ch_nnote `nnote'
        if `svynote' c_local ch_svynote "data are svyset with pweight `svywvar', but no weights were specified"
        c_local ch_nterms `nterms'
        c_local ch_nzw `nzw'
        foreach m in bsreps bsclust bsseed bsbca bsties bsgen bsgenreplace bsdots bsnodots bsseq ///
            bsfilereplace bstypes level ccloop ccrep ccpfx ccname o_force recut {
            c_local ch_`m' `"``m''"'
        }
        c_local ch_bsfile `"`bsfile'"'
        c_local ch_dispopts `diagnostics' `cformat' `pformat' `sformat' ///
            `citable' `varlabels'
        c_local ch_makeplot `makeplot'
        c_local ch_makeplot2 `"`macval(makeplot2)'"'
        qui gen byte `chainprobe' = `touse'
        exit
    }
    if "`chainprobe'" != "" & !`ccloop' {
        di as err "{p 0 4 2}by(`by') takes one value (`glev'); chain requires two or more{p_end}"
        exit 459
    }
    if `ng' != 2 {
        // a replicate of a prefix (a resample or deletion without one of the groups): the
        // cause recorded, no message (Misreps_warn)
        if `ccloop' {
            mata: seldecomp_failcache("add", "group", "`by'")
            exit 459
        }
        // one value left: say so when the other group's observations all have zero weight
        local zwhy
        if `ng' == 1 & `nzw' > 0 {
            local zopt none
            if `"`macval(focalgroup)'"' != "" local zopt focalgroup
            tempname ztxt
            mata: st_strscalar("`ztxt'", "`zopt'" == "none" ? "" : st_local("`zopt'"))
            Zero_group zwhy `by' `touse' `zw' `zopt' `ztxt'
        }
        if `ng' == 1 {
            // a string value in plain quotes, as Group_values shows it (levelsof adds compound quotes)
            local gvl `"`macval(glev)'"'
            if `bystr' Group_values gvl `by' `touse'
            if `"`macval(zwhy)'"' == "" {
                local zwhy "the other group has no observations with nonmissing data and a positive weight"
            }
            di as err `"{p 0 4 2}by(`by') takes one value (`macval(gvl)'): `macval(zwhy)'{p_end}"'
        }
        else {
            Group_values gvl `by' `touse'
            di as err `"{p 0 4 2}by(`by') takes `ng' values (`macval(gvl)'); it must take two{p_end}"'
        }
        exit 459
    }
    gettoken gv1 grest : glev
    gettoken gv2 grest : grest
    tempname gs1 gs2
    if !`bystr' {
        // levelsof (and macro text) store non-integer values with limited precision: the
        // comparisons use exact scalars; the macro text serves display and e(group_#)
        qui sum `by' if `touse', meanonly
        scalar `gs1' = r(min)
        scalar `gs2' = r(max)
        local gv1 = r(min)
        local gv2 = r(max)
    }
    else {
        // levelsof re-expands macros inside string values ($ and backquotes would vanish):
        // read the two values directly
        tempvar gcode
        qui egen long `gcode' = group(`by') if `touse'
        mata: st_local("gv1", st_sdata(selectindex(st_data(., "`gcode'") :== 1)[1], "`by'"))
        mata: st_local("gv2", st_sdata(selectindex(st_data(., "`gcode'") :== 2)[1], "`by'"))
    }
    local kfoc 1
    if `"`macval(focalgroup)'"' == "" {
        local gfocal `"`macval(gv1)'"'
    }
    else {
        if `bystr' {
            if `"`macval(focalgroup)'"' != `"`macval(gv1)'"' & `"`macval(focalgroup)'"' != `"`macval(gv2)'"' {
                Group_values gvl `by' `touse'
                tempname gtxt
                mata: st_strscalar("`gtxt'", st_local("focalgroup"))
                Zero_named zwhy `by' `gtxt' `zw'
                di as err `"{p 0 4 2}focalgroup(`macval(focalgroup)') is not a value of `by'`macval(zwhy)'; its values are `macval(gvl)'{p_end}"'
                exit 459
            }
            local gfocal `"`macval(focalgroup)'"'
        }
        else {
            tempname gtxt
            mata: st_strscalar("`gtxt'", st_local("focalgroup"))
            Group_pick kfoc focalgroup `by' `touse' `gs1' `gs2' `gtxt' `zw'
            local gfocal `gv`kfoc''
        }
    }
    // gsfocal/gsref: the scalars holding the exact focal and reference values of numeric by()
    if (!`bystr' & `kfoc' == 1) | (`bystr' & `"`macval(gfocal)'"' == `"`macval(gv1)'"') {
        local gref `"`macval(gv2)'"'
        local gsfocal `gs1'
        local gsref   `gs2'
    }
    else {
        local gref `"`macval(gv1)'"'
        local gsfocal `gs2'
        local gsref   `gs1'
    }
    tempvar g1
    if `bystr' {
        local gcf = cond(`"`macval(gfocal)'"' == `"`macval(gv1)'"', 1, 2)
        qui gen byte `g1' = (`gcode' == `gcf') if `touse'
    }
    else {
        qui gen byte `g1' = (`by' == `gsfocal') if `touse'
    }
    qui count if `touse' & `g1' == 1
    local N1 = r(N)
    qui count if `touse' & `g1' == 0
    local N2 = r(N)
    if `N1' == 0 | `N2' == 0 {
        if `ccloop' {
            mata: seldecomp_failcache("add", "group", "`by'")
            exit 2000
        }
        di as err "one of the groups has no observations in the estimation sample"
        exit 2000
    }
    // zero weights: the observations left out, by group (a third by() value is in neither);
    // no note in the replicate calls of a replication prefix, whose replicate weights are 0
    // for the units left out of a replicate (svy brr, jackknife, bootstrap, sdr)
    local zwnote
    if `nzw' > 0 {
        if `bystr' {
            // the values compared as data (macro text would re-expand $ and backquotes)
            mata: st_local("nzw1", strofreal(sum(st_sdata(., "`by'", "`zw'") :== st_local("gfocal"))))
            mata: st_local("nzw2", strofreal(sum(st_sdata(., "`by'", "`zw'") :== st_local("gref"))))
        }
        else {
            qui count if `zw' & `by' == `gsfocal'
            local nzw1 = r(N)
            qui count if `zw' & `by' == `gsref'
            local nzw2 = r(N)
        }
        local zwo
        if `nzw' - `nzw1' - `nzw2' > 0 local zwo ", `=strtrim(string(`nzw' - `nzw1' - `nzw2', "%12.0fc"))' with other values of `by'"
        if !`ccloop' {
            local zwnote "`=strtrim(string(`nzw', "%12.0fc"))' `=plural(`nzw', "observation")' with zero weight excluded (`=strtrim(string(`nzw1', "%12.0fc"))' in the focal group, `=strtrim(string(`nzw2', "%12.0fc"))' in the reference group`zwo')"
        }
    }
    // the other exclusions (xlw, xmw, xmc above), noted in the same way, before the zero weights
    local xnotes 0
    foreach k in xlw xmw xmc {
        if `ccloop' continue
        qui count if ``k''
        local nx = r(N)
        if `nx' == 0 continue
        if `bystr' {
            mata: st_local("nx1", strofreal(sum(st_sdata(., "`by'", "``k''") :== st_local("gfocal"))))
            mata: st_local("nx2", strofreal(sum(st_sdata(., "`by'", "``k''") :== st_local("gref"))))
        }
        else {
            qui count if ``k'' & `by' == `gsfocal'
            local nx1 = r(N)
            qui count if ``k'' & `by' == `gsref'
            local nx2 = r(N)
        }
        local xo
        if `nx' - `nx1' - `nx2' > 0 local xo ", `=strtrim(string(`nx' - `nx1' - `nx2', "%12.0fc"))' with other values of `by'"
        local xtail
        if "`k'" == "xlw" {
            if "`outcomes'" != "" & "`ratio'" != "" local what "a missing value of an outcomes() or the ratio() variable"
            else if "`outcomes'" != "" local what "a missing value of an outcomes() variable"
            else local what "a missing value of `ratio' (the ratio() variable)"
            if "`outcomes'" != "" local xtail "; every outcome uses the same observations"
        }
        else if "`k'" == "xmw" local what "a missing weight"
        else {
            local what "a missing bayesboot() cluster id"
        }
        local ++xnotes
        local xnote`xnotes' "`=strtrim(string(`nx', "%12.0fc"))' `=plural(`nx', "observation")' with `what' excluded (`=strtrim(string(`nx1', "%12.0fc"))' in the focal group, `=strtrim(string(`nx2', "%12.0fc"))' in the reference group`xo')`xtail'"
    }
    // group value labels for display
    local gfocallbl `"`macval(gfocal)'"'
    local greflbl   `"`macval(gref)'"'
    if !`bystr' {
        local vl : value label `by'
        if "`vl'" != "" & `gfocal' == int(`gfocal') & `gref' == int(`gref') {
            local l1 : label `vl' `gfocal', strict
            local l2 : label `vl' `gref', strict
            if `"`macval(l1)'"' != "" {
                local gfocallbl `"`gfocal' (`macval(l1)')"'
            }
            if `"`macval(l2)'"' != "" {
                local greflbl `"`gref' (`macval(l2)')"'
            }
        }
    }
    // the same for the error messages: user text, SMCL braces escaped
    foreach x in greflbl gfocallbl {
        local e`x' `"`macval(`x')'"'
        local e`x' : subinstr local e`x' "{" "`=char(1)'", all
        local e`x' : subinstr local e`x' "}" "{c )-}", all
        local e`x' : subinstr local e`x' "`=char(1)'" "{c -(}", all
    }

    // ---- build internal categorical raking variables ------------------------
    // each term t -> tempvar tv`t' coded 1..nlev`t' on the estimation sample,
    // level labels lev`t'_`l'
    // one: the variable of ones named in the column equations of the control matrices
    // (seldecomp_ctrlmat); cosmetic: the raking (seldecomp_rake) reads the column names and
    // the row name only
    tempvar one
    qui gen byte `one' = 1
    forvalues t = 1/`nterms' {
        local comps
        local j 0
        foreach v of local term`t'vars {
            local ++j
            local kind : word `j' of `term`t'kinds'
            capture confirm string variable `v'
            if _rc == 0 {
                if "`kind'" == "c" {
                    di as err "c.`v': `v' is a string variable"
                    exit 109
                }
                tempvar s`t'_`j'
                qui egen long `s`t'_`j'' = group(`v') if `touse'
                local comps `comps' `s`t'_`j''
                local ckind`t'_`j' s
            }
            else if "`kind'" == "c" {
                // continuous: bin by quantiles
                local ci : list posof "`v'" in contlist
                local q `cq`ci''
                tempvar bin`t'_`j'
                tempname cm`t'_`j'
                local fixopt
                if `ccfixed' local fixopt fixed(`ccm`ci'')
                // a replicate of a prefix (recut): the cause recorded, no message (Misreps_warn)
                if `ccloop' {
                    capture Bin_quantiles `v', touse(`touse') g1(`g1') w(`w') q(`q') sample(`contsample') ///
                        gen(`bin`t'_`j'') cutmat(`cm`t'_`j'') `fixopt'
                    if _rc == 2001 {
                        mata: seldecomp_failcache("add", "other", "`chp'contvars() found too few observations for the quantile groups of `v'")
                        exit 2001
                    }
                    else if _rc error _rc
                }
                else Bin_quantiles `v', touse(`touse') g1(`g1') w(`w') q(`q') sample(`contsample') ///
                    gen(`bin`t'_`j'') cutmat(`cm`t'_`j'') `fixopt'
                local cuts`t'_`j' `s_cuts'
                local ncuts`t'_`j' `s_ncuts'
                if !`: list v in cutvars' {
                    local cutvars `cutvars' `v'
                    local cutmats `cutmats' `cm`t'_`j''
                    local cutq    `cutq' `q'
                    // the nonempty bins (fewer than q with ties), for the header line
                    local cutnb   `cutnb' `s_nbins'
                }
                local comps `comps' `bin`t'_`j''
                local ckind`t'_`j' c
                local rcany`t' 1
            }
            else {
                local comps `comps' `v'
                local ckind`t'_`j' i
            }
        }
        // bayesboot() with recut: each term with a c. component, as its number and its
        // components, c.name for a binned one (Bootstrap bins it again in every replicate)
        if "`rcany`t''" != "" & "`recut'" != "" & `bsreps' > 0 {
            local rcc
            local j 0
            foreach cv of local comps {
                local ++j
                if "`ckind`t'_`j''" == "c" local rcc `rcc' c.`: word `j' of `term`t'vars''
                else local rcc `rcc' `cv'
            }
            if "`rcterms'" != "" local rcterms `rcterms' |
            local rcterms `rcterms' `t' `rcc'
        }
        // the levels numbered 1, 2, ... in the ascending order of the components, as
        // egen group() numbers them; row l of LV`t' holds level l's component values
        tempvar tv`t'
        tempname LV`t'
        qui gen long `tv`t'' = .
        mata: seldecomp_group("`comps'", "`touse'", "`tv`t''", "`LV`t''")
        local nlev`t' = rowsof(`LV`t'')
        if `nlev`t'' == 1 {
            di as txt "note: `chp'`term`t'name' is constant in the estimation sample"
        }
        // level labels
        local nv : word count `term`t'vars'
        forvalues l = 1/`nlev`t'' {
            local lab
            local j 0
            foreach v of local term`t'vars {
                local ++j
                local cv : word `j' of `comps'
                local val = `LV`t''[`l', `j']
                if "`ckind`t'_`j''" == "s" {
                    // the string value itself, read without macro expansion
                    mata: st_local("piece", st_sdata(selectindex(st_data(., "`cv'") :== `val')[1], "`v'"))
                }
                else if "`ckind`t'_`j''" == "c" {
                    Bin_label, bin(`val') ncuts(`ncuts`t'_`j'') cuts(`cuts`t'_`j'')
                    local piece `"`s_label'"'
                }
                else {
                    local vl : value label `v'
                    local piece `"`val'"'
                    // a non-integer level of a float variable: the shortest text of 7 to 9
                    // significant digits that is that float (.1, not .1000000014901161)
                    if "`: type `v''" == "float" & `val' != int(`val') {
                        foreach d in 7 8 9 {
                            local fv = strtrim(string(`LV`t''[`l', `j'], "%12.`d'g"))
                            if float(real("`fv'")) == `LV`t''[`l', `j'] {
                                local piece `fv'
                                continue, break
                            }
                        }
                    }
                    if "`vl'" != "" & `val' == int(`val') {
                        local vlab : label `vl' `val', strict
                        if `"`macval(vlab)'"' != "" {
                            local piece `"`val' (`macval(vlab)')"'
                        }
                    }
                }
                if `nv' > 1 {
                    local piece `"`v'=`macval(piece)'"'
                }
                if `"`macval(lab)'"' == "" {
                    local lab `"`macval(piece)'"'
                }
                else {
                    local lab `"`macval(lab)' `macval(piece)'"'
                }
            }
            // e(levels#) joins the labels with " | ", so a | inside a label is shown as /
            local lab : subinstr local lab "|" "/", all
            local lev`t'_`l' `"`macval(lab)'"'
        }
    }
    // insupport: every Focal Group observation with a level of any raking term that has no Reference Group
    // weight leaves the estimation sample, once, before anything is computed (Reference Group does not
    // change, so one pass removes them all); the levels are then renumbered on the screened
    // sample, the surviving ones keeping their order and labels (seldecomp_screen, _recode)
    local insnnote 0
    if "`o_insupport'" != "" {
        local tvl0
        local nll0
        forvalues t = 1/`nterms' {
            local tvl0 `tvl0' `tv`t''
            local nll0 `nll0' `nlev`t''
        }
        // Focal Group before the exclusion: count, weight, and each outcome's mean (a ratio with ratio())
        qui count if `touse' & `g1' == 1
        local N1all = r(N)
        qui sum `w' if `touse' & `g1' == 1, meanonly
        local W1all = r(sum)
        if "`ratio'" != "" {
            qui sum `ratio' [iw=`w'] if `touse' & `g1' == 1, meanonly
            local d1all = r(mean)
        }
        local j 0
        foreach v of local ylist {
            local ++j
            qui sum `v' [iw=`w'] if `touse' & `g1' == 1, meanonly
            local m1all_`j' = r(mean)
            if "`ratio'" != "" local m1all_`j' = cond(`d1all' != 0, `m1all_`j'' / `d1all', .)
        }
        tempvar insout
        qui gen byte `insout' = 0
        mata: seldecomp_screen("`tvl0'", "`nll0'", "`touse'", "`g1'", "`insout'")
        // a replicate of a prefix: the cause recorded, no message (Misreps_warn); a level that
        // no longer occurs in the Focal Group as a support failure
        if `ccloop' & `nexcl' == `N1all' {
            mata: seldecomp_failcache("add", "other", "insupport left no focal group observation")
            exit 2000
        }
        if `ccloop' & "`insfail'" != "" {
            gettoken ff : insfail
            gettoken ft fl : ff, parse(":")
            gettoken colon fl : fl, parse(":")
            local flab `"`macval(lev`ft'_`fl')'"'
            local flab : subinstr local flab "{" "`=char(1)'", all
            local flab : subinstr local flab "}" "{c )-}", all
            local flab : subinstr local flab "`=char(1)'" "{c -(}", all
            local flab `"`term`ft'name' = `macval(flab)'"'
            if `chainpair' > 0 mata: st_local("flab", st_local("flab") + " (pair `chainpair')")
            mata: seldecomp_failcache("add", "support", st_local("flab"))
            exit 459
        }
        if `nexcl' == `N1all' {
            // the culprit: the terms with such levels, with their Focal Group counts; a term
            // that excludes every Focal Group observation alone has no level found in both groups
            local tnm
            forvalues t = 1/`nterms' {
                local tnm `tnm' `term`t'name'
            }
            Zero_terms "`inslist'" `N1all' "`tnm'"
            if "`zsole'" != "" {
                di as err "{p 0 4 2}insupport excludes every focal group observation (levels of `term`zsole'name' absent from the reference group){p_end}"
            }
            else {
                di as err `"{p 0 4 2}insupport excludes every focal group observation (levels absent from the reference group: `zctxt'){p_end}"'
            }
            exit 2000
        }
        if "`insfail'" != "" {
            di as err `"{p 0 4 2}insupport: after the exclusion, these levels of the reference group (`by' = `macval(egreflbl)') have no focal group observations:{p_end}"'
            local nf 0
            foreach ff of local insfail {
                local ++nf
                if `nf' > 10 {
                    di as err "    ... and `=`: word count `insfail'' - 10' more level(s)"
                    continue, break
                }
                gettoken ft fl : ff, parse(":")
                gettoken colon fl : fl, parse(":")
                local flab `"`macval(lev`ft'_`fl')'"'
                local flab : subinstr local flab "{" "`=char(1)'", all
                local flab : subinstr local flab "}" "{c )-}", all
                local flab : subinstr local flab "`=char(1)'" "{c -(}", all
                di as err `"    `term`ft'name' = `macval(flab)'"'
            }
            exit 459
        }
        if `nexcl' > 0 {
            // the note: the levels absent from the Reference Group, each with its Focal Group count, at most 10
            // (the header line gives the count and the reason); one level on the note's line
            local nins : word count `inslist'
            local ++insnnote
            local insnote`insnnote' "insupport excluded:"
            local shown 0
            foreach zz of local inslist {
                local ++shown
                if `shown' > 10 {
                    local ++insnnote
                    local insnote`insnnote' "    ... (`=`nins' - 10' more)"
                    continue, break
                }
                gettoken zt zrest : zz, parse(":")
                gettoken colon zrest : zrest, parse(":")
                gettoken zl zrest : zrest, parse(":")
                gettoken colon zcount : zrest, parse(":")
                if `nins' == 1 {
                    local insnote`insnnote' `"insupport excluded `term`zt'name' = `macval(lev`zt'_`zl')' (`=strtrim(string(`zcount', "%12.0fc"))' obs)"'
                    continue
                }
                local ++insnnote
                local insnote`insnnote' `"    `term`zt'name' = `macval(lev`zt'_`zl')' (`=strtrim(string(`zcount', "%12.0fc"))' obs)"'
            }
            qui replace `touse' = 0 if `insout'
            local N1 = `N1all' - `nexcl'
            forvalues t = 1/`nterms' {
                local nl0 `nlev`t''
                mata: seldecomp_recode("`tv`t''", "`touse'", "`LV`t''")
                local nlev`t' = rowsof(`LV`t'')
                local l 0
                foreach l0 of local rkeep {
                    local ++l
                    local nlab`l' `"`macval(lev`t'_`l0')'"'
                }
                forvalues l = 1/`nlev`t'' {
                    local lev`t'_`l' `"`macval(nlab`l')'"'
                }
                if `nlev`t'' == 1 & `nl0' > 1 {
                    di as txt "note: `chp'`term`t'name' is constant in the estimation sample"
                }
            }
        }
    }
    // term lists per block, stage definitions
    forvalues b = 1/`nblocks' {
        local blk`b'desc
        foreach t of local blk`b'terms {
            local blk`b'desc `blk`b'desc' `term`t'name'
        }
    }
    if "`sequential'" != "" & "`o_variant'" != "" {
        // individually: raking k rakes block k alone to the Reference Group (Engine's holdothers holds the
        // other terms at the Focal Group shares with fixothermargins), raking K + 1 all blocks (the overall rows)
        local nstages = `nblocks' + 1
        local stg`nstages'terms
        forvalues s = 1/`nblocks' {
            local stg`s'terms `blk`s'terms'
            local stg`s'name  `blk`s'name'
            local stg`nstages'terms `stg`nstages'terms' `blk`s'terms'
        }
        local stg`nstages'name overall
    }
    else if "`sequential'" != "" {
        local nstages `nblocks'
        local cum
        forvalues s = 1/`nstages' {
            local cum `cum' `blk`s'terms'
            local stg`s'terms `cum'
            local stg`s'name  `blk`s'name'
        }
    }
    else {
        local nstages 1
        local stg1terms
        forvalues b = 1/`nblocks' {
            local stg1terms `stg1terms' `blk`b'terms'
        }
        local stg1name overall
    }
    forvalues s = 1/`nstages' {
        local stg`s'desc
        foreach t of local stg`s'terms {
            local stg`s'desc `stg`s'desc' `term`t'name'
        }
    }
    local tvlist
    local nlevlist
    forvalues t = 1/`nterms' {
        local tvlist `tvlist' `tv`t''
        local nlevlist `nlevlist' `nlev`t''
    }
    local stagedefs
    forvalues s = 1/`nstages' {
        if `s' > 1 local stagedefs `stagedefs' |
        local stagedefs `stagedefs' `stg`s'terms'
    }
    // reverse: the blocks added last to first. The terms are renumbered block-wise (blocks
    // K, ..., 1, each in its own order), so that the second Engine call does exactly what a
    // run with the blocks typed in reverse order does (targets and held margins in the same
    // order); rstg#name: the block added at stage # of the reverse order
    if "`reverse'" != "" {
        local rtvlist
        local rnlevlist
        local rtnames
        local rstagedefs
        // rtorig: the original number of each term of the reverse order; rcum#: the number of
        // terms raked to the Reference Group at stage # of the reverse order (its first rcum# terms)
        local rtorig
        local rcum 0
        forvalues s = 1/`nstages' {
            local bb = `nblocks' + 1 - `s'
            foreach t of local blk`bb'terms {
                local rtvlist `rtvlist' `tv`t''
                local rnlevlist `rnlevlist' `nlev`t''
                local rtnames `rtnames' `term`t'name'
                local rtorig `rtorig' `t'
            }
            local rcum = `rcum' + `: word count `blk`bb'terms''
            local rcum`s' `rcum'
            if `s' > 1 local rstagedefs `rstagedefs' |
            forvalues t = 1/`rcum' {
                local rstagedefs `rstagedefs' `t'
            }
            local rstg`s'name `blk`bb'name'
        }
    }
    // detail: the term chain, stage j raking terms 1 ... j (the terms are numbered in the order
    // typed, block by block); dends: the term stage at each block's last term, which is that
    // block's stage; dblk#: the block of term #; dname#: term #'s name in e(b), the term's name
    // with # replaced by _ (as default block names are made), else term#
    if "`detail'" != "" {
        local dstagedefs
        forvalues t = 1/`nterms' {
            if `t' > 1 local dstagedefs `dstagedefs' |
            forvalues u = 1/`t' {
                local dstagedefs `dstagedefs' `u'
            }
        }
        local dends
        forvalues b = 1/`nblocks' {
            foreach t of local blk`b'terms {
                local dblk`t' `b'
                local dlast `t'
            }
            local dends `dends' `dlast'
        }
        local dnames
        forvalues t = 1/`nterms' {
            local dn : subinstr local term`t'name "#" "_", all
            capture confirm name `dn'
            if _rc | strlen("`dn'") > 32 | `: list dn in dnames' {
                if !`ccloop' {
                    di as txt "{p 0 6 2 78}note: term `t' (`term`t'name') is shown as term`t' in detail{p_end}"
                }
                local dn term`t'
                while `: list dn in dnames' {
                    local dn `dn'_
                }
            }
            local dname`t' `dn'
            local dnames `dnames' `dn'
        }
    }
    // outcomes(): the equations of e(b) are prefixed by the outcome and must be unique
    if `nout' > 1 {
        local eqs
        foreach v of local ylist {
            local veqs `v'
            if "`o_variant'" != "" local veqs `veqs' `v'_marg_expl `v'_raked_focal
            else if "`sequential'" != "" {
                local veqs `veqs' `v'_marg_expl `v'_cum_expl `v'_cum_unexpl `v'_raked_focal
            }
            if "`detail'" != "" local veqs `veqs' `v'_detail
            if "`reverse'" != "" {
                local veqs `veqs' `v'_rev_marg_expl `v'_rev_cum_expl `v'_rev_cum_unexpl `v'_rev_raked_focal
            }
            foreach q of local veqs {
                if `: list q in eqs' {
                    di as err "{p 0 4 2}outcomes(): two outcomes give the equation name `q'; rename one{p_end}"
                    exit 198
                }
            }
            local eqs `eqs' `veqs'
        }
    }

    // ---- generate() names ---------------------------------------------------
    // under a replication prefix only the full-sample call makes the variables: the replicate
    // calls skip generate(), so what is left is the full-sample weights
    if "`generate'" != "" & !`ccloop' {
        local gennames `generate'
        if "`sequential'" != "" {
            forvalues s = 1/`=`nstages'-1' {
                local gennames `gennames' `generate'_s`s'
            }
        }
        local gnsvy
        if "`: word 1 of `ccpfx''" == "svy" {
            qui svyset
            local gnl `r(wvar)' `r(brrweight)' `r(jkrweight)' `r(sdrweight)' `r(bsrweight)'
            local gnl `gnl' `r(poststrata)' `r(postweight)'
            forvalues k = 1/`=max(1, r(stages))' {
                local gnl `gnl' `r(su`k')' `r(strata`k')' `r(fpc`k')'
            }
            local gncal `"`r(rake)' `r(regress)'"'
            gettoken gncal : gncal, parse(",")
            local gncal = ustrregexra(`"`gncal'"', "[#()]", " ")
            foreach tk of local gncal {
                local gnl `gnl' `= ustrregexra("`tk'", "^.*[.]", "")'
            }
            foreach tk of local gnl {
                capture unab tv : `tk'
                if !_rc local gnsvy `gnsvy' `tv'
            }
        }
        // variables made by an earlier full-sample call of the same prefix run, with the same
        // specification and estimation sample, are replaced (Stata 19's bayesboot prefix makes
        // a dry run first, and bootstrap's bca a jackknife: seldecomp_gencache); anything else
        // needs replace. Under a prefix, replace overwrites only a variable that generate()
        // made (characteristic seldecomp_cmd): the prefix's own variables (the cluster() and
        // strata() of bootstrap and jackknife, svy's subpop()) cannot be seen here, and the
        // prefix reads them after this call
        local gcok 0
        if `ccrep' mata: seldecomp_gencache("check", st_local("gennames"))
        if `gcok' local replace replace
        foreach gn of local gennames {
            capture confirm name `gn'
            if _rc {
                di as err "generate(): `gn' is not a valid variable name"
                exit 198
            }
            local role
            if "`gn'" == "`depvar'" local role "the dependent variable"
            else if `: list gn in outcomes' local role "an outcomes() variable"
            else if "`gn'" == "`ratio'" local role "the ratio() variable"
            else if "`gn'" == "`by'" local role "the by() variable"
            else if `: list gn in contlist' local role "a contvars() variable"
            else if `: list gn in gnrk' local role "a raking variable"
            else if `: list gn in gnw' local role "the weight variable"
            else if "`gn'" == "`bsclust'" local role "the bayesboot(cluster()) variable"
            else if `: list gn in gnsvy' local role "a variable of the svyset design"
            if "`role'" != "" {
                di as err "generate(): `gn' is `role'; choose another name"
                exit 198
            }
            if "`replace'" == "" confirm new variable `gn'
            else if `ccrep' & !`gcok' {
                capture confirm variable `gn', exact
                if !_rc {
                    mata: st_local("gnmade", strofreal(st_global("`gn'[seldecomp_cmd]") != ""))
                    if !`gnmade' {
                        di as err "{p 0 4 2}generate(..., replace): under the `ccname' prefix, replace overwrites only variables made by generate(); `gn' was not{p_end}"
                        exit 198
                    }
                }
            }
        }
    }
    if "`bsgen'" != "" {
        forvalues r = 1/`bsreps' {
            if "`bsgenreplace'" == "" confirm new variable `bsgen'_`r'
            else Gen_role `bsgen'_`r' `"`ginputs'"' `"bayesboot(generate(`bsgen', replace))"'
            if `: list posof "`bsgen'_`r'" in gennames' {
                di as err "{p 0 4 2}generate() and bayesboot(generate()) would both create `bsgen'_`r'{p_end}"
                exit 198
            }
        }
    }

    // ---- point estimates ----------------------------------------------------
    local tnames
    forvalues t = 1/`nterms' {
        local tnames `tnames' `term`t'name'
    }
    tempname RES RCHK RMISS
    forvalues s = 1/`nstages' {
        tempvar rw`s'
        qui gen double `rw`s'' = .
        local rwlist `rwlist' `rw`s''
    }
    local engopts trimhiabs(`trimhiabs') trimloabs(`trimloabs') mincell(`mincell')
    if "`o_tso'" != "" local engopts `engopts' trimsupportoverride
    // outcomes() and ratio() reach every Engine call (bootstrap replicates too) through engopts;
    // RESX#: the stage matrix of outcome # (# > 1), RESR/RESRX#: the same in reverse order
    local yextra : list ylist - depvar
    if `nout' > 1 local engopts `engopts' ylist(`yextra')
    if "`ratio'" != "" local engopts `engopts' den(`ratio')
    // individually: column 5 holds each raking's own explained part, the remainder last
    if "`o_variant'" != "" local engopts `engopts' contrib
    local resx
    local resrx
    forvalues j = 2/`nout' {
        tempname RESX`j' RESRX`j'
        local resx `resx' `RESX`j''
        local resrx `resrx' `RESRX`j''
    }
    // hold the user's previous estimation results while the decomposition runs: on any error
    // path below the held results are restored when Estimate exits
    tempname ehold
    _estimates hold `ehold', restore nullok
    Engine, touse(`touse') g1(`g1') w(`w') y(`depvar') one(`one') nterms(`nterms') ///
        tvlist(`tvlist') nlevlist(`nlevlist') nstages(`nstages') stagedefs(`stagedefs') ///
        res(`RES') rwlist(`rwlist') ipfopts(`ipfopts') tnames(`tnames') `engopts' `holdothers' `noi' ///
        resx(`resx') checks(`RCHK') misses(`RMISS')
    local rc = r(rc)
    if `rc' {
        if `rc' == 459 & !missing(r(nfail)) & r(nfail) > 0 {
            local nfail = r(nfail)
            local fstage = r(stage)
            local fzero = r(nzero)
            // a replicate of a prefix: no message; the
            // first failing level is recorded for the one message after the prefix's loop
            // (Misreps_warn), and the prefix drops the replicate
            if `ccloop' {
                local ft = r(fail1_t)
                local fl = r(fail1_l)
                // the label is user text: SMCL braces escaped
                local flab `"`macval(lev`ft'_`fl')'"'
                local flab : subinstr local flab "{" "`=char(1)'", all
                local flab : subinstr local flab "}" "{c )-}", all
                local flab : subinstr local flab "`=char(1)'" "{c -(}", all
                local flab `"`term`ft'name' = `macval(flab)'"'
                if `chainpair' > 0 mata: st_local("flab", st_local("flab") + " (pair `chainpair')")
                mata: seldecomp_failcache("add", "support", st_local("flab"))
                exit `rc'
            }
            local fzlist `r(zerolist)'
            // the failing terms (fterms) and those of the set-aside observations (zterms)
            local fterms
            forvalues f = 1/`nfail' {
                local ft`f' = r(fail`f'_t)
                local fl`f' = r(fail`f'_l)
                local fz`f' = r(fail`f'_z)
                if !`: list ft`f' in fterms' local fterms `fterms' `ft`f''
            }
            local zterms
            foreach zz of local fzlist {
                gettoken zt : zz, parse(":")
                if !`: list zt in zterms' local zterms `zterms' `zt'
            }
            local nshow 10
            // the labels are user text: SMCL braces escaped (egreflbl and egfocallbl above)
            di as err `"{p 0 4 2}raking support failure"' cond(`nstages' > 1, " at stage `fstage'", "") ///
                `": these levels of the reference group (`by' = `macval(egreflbl)') have no focal group observations (`by' = `macval(egfocallbl)'):{p_end}"'
            forvalues f = 1/`nfail' {
                if `f' > `nshow' {
                    di as err "    ... and `=`nfail' - `nshow'' more " plural(`nfail' - `nshow', "level")
                    continue, break
                }
                local flab `"`macval(lev`ft`f''_`fl`f'')'"'
                local flab : subinstr local flab "{" "`=char(1)'", all
                local flab : subinstr local flab "}" "{c )-}", all
                local flab : subinstr local flab "`=char(1)'" "{c -(}", all
                if `fz`f'' == 1 {
                    di as err `"    `term`ft`f''name' = `macval(flab)' (only in set-aside focal group observations)"'
                }
                else {
                    di as err `"    `term`ft`f''name' = `macval(flab)'"'
                }
            }
            if `fzero' > 0 {
                // "another" raking variable only when none of them is a failing one
                local zother : list zterms & fterms
                local zother = cond("`zother'" == "", "another", "a")
                di as err "{p 0 4 2}`fzero' focal group " plural(`fzero', "observation") " set aside (levels absent from the reference group):{p_end}"
                local nz 0
                foreach zz of local fzlist {
                    local ++nz
                    if `nz' > `nshow' {
                        di as err "      ..."
                        continue, break
                    }
                    gettoken zt zrest : zz, parse(":")
                    gettoken colon zrest : zrest, parse(":")
                    gettoken zl zrest : zrest, parse(":")
                    gettoken colon zcount : zrest, parse(":")
                    local flab `"`macval(lev`zt'_`zl')'"'
                    local flab : subinstr local flab "{" "`=char(1)'", all
                    local flab : subinstr local flab "}" "{c )-}", all
                    local flab : subinstr local flab "`=char(1)'" "{c -(}", all
                    di as err `"      `term`zt'name' = `macval(flab)' (`zcount' obs)"'
                }
            }
            // the remedies: collapsing or dropping the levels, and raking the other group when it
            // has support (seldecomp_other: the check with the groups
            // swapped, on all terms; with insupport, excluded Focal Group observations have levels
            // that Reference Group lacks, so raking Reference Group would fail)
            if 1 {
                local rem "collapse or drop these levels"
                // the other group: focalgroup() with the value of the Reference Group
                local oopt `"focalgroup(`macval(gref)')"'
                local oin
                // user text: a brace block (a one-line if would expand it again)
                if `"`macval(o_focalgroup)'"' != "" {
                    local oin `" in place of focalgroup(`macval(o_focalgroup)')"'
                }
                // (chain: each pair's Focal Group is its lower value; nothing is said)
                local otherok = cond(`chainpair' > 0, ".", "0")
                if !("`o_insupport'" != "" & 0`nexcl' > 0) & `chainpair' == 0 {
                    // cut points from one group's sample (contvars(, focalgroup or refgroup)):
                    // raking the other group cuts at the other group's quantiles, so its binned
                    // terms are coded again as that run would code them; in a replicate call
                    // with the full-sample call's cut points, those of the other run are not
                    // known here, and nothing is said about it (otherok missing)
                    local otvl `tvlist'
                    local onll `nlevlist'
                    if "`contsample'" != "pooled" & "`cutvars'" != "" {
                        if `ccfixed' local otherok .
                        else {
                            local osample = cond("`contsample'" == "refgroup", "focalgroup", "refgroup")
                            local otvl
                            local onll
                            forvalues t = 1/`nterms' {
                                if "`rcany`t''" == "" {
                                    local otvl `otvl' `tv`t''
                                    local onll `onll' `nlev`t''
                                    continue
                                }
                                local ocomps
                                local j 0
                                foreach v of local term`t'vars {
                                    local ++j
                                    if "`ckind`t'_`j''" == "c" {
                                        local ci : list posof "`v'" in contlist
                                        tempvar ob`t'_`j'
                                        // too few observations in the other group for the
                                        // quantile groups: that run fails (2001; otherok .z)
                                        capture Bin_quantiles `v', touse(`touse') g1(`g1') w(`w') ///
                                            q(`cq`ci'') sample(`osample') gen(`ob`t'_`j'')
                                        if _rc local otherok .z
                                        local ocomps `ocomps' `ob`t'_`j''
                                    }
                                    else if "`ckind`t'_`j''" == "s" local ocomps `ocomps' `s`t'_`j''
                                    else local ocomps `ocomps' `v'
                                }
                                if "`otherok'" == ".z" continue, break
                                tempvar otv`t'
                                tempname OLV`t'
                                qui gen long `otv`t'' = .
                                mata: seldecomp_group("`ocomps'", "`touse'", "`otv`t''", "`OLV`t''")
                                local otvl `otvl' `otv`t''
                                local onll `onll' `=rowsof(`OLV`t'')'
                            }
                        }
                    }
                    if "`otherok'" == "0" {
                        mata: seldecomp_other("`otvl'", "`onll'", "`touse'", "`g1'")
                    }
                }
                if "`otherok'" == "1" {
                    local rem `"`rem', or rake the other group with `macval(oopt)'`macval(oin)'"'
                }
            }
            di as err `"{p 0 4 2}`macval(rem)'{p_end}"'
        }
        else if `"`r(msg1)'"' != "" {
            // a replicate of a prefix: recorded, no message (as above)
            if `ccloop' {
                mata: seldecomp_failcache("add", "other", st_local("chp") + st_global("r(msg1)"))
                exit `rc'
            }
            di as err `"{p 0 4 2}`chp'`r(msg1)'{p_end}"'
        }
        exit `rc'
    }
    local m1  = r(m1)
    local m2  = r(m2)
    local W1  = r(W1)
    local W2  = r(W2)
    local N   = `N1' + `N2'
    local diff = `m1' - `m2'
    // the further outcomes, and the numerators and denominators of ratio()
    forvalues j = 2/`nout' {
        local m1_`j' = r(m1_`j')
        local m2_`j' = r(m2_`j')
        local diff_`j' = `m1_`j'' - `m2_`j''
    }
    if "`ratio'" != "" {
        tempname RP
        matrix `RP' = J(3 * `nout', 2, .)
        local rpnames
        forvalues j = 1/`nout' {
            local yo : word `j' of `ylist'
            matrix `RP'[3 * `j' - 2, 1] = r(n1_`j')
            matrix `RP'[3 * `j' - 2, 2] = r(d1)
            matrix `RP'[3 * `j' - 1, 1] = r(n1r_`j')
            matrix `RP'[3 * `j' - 1, 2] = r(d1r)
            matrix `RP'[3 * `j', 1] = r(n2_`j')
            matrix `RP'[3 * `j', 2] = r(d2)
            local rpnames `rpnames' `yo':focal_group `yo':raked_focal `yo':ref_group
        }
        matrix rownames `RP' = `rpnames'
        matrix colnames `RP' = numerator denominator
    }
    // support and lower bounds (o_lobind): the Focal Group observations with a level absent from the
    // Reference Group stay in the raking, and the lower bounds keep them above zero weight, so their
    // margins (target 0) will likely be missed, which fails the stage: a warning, with the
    // levels (Engine's r(held#_*); the stage with the most observations has them all), unless
    // force keeps the results (svy shows only error-class output from the full-sample call)
    if r(nhs) > 0 & r(nhs) < . & "`o_force'" == "" & !`ccloop' {
        local hn 0
        local hl
        forvalues h = 1/`=r(nhs)' {
            if r(held`h'_n) > `hn' {
                local hn = r(held`h'_n)
                local hl `r(held`h'_list)'
            }
        }
        local htx
        local nhl : word count `hl'
        local k 0
        foreach zz of local hl {
            local ++k
            if `k' > 5 {
                local htx `"`macval(htx)'; and `=`nhl' - 5' more"'
                continue, break
            }
            gettoken zt zrest : zz, parse(":")
            gettoken colon zrest : zrest, parse(":")
            gettoken zl zrest : zrest, parse(":")
            gettoken colon zcount : zrest, parse(":")
            local sep = cond(`k' > 1, "; ", "")
            local htx `"`macval(htx)'`sep'`term`zt'name' = `macval(lev`zt'_`zl')': `=strtrim(string(`zcount', "%12.0fc"))' obs"'
        }
        local htx : subinstr local htx "{" "`=char(1)'", all
        local htx : subinstr local htx "}" "{c )-}", all
        local htx : subinstr local htx "`=char(1)'" "{c -(}", all
        local hcls = cond(`ccrep' & "`: word 1 of `ccpfx''" == "svy", "err", "txt")
        di as `hcls' `"{p 0 4 2 78}warning: `chp'`=strtrim(string(`hn', "%12.0fc"))' focal group "' ///
            plural(`hn', "observation") cond(`hn' == 1, " has", " have") ///
            `" a level absent from the reference group (`macval(htx)'). `o_lbtxt' "' ///
            cond(strpos("`o_lbtxt'", " and "), "keep ", "keeps ") cond(`hn' == 1, "its", "their") ///
            " raked " plural(`hn', "weight") " above 0, so raking will likely fail. trimsupportoverride sets " ///
            cond(`hn' == 1, "it", "them") " to 0.{p_end}"
    }
    // unmet targets: the raking did not produce the decomposition asked for, and only force
    // keeps the results, whatever kept the margins from their targets. tchg#: trimming changed
    // the stage's raked weights (Engine's r(trimchg#)), so the bounds are the cause, which the
    // error and the note say; supm#: a level absent from the Reference Group was missed
    // (r(supmiss#)). A replicate fails by the same rule (Bootstrap)
    local badstage 0
    forvalues s = 1/`nstages' {
        local tchg`s' = (r(trimchg`s') == 1)
        local supm`s' = (r(supmiss`s') > 0 & r(supmiss`s') < .)
        if `RES'[`s', 7] > 0 local badstage `s'
    }
    if `badstage' > 0 {
        // a replicate of a prefix: recorded, no message (as a support failure above); the
        // raking stopped at the iterate() limit (the remedy differs) with the stage, unless
        // the trimming bounds are the cause
        if "`o_force'" == "" & `ccloop' {
            if `RES'[`badstage', 6] == 0 & !`tchg`badstage'' mata: seldecomp_failcache("add", "iterate", "`badstage'")
            else mata: seldecomp_failcache("add", "targets", "")
            exit 430
        }
        if "`o_force'" == "" {
            local nbad = `RES'[`badstage', 7]
            // the worst margin of the stage (Engine's r(worst#_t), r(worst#_l)), labelled; for a
            // term raked to the Reference Group at the stage, with its raked Focal Group share and its target,
            // Reference Group's share (the level shares of the Reference Group sum to 1 over the levels raked: a level
            // without Focal Group observations is a support failure). A term held at the Focal Group's shares
            // (fixothermargins) is raked last in each cycle, so its margins are met unless trimming moved them
            local wt = r(worst`badstage'_t)
            local wl = r(worst`badstage'_l)
            local wtxt
            if `wt' < . & `wt' > 0 & `wl' < . {
                local wtxt `"; worst: `term`wt'name' = `macval(lev`wt'_`wl')'"'
                if `: list wt in stg`badstage'terms' {
                    local rwb : word `badstage' of `rwlist'
                    Worst_share `rwb' `tv`wt'' `wl' `touse' `g1' `w' `W2'
                    local wtxt `"`macval(wtxt)' (`wshare')"'
                }
            }
            // raising iterate() helps only when the raking stopped at its limit before the weights converged;
            // a missed level absent from the Reference Group: the lower bounds are the cause;
            // trimming changed the weights: the bounds are
            local sup `supm`badstage''
            Miss_error "`chp'" `nstages' `badstage' `nbad' `=`RES'[`badstage', 6] == 0' `"`macval(wtxt)'"' `sup' `tchg`badstage''
        }
        // with force, the results are kept and a note says so (below: one message, no separate
        // warning)
    }
    // converged: the weights converged and the targets were met; a miss kept by force makes
    // it 0. targets_met (e(targets_met)): every margin of every stage met.
    // rconv#: the raking's own convergence of stage #, for the notes
    local tmet 1
    forvalues s = 1/`nstages' {
        local rconv`s' = `RES'[`s', 6]
        if `RES'[`s', 7] > 0 local tmet 0
        if `RES'[`s', 7] > 0 {
            matrix `RES'[`s', 6] = 0
            foreach X of local resx {
                matrix `X'[`s', 6] = 0
            }
        }
    }
    local m1r   = `RES'[`nstages', 1]
    local expl  = `RES'[`nstages', 2]
    local unexpl = `RES'[`nstages', 3]
    local pct   = `RES'[`nstages', 4]
    forvalues j = 2/`nout' {
        local m1r_`j'    = `RESX`j''[`nstages', 1]
        local expl_`j'   = `RESX`j''[`nstages', 2]
        local unexpl_`j' = `RESX`j''[`nstages', 3]
        local pct_`j'    = `RESX`j''[`nstages', 4]
    }
    // raking notes from the point-estimate run (the cell lists are labelled here). A cell
    // (term, level, count) with zero raked weight or fewer than mincell() observations is
    // listed once, and the cells found at the same stages share one note. r() is read first.
    foreach kind in zero thin {
        local `kind'keys
        local nk = cond("`kind'" == "zero", r(nzs), r(nts))
        forvalues z = 1/`nk' {
            local zs = r(`kind'`z'_s)
            foreach zz in `r(`kind'`z'_list)' {
                local kk = subinstr("`zz'", ":", "_", .)
                if !`: list zz in `kind'keys' local `kind'keys ``kind'keys' `zz'
                local `kind'st_`kk' ``kind'st_`kk'' `zs'
            }
        }
    }
    // worst-fitting control of each stage, labelled; the stage's note (one note per cause,
    // saying what to do), mt# with the parts mt#a and mt#b of a
    // force note, between which "; also in reverse order" goes (reverse, below):
    //   margins missed and kept by force (those that trimming left unmet are said so;
    //   trimming that met every margin gets no note), or the raking stopped before the
    //   weights converged (the iterate() limit)
    local stoptxt "at the iterate() limit"
    forvalues s = 1/`nstages' {
        local wt = r(worst`s'_t)
        local wl = r(worst`s'_l)
        local worstvar`s'
        local worstcat`s'
        if `wt' < . & `wt' > 0 {
            local worstvar`s' `term`wt'name'
            local worstcat`s' `"`macval(lev`wt'_`wl')'"'
        }
        local mt`s'
        local mt`s'a
        local mt`s'b
        local nbad = `RES'[`s', 7]
        local plm = plural(`nbad', "margin")
        if `nbad' > 0 & `tchg`s'' {
            local ntr = `RES'[`s', 22]
            local ntx = cond(`ntr' > 0, "; " + strtrim(string(`ntr', "%12.0fc")) + " " + plural(`ntr', "weight") + " at a bound", "")
            local mt`s'a `"trimming left `nbad' raking `plm' unmet (worst: `worstvar`s'' = `macval(worstcat`s')'`ntx'"'
            local mt`s'b "); results kept because force was specified"
            local mt`s' `"`macval(mt`s'a)'`mt`s'b'"'
        }
        else if `nbad' > 0 {
            local stp = cond(`rconv`s'' == 0, "raking stopped `stoptxt' with ", "")
            local mt`s'a `"`stp'`nbad' raking `plm' not met (worst: `worstvar`s'' = `macval(worstcat`s')'"'
            local mt`s'b = cond(`rconv`s'' == 0, "), kept by force", "); results kept because force was specified")
            local mt`s' `"`macval(mt`s'a)'`mt`s'b'"'
        }
        else if `rconv`s'' == 0 {
            local mt`s' "raking stopped `stoptxt' before the weights converged; raise iterate()"
            local itstages `itstages' `s'
        }
    }
    local nnote 0
    // the observations left out for a missing value, then for a zero weight, first: the
    // earliest steps of the sample
    forvalues m = 1/`xnotes' {
        local ++nnote
        local note`nnote' "`xnote`m''"
    }
    if "`zwnote'" != "" {
        local ++nnote
        local note`nnote' "`zwnote'"
    }
    // insupport: the excluded levels, listed once
    forvalues m = 1/`insnnote' {
        local ++nnote
        local note`nnote' `"`macval(insnote`m')'"'
    }
    foreach kind in zero thin {
        local ng 0
        foreach zz of local `kind'keys {
            local kk = subinstr("`zz'", ":", "_", .)
            Note_stages sttxt `nstages' ``kind'st_`kk''
            local gi 0
            forvalues h = 1/`ng' {
                if "`gst`h''" == "`sttxt'" local gi `h'
            }
            if `gi' == 0 {
                local ++ng
                local gi `ng'
                local gst`gi' `sttxt'
                local gcells`gi'
            }
            local gcells`gi' `gcells`gi'' `zz'
        }
        forvalues h = 1/`ng' {
            local nc : word count `gcells`h''
            local stx = cond("`gst`h''" == "", "", ", `gst`h''")
            local ++nnote
            local nct = cond(`nc' > 1, "`nc' ", "")
            if "`kind'" == "zero" {
                local pll = plural(`nc', "level")
                local note`nnote' "`nct'`pll' absent from the reference group (zero raked weight)`stx':"
            }
            else {
                local plc = plural(`nc', "target level")
                local note`nnote' "`nct'thin `plc'`stx':"
            }
            local shown 0
            foreach zz of local gcells`h' {
                local ++shown
                if "`kind'" == "thin" & `shown' > 8 {
                    local ++nnote
                    local note`nnote' "    ... (`=`nc' - 8' more)"
                    continue, break
                }
                gettoken zt zrest : zz, parse(":")
                gettoken colon zrest : zrest, parse(":")
                gettoken zl zrest : zrest, parse(":")
                gettoken colon zcount : zrest, parse(":")
                local ++nnote
                local note`nnote' `"    `term`zt'name' = `macval(lev`zt'_`zl')' (`=strtrim(string(`zcount', "%12.0fc"))' obs)"'
            }
        }
    }
    // the stages' notes: the stages with the same note share it; fnote, the number of the last
    // force note, and its parts (fnotea, fnoteb)
    local fnote 0
    forvalues s = 1/`nstages' {
        if `"`macval(mt`s')'"' != "" {
            local sl `s'
            forvalues s2 = `=`s' + 1'/`nstages' {
                if `"`macval(mt`s2')'"' == `"`macval(mt`s')'"' {
                    local sl `sl' `s2'
                    local mt`s2'
                }
            }
            Note_stages sttxt `nstages' `sl'
            local stx = cond("`sttxt'" == "", "", "`sttxt': ")
            local ++nnote
            local note`nnote' `"`stx'`macval(mt`s')'"'
            // the note of the iterate() limit (its stages: itstages), for the bayesboot()
            // replications that stopped there too (below)
            if `: list s in itstages' {
                local inote `nnote'
                local istx "`stx'"
            }
            if `"`macval(mt`s'a)'"' != "" {
                local fnote `nnote'
                local fnotea `"`stx'`macval(mt`s'a)'"'
                local fnoteb `"`mt`s'b'"'
                local fnlist `fnlist' `nnote'
            }
        }
    }
    if `svynote' {
        local ++nnote
        local note`nnote' "data are svyset with pweight `svywvar', but no weights were specified"
    }

    // ---- reverse order ------------------------------------------------------
    // the blocks added last to first, on the same weights (the terms renumbered above); after
    // the notes, because this second Engine call replaces r()
    if "`reverse'" != "" {
        tempname RESR RCHKR RMISSR
        // the raked weights of each stage, for the message of a missed margin
        local rrwlist
        forvalues s = 1/`nstages' {
            tempvar rrw`s'
            qui gen double `rrw`s'' = .
            local rrwlist `rrwlist' `rrw`s''
        }
        Engine, touse(`touse') g1(`g1') w(`w') y(`depvar') one(`one') nterms(`nterms') ///
            tvlist(`rtvlist') nlevlist(`rnlevlist') nstages(`nstages') stagedefs(`rstagedefs') ///
            res(`RESR') rwlist(`rrwlist') ipfopts(`ipfopts') tnames(`rtnames') `engopts' `holdothers' `noi' ///
            resx(`resrx') checks(`RCHKR') misses(`RMISSR')
        local rc = r(rc)
        if `rc' {
            // a Break: no message, as in the forward order and the detail chain
            if `rc' == 1 exit 1
            // a replicate of a prefix: recorded, no message (as for the forward order)
            if `ccloop' {
                mata: seldecomp_failcache("add", "other", "reverse order: " + st_global("r(msg1)"))
                exit `rc'
            }
            if `"`r(msg1)'"' != "" di as err `"reverse order: `r(msg1)'"'
            else di as err "reverse order: raking failed (return code `rc')"
            exit `rc'
        }
        local rbad 0
        forvalues s = 1/`nstages' {
            local rtchg`s' = (r(trimchg`s') == 1)
            if `RESR'[`s', 7] > 0 local rbad `s'
        }
        if `rbad' > 0 {
            if "`o_force'" == "" & `ccloop' {
                // its stage, marked r (not one of the forward order's notes)
                if `RESR'[`rbad', 6] == 0 & !`rtchg`rbad'' mata: seldecomp_failcache("add", "iterate", "r`rbad'")
                else mata: seldecomp_failcache("add", "targets", "")
                exit 430
            }
            if "`o_force'" == "" {
                // as for the forward order; Engine's r(worst#_t) numbers the terms of the
                // reverse order (rtorig gives the original number)
                local wt = r(worst`rbad'_t)
                local wl = r(worst`rbad'_l)
                local wtxt
                if `wt' < . & `wt' > 0 & `wl' < . {
                    local wo : word `wt' of `rtorig'
                    local wtxt `"; worst: `term`wo'name' = `macval(lev`wo'_`wl')'"'
                    if `wt' <= `rcum`rbad'' {
                        local rwb : word `rbad' of `rrwlist'
                        Worst_share `rwb' `tv`wo'' `wl' `touse' `g1' `w' `W2'
                        local wtxt `"`macval(wtxt)' (`wshare')"'
                    }
                }
                Miss_error "reverse order: " `nstages' `rbad' `=`RESR'[`rbad', 7]' ///
                    `=`RESR'[`rbad', 6] == 0' `"`macval(wtxt)'"' 0 `rtchg`rbad''
            }
            // kept by force: said in the forward order's force note, else in a note of its own
            if `fnote' > 0 {
                local note`fnote' `"`macval(fnotea)'; also in reverse order`fnoteb'"'
                local fnotea `"`macval(fnotea)'; also in reverse order"'
            }
            else {
                local nrb = `RESR'[`rbad', 7]
                local wt = r(worst`rbad'_t)
                local wl = r(worst`rbad'_l)
                local wtxt
                if `wt' < . & `wt' > 0 & `wl' < . {
                    local wo : word `wt' of `rtorig'
                    local wtxt `"worst: `term`wo'name' = `macval(lev`wo'_`wl')'"'
                }
                // margins that trimming left unmet: said so, with the weights at a bound (as in
                // the forward order's note)
                local rstem `"`nrb' raking `=plural(`nrb', "margin")' not met"'
                if `rtchg`rbad'' {
                    local rstem `"trimming left `nrb' raking `=plural(`nrb', "margin")' unmet"'
                    local ntr = `RESR'[`rbad', 22]
                    if `ntr' > 0 & `ntr' < . {
                        if `"`macval(wtxt)'"' != "" {
                            local wtxt `"`macval(wtxt)'; "'
                        }
                        local wtxt `"`macval(wtxt)'`=strtrim(string(`ntr', "%12.0fc"))' `=plural(`ntr', "weight")' at a bound"'
                    }
                }
                if `"`macval(wtxt)'"' != "" {
                    local wtxt `" (`macval(wtxt)')"'
                }
                local ++nnote
                local note`nnote' `"reverse order: stage `rbad': `rstem'`macval(wtxt)'; results kept because force was specified"'
                local fnlist `fnlist' `nnote'
            }
        }
        forvalues s = 1/`nstages' {
            if `RESR'[`s', 7] > 0 {
                matrix `RESR'[`s', 6] = 0
                foreach X of local resrx {
                    matrix `X'[`s', 6] = 0
                }
            }
        }
        // Without trimming both final stages rake to the same targets in another order
        // and agree to tolerance. Intermediate clipping can make their final fits differ.
        if "`anytrim'" == "" & `badstage' == 0 & `rbad' == 0 {
            local tolchk 1e-8
            if "`tolerance'" != "" local tolchk `tolerance'
            local fd = reldif(`RES'[`nstages', 1], `RESR'[`nstages', 1])
            forvalues j = 2/`nout' {
                local fd = max(`fd', reldif(`RESX`j''[`nstages', 1], `RESRX`j''[`nstages', 1]))
            }
            if `fd' > 1000 * `tolchk' {
                local fdt = strtrim(string(`fd', "%8.2g"))
                local ++nnote
                local note`nnote' "the final stages of the two orders differ by `fdt' (relative); tighten tolerance()"
            }
        }
    }

    // ---- detail: the term chain ----------------------------------------------
    // the terms added one at a time, on the same weights; its stage at
    // each block's last term is that block's stage (the same raking), so only the other term
    // stages can fail where the block stages did not, and the notes are the block stages'. When
    // every block has one term, the chain is the block stages. After the notes: r() is replaced
    if "`detail'" != "" {
        tempname RESD RCHKD RMISSD
        local resdx
        forvalues j = 2/`nout' {
            tempname RESDX`j'
            local resdx `resdx' `RESDX`j''
        }
        if `nterms' == `nblocks' {
            matrix `RESD' = `RES'
            matrix `RCHKD' = `RCHK'
            capture matrix `RMISSD' = `RMISS'
            forvalues j = 2/`nout' {
                matrix `RESDX`j'' = `RESX`j''
            }
        }
        else {
            // the raked weights of each term stage, for the message of a missed margin
            local drwlist
            forvalues s = 1/`nterms' {
                tempvar drw`s'
                qui gen double `drw`s'' = .
                local drwlist `drwlist' `drw`s''
            }
            Engine, touse(`touse') g1(`g1') w(`w') y(`depvar') one(`one') nterms(`nterms') ///
                tvlist(`tvlist') nlevlist(`nlevlist') nstages(`nterms') stagedefs(`dstagedefs') ///
                res(`RESD') rwlist(`drwlist') ipfopts(`ipfopts') tnames(`tnames') `engopts' `holdothers' `noi' ///
                resx(`resdx') checks(`RCHKD') misses(`RMISSD')
            local rc = r(rc)
            if `rc' == 1 exit 1
            if `rc' {
                // Engine's message names the term stage: the term added there instead
                local dmsg `"`r(msg1)'"'
                mata: st_local("dsn", ustrregexm(st_local("dmsg"), "stage ([0-9]+)") ? ustrregexs(1) : "")
                if "`dsn'" != "" {
                    local dmsg : subinstr local dmsg "stage `dsn'" "adding `term`dsn'name'"
                }
                if `"`dmsg'"' == "" local dmsg "raking failed (return code `rc')"
                // a replicate of a prefix: recorded, no message (as for the reverse order)
                if `ccloop' {
                    mata: seldecomp_failcache("add", "other", "detail: " + st_local("dmsg"))
                    exit `rc'
                }
                di as err `"{p 0 4 2}detail: `dmsg'{p_end}"'
                exit `rc'
            }
            // a missed margin at a term stage that is not the last of its block (those are the
            // block stages, checked above); dtchg#: trimming changed the weights (as above)
            local dbad 0
            forvalues s = 1/`nterms' {
                local dtchg`s' = (r(trimchg`s') == 1)
                if !`: list s in dends' & `RESD'[`s', 7] > 0 local dbad `s'
            }
            if `dbad' > 0 {
                if "`o_force'" == "" & `ccloop' {
                    if `RESD'[`dbad', 6] == 0 & !`dtchg`dbad'' mata: seldecomp_failcache("add", "iterate", "d`dbad'")
                    else mata: seldecomp_failcache("add", "targets", "")
                    exit 430
                }
                local db `dblk`dbad''
                local wt = r(worst`dbad'_t)
                local wl = r(worst`dbad'_l)
                if "`o_force'" == "" {
                    // one message of at most two lines: the term, the block, the
                    // worst margin, the cause (with fixothermargins the block's later terms are held
                    // at the Focal Group's distributions) and the remedies; with fixothermargins dropping it
                    // first, since reordering cannot help when a later term is a function of the term
                    local nbd = `RESD'[`dbad', 7]
                    local wtxt
                    if `wt' < . & `wt' > 0 & `wl' < . {
                        local wtxt `" (worst: `term`wt'name' = `macval(lev`wt'_`wl')')"'
                    }
                    local stx = cond(`RESD'[`dbad', 6] == 0, " (raking stopped at iterate())", "")
                    local dfix "specify force to keep the results"
                    // trimming changed the weights: the bounds are the cause (as in Miss_error)
                    if `dtchg`dbad'' {
                        local stx
                        local dfix "the trimming bounds are binding: relax them, or `dfix'"
                    }
                    local dm `"detail: adding `term`dbad'name' (block `blk`db'name') misses `nbd' raking `=plural(`nbd', "margin")'`macval(wtxt)'`stx'; `dfix'"'
                    local dm : subinstr local dm "{" "`=char(1)'", all
                    local dm : subinstr local dm "}" "{c )-}", all
                    local dm : subinstr local dm "`=char(1)'" "{c -(}", all
                    di as err `"{p 0 4 2}`macval(dm)'{p_end}"'
                    exit 430
                }
                // kept by force: said in the forward order's force note, else in a note of its own
                if `fnote' > 0 {
                    local note`fnote' `"`macval(fnotea)'; also in the detail, adding `term`dbad'name'`fnoteb'"'
                }
                else {
                    local ndb = `RESD'[`dbad', 7]
                    local wtxt
                    if `wt' < . & `wt' > 0 & `wl' < . {
                        local wtxt `"worst: `term`wt'name' = `macval(lev`wt'_`wl')'"'
                    }
                    // margins that trimming left unmet: said so, with the weights at a bound
                    local dstem `"`ndb' raking `=plural(`ndb', "margin")' not met"'
                    if `dtchg`dbad'' {
                        local dstem `"trimming left `ndb' raking `=plural(`ndb', "margin")' unmet"'
                        local ntr = `RESD'[`dbad', 22]
                        if `ntr' > 0 & `ntr' < . {
                            if `"`macval(wtxt)'"' != "" {
                                local wtxt `"`macval(wtxt)'; "'
                            }
                            local wtxt `"`macval(wtxt)'`=strtrim(string(`ntr', "%12.0fc"))' `=plural(`ntr', "weight")' at a bound"'
                        }
                    }
                    if `"`macval(wtxt)'"' != "" {
                        local wtxt `" (`macval(wtxt)')"'
                    }
                    local ++nnote
                    local note`nnote' `"detail, adding `term`dbad'name' (block `blk`db'name'): `dstem'`macval(wtxt)'; results kept because force was specified"'
                    local fnlist `fnlist' `nnote'
                }
            }
            forvalues s = 1/`nterms' {
                if `RESD'[`s', 7] > 0 {
                    matrix `RESD'[`s', 6] = 0
                    foreach X of local resdx {
                        matrix `X'[`s', 6] = 0
                    }
                }
            }
        }
    }
    // margins kept by force under svy, whose output shows neither the notes nor the
    // diagnostics: the force note once, in the full-sample call (error class, the only output
    // svy shows from it); the other prefixes display the notes through seldecomp's replay
    if "`fnlist'" != "" & `ccrep' & !`ccloop' & "`: word 1 of `ccpfx''" == "svy" {
        foreach m of local fnlist {
            local nt `"`macval(note`m')'"'
            local nt : subinstr local nt "{" "`=char(1)'", all
            local nt : subinstr local nt "}" "{c )-}", all
            local nt : subinstr local nt "`=char(1)'" "{c -(}", all
            di as err `"{p 0 4 2 78}warning: `chp'`macval(nt)'{p_end}"'
        }
    }

    // ---- pre-flight check of an svy prefix's replicates ----------------------
    // in the full-sample call of svy brr, jackknife, sdr or bootstrap
    // (not in the replicates), the replicates that will fail with a support error are found
    // from the svyset replicate weights (or, for svy jackknife without jkrweight(), from the
    // PSUs) before the prefix's loop, and the levels named, once; e(N_prefail). svy brr with
    // hadamard() builds its half-samples from a matrix this call cannot see: not checked. The
    // output is error-class, the only kind the prefix shows from this call; it never stops.
    // chain: made once, over the pairs, by Chain
    local prefail
    if `ccrep' & !`ccloop' & "`: word 1 of `ccpfx''" == "svy" & `chainpair' == 0 {
        local pfm : word 2 of `ccpfx'
        local pfrw
        local pfpsu
        if inlist("`pfm'", "brr", "jackknife", "sdr", "bootstrap") {
            local pfc = cond("`pfm'" == "jackknife", "jkr", cond("`pfm'" == "bootstrap", "bsr", "`pfm'"))
            local pfrw : char _dta[_svy_`pfc'weight]
            if "`pfrw'" != "" {
                capture confirm numeric variable `pfrw', exact
                local pfok = (_rc == 0)
                if !`pfok' local pfrw
                if `pfok' & "`pfm'" == "brr" {
                    // svy brr, hadamard() typed while brrweight() variables are svyset: the
                    // prefix uses the matrix and not the weights (its note says so), so
                    // checking them would describe replicates that are never run. That route
                    // is recognised from what the prefix has done before this call: sorted the
                    // data by the svyset strata and PSU variables and made two working
                    // variables, the stratum's index in that order and the PSU's number (1 or
                    // 2) within it; then nothing is checked
                    capture mata: seldecomp_hadroute("`: char _dta[_svy_strata1]'", ///
                        "`: char _dta[_svy_su1]'", "`ccsort'")
                    if !_rc & "`pf_had'" == "1" local pfrw
                }
            }
            else if "`pfm'" == "jackknife" {
                // the delete-one-PSU jackknife: the first-stage PSUs within strata
                local pfsu : char _dta[_svy_su1]
                local pfst : char _dta[_svy_strata1]
                tempvar pfpsu
                if inlist("`pfsu'", "", "_n") capture gen double `pfpsu' = _n if `touse'
                else capture egen double `pfpsu' = group(`pfst' `pfsu') if `touse'
                if _rc local pfpsu
            }
        }
        // the terms (with recut, those without a binned component: they are cut again in each
        // replicate)
        local pftv
        local pfnl
        local pfts
        forvalues t = 1/`nterms' {
            if "`recut'" != "" & "`rcany`t''" != "" continue
            local pftv `pftv' `tv`t''
            local pfnl `pfnl' `nlev`t''
            local pfts `pfts' `t'
        }
        if "`pfrw'`pfpsu'" != "" & "`pftv'" != "" {
            capture mata: seldecomp_prefail("`pftv'", "`pfnl'", "`touse'", "`g1'", "`pfrw'", "`pfpsu'")
            if !_rc local prefail `pf_nfail'
            // svy jackknife skips a replicate-weight variable equal to the full-sample weight on
            // the estimation sample (the subpopulation's, with subpop()): the replicates it runs
            if "`prefail'" != "" & "`pfrw'" != "" & "`pfm'" == "jackknife" {
                mata: st_local("pf_R", strofreal(sum(colsum(st_data(., tokens("`pfrw'"), "`touse'") ///
                    :!= st_data(., "`w'", "`touse'")) :> 0)))
            }
        }
        if "`prefail'" == "" local pf_nfail 0
        if "`prefail'" == "" local pf_nvan 0
        // one warning: what fails the replicates and what to
        // do. A level that fails replicates (no Focal Group weight) or
        // vanishes from the Reference Group in replicates that run (its Focal Group observations get zero raked
        // weight there) is listed once, with both counts; one level that only fails replicates
        // is said in one sentence. What a failed replicate does to the standard errors is
        // svy's to say, in the notes under its table (so the warning says "will fail" and
        // nothing on what follows): they come from the other replicates, or are missing when the
        // PSU jackknife is left with one replicate in a stratum
        if `pf_nfail' > 0 | `pf_nvan' > 0 {
            local pfwhat = cond("`pfm'" == "brr", "BRR", cond("`pfm'" == "sdr", "SDR", "`pfm'"))
            local pffay = cond("`pfm'" == "brr", ", or use Fay's BRR (fay())", "")
            local lkeys
            foreach k in fail van {
                foreach ff of local pf_`k' {
                    tokenize "`ff'", parse(":")
                    local key `1'_`3'
                    if !`: list key in lkeys' local lkeys `lkeys' `key'
                    local pfn`k'_`key' `5'
                    local pfo`k'_`key' `7'
                }
            }
            local nkeys : word count `lkeys'
            // the text of each level, and with the PSU jackknife the PSU deleted in the first
            // replicate it fails or leaves without Reference Group weight (by its svyset values)
            foreach key of local lkeys {
                gettoken t l : key, parse("_")
                gettoken us l : l, parse("_")
                local t : word `t' of `pfts'
                local pflev`key' `"`term`t'name' = `macval(lev`t'_`l')'"'
                foreach k in fail van {
                    local pfw`k'_`key'
                    if "`pfn`k'_`key''" == "" continue
                    local nk = `pfn`k'_`key''
                    if "`pfpsu'" != "" {
                        local pfo `pfo`k'_`key''
                        if inlist("`pfsu'", "", "_n") local pfp "observation `pfo'"
                        else {
                            local pfv = `pfsu'[`pfo']
                            local pfp `"`pfsu' `pfv'"'
                            if "`pfst'" != "" {
                                local pfv = `pfst'[`pfo']
                                local pfp `"`macval(pfp)' of `pfst' `pfv'"'
                            }
                        }
                        local pfw`k'_`key' `"when `macval(pfp)' is deleted"'
                        if `nk' > 1 {
                            local pfw`k'_`key' `"`macval(pfw`k'_`key')', and in `=`nk' - 1' more `=plural(`nk' - 1, "replicate")'"'
                        }
                    }
                    else local pfw`k'_`key' "in `nk' `=plural(`nk', "replicate")'"
                }
            }
            local key : word 1 of `lkeys'
            if `nkeys' == 1 & `pf_nvan' == 0 & `pf_ng2' == 0 & ("`pfpsu'" == "" | `pf_nfail' == 1) {
                if "`pfpsu'" != "" {
                    local pfl `"`macval(pflev`key')' has no focal group weight `macval(pfwfail_`key')': that jackknife replicate will fail. Collapse that level."'
                }
                else {
                    local pfl `"`macval(pflev`key')' has no focal group weight in `pf_nfail' of `pf_R' `pfwhat' replicates, which will fail. Collapse that level`pffay'."'
                }
                local pflines 0
            }
            else {
                // with the PSU (or observation) jackknife, the number of replicates svy runs is
                // not known here (under subpop(), it has the PSUs of the whole strata): the
                // count alone; svy's note after the loop gives both
                local pfl
                if "`pfpsu'" != "" {
                    if `pf_nfail' > 0 {
                        local pfl "`pf_nfail' jackknife `=plural(`pf_nfail', "replicate")' will fail"
                    }
                }
                else if `pf_nfail' > 0 local pfl "`pf_nfail' of `pf_R' `pfwhat' replicates will fail"
                if `pf_nvan' > 0 {
                    // what such a replicate does to its Focal Group observations: zero raked
                    // weight; with the lower bounds (o_lobind) it fails, unless force keeps it
                    local vw1 "gives focal group obs zero weight"
                    local vw2 "give focal group obs zero weight"
                    if `o_lobind' & "`o_force'" == "" {
                        local vw1 "will fail (the lower trimming bounds keep focal group obs above zero weight)"
                        local vw2 `vw1'
                    }
                    else if `o_lobind' {
                        local vw1 "holds focal group obs at the lower trimming bound"
                        local vw2 "hold focal group obs at the lower trimming bound"
                    }
                    if `pf_nfail' > 0 {
                        local pfl "`pfl', and `pf_nvan' `=cond(`pf_nvan' == 1, "other `vw1'", "others `vw2'")'"
                    }
                    else if "`pfpsu'" != "" {
                        local pfl "`pf_nvan' jackknife `=cond(`pf_nvan' == 1, "replicate `vw1'", "replicates `vw2'")'"
                    }
                    else local pfl "`pf_nvan' of `pf_R' `pfwhat' replicates `vw2'"
                }
                if `nkeys' > 0 local pfl = "`pfl'; collapse " + cond(`nkeys' == 1, "this level", "these levels") + "`pffay'"
                local pfl "`pfl':"
                local pflines 1
            }
            // the labels and the PSU values are user text: SMCL braces escaped
            local pfl : subinstr local pfl "{" "`=char(1)'", all
            local pfl : subinstr local pfl "}" "{c )-}", all
            local pfl : subinstr local pfl "`=char(1)'" "{c -(}", all
            di as err `"{p 0 4 2 78}warning: `macval(pfl)'{p_end}"'
            if `pflines' {
                local nsh 0
                foreach key of local lkeys {
                    local ++nsh
                    if `nsh' > 10 {
                        di as err "{p 8 12 2 78}... and `=`nkeys' - 10' more{p_end}"
                        continue, break
                    }
                    // user text (labels, PSU values): brace blocks, not one-line ifs
                    local pfl `"`macval(pflev`key')':"'
                    local hasf = ("`pfnfail_`key''" != "")
                    if `hasf' {
                        local pfl `"`macval(pfl)' no focal group weight `macval(pfwfail_`key')'"'
                    }
                    if "`pfnvan_`key''" != "" {
                        local pfv `"`macval(pfwvan_`key')'"'
                        if `hasf' & "`pfpsu'" == "" local pfv "in `pfnvan_`key'' others"
                        if `hasf' {
                            local pfl `"`macval(pfl)',"'
                        }
                        local pfl `"`macval(pfl)' no reference group weight `macval(pfv)'"'
                    }
                    local pfl : subinstr local pfl "{" "`=char(1)'", all
                    local pfl : subinstr local pfl "}" "{c )-}", all
                    local pfl : subinstr local pfl "`=char(1)'" "{c -(}", all
                    di as err `"{p 8 12 2 78}`macval(pfl)'{p_end}"'
                }
                if `pf_ng2' > 0 {
                    di as err "{p 8 12 2 78}Reference group: no weight in `pf_ng2' `=plural(`pf_ng2', "replicate")'{p_end}"
                }
            }
        }
    }

    // ---- e(b) layout --------------------------------------------------------
    // per outcome: the overall equation (named after the outcome with outcomes()), then with
    // a staged mode the marginal explained parts (marg_expl), cumulative explained and unexplained
    // parts and raked Focal Group means (raked_focal) by stage, then with detail each term's
    // marginal explained part, then with reverse the stage equations in the reverse order (rev_);
    // outcomes one after the other
    tempname b bj
    local bnames
    forvalues j = 1/`nout' {
        local eq0 overall
        local pfx
        if `nout' > 1 {
            local eq0 : word `j' of `ylist'
            local pfx `eq0'_
        }
        if `j' == 1 {
            local R `RES'
            local RR `RESR'
            matrix `bj' = (`m1', `m1r', `m2', `diff', `expl', `unexpl')
        }
        else {
            local R `RESX`j''
            local RR `RESRX`j''
            matrix `bj' = (`m1_`j'', `m1r_`j'', `m2_`j'', `diff_`j'', `expl_`j'', `unexpl_`j'')
        }
        local bnames `bnames' `eq0':focal_group `eq0':raked_focal `eq0':ref_group `eq0':difference ///
            `eq0':explained `eq0':unexplained
        if "`sequential'" != "" {
            local ords 1
            if "`reverse'" != "" local ords 1 2
            foreach o of local ords {
                local M = cond(`o' == 1, "`R'", "`RR'")
                local rp = cond(`o' == 1, "", "rev_")
                local sn = cond(`o' == 1, "stg", "rstg")
                // the marginal explained parts, cumulative explained and unexplained parts and
                // raked Focal Group means by stage; individually: the contributions and the
                // remainder, and the raked mean of each block's own raking (the last raking, all
                // blocks, is the overall raked_focal)
                local cols 5 2 3 1
                if "`o_variant'" != "" local cols 5 1
                foreach c of local cols {
                    local eqn = cond(`c' == 5, "marg_expl", cond(`c' == 2, "cum_expl", ///
                        cond(`c' == 3, "cum_unexpl", "raked_focal")))
                    forvalues s = 1/`nstages' {
                        if "`o_variant'" != "" & `c' == 1 & `s' == `nstages' continue
                        matrix `bj' = (`bj', `M'[`s', `c'])
                        local sname ``sn'`s'name'
                        if "`o_variant'" != "" & `s' == `nstages' local sname _remainder
                        local bnames `bnames' `pfx'`rp'`eqn':`sname'
                    }
                }
                // detail: each term's increment, after the forward equations
                if `o' == 1 & "`detail'" != "" {
                    local D = cond(`j' == 1, "`RESD'", "`RESDX`j''")
                    forvalues t = 1/`nterms' {
                        matrix `bj' = (`bj', `D'[`t', 5])
                        local bnames `bnames' `pfx'detail:`dname`t''
                    }
                }
            }
        }
        if `j' == 1 matrix `b' = `bj'
        else matrix `b' = (`b', `bj')
    }
    matrix colnames `b' = `bnames'

    // ---- covariate detail ---------------------------------------------------
    if "`covariatedetail'" != "" {
        tempname COV
        Covariate_detail, touse(`touse') g1(`g1') w(`w') rw(`rw`nstages'') nterms(`nterms') ///
            tvlist(`tvlist') nlevlist(`nlevlist') res(`COV')
        local covrows `s_rows'
        matrix rownames `COV' = `covrows'
        matrix colnames `COV' = focal_group ref_group raked_focal N_focal N_ref
    }

    // ---- Bayesian bootstrap -------------------------------------------------
    if `bsreps' > 0 {
        tempname V BS
        local bsopts
        if "`bsclust'" != ""    local bsopts `bsopts' cluster(`bsclust')
        if `"`bsseed'"' != ""   local bsopts `bsopts' seed(`bsseed')
        if "`bsgen'" != ""      local bsopts `bsopts' bsgen(`bsgen')
        if "`bsgenreplace'" != "" local bsopts `bsopts' bsgenreplace
        if `"`bsfile'"' != ""   local bsopts `bsopts' bsfile(`"`bsfile'"') `bsfilereplace'
        // outcomes(): each outcome's observed explained share; reverse: the renumbered terms
        local pctlist `pct'
        forvalues j = 2/`nout' {
            local pctlist `pctlist' `pct_`j''
        }
        if `nout' > 1 local bsopts `bsopts' outcomes(`ylist')
        if "`reverse'" != "" {
            local bsopts `bsopts' reverse rtvlist(`rtvlist') rnlevlist(`rnlevlist') ///
                rtnames(`rtnames') rstagedefs(`rstagedefs')
        }
        // recut: the terms to re-bin, and the rule of contvars()
        if "`recut'" != "" {
            local rcqs
            forvalues i = 1/`ncont' {
                local rcqs `rcqs' `cq`i''
            }
            local bsopts `bsopts' recut(`rcterms') rcvars(`contlist') rcq(`rcqs') rcsample(`contsample')
        }
        if "`o_force'" != "" local bsopts `bsopts' force
        // detail: the replicates run the term chain
        if "`detail'" != "" & "`bsseq'" != "" {
            local bsopts `bsopts' dstagedefs(`dstagedefs') dends(`dends')
        }
        Bootstrap, touse(`touse') g1(`g1') w(`w') y(`depvar') one(`one') nterms(`nterms')    ///
            tvlist(`tvlist') nlevlist(`nlevlist') nstages(`nstages') stagedefs(`stagedefs')  ///
            b(`b') v(`V') bs(`BS') reps(`bsreps') dots(`bsdots') `bsnodots' ipfopts(`ipfopts') ///
            tnames(`tnames') level(`level') engopts(`engopts') pct(`pctlist')                 ///
            `holdothers' `bsseq' `bsbca' `bsties' `bsopts'
        local bsk       = r(k)
        local bsseed    = r(seed)
        local N_reps    = r(N_reps)
        local N_fail    = r(N_fail)
        local N_iter    = r(N_iter)
        local N_force   = r(N_force)
        local N_jkforce = r(N_jkforce)
        local N_jkfail  = r(N_jkfail)
        local pct_se    = r(pct_se)
        local pct_ll    = r(pct_ll)
        local pct_ul    = r(pct_ul)
        forvalues j = 2/`nout' {
            local pct_se_`j' = r(pct_se_`j')
            local pct_ll_`j' = r(pct_ll_`j')
            local pct_ul_`j' = r(pct_ul_`j')
        }
        local N_clust   = r(N_clust)
        local df_bs     = r(df_bs)
        local rngstate  `"`r(rngstate)'"'
        local bsnames   `r(bsnames)'
        // failed replications, and those kept at the iterate() limit (one message per cause).
        // When the point estimate stopped at that limit too,
        // at the stages where the replications did, and every failed replication stopped
        // there, its note counts them; else one note here. At replay a one-line footnote
        // (Display_table), none when the note counts them
        local nstop = r(N_failit) + `N_iter'
        local bsfold 0
        if "`inote'" != "" & `N_fail' == r(N_failit) & `nstop' > 0 & r(itrev) == 0 {
            local itbs `r(itstages)'
            local bsfold : list itbs in itstages
        }
        if `bsfold' {
            local dtx = cond(`N_fail' > 0, " (`N_fail' short of the targets, dropped)", "")
            local caution = cond(`N_fail' > 0, "; bootstrap results may be incorrect.", "")
            local note`inote' `"`istx'raking stopped `stoptxt', also in `nstop' of `bsreps' replications`dtx'; raise iterate()`caution'"'
        }
        else if `N_fail' > 0 & `N_fail' == r(N_failit) & `N_iter' > 0 {
            di as txt "{p 0 6 2 78}note: raking stopped `stoptxt' in `nstop' of `bsreps' replications; " ///
                "`N_fail' of them missed their targets and " cond(`N_fail' == 1, "was", "were") ///
                " dropped; raise iterate(); bootstrap results may be incorrect.{p_end}"
        }
        else if `N_fail' > 0 {
            di as txt "{p 0 6 2 78}note: `N_fail' of `bsreps' replications failed and " ///
                cond(`N_fail' == 1, "was", "were") " dropped: in " cond(`N_fail' == 1, "it", "them") ///
                `", `r(bswhy)'`r(bsfix)'; bootstrap results may be incorrect.{p_end}"'
        }
        else if `N_iter' > 0 {
            di as txt "{p 0 6 2 78}note: `N_iter' of `bsreps' replications stopped `stoptxt' with the " ///
                "raking targets met and " cond(`N_iter' == 1, "was", "were") " kept; raise iterate(){p_end}"
        }
        // force: the replications (and BCa jackknife deletions) kept although a margin was missed
        if `N_force' + `N_jkforce' > 0 {
            local jkt = cond(`N_jkforce' > 0, " and `N_jkforce' jackknife " + plural(`N_jkforce', "deletion"), "")
            di as txt "{p 0 6 2 78}note: force kept `N_force' of `bsreps' " plural(`bsreps', "replication") ///
                "`jkt' that missed raking targets{p_end}"
        }
        // few clusters: no note
    }

    // ---- generate() ---------------------------------------------------------
    if "`generate'" != "" & !`ccloop' {
        forvalues s = 1/`nstages' {
            if `s' == `nstages' local gn `generate'
            else local gn `generate'_s`s'
            // the exact name: with varabbrev on, drop takes an abbreviation (another variable)
            if "`replace'" != "" {
                capture confirm variable `gn', exact
                if !_rc drop `gn'
            }
            qui gen double `gn' = cond(`g1' == 1, `rw`s'', `w') if `touse'
            if "`o_variant'" != "" {
                if `s' < `nstages' {
                    local lbl = usubstr("Raked wts (Focal)/base wts (Ref.), block `s' alone (`o_variant'): `stg`s'name'", 1, 80)
                }
                else local lbl "Raked weights (focal group), base weights (reference group)"
            }
            else if "`sequential'" != "" {
                local lbl = usubstr("Raked wts (Focal)/base wts (Ref.), stage `s': `stg`s'name'", 1, 80)
            }
            else local lbl "Raked weights (focal group), base weights (reference group)"
            label variable `gn' "`lbl'"
            char `gn'[seldecomp_cmd] `"seldecomp `macval(cmdline0)'"'
            char `gn'[seldecomp_stage] `s'
            char `gn'[seldecomp_rakevars] `stg`s'desc'
            char `gn'[seldecomp_converged] `=`RES'[`s', 6]'
        }
        if `ccrep' mata: seldecomp_gencache("put", st_local("gennames"))
    }

    // chain, pair mode: the pair's estimation sample, Focal Group flag and coded terms, into the
    // variables Chain named (they live in its scope, for the replicates), and what Engine is
    // given besides; everything else Chain reads from the pair's e()
    if `chainpair' > 0 {
        gettoken cvt cvrest : chainvars
        gettoken cvg cvrest : cvrest
        qui gen byte `cvt' = `touse'
        qui gen byte `cvg' = `g1'
        forvalues t = 1/`nterms' {
            gettoken cv cvrest : cvrest
            qui gen long `cv' = `tv`t''
        }
        foreach m in nlevlist nstages stagedefs tnames engopts ipfopts holdothers {
            c_local ch_`m' `"``m''"'
        }
    }

    // ---- post ---------------------------------------------------------------
    if `bsreps' > 0 {
        _estimates unhold `ehold', not
        ereturn post `b' `V', depname(`depvar') obs(`N') esample(`touse')
    }
    else {
        _estimates unhold `ehold', not
        ereturn post `b', depname(`depvar') obs(`N') esample(`touse')
    }
    ereturn scalar N   = `N'
    ereturn scalar N_focal = `N1'
    ereturn scalar N_ref = `N2'
    ereturn scalar W_focal = `W1'
    ereturn scalar W_ref = `W2'
    ereturn scalar N_zerowt = `nzw'
    ereturn scalar mean_focal       = `m1'
    ereturn scalar mean_raked_focal = `m1r'
    ereturn scalar mean_ref       = `m2'
    ereturn scalar difference   = `diff'
    ereturn scalar explained    = `expl'
    ereturn scalar unexplained  = `unexpl'
    ereturn scalar pct_explained = `pct'
    ereturn local diffsign "focal_group - ref_group"
    ereturn scalar N_stages = `nstages'
    ereturn scalar k_eq = `nout' * (1 + cond("`o_variant'" != "", 2, 4) * ("`sequential'" != "") + 4 * ("`reverse'" != "") + ("`detail'" != ""))
    ereturn scalar level = `level'
    local diagnames converged badcontrols maxreldif maxctrl raked_mean raked_min raked_max     ///
        raked_sd raked_cv factor_mean factor_min factor_max factor_sd factor_cv N_raked N_focal_zero  ///
        N_at_bound N_eff_base N_eff_raked ess_ratio N_focal_mincell pct_increment pct_W_focal_zero       ///
        mean_focal_zero mean_focal_cs
    local c 5
    foreach d of local diagnames {
        local ++c
        ereturn scalar `d' = `RES'[`nstages', `c']
    }
    // converged: the minimum over the stages, and with reverse over the reverse order's too, with
    // detail over the term stages' (1 only when every raking converged; svy warns "Convergence
    // not achieved" at 0); column converged of e(bystage), e(bystage_rev) and
    // e(bystage_detail) keeps each stage's value
    local cmin .
    local cms `RES'
    if "`reverse'" != "" local cms `RES' `RESR'
    if "`detail'" != "" local cms `cms' `RESD'
    foreach M of local cms {
        forvalues s = 1/`=rowsof(`M')' {
            if `M'[`s', 6] < `cmin' local cmin = `M'[`s', 6]
        }
    }
    if `cmin' < . ereturn scalar converged = `cmin'
    ereturn scalar targets_met = `tmet'
    // individually: columns 2, 3, 5, 27 are each raking's own explained and unexplained
    // parts, the block's contribution (the remainder on the overall row) and its percent
    local c5names raked_focal cum_expl cum_unexpl pct_explained increment
    local dgnames `diagnames'
    if "`o_variant'" != "" {
        local c5names raked_focal explained unexplained pct_explained contribution
        local dgnames : subinstr local diagnames "pct_increment" "pct_contribution"
    }
    matrix colnames `RES' = `c5names' `dgnames'
    local stgnames
    forvalues s = 1/`nstages' {
        local stgnames `stgnames' `stg`s'name'
    }
    matrix rownames `RES' = `stgnames'
    ereturn matrix bystage = `RES', copy       // not e(stages): svy reserves that name for a scalar
    // outcomes(): e(bystage#) for outcome # > 1; reverse: e(bystage_rev), e(bystage_rev#)
    forvalues j = 2/`nout' {
        matrix colnames `RESX`j'' = `c5names' `dgnames'
        matrix rownames `RESX`j'' = `stgnames'
        ereturn matrix bystage`j' = `RESX`j''
    }
    if "`reverse'" != "" {
        local rstgnames
        forvalues s = 1/`nstages' {
            local rstgnames `rstgnames' `rstg`s'name'
        }
        local RR `RESR' `resrx'
        local j 0
        foreach M of local RR {
            local ++j
            matrix colnames `M' = raked_focal cum_expl cum_unexpl pct_explained increment `diagnames'
            matrix rownames `M' = `rstgnames'
            if `j' == 1 ereturn matrix bystage_rev = `M'
            else ereturn matrix bystage_rev`j' = `M'
        }
    }
    // detail: e(bystage_detail), e(bystage_detail#), one row per term stage
    if "`detail'" != "" {
        local RD `RESD' `resdx'
        local j 0
        foreach M of local RD {
            local ++j
            matrix colnames `M' = `c5names' `dgnames'
            matrix rownames `M' = `dnames'
            if `j' == 1 ereturn matrix bystage_detail = `M'
            else ereturn matrix bystage_detail`j' = `M'
        }
    }
    // Diagnostic sidecars preserve raw iteration convergence and level-target counts without
    // changing e(bystage)'s established columns. Miss labels are stored separately as text,
    // so replay and estimates save do not depend on temporary term variables or value labels.
    foreach order in forward reverse detail {
        if "`order'" != "forward" & "``order''" == "" continue
        local suffix
        local CK `RCHK'
        local MI `RMISS'
        local rn `stgnames'
        if "`order'" == "reverse" {
            local suffix _rev
            local CK `RCHKR'
            local MI `RMISSR'
            local rn `rstgnames'
        }
        else if "`order'" == "detail" {
            local suffix _detail
            local CK `RCHKD'
            local MI `RMISSD'
            local rn `dnames'
        }
        matrix colnames `CK' = iterations_converged missed_variables missed_targets targets_checked control_scale
        matrix rownames `CK' = `rn'
        ereturn matrix raking_checks`suffix' = `CK', copy
        capture confirm matrix `MI'
        if !_rc {
            forvalues m = 1/`=rowsof(`MI')' {
                local t = `MI'[`m', 2]
                if "`order'" == "reverse" {
                    local t : word `t' of `rtorig'
                    matrix `MI'[`m', 2] = `t'
                }
                local l = `MI'[`m', 3]
                ereturn local raking_miss`suffix'`m'_var `term`t'name'
                ereturn local raking_miss`suffix'`m'_level `"`macval(lev`t'_`l')'"'
            }
            matrix colnames `MI' = stage term level target achieved discrepancy
            ereturn matrix raking_misses`suffix' = `MI', copy
        }
    }
    if `nout' > 1 {
        tempname BO
        // insupport: an 11th column, the Focal Group mean before the exclusion
        local ins = ("`o_insupport'" != "")
        matrix `BO' = J(`nout', 10 + `ins', .)
        forvalues j = 1/`nout' {
            if `ins' matrix `BO'[`j', 11] = `m1all_`j''
            local sfx = cond(`j' == 1, "", "_`j'")
            local c 0
            foreach q in m1 m1r m2 diff expl unexpl pct {
                matrix `BO'[`j', `++c'] = ``q'`sfx''
            }
            if `bsreps' > 0 {
                matrix `BO'[`j', 8]  = `pct_se`sfx''
                matrix `BO'[`j', 9]  = `pct_ll`sfx''
                matrix `BO'[`j', 10] = `pct_ul`sfx''
            }
        }
        matrix rownames `BO' = `ylist'
        matrix colnames `BO' = mean_focal mean_raked_focal mean_ref difference explained unexplained ///
            pct_explained pct_explained_se pct_explained_ll pct_explained_ul `=cond(`ins', "mean_focal_all", "")'
        ereturn matrix by_outcome = `BO'
    }
    if "`ratio'" != "" ereturn matrix ratio_parts = `RP'
    if "`covariatedetail'" != "" ereturn matrix covariates = `COV'
    if `bsreps' > 0 {
        ereturn scalar N_reps = `N_reps'
        ereturn scalar N_misreps = `N_fail'
        ereturn scalar N_iterlimit = `N_iter'
        if "`bsbca'" != "" ereturn scalar N_jkfail = `N_jkfail'
        if "`o_force'`ch_o_force'" != "" ereturn scalar N_forced = `N_force'
        if `N_fail' > 0 ereturn local missing missing
        ereturn scalar pct_explained_se = `pct_se'
        tempname PCI
        matrix `PCI' = (`pct_ll', `pct_ul')
        matrix colnames `PCI' = ll ul
        ereturn matrix pct_explained_ci = `PCI'
        ereturn scalar seed = `bsseed'
        ereturn scalar bca  = ("`bsbca'" != "")
        ereturn scalar ties = ("`bsties'" != "")
        if "`bsclust'" != "" ereturn scalar N_clust = `N_clust'
        ereturn local rngstate `"`rngstate'"'
        ereturn local clustvar `bsclust'
        ereturn local vce bayesboot
        ereturn local vcetype Bootstrap
        ereturn local bsnames `bsnames'
        // the interval types displayed, and the t intervals' degrees of freedom: e(df_r) too
        // when t is the type of the coefficient table's tests (the first of normal and t
        // displayed, else the default), so that test, lincom and the graphs use it
        ereturn local citypes `bstypes'
        ereturn scalar df_bs = `df_bs'
        local tnorm normal t
        local tfam : list bstypes & tnorm
        local tfam : word 1 of `tfam'
        if "`tfam'" == "" local tfam = cond("`bsclust'" != "", "t", "normal")
        if "`tfam'" == "t" ereturn scalar df_r = `df_bs'
        foreach m in b_bs bias se ci_normal ci_t ci_percentile ci_bc ci_bca z0 accel reps bsdraws bsdraws_pct {
            capture confirm matrix `BS'_`m'
            if _rc == 0 {
                ereturn matrix `m' = `BS'_`m'
            }
        }
        if `"`bsfile'"' != "" ereturn local bsfile `"`bsfile'"'
        if "`bsgen'" != "" ereturn local bsgenerate `bsgen'
    }
    if `RES'[`nstages', 21] > 0 {
        ereturn scalar explained_support = `m1' - `RES'[`nstages', 30]
        ereturn scalar explained_cs = `RES'[`nstages', 30] - `m1r'
    }
    ereturn scalar mincell = `mincell'
    if "`prefail'" != "" ereturn scalar N_prefail = `prefail'
    if `RES'[`nstages', 7] > 0 {
        ereturn local worstvar `worstvar`nstages''
        ereturn local worstcat `"`worstcat`nstages''"'
    }
    // e(stagetype): sequentially or individually when staged; e(fixothermargins) when in effect
    if "`o_variant'" != "" ereturn local stagetype `o_variant'
    else if "`sequential'" != "" ereturn local stagetype sequentially
    if "`holdothers'" != "" ereturn local fixothermargins fixothermargins
    if "`o_insupport'" != "" {
        ereturn local insupport insupport
        ereturn scalar N_focal_all = `N1all'
        ereturn scalar W_focal_all = `W1all'
        ereturn scalar mean_focal_all = `m1all_1'
    }
    local ncv : word count `cutvars'
    forvalues c = 1/`ncv' {
        local cmn : word `c' of `cutmats'
        ereturn matrix cuts`c' = `cmn'
    }
    if `ncv' > 0 {
        ereturn local contsample `contsample'
        ereturn local cutvars `cutvars'
        ereturn local cutq `cutq'
        ereturn local cutbins `cutnb'
    }
    if "`recut'" != "" ereturn local recut recut
    ereturn local generate `generate'
    local ipfrec `ipfopts'
    if `trimhiabs' != -1 local ipfrec `ipfrec' trimhiabs(`trimhiabs')
    if `trimloabs' != -1 local ipfrec `ipfrec' trimloabs(`trimloabs')
    if "`o_tso'" != "" local ipfrec `ipfrec' trimsupportoverride
    ereturn local rakeopts `ipfrec'
    ereturn local contvars `"`contvars'"'
    ereturn scalar N_blocks = `nblocks'
    forvalues s = 1/`nstages' {
        ereturn local stage`s' `stg`s'desc'
    }
    forvalues b = 1/`nblocks' {
        ereturn local block`b' `blk`b'desc'
    }
    forvalues t = 1/`nterms' {
        local ll
        forvalues l = 1/`nlev`t'' {
            if `l' > 1 {
                local ll `"`macval(ll)' | `macval(lev`t'_`l')'"'
            }
            else {
                local ll `"`macval(lev`t'_`l')'"'
            }
        }
        ereturn local levels`t' `"`macval(ll)'"'
        ereturn local term`t' `term`t'name'
    }
    ereturn scalar N_terms = `nterms'
    local allterms
    forvalues t = 1/`nterms' {
        local allterms `allterms' `term`t'name'
    }
    ereturn local rakevars `allterms'
    local blknames
    forvalues b = 1/`nblocks' {
        local blknames `blknames' `blk`b'name'
    }
    ereturn local blocks `blknames'
    ereturn local stagenames `stgnames'
    ereturn local rakespec `"`rakespec'"'
    ereturn local wexp `"`o_exp'"'
    ereturn local wtype `o_weight'
    ereturn local ref_group `"`macval(gref)'"'
    ereturn local focal_group `"`macval(gfocal)'"'
    ereturn local ref_group_label `"`macval(greflbl)'"'
    ereturn local focal_group_label `"`macval(gfocallbl)'"'
    ereturn local by `by'
    ereturn local depvar `ylist'
    if `nout' > 1 {
        ereturn local outcomes `ylist'
        ereturn scalar N_outcomes = `nout'
    }
    if "`ratio'" != "" ereturn local ratio `ratio'
    if "`reverse'" != "" ereturn local reverse reverse
    if "`detail'" != "" ereturn local detail detail
    ereturn local title "Raking decomposition (Pylypchuk-Selden)"
    ereturn scalar N_notes = `nnote'
    forvalues m = 1/`nnote' {
        ereturn local note`m' `"`macval(note`m')'"'
    }
    if `bsreps' > 0 ereturn local properties "b V"
    else ereturn local properties "b"
    ereturn local predict _no_predict
    // the full-sample call of an svy prefix: the prefix runs this command line after its
    // replicate loop, once its own results are posted (Prefix_epilog; the hook that mixed
    // uses, in brr, jackknife, sdr and svy bootstrap of Stata 18 and 19)
    if `ccrep' & !`ccloop' & "`: word 1 of `ccpfx''" == "svy" {
        ereturn hidden local prefix_epilog "seldecomp, prefixepilog"
    }
    // the full-sample call under a replication prefix: the display options typed, for the
    // display that bootstrap: and jackknife: run after their loop (the replay above)
    if `ccrep' & !`ccloop' {
        local pdo `diagnostics' `covariatedetail' `varlabels' `cformat' `pformat' `sformat' `citable'
        if "`pdo'" != "" ereturn hidden local prefix_dispopts `pdo'
    }
    ereturn local cmd seldecomp
    // the cut points of a full-sample call under a replication prefix, for its replicates; a
    // call without a prefix drops any record, so it can never be read by a later call
    if !`ccloop' {
        if `ccrep' {
            if `ncont' > 0 & "`recut'" == "" mata: seldecomp_cutcache("put", st_global("e(cutvars)"), "")
        }
        else mata: seldecomp_cutcache("reset", "", "")
    }

    c_local dispopts `diagnostics' `covariatedetail' `varlabels' `cformat' ///
        `pformat' `sformat' `citable'
    c_local makeplot `makeplot'
    c_local makeplot2 `"`macval(makeplot2)'"'
end


* ============================================================================
* Chain: chain, the complete decomposition of each pair of adjacent values of a
*   numeric by(), the lower one Focal Group, and their total. Estimate runs once as a probe
*   (chainprobe(): the options checked, the estimation sample of the whole chain marked and
*   its exclusions noted) and once per pair in pair mode (chainpair(): the ordinary code path
*   on the pair's two values, with the pair's coded terms left in variables named here); then
*   the svy pre-flight check over the pairs, Bootstrap over the pairs (one replicate loop),
*   and the joint ereturn. The user's results are restored if anything fails
* ============================================================================
program define Chain, eclass
    version 16
    // bayesboot() and makeplot() typed with empty parentheses, as in Estimate
    foreach a in bayes bayesb bayesbo bayesboo bayesboot {
        local 0 : subinstr local 0 " `a'()" " bayesboot", all
        local 0 : subinstr local 0 ",`a'()" ", bayesboot", all
    }
    foreach a in make makep makepl makeplo makeplot {
        local 0 : subinstr local 0 " `a'()" " makeplot", all
        local 0 : subinstr local 0 ",`a'()" ", makeplot", all
    }
    local cmd0 : copy local 0
    // the caller's sort order (the pre-flight check below reads it)
    local ccsort : sortedby
    syntax anything(equalok) [if] [in] [aw pw iw fw] , BY(varname) [ * ]
    local wgt
    if "`weight'" != "" {
        local wgt `"[`weight'`macval(exp)']"'
    }
    tempname ehold
    _estimates hold `ehold', restore nullok

    // ---- the probe: the options, the estimation sample of the chain ---------
    tempvar CT
    Estimate `macval(cmd0)' chainprobe(`CT')
    tempname LV
    mata: st_matrix("`LV'", uniqrows(st_data(., "`by'", "`CT'")))
    local K = rowsof(`LV')
    local J = `K' - 1
    // the values in ascending order, and their texts as Estimate makes them (value (label))
    local vl : value label `by'
    local blevs
    forvalues k = 1/`K' {
        local v`k' = `LV'[`k', 1]
        local blevs `blevs' `v`k''
        local lb`k' `v`k''
        if "`vl'" != "" & `v`k'' == int(`v`k'') {
            local l : label `vl' `v`k'', strict
            if `"`macval(l)'"' != "" {
                local lb`k' `"`v`k'' (`macval(l)')"'
            }
        }
        local lbx : subinstr local lb`k' "|" "/", all
        if `k' == 1 {
            local blabs `"`macval(lbx)'"'
        }
        else {
            local blabs `"`macval(blabs)' | `macval(lbx)'"'
        }
    }
    local nterms `ch_nterms'
    local ccloop `ch_ccloop'

    // ---- the pairs ------------------------------------------------------------
    // pair j: the values j and j + 1; Estimate in pair mode leaves its sample, Focal Group flag and
    // coded terms in PTj, PGj and TVj_#, and its results in e(), read here before the next pair
    mata: seldecomp_chainnotes("reset", 0)
    tempname BP PR T b
    local ptl
    local pgl
    local ptvl
    local pnll
    local cmin 1
    local tmin 1
    local nexcl 0
    forvalues j = 1/`J' {
        tempvar PV`j' PT`j' PG`j'
        tempname B`j' S`j' RP`j' CK`j' MI`j'
        qui gen byte `PV`j'' = `CT' & (`by' == `LV'[`j', 1] | `by' == `LV'[`j' + 1, 1])
        local tvn`j'
        forvalues t = 1/`nterms' {
            tempvar TV`j'_`t'
            local tvn`j' `tvn`j'' `TV`j'_`t''
        }
        // a pair that fails stops the chain with its own message; the pairs' notes collected
        // so far (a Mata external) are dropped first. Its output as without capture: noisily
        // alone would show it under quietly, which shows only the error messages
        if c(noisily) {
            capture noisily Estimate `macval(anything)' if `PV`j'' `macval(wgt)', by(`by') ///
                `macval(options)' chainpair(`j') chainvars(`PT`j'' `PG`j'' `tvn`j'')
        }
        else {
            capture noisily quietly Estimate `macval(anything)' if `PV`j'' `macval(wgt)', ///
                by(`by') `macval(options)' chainpair(`j') chainvars(`PT`j'' `PG`j'' `tvn`j'')
        }
        if _rc {
            local rc = _rc
            mata: seldecomp_chainnotes("reset", 0)
            exit `rc'
        }
        local ptl `ptl' `PT`j''
        local pgl `pgl' `PG`j''
        local ptvl `ptvl' `tvn`j''
        local nl`j' `ch_nlevlist'
        local pnll `pnll' `ch_nlevlist'
        matrix `B`j'' = e(b)
        matrix `S`j'' = e(bystage)
        matrix `CK`j'' = e(raking_checks)
        local nm`j' 0
        capture matrix `MI`j'' = e(raking_misses)
        if !_rc {
            local nm`j' = rowsof(`MI`j'')
            forvalues m = 1/`nm`j'' {
                mata: st_local("rmv`j'_`m'", st_global("e(raking_miss`m'_var)"))
                mata: st_local("rml`j'_`m'", st_global("e(raking_miss`m'_level)"))
            }
        }
        if "`e(ratio)'" != "" matrix `RP`j'' = e(ratio_parts)
        // e(by_pair): the pair's means, parts, share, counts and weights (and insupport's)
        local ins = ("`e(insupport)'" != "")
        if `j' == 1 matrix `BP' = J(`J', 11 + 2 * `ins', .)
        local c 0
        foreach s in mean_focal mean_raked_focal mean_ref difference explained unexplained pct_explained ///
            N_focal N_ref W_focal W_ref {
            matrix `BP'[`j', `++c'] = e(`s')
        }
        if `ins' {
            matrix `BP'[`j', 12] = e(N_focal_all)
            matrix `BP'[`j', 13] = e(mean_focal_all)
            local nexcl = `nexcl' + e(N_focal_all) - e(N_focal)
        }
        if e(converged) < `cmin' local cmin = e(converged)
        if e(targets_met) < `tmin' local tmin = e(targets_met)
        // the cut points of the pair's c. variables
        local ncv : word count `e(cutvars)'
        forvalues c = 1/`ncv' {
            tempname C`c'_`j'
            matrix `C`c'_`j'' = e(cuts`c')
        }
        // the labels of the levels of the pair's terms, for the pre-flight check
        forvalues t = 1/`nterms' {
            mata: st_local("lv`j'_`t'", st_global("e(levels`t')"))
        }
        mata: seldecomp_chainnotes("add", `j')
        if `j' == `J' {
            // what does not depend on the pair
            local nstages `ch_nstages'
            local stagedefs `ch_stagedefs'
            local tnames `ch_tnames'
            local engopts `ch_engopts'
            local ipfopts `ch_ipfopts'
            local holdothers `ch_holdothers'
            foreach m in depvar blocks stagenames rakevars rakespec stagetype fixothermargins ///
                diffsign wtype wexp rakeopts ratio insupport contvars cutvars cutq       ///
                contsample recut title {
                mata: st_local("em_`m'", st_global("e(`m')"))
            }
            foreach m in N_stages N_blocks N_terms mincell {
                local es_`m' = e(`m')
            }
            forvalues s = 1/`es_N_stages' {
                local em_stage`s' `e(stage`s')'
            }
            forvalues s = 1/`es_N_blocks' {
                local em_block`s' `e(block`s')'
            }
            forvalues t = 1/`es_N_terms' {
                local em_term`t' `e(term`t')'
            }
        }
    }
    // the estimation sample: the union of the pairs' samples (with insupport, the first
    // value's observations that pair 1 excludes are in no pair)
    tempvar CU
    qui gen byte `CU' = 0
    foreach pt of local ptl {
        qui replace `CU' = 1 if `pt'
    }
    local staged = ("`em_stagetype'" != "")
    local nb `es_N_blocks'
    local prn
    forvalues j = 1/`J' {
        local prn `prn' p`j'
    }
    matrix rownames `BP' = `prn'
    local bpc mean_focal mean_raked_focal mean_ref difference explained unexplained pct_explained N_focal N_ref W_focal W_ref
    if `ins' local bpc `bpc' N_focal_all mean_focal_all
    matrix colnames `BP' = `bpc'
    matrix `PR' = `LV'[1..`J', 1], `LV'[2..`K', 1]
    matrix rownames `PR' = `prn'
    matrix colnames `PR' = focal_group ref_group

    // ---- e(b): each pair's equations prefixed pj, then the total -------------
    local bn
    forvalues j = 1/`J' {
        local cn : colfullnames `B`j''
        foreach c of local cn {
            gettoken eq nm : c, parse(":")
            gettoken colon nm : nm, parse(":")
            if "`eq'" == "overall" local bn `bn' p`j':`nm'
            else local bn `bn' p`j'_`eq':`nm'
        }
        if `j' == 1 matrix `b' = `B`j''
        else matrix `b' = (`b', `B`j'')
    }
    // the total: the difference, explained and unexplained parts and each block's contribution,
    // summed over the pairs in their order
    local kt = 3 + `nb' * `staged'
    matrix `T' = J(1, `kt', 0)
    forvalues j = 1/`J' {
        forvalues c = 1/3 {
            matrix `T'[1, `c'] = `T'[1, `c'] + `B`j''[1, 3 + `c']
        }
        forvalues s = 1/`=`nb' * `staged'' {
            matrix `T'[1, 3 + `s'] = `T'[1, 3 + `s'] + `B`j''[1, 6 + `s']
        }
    }
    local bn `bn' total:difference total:explained total:unexplained
    if `staged' {
        foreach bl of local em_stagenames {
            local bn `bn' total_marg_expl:`bl'
        }
    }
    matrix `b' = (`b', `T')
    matrix colnames `b' = `bn'
    local tpct = cond(`T'[1, 1] != 0, 100 * `T'[1, 2] / `T'[1, 1], .)

    // ---- the notes: the chain's exclusions, then the pairs' (merged) ------------
    local nnote 0
    forvalues m = 1/`ch_nnote' {
        local ++nnote
        local note`nnote' `"`macval(ch_note`m')'"'
    }
    mata: seldecomp_chainnotes("get", 0)
    forvalues m = 1/`cnote_n' {
        local ++nnote
        local note`nnote' `"`macval(cnote`m')'"'
    }
    if `"`ch_svynote'"' != "" {
        local ++nnote
        local note`nnote' `"`ch_svynote'"'
    }

    // ---- the pre-flight check of an svy prefix's replicates, over the pairs ----
    // (Estimate's, pair by pair: seldecomp_prefail on the chain's replicate weights or PSUs;
    // a replicate that fails in two pairs counts once; one warning). Not with recut, whose
    // binned terms are cut again in each replicate
    local prefail
    if `ch_ccrep' & !`ccloop' & "`: word 1 of `ch_ccpfx''" == "svy" & "`ch_recut'" == "" {
        local pfm : word 2 of `ch_ccpfx'
        local pfrw
        local pfpsu
        if inlist("`pfm'", "brr", "jackknife", "sdr", "bootstrap") {
            local pfc = cond("`pfm'" == "jackknife", "jkr", cond("`pfm'" == "bootstrap", "bsr", "`pfm'"))
            local pfrw : char _dta[_svy_`pfc'weight]
            if "`pfrw'" != "" {
                capture confirm numeric variable `pfrw', exact
                if _rc local pfrw
                else if "`pfm'" == "brr" {
                    capture mata: seldecomp_hadroute("`: char _dta[_svy_strata1]'", ///
                        "`: char _dta[_svy_su1]'", "`ccsort'")
                    if !_rc & "`pf_had'" == "1" local pfrw
                }
            }
            else if "`pfm'" == "jackknife" {
                local pfsu : char _dta[_svy_su1]
                local pfst : char _dta[_svy_strata1]
                tempvar pfpsu
                if inlist("`pfsu'", "", "_n") capture gen double `pfpsu' = _n if `CT'
                else capture egen double `pfpsu' = group(`pfst' `pfsu') if `CT'
                if _rc local pfpsu
            }
        }
        if "`pfrw'`pfpsu'" != "" {
            local fids
            local vids
            local pfitems
            local pfok 1
            // the replicates in which value k of by() has no weight (cv1 ... cvK): pair j
            // without Focal Group (value j) or without Reference Group (value j + 1)
            forvalues k = 1/`K' {
                local cv`k'
            }
            forvalues j = 1/`J' {
                capture mata: seldecomp_prefail("`tvn`j''", "`nl`j''", "`PT`j''", "`PG`j''", ///
                    "`pfrw'", "`pfpsu'", 1)
                if _rc {
                    local pfok 0
                    continue, break
                }
                local fids : list fids | pf_fids
                local vids : list vids | pf_vids
                local cv`j' : list cv`j' | pf_g1ids
                local k = `j' + 1
                local cv`k' : list cv`k' | pf_g2ids
                foreach k in fail van {
                    foreach ff of local pf_`k' {
                        tokenize "`ff'", parse(":")
                        local key `j'_`1'_`3'
                        if !`: list key in pfitems' local pfitems `pfitems' `key'
                        local pfn`k'_`key' `5'
                    }
                }
            }
            if `pfok' {
                // the replicates svy runs: with the PSU (or observation) jackknife not known here
                // (under subpop(), it has the PSUs of the whole strata), so the count alone is
                // given; svy jackknife skips a replicate-weight variable equal to the full-sample
                // weight on the estimation sample
                if "`pfpsu'" != "" local pfR
                else if "`pfm'" == "jackknife" {
                    tempvar pfw
                    if "`weight'" == "" qui gen double `pfw' = 1
                    else qui gen double `pfw' `macval(exp)' if `CU'
                    mata: st_local("pfR", strofreal(sum(colsum(st_data(., tokens("`pfrw'"), "`CU'") ///
                        :!= st_data(., "`pfw'", "`CU'")) :> 0)))
                }
                else local pfR : word count `pfrw'
                local vids : list vids - fids
                local nf : word count `fids'
                local nv : word count `vids'
                local prefail `nf'
                local nkeys : word count `pfitems'
                // the values of by() without weight in some replicates (at most 3 named)
                local cvids
                local ncvv 0
                forvalues k = 1/`K' {
                    if "`cv`k''" == "" continue
                    local cvids : list cvids | cv`k'
                    local ++ncvv
                    if `ncvv' <= 3 local cvl`ncvv' "`by' = `v`k''"
                }
                local ncvr : word count `cvids'
                if `ncvr' > 0 {
                    local cvw `"`macval(cvl1)'"'
                    if `ncvv' == 2 {
                        local cvw `"`macval(cvw)' or `macval(cvl2)'"'
                    }
                    else if `ncvv' == 3 {
                        local cvw `"`macval(cvw)', `macval(cvl2)' or `macval(cvl3)'"'
                    }
                    else if `ncvv' > 3 {
                        local cvw `"`macval(cvw)', `macval(cvl2)', `macval(cvl3)' or `=`ncvv' - 3' more values"'
                    }
                }
                if `nkeys' > 0 | `ncvr' > 0 {
                    local pfwhat = cond("`pfm'" == "brr", "BRR", cond("`pfm'" == "sdr", "SDR", "`pfm'"))
                    local pffay = cond("`pfm'" == "brr", ", or use Fay's BRR (fay())", "")
                    // the text of each (pair, term, level)
                    foreach key of local pfitems {
                        tokenize "`key'", parse("_")
                        local pj `1'
                        local pt `3'
                        local pl `5'
                        mata: st_local("pflab", strtrim(ustrsplit(st_local("lv`pj'_`pt'"), "[|]")[`pl']))
                        local pflev`key' `"p`pj': `em_term`pt'' = `macval(pflab)'"'
                    }
                    local key : word 1 of `pfitems'
                    // "4 of 8 BRR replicates", or with the PSU jackknife "1 jackknife replicate"
                    if "`pfR'" != "" {
                        local pfnf "`nf' of `pfR' `pfwhat' replicates"
                        local pfnv "`nv' of `pfR' `pfwhat' replicates give"
                    }
                    else {
                        local pfnf "`nf' jackknife `=plural(`nf', "replicate")'"
                        local pfnv "`nv' jackknife `=cond(`nv' == 1, "replicate gives", "replicates give")'"
                    }
                    if `nkeys' == 0 {
                        // only values of by() without weight: they are the cause
                        local pfl `"`pfnf' will fail: in `=cond(`nf' == 1, "it", "them")', `macval(cvw)' has no weight."'
                        if "`pfm'" == "brr" {
                            local pfl `"`macval(pfl)' Use Fay's BRR (fay())."'
                        }
                        local pflines 0
                    }
                    else if `nkeys' == 1 & `nv' == 0 & "`pfnfail_`key''" != "" & `ncvr' == 0 {
                        local pfl `"`macval(pflev`key')' has no focal group weight in `pfnf', which will fail. Collapse that level`pffay'."'
                        local pflines 0
                    }
                    else {
                        local pfl
                        if `nf' > 0 local pfl "`pfnf' will fail"
                        if `ncvr' > 0 {
                            local pfl `"`macval(pfl)' (in `ncvr', `macval(cvw)' has no weight)"'
                        }
                        if `nv' > 0 {
                            if `nf' > 0 {
                                local pfl `"`macval(pfl)', and `nv' `=cond(`nv' == 1, "other gives", "others give")' focal group obs zero weight"'
                            }
                            else local pfl "`pfnv' focal group obs zero weight"
                        }
                        local pfcl = cond(`nkeys' == 1, "this level", "these levels")
                        local pfl `"`macval(pfl)'; collapse `pfcl'`pffay':"'
                        local pflines 1
                    }
                    // user text (labels): SMCL braces escaped
                    local pfl : subinstr local pfl "{" "`=char(1)'", all
                    local pfl : subinstr local pfl "}" "{c )-}", all
                    local pfl : subinstr local pfl "`=char(1)'" "{c -(}", all
                    di as err `"{p 0 4 2 78}warning: `macval(pfl)'{p_end}"'
                    if `pflines' {
                        local nsh 0
                        foreach key of local pfitems {
                            local ++nsh
                            if `nsh' > 10 {
                                di as err "{p 8 12 2 78}... and `=`nkeys' - 10' more{p_end}"
                                continue, break
                            }
                            local pfl `"`macval(pflev`key')':"'
                            if "`pfnfail_`key''" != "" {
                                local pfl `"`macval(pfl)' no focal group weight in `pfnfail_`key'' `=plural(`pfnfail_`key'', "replicate")'"'
                            }
                            if "`pfnvan_`key''" != "" {
                                if "`pfnfail_`key''" != "" {
                                    local pfl `"`macval(pfl)',"'
                                }
                                local pfl `"`macval(pfl)' no reference group weight in `pfnvan_`key'' `=plural(`pfnvan_`key'', "replicate")'"'
                            }
                            local pfl : subinstr local pfl "{" "`=char(1)'", all
                            local pfl : subinstr local pfl "}" "{c )-}", all
                            local pfl : subinstr local pfl "`=char(1)'" "{c -(}", all
                            di as err `"{p 8 12 2 78}`macval(pfl)'{p_end}"'
                        }
                    }
                }
            }
        }
    }

    // ---- bayesboot(): one replicate loop over the pairs ---------------------
    local bsreps `ch_bsreps'
    if `bsreps' > 0 {
        tempvar W ONE
        if "`weight'" == "" qui gen double `W' = 1
        else qui gen double `W' `macval(exp)' if `CT'
        qui gen byte `ONE' = 1
        tempname V BS
        local bsopts
        if "`ch_bsclust'" != "" local bsopts `bsopts' cluster(`ch_bsclust')
        if `"`ch_bsseed'"' != "" local bsopts `bsopts' seed(`ch_bsseed')
        if "`ch_bsgen'" != "" local bsopts `bsopts' bsgen(`ch_bsgen')
        if "`ch_bsgenreplace'" != "" local bsopts `bsopts' bsgenreplace
        if `"`ch_bsfile'"' != "" local bsopts `bsopts' bsfile(`"`ch_bsfile'"') `ch_bsfilereplace'
        if "`ch_o_force'" != "" local bsopts `bsopts' force
        // the observed shares: each pair's, then the total's
        local pctlist
        forvalues j = 1/`J' {
            local pctlist `pctlist' `=`BP'[`j', 7]'
        }
        local pctlist `pctlist' `tpct'
        local y1 : word 1 of `em_depvar'
        // the draws on the chain's estimation sample CU (with insupport, not on the first
        // value's observations that no pair uses)
        Bootstrap, touse(`CU') g1(`PG1') w(`W') y(`y1') one(`ONE') nterms(`nterms')          ///
            tvlist(`tvn1') nlevlist(`nl1') nstages(`nstages') stagedefs(`stagedefs')          ///
            b(`b') v(`V') bs(`BS') reps(`bsreps') dots(`ch_bsdots') `ch_bsnodots'             ///
            ipfopts(`ipfopts') tnames(`tnames') level(`ch_level') engopts(`engopts')          ///
            pct(`pctlist') `holdothers' `ch_bsseq' `ch_bsbca' `ch_bsties' `bsopts'            ///
            chain(`J') ptouse(`ptl') pg1(`pgl') ptv(`ptvl') pnl(`pnll')
        local N_reps    = r(N_reps)
        local N_fail    = r(N_fail)
        local N_iter    = r(N_iter)
        local N_force   = r(N_force)
        local N_jkforce = r(N_jkforce)
        local N_jkfail  = r(N_jkfail)
        local N_clust   = r(N_clust)
        local df_bs     = r(df_bs)
        local bsseed    = r(seed)
        local rngstate  `"`r(rngstate)'"'
        local bsnames   `r(bsnames)'
        forvalues j = 1/`K' {
            local pct_se_`j' = r(pct_se_`j')
            local pct_ll_`j' = r(pct_ll_`j')
            local pct_ul_`j' = r(pct_ul_`j')
        }
        // failed replications and those kept at the iterate() limit: one note (the point
        // estimate's notes are the pairs'; they are not folded together here)
        local stoptxt "at the iterate() limit"
            local nstop = r(N_failit) + `N_iter'
        if `N_fail' > 0 & `N_fail' == r(N_failit) & `N_iter' > 0 {
            di as txt "{p 0 6 2 78}note: raking stopped `stoptxt' in `nstop' of `bsreps' replications; " ///
                "`N_fail' of them missed their targets and " cond(`N_fail' == 1, "was", "were") ///
                " dropped; raise iterate(); bootstrap results may be incorrect.{p_end}"
        }
        else if `N_fail' > 0 {
            di as txt "{p 0 6 2 78}note: `N_fail' of `bsreps' replications failed and " ///
                cond(`N_fail' == 1, "was", "were") " dropped: in " cond(`N_fail' == 1, "it", "them") ///
                `", `r(bswhy)'`r(bsfix)'; bootstrap results may be incorrect.{p_end}"'
        }
        else if `N_iter' > 0 {
            di as txt "{p 0 6 2 78}note: `N_iter' of `bsreps' replications stopped `stoptxt' with the " ///
                "raking targets met and " cond(`N_iter' == 1, "was", "were") " kept; raise iterate(){p_end}"
        }
        // force: the replications (and BCa jackknife deletions) kept although a margin was missed
        if `N_force' + `N_jkforce' > 0 {
            local jkt = cond(`N_jkforce' > 0, " and `N_jkforce' jackknife " + plural(`N_jkforce', "deletion"), "")
            di as txt "{p 0 6 2 78}note: force kept `N_force' of `bsreps' " plural(`bsreps', "replication") ///
                "`jkt' that missed raking targets{p_end}"
        }
        // e(by_pair): the pairs' shares' std. err. and percentile CIs
        matrix `BP' = `BP', J(`J', 3, .)
        forvalues j = 1/`J' {
            matrix `BP'[`j', colsof(`BP') - 2] = `pct_se_`j''
            matrix `BP'[`j', colsof(`BP') - 1] = `pct_ll_`j''
            matrix `BP'[`j', colsof(`BP')] = `pct_ul_`j''
        }
        matrix colnames `BP' = `bpc' pct_explained_se pct_explained_ll pct_explained_ul
    }

    // ---- post ---------------------------------------------------------------
    // the estimation sample: CU, the union of the pairs' samples (above)
    qui count if `CU'
    local N = r(N)
    _estimates unhold `ehold', not
    if `bsreps' > 0 {
        ereturn post `b' `V', depname(`y1') obs(`N') esample(`CU')
    }
    else {
        local y1 : word 1 of `em_depvar'
        ereturn post `b', depname(`y1') obs(`N') esample(`CU')
    }
    ereturn local chain chain
    ereturn scalar N = `N'
    ereturn scalar N_pairs = `J'
    ereturn scalar N_zerowt = `ch_nzw'
    ereturn scalar mean_focal = el(`BP', 1, 1)
    ereturn scalar mean_ref = el(`BP', `J', 3)
    ereturn scalar difference = el(`T', 1, 1)
    ereturn scalar explained = el(`T', 1, 2)
    ereturn scalar unexplained = el(`T', 1, 3)
    ereturn scalar pct_explained = `tpct'
    ereturn scalar N_stages = `es_N_stages'
    ereturn scalar k_eq = `J' * (1 + 4 * `staged') + 1 + `staged'
    ereturn scalar level = `ch_level'
    ereturn scalar converged = `cmin'
    ereturn scalar targets_met = `tmin'
    ereturn scalar mincell = `es_mincell'
    ereturn scalar N_blocks = `es_N_blocks'
    ereturn scalar N_terms = `es_N_terms'
    if "`prefail'" != "" ereturn scalar N_prefail = `prefail'
    ereturn matrix pairs = `PR'
    ereturn matrix by_pair = `BP'
    forvalues j = 1/`J' {
        if `j' == 1 ereturn matrix bystage = `S`j''
        else ereturn matrix bystage`j' = `S`j''
        local js = cond(`j' == 1, "", "`j'")
        ereturn matrix raking_checks`js' = `CK`j''
        if `nm`j'' > 0 {
            ereturn matrix raking_misses`js' = `MI`j''
            local mp = cond(`j' == 1, "raking_miss", "raking_miss`j'_")
            forvalues m = 1/`nm`j'' {
                ereturn local `mp'`m'_var `"`macval(rmv`j'_`m')'"'
                ereturn local `mp'`m'_level `"`macval(rml`j'_`m')'"'
            }
        }
    }
    if "`em_ratio'" != "" {
        tempname RPA
        local rpn
        forvalues j = 1/`J' {
            if `j' == 1 matrix `RPA' = `RP`j''
            else matrix `RPA' = `RPA' \ `RP`j''
            local rpn `rpn' p`j':focal_group p`j':raked_focal p`j':ref_group
        }
        matrix rownames `RPA' = `rpn'
        ereturn matrix ratio_parts = `RPA'
    }
    forvalues j = 1/`J' {
        forvalues c = 1/`ncv' {
            ereturn matrix cuts`c'_p`j' = `C`c'_`j''
        }
    }
    if `bsreps' > 0 {
        ereturn scalar N_reps = `N_reps'
        ereturn scalar N_misreps = `N_fail'
        ereturn scalar N_iterlimit = `N_iter'
        if "`ch_bsbca'" != "" ereturn scalar N_jkfail = `N_jkfail'
        if "`o_force'`ch_o_force'" != "" ereturn scalar N_forced = `N_force'
        if `N_fail' > 0 ereturn local missing missing
        ereturn scalar pct_explained_se = `pct_se_`K''
        tempname PCI
        matrix `PCI' = (`pct_ll_`K'', `pct_ul_`K'')
        matrix colnames `PCI' = ll ul
        ereturn matrix pct_explained_ci = `PCI'
        ereturn scalar seed = `bsseed'
        ereturn scalar bca = ("`ch_bsbca'" != "")
        ereturn scalar ties = ("`ch_bsties'" != "")
        if "`ch_bsclust'" != "" ereturn scalar N_clust = `N_clust'
        ereturn local rngstate `"`rngstate'"'
        ereturn local clustvar `ch_bsclust'
        ereturn local vce bayesboot
        ereturn local vcetype Bootstrap
        ereturn local bsnames `bsnames'
        // the interval types displayed, and the t intervals' degrees of freedom: e(df_r) too
        // when t is the type of the coefficient table's tests (the first of normal and t
        // displayed, else the default), so that test, lincom and the graphs use it
        ereturn local citypes `ch_bstypes'
        ereturn scalar df_bs = `df_bs'
        local tnorm normal t
        local tfam : list ch_bstypes & tnorm
        local tfam : word 1 of `tfam'
        if "`tfam'" == "" local tfam = cond("`ch_bsclust'" != "", "t", "normal")
        if "`tfam'" == "t" ereturn scalar df_r = `df_bs'
        foreach m in b_bs bias se ci_normal ci_t ci_percentile ci_bc ci_bca z0 accel reps bsdraws bsdraws_pct {
            capture confirm matrix `BS'_`m'
            if _rc == 0 {
                ereturn matrix `m' = `BS'_`m'
            }
        }
        if `"`ch_bsfile'"' != "" ereturn local bsfile `"`ch_bsfile'"'
        if "`ch_bsgen'" != "" ereturn local bsgenerate `ch_bsgen'
    }
    foreach m in cutvars cutq contsample {
        if "`em_`m''" != "" ereturn local `m' `em_`m''
    }
    if "`em_recut'" != "" ereturn local recut recut
    if "`em_insupport'" != "" ereturn local insupport insupport
    if "`em_ratio'" != "" ereturn local ratio `em_ratio'
    ereturn local stagetype `em_stagetype'
    ereturn local fixothermargins `em_fixothermargins'
    ereturn local rakeopts `em_rakeopts'
    ereturn local contvars `"`em_contvars'"'
    forvalues s = 1/`es_N_stages' {
        ereturn local stage`s' `em_stage`s''
    }
    forvalues bl = 1/`es_N_blocks' {
        ereturn local block`bl' `em_block`bl''
    }
    forvalues t = 1/`es_N_terms' {
        ereturn local term`t' `em_term`t''
    }
    ereturn local rakevars `em_rakevars'
    ereturn local blocks `em_blocks'
    ereturn local stagenames `em_stagenames'
    ereturn local rakespec `"`em_rakespec'"'
    ereturn local wexp `"`em_wexp'"'
    ereturn local wtype `em_wtype'
    ereturn local diffsign `em_diffsign'
    ereturn local by_levels `blevs'
    ereturn local by_labels `"`macval(blabs)'"'
    ereturn local ref_group `v`K''
    ereturn local focal_group `v1'
    ereturn local ref_group_label `"`macval(lb`K')'"'
    ereturn local focal_group_label `"`macval(lb1)'"'
    ereturn local by `by'
    ereturn local depvar `em_depvar'
    ereturn local title `"`em_title'"'
    ereturn scalar N_notes = `nnote'
    forvalues m = 1/`nnote' {
        ereturn local note`m' `"`macval(note`m')'"'
    }
    if `bsreps' > 0 ereturn local properties "b V"
    else ereturn local properties "b"
    ereturn local predict _no_predict
    if `ch_ccrep' & !`ccloop' & "`: word 1 of `ch_ccpfx''" == "svy" {
        ereturn hidden local prefix_epilog "seldecomp, prefixepilog"
    }
    if `ch_ccrep' & !`ccloop' & "`ch_dispopts'" != "" ereturn hidden local prefix_dispopts `ch_dispopts'
    ereturn local cmd seldecomp
    c_local dispopts `ch_dispopts'
    c_local makeplot `ch_makeplot'
    c_local makeplot2 `"`macval(ch_makeplot2)'"'
end


* ============================================================================
* Parse_term: one token of the raking list -> variables, kinds (i/c), name
*   returns (c_local) s_nnew, s_vars#, s_kinds#, s_name#
* ============================================================================
program define Parse_term
    args tok
    if strpos("`tok'", "(") | strpos("`tok'", ")") | strpos("`tok'", ",") {
        di as err "invalid term `tok' in the raking variable list: i.(varlist), nested " ///
            "parentheses, and commas are not allowed"
        exit 198
    }
    local t = subinstr("`tok'", "##", "#", .)
    if substr("`t'", 1, 1) == "#" | substr("`t'", -1, 1) == "#" | strpos("`t'", "##") {
        di as err "invalid interaction `tok' in the raking variable list"
        exit 198
    }
    local inter = strpos("`t'", "#") > 0
    local vars
    local kinds
    tokenize "`t'", parse("#")
    local i 1
    while "``i''" != "" {
        local piece ``i''
        if "`piece'" != "#" {
            if substr("`piece'", 1, 2) == "c." {
                local kind c
                local v = substr("`piece'", 3, .)
            }
            else if substr("`piece'", 1, 2) == "i." {
                local kind i
                local v = substr("`piece'", 3, .)
            }
            else if strpos("`piece'", ".") {
                di as err "{p 0 4 2}`piece': only the i. and c. operators are allowed in the raking variable list{p_end}"
                exit 198
            }
            else {
                local kind i
                local v `piece'
            }
            if "`v'" == "" {
                di as err "invalid term `tok' in the raking variable list"
                exit 198
            }
            capture unab vl : `v'
            if _rc {
                di as err "variable `v' not found"
                exit 111
            }
            if `inter' & `: word count `vl'' > 1 {
                di as err "{p 0 4 2}`tok': wildcards and variable lists are not allowed inside interaction terms{p_end}"
                exit 198
            }
            foreach vv of local vl {
                local vars `vars' `vv'
                local kinds `kinds' `kind'
            }
        }
        local ++i
    }
    if "`vars'" == "" {
        di as err "invalid term `tok' in the raking variable list"
        exit 198
    }
    if `inter' {
        local uvars : list uniq vars
        if `: list sizeof uvars' < `: list sizeof vars' {
            di as err "`tok': a variable may appear only once in an interaction"
            exit 198
        }
        c_local s_nnew 1
        c_local s_vars1 `vars'
        c_local s_kinds1 `kinds'
        local nm : subinstr local vars " " "#", all
        c_local s_name1 `nm'
    }
    else {
        local n : word count `vars'
        c_local s_nnew `n'
        forvalues k = 1/`n' {
            c_local s_vars`k' `: word `k' of `vars''
            c_local s_kinds`k' `: word `k' of `kinds''
            c_local s_name`k' `: word `k' of `vars''
        }
    }
end


* ============================================================================
* Subopt_error <option> <spec> <contents> [<pre>]: -syntax- rejected <contents> (the
*   caller's local 0 for the list inside option(); <pre> is what precedes its comma) with
*   the options <spec>; say which option list did so. An unknown suboption is 198; a known
*   one with a bad argument keeps -syntax-'s message and return code
* ============================================================================
program define Subopt_error
    args opt spec zero pre
    // every suboption declared as free text: only unknown names are left over in options
    local loose = ustrregexra(`"`spec'"', "\([^)]*\)", "(string asis)")
    local 0 `"`macval(zero)'"'
    capture syntax `pre' [, `loose' * ]
    if _rc == 0 & `"`macval(options)'"' != "" {
        gettoken bad : options, bind
        local nm = ustrregexra(`"`macval(bad)'"', "\(.*$", "()")
        local hint
        if "`opt'" == "bayesboot" & inlist(`"`macval(nm)'"', "seed", "seed()") {
            local hint "; use rseed()"
        }
        // strata() is not a suboption: the hint says where a stratified design goes (not the
        // bootstrap prefix with strata(): it refuses weights without force and, resampling n of n
        // PSUs, gives variances (n - 1)/n of the design-based ones)
        if "`opt'" == "bayesboot" & ustrregexm(`"`macval(nm)'"', "^str(a(t(a)?)?)?(\(\))?$") {
            local hint "; for a stratified design use svy jackknife: or svy brr:"
        }
        di as err `"{p 0 4 2}`opt'(): option `macval(nm)' not allowed`hint'{p_end}"'
        exit 198
    }
    local 0 `"`macval(zero)'"'
    capture noisily syntax `pre' [, `spec' ]
    local rc = cond(_rc, _rc, 198)
    di as err "(error in option `opt'())"
    exit `rc'
end



* ============================================================================
* Comma_check <spec> <command line>: a word directly followed by a parenthesis among the
*   raking variables (by(married), focalgroup(1)) is an option typed before the comma: say
*   so (error 198) and show the command with it after the comma. No valid raking term has
*   that form (a block begins with the parenthesis; Parse_term refuses any other). A plain
*   word that is no variable goes after the comma too when it follows such an option or is
*   one of seldecomp's options without argument in <spec> (sequentially, force). by() after an
*   if or in qualifier (which takes it) is found in the text before the comma, and the
*   option words before the qualifier are moved with it. Without an option typed before
*   the comma, option words that are no variables get "sequentially is an option: put it after
*   the comma" (198; the command would fail with 111 in Parse_term). The outcome (the
*   first word) is never moved; a valid command has no such word and passes unchanged
* ============================================================================
program define Comma_check
    args spec zero
    local 0 `"`macval(zero)'"'
    local bare
    foreach o of local spec {
        if !strpos("`o'", "(") local bare `bare' `o'
    }
    local first
    local moved
    capture syntax [anything(equalok)] [if] [in] [aw pw iw fw] [, * ]
    if !_rc {
        local keep
        local rest `"`macval(anything)'"'
        local k 0
        while `"`macval(rest)'"' != "" {
            gettoken tok rest : rest, bind
            local ++k
            mata: st_local("nm", ustrregexm(st_local("tok"), "^([A-Za-z_][A-Za-z0-9_]*)[(]") ? ustrregexs(1) : "")
            // a factor-variable operator with parentheses (ib(first).x, i(1 2).x) is no option:
            // Estimate refuses it as an operator
            mata: st_local("nm", ustrregexm(st_local("tok"), "^[A-Za-z_][A-Za-z0-9_]*[(][^()]*[)][.]") ? "" : st_local("nm"))
            // the options that the by prefix appends follow a command line without a comma
            if inlist("`nm'", "seldbyvar", "seldbyidx", "seldbyvars") continue
            local mv 0
            if "`nm'" != "" {
                if "`first'" == "" local first `nm'
                local mv 1
            }
            else {
                // the first word is the outcome: it stays where it is
                mata: st_local("pw", strofreal(ustrregexm(st_local("tok"), "^[A-Za-z_][A-Za-z0-9_]*$")))
                if `pw' & `k' > 1 {
                    capture confirm variable `tok'
                    if _rc {
                        mata: st_local("mv", strofreal(st_local("first") != "" | seldecomp_isopt(st_local("tok"), st_local("bare"))))
                    }
                }
            }
            if `mv' {
                local moved `"`macval(moved)' `macval(tok)'"'
            }
            else {
                local keep `"`macval(keep)' `macval(tok)'"'
            }
        }
        // the command with those words after the comma
        local wt
        if "`weight'" != "" local wt `"[`weight'`macval(exp)']"'
        local pre `"`macval(keep)' `macval(if)' `macval(in)' `macval(wt)'"'
        local post `"`macval(moved)' `macval(options)'"'
    }
    local also
    if "`first'" == "" {
        // by() typed after an if or in qualifier, which then takes it: by( in the text before
        // the first comma outside parentheses and quotes (seldecomp_commaby)
        mata: seldecomp_commaby(st_local("zero"))
        if "`cbfound'" != "1" {
            // no option typed before the comma, but seldecomp's options without argument
            // that are not variables (sequentially, force): they go after the comma, after the
            // options typed there (a command without a comma lacks by(): -syntax- says so)
            mata: st_local("opts", stritrim(strtrim(ustrregexra(st_local("options"), "seldby(var|idx|vars)[(][^)]*[)]", ""))))
            if `"`macval(moved)'"' == "" | `"`macval(opts)'"' == "" exit
            local post `"`macval(options)' `macval(moved)'"'
            mata: st_local("pre", stritrim(strtrim(st_local("pre"))))
            mata: st_local("post", stritrim(strtrim(ustrregexra(st_local("post"), "seldby(var|idx|vars)[(][^)]*[)]", ""))))
            local nmv : word count `moved'
            And_list moved `"`moved'"'
            di as err `"{p 0 4 2}`moved' must follow the comma:{p_end}"'
            di as err `"{p 8 8 2}seldecomp `macval(pre)', `macval(post)'{p_end}"'
            exit 198
        }
        local first by
        local pre `"`macval(cbpre)'"'
        local post `"`macval(cbpost)'"'
        // the option words typed before the if or in qualifier go after the comma too (the
        // text before by( parsed again: an in range followed by by() fails the parse above)
        local 0 `"`macval(cbpre)'"'
        capture syntax [anything(equalok)] [if] [in] [aw pw iw fw]
        if !_rc {
            local keep
            local moved
            local rest `"`macval(anything)'"'
            local k 0
            while `"`macval(rest)'"' != "" {
                gettoken tok rest : rest, bind
                local ++k
                local mv 0
                mata: st_local("pw", strofreal(ustrregexm(st_local("tok"), "^[A-Za-z_][A-Za-z0-9_]*$")))
                if `pw' & `k' > 1 {
                    capture confirm variable `tok'
                    if _rc mata: st_local("mv", strofreal(seldecomp_isopt(st_local("tok"), st_local("bare"))))
                }
                if `mv' {
                    local moved `"`macval(moved)' `macval(tok)'"'
                }
                else {
                    local keep `"`macval(keep)' `macval(tok)'"'
                }
            }
            if `"`macval(moved)'"' != "" {
                local wt
                if "`weight'" != "" local wt `"[`weight'`macval(exp)']"'
                local pre `"`macval(keep)' `macval(if)' `macval(in)' `macval(wt)'"'
                local post `"`macval(cbpost)' `macval(moved)'"'
                local nmv : word count `moved'
                And_list moved `"`moved'"'
                local also = " (and `moved')"
            }
        }
    }
    // (the options that the by prefix appends are left out)
    mata: st_local("pre", stritrim(strtrim(st_local("pre"))))
    mata: st_local("post", stritrim(strtrim(ustrregexra(st_local("post"), "seldby(var|idx|vars)[(][^)]*[)]", ""))))
    di as err `"{p 0 4 2}option `first'() must follow the comma`also':{p_end}"'
    di as err `"{p 8 8 2}seldecomp `macval(pre)', `macval(post)'{p_end}"'
    exit 198
end


* ============================================================================
* Top_error <spec> <command line>: -syntax- rejected the command line, whose options are
*   <spec>. When the first option it does not know is a suboption of bayesboot(), contvars()
*   or generate() typed outside their parentheses, or an option of oaxaca, the
*   message names seldecomp's equivalent (error 198); anything else gets -syntax-'s own
*   message and code
* ============================================================================
program define Top_error
    args spec zero
    local loose = ustrregexra(`"`spec'"', "\([^)]*\)", "(string asis)")
    local 0 `"`macval(zero)'"'
    capture syntax [anything(equalok)] [if] [in] [aw pw iw fw] [, BY(string asis) `loose' * ]
    local rc0 = _rc
    // recut typed outside contvars(): a suboption of contvars(), shown inside it
    if `rc0' == 0 & `"`macval(options)'"' != "" {
        gettoken w1 : options, bind
        if strlower(`"`macval(w1)'"') == "recut" {
            if `"`contvars'"' == "" {
                di as err "recut is a suboption of contvars(): type contvars({it:varname} {it:#}, recut)"
            }
            else if strpos(`"`contvars'"', ",") {
                di as err `"{p 0 4 2}recut is a suboption of contvars(): type contvars(`contvars' recut){p_end}"'
            }
            else di as err `"{p 0 4 2}recut is a suboption of contvars(): type contvars(`contvars', recut){p_end}"'
            exit 198
        }
        // replace typed outside generate(): a suboption of generate(), shown inside it
        if strlower(`"`macval(w1)'"') == "replace" {
            if `"`generate'"' == "" {
                di as err "replace is a suboption of generate(): type generate({it:newvar}, replace)"
            }
            else if strpos(`"`generate'"', ",") {
                di as err `"{p 0 4 2}replace is a suboption of generate(): type generate(`generate' replace){p_end}"'
            }
            else di as err `"{p 0 4 2}replace is a suboption of generate(): type generate(`generate', replace){p_end}"'
            exit 198
        }
    }
    local class
    if `rc0' == 0 & `"`macval(options)'"' != "" {
        // every unknown option of the class of the first one: bayesboot() suboptions are
        // collected (names and the words as typed), an oaxaca option is taken alone
        local names
        local toks
        local seed 0
        local k 0
        local rest `"`macval(options)'"'
        while `"`macval(rest)'"' != "" {
            gettoken tok rest : rest, bind
            local ++k
            mata: st_local("nm", strlower(ustrregexra(st_local("tok"), "[(].*$", "")))
            Habit_class `nm'
            if `k' == 1 {
                if "`hclass'" == "" continue, break
                local class `hclass'
                local oword `hfull'
                mata: st_local("oarg", ustrregexm(st_local("tok"), "^[^(]*[(](.*)[)]$") ? strtrim(ustrregexs(1)) : "")
                // the suboptions already typed inside the option (gname: their full names; gtok_<name>:
                // the word as typed), so that one typed outside as well is shown once
                local gname
                if "`class'" != "oax" {
                    local grest `"`macval(bayesboot2)'"'
                    while `"`macval(grest)'"' != "" {
                        gettoken gtok grest : grest, bind
                        mata: st_local("gnm", strlower(ustrregexra(st_local("gtok"), "[(].*$", "")))
                        Habit_class `gnm'
                        if "`hclass'" == "`class'" {
                            if "`hfull'" == "seed" local hfull rseed
                            if !`: list hfull in gname' {
                                local gname `gname' `hfull'
                                local gtok_`hfull' `"`macval(gtok)'"'
                            }
                        }
                    }
                    // the classification of the word typed outside again
                    Habit_class `nm'
                }
                local dtoks
                local dins
                local nd 0
            }
            if "`hclass'" != "`class'" | "`class'" == "oax" continue
            local ttok `"`macval(tok)'"'
            local par = cond(strpos(`"`macval(tok)'"', "("), "()", "")
            if "`hfull'" == "seed" {
                // seed() is rseed() in bayesboot()
                local seed 1
                local hfull rseed
                mata: st_local("tok", "r" + st_local("tok"))
            }
            // a suboption that takes an argument typed without one: shown with it (dots(1), the
            // default of dots(); a placeholder for the others)
            if "`par'" == "" & "`harg'" != "" {
                local par ()
                local tok = cond("`hfull'" == "dots", "dots(1)", "`hfull'({it:`harg'})")
            }
            local hn `hfull'`par'
            if !`: list hn in names' local names `names' `hn'
            if `: list hfull in gname' {
                // given already (inside the option, or earlier outside it): left out of the
                // suggestion, and named as typed
                local ++nd
                local dtoks `"`macval(dtoks)' `macval(ttok)'"'
                local gt `"`macval(gtok_`hfull')'"'
                if !`: list gt in dins' local dins `"`macval(dins)' `macval(gt)'"'
                continue
            }
            local gname `gname' `hfull'
            local gtok_`hfull' `"`macval(tok)'"'
            local toks `"`macval(toks)' `macval(tok)'"'
        }
    }
    if "`class'" == "bb" {
        local opt bayesboot
        local n : word count `names'
        local lst
        forvalues i = 1/`n' {
            local x : word `i' of `names'
            if `i' == 1 local lst `x'
            else if `i' == `n' local lst `lst' and `x'
            else local lst `lst', `x'
        }
        // the suggestion: what was typed inside the option, then the words typed outside
        local sug `"`macval(bayesboot2)' `macval(toks)'"'
        mata: st_local("sug", stritrim(strtrim(st_local("sug"))))
        local sh = cond(`seed', "; use rseed(), not seed()", "")
        // level(): a replication prefix takes its own
        if `: list posof "level()" in names' {
            local sh "`sh'; with a replication prefix, give level() to the prefix"
        }
        local dh
        if `nd' > 0 {
            And_list dtoks `"`macval(dtoks)'"'
            And_list dins `"`macval(dins)'"'
            local dh `"; `macval(dtoks)' `=cond(`nd' == 1, "is", "are")' left out, as `opt'() has `macval(dins)' already"'
        }
        di as err `"{p 0 4 2}`lst' "' cond(`n' == 1, "is a suboption", "are suboptions") ///
            `" of `opt'(): type `opt'(`macval(sug)')`macval(dh)'`sh'{p_end}"'
        exit 198
    }
    if "`class'" == "oax" {
        // the group oaxaca's swap and weight(0) make the Focal Group here: the higher value of a numeric
        // by(), or with focalgroup() typed the other of its two values, named in focalgroup().
        // The typed text is read as Group_pick reads it: a number equal to one of the values
        // or within a relative 1e-7 of it, or the value label of one; a string by()'s two
        // values are read directly (levelsof would expand a $ or a backquote in them). Text
        // that names neither value leaves the placeholder
        local other "<the other group>"
        if `"`by'"' != "" {
            capture unab byv : `by', max(1)
            if !_rc {
                capture confirm numeric variable `byv'
                if !_rc {
                    capture quietly summarize `byv' `if' `in', meanonly
                    if !_rc & r(N) > 0 & r(min) < r(max) {
                        tempname vmin vmax
                        scalar `vmin' = r(min)
                        scalar `vmax' = r(max)
                        local kf 0
                        // user text: brace blocks (a one-line if would expand it again)
                        if `"`macval(focalgroup)'"' != "" {
                            mata: _sd_x = strtoreal(st_local("focalgroup"))
                            mata: _sd_v = (st_numscalar("`vmin'"), st_numscalar("`vmax'"))
                            mata: _sd_d = (reldif(_sd_x, _sd_v[1]) < 1e-7, reldif(_sd_x, _sd_v[2]) < 1e-7)
                            mata: st_local("kf", strofreal(missing(_sd_x) ? 0 : (_sd_x == _sd_v[1] ? 1 :  ///
                                (_sd_x == _sd_v[2] ? 2 : (_sd_d[1] & !_sd_d[2] ? 1 : (_sd_d[2] & !_sd_d[1] ? 2 : 0))))))
                            local vl : value label `byv'
                            local haslab 0
                            if "`vl'" != "" {
                                mata: st_local("haslab", strofreal(st_vlexists("`vl'")))
                            }
                            if `kf' == 0 & `haslab' {
                                mata: _sd_m = (_sd_v :== trunc(_sd_v)) :* (st_vlmap("`vl'", trunc(_sd_v)) :== st_local("focalgroup"))
                                mata: st_local("kf", strofreal(_sd_m[1] & !_sd_m[2] ? 1 : (_sd_m[2] & !_sd_m[1] ? 2 : 0)))
                            }
                            capture mata: mata drop _sd_x _sd_v _sd_d
                            capture mata: mata drop _sd_m
                        }
                        if `kf' == 1 local other = strtrim(string(`vmax', "%12.0g"))
                        else if `kf' == 2 local other = strtrim(string(`vmin', "%12.0g"))
                        else if `"`macval(focalgroup)'"' == "" {
                            local other = strtrim(string(`vmax', "%12.0g"))
                        }
                    }
                }
                else {
                    capture confirm string variable `byv'
                    if !_rc & `"`macval(focalgroup)'"' != "" {
                        tempvar gcode
                        capture quietly egen long `gcode' = group(`byv') `if' `in'
                        if !_rc {
                            capture quietly summarize `gcode', meanonly
                            if !_rc & r(N) > 0 & r(max) == 2 {
                                mata: st_local("sv1", st_sdata(selectindex(st_data(., "`gcode'") :== 1)[1], "`byv'"))
                                mata: st_local("sv2", st_sdata(selectindex(st_data(., "`gcode'") :== 2)[1], "`byv'"))
                                if `"`macval(focalgroup)'"' == `"`macval(sv1)'"' {
                                    local other `"`macval(sv2)'"'
                                }
                                else if `"`macval(focalgroup)'"' == `"`macval(sv2)'"' {
                                    local other `"`macval(sv1)'"'
                                }
                            }
                        }
                    }
                }
            }
        }
        local sw `"focalgroup(`macval(other)')"'
        if `"`macval(focalgroup)'"' != "" {
            local sw `"focalgroup(`macval(other)') in place of focalgroup(`macval(focalgroup)')"'
        }
        if "`oword'" == "swap" {
            local msg `"swap not allowed; to rake the other group, use `macval(sw)'"'
        }
        else if "`oword'" == "weight" & inlist(`"`oarg'"', "0", "0.0") {
            local msg `"weight(0) not allowed; to rake the other group, use `macval(sw)'"'
        }
        else if "`oword'" == "weight" & inlist(`"`oarg'"', "1", "1.0") {
            local msg "weight(1) not allowed; it is the default"
        }
        else if "`oword'" == "weight" {
            local msg `"weight() not allowed; the default corresponds to weight(1), and `macval(sw)' to weight(0)"'
        }
        else if "`oword'" == "pooled" {
            local msg "pooled not allowed; pooled is a suboption of contvars()"
        }
        else {
            // vce(): the vcetype decides
            mata: st_local("vt", strlower(ustrregexra(st_local("oarg"), "[ ,].*$", "")))
            mata: st_local("cv", ustrregexm(st_local("oarg"), "^[^ ,]+[ ]+([^ ,]+)") ? ustrregexs(1) : "")
            if "`vt'" != "" & strpos("bootstrap", "`vt'") == 1 {
                local msg "vce(bootstrap) not allowed; use bayesboot() or the bootstrap: prefix"
            }
            else if "`vt'" != "" & strpos("jackknife", "`vt'") == 1 {
                local msg "vce(jackknife) not allowed; use the jackknife: prefix"
            }
            else if "`vt'" != "" & strpos("cluster", "`vt'") == 1 & "`cv'" != "" {
                local msg "vce(cluster) not allowed; use bayesboot(cluster(`cv'))"
            }
            else {
                local msg "vce() not allowed; use bayesboot() or a replication prefix"
            }
        }
        di as err `"{p 0 4 2}`macval(msg)'{p_end}"'
        exit 198
    }
    // anything else: -syntax-'s own message and return code
    local 0 `"`macval(zero)'"'
    syntax anything(id="varlist" equalok) [if] [in] [aw pw iw fw] , BY(varname) [ `spec' ]
    exit 198
end


* ============================================================================
* And_list <lmacname> "<words>": c_local <lmacname>, the words (options, bound as gettoken
*   binds them) joined as "a", "a and b", "a, b and c"
* ============================================================================
program define And_list
    args lm txt
    local n 0
    local rest `"`macval(txt)'"'
    while `"`macval(rest)'"' != "" {
        gettoken tok rest : rest, bind
        local ++n
        local t`n' `"`macval(tok)'"'
    }
    local out
    forvalues i = 1/`n' {
        if `i' == 1 local out `"`macval(t`i')'"'
        else if `i' == `n' local out `"`macval(out)' and `macval(t`i')'"'
        else local out `"`macval(out)', `macval(t`i')'"'
    }
    c_local `lm' `"`macval(out)'"'
end


* ============================================================================
* Gen_role <varname> "<inputs>" "<option>": refuse (198) to overwrite an input of the command:
*   <inputs> are pairs of a role and its variables ("the by() variable" "group" ...), and the
*   message names <option>, the option that would overwrite <varname>
* ============================================================================
program define Gen_role
    args gn inputs opt
    while `"`inputs'"' != "" {
        gettoken role inputs : inputs
        gettoken vars inputs : inputs
        if `: list gn in vars' {
            di as err "{p 0 4 2}`opt' would overwrite `gn', `role'; choose another name{p_end}"
            exit 198
        }
    }
end


* ============================================================================
* Habit_class <word>: c_local hclass (bb: a suboption of bayesboot(),
*   oax: an option of oaxaca; empty otherwise), hfull, the full name, and harg, the kind of
*   its argument (#, varname, filename, keyword; empty when it takes none); a word matches
*   when it is the name abbreviated to no fewer letters than the minimum abbreviation (in the
*   option list of bayesboot(); oaxaca's weight() and pooled by analogy)
* ============================================================================
program define Habit_class
    args w
    c_local hclass
    c_local hfull
    c_local harg
    foreach e in "bb reps r #" "bb rseed rseed #" "bb seed seed #" "bb cluster cl varname"      ///
        "bb saving sav filename" "bb dots dots #" "bb nodots nodots"                             ///
        "bb bca bca" "bb ties ties" "bb level l #" "bb citable cit" "bb normal norm"             ///
        "bb tinterval tint" "bb percentile perc" "bb bc bc"                                     ///
        "oax swap swap" "oax weight w" "oax pooled pool" "oax vce vce" {
        gettoken c e : e
        gettoken f e : e
        gettoken m a : e
        local a = strtrim("`a'")
        if strlen("`w'") >= strlen("`m'") & substr("`f'", 1, strlen("`w'")) == "`w'" {
            c_local hclass `c'
            c_local hfull `f'
            c_local harg `a'
            exit
        }
    }
end


* ============================================================================
* Replay_error <spec> <options>: -syntax- rejected the options of a replay, whose options
*   are <spec>. An option replay does not take (an estimation option, a misspelling) gets a
*   message listing those it takes (error 198); a bad argument gets -syntax-'s own message
* ============================================================================
program define Reject_nodiagnostics
    // -syntax- silently consumes the no-form of a switch. Inspect the original
    // top-level option tokens, binding parentheses/quotes, before that happens.
    capture syntax [anything(equalok)] [if] [in] [aw pw iw fw] [, * ]
    if _rc exit
    local rest `"`macval(options)'"'
    while `"`macval(rest)'"' != "" {
        gettoken tok rest : rest, bind
        mata: st_local("bad", strofreal(strlen(st_local("tok")) >= 6 & substr("nodiagnostics", 1, strlen(st_local("tok"))) == st_local("tok")))
        if `bad' {
            di as err "option nodiagnostics not allowed"
            exit 198
        }
    }
end

program define Replay_error
    args spec zero
    local loose = ustrregexra(`"`spec'"', "\([^)]*\)", "(string asis)")
    local 0 `"`macval(zero)'"'
    capture syntax [, `loose' * ]
    if _rc == 0 & `"`macval(options)'"' != "" {
        gettoken bad : options, bind
        mata: st_local("nm", ustrregexra(st_local("bad"), "[(].*$", "()"))
        di as err `"{p 0 4 2}option `macval(nm)' not allowed at replay; replay takes level(), diagnostics, covariatedetail, citable, varlabels, cformat(), pformat(), sformat(), and makeplot(){p_end}"'
        exit 198
    }
    local 0 `"`macval(zero)'"'
    syntax [, `spec' ]
    exit 198
end


* ============================================================================
* Worst_share <raked w> <term codes> <level> <touse> <g1> <w> <W2>: c_local wshare, the
*   raked Focal Group share of the level against its target, Reference Group's share (the level shares
*   of the Reference Group sum to 1 over the levels raked: a level without Focal Group observations is a
*   support failure), at two decimals, or six when two do not tell them apart
* ============================================================================
program define Worst_share
    args rw tv l touse g1 w W2
    qui sum `rw' if `touse' & `g1' == 1, meanonly
    local rtot = r(sum)
    qui sum `rw' if `touse' & `g1' == 1 & `tv' == `l', meanonly
    local rsh = r(sum) / `rtot'
    qui sum `w' if `touse' & `g1' == 0 & `tv' == `l', meanonly
    local tsh = r(sum) / `W2'
    local pf %9.2f
    if string(100 * `rsh', "`pf'") == string(100 * `tsh', "`pf'") local pf %12.6f
    local rsh = strtrim(string(100 * `rsh', "`pf'"))
    local tsh = strtrim(string(100 * `tsh', "`pf'"))
    c_local wshare "raked share `rsh'%, target `tsh'%"
end


* ============================================================================
* Miss_error "<prefix>" <nstages> <stage> <nbad> <itstop> "<wtxt>" [<sup>]: the message of
*   a raking that did not reach its targets (error 430), for the forward order and the
*   reverse one (<prefix> "reverse order: "); <itstop> 1 when the raking stopped at its
*   iterate() limit before the weights converged, the one case where raising it helps
* ============================================================================
program define Miss_error
    args pre nst s nbad itstop wtxt sup trim
    // the worst level's label is user text: SMCL braces escaped
    local wtxt : subinstr local wtxt "{" "`=char(1)'", all
    local wtxt : subinstr local wtxt "}" "{c )-}", all
    local wtxt : subinstr local wtxt "`=char(1)'" "{c -(}", all
    local itxt = cond(`itstop', " (raking stopped at iterate())", "")
    // sup: a level absent from the Reference Group was missed: the lower trimming bounds kept its
    // Focal Group observations above zero weight
    if 0`sup' {
        di as err `"{p 0 4 2}`pre'raking did not reach its targets"' cond(`nst' > 1, " at stage `s'", "") ///
            ": `nbad' " plural(`nbad', "margin") " missed by more than ctrltolerance()" ///
            `"`macval(wtxt)'; the lower trimming bounds keep focal group observations with levels "' ///
            "absent from the reference group above 0; specify trimsupportoverride, insupport, or force{p_end}"
        exit 430
    }
    // trim: trimming changed the stage's raked weights: the bounds keep the margins from their
    // targets, whether or not the raking then stopped at iterate()
    if 0`trim' {
        di as err `"{p 0 4 2}`pre'raking did not reach its targets"' cond(`nst' > 1, " at stage `s'", "") ///
            ": `nbad' " plural(`nbad', "margin") " missed by more than ctrltolerance()" ///
            `"`macval(wtxt)'; the trimming bounds are binding: relax them, or specify force to keep the results{p_end}"'
        exit 430
    }
    di as err `"{p 0 4 2}`pre'raking did not reach its targets"' cond(`nst' > 1, " at stage `s'", "") ///
        ": `nbad' " plural(`nbad', "margin") " missed by more than ctrltolerance()`itxt'" ///
        `"`macval(wtxt)'; specify force to keep the results{p_end}"'
    exit 430
end


* ============================================================================
* Zero_terms "<list>" <total> "<names>": <list> holds t:l:count for each level l of term t
*   absent from the Reference Group, with its number of the Focal Group observations; <names> the term names.
*   c_local zsole: the term whose levels cover all <total> Focal Group observations (it has no
*   level found in both groups), else empty; zctxt: each term with its count, as
*   "x6: 500 obs; x7: 500 obs"
* ============================================================================
program define Zero_terms
    args zlist tot names
    local ts
    foreach zz of local zlist {
        gettoken t zr : zz, parse(":")
        gettoken colon zr : zr, parse(":")
        gettoken l zr : zr, parse(":")
        gettoken colon c : zr, parse(":")
        if !`: list t in ts' {
            local ts `ts' `t'
            local n`t' 0
        }
        local n`t' = `n`t'' + `c'
    }
    local sole
    local txt
    foreach t of local ts {
        if `n`t'' == `tot' & "`sole'" == "" local sole `t'
        local sepr = cond("`txt'" == "", "", "; ")
        local txt `"`txt'`sepr'`: word `t' of `names'': `=strtrim(string(`n`t'', "%12.0fc"))' obs"'
    }
    c_local zsole `sole'
    c_local zctxt `"`txt'"'
end


* ============================================================================
* Group_values <lmacname> <byvar> <touse>: the values of the by() variable in the
*   estimation sample, as the group-value error messages list them: 0 (Single) and
*   1 (Married); "a" and "b"; 1, 2, and 3; at most 6, then "...". <lmacname>_ex holds
*   the first two as inlist() arguments
* ============================================================================
program define Group_values
    args lm by touse
    capture confirm string variable `by'
    local isstr = (_rc == 0)
    if `isstr' {
        // read the strings directly (macro text would re-expand $ and backquotes)
        mata: _sd_gv = uniqrows(st_sdata(., "`by'", "`touse'"))
        mata: st_local("n", strofreal(rows(_sd_gv)))
    }
    else {
        qui levelsof `by' if `touse', local(vals)
        local n : word count `vals'
        local vl : value label `by'
    }
    local txt
    local ex
    forvalues j = 1/`=min(`n', 6)' {
        if `isstr' {
            mata: st_local("v", _sd_gv[`j'])
            local v `""`macval(v)'""'
            local piece `"`macval(v)'"'
        }
        else {
            local v : word `j' of `vals'
            local piece `v'
            if "`vl'" != "" & `v' == int(`v') {
                local lab : label `vl' `v', strict
                if `"`macval(lab)'"' != "" {
                    local piece `"`v' (`macval(lab)')"'
                }
            }
        }
        if `j' == 1 {
            local txt `"`macval(piece)'"'
            local ex `"`macval(v)'"'
        }
        else {
            if `j' == 2 {
                local ex `"`macval(ex)', `macval(v)'"'
            }
            local sep ", "
            if `j' == `n' local sep = cond(`n' == 2, " and ", ", and ")
            local txt `"`macval(txt)'`sep'`macval(piece)'"'
        }
    }
    if `n' > 6 {
        local txt `"`macval(txt)', ..."'
    }
    capture mata: mata drop _sd_gv
    // labels and string values are user text, shown in messages: SMCL braces escaped
    foreach x in txt ex {
        local `x' : subinstr local `x' "{" "`=char(1)'", all
        local `x' : subinstr local `x' "}" "{c )-}", all
        local `x' : subinstr local `x' "`=char(1)'" "{c -(}", all
    }
    c_local `lm' `"`macval(txt)'"'
    c_local `lm'_ex `"`macval(ex)'"'
end


* ============================================================================
* Group_pick <lmacname> <opt> <byvar> <touse> <gs1> <gs2> <txt> [<zw>]: which of the two
*   values of a numeric by(), 1 (scalar <gs1>) or 2 (<gs2>), the text of <opt>() names (held
*   in the string scalar <txt>): a value first; if none, a value label (exact,
*   case-sensitive). Text that is a value and the label of the other value, or the label of
*   both, is ambiguous (198); text that is neither is 459 (198 for a text that is not a
*   number when by() has no value label); each message lists the values and their labels,
*   and the 459 says so when the text names observations left out for a zero weight alone
*   (marked in <zw>)
* ============================================================================
program define Group_pick
    args lm opt by touse gs1 gs2 txt zw
    mata: st_local("t", st_strscalar("`txt'"))
    // a value: the text read as a number equals one of the two values; an exact match first,
    // else the one value within a relative 1e-7 (a float by() typed as a decimal). When both
    // values are within it (large codes such as yyyymmdd dates), only an exact match counts
    mata: _sd_x = strtoreal(st_strscalar("`txt'"))
    mata: st_local("isnum", strofreal(!missing(_sd_x)))
    mata: _sd_v = (st_numscalar("`gs1'"), st_numscalar("`gs2'"))
    mata: _sd_d = (reldif(_sd_x, _sd_v[1]) < 1e-7, reldif(_sd_x, _sd_v[2]) < 1e-7)
    mata: st_local("kv", strofreal(missing(_sd_x) ? 0 : (_sd_x == _sd_v[1] ? 1 :            ///
        (_sd_x == _sd_v[2] ? 2 : (_sd_d[1] & !_sd_d[2] ? 1 : (_sd_d[2] & !_sd_d[1] ? 2 : 0))))))
    capture mata: mata drop _sd_x _sd_v _sd_d
    // a value label: the labels of the two values (integers only) that equal the text
    local vl : value label `by'
    local haslab 0
    local nl 0
    local kl 0
    if "`vl'" != "" {
        mata: st_local("haslab", strofreal(st_vlexists("`vl'")))
    }
    if `haslab' {
        mata: _sd_v = (st_numscalar("`gs1'"), st_numscalar("`gs2'"))
        mata: _sd_m = (_sd_v :== trunc(_sd_v)) :* (st_vlmap("`vl'", trunc(_sd_v)) :== st_strscalar("`txt'"))
        mata: st_local("nl", strofreal(sum(_sd_m)))
        mata: st_local("kl", strofreal(_sd_m[1] ? 1 : (_sd_m[2] ? 2 : 0)))
        capture mata: mata drop _sd_v _sd_m
    }
    if `nl' == 2 | (`kv' & `nl' == 1 & `kl' != `kv') {
        Group_values gvl `by' `touse'
        local why = cond(`nl' == 2, "the value label of both values", ///
            "a value of `by' and the value label of the other value")
        di as err `"{p 0 4 2}`opt'(`macval(t)') is ambiguous: `why'; its values are `macval(gvl)'{p_end}"'
        exit 198
    }
    if `kv' {
        c_local `lm' `kv'
        exit
    }
    if `kl' {
        if `isnum' {
            local v = strofreal(`gs`kl'', "%21.0g")
            di as txt `"{p 0 6 2 78}note: `opt'(`macval(t)'): `macval(t)' is not a value of `by'; "' ///
                `"it is taken as the value label of `v'{p_end}"'
        }
        c_local `lm' `kl'
        exit
    }
    Group_values gvl `by' `touse'
    // the text may name a value whose observations all have zero weight (Zero_named)
    Zero_named zwhy `by' `txt' `zw'
    if `haslab' {
        di as err `"{p 0 4 2}`opt'(`macval(t)') is neither a value nor a value label of `by'`macval(zwhy)'; its values are `macval(gvl)'{p_end}"'
        exit 459
    }
    if !`isnum' {
        di as err `"{p 0 4 2}`opt'(`macval(t)') must be a number; the values of `by' are `macval(gvl)'{p_end}"'
        exit 198
    }
    di as err `"{p 0 4 2}`opt'(`macval(t)') is not a value of `by'`macval(zwhy)'; its values are `macval(gvl)'{p_end}"'
    exit 459
end


* ============================================================================
* Zero_group <lmacname> <byvar> <touse> <zw> <opt> <txt>: by() takes one value on <touse>;
*   when the observations of the other value(s) all have zero weight (marked in <zw>), the
*   text of that reason, naming the group (Focal Group or 2) when it can be told, else empty.
*   <opt> is focalgroup or none; the string scalar <txt> holds its text
* ============================================================================
program define Zero_group
    args lm by touse zw opt txt
    c_local `lm'
    capture confirm string variable `by'
    local isstr = (_rc == 0)
    tempvar zo
    tempname gk gz
    if `isstr' {
        // the values read and compared as data (macro text would re-expand $ and backquotes)
        qui gen byte `zo' = 0
        mata: _sd_k = st_sdata(selectindex(st_data(., "`touse'"))[1], "`by'")
        mata: st_store(., "`zo'", st_data(., "`zw'") :& (st_sdata(., "`by'") :!= _sd_k))
    }
    else {
        qui sum `by' if `touse', meanonly
        scalar `gk' = r(min)
        qui gen byte `zo' = `zw' & `by' != `gk'
    }
    qui count if `zo'
    local nz = r(N)
    if `nz' == 0 {
        capture mata: mata drop _sd_k
        exit
    }
    Group_values gvz `by' `zo'
    if `isstr' {
        mata: _sd_z = uniqrows(st_sdata(., "`by'", "`zo'"))
        mata: st_local("nv", strofreal(rows(_sd_z)))
    }
    else {
        mata: st_local("nv", strofreal(rows(uniqrows(st_data(., "`by'", "`zo'")))))
    }
    if `nv' > 1 {
        c_local `lm' `"the observations with the other values of `by' (`macval(gvz)') and nonmissing data all have zero weight"'
        capture mata: mata drop _sd_k _sd_z
        exit
    }
    // the role of the group left out: Focal Group is the lower value unless focalgroup() names
    // one (a text that names neither value leaves the role unsaid)
    local role 0
    if `isstr' {
        if "`opt'" == "none" {
            mata: st_local("zlow", strofreal(_sd_z[1] < _sd_k))
            local role = cond(`zlow', 1, 2)
        }
        else {
            mata: st_local("zhit", strofreal(st_strscalar("`txt'") == _sd_z[1] ? 2 : ///
                (st_strscalar("`txt'") == _sd_k ? 1 : 0)))
            if `zhit' local role = 3 - `zhit'
        }
    }
    else {
        qui sum `by' if `zo', meanonly
        scalar `gz' = r(min)
        if "`opt'" == "none" local role = cond(`gz' < `gk', 1, 2)
        else {
            // as Group_pick reads the option: a value, else a value label
            capture Group_pick kz `opt' `by' `zo' `gk' `gz' `txt'
            if _rc == 0 local role = 3 - `kz'
        }
    }
    capture mata: mata drop _sd_k _sd_z
    // "with nonmissing data": the group's observations out for a missing value are not counted
    local lead = cond(`nz' == 1, "the one", "all `nz'")
    local obs = plural(`nz', "observation")
    local have = cond(`nz' == 1, "has", "have")
    if `role' {
        local rname = cond(`role' == 1, "focal group", "reference group")
        c_local `lm' `"`lead' `rname' `obs' (`by' = `macval(gvz)') with nonmissing data `have' zero weight"'
    }
    else {
        c_local `lm' `"`lead' `obs' of the other group (`by' = `macval(gvz)') with nonmissing data `have' zero weight"'
    }
end


* ============================================================================
* Zero_named <lmacname> <byvar> <txt> [<zw>]: the text of focalgroup() (the
*   string scalar <txt>) names no value of <byvar> in the estimation sample; when it names
*   the value of observations left out for a zero weight alone (marked in <zw>; a value,
*   or for numeric by() also a value label, read as Group_pick reads them), ": " and the
*   text of that reason, else empty
* ============================================================================
program define Zero_named
    args lm by txt zw
    c_local `lm'
    if "`zw'" == "" exit
    qui count if `zw'
    if r(N) == 0 exit
    capture confirm string variable `by'
    local isstr = (_rc == 0)
    tempvar zn
    qui gen byte `zn' = 0
    if `isstr' {
        // compared as data (macro text would re-expand $ and backquotes)
        mata: st_store(., "`zn'", st_data(., "`zw'") :& (st_sdata(., "`by'") :== st_strscalar("`txt'")))
    }
    else {
        tempname x
        mata: st_numscalar("`x'", strtoreal(st_strscalar("`txt'")))
        if !missing(`x') {
            // as Group_pick: the value equal to the text, else the one value within 1e-7
            mata: _sd_v = uniqrows(st_data(., "`by'", "`zw'"))
            mata: _sd_e = select(_sd_v, _sd_v :== st_numscalar("`x'"))
            mata: _sd_e = (rows(_sd_e) ? _sd_e : select(_sd_v, reldif(_sd_v, J(rows(_sd_v), 1, st_numscalar("`x'"))) :< 1e-7))
            mata: st_local("nxv", strofreal(rows(_sd_e)))
            if `nxv' == 1 {
                mata: st_numscalar("`x'", _sd_e[1])
                qui replace `zn' = 1 if `zw' & `by' == `x'
            }
            capture mata: mata drop _sd_v _sd_e
        }
        local vl : value label `by'
        local haslab 0
        if "`vl'" != "" {
            mata: st_local("haslab", strofreal(st_vlexists("`vl'")))
        }
        if `haslab' {
            // the integer values of the left-out observations whose value label is the text
            mata: _sd_v = uniqrows(st_data(., "`by'", "`zw'"))
            mata: _sd_v = select(_sd_v, _sd_v :== trunc(_sd_v))
            mata: _sd_v = select(_sd_v, st_vlmap("`vl'", _sd_v) :== st_strscalar("`txt'"))
            mata: st_local("nlv", strofreal(rows(_sd_v)))
            forvalues j = 1/`nlv' {
                mata: st_numscalar("`x'", _sd_v[`j'])
                qui replace `zn' = 1 if `zw' & `by' == `x'
            }
            capture mata: mata drop _sd_v
        }
    }
    qui count if `zn'
    local nz = r(N)
    if `nz' == 0 exit
    Group_values gvz `by' `zn'
    if `nz' == 1 {
        c_local `lm' `": the one observation with `by' = `macval(gvz)' and nonmissing data has zero weight"'
    }
    else {
        c_local `lm' `": all `nz' observations with `by' = `macval(gvz)' and nonmissing data have zero weight"'
    }
end


* ============================================================================
* Bin_quantiles: quantile-group variable for a continuous raking variable
*   cut points from _pctile on the requested sample (pooled/focalgroup/refgroup), or those
*   of the 1 x (q - 1) matrix fixed() (a replicate of a replication prefix)
*   returns (c_local) s_ncuts, s_cuts (formatted), s_nbins
* ============================================================================
program define Bin_quantiles
    syntax varname(numeric), touse(varname) g1(varname) w(varname) q(integer) sample(string) ///
        gen(name) [ cutmat(name) fixed(name) ]
    local x `varlist'
    if "`sample'" == "pooled"          local cond `touse'
    else if "`sample'" == "focalgroup" local cond `touse' & `g1' == 1
    else                               local cond `touse' & `g1' == 0
    if "`fixed'" == "" {
        qui count if `cond'
        if r(N) < `q' {
            di as err "{p 0 4 2}contvars(): too few observations to form `q' quantile groups of `x'{p_end}"
            exit 2001
        }
        qui _pctile `x' [aw=`w'] if `cond', nq(`q')
    }
    local cuts
    local cnames
    forvalues j = 1/`=`q'-1' {
        tempname c`j'
        if "`fixed'" == "" scalar `c`j'' = r(r`j')
        else scalar `c`j'' = `fixed'[1, `j']
        local cuts `cuts' `: di %9.0g scalar(`c`j'')'
        local cnames `cnames' cut`j'
    }
    if "`cutmat'" != "" {
        matrix `cutmat' = J(1, `q' - 1, .)
        forvalues j = 1/`=`q'-1' {
            matrix `cutmat'[1, `j'] = scalar(`c`j'')
        }
        matrix colnames `cutmat' = `cnames'
        matrix rownames `cutmat' = `x'
    }
    qui gen int `gen' = `q' if `touse'
    forvalues j = `=`q'-1'(-1)1 {
        qui replace `gen' = `j' if `touse' & `x' <= scalar(`c`j'')
    }
    qui tab `gen' if `touse'
    c_local s_nbins `r(r)'
    c_local s_ncuts `=`q'-1'
    c_local s_cuts `cuts'
end


* ============================================================================
* Bin_label: text label for quantile group -bin- given the cut points
* ============================================================================
program define Bin_label
    syntax, bin(integer) ncuts(integer) [ cuts(string) ]
    local q = `ncuts' + 1
    if `bin' == 1 {
        local c1 : word 1 of `cuts'
        local lab "<= `c1'"
    }
    else if `bin' == `q' {
        local ck : word `ncuts' of `cuts'
        local lab "> `ck'"
    }
    else {
        local ca : word `=`bin'-1' of `cuts'
        local cb : word `bin' of `cuts'
        local lab "(`ca', `cb']"
    }
    c_local s_label `"`lab'"'
end


* ============================================================================
* Engine: raking decomposition for one weight vector, all stages
*   res(): nstages x 30 matrix, columns
*     1 raked_focal 2 explained 3 unexplained 4 pct_explained 5 increment
*     6 converged 7 badcontrols 8 maxreldif 9 maxctrl
*    10 raked_mean 11 raked_min 12 raked_max 13 raked_sd 14 raked_cv
*    15 factor_mean 16 factor_min 17 factor_max 18 factor_sd 19 factor_cv
*    20 N_raked 21 N_focal_zero 22 N_at_bound 23 N_eff_base 24 N_eff_raked 25 ess_ratio
*    26 N_focal_mincell 27 pct_increment 28 pct_W_focal_zero 29 mean_focal_zero 30 mean_focal_cs
*   returns r(rc) (0 ok, 1 Break inside the raking, 459 support failure, 2000 zero weight or
*   nothing left to rake, other = raking error),
*   r(m1) r(m2) r(W1) r(W2), r(nfail) r(fail#_t) r(fail#_l), r(msg1),
*   r(zero#_*) r(thin#_*) (cells for the notes: zero raked weight, thin), r(held#_*) (the cells
*   held at a lower trimming bound, for the warning), r(supmiss#) (their margins missed),
*   r(trimchg#) (1 when stage # missed a margin with a trimming option
*   and trimming changed its raked weights: the bounds are then the cause of the miss)
*   checks(): optional stage x 5 sidecar (raw iteration convergence, missed terms,
*   missed level targets, checked level targets, control scale); misses(): optional rows
*   stage, term, level, target, achieved, discrepancy, the worst missed level per term.
*   Target and achieved totals are in original weight units; discrepancy uses the scale.
*   ylist(): further outcomes, raked with the same weights: outcome j (2, 3, ...) fills the
*   matrix named by word j-1 of resx(), with the columns of res() (1-5, 27, 29, 30 its own,
*   the others copied), and returns r(m1_j) r(m2_j) (also for j = 1). den(): every mean is
*   the ratio of the means of the outcome and den under the same weights; r(n1_j) r(n2_j)
*   r(n1r_j) (numerators: groups, raked at the last stage) and r(d1) r(d2) r(d1r)
*   contrib (individually): column 5 of stage s < nstages is the stage's own explained
*   part (its block's contribution), of the last stage explained minus their sum (remainder)
* ============================================================================
program define Engine, rclass
    syntax, touse(varname) g1(varname) w(varname) y(varname) one(varname) nterms(integer) ///
        tvlist(string) nlevlist(string) nstages(integer) stagedefs(string) res(name)     ///
        [ rwlist(string) ipfopts(string asis) tnames(string) trimhiabs(real -1) ///
          trimloabs(real -1) mincell(integer 10) holdothers contrib NOIsily              ///
          ylist(string) resx(string) den(string) skipdiag trimsupportoverride ///
          checks(name) misses(name) ]
    local ny = 1 + `: word count `ylist''

    forvalues t = 1/`nterms' {
        local tv`t'   : word `t' of `tvlist'
        local nlev`t' : word `t' of `nlevlist'
        local tn`t'   : word `t' of `tnames'
        if "`tn`t''" == "" local tn`t' term`t'
    }
    local s 0
    local rest `stagedefs'
    while `"`rest'"' != "" {
        gettoken piece rest : rest, parse("|")
        if "`piece'" == "|" continue
        local ++s
        local stg`s' `piece'
    }
    if `s' != `nstages' {
        di as err "seldecomp: internal error (stage definitions)"
        exit 498
    }
    if "`noisily'" != "" local noi noisily
    // trimming options: the absolute bounds are Engine options, the relative
    // bounds and trimfrequency() are in ipfopts. ipfnotrim: ipfopts without them, for the
    // untrimmed raking that tells whether trimming changed the weights (below); tolchg: the
    // change of an adjustment factor that counts, the convergence tolerance
    local hastrim = (`trimhiabs' != -1 | `trimloabs' != -1)
    local ipfnotrim
    local tolchg 1e-6
    local haslorel 0
    foreach o of local ipfopts {
        if substr(`"`o'"', 1, 4) == "trim" local hastrim 1
        else local ipfnotrim `ipfnotrim' `o'
        if regexm(`"`o'"', "^tolerance\((.+)\)$") local tolchg = regexs(1)
        if substr(`"`o'"', 1, 10) == "trimlorel(" local haslorel 1
    }
    // a lower trimming bound binds on the Focal Group observations with a level absent from the
    // Reference Group as on any other: they stay in the raking (incl), their level's target 0,
    // and the bound holds them; the margin it leaves unmet fails the stage (r(supmiss#)).
    // trimsupportoverride, or no lower bound: they are set aside with zero raked weight
    local incl = ("`trimsupportoverride'" == "" & (`trimloabs' != -1 | `haslorel'))

    // group means and totals (zero weights contribute nothing)
    qui sum `y' [iw=`w'] if `touse' & `g1' == 1, meanonly
    local m1 = r(mean)
    local W1 = r(sum_w)
    qui sum `y' [iw=`w'] if `touse' & `g1' == 0, meanonly
    local m2 = r(mean)
    local W2 = r(sum_w)
    if `W1' <= 0 | `W2' <= 0 | missing(`m1') | missing(`m2') {
        return scalar rc = 2000
        return local msg1 "one of the groups has zero total weight"
        exit
    }
    local m1_1 `m1'
    local m2_1 `m2'
    forvalues j = 2/`ny' {
        local yv : word `=`j' - 1' of `ylist'
        qui sum `yv' [iw=`w'] if `touse' & `g1' == 1, meanonly
        local m1_`j' = r(mean)
        qui sum `yv' [iw=`w'] if `touse' & `g1' == 0, meanonly
        local m2_`j' = r(mean)
    }
    // ratio of means: the numerators above over the means of den under the same weights
    if "`den'" != "" {
        qui sum `den' [iw=`w'] if `touse' & `g1' == 1, meanonly
        local d1 = r(mean)
        qui sum `den' [iw=`w'] if `touse' & `g1' == 0, meanonly
        local d2 = r(mean)
        if `d1' == 0 | `d2' == 0 | missing(`d1') | missing(`d2') {
            return scalar rc = 459
            return local msg1 "ratio(): the mean of `den' is zero in the `=cond(`d1' == 0 | missing(`d1'), "Focal", "Reference")' Group"
            exit
        }
        forvalues j = 1/`ny' {
            local n1_`j' `m1_`j''
            local n2_`j' `m2_`j''
            local m1_`j' = `n1_`j'' / `d1'
            local m2_`j' = `n2_`j'' / `d2'
        }
        local m1 `m1_1'
        local m2 `m2_1'
    }
    local diff = `m1' - `m2'
    local diff_1 `diff'
    forvalues j = 2/`ny' {
        local diff_`j' = `m1_`j'' - `m2_`j''
        local prevexpl_`j' 0
        local X`j' : word `=`j' - 1' of `resx'
        matrix `X`j'' = J(`nstages', 30, .)
    }

    tempvar rs zero rss rw rwtmp
    qui gen byte `rs' = `touse' & `g1' == 1 & `w' > 0
    qui gen byte `zero' = 0
    qui gen byte `rss' = 0
    qui gen double `rw' = .
    matrix `res' = J(`nstages', 30, .)
    if "`checks'" != "" matrix `checks' = J(`nstages', 5, .)
    if "`misses'" != "" capture matrix drop `misses'
    local nmissrows 0
    tempname RMISS MROW
    local nzs 0
    local nhs 0
    local nts 0
    local prevexpl 0
    // every term's Reference Group shares, once per call (they do not depend on the stage); one
    // control-matrix name per term, reused by every stage
    tempvar s2
    tempname P2 F
    qui gen byte `s2' = `touse' & `g1' == 0 & `w' > 0
    mata: seldecomp_g2shares("`tvlist'", "`nlevlist'", "`w'", "`s2'", `W2', "`P2'")
    drop `s2'
    local mnames
    forvalues t = 1/`nterms' {
        tempname MT`t'
        local mnames `mnames' `MT`t''
    }

    forvalues s = 1/`nstages' {
        local sterms `stg`s''
        local oterms
        if "`holdothers'" != "" {
            forvalues t = 1/`nterms' {
                if !`: list t in sterms' local oterms `oterms' `t'
            }
        }
        // Reference Group targets of the stage's terms (seldecomp_stage1, from the shares in P2): a
        // level absent from the Reference Group gets zero target -> zero raked weight (the zero flag); the
        // raking sample rss, the thin target levels, the support failures, the control matrices
        mata: seldecomp_stage1("`tvlist'", "`nlevlist'", "`sterms'", "`rs'", "`zero'", "`rss'", ///
            "`P2'", `W1', `mincell', "`one'", "`mnames'", "`F'", `incl')
        if `Nin' == 0 {
            return scalar rc = 2000
            // the culprit: the terms whose levels absent from the Reference Group set Focal Group observations
            // aside, with their counts; r(septerm): a term that sets them all aside alone (it has
            // no level found in both groups: the caller asks whether it is continuous or a copy of by())
            local tnl
            forvalues t = 1/`nterms' {
                local tnl `tnl' `tn`t''
            }
            Zero_terms "`zerolist'" `nsup' "`tnl'"
            if "`zsole'" != "" {
                return scalar septerm = `zsole'
                return local msg1 "stage `s': every focal group observation has a level of `tn`zsole'' absent from the reference group"
            }
            else {
                return local msg1 "stage `s': every focal group observation has a level absent from the reference group (`zctxt')"
            }
            local engstop stop
            continue, break
        }
        if `nfail' > 0 {
            forvalues f = 1/`nfail' {
                return scalar fail`f'_t = `F'[`f', 1]
                return scalar fail`f'_l = `F'[`f', 2]
                return scalar fail`f'_z = `F'[`f', 3]
            }
            return scalar rc = 459
            return scalar nfail = `nfail'
            return scalar stage = `s'
            return scalar nzero = `nzero'
            return local zerolist `zerolist'
            local engstop stop
            continue, break
        }
        // holdothers (fixothermargins): the other raking variables pinned at the Focal Group's own shares,
        // computed on the raking sample (observations with zero raked weight excluded)
        qui sum `w' if `rss', meanonly
        local Wrss = r(sum)
        if "`oterms'" != "" {
            mata: seldecomp_stage2("`tvlist'", "`nlevlist'", "`oterms'", "`w'", "`rss'", `Wrss', ///
                `W1', "`one'", "`mnames'")
        }
        // rake; the base weights are rescaled to mean 1 on the raking sample (the absolute
        // trimming bounds with them) and scaled back afterwards
        qui sum `w' if `rss', meanonly
        local Wrss = r(sum)
        local sc = `Nrake' / `Wrss'
        local ntrim 0
        local tchg 0
        local badtargets 0
        local ntargets 0
        // Constant controls still constrain the total weight. With trimming, send that
        // constraint through the raker so bounds, conflicts and missed totals are checked.
        // Use the first stage term's retained level to identify it in the diagnostics.
        if "`mats'" == "" & `hastrim' {
            tempname Mtotal
            local tt : word 1 of `sterms'
            qui sum `tv`tt'' if `rss', meanonly
            mata: seldecomp_ctrlmat("`Mtotal'", 1, `r(min)', `W1', "`one'", "`tv`tt''")
            local mats `Mtotal'
        }
        if "`mats'" == "" {
            qui replace `rw' = `w' * `W1' / `Wrss' if `rss'
            local conv 1
            local badc 0
            local mrd 0
            local mctrl 0
        }
        else {
            tempvar ws
            qui gen double `ws' = `w' * `sc' if `rss'
            local nmats : word count `mats'
            forvalues m = 1/`nmats' {
                local Mm : word `m' of `mats'
                matrix `Mm' = `Mm' * `sc'
            }
            local trimopts
            if `trimhiabs' != -1 local trimopts `trimopts' trimhiabs(`=`trimhiabs' * `sc'')
            if `trimloabs' != -1 local trimopts `trimopts' trimloabs(`=`trimloabs' * `sc'')
            capture drop `rwtmp'
            local trace
            local missopts
            if "`misses'" != "" local missopts misses(`RMISS')
            if "`noi'" != "" {
                local trace trace
                if `nstages' > 1 di as txt _n "Stage `s'"
            }
            capture `noi' Rake_mata `ws', touse(`rss') ctotal(`mats') generate(`rwtmp') ///
                `ipfopts' `trimopts' `trace' `missopts'
            local rc = _rc
            // a Break (rc 1) stops Engine; the untrimmed re-raking below handles it the same way
            // (continue, break must stay in this loop, so the three lines cannot move to a helper)
            if `rc' == 1 {
                return scalar rc = 1
                local engstop stop
                continue, break
            }
            if `rc' {
                return scalar rc = `rc'
                return local msg1 "stage `s': raking failed (error `rc')"
                local engstop stop
                continue, break
            }
            // trimming bounds that contradict for some observations (a lower bound above an upper
            // bound): every bound must bind, so the stage fails, at estimation with the error
            // (Estimate) and in a replicate as a failed replicate (Bootstrap)
            if r(tclash_lo) + r(tclash_hi) > 0 {
                local tcl
                if r(tclash_lo) > 0 {
                    local tcl "trimloabs() is above trimhirel() times the base weight of `r(tclash_lo)' focal group `=plural(r(tclash_lo), "observation")'"
                }
                if r(tclash_hi) > 0 {
                    if "`tcl'" != "" local tcl "`tcl', and "
                    local tcl "`tcl'trimlorel() times the base weight of `r(tclash_hi)' focal group `=plural(r(tclash_hi), "observation")' is above trimhiabs()"
                }
                return scalar rc = 198
                return scalar trimclash = 1
                return local msg1 "stage `s': the trimming bounds contradict each other: `tcl'"
                local engstop stop
                continue, break
            }
            local conv  = r(converged)
            local badc  = r(badcontrols)
            local bad0  = r(bad0)
            local mrd   = r(maxreldif)
            local mctrl = r(maxctrl)
            local badtargets = r(badtargets)
            local ntargets = r(ntargets)
            // One worst missed level per control, in original weight units. Control-matrix
            // order includes held margins; map it back to Engine's term numbering.
            if "`misses'" != "" & `badc' > 0 {
                forvalues m = 1/`=rowsof(`RMISS')' {
                    local km = `RMISS'[`m', 1]
                    local cm : word `km' of `mats'
                    local mv : rownames `cm'
                    local mt : list posof "`mv'" in tvlist
                    matrix `MROW' = (`s', `mt', `RMISS'[`m', 2], ///
                        `RMISS'[`m', 3] / `sc', `RMISS'[`m', 4] / `sc', `RMISS'[`m', 5])
                    if `nmissrows' == 0 matrix `misses' = `MROW'
                    else matrix `misses' = `misses' \ `MROW'
                    local ++nmissrows
                }
            }
            local ntrim = r(n_at_bound)
            if missing(`ntrim') local ntrim 0
            local wv `r(worstvar)'
            local wc = r(worstcat)
            local worst_t 0
            forvalues t = 1/`nterms' {
                if "`tv`t''" == "`wv'" local worst_t `t'
            }
            return scalar worst`s'_t = `worst_t'
            return scalar worst`s'_l = `wc'
            // a missed margin is put down to the trimming bounds (r(trimchg#): the error and the
            // notes of Estimate say so) only when trimming changed the raked weights: they are
            // compared with an untrimmed raking of the stage (Rake_mata without the trimming
            // options, which follows the same path when no bound binds) and count as changed
            // when an adjustment factor differs by more than the convergence tolerance
            // a missed level absent from the Reference Group (bad0: a lower bound kept its observations
            // above zero weight) has its own message (r(supmiss#)) and is not compared
            return scalar supmiss`s' = `bad0'
            if `badc' > 0 & `hastrim' & `bad0' == 0 {
                tempvar rwu fd
                capture Rake_mata `ws', touse(`rss') ctotal(`mats') generate(`rwu') `ipfnotrim'
                local rcu = _rc
                // a Break: as after the raking call above
                if `rcu' == 1 {
                    return scalar rc = 1
                    local engstop stop
                    continue, break
                }
                if `rcu' local tchg = (`ntrim' > 0)
                else {
                    qui gen double `fd' = abs(`rwtmp' - `rwu') / `ws' if `rss'
                    qui sum `fd', meanonly
                    local tchg = (r(max) > `tolchg')
                    if "`noi'" != "" {
                        di as txt "{p 0 4 2 78}stage `s': trimming changed the raked weights: " ///
                            cond(`tchg', "yes", "no") " (largest change of a factor " %9.2e r(max) "){p_end}"
                    }
                    drop `fd'
                }
                capture drop `rwu'
            }
            qui replace `rw' = `rwtmp' / `sc' if `rss'
        }
        return scalar trimchg`s' = `tchg'
        if "`checks'" != "" {
            matrix `checks'[`s', 1] = (`conv', `badc', `badtargets', `ntargets', `sc')
        }
        qui replace `rw' = 0 if `rs' & `zero'
        qui replace `rw' = . if !`rs'
        // raked mean and decomposition
        if "`den'" != "" {
            qui sum `den' [iw=`rw'] if `rss', meanonly
            local d1r = r(mean)
            if `d1r' == 0 | missing(`d1r') {
                return scalar rc = 459
                return local msg1 "ratio(): stage `s': the mean of `den' under the raked focal group weights is zero"
                local engstop stop
                continue, break
            }
        }
        qui sum `y' [iw=`rw'] if `rss', meanonly
        local m1r = r(mean)
        if "`den'" != "" {
            local n1r_1 `m1r'
            local m1r = `n1r_1' / `d1r'
        }
        local expl = `m1' - `m1r'
        local unexpl = `m1r' - `m2'
        local pct = cond(`diff' != 0, 100 * `expl' / `diff', .)
        // contrib (individually): a block's raking gives its contribution, its own explained
        // part; the last raking (all blocks) the remainder, explained minus their sum
        if "`contrib'" == "" {
            local incr = `expl' - `prevexpl'
            local prevexpl `expl'
        }
        else if `s' < `nstages' {
            local incr `expl'
            local prevexpl = `prevexpl' + `expl'
        }
        else local incr = `expl' - `prevexpl'
        // weight diagnostics on the raking sample (seldecomp_diag: the raked weights, the
        // factors raked / base, the Kish effective sizes); skipdiag (bootstrap replicates, which
        // use columns 1-7 only) leaves them missing
        tempname DG
        if "`skipdiag'" == "" {
            mata: seldecomp_diag("`rw'", "`w'", "`rss'", `Wrss', "`DG'")
        }
        else matrix `DG' = J(1, 12, .)
        matrix `res'[`s', 1]  = `m1r'
        matrix `res'[`s', 2]  = `expl'
        matrix `res'[`s', 3]  = `unexpl'
        matrix `res'[`s', 4]  = `pct'
        matrix `res'[`s', 5]  = `incr'
        matrix `res'[`s', 6]  = `conv'
        matrix `res'[`s', 7]  = `badc'
        matrix `res'[`s', 8]  = `mrd'
        matrix `res'[`s', 9]  = `mctrl'
        // raked weights: mean, min, max, sd, cv; factors: the same; Kish effective sizes
        forvalues c = 10/19 {
            matrix `res'[`s', `c'] = `DG'[1, `c' - 9]
        }
        matrix `res'[`s', 20] = `Nrake'
        matrix `res'[`s', 21] = `nzero'
        matrix `res'[`s', 22] = `ntrim'
        matrix `res'[`s', 23] = `DG'[1, 11]
        matrix `res'[`s', 24] = `DG'[1, 12]
        matrix `res'[`s', 25] = `DG'[1, 11] / `DG'[1, 12]
        matrix `res'[`s', 26] = `nthin'
        matrix `res'[`s', 27] = cond(`diff' != 0, 100 * `incr' / `diff', .)
        // out-of-support accounting: weight share and mean of the zeroed observations,
        // and the Focal Group mean on the raking sample under the base weights
        if "`skipdiag'" != "" {
        }
        else if `nzero' > 0 {
            matrix `res'[`s', 28] = 100 * (`W1' - `Wrss') / `W1'
            if "`den'" != "" {
                qui sum `den' [iw=`w'] if `rs' & `zero', meanonly
                local dz = r(mean)
                qui sum `den' [iw=`w'] if `rss', meanonly
                local dcs = r(mean)
            }
            qui sum `y' [iw=`w'] if `rs' & `zero', meanonly
            if "`den'" != "" matrix `res'[`s', 29] = cond(`dz' != 0, r(mean) / `dz', .)
            else matrix `res'[`s', 29] = r(mean)
            qui sum `y' [iw=`w'] if `rss', meanonly
            if "`den'" != "" matrix `res'[`s', 30] = cond(`dcs' != 0, r(mean) / `dcs', .)
            else matrix `res'[`s', 30] = r(mean)
        }
        else {
            matrix `res'[`s', 28] = 0
            matrix `res'[`s', 30] = `m1'
        }
        // the further outcomes, from the same raked weights
        forvalues j = 2/`ny' {
            local yv : word `=`j' - 1' of `ylist'
            qui sum `yv' [iw=`rw'] if `rss', meanonly
            local m1r_j = r(mean)
            if "`den'" != "" {
                local n1r_`j' `m1r_j'
                local m1r_j = `n1r_`j'' / `d1r'
            }
            // the same steps (through the same macros) as for the first outcome above
            local expl_j = `m1_`j'' - `m1r_j'
            local unexpl_j = `m1r_j' - `m2_`j''
            local pct_j = cond(`diff_`j'' != 0, 100 * `expl_j' / `diff_`j'', .)
            if "`contrib'" == "" {
                local incr_j = `expl_j' - `prevexpl_`j''
                local prevexpl_`j' `expl_j'
            }
            else if `s' < `nstages' {
                local incr_j `expl_j'
                local prevexpl_`j' = `prevexpl_`j'' + `expl_j'
            }
            else local incr_j = `expl_j' - `prevexpl_`j''
            matrix `X`j''[`s', 1] = `m1r_j'
            matrix `X`j''[`s', 2] = `expl_j'
            matrix `X`j''[`s', 3] = `unexpl_j'
            matrix `X`j''[`s', 4] = `pct_j'
            matrix `X`j''[`s', 5] = `incr_j'
            matrix `X`j''[`s', 27] = cond(`diff_`j'' != 0, 100 * `incr_j' / `diff_`j'', .)
            if `nzero' > 0 & "`skipdiag'" == "" {
                qui sum `yv' [iw=`w'] if `rs' & `zero', meanonly
                if "`den'" != "" matrix `X`j''[`s', 29] = cond(`dz' != 0, r(mean) / `dz', .)
                else matrix `X`j''[`s', 29] = r(mean)
                qui sum `yv' [iw=`w'] if `rss', meanonly
                if "`den'" != "" matrix `X`j''[`s', 30] = cond(`dcs' != 0, r(mean) / `dcs', .)
                else matrix `X`j''[`s', 30] = r(mean)
            }
            else matrix `X`j''[`s', 30] = `m1_`j''
        }
        if `nthin' > 0 {
            local ++nts
            return scalar thin`nts'_s = `s'
            return scalar thin`nts'_n = `nthin'
            return local  thin`nts'_list `thinlist'
        }
        if "`rwlist'" != "" {
            local rwv : word `s' of `rwlist'
            qui replace `rwv' = `rw' if `rs'
        }
        // for the notes (zero-weight lists are returned as data; the caller attaches the
        // labels, and words the notes of convergence and trimming from res())
        if `nzero' > 0 {
            local ++nzs
            return scalar zero`nzs'_s = `s'
            return scalar zero`nzs'_n = `nzero'
            return local  zero`nzs'_list `zerolist'
        }
        // the out-of-support observations held by the lower bounds (incl), for the warning
        else if `incl' & `nsup' > 0 {
            local ++nhs
            return scalar held`nhs'_s = `s'
            return scalar held`nhs'_n = `nsup'
            return local  held`nhs'_list `zerolist'
        }
    }
    // the error paths above leave the stage loop with -continue, break- and exit here: under
    // version 16 a plain -exit- inside forvalues only ends the current iteration ([P] version)
    if "`engstop'" != "" exit
    // further outcomes: the weight diagnostics, which do not depend on the outcome
    forvalues j = 2/`ny' {
        mata: _sd_X = st_matrix("`X`j''"); _sd_R = st_matrix("`res'")
        mata: _sd_X[., (6..26, 28)] = _sd_R[., (6..26, 28)]; st_matrix("`X`j''", _sd_X)
        mata: mata drop _sd_X _sd_R
    }
    return scalar rc = 0
    return scalar nfail = 0
    return scalar nzs = `nzs'
    return scalar nhs = `nhs'
    return scalar nts = `nts'
    return scalar m1 = `m1'
    return scalar m2 = `m2'
    forvalues j = 1/`ny' {
        return scalar m1_`j' = `m1_`j''
        return scalar m2_`j' = `m2_`j''
        if "`den'" != "" {
            return scalar n1_`j' = `n1_`j''
            return scalar n2_`j' = `n2_`j''
            return scalar n1r_`j' = `n1r_`j''
        }
    }
    if "`den'" != "" {
        return scalar d1 = `d1'
        return scalar d2 = `d2'
        return scalar d1r = `d1r'
    }
    return scalar W1 = `W1'
    return scalar W2 = `W2'
end


* ============================================================================
* Rake_mata: Engine's raking, done by seldecomp_rake below (iterative proportional fitting,
*   trimming included): the r() results Engine reads; r(cycles) is the number of cycles run,
*   r(tclash_lo) and r(tclash_hi) the observations whose trimming bounds contradict
* ============================================================================
program define Rake_mata, rclass
    // the trimming bounds are real options (a number with 16 digits: the same bound to the
    // last digit); -1 marks an absent bound
    syntax varname(numeric), touse(varname) ctotal(namelist) generate(name) ///
        [ ITERate(integer 2000) TOLerance(string) CTRLTOLerance(string) trace             ///
          trimhirel(real -1) trimlorel(real -1) trimhiabs(real -1) trimloabs(real -1)    ///
          TRIMFREQuency(string) misses(name) ]
    // the defaults when they are not passed (Engine always passes all three)
    if "`tolerance'" == "" local tolerance 1e-6
    if "`ctrltolerance'" == "" local ctrltolerance 1e-6
    local trace = ("`trace'" != "")
    // when to trim: 0 no bound, 1 once (after the cycles), 2 sometimes (after each cycle;
    // the default), 3 often (after each margin)
    local tfreq 0
    if `trimhirel' != -1 | `trimlorel' != -1 | `trimhiabs' != -1 | `trimloabs' != -1 {
        local tfreq = cond("`trimfrequency'" == "once", 1, cond("`trimfrequency'" == "often", 3, 2))
    }
    tempname R
    qui gen double `generate' = .
    mata: seldecomp_rake("`varlist'", "`touse'", "`ctotal'", `iterate', `tolerance', ///
        `ctrltolerance', (`trimhirel', `trimlorel', `trimhiabs', `trimloabs'), ///
        `tfreq', "`generate'", `trace', "`R'", "`misses'")
    return scalar converged   = `R'[1, 1]
    return scalar badcontrols = `R'[1, 2]
    // with trimfrequency(once) no r(maxreldif) is returned: the trimming after the cycles
    // changes the weights whose largest change the cycles measured
    if `tfreq' != 1 return scalar maxreldif = `R'[1, 3]
    return scalar maxctrl     = `R'[1, 4]
    return local  worstvar `worstvar'
    return scalar worstcat    = `R'[1, 6]
    return scalar cycles      = `R'[1, 7]
    local c 7
    foreach t in hiabs loabs hirel lorel {
        return scalar n_trim`t' = `R'[1, `++c']
    }
    // the observations whose bounds contradict (no raking was done): trimloabs() above
    // trimhirel() times the base weight, trimlorel() times the base weight above trimhiabs()
    return scalar tclash_lo = `R'[1, 12]
    return scalar tclash_hi = `R'[1, 13]
    // the levels absent from the Reference Group (target 0) missed: a lower bound kept weight on
    // their observations
    return scalar bad0 = `R'[1, 14]
    return scalar badtargets = `R'[1, 15]
    return scalar ntargets = `R'[1, 16]
    return scalar n_at_bound = `R'[1, 17]
end


* ============================================================================
* Covariate_detail: shares of each level by group and under the final raked weights
* ============================================================================
program define Covariate_detail
    syntax, touse(varname) g1(varname) w(varname) rw(varname) nterms(integer) ///
        tvlist(string) nlevlist(string) res(name)
    local rows 0
    forvalues t = 1/`nterms' {
        local tv`t'   : word `t' of `tvlist'
        local nlev`t' : word `t' of `nlevlist'
        local rows = `rows' + `nlev`t''
    }
    matrix `res' = J(`rows', 5, .)
    qui sum `w' if `touse' & `g1' == 1, meanonly
    local W1 = r(sum)
    qui sum `w' if `touse' & `g1' == 0, meanonly
    local W2 = r(sum)
    qui sum `rw' if `touse' & `g1' == 1, meanonly
    local WR = r(sum)
    local row 0
    local rnames
    forvalues t = 1/`nterms' {
        forvalues l = 1/`nlev`t'' {
            local ++row
            qui sum `w' if `touse' & `g1' == 1 & `tv`t'' == `l', meanonly
            matrix `res'[`row', 1] = cond(r(N) > 0, r(sum) / `W1', 0)
            qui sum `w' if `touse' & `g1' == 0 & `tv`t'' == `l', meanonly
            matrix `res'[`row', 2] = cond(r(N) > 0, r(sum) / `W2', 0)
            qui sum `rw' if `touse' & `g1' == 1 & `tv`t'' == `l', meanonly
            matrix `res'[`row', 3] = cond(r(N) > 0 & `WR' > 0, r(sum) / `WR', 0)
            qui count if `touse' & `g1' == 1 & `tv`t'' == `l'
            matrix `res'[`row', 4] = r(N)
            qui count if `touse' & `g1' == 0 & `tv`t'' == `l'
            matrix `res'[`row', 5] = r(N)
            local rnames `rnames' `t':`l'
        }
    }
    c_local s_rows `rnames'
end


* ============================================================================
* Draw_weights: Bayesian bootstrap (Dirichlet) replicate weights
*   bw = w * g / gbar, g ~ Gamma((m-2)/m, 1) per unit (or per cluster), m the number of units,
*   gbar = w-weighted mean of g over the sample
*   prep(shp den first): the variables made by Draw_prep, which do not change between
*   replicates (Bootstrap makes them once); without prep() they are made here
* ============================================================================
program define Draw_weights
    syntax, w(varname) bw(varname) touse(varname) obsid(varname) [ cluster(varname) prep(string) ]
    if "`prep'" == "" {
        tempvar shp den first
        local dpopts
        if "`cluster'" != "" local dpopts `dpopts' cluster(`cluster')
        Draw_prep, w(`w') touse(`touse') obsid(`obsid') shp(`shp') den(`den') ///
            first(`first') `dpopts'
    }
    else {
        tokenize `prep'
        local shp `1'
        local den `2'
        local first `3'
    }
    tempvar g
    // sort on a unique key before drawing so that the draws are reproducible from the seed
    if "`cluster'" == "" {
        sort `obsid'
        qui gen double `g' = rgamma(`shp', 1) if `touse'
    }
    else {
        sort `touse' `cluster' `obsid'
        qui gen double `g' = rgamma(`shp', 1) if `first'
        qui by `touse' `cluster': replace `g' = `g'[1] if `touse'
    }
    // bw = w * g * (total of w) / (total of w * g)
    mata: seldecomp_bw("`w'", "`g'", "`den'", "`touse'", "`bw'")
end

* ============================================================================
* Draw_prep: the parts of Draw_weights that do not change between replicates: the shape
*   (m - 2)/m of the Gamma draws (shp), the total of the base weights (den), and the first
*   observation of each sampling unit (first)
* ============================================================================
program define Draw_prep
    syntax, w(varname) touse(varname) obsid(varname) shp(name) den(name) first(name) ///
        [ cluster(varname) ]
    tempvar msz
    // Gamma((m - 2)/m, 1) draws, m = number of sampling units (observations or clusters) in the
    // sample: with Gamma(1,1) the Bayesian bootstrap variance of a mean is only (m - 1)/(m + 1)
    // of its unbiased design-based variance; the shape (m - 2)/m makes the two equal. Bootstrap
    // refuses m < 3.
    if "`cluster'" == "" {
        qui gen byte `first' = `touse'
    }
    else {
        sort `touse' `cluster' `obsid'
        qui by `touse' `cluster': gen byte `first' = (_n == 1) & `touse'
    }
    qui bysort `touse': egen double `msz' = total(`first') if `touse'
    qui gen double `shp' = (`msz' - 2) / `msz' if `touse'
    qui by `touse': egen double `den' = total(`w') if `touse'
end


* ============================================================================
* Bootstrap: replicate loop, SEs, CIs (normal, percentile, BC, BCa)
*   fills v (k x k, zero outside the bootstrapped block) and matrices
*   `bs'_b_bs `bs'_bias `bs'_se `bs'_ci_normal `bs'_ci_percentile `bs'_ci_bc [`bs'_ci_bca `bs'_accel] `bs'_z0
* ============================================================================
program define Bootstrap, rclass
    syntax, touse(varname) g1(varname) w(varname) y(varname) one(varname) nterms(integer)  ///
        tvlist(string) nlevlist(string) nstages(integer) stagedefs(string) b(name) v(name) ///
        bs(name) reps(integer) dots(integer) level(cilevel)                                 ///
        [ ipfopts(string asis) engopts(string asis) tnames(string) pct(string) NODOTS       ///
          holdothers sequential bca ties cluster(varname) seed(string)                         ///
          bsgen(name) bsgenreplace bsfile(string) replace outcomes(string) reverse           ///
          rtvlist(string) rnlevlist(string) rtnames(string) rstagedefs(string)               ///
          recut(string) rcvars(string) rcq(string) rcsample(string) force                   ///
          dstagedefs(string) dends(string)                                              ///
          chain(integer 0) ptouse(string) pg1(string) ptv(string) pnl(string) ]
    // chain(J): J pairs, each with its sample, Focal Group flag, coded terms and
    // their level counts in ptouse() pg1() ptv() pnl() (pair-major); the draws are made on
    // touse(), the chain's sample, and Chain_rep computes every pair and the total from them

    // individually (Engine's contrib): a replicate's stage values are columns 5 and 1 (the
    // contributions and the remainder; the raked means of the blocks' own rakings)
    local contrib = (strpos(" `engopts' ", " contrib ") > 0)
    // outcomes(): all outcomes (Engine gets the further ones through engopts); RB# holds the
    // replicate's stage results of outcome #, RBR# those of the reverse order
    local ny = max(1, `: word count `outcomes'')
    local rev = ("`reverse'" != "" & "`sequential'" != "")
    tempname RB1 RBR1
    local rbl `RB1'
    local rbrl `RBR1'
    forvalues j = 2/`ny' {
        tempname RB`j' RBR`j'
        local rbl `rbl' `RB`j''
        local rbrl `rbrl' `RBR`j''
    }
    local rbx : list rbl - RB1
    local rbrx : list rbrl - RBR1
    // ---- seed ---------------------------------------------------------------
    if `"`seed'"' == "" local seed = runiformint(1, 2147483647)
    set seed `seed'
    local rngstate = c(rngstate)
    // unique observation id in the current order (sortpreserve restores the order at exit)
    tempvar obsid
    qui gen long `obsid' = _n

    // ---- cluster checks -------------------------------------------------------
    local N_clust .
    if "`cluster'" != "" {
        tempvar first
        qui bysort `touse' `cluster': gen byte `first' = (_n == 1) if `touse'
        qui count if `first' == 1
        local N_clust = r(N)
        if `N_clust' < 3 {
            di as err "bayesboot(cluster()): fewer than 3 clusters in the estimation sample"
            exit 459
        }
    }
    else {
        qui count if `touse'
        if r(N) < 3 {
            di as err "bayesboot(): fewer than 3 observations"
            exit 459
        }
    }

    // ---- which stages to replicate ------------------------------------------
    // detail (dstagedefs(), with the stages): the term chain replaces the block stages, whose
    // values are its rows dends() (the same rakings), and its increments are replicated too
    local det = ("`dstagedefs'" != "" & "`sequential'" != "")
    local k = colsof(`b')
    if "`sequential'" != "" {
        local nsb `nstages'
        local sdb `stagedefs'
        local kb = `ny' * (6 + 4 * `nstages' * (1 + `rev'))
        if `contrib' local kb = `ny' * (6 + 2 * `nstages' - 1)
        if `det' {
            local nsb `nterms'
            local sdb `dstagedefs'
            local kb = `kb' + `ny' * `nterms'
        }
    }
    else {
        local nsb 1
        local rest `stagedefs'
        while `"`rest'"' != "" {
            gettoken piece rest : rest, parse("|")
            if "`piece'" != "|" local sdb `piece'
        }
        local kb = 6 * `ny'
    }
    // names and positions in e(b), in the order of Rep_vals: per outcome (e(b) has k / ny
    // columns each) the overall six, the marginal explained parts, the cumulative explained
    // and unexplained parts and the raked means by stage, the terms' marginal explained parts
    // (detail), and the stage values in reverse order; outcome # > 1 is prefixed o#_ in the
    // saving() file
    local names : colfullnames `b'
    local kper = `k' / `ny'
    local vnames
    local pos
    local pnames
    forvalues j = 1/`ny' {
        local vp = cond(`j' == 1, "", "o`j'_")
        local off = (`j' - 1) * `kper'
        foreach vv in focal_group raked_focal ref_group difference explained unexplained {
            local vnames `vnames' `vp'`vv'
        }
        forvalues p = 1/6 {
            local pos `pos' `=`off' + `p''
        }
        if "`sequential'" != "" {
            local c 6
            foreach vv in marg_s cumexpl_s cumunexpl_s raked_s det_t rmarg_s rcumexpl_s ///
                rcumunexpl_s rraked_s {
                if !`rev' & inlist("`vv'", "rmarg_s", "rcumexpl_s", "rcumunexpl_s", "rraked_s") continue
                if `contrib' & !inlist("`vv'", "marg_s", "raked_s") continue
                if !`det' & "`vv'" == "det_t" continue
                local nv = cond("`vv'" == "det_t", `nterms', `nstages')
                if `contrib' & "`vv'" == "raked_s" local nv = `nstages' - 1
                forvalues s = 1/`nv' {
                    local vnames `vnames' `vp'`vv'`s'
                    local pos `pos' `=`off' + `c' + `s''
                }
                local c = `c' + `nv'
            }
        }
        local pnames `pnames' `vp'pct_explained
    }
    // chain: each pair's block of e(b) (the six, then with a staged mode its marginal explained
    // and cumulative parts and raked means by stage), replicated as an outcome's is, prefixed
    // p#_; then the total (difference, explained, unexplained, and the summed marginal explained
    // parts), prefixed t_; then the
    // explained share of each pair and of the total (npc of them). OBC: the observed row
    local npc `ny'
    if `chain' > 0 {
        local sq = ("`sequential'" != "")
        local kpair = 6 + 4 * `nstages' * (`k' > 6 * `chain' + 3)
        local vnames
        local pos
        local pnames
        forvalues j = 1/`chain' {
            local off = (`j' - 1) * `kpair'
            foreach vv in focal_group raked_focal ref_group difference explained unexplained {
                local vnames `vnames' p`j'_`vv'
            }
            forvalues p = 1/6 {
                local pos `pos' `=`off' + `p''
            }
            if `sq' {
                local c 6
                foreach vv in marg_s cumexpl_s cumunexpl_s raked_s {
                    forvalues s = 1/`nstages' {
                        local vnames `vnames' p`j'_`vv'`s'
                        local pos `pos' `=`off' + `c' + `s''
                    }
                    local c = `c' + `nstages'
                }
            }
            local pnames `pnames' p`j'_pct_explained
        }
        local off = `chain' * `kpair'
        foreach vv in difference explained unexplained {
            local vnames `vnames' t_`vv'
        }
        forvalues p = 1/3 {
            local pos `pos' `=`off' + `p''
        }
        if `sq' {
            forvalues s = 1/`nstages' {
                local vnames `vnames' t_marg_s`s'
                local pos `pos' `=`off' + 3 + `s''
            }
        }
        local pnames `pnames' t_pct_explained
        local kb : word count `pos'
        local npc = `chain' + 1
        tempname OBC CR
        matrix `OBC' = J(1, `kb' + `npc', .)
        local j 0
        foreach p of local pos {
            matrix `OBC'[1, `++j'] = `b'[1, `p']
        }
        foreach p of local pct {
            matrix `OBC'[1, `++j'] = `p'
        }
        local cvals
        forvalues j = 1/`=`kb' + `npc'' {
            local cvals `cvals' (`CR'[1, `j'])
        }
        local crop npairs(`chain') ptouse(`ptouse') pg1(`pg1') ptv(`ptv') pnl(`pnl') y(`y') ///
            one(`one') nterms(`nterms') nstages(`nstages') nsb(`nsb') sdb(`sdb') res(`RB1')     ///
            row(`CR') tnames(`tnames') `holdothers' stages(`=`sq' * `nstages'') `force'
        local fpairs
    }
    local bsnames
    foreach p of local pos {
        local bsnames `bsnames' `: word `p' of `names''
    }

    qui count if `touse'
    local Nobs = r(N)

    // ---- replicate loop -----------------------------------------------------
    local dwopts
    if "`cluster'" != "" local dwopts `dwopts' cluster(`cluster')
    // the parts of the draws that do not change between replicates, made once; the replicates
    // and jackknife deletions skip Engine's weight diagnostics (columns 10-30 are not used)
    tempvar dshp dden dfirst
    Draw_prep, w(`w') touse(`touse') obsid(`obsid') shp(`dshp') den(`dden') ///
        first(`dfirst') `dwopts'
    local dwopts `dwopts' prep(`dshp' `dden' `dfirst')
    local engopts `engopts' skipdiag
    tempvar bw
    qui gen double `bw' = .
    // recut: every replicate (and BCa deletion) bins the c. variables again, from its own
    // weights by the rule of contvars(), into rcb#; each term that uses them is coded again
    // into rtv# (Recut), which replaces its variable in the term lists given to Engine, and
    // its level count in nll/rnll. Without recut the lists are those of the point estimate
    local tvl `tvlist'
    local nll `nlevlist'
    local rtvl `rtvlist'
    local rnll `rnlevlist'
    if `"`recut'"' != "" {
        local rcbins
        foreach rcv of local rcvars {
            tempvar rcb
            qui gen int `rcb' = .
            local rcbins `rcbins' `rcb'
        }
        local rcspec
        local rct
        local rest `recut'
        while `"`rest'"' != "" {
            gettoken piece rest : rest, parse("|")
            if `"`piece'"' == "|" continue
            gettoken t comps : piece
            local cc
            foreach c of local comps {
                if substr("`c'", 1, 2) == "c." {
                    local i : list posof "`=substr("`c'", 3, .)'" in rcvars
                    local cc `cc' `: word `i' of `rcbins''
                }
                else local cc `cc' `c'
            }
            tempvar rtv`t'
            qui gen long `rtv`t'' = .
            if "`rcspec'" != "" local rcspec `rcspec' |
            local rcspec `rcspec' `t' `rtv`t'' `cc'
            local rct `rct' `t'
        }
        local tvl
        local j 0
        foreach rcv of local tvlist {
            local ++j
            if `: list j in rct' local tvl `tvl' `rtv`j''
            else local tvl `tvl' `rcv'
        }
        // the reverse order: the position of each of its terms in the forward order
        local rtix
        local rtvl
        foreach rcv of local rtvlist {
            local j : list posof "`rcv'" in tvlist
            local rtix `rtix' `j'
            local rtvl `rtvl' `: word `j' of `tvl''
        }
        local rcopts touse(`touse') g1(`g1') sample(`rcsample') vars(`rcvars') q(`rcq') ///
            bins(`rcbins') terms(`rcspec') nlevlist(`nlevlist') rtix(`rtix')
    }
    tempname ph
    local RB `RB1'
    tempfile bsdta
    postfile `ph' double(`vnames' `pnames') using `"`bsdta'"'
    local rvopts nsb(`nsb')
    if "`sequential'" != "" local rvopts `rvopts' stages(`nstages')
    if `rev' local rvopts `rvopts' rres(`rbrl')
    if `contrib' local rvopts `rvopts' contrib
    if `det' local rvopts `rvopts' dends(`dends')
    if "`nodots'" == "" {
        di as txt _n "Bayesian bootstrap replications (" as res `reps' as txt ")"
        di as txt "{hline 4}{c +}{hline 3} 1 {hline 3}{c +}{hline 3} 2 {hline 3}{c +}{hline 3} 3 " ///
                  "{hline 3}{c +}{hline 3} 4 {hline 3}{c +}{hline 3} 5 "
    }
    local nfail 0
    local niter 0
    local ndots 0
    // the causes of the failed replications: the raking stopped at the iterate() limit short
    // of its targets (it), missed them otherwise (tg), the trimming bounds contradicted each
    // other on the replicate's weights (tb), or the raking failed (rk); itst, the stages at which
    // replications stopped at that limit (numbered as the point estimate's), itrev 1 when the
    // reverse order did
    local nf_it 0
    local nf_tg 0
    local nf_rk 0
    local nf_tb 0
    // force: the replicates kept although a margin was missed (nforce; e(N_forced))
    local nforce 0
    local itst
    local itrev 0
    // detail: a term stage that is not the last of its block has no number of the point
    // estimate's stages; it is marked as the reverse order is (the note is then not folded)
    local itdet 0
    // Track only variables actually written: an interruption can leave later names unused,
    // including pre-existing variables that generate(..., replace) has not reached yet.
    local bscreated
    forvalues r = 1/`reps' {
        Draw_weights, w(`w') bw(`bw') touse(`touse') obsid(`obsid') `dwopts'
        if "`bsgen'" != "" {
            // (the exact name: with varabbrev on, drop takes an abbreviation, another variable)
            if "`bsgenreplace'" != "" {
                capture confirm variable `bsgen'_`r', exact
                if !_rc drop `bsgen'_`r'
            }
            qui gen double `bsgen'_`r' = `bw' if `touse'
            local bscreated `bscreated' `bsgen'_`r'
            label variable `bsgen'_`r' "Bayesian bootstrap weights, replicate `r'"
        }
        local rc 0
        local tbc 0
        if `"`recut'"' != "" {
            Recut, w(`bw') `rcopts'
            local nll `s_nll'
            local rnll `s_rnll'
            local rc `s_rcrc'
        }
        // chain: every pair on the replicate's weights, and the total (Chain_rep)
        if `chain' > 0 {
            Chain_rep, w(`bw') `crop' ipfopts(`ipfopts') engopts(`engopts')
            local rc `ch_rc'
        }
        else if `rc' == 0 {
            Engine, touse(`touse') g1(`g1') w(`bw') y(`y') one(`one') nterms(`nterms')   ///
                tvlist(`tvl') nlevlist(`nll') nstages(`nsb') stagedefs(`sdb') res(`RB') ///
                ipfopts(`ipfopts') tnames(`tnames') `engopts' `holdothers' resx(`rbx')
            local rc = r(rc)
            local tbc = (r(trimclash) == 1)
        }
        if `rc' == 1 {                       // Break inside Engine: stop, do not count it as a failed replicate
            local bsbreak 1
            continue, break
        }
        // a replicate fails when raking fails or a margin is missed, as the point estimate
        // would give 430 (force keeps it: frc); one whose raking stopped before the weights
        // converged (the iterate() limit) with every margin met is kept, and counted in niter
        // (e(N_iterlimit)). A margin that trimming left unmet (Engine's r(trimchg#)) is the
        // trimming's, as in the point estimate's message, not the limit's (itb)
        local itl 0
        local frc 0
        local why rk
        if `chain' == 0 & `rc' {
            if `tbc' local why tb
        }
        if `chain' > 0 {
            // (Chain_rep applied that rule to every pair; the pair that failed, for the note)
            local itl `ch_itl'
            local frc `ch_forced'
            if `rc' local why `ch_why'
            local itst : list itst | ch_itst
            if `rc' & !`: list ch_fpair in fpairs' local fpairs `fpairs' `ch_fpair'
        }
        else if `rc' == 0 {
            forvalues s = 1/`nsb' {
                local st = cond(`nsb' == 1, `nstages', `s')
                if `det' {
                    local st : list posof "`s'" in dends
                    if `st' == 0 & `RB'[`s', 6] == 0 local itdet 1
                }
                if `RB'[`s', 7] > 0 & "`force'" != "" local frc 1
                else if `RB'[`s', 7] > 0 {
                    local itb = (`RB'[`s', 6] == 0 & r(trimchg`s') != 1)
                    if `rc' == 0 local why = cond(`itb', "it", "tg")
                    if `itb' & `st' > 0 & !`: list st in itst' local itst `itst' `st'
                    local rc 430
                }
                else if `RB'[`s', 6] == 0 & `RB'[`s', 7] == 0 {
                    local itl 1
                    if `st' > 0 & !`: list st in itst' local itst `itst' `st'
                }
            }
        }
        if `rc' == 0 & `chain' == 0 {
            local m1l
            local m2l
            forvalues j = 1/`ny' {
                local m1 = r(m1_`j')
                local m2 = r(m2_`j')
                local m1l `m1l' `m1'
                local m2l `m2l' `m2'
            }
            // the reverse order on the same draws
            if `rev' {
                Engine, touse(`touse') g1(`g1') w(`bw') y(`y') one(`one') nterms(`nterms')  ///
                    tvlist(`rtvl') nlevlist(`rnll') nstages(`nstages')             ///
                    stagedefs(`rstagedefs') res(`RBR1') ipfopts(`ipfopts') tnames(`rtnames') ///
                    `engopts' `holdothers' resx(`rbrx')
                local rc = r(rc)
                if `rc' == 1 {
                    local bsbreak 1
                    continue, break
                }
                if `rc' & r(trimclash) == 1 local why tb
                if `rc' == 0 {
                    forvalues s = 1/`nstages' {
                        if `RBR1'[`s', 7] > 0 & "`force'" != "" local frc 1
                        else if `RBR1'[`s', 7] > 0 {
                            local itb = (`RBR1'[`s', 6] == 0 & r(trimchg`s') != 1)
                            if `rc' == 0 local why = cond(`itb', "it", "tg")
                            if `itb' local itrev 1
                            local rc 430
                        }
                        else if `RBR1'[`s', 6] == 0 & `RBR1'[`s', 7] == 0 {
                            local itl 1
                            local itrev 1
                        }
                    }
                }
            }
        }
        if `rc' == 0 {
            if `chain' > 0 local vals `cvals'
            else Rep_vals, res(`rbl') m1(`m1l') m2(`m2l') `rvopts'
            post `ph' `vals'
            if `itl' local ++niter
            if `frc' local ++nforce
        }
        else {
            local ++nfail
            local ++nf_`why'
            // a row of missing values, as -bootstrap- writes for a failed replicate, so that
            // -bstat using- the saving() file reports e(N_misreps)
            local nvn : word count `vnames' `pnames'
            local vals
            forvalues j = 1/`nvn' {
                local vals `vals' (.)
            }
            post `ph' `vals'
        }
        if "`nodots'" == "" & mod(`r', `dots') == 0 {
            local ++ndots
            if `rc' == 0 di as txt "." _c
            else di as err "x" _c
            if mod(`ndots', 50) == 0 di as txt " " %5.0f `r'
        }
    }
    if "`nodots'" == "" & mod(`ndots', 50) != 0 di ""
    if "`bsbreak'" != "" {
        // do not leave partial bayesboot(generate()) variables behind
        foreach gn of local bscreated {
            capture confirm variable `gn', exact
            if !_rc drop `gn'
        }
        exit 1
    }
    postclose `ph'
    // what failed the replications, and what to do: given by
    // Estimate in one note, or in the note of the point estimate's iterate() limit; here in
    // the error when fewer than 2 succeeded. With force, a replication that missed a margin is
    // kept (nforce), so these causes are then the others
    local stoptxt "at the iterate() limit"
    // (a limit reached short of the targets can be too few iterations or margins that cannot
    // be met, as for the point estimate's error, Miss_error: both remedies)
    local bswhy
    local bsfix
    if `nf_it' > 0 & `nf_tg' == 0 local bswhy "raking stopped at iterate() short of its targets"
    else if `nf_tg' > 0 {
        local bswhy "raking did not reach its targets"
        if `nf_it' > 0 local bswhy "`bswhy' (`nf_it' stopped at iterate())"
    }
    if `nf_tb' > 0 {
        local bswhy = "`bswhy'" + cond("`bswhy'" != "", ", or ", "") + ///
            "the trimming bounds contradicted each other"
    }
    if `nf_rk' > 0 local bswhy = "`bswhy'" + cond("`bswhy'" != "", ", or ", "") + "raking failed"
    // chain: the pairs whose raking failed them
    if `chain' > 0 & "`fpairs'" != "" {
        numlist "`fpairs'", sort
        local fpl
        foreach f in `r(numlist)' {
            local fpl = "`fpl'" + cond("`fpl'" != "", ", ", "") + "p`f'"
        }
        local bswhy "`bswhy' (`=plural(`: word count `fpairs'', "pair")' `fpl')"
    }
    local N_reps = `reps' - `nfail'
    if `N_reps' < 2 {
        foreach gn of local bscreated {
            capture confirm variable `gn', exact
            if !_rc drop `gn'
        }
        di as err "{p 0 4 2}bayesboot(): `N_reps' of `reps' replications succeeded; 2 or more are required. Failures: `bswhy'{p_end}"
        exit 430
    }

    // ---- jackknife acceleration for BCa -------------------------------------
    if "`bca'" != "" {
        tempvar unit jw
        sort `obsid'
        if "`cluster'" == "" {
            qui gen long `unit' = sum(`touse')
            qui replace `unit' = . if !`touse'
        }
        else qui egen long `unit' = group(`cluster') if `touse'
        qui sum `unit' if `touse', meanonly
        local nunits = r(max)
        tempname ACC R
        mata: seldecomp_JK = J(`nunits', `kb' + `npc', .)
        qui gen double `jw' = .
        // chain: the pairs that contain each unit; a deletion recomputes only those, the other
        // pairs keep their point estimates (Chain_rep's only() and obs())
        if `chain' > 0 {
            mata: _sd_INC = J(`nunits', `chain', 0)
            forvalues j = 1/`chain' {
                local ptj : word `j' of `ptouse'
                mata: _sd_u = st_data(., "`unit'", "`ptj'"); _sd_INC[_sd_u, `j'] = J(rows(_sd_u), 1, 1)
            }
        }
        if "`nodots'" == "" {
            di as txt _n "Jackknife for BCa acceleration (" as res `nunits' as txt " deletions)"
            di as txt "{hline 4}{c +}{hline 3} 1 {hline 3}{c +}{hline 3} 2 {hline 3}{c +}{hline 3} 3 " ///
                      "{hline 3}{c +}{hline 3} 4 {hline 3}{c +}{hline 3} 5 "
        }
        local njfail 0
        local njiter 0
        local njforce 0
        local jdots 0
        forvalues u = 1/`nunits' {
            qui replace `jw' = cond(`unit' == `u', 0, `w') if `touse'
            local rc 0
            if `"`recut'"' != "" {
                Recut, w(`jw') `rcopts'
                local nll `s_nll'
                local rnll `s_rnll'
                local rc `s_rcrc'
            }
            if `chain' > 0 {
                mata: st_local("tp", invtokens(strofreal(selectindex(_sd_INC[`u', .]))))
                Chain_rep, w(`jw') `crop' ipfopts(`ipfopts') engopts(`engopts') only(`tp') obs(`OBC')
                local rc `ch_rc'
            }
            else if `rc' == 0 {
                Engine, touse(`touse') g1(`g1') w(`jw') y(`y') one(`one') nterms(`nterms')   ///
                    tvlist(`tvl') nlevlist(`nll') nstages(`nsb') stagedefs(`sdb') res(`RB') ///
                    ipfopts(`ipfopts') tnames(`tnames') `engopts' `holdothers' resx(`rbx')
                local rc = r(rc)
            }
            if `rc' == 1 {                   // Break inside Engine: stop the jackknife
                local jkbreak 1
                continue, break
            }
            // the rule of the replicates: a missed margin fails the deletion (unless force keeps
            // it: njforce), the iterate() limit alone does not (counted in njiter)
            local itl 0
            local frc 0
            if `chain' > 0 {
                local itl `ch_itl'
                local frc `ch_forced'
            }
            else if `rc' == 0 {
                forvalues s = 1/`nsb' {
                    if `RB'[`s', 7] > 0 & "`force'" != "" local frc 1
                    else if `RB'[`s', 7] > 0 local rc 430
                    else if `RB'[`s', 6] == 0 local itl 1
                }
            }
            if `rc' == 0 & `chain' == 0 {
                local m1l
                local m2l
                forvalues j = 1/`ny' {
                    local m1 = r(m1_`j')
                    local m2 = r(m2_`j')
                    local m1l `m1l' `m1'
                    local m2l `m2l' `m2'
                }
                if `rev' {
                    Engine, touse(`touse') g1(`g1') w(`jw') y(`y') one(`one') nterms(`nterms') ///
                        tvlist(`rtvl') nlevlist(`rnll') nstages(`nstages')            ///
                        stagedefs(`rstagedefs') res(`RBR1') ipfopts(`ipfopts')                ///
                        tnames(`rtnames') `engopts' `holdothers' resx(`rbrx')
                    local rc = r(rc)
                    if `rc' == 1 {
                        local jkbreak 1
                        continue, break
                    }
                    if `rc' == 0 {
                        forvalues s = 1/`nstages' {
                            if `RBR1'[`s', 7] > 0 & "`force'" != "" local frc 1
                            else if `RBR1'[`s', 7] > 0 local rc 430
                            else if `RBR1'[`s', 6] == 0 local itl 1
                        }
                    }
                }
            }
            if `rc' == 0 {
                // the replicate's row, in the order of the replicate loop (Rep_vals)
                if `chain' > 0 local vals `cvals'
                else Rep_vals, res(`rbl') m1(`m1l') m2(`m2l') `rvopts'
                matrix `R' = J(1, `kb' + `npc', .)
                local c 0
                while `"`vals'"' != "" {
                    gettoken e vals : vals, match(par)
                    matrix `R'[1, `++c'] = `e'
                }
                mata: seldecomp_JK[`u', .] = st_matrix("`R'")
                if `itl' local ++njiter
                if `frc' local ++njforce
            }
            else local ++njfail
            if "`nodots'" == "" & mod(`u', `dots') == 0 {
                local ++jdots
                if `rc' == 0 di as txt "." _c
                else di as err "x" _c
                if mod(`jdots', 50) == 0 di as txt " " %5.0f `u'
            }
        }
        if "`nodots'" == "" & mod(`jdots', 50) != 0 di ""
        if "`jkbreak'" != "" {
            capture mata: mata drop seldecomp_JK
            capture mata: mata drop _sd_INC _sd_u
            // nor the bayesboot(generate()) variables made by the replicate loop
            foreach gn of local bscreated {
                capture confirm variable `gn', exact
                if !_rc drop `gn'
            }
            exit 1
        }
        // Keep the acceleration from the deletions that ran. Return the failed-deletion count
        // for the caution at estimation and replay.
        mata: st_matrix("`ACC'", seldecomp_accel(seldecomp_JK))
        mata: mata drop seldecomp_JK
        capture mata: mata drop _sd_INC _sd_u
    }

    // ---- summarize the replicates -------------------------------------------
    tempname C OB
    matrix `v' = J(`k', `k', 0)
    if "`bca'" != "" matrix `bs'_accel = `ACC'[1, 1..`kb']
    matrix `OB' = J(1, `kb', .)
    local j 0
    foreach p of local pos {
        local ++j
        matrix `OB'[1, `j'] = `b'[1, `p']
    }
    preserve
    qui use `"`bsdta'"', clear
    qui corr `vnames', covariance
    matrix `C' = r(C)
    // explained share: SE and percentile CI (not part of e(b)); none when the observed share is
    // undefined (a difference of 0); pct_se_# etc. for outcome # (from Bootstrap_ci)
    tempname PO
    matrix `PO' = J(1, `npc', .)
    forvalues j = 1/`npc' {
        matrix `PO'[1, `j'] = `: word `j' of `pct''
        local pct_se_`j' .
        local pct_ll_`j' .
        local pct_ul_`j' .
    }
    // the t intervals: the observations, or the clusters, less 1 degrees of freedom
    local df_bs = cond("`cluster'" != "", `N_clust', `Nobs') - 1
    local ciopts `ties' pct(`pnames') pobs(`PO') df(`df_bs')
    if "`bca'" != "" local ciopts `ciopts' accel(`bs'_accel)
    Bootstrap_ci `vnames', obs(`OB') level(`level') nreps(`N_reps') prefix(`bs') `ciopts'
    local pct_se `pct_se_1'
    local pct_ll `pct_ll_1'
    local pct_ul `pct_ul_1'
    matrix `bs'_reps = J(1, `kb', `N_reps')
    // covariance block into V
    forvalues j = 1/`kb' {
        local pj : word `j' of `pos'
        forvalues i = 1/`kb' {
            local pi : word `i' of `pos'
            matrix `v'[`pi', `pj'] = `C'[`i', `j']
        }
    }
    // the replicates (a row of missing values for a failed one), stored in e(bsdraws) and
    // e(bsdraws_pct) so that a replay at another level recomputes the intervals from them;
    // more than c(max_matdim) are not stored, which the replay at another level says
    // (Display_bootci)
    if `reps' <= c(max_matdim) {
        mkmat `vnames', matrix(`bs'_bsdraws)
        // outcomes(): one column per outcome (chain: per pair and the total)
        if `npc' > 1 mkmat `pnames', matrix(`bs'_bsdraws_pct)
        else if !missing(`pct') mkmat pct_explained, matrix(`bs'_bsdraws_pct)
    }
    // saving(): a version-3 bootstrap dataset readable by -bstat using-
    if `"`bsfile'"' != "" {
        local j 0
        foreach vn of local vnames {
            local ++j
            local full : word `j' of `bsnames'
            gettoken eq nm : full, parse(":")
            gettoken colon nm : nm, parse(":")
            local pj : word `j' of `pos'
            label variable `vn' "`full'"
            char `vn'[observed] `=`b'[1, `pj']'
            char `vn'[expression] "e(b)[`full']"
            char `vn'[coleq] `eq'
            char `vn'[colname] `nm'
            if "`bca'" != "" {
                char `vn'[acceleration] `=`ACC'[1, `j']'
            }
        }
        // bstat needs an observed value for every variable: no share column when it is undefined
        forvalues j = 1/`npc' {
            local pn : word `j' of `pnames'
            local pj : word `j' of `pct'
            if missing(`pj') drop `pn'
            else {
                label variable `pn' "explained share of the difference (%)"
                char `pn'[observed] `pj'
                if `ny' == 1 & `chain' == 0 {
                    char `pn'[expression] "100*e(explained)/e(difference)"
                    char `pn'[coleq] overall
                }
                else {
                    local yo = cond(`chain' > 0, cond(`j' <= `chain', "p`j'", "total"), word("`outcomes'", `j'))
                    char `pn'[expression] "100*_b[`yo':explained]/_b[`yo':difference]"
                    char `pn'[coleq] `yo'
                }
                char `pn'[colname] pct_explained
                if "`bca'" != "" {
                    char `pn'[acceleration] `=`ACC'[1, `kb' + `j']'
                }
            }
        }
        char _dta[bs_version] 3
        char _dta[N] `Nobs'
        char _dta[N_reps] `N_reps'
        char _dta[command] seldecomp
        char _dta[seed] `seed'
        char _dta[rngstate] `"`rngstate'"'
        if "`cluster'" != "" {
            char _dta[N_cluster] `N_clust'
            char _dta[cluster] `cluster'
        }
        label data "seldecomp: Bayesian bootstrap replicates"
        capture noisily qui save `"`bsfile'"', `replace'
        if _rc {
            // the file could not be written: no bayesboot(generate()) variables left behind
            local rc = _rc
            restore
            foreach gn of local bscreated {
                capture confirm variable `gn', exact
                if !_rc drop `gn'
            }
            exit `rc'
        }
        di as txt `"(replicate results saved in `bsfile')"'
    }
    restore

    // stripes
    matrix rownames `v' = `names'
    matrix colnames `v' = `names'
    foreach m in b_bs bias se z0 ci_normal ci_t ci_percentile ci_bc reps {
        matrix colnames `bs'_`m' = `bsnames'
    }
    if "`bca'" != "" {
        matrix colnames `bs'_ci_bca = `bsnames'
        matrix colnames `bs'_accel  = `bsnames'
        matrix rownames `bs'_ci_bca = ll ul
    }
    matrix rownames `bs'_ci_normal     = ll ul
    matrix rownames `bs'_ci_t          = ll ul
    matrix rownames `bs'_ci_percentile = ll ul
    matrix rownames `bs'_ci_bc         = ll ul
    capture confirm matrix `bs'_bsdraws
    if !_rc matrix colnames `bs'_bsdraws = `bsnames'
    capture confirm matrix `bs'_bsdraws_pct
    if !_rc {
        local pcn pct_explained
        if `ny' > 1 {
            local pcn
            foreach yo of local outcomes {
                local pcn `pcn' `yo':pct_explained
            }
        }
        if `chain' > 0 {
            local pcn
            forvalues j = 1/`chain' {
                local pcn `pcn' p`j':pct_explained
            }
            local pcn `pcn' total:pct_explained
        }
        matrix colnames `bs'_bsdraws_pct = `pcn'
    }

    return scalar k = `kb'
    return scalar seed = `seed'
    return scalar N_reps = `N_reps'
    return scalar N_fail = `nfail'
    return scalar N_iter = `niter'
    return scalar N_force = `nforce'
    return scalar N_jkforce = cond("`bca'" != "", 0`njforce', 0)
    return scalar N_jkfail = cond("`bca'" != "", 0`njfail', 0)
    // for the note on the replications (Estimate): the failures at the iterate() limit, the
    // stages that stopped there, the causes and what to do
    return scalar N_failit = `nf_it'
    return scalar itrev = (`itrev' | `itdet')
    return local itstages `itst'
    return local bswhy `"`bswhy'"'
    return local bsfix `"`bsfix'"'
    return scalar pct_se = `pct_se'
    return scalar pct_ll = `pct_ll'
    return scalar pct_ul = `pct_ul'
    forvalues j = 1/`npc' {
        return scalar pct_se_`j' = `pct_se_`j''
        return scalar pct_ll_`j' = `pct_ll_`j''
        return scalar pct_ul_`j' = `pct_ul_`j''
    }
    return scalar N_clust = `N_clust'
    return scalar df_bs = `df_bs'
    return local rngstate `"`rngstate'"'
    return local bsnames `bsnames'
end

* ============================================================================
* Bootstrap_ci: summaries of the replicates in memory (one variable per statistic, missing in
*   a failed replicate) at level(): `prefix'_b_bs _bias _se _z0 _ci_normal _ci_t (with df(),
*   its degrees of freedom) _ci_percentile _ci_bc [_ci_bca with accel()]; with pct(), the explained share's SE and percentile CI in
*   the caller's locals pct_se pct_ll pct_ul, and for the #th variable of pct() (outcomes())
*   in pct_se_# pct_ll_# pct_ul_#, skipped where pobs() (the observed shares) is missing.
*   Used by Bootstrap at estimation and by a replay at another level (Display_ci), so that
*   both compute the intervals in the same way.
* ============================================================================
program define Bootstrap_ci
    syntax varlist, obs(name) level(cilevel) nreps(integer) prefix(name) ///
        [ ties accel(name) pct(varlist) pobs(name) df(real -1) ]
    local bs `prefix'
    local N_reps `nreps'
    local kb : word count `varlist'
    tempname ob
    local zz = invnormal(1 - (100 - `level') / 200)
    local tt = cond(`df' > 0, invttail(`df', (100 - `level') / 200), .)
    local alo = (100 - `level') / 200
    local ahi = 1 - `alo'
    matrix `bs'_b_bs = J(1, `kb', .)
    matrix `bs'_bias = J(1, `kb', .)
    matrix `bs'_se   = J(1, `kb', .)
    matrix `bs'_z0   = J(1, `kb', .)
    matrix `bs'_ci_normal     = J(2, `kb', .)
    matrix `bs'_ci_t          = J(2, `kb', .)
    matrix `bs'_ci_percentile = J(2, `kb', .)
    matrix `bs'_ci_bc         = J(2, `kb', .)
    if "`accel'" != "" matrix `bs'_ci_bca = J(2, `kb', .)
    local k 0
    foreach pv of local pct {
        local ++k
        if "`pobs'" != "" {
            if missing(`pobs'[1, `k']) continue
        }
        qui sum `pv'
        local pct_se = r(sd)
        qui _pctile `pv', p(`=100*`alo'' `=100*`ahi'')
        local pct_ll = r(r1)
        local pct_ul = r(r2)
        c_local pct_se_`k' `pct_se'
        c_local pct_ll_`k' `pct_ll'
        c_local pct_ul_`k' `pct_ul'
        if `k' == 1 {
            c_local pct_se `pct_se'
            c_local pct_ll `pct_ll'
            c_local pct_ul `pct_ul'
        }
    }
    local j 0
    foreach vn of local varlist {
        local ++j
        scalar `ob' = `obs'[1, `j']
        qui sum `vn'
        local se = r(sd)
        matrix `bs'_b_bs[1, `j'] = r(mean)
        matrix `bs'_bias[1, `j'] = r(mean) - `ob'
        matrix `bs'_se[1, `j']   = `se'
        matrix `bs'_ci_normal[1, `j'] = `ob' - `zz' * `se'
        matrix `bs'_ci_normal[2, `j'] = `ob' + `zz' * `se'
        matrix `bs'_ci_t[1, `j'] = `ob' - `tt' * `se'
        matrix `bs'_ci_t[2, `j'] = `ob' + `tt' * `se'
        qui _pctile `vn', p(`=100*`alo'' `=100*`ahi'')
        matrix `bs'_ci_percentile[1, `j'] = r(r1)
        matrix `bs'_ci_percentile[2, `j'] = r(r2)
        // bias correction as in Stata's bstat: count replicates <= observed, minus half the ties
        qui count if `vn' <= `ob'
        local nle = r(N)
        qui count if `vn' == `ob'
        local neq = r(N)
        if "`ties'" != "" local p0 = (`nle' - 0.5 * `neq') / `N_reps'
        else local p0 = `nle' / `N_reps'
        if `p0' > 0 & `p0' < 1 {
            local z0 = invnormal(`p0')
            matrix `bs'_z0[1, `j'] = `z0'
            local pl = normal(2 * `z0' - `zz')
            local pu = normal(2 * `z0' + `zz')
            // _pctile rejects p(0) and p(100), which normal() reaches for extreme z0: leave missing
            if `pl' > 0 & `pl' < 1 & `pu' > 0 & `pu' < 1 & `pl' < `pu' {
                qui _pctile `vn', p(`=100*`pl'' `=100*`pu'')
                matrix `bs'_ci_bc[1, `j'] = r(r1)
                matrix `bs'_ci_bc[2, `j'] = r(r2)
            }
            if "`accel'" != "" {
                local a = `accel'[1, `j']
                if `a' < . {
                    local pl = normal(`z0' + (`z0' - `zz') / (1 - `a' * (`z0' - `zz')))
                    local pu = normal(`z0' + (`z0' + `zz') / (1 - `a' * (`z0' + `zz')))
                    if `pl' > 0 & `pl' < 1 & `pu' > 0 & `pu' < 1 & `pl' < `pu' {
                        qui _pctile `vn', p(`=100*`pl'' `=100*`pu'')
                        matrix `bs'_ci_bca[1, `j'] = r(r1)
                        matrix `bs'_ci_bca[2, `j'] = r(r2)
                    }
                }
            }
        }
    }
end

* Rep_vals: the row of one bootstrap replicate or jackknife deletion, as parenthesized
*   expressions in c_local vals: for each outcome (res(): its Engine matrices; m1(), m2():
*   its group means) focal_group raked_focal ref_group difference explained unexplained, then with
*   stages() the marginal explained parts, cumulative explained and unexplained parts and raked
*   means by stage, then with rres() the same in reverse order; last, each outcome's explained share
program define Rep_vals
    syntax, res(string) m1(string) m2(string) nsb(integer) [ stages(integer 0) rres(string) ///
        contrib dends(string) ]
    // the columns of e(b)'s stage equations: marginal explained, cumulative explained and
    // unexplained, raked mean; contrib (individually): the contributions and the remainder
    // (column 5), then the raked means of the blocks' own rakings (column 1, not the last)
    local cols 5 2 3 1
    if "`contrib'" != "" local cols 5 1
    local vals
    local j 0
    foreach R of local res {
        local ++j
        local a : word `j' of `m1'
        local c : word `j' of `m2'
        local vals `vals' (`a') (`R'[`nsb', 1]) (`c') (`a' - `c') (`R'[`nsb', 2]) (`R'[`nsb', 3])
        local mats `R'
        if "`rres'" != "" local mats `mats' `: word `j' of `rres''
        if `stages' > 0 {
            foreach M of local mats {
                // detail (dends()): R holds the term chain; a block's values are those of the term
                // stage at its last term, and its increment the difference of the explained parts
                // at the last terms of it and of the block before (as Engine forms increments);
                // then each term's increment
                if "`dends'" != "" & "`M'" == "`R'" {
                    foreach col of local cols {
                        local e0
                        foreach e of local dends {
                            if `col' != 5 local vals `vals' (`M'[`e', `col'])
                            else if "`e0'" == "" local vals `vals' (`M'[`e', 2])
                            else local vals `vals' (`M'[`e', 2] - `M'[`e0', 2])
                            local e0 `e'
                        }
                    }
                    forvalues s = 1/`nsb' {
                        local vals `vals' (`M'[`s', 5])
                    }
                    continue
                }
                foreach col of local cols {
                    forvalues s = 1/`stages' {
                        if "`contrib'" != "" & `col' == 1 & `s' == `stages' continue
                        local vals `vals' (`M'[`s', `col'])
                    }
                }
            }
        }
    }
    foreach R of local res {
        local vals `vals' (`R'[`nsb', 4])
    }
    c_local vals `vals'
end

* ============================================================================
* Chain_rep: chain, one bootstrap replicate or BCa deletion: Engine for every
*   pair on the weights w() (the pair's sample, Focal Group flag and coded terms from ptouse()
*   pg1() ptv(), its level counts from pnl(), pair-major), Bootstrap's rule for a failed
*   raking applied to each, and the replicate's row in matrix row(): for each pair focal_group
*   raked_focal ref_group difference explained unexplained (with stages() also the marginal
*   explained, cumulative explained and cumulative unexplained parts and the raked means by
*   stage), then the total's difference, explained and unexplained parts (and the summed
*   marginal explained parts), then the
*   explained share of each pair and of the total. only(): the pairs to compute (a BCa
*   deletion touches the pairs that contain the unit), the others copied from obs(), the
*   observed row. c_local ch_rc (0; 1, a Break; else the failure's code), ch_why (it, tg,
*   tb, rk: as Bootstrap), ch_fpair (the pair that failed), ch_itl (1 when a raking stopped at
*   the iterate() limit with its margins met) and ch_itst (the stages that stopped there)
* ============================================================================
program define Chain_rep
    syntax, w(varname) npairs(integer) ptouse(varlist) pg1(varlist) ptv(varlist) pnl(string) ///
        y(varname) one(varname) nterms(integer) nstages(integer) nsb(integer) sdb(string)     ///
        res(name) row(name) [ ipfopts(string asis) engopts(string asis)                       ///
        tnames(string) holdothers stages(integer 0) only(numlist) obs(name) force ]
    local kp = 6 + 4 * `stages'
    local kt = 3 + `stages'
    local toff = `npairs' * `kp'
    local poff = `toff' + `kt'
    matrix `row' = J(1, `poff' + `npairs' + 1, .)
    forvalues c = 1/`kt' {
        matrix `row'[1, `toff' + `c'] = 0
    }
    local rc 0
    local why rk
    local fpair
    local itl 0
    local itst
    local frc 0
    forvalues j = 1/`npairs' {
        local off = (`j' - 1) * `kp'
        if "`only'" != "" & !`: list j in only' {
            // a deletion that does not touch the pair: its point estimates
            forvalues c = 1/`kp' {
                matrix `row'[1, `off' + `c'] = `obs'[1, `off' + `c']
            }
            matrix `row'[1, `poff' + `j'] = `obs'[1, `poff' + `j']
        }
        else {
            local pt : word `j' of `ptouse'
            local pg : word `j' of `pg1'
            local tvl
            local nll
            forvalues t = 1/`nterms' {
                local tvl `tvl' `: word `=(`j' - 1) * `nterms' + `t'' of `ptv''
                local nll `nll' `: word `=(`j' - 1) * `nterms' + `t'' of `pnl''
            }
            Engine, touse(`pt') g1(`pg') w(`w') y(`y') one(`one') nterms(`nterms') tvlist(`tvl') ///
                nlevlist(`nll') nstages(`nsb') stagedefs(`sdb') res(`res') ipfopts(`ipfopts')    ///
                tnames(`tnames') `engopts' `holdothers'
            local rc = r(rc)
            if `rc' == 1 continue, break
            if `rc' & r(trimclash) == 1 local why tb
            if `rc' == 0 {
                forvalues s = 1/`nsb' {
                    local st = cond(`nsb' == 1, `nstages', `s')
                    // a missed margin fails the pair, unless force keeps it (frc); itb: at the
                    // iterate() limit, and not for the trimming bounds
                    if `res'[`s', 7] > 0 & "`force'" != "" local frc 1
                    else if `res'[`s', 7] > 0 {
                        local itb = (`res'[`s', 6] == 0 & r(trimchg`s') != 1)
                        if `rc' == 0 local why = cond(`itb', "it", "tg")
                        if `itb' & !`: list st in itst' local itst `itst' `st'
                        local rc 430
                    }
                    else if `res'[`s', 6] == 0 & `res'[`s', 7] == 0 {
                        local itl 1
                        if !`: list st in itst' local itst `itst' `st'
                    }
                }
            }
            if `rc' {
                local fpair `j'
                continue, break
            }
            matrix `row'[1, `off' + 1] = r(m1_1)
            matrix `row'[1, `off' + 2] = `res'[`nsb', 1]
            matrix `row'[1, `off' + 3] = r(m2_1)
            matrix `row'[1, `off' + 4] = r(m1_1) - r(m2_1)
            matrix `row'[1, `off' + 5] = `res'[`nsb', 2]
            matrix `row'[1, `off' + 6] = `res'[`nsb', 3]
            local c 6
            foreach col in 5 2 3 1 {
                forvalues s = 1/`stages' {
                    matrix `row'[1, `off' + `c' + `s'] = `res'[`s', `col']
                }
                local c = `c' + `stages'
            }
            matrix `row'[1, `poff' + `j'] = `res'[`nsb', 4]
        }
        // the total: the sums over the pairs, in the order of the pairs (as Chain sums them)
        forvalues c = 1/3 {
            matrix `row'[1, `toff' + `c'] = `row'[1, `toff' + `c'] + `row'[1, `off' + 3 + `c']
        }
        forvalues s = 1/`stages' {
            matrix `row'[1, `toff' + 3 + `s'] = `row'[1, `toff' + 3 + `s'] + `row'[1, `off' + 6 + `s']
        }
    }
    if `rc' == 0 {
        matrix `row'[1, `poff' + `npairs' + 1] = cond(`row'[1, `toff' + 1] != 0, ///
            100 * `row'[1, `toff' + 2] / `row'[1, `toff' + 1], .)
    }
    c_local ch_rc `rc'
    c_local ch_why `why'
    c_local ch_fpair `fpair'
    c_local ch_itl `itl'
    c_local ch_itst `itst'
    c_local ch_forced `frc'
end

* ============================================================================
* Recut: bayesboot() with recut, before a replicate or a BCa deletion: the bins of the c.
*   variables vars() (q() groups each, into the variables bins()) from the weights w() by the
*   rule of contvars() (weighted quantiles on sample(), over the observations with positive
*   weight), then the levels of each term that uses them (terms(): its number, its variable
*   and its components, | between terms), numbered as seldecomp_group numbers them. The level
*   counts go to c_local s_nll (nlevlist() with those terms' counts replaced) and s_rnll (the
*   reverse order: rtix() holds each of its terms' positions in the forward order); s_rcrc is
*   0, or _pctile's return code when a deletion leaves no observation to take quantiles of
* ============================================================================
program define Recut
    syntax, w(varname) touse(varname) g1(varname) sample(string) vars(string) q(string) ///
        bins(string) terms(string) nlevlist(string) [ rtix(string) ]
    local cond `touse' & `w' > 0
    if "`sample'" == "focalgroup" local cond `cond' & `g1' == 1
    else if "`sample'" == "refgroup" local cond `cond' & `g1' == 0
    local i 0
    foreach x of local vars {
        local ++i
        local qi : word `i' of `q'
        local b : word `i' of `bins'
        capture _pctile `x' [aw=`w'] if `cond', nq(`qi')
        if _rc {
            c_local s_rcrc = _rc
            exit
        }
        forvalues j = 1/`=`qi'-1' {
            tempname c`j'
            scalar `c`j'' = r(r`j')
        }
        qui replace `b' = `qi' if `touse'
        forvalues j = `=`qi'-1'(-1)1 {
            qui replace `b' = `j' if `touse' & `x' <= scalar(`c`j'')
        }
    }
    tempname LV
    local nll `nlevlist'
    local rest `terms'
    while `"`rest'"' != "" {
        gettoken piece rest : rest, parse("|")
        if `"`piece'"' == "|" continue
        gettoken t piece : piece
        gettoken tv comps : piece
        mata: seldecomp_group(st_local("comps"), "`touse'", "`tv'", "`LV'")
        local nl = rowsof(`LV')
        local new
        local j 0
        foreach n of local nll {
            local ++j
            local new `new' `=cond(`j' == `t', `nl', `n')'
        }
        local nll `new'
    }
    local rnll
    foreach j of local rtix {
        local rnll `rnll' `: word `j' of `nll''
    }
    c_local s_nll `nll'
    c_local s_rnll `rnll'
    c_local s_rcrc 0
end


* ============================================================================
* Display: the output in sections, each after a rule (Display_sep): the header, with the
*   definitions and the notes, and the warning on the replicates a prefix dropped; the raking
*   diagnostics; the covariate detail (covariatedetail); the decomposition, for each outcome
*   its stage table(s) and then its main table, the last table of the output (with citable
*   its bootstrap interval table, the same rows with intervals, follows as an appendix).
*   Chain: the same sections, the pairs' stage tables before Display_chain
* ============================================================================
program define Display, rclass
    syntax [, Level(cilevel) DIAGnostics COVariatedetail cformat(string) ///
              pformat(string) sformat(string) CITable VARLabels replay ]
    // varlabels: the variable labels in place of the variable names, in every routine below
    mata: seldecomp_vlset(`=("`varlabels'" != "")')
    local fopts
    foreach f in cformat pformat sformat {
        if `"``f''"' != "" local fopts `fopts' `f'(``f'')
    }
    Display_fmts, `fopts'
    local fmts cf(`cf') pf(`pf') sf(`sf')
    // the display of the estimation: seldecomp's own, or the replay that bootstrap: and
    // jackknife: run after their loop, whose c(prefix) holds the prefix (a replay typed by the
    // user has none); failed replicates get one message there, a one-line footnote at replay
    local atest = ("`replay'" == "" | strtrim(c(prefix)) != "")
    local ownbs = ("`e(bsnames)'" != "")           // seldecomp's own bayesboot() results
    // its intervals at the displayed level: e() at the estimation level, else recomputed
    tempname ci
    local ciopts
    if `ownbs' {
        Display_ci, level(`level') prefix(`ci')
        local ciopts ci(`ci') cilevel(`cilevel')
    }
    // the interval types shown (e(citypes)); citable adds percentile and bc, the default of
    // normal and t (t with bayesboot(cluster())) unless one of them is shown, and bca when it
    // was computed. Two types or more: the coefficient table shows the first, the bootstrap
    // interval table all of them
    local citypes
    if `ownbs' {
        local citypes `e(citypes)'
        if "`citypes'" == "" local citypes normal
        if "`citable'" != "" {
            local add percentile bc
            local tnorm normal t
            local has : list citypes & tnorm
            if "`has'" == "" local add `add' `=cond("`e(clustvar)'" != "", "t", "normal")'
            if e(bca) == 1 local add `add' bca
            local citypes : list citypes | add
        }
        local ciopts `ciopts' citypes(`citypes')
    }
    // (a replay shows the interval table only when citable is typed, as at estimation it shows
    // what bayesboot() requested; its coefficient table keeps the first type)
    local citab = (`: word count `citypes'' > 1)
    if "`replay'" != "" & "`citable'" == "" local citab 0
    // chain: its own table layout, after the pairs' stage tables
    local chain = ("`e(chain)'" != "")
    // 1. the header, the definitions and the notes included; 2. the replicates that bootstrap:
    // or jackknife: dropped: one warning, in the display that follows their loop (Misreps_warn;
    // at replay Display_table's footnote), none when the note of the iterate() limit in the
    // header counts them (Misreps_fold)
    Display_header
    if `atest' & !`ownbs' Misreps_warn, shownotes
    if `ownbs' Display_jkfail
    // 3. detailed diagnostics appear on request or when a raking check fails.
    Display_diag_scan
    local showdiag = ("`diagnostics'" != "" | `diag_issue')
    if `showdiag' {
        Display_sep
        Display_diag
    }
    // 4. the covariate detail, computed at estimation with covariatedetail (never under chain)
    if "`covariatedetail'" != "" {
        Display_sep
        capture confirm matrix e(covariates)
        if _rc di as txt _n "(covariatedetail was not specified at estimation; nothing to display)"
        else Display_covdetail
    }
    // 5. the decomposition: titled "Decomposition", or by the stage table that comes first
    Display_sep
    local stg = (`e(N_stages)' > 1) | (!`chain' & "`e(detail)'" != "" & e(N_terms) > e(N_blocks))
    if !`stg' di as txt _n "Decomposition"
    if `chain' {
        // the pairs' stage tables (a staged run of two or more stages), then the main table
        // and, with citable, its bootstrap interval table
        if `e(N_stages)' > 1 {
            forvalues p = 1/`e(N_pairs)' {
                Display_stages, `fmts' pair(`p')
            }
        }
        Display_chain, level(`level') `fmts' `ciopts' atest(`atest')
        if `citab' {
            Display_chain, level(`level') `fmts' `ciopts' bootci
        }
    }
    else {
        // outcomes(): for each outcome its stage table(s) (detail: also the one stage of a
        // single block of two or more terms, its term rows; reverse: the reverse-order table
        // too), then its main table and, with citable, its bootstrap interval table
        local nout = max(1, wordcount("`e(outcomes)'"))
        forvalues o = 1/`nout' {
            if `e(N_stages)' > 1 | ("`e(detail)'" != "" & e(N_terms) > e(N_blocks)) {
                Display_stages, `fmts' outcome(`o')
            }
            if `e(N_stages)' > 1 & "`e(reverse)'" != "" Display_stages, `fmts' outcome(`o') reverse
            Display_table, level(`level') `fmts' `ciopts' outcome(`o') atest(`atest')
            if `citab' {
                Display_bootci, level(`level') `fmts' `ciopts' outcome(`o')
            }
        }
    }
    // citable with nothing to show says so
    if "`citable'" != "" & !`ownbs' {
        capture confirm matrix e(V)
        if _rc di as txt _n "(citable: no confidence intervals without bayesboot() or a prefix)"
        else {
            Display_vce
            di as txt _n "(citable: no bayesboot() replicates; the CIs above are `vcelab')"
        }
    }
    Display_raking_footer, shown(`showdiag') issue(`diag_issue')
    // at replay the intervals displayed are left in r(); e() keeps those of the estimation
    if `ownbs' {
        foreach m in ci_normal ci_t ci_percentile ci_bc ci_bca pct_explained_ci {
            capture confirm matrix `ci'_`m'
            if _rc continue
            if "`replay'" != "" return matrix `m' = `ci'_`m'
            else matrix drop `ci'_`m'
        }
        if "`replay'" != "" return scalar level = `cilevel'
    }
end

program define Display_jkfail
    // The BCa calculation retains completed deletions; make an incomplete jackknife explicit.
    if e(N_jkfail) > 0 & e(N_jkfail) < . {
        di as txt "{p 0 4 2 78}Caution: `=e(N_jkfail)' jackknife replicates failed; BCa intervals may be incorrect.{p_end}"
    }
end

program define Display_sep
    // the rule between the sections of the output: a blank line and a 78-column rule of
    // "="; the blank line under the rule is the _n of each section's first line
    di ""
    di as txt _dup(78) "="
end

program define Display_ci
    // seldecomp's own bootstrap intervals at level(), in `prefix'_ci_normal, _ci_t,
    // _ci_percentile, _ci_bc, _ci_bca (bca) and _pct_explained_ci (the explained share): those of e() at the
    // estimation level; at another level recomputed from the replicates in e(bsdraws) by
    // Bootstrap_ci, as the estimation at that level computes them. cilevel: their level.
    // outcomes(): _pct_explained_ci has one row per outcome (e(by_outcome), e(bsdraws_pct))
    syntax, level(cilevel) prefix(name)
    local bca = (e(bca) == 1)
    local nout = max(1, wordcount("`e(outcomes)'"))
    // chain: a share for each pair and the total (e(by_pair), e(pct_explained_ci))
    local ch = ("`e(chain)'" != "")
    if `ch' local nout = e(N_pairs) + 1
    capture confirm matrix e(bsdraws)
    if `level' == e(level) | _rc {
        foreach m in ci_normal ci_t ci_percentile ci_bc {
            matrix `prefix'_`m' = e(`m')
        }
        if `bca' matrix `prefix'_ci_bca = e(ci_bca)
        if `ch' {
            tempname BP
            matrix `BP' = e(by_pair)
            local cl = colnumb(`BP', "pct_explained_ll")
            matrix `prefix'_pct_explained_ci = `BP'[1..., `cl'..`cl' + 1] \ e(pct_explained_ci)
            local rn : rownames `BP'
            matrix rownames `prefix'_pct_explained_ci = `rn' total
            matrix colnames `prefix'_pct_explained_ci = ll ul
        }
        else if `nout' > 1 {
            matrix `prefix'_pct_explained_ci = e(by_outcome)
            matrix `prefix'_pct_explained_ci = `prefix'_pct_explained_ci[1..., 9..10]
            matrix colnames `prefix'_pct_explained_ci = ll ul
        }
        else {
            capture confirm matrix e(pct_explained_ci)
            if !_rc matrix `prefix'_pct_explained_ci = e(pct_explained_ci)
        }
        local lv = string(e(level), "%12.0g")
        c_local cilevel `lv'
        exit
    }
    tempname b D DP OB A PO
    matrix `b' = e(b)
    matrix `D' = e(bsdraws)
    local bsnames `e(bsnames)'
    local kb : word count `bsnames'
    matrix `OB' = J(1, `kb', .)
    forvalues j = 1/`kb' {
        local c = colnumb(`b', "`: word `j' of `bsnames''")
        matrix `OB'[1, `j'] = `b'[1, `c']
    }
    local opts df(`=e(df_bs)')
    if e(ties) == 1 local opts `opts' ties
    if `bca' {
        matrix `A' = e(accel)
        local opts `opts' accel(`A')
    }
    capture confirm matrix e(bsdraws_pct)
    local haspct = (_rc == 0)
    if `haspct' matrix `DP' = e(bsdraws_pct)
    local nreps = e(N_reps)
    preserve
    qui drop _all
    qui svmat double `D', names(`prefix'_v)
    local vl
    forvalues j = 1/`kb' {
        local vl `vl' `prefix'_v`j'
    }
    if `haspct' {
        qui svmat double `DP', names(`prefix'_p)
        local pl
        forvalues k = 1/`=colsof(`DP')' {
            local pl `pl' `prefix'_p`k'
        }
        local opts `opts' pct(`pl')
        if `nout' > 1 {
            if `ch' {
                matrix `PO' = e(by_pair)
                matrix `PO' = (`PO'[1..., 7] \ e(pct_explained))'
            }
            else {
                matrix `PO' = e(by_outcome)
                matrix `PO' = (`PO'[1..., 7])'
            }
            local opts `opts' pobs(`PO')
            forvalues k = 1/`nout' {
                local pct_ll_`k' .
                local pct_ul_`k' .
            }
        }
    }
    Bootstrap_ci `vl', obs(`OB') level(`level') nreps(`nreps') prefix(`prefix') `opts'
    restore
    matrix drop `prefix'_b_bs `prefix'_bias `prefix'_se `prefix'_z0
    local ms ci_normal ci_t ci_percentile ci_bc
    if `bca' local ms `ms' ci_bca
    foreach m of local ms {
        matrix colnames `prefix'_`m' = `bsnames'
        matrix rownames `prefix'_`m' = ll ul
    }
    if `haspct' {
        if `nout' > 1 {
            matrix `prefix'_pct_explained_ci = J(`nout', 2, .)
            forvalues k = 1/`nout' {
                matrix `prefix'_pct_explained_ci[`k', 1] = `pct_ll_`k''
                matrix `prefix'_pct_explained_ci[`k', 2] = `pct_ul_`k''
            }
            local rn `e(outcomes)'
            if `ch' {
                local rn
                forvalues k = 1/`=`nout' - 1' {
                    local rn `rn' p`k'
                }
                local rn `rn' total
            }
            matrix rownames `prefix'_pct_explained_ci = `rn'
        }
        else matrix `prefix'_pct_explained_ci = (`pct_ll', `pct_ul')
        matrix colnames `prefix'_pct_explained_ci = ll ul
    }
    c_local cilevel `level'
end

program define Display_pct
    // Display_pct lmacname value: the value with one decimal (%9.1f), trimmed, never "-0.0"
    args lname x
    local s = strtrim(string(`x', "%9.1f"))
    if real("`s'") == 0 local s 0.0
    c_local `lname' `s'
end

program define Display_pct6
    // Display_pct6 lmacname value: a percentage of the stage table in at most 6 characters,
    // so that it never runs into the next column: one decimal below 1,000, none below
    // 100,000, else a short power of ten (-1.4e5, 2e12); missing is "."
    args lname x
    if missing(`x') {
        c_local `lname' .
        exit
    }
    local s = string(`x', "%6.1f")
    if abs(`x') >= 1000 | strlen("`s'") > 6 local s = string(`x', "%6.0f")
    if abs(`x') >= 1e5 | strlen("`s'") > 6 {
        local sg = cond(`x' < 0, "-", "")
        local k = floor(log10(abs(`x')))
        local m = round(abs(`x') / 10^`k', .1)
        if `m' >= 10 {
            local m = `m' / 10
            local ++k
        }
        local s = "`sg'" + string(`m', "%3.1f") + "e" + string(`k')
        if strlen("`s'") > 6 {
            local m = round(abs(`x') / 10^`k')
            if `m' >= 10 {
                local m 1
                local ++k
            }
            local s = "`sg'" + string(`m') + "e" + string(`k')
        }
    }
    c_local `lname' `s'
end

program define Display_eps
    // Display_eps lmacname [#]: the rounding threshold of outcome #, 1e-13 times the
    // largest magnitude of its focal_group, raked_focal and ref_group in e(b). Every value of the
    // decomposition is a mean or a difference of such means, computed in double precision
    // (about 16 significant digits); a value smaller in magnitude than the threshold is rounding
    // error of that subtraction and is displayed as 0. The display only: e() keeps the values
    args lname o
    if "`o'" == "" local o 1
    local eq0 = cond("`e(outcomes)'" != "", word("`e(outcomes)'", `o'), "overall")
    // chain: # is the pair (its means), 0 the total (the means of every pair)
    if "`e(chain)'" != "" {
        local eq0 p`o'
        if `o' == 0 {
            local eq0
            forvalues j = 1/`=e(N_pairs)' {
                local eq0 `eq0' p`j'
            }
        }
    }
    tempname b
    matrix `b' = e(b)
    local s 0
    foreach eq of local eq0 {
        foreach r in focal_group raked_focal ref_group {
            local c = colnumb(`b', "`eq':`r'")
            if `c' < . local s = max(`s', abs(`b'[1, `c']))
        }
    }
    c_local `lname' = 1e-13 * `s'
end

program define Display_fmts
    // the formats of the tables: the option, else -set cformat- (pformat, sformat), else
    // Stata's defaults; returned in cf, pf and sf. check: validate only, without the note
    syntax [, cformat(string) pformat(string) sformat(string) check ]
    foreach f in c p s {
        local maxw = cond("`f'" == "c", 9, cond("`f'" == "p", 5, 8))
        local fmt `"``f'format'"'
        if `"`fmt'"' != "" {
            confirm numeric format `fmt'
            if fmtwidth(`"`fmt'"') > `maxw' {
                if "`check'" == "" di as txt "note: invalid `f'format(), using default"
                local fmt
            }
        }
        if `"`fmt'"' == "" local fmt `c(`f'format)'
        if `"`fmt'"' == "" local fmt = cond("`f'" == "c", "%9.0g", cond("`f'" == "p", "%5.3f", "%8.2f"))
        c_local `f'f `fmt'
    }
end

program define Display_vce
    // one source for the inference label of the header (vcelab) and the column heading
    // (vcehead), at estimation and at replay; vcebs = 1 for the bootstrap layout
    local lab `"`e(vcetype)'"'
    local head `"`e(vcetype)'"'
    local bs = inlist("`e(vce)'", "bootstrap", "bayesboot")
    if "`e(bsnames)'" != "" {
        local lab "Bayesian bootstrap"
        local head Bootstrap
        local bs 1
    }
    else if "`e(prefix)'" == "bayesboot" | "`e(vce)'" == "bayesboot" {
        // Stata 19's bayesboot prefix sets e(vcetype) only after it has displayed the results
        local lab "Bayesian bootstrap"
        local head "Bayesian bootstrap"
    }
    c_local vcelab `"`lab'"'
    c_local vcehead `"`head'"'
    c_local vcebs `bs'
end

program define Prefix_epilog, eclass
    // run by a replication prefix once its results are posted (e(prefix_epilog)). After svy,
    // whose table is its own: the message on failed replicates (Misreps_warn), and the design
    // effects as missing, which estat effects would otherwise show as 0: they need a
    // linearized (SRS) variance, which seldecomp does not compute.
    // bootstrap: and jackknife: display through seldecomp's replay, which gives the message
    if "`e(prefix)'" != "svy" exit
    tempname D
    matrix `D' = J(1, colsof(e(b)), .)
    matrix colnames `D' = `: colfullnames e(b)'
    local dl deff deft
    if "`e(subpop)'" != "" local dl `dl' deffsub deftsub
    foreach d of local dl {
        ereturn matrix `d' = `D', copy
    }
    Misreps_fold
    Misreps_warn
    // svy shows its own table, with e(b)'s names: before it, what that table does not say:
    // the two groups (or the chain's pairs), the mode lines, the order of the
    // difference and the raking status. seldecomp typed afterwards redisplays everything
    mata: seldecomp_vlset(ustrregexm(" " + st_global("e(prefix_dispopts)") + " ", " varl[a-z]* "))
    Display_header brief
    Display_table_note
    Display_raking_status
    di as txt "(type {stata seldecomp} to redisplay these results in seldecomp's own tables)"
end

program define Misreps_fold, eclass
    // replicates of a prefix that failed at the iterate() limit where the point estimate
    // stopped at that limit too (one message per cause): its
    // note counts them ("raking stopped at the iterate() limit, also in 5 of 8 BRR replicates
    // (short of the targets, dropped); raise iterate()"), when they stopped at its
    // stages and every replicate beyond those the pre-flight check foresaw failed so. That
    // note is then the message where the notes are shown (Misreps_warn) and at replay
    // (Display_table). Run after the loop: by Prefix_epilog, and by the replay that
    // bootstrap: and jackknife: run. Not under chain, whose pairs' stages are not one note's
    if !(e(N_misreps) > 0 & e(N_misreps) < .) | "`e(bsnames)'" != "" | "`e(chain)'" != "" exit
    local nf = e(N_misreps)
    local nr = e(N_reps) + e(N_misreps)
    local vce `e(vce)'
    local what = cond("`vce'" == "brr", "BRR", cond("`vce'" == "sdr", "SDR", "`vce'"))
    local pre = cond(e(N_prefail) < ., e(N_prefail), 0)
    mata: seldecomp_failcache("get", "", "")
    local nit 0
    local itst
    forvalues i = 1/`fc_n' {
        if "`fc_kind`i''" == "iterate" {
            local nit = `nit' + `fc_count`i''
            local itst `itst' `fc_text`i''
        }
    }
    if `nit' == 0 | `nit' != `nf' - `pre' exit
    // the point estimate's stages that stopped at the limit with every margin met (the note
    // "raking stopped at the iterate() limit before the weights converged")
    capture confirm matrix e(bystage)
    if _rc exit
    tempname BS
    matrix `BS' = e(bystage)
    local pest
    forvalues s = 1/`=rowsof(`BS')' {
        if `BS'[`s', 6] == 0 & `BS'[`s', 7] == 0 local pest `pest' `s'
    }
    if !`: list itst in pest' exit
    local suf " before the weights converged; raise iterate()"
    local nn = cond(e(N_notes) < ., e(N_notes), 0)
    forvalues m = 1/`nn' {
        mata: st_local("nt", st_global("e(note`m')"))
        if strpos(`"`macval(nt)'"', "raking stopped ") & !strpos(`"`macval(nt)'"', ", also in ") & ///
            substr(`"`macval(nt)'"', -strlen("`suf'"), .) == "`suf'" {
            local nt = substr(`"`macval(nt)'"', 1, strlen(`"`macval(nt)'"') - strlen("`suf'"))
            // (after svy without "dropped", as in Misreps_warn)
            local dr = cond("`e(prefix)'" == "svy", "", ", dropped")
            ereturn local note`m' `"`macval(nt)', also in `nit' of `nr' `what' replicates (short of the targets`dr'); raise iterate()"'
            continue, break
        }
    }
end

program define Fold_has
    // Fold_has <lmacname> "<text>": 1 when a note of the results (e(note#)) contains text
    // (the replicates counted in the note of the iterate() limit: Misreps_fold, Estimate)
    args lname txt
    local f 0
    local nn = cond(e(N_notes) < ., e(N_notes), 0)
    forvalues m = 1/`nn' {
        mata: st_local("nt", st_global("e(note`m')"))
        if strpos(`"`macval(nt)'"', `"`txt'"') local f 1
    }
    c_local `lname' `f'
end

program define Misreps_warn
    // replicates that failed under a replication prefix: after
    // the prefix's loop (Prefix_epilog after svy, Display in the replay that bootstrap: and
    // jackknife: run), one warning that names what failed them, from the record of the
    // replicate calls (seldecomp_failcache), and what to do. When the full-sample call's
    // pre-flight check (svy with replicate weights, or the PSU jackknife) foresaw failures
    // (e(N_prefail) > 0), its warning was the message, and this one is given only for the
    // replicates it did not foresee. With shownotes (Display, whose header shows the notes;
    // Prefix_epilog after svy, whose table shows none, gives none), none when the
    // note of the iterate() limit counts them (Misreps_fold). The levels are named (at most 3).
    // After svy the message says that the replicates failed and no more: svy's own
    // notes give the consequence for the standard errors (computed from the other replicates,
    // or missing when the PSU jackknife is left with one replicate in a stratum). bootstrap:
    // and jackknife: drop a failed replicate, and the message says so. At replay, Display_table
    // gives a one-line footnote, by the same rule
    syntax [, shownotes ]
    if !(e(N_misreps) > 0 & e(N_misreps) < .) | "`e(bsnames)'" != "" exit
    local vce `e(vce)'
    local nf = e(N_misreps)
    local nr = e(N_reps) + e(N_misreps)
    local what = cond("`vce'" == "brr", "BRR", cond("`vce'" == "sdr", "SDR", "`vce'"))
    if "`shownotes'" != "" {
        Fold_has fin ", also in `nf' of `nr' `what' replicates ("
        if `fin' exit
    }
    // the causes recorded: the levels without focal group weight (at most 3 named), the raking
    // stopped at the iterate() limit short of its targets, missed targets, by() with one
    // value, chain's values of by() missing (at most 3 named), another error (the first)
    mata: seldecomp_failcache("get", "", "")
    local nlev 0
    local nsup 0
    local tgt 0
    local itl 0
    local grp
    local oth
    local nch 0
    forvalues i = 1/`fc_n' {
        if "`fc_kind`i''" == "support" {
            local nsup = `nsup' + `fc_count`i''
            local ++nlev
            if `nlev' <= 3 {
                local lev`nlev' `"`macval(fc_text`i')'"'
            }
        }
        else if "`fc_kind`i''" == "targets" local tgt 1
        else if "`fc_kind`i''" == "iterate" local itl = `itl' + `fc_count`i''
        else if "`fc_kind`i''" == "group" local grp `fc_text`i''
        else if "`fc_kind`i''" == "chain" {
            // chain: a value of by() missing from the replicate
            local ++nch
            if `nch' <= 3 {
                local chv`nch' `"`macval(fc_text`i')'"'
            }
        }
        else if `"`macval(oth)'"' == "" {
            local oth `"`macval(fc_text`i')'"'
            local oth : subinstr local oth "{" "`=char(1)'", all
            local oth : subinstr local oth "}" "{c )-}", all
            local oth : subinstr local oth "`=char(1)'" "{c -(}", all
        }
    }
    local pre = e(N_prefail)
    if `pre' < . & `pre' > 0 {
        if `nf' <= `pre' exit
        local nf = `nf' - `pre'
        if `nsup' <= `pre' local nlev 0
        local lead = "`nf' further `what' " + plural(`nf', "replicate")
    }
    else local lead "`nf' of `nr' `what' replicates"
    local lead "`lead' failed"
    if "`e(prefix)'" != "svy" local lead = "`lead' and " + cond(`nf' == 1, "was", "were") + " dropped"
    local why
    local sep
    if `nlev' > 0 {
        local why `"`macval(lev1)'"'
        if `nlev' == 2 {
            local why `"`macval(why)' or `macval(lev2)'"'
        }
        else if `nlev' >= 3 {
            local why `"`macval(why)', `macval(lev2)'"'
            if `nlev' == 3 {
                local why `"`macval(why)' or `macval(lev3)'"'
            }
            else {
                local why `"`macval(why)', `macval(lev3)' or `=`nlev' - 3' more levels"'
            }
        }
        local why `"`macval(why)' had no focal group weight"'
        local sep ", or "
    }
    if `itl' & !`tgt' {
        local why `"`macval(why)'`sep'raking stopped at the iterate() limit short of its targets"'
        local sep ", or "
    }
    else if `tgt' {
        local why `"`macval(why)'`sep'raking did not reach its targets"'
        if `itl' local why `"`macval(why)' (in some, at the iterate() limit)"'
        local sep ", or "
    }
    if "`grp'" != "" {
        local why `"`macval(why)'`sep'by(`grp') took only one value"'
        local sep ", or "
    }
    if `nch' > 0 {
        local chw `"`macval(chv1)'"'
        if `nch' == 2 {
            local chw `"`macval(chw)' or `macval(chv2)'"'
        }
        else if `nch' == 3 {
            local chw `"`macval(chw)', `macval(chv2)' or `macval(chv3)'"'
        }
        else if `nch' > 3 {
            local chw `"`macval(chw)', `macval(chv2)', `macval(chv3)' or `=`nch' - 3' more values"'
        }
        local why `"`macval(why)'`sep'`macval(chw)' had no weight"'
        local sep ", or "
    }
    if `"`macval(oth)'"' != "" {
        local why `"`macval(why)'`sep'`macval(oth)'"'
    }
    if `"`macval(why)'"' != "" {
        local why `": in `=cond(`nf' == 1, "it", "them")', `macval(why)'"'
    }
    // what to do, when it is certain: collapse a level without Focal Group weight; resample
    // bootstrap: within the groups when a resample lacks one
    local fix
    if `nlev' > 0 {
        local fix = " Collapse " + cond(`nlev' == 1, "that level", "those levels") + ///
            cond("`vce'" == "brr", ", or use Fay's BRR (fay())", "") + "."
    }
    if "`grp'" != "" & "`e(prefix)'" == "bootstrap" local fix "`fix' Resample within the groups: bootstrap, strata(`grp')."
    else if `nch' > 0 & "`e(prefix)'" == "bootstrap" local fix "`fix' Use bootstrap, strata(`e(by)')."
    di as txt `"{p 0 4 2 78}warning: `lead'`macval(why)'.`fix'{p_end}"'
end

program define Display_header
    // the title and counts, the groups, weights, inference, the mode lines,
    // and the notes of the results.
    // two columns: a left entry must end by col 47 (one blank before _col(49)); one that does
    // not fit leaves its slot blank and is printed in full below the rows, as a {p} paragraph
    // when longer than 78 columns (p0, p1), so that nothing is truncated; by() as shown
    // (varlabels: its label)
    // brief (Prefix_epilog, after svy): without the title and counts, the weights and the
    // inference lines, which svy's own header gives
    args brief
    local brief = ("`brief'" == "brief")
    Vlab by "`e(by)'"
    // chain: the pairs in a table below, in place of the group lines
    local ch = ("`e(chain)'" != "")
    // group labels are user text: copy without macro expansion, measure, escape SMCL braces;
    // the two group lines stay together
    local gfit 1
    foreach g in 1 2 {
        mata: st_local("g`g'", st_global("e(" + (`g' == 1 ? "focal" : "ref") + "_group_label)"))
        local gw`g' = 24 + udstrlen(`"`by'"') + udstrlen(`"`macval(g`g')'"')
        if `gw`g'' > 47 local gfit 0
        local g`g' : subinstr local g`g' "{" "`=char(1)'", all
        local g`g' : subinstr local g`g' "}" "{c )-}", all
        local g`g' : subinstr local g`g' "`=char(1)'" "{c -(}", all
    }
    if `brief' di ""
    else di as txt _n "`e(title)'" _col(49) "Number of obs" _col(67) "= " as res %10.0fc e(N)
    // svy with subpop(): the header N is svy's, all of its estimation sample; the
    // subpopulation's count below it, as svy shows them
    if e(N_sub) < . & !`brief' {
        di as txt _col(49) "Subpop. no. obs" _col(67) "= " as res %10.0fc e(N_sub)
    }
    if `ch' {
        // a long by() name: the line follows the Pairs row, in full
        if udstrlen(`"Chain of adjacent pairs of `by'"') <= 47 {
            di as txt `"Chain of adjacent pairs of `by'"' _col(49) "Pairs" _col(67) "= " as res %10.0fc e(N_pairs)
        }
        else {
            di as txt _col(49) "Pairs" _col(67) "= " as res %10.0fc e(N_pairs)
            di as txt `"Chain of adjacent pairs of `by'"'
        }
    }
    else if `gfit' {
        di as txt `"Focal Group (raked): `by' = "' as res `"`macval(g1)'"' ///
            as txt _col(49) "Focal obs" _col(67) "= " as res %10.0fc e(N_focal)
        di as txt `"Reference Group:     `by' = "' as res `"`macval(g2)'"' ///
            as txt _col(49) "Reference obs" _col(67) "= " as res %10.0fc e(N_ref)
    }
    else {
        di as txt _col(49) "Focal obs" _col(67) "= " as res %10.0fc e(N_focal)
        di as txt _col(49) "Reference obs" _col(67) "= " as res %10.0fc e(N_ref)
    }
    local wfit 1
    local wtp = cond(`brief', "", "`e(wtype)'")
    if "`wtp'" != "" & `ch' {
        // chain: the weights of each pair are in e(by_pair)
        local wtxt `"Weights: `e(wtype)' `e(wexp)'"'
        local ww = udstrlen(`"`wtxt'"')
        local wfit 0
    }
    else if "`wtp'" != "" {
        local wtxt `"Weights: `e(wtype)' `e(wexp)'"'
        local ww = udstrlen(`"`wtxt'"')
        local wfit = (`ww' <= 47)
        local wcol 69
        local wf %10.0gc
        if max(e(W_focal), e(W_ref)) >= 1e7 {
            local wcol 68
            local wf %11.0fc
        }
        if `wfit' {
            di as txt `"`wtxt'"' _col(49) "Focal wgt sum" _col(67) "=" as res _col(`wcol') `wf' e(W_focal)
        }
        else di as txt _col(49) "Focal wgt sum" _col(67) "=" as res _col(`wcol') `wf' e(W_focal)
        di as txt _col(49) "Reference wgt sum" _col(67) "=" as res _col(`wcol') `wf' e(W_ref)
    }
    if !`gfit' & !`ch' {
        // hanging indent 21: continuation lines start under the by() variable
        foreach g in 1 2 {
            local gt = cond(`g' == 1, "Focal Group (raked): ", "Reference Group:{space 5}")
            local p0 = cond(`gw`g'' > 78, "{p 0 21 2 78}", "")
            local p1 = cond(`gw`g'' > 78, "{p_end}", "")
            di as txt `"`p0'`gt'`by' = "' as res `"`macval(g`g')'"' as txt "`p1'"
        }
    }
    if !`wfit' {
        local p0 = cond(`ww' > 78, "{p 0 9 2 78}", "")
        local p1 = cond(`ww' > 78, "{p_end}", "")
        di as txt `"`p0'`wtxt'`p1'"'
    }
    Display_vce
    if `"`vcelab'"' != "" & !`brief' {
        local itxt "Inference: `vcelab'"
        local iw = udstrlen("`itxt'")
        local p0 = cond(`iw' > 78, "{p 0 11 2 78}", "")
        local p1 = cond(`iw' > 78, "{p_end}", "")
        if e(N_reps) >= . di as txt "`p0'`itxt'`p1'"
        else if `iw' <= 47 di as txt "`itxt'" _col(49) "Replications" _col(67) "= " as res %10.0fc e(N_reps)
        else di as txt _col(49) "Replications" _col(67) "= " as res %10.0fc e(N_reps)
        if "`e(bsnames)'" != "" & e(N_clust) < . {
            Vlab cvl "`e(clustvar)'"
            di as txt "{ralign 78:(Replications based on {res:`=strtrim(string(e(N_clust), "%12.0fc"))'} clusters in `cvl')}"
        }
        if e(N_reps) < . & `iw' > 47 di as txt "`p0'`itxt'`p1'"
    }
    if "`e(outcomes)'" != "" {
        Vlab otl "`e(outcomes)'"
        di as txt "{p 0 4 2 78}Outcomes (one raking): {res:`otl'}{p_end}"
    }
    if "`e(ratio)'" != "" {
        Vlab num "`e(depvar)'"
        if "`e(outcomes)'" != "" local num each outcome
        Vlab rvl "`e(ratio)'"
        di as txt "{p 0 4 2 78}Ratio of means: {res:`num'} / {res:`rvl'}{p_end}"
    }
    // the mode lines: the mode word, and whether the other blocks are held (fixothermargins)
    local how = cond("`e(fixothermargins)'" != "", "others held at focal group margins", "others unconstrained")
    if "`e(stagetype)'" == "individually" {
        // each block raked alone from the focal group's starting point, then all blocks
        local nb = e(N_blocks)
        di as txt "{p 0 4 2 78}Blocks raked " as res "individually" as txt " (`how'):{p_end}"
        forvalues s = 1/`nb' {
            local bn : word `s' of `e(blocks)'
            local bd `"`e(block`s')'"'
            Block_lab bnl `bn'
            Vlab bdl "`bd'"
            if "`bn'" == "`bd'" | `"`bnl'"' == `"`bdl'"' di as txt "{p 2 14 2 78}block " as res %2.0f `s' as txt ": " as res `"`bnl'"' as txt "{p_end}"
            else di as txt "{p 2 14 2 78}block " as res %2.0f `s' as txt ": " as res `"`bnl'"' as txt `" (`bdl'){p_end}"'
        }
    }
    else if "`e(stagetype)'" != "" {
        local ns = e(N_stages)
        di as txt "{p 0 4 2 78}Blocks added " as res "sequentially" as txt " (`how'):{p_end}"
        forvalues s = 1/`ns' {
            local bn : word `s' of `e(blocks)'
            local bd `"`e(block`s')'"'
            Block_lab bnl `bn'
            Vlab bdl "`bd'"
            if "`bn'" == "`bd'" | `"`bnl'"' == `"`bdl'"' di as txt "{p 2 14 2 78}stage " as res %2.0f `s' as txt ": + " as res `"`bnl'"' as txt "{p_end}"
            else di as txt "{p 2 14 2 78}stage " as res %2.0f `s' as txt ": + " as res `"`bnl'"' as txt `" (`bdl'){p_end}"'
        }
        if "`e(reverse)'" != "" {
            local rb
            forvalues s = `ns'(-1)1 {
                Block_lab bnl `: word `s' of `e(blocks)''
                local rb `"`rb'`=cond(`s' < `ns', " ", "")'+ {res:`bnl'}"'
            }
            di as txt "{p 2 4 2 78}and in reverse order: `rb'{p_end}"
        }
        if "`e(detail)'" != "" di as txt "{p 2 4 2 78}and term by term within each block (detail){p_end}"
    }
    else {
        Vlab rkl "`e(rakevars)'"
        di as txt "{p 0 18 2 78}Raking variables: {res:`rkl'}{p_end}"
    }
    if "`e(insupport)'" != "" & `ch' {
        // chain: the Focal Group observations excluded over the pairs (e(by_pair)); the notes list them
        tempname BP
        matrix `BP' = e(by_pair)
        local nx 0
        forvalues j = 1/`=rowsof(`BP')' {
            local nx = `nx' + `BP'[`j', colnumb(`BP', "N_focal_all")] - `BP'[`j', colnumb(`BP', "N_focal")]
        }
        if `nx' > 0 {
            di as txt "{p 0 4 2 78}Support (insupport): {res:`=strtrim(string(`nx', "%12.0fc"))'} " ///
                "focal group obs excluded over the pairs (levels absent from the pair's reference group){p_end}"
        }
        else {
            di as txt "{p 0 4 2 78}Support (insupport): no observation excluded{p_end}"
        }
    }
    else if "`e(insupport)'" != "" {
        // insupport: how many Focal Group observations the screen excluded, and the mean before it
        local nx = e(N_focal_all) - e(N_focal)
        if `nx' > 0 {
            local nxs = strtrim(string(`nx', "%12.0fc"))
            local nas = strtrim(string(e(N_focal_all), "%12.0fc"))
            local mas = strtrim(string(e(mean_focal_all), "%9.0g"))
            Vlab dvl "`=word("`e(depvar)'", 1)'"
            local mtxt = cond("`e(ratio)'" != "", "ratio", "mean") + " of `dvl'"
            di as txt "{p 0 4 2 78}Support (insupport): {res:`nxs'} of {res:`nas'} focal group obs " ///
                "excluded (levels absent from the reference group); focal group `mtxt' before the " ///
                "exclusion: {res:`mas'}{p_end}"
        }
        else {
            di as txt "{p 0 4 2 78}Support (insupport): no observation excluded{p_end}"
        }
    }
    if "`e(cutvars)'" != "" {
        // binned c. variables: the number of bins (the nonempty ones when ties leave fewer,
        // e(cutbins)) and the sample of the cut points
        local ctxt
        local c 0
        foreach cv in `e(cutvars)' {
            local ++c
            local q : word `c' of `e(cutq)'
            local nb : word `c' of `e(cutbins)'
            local qtxt "`q' quantile bins"
            if "`nb'" != "" {
                if `nb' < `q' local qtxt "`nb' nonempty of `q' quantile bins (ties)"
            }
            Vlab cvl "`cv'"
            local ctxt = "`ctxt'" + cond(`c' > 1, ", ", "") + "{res:`cvl'}: `qtxt'"
        }
        local sep = cond(`c' > 1, ";", ",")
        local csx = cond(`ch', "each pair's `e(contsample)' sample", "the `e(contsample)' sample")
        di as txt "{p 2 3 2 78}(`ctxt'`sep' cut points from `csx'){p_end}"
    }
    // chain: the pairs, each with its groups and their observations; then the total
    if `ch' {
        tempname BP PR
        matrix `BP' = e(by_pair)
        matrix `PR' = e(pairs)
        local J = rowsof(`PR')
        mata: _sd_lab = strtrim(ustrsplit(st_global("e(by_labels)"), "[|]"))
        // "<by> = " before each value when every cell fits 26 characters; else the values alone
        // (the header names by())
        mata: st_local("bypfx", udstrlen(st_local("by")) + 3 + max(udstrlen(_sd_lab)) <= 26 ? st_local("by") + " = " : "")
        di ""
        di as txt "Pair  Focal Group (raked)" _col(34) "Reference Group" _col(61) "Focal obs" _col(71) "Ref. obs"
        forvalues j = 1/`J' {
            // the texts of the two values (user text: cut to 26 characters, SMCL braces escaped)
            foreach q in 1 2 {
                mata: st_local("t`q'", st_local("bypfx") + _sd_lab[`j' + `q' - 1])
                mata: st_local("t`q'", udstrlen(st_local("t`q'")) > 26 ? usubstr(st_local("t`q'"), 1, 23) + "..." : st_local("t`q'"))
                local t`q' : subinstr local t`q' "{" "`=char(1)'", all
                local t`q' : subinstr local t`q' "}" "{c )-}", all
                local t`q' : subinstr local t`q' "`=char(1)'" "{c -(}", all
            }
            di as txt %4s "p`j'" "  " as res `"`macval(t1)'"' _col(34) `"`macval(t2)'"' ///
                _col(61) %9.0fc `BP'[`j', colnumb(`BP', "N_focal")] _col(70) %9.0fc `BP'[`j', colnumb(`BP', "N_ref")]
        }
        foreach q in 1 2 {
            mata: st_local("t`q'", st_local("by") + " = " + _sd_lab[`q' == 1 ? 1 : cols(_sd_lab)])
            local t`q' : subinstr local t`q' "{" "`=char(1)'", all
            local t`q' : subinstr local t`q' "}" "{c )-}", all
            local t`q' : subinstr local t`q' "`=char(1)'" "{c -(}", all
        }
        mata: mata drop _sd_lab
        di as txt `"{p 0 6 2 78}total: the sum of the pairs, from {res:`macval(t1)'} to {res:`macval(t2)'}{p_end}"'
    }
    // Thin-cell and target-failure details belong in the raking diagnostics block.
    // Keep the stored notes intact for callers and replicate-warning accounting.
    local thin 0
    local missnote 0
    local shown 0
    // (brief, before svy's table: no notes; the failed replicates have their one message,
    // Misreps_warn, and seldecomp typed afterwards shows the notes)
    if `brief' exit
    forvalues m = 1/`e(N_notes)' {
        mata: st_local("nt", seldecomp_vlnote(st_global("e(note`m')")))
        // continuation lines (levels listed under a note) start with four blanks
        local cont = (substr(`"`macval(nt)'"', 1, 4) == "    ")
        if !`cont' {
            Thin_is thin `"`macval(nt)'"'
            mata: st_local("missnote", strofreal(regexm(st_local("nt"), "(^|: )(trimming left [0-9,]+ raking margins? unmet|[0-9,]+ raking margins? not met|raking stopped .* with [0-9,]+ raking margins? not met)")))
        }
        if `thin' | `missnote' continue
        mata: st_local("nt", seldecomp_smcl(st_local("nt")))
        if !`shown' di ""
        local shown 1
        if `cont' {
            di as txt `"{p 6 8 2 78}`macval(nt)'{p_end}"'
        }
        else {
            di as txt `"{p 0 6 2 78}note: `macval(nt)'{p_end}"'
        }
    }
end

* Thin_is <lmacname> "<note>": c_local lmacname 1 if the note is the head of a thin-cell note
*   ("thin target level, stages 1-3:", "3 thin target levels:", "p1, p3: thin target level:"), else 0
program define Thin_is
    args lname nt
    mata: st_local("t", strofreal(ustrregexm(st_local("nt"), "^(p[0-9]+(, p[0-9]+)*: )?([0-9,]+ )?thin target levels?[,:]")))
    c_local `lname' `t'
end


program define Display_table
    syntax [, Level(cilevel) cf(string) pf(string) sf(string) ci(name) cilevel(string) ///
              outcome(integer 1) atest(integer 1) citypes(string) ]
    capture confirm matrix e(V)
    local hasV = (_rc == 0)
    Display_vce
    tempname b V P PSE MZ ES ECS
    matrix `b' = e(b)
    local k = colsof(`b')
    local names : colfullnames `b'
    if `hasV' matrix `V' = e(V)
    Vlab depvar "`e(depvar)'" 16
    // the columns of this outcome: all of e(b), or with outcomes() the outcome's block, whose
    // equations are shown without the outcome prefix (overall, marg_expl, cum_expl, ...)
    local nout = max(1, wordcount("`e(outcomes)'"))
    local j0 = (`outcome' - 1) * `k' / `nout' + 1
    local j1 = `outcome' * `k' / `nout'
    scalar `P' = e(pct_explained)
    scalar `PSE' = e(pct_explained_se)
    scalar `MZ' = e(mean_focal_zero)
    scalar `ES' = e(explained_support)
    scalar `ECS' = e(explained_cs)
    local ylab `e(depvar)'
    local eq0 overall
    if `nout' > 1 {
        local yo : word `outcome' of `e(outcomes)'
        local ylab `yo'
        local eq0 `yo'
        Vlab depvar "`yo'" 16
        local dfn = colnumb(`b', "`yo':difference")
        scalar `P' = cond(`b'[1, `dfn'] != 0, 100 * `b'[1, `dfn' + 1] / `b'[1, `dfn'], .)
        scalar `PSE' = .
        capture confirm matrix e(by_outcome)
        if !_rc {
            scalar `P' = el(e(by_outcome), `outcome', 7)
            scalar `PSE' = el(e(by_outcome), `outcome', 8)
        }
        if `outcome' > 1 & e(N_focal_zero) > 0 {
            // change of support: the outcome's own Focal Group means (e(bystage#), last stage)
            scalar `MZ' = el(e(bystage`outcome'), e(N_stages), 29)
            local mcs = el(e(bystage`outcome'), e(N_stages), 30)
            scalar `ES' = `b'[1, `dfn' - 3] - `mcs'
            scalar `ECS' = `mcs' - `b'[1, `dfn' - 2]
        }
    }
    if "`e(ratio)'" != "" local ylab `ylab'/`e(ratio)'
    // values below the rounding threshold of the outcome are displayed as 0 (Display_eps)
    Display_eps eps `outcome'
    // t with the residual (design) df when the prefix set one (svy, jackknife), else z
    local df = e(df_r)
    local st = cond(`df' < ., "t", "z")
    if `df' < . local crit = invttail(`df', (100 - `level') / 200)
    else local crit = invnormal(1 - (100 - `level') / 200)
    // seldecomp's own bootstrap: the interval column holds the first of the interval types
    // shown (citypes), normal and t from the std. err. (t with e(df_bs) degrees of freedom), the
    // others the intervals of Display_ci (`ci'_ci_<type>, at cilevel)
    local ownbs = ("`e(bsnames)'" != "")
    local mt normal
    local ccrit `crit'
    if `ownbs' {
        local mt : word 1 of `citypes'
        if "`mt'" == "normal" local ccrit = invnormal(1 - (100 - `level') / 200)
        else if "`mt'" == "t" local ccrit = invttail(e(df_bs), (100 - `level') / 200)
        else {
            tempname CM
            matrix `CM' = `ci'_ci_`mt'
        }
    }
    // A 16-character row label; numeric columns end at 29, 40, 49, 55, 66 and 78.
    local cw = fmtwidth("`cf'")
    local pw = fmtwidth("`pf'")
    local sw = fmtwidth("`sf'")
    di ""
    Table_headings
    di as txt "{hline 17}{c TT}{hline 60}"
    if `hasV' {
        // the variance type is never abbreviated: all but its last word go on a line above
        local nw : word count `vcehead'
        local vlast : word `nw' of `vcehead'
        local vfirst
        forvalues i = 1/`=`nw' - 1' {
            local vfirst `vfirst' `: word `i' of `vcehead''
        }
        local cilab "[`level'% `h_ci']"
        if !inlist("`mt'", "normal", "t") local cilab "[`cilevel'% `h_ci']"
        local cicol = 79 - udstrlen("`cilab'")
        // centred over the std. err. column (32-40)
        if "`vfirst'" != "" {
            di as txt "{space 17}{c |}" _col(`=32 + floor((9 - udstrlen("`vfirst'")) / 2)') "`vfirst'"
        }
        local vcol = 32 + floor((9 - udstrlen("`vlast'")) / 2)
        if `vcebs' {
            // the type of the intervals, centred over them
            Ci_label cihd `mt'
            di as txt "{space 17}{c |}" _col(`vcol') "`vlast'" _col(`=69 - floor(udstrlen("`cihd'") / 2)') "`cihd'"
            local hc `"`h_coefl'"'
        }
        else {
            di as txt "{space 17}{c |}" _col(`vcol') "`vlast'"
            local hc `"`h_coef'"'
        }
        di as txt `"{ralign 16:`macval(depvar)'} {c |}"' _col(19) "`hc'" _col(32) "`h_sel'" _col(49) "`st'" ///
            _col(51) "P>|`st'|" _col(`cicol') "`cilab'"
    }
    else {
        di as txt `"{ralign 16:`macval(depvar)'} {c |}"' _col(19) "`h_coef'"
    }
    di as txt "{hline 17}{c +}{hline 60}"
    local preveq
    forvalues j = `j0'/`j1' {
        local fn : word `j' of `names'
        gettoken eq nm : fn, parse(":")
        gettoken colon nm : nm, parse(":")
        if `nout' > 1 {
            if "`eq'" == "`yo'" local eq overall
            else local eq = substr("`eq'", strlen("`yo'") + 2, .)
        }
        // without inference the stage equations (the terms' and the reverse order's too) repeat
        // the stage tables above: the overall equation only; with it every equation, titled
        if !`hasV' & "`eq'" != "overall" continue
        if "`eq'" != "`preveq'" {
            if "`preveq'" != "" di as txt "{hline 17}{c +}{hline 60}"
            Eq_title ttl `eq'
            Eq_lab eql `eq'
            di as res %-17s "`eql'" as txt "{c |}`ttl'"
            local preveq `eq'
        }
        local coef = `b'[1, `j']
        if abs(`coef') < `eps' local coef 0
        // a stage row: its block or term as shown (varlabels: their labels)
        Row_lab nml `eq' `nm'
        di as txt `"{ralign 16:`macval(nml)'} {c |}"' as res _col(`=30 - `cw'') `cf' `coef' _c
        if `hasV' {
            // a std. err. below the threshold is 0: the row shows its coefficient only
            local se = sqrt(`V'[`j', `j'])
            if `se' < `eps' local se 0
            if `se' < . & `se' > 0 {
                local z = `coef' / `se'
                if `df' < . local p = 2 * ttail(`df', abs(`z'))
                else local p = 2 * normal(-abs(`z'))
                if inlist("`mt'", "normal", "t") {
                    local ll = `coef' - `ccrit' * `se'
                    local ul = `coef' + `ccrit' * `se'
                }
                else {
                    local cc = colnumb(`CM', "`fn'")
                    local ll = `CM'[1, `cc']
                    local ul = `CM'[2, `cc']
                }
                if abs(`ll') < `eps' local ll 0
                if abs(`ul') < `eps' local ul 0
                di as res _col(`=41 - `cw'') `cf' `se' _col(`=50 - `sw'') `sf' `z' ///
                    _col(`=56 - `pw'') `pf' `p' _col(`=67 - `cw'') `cf' `ll' _col(`=79 - `cw'') `cf' `ul'
            }
            else di ""
        }
        else di ""
    }
    di as txt "{hline 17}{c BT}{hline 60}"
    Display_table_note
    // t intervals or tests: their degrees of freedom
    if `ownbs' & `hasV' & ("`mt'" == "t" | `df' < .) {
        local dfw = cond("`e(clustvar)'" != "", "clusters", "observations")
        di as txt "{p 0 1 2 78}(t: `=e(df_bs)' degrees of freedom, the number of `dfw' less 1){p_end}"
    }
    // explained share, its bootstrap standard error and percentile interval when available;
    // outcomes(): this outcome's row. None when the difference is displayed as 0 (Display_eps).
    local cd = colnumb(`b', "`eq0':difference")
    if `cd' < . {
        if abs(`b'[1, `cd']) < `eps' scalar `P' = .
    }
    if `P' < . {
        Display_pct sh `P'
        local shtxt "{res:`sh'%}"
        if `hasV' & `PSE' < . {
            capture confirm matrix `ci'_pct_explained_ci
            local haspci = (_rc == 0)
            Display_pct sse `PSE'
            di as txt "{p 0 4 2 78}Explained share of the difference: `shtxt', bootstrap " ///
                "std. err. {res:`sse'}{p_end}"
            if `haspci' {
                Display_pct pl el(`ci'_pct_explained_ci,`outcome',1)
                Display_pct pu el(`ci'_pct_explained_ci,`outcome',2)
                di as txt "{p 4 4 2 78}`cilevel'% percentile CI: {res:`pl'%} to {res:`pu'%}{p_end}"
            }
        }
        else {
            di as txt "{p 0 4 2 78}Explained share of the difference: `shtxt'{p_end}"
        }
    }
    else di as txt "Explained share of the difference: n/a (the difference is 0)"
    if e(N_focal_zero) > 0 {
        local mtxt = cond("`e(ratio)'" != "", "ratio", "mean")
        foreach q in ES ECS MZ {
            if abs(``q'') < `eps' scalar ``q'' = 0
        }
        di as txt "{p 0 4 2 78}explained = " as res %9.0g `ES' as txt " (change of support: " ///
            as res "`e(N_focal_zero)'" as txt " focal group obs with " as res %4.1f e(pct_W_focal_zero) as txt ///
            "% of its weight, `mtxt' `ylab' " as res %9.0g `MZ' as txt ", get zero weight) + " ///
            as res %9.0g `ECS' as txt " (reweighting within the common support){p_end}"
    }
    if `hasV' & `outcome' == 1 {
        // failed replicates: at estimation one message,
        // given by Display (Misreps_warn) in the display that follows the loop of bootstrap: or
        // jackknife:, by Prefix_epilog after svy, by Estimate for seldecomp's own bootstrap. At
        // replay, a one-line footnote, for the bootstrap also when replications were kept at
        // the iterate() limit (e(N_iterlimit)); none when the note of that limit counts them
        // (Estimate, Misreps_fold)
        local nf = cond(e(N_misreps) < ., e(N_misreps), 0)
        local ni = cond(e(N_iterlimit) < ., e(N_iterlimit), 0)
        local nr = e(N_reps) + `nf'
        if !`atest' & "`e(bsnames)'" != "" & (`nf' > 0 | `ni' > 0) {
            Fold_has fin ", also in `=`nf' + `ni'' of `nr' replications"
            if `fin' {
            }
            else if `nf' > 0 {
                di as txt "{p 0 1 2 78}(`nf' of `nr' replications failed and " cond(`nf' == 1, "was", "were") ///
                    " dropped; bootstrap results may be incorrect){p_end}"
            }
            else {
                di as txt "{p 0 1 2 78}(`ni' of `nr' replications stopped before the weights converged){p_end}"
            }
        }
        else if !`atest' & "`e(bsnames)'" == "" & `nf' > 0 {
            local vce `e(vce)'
            local what = cond("`vce'" == "brr", "BRR", cond("`vce'" == "sdr", "SDR", "`vce'"))
            Fold_has fin ", also in `nf' of `nr' `what' replicates ("
            if !`fin' {
                // after svy "failed" alone (Misreps_warn)
                local dr = cond("`e(prefix)'" == "svy", "", " and " + cond(`nf' == 1, "was", "were") + " dropped")
                di as txt "{p 0 1 2 78}(`nf' of `nr' `what' replicates failed`dr'){p_end}"
            }
        }
        // force: the replications kept although they missed raking targets, at replay
        if !`atest' & e(N_forced) > 0 & e(N_forced) < . {
            di as txt "{p 0 1 2 78}(force kept `=e(N_forced)' of `nr' replications that missed raking targets){p_end}"
        }
    }
end

* Vlab lname "names" [width]: c_local lname, the variable names (or terms) as shown: with
*   varlabels their variable labels (seldecomp_vlab), braces escaped for display; with width,
*   a label shortened to width characters with "~", and a name abbreviated (abbrev())
program define Vlab
    args lname names w
    mata: st_local("t", seldecomp_vlab(st_local("names")))
    if "`w'" != "" {
        mata: st_local("same", strofreal(st_local("t") == st_local("names")))
        if `same' local t = abbrev(`"`names'"', `w')
        else mata: st_local("t", seldecomp_vltrim(st_local("t"), `w'))
    }
    mata: st_local("t", seldecomp_smcl(st_local("t")))
    c_local `lname' `"`macval(t)'"'
end

* Block_lab lname name [width]: c_local lname, block name as shown: with varlabels, a block that
*   is named after its terms (an unnamed block, whose name is its term's with # as _) shows
*   their variable labels; otherwise its name. Shortened and escaped as by Vlab
program define Block_lab
    args lname bn w
    local blks `e(blocks)'
    local b : list posof "`bn'" in blks
    local t `bn'
    mata: st_local("on", strofreal(seldecomp_vlon()))
    if `b' > 0 & `on' {
        local terms `"`e(block`b')'"'
        if "`bn'" == "`terms'" | "`bn'" == subinstr("`terms'", "#", "_", .) local t `terms'
    }
    Vlab out "`t'" `w'
    c_local `lname' `"`macval(out)'"'
end

* Row_lab lname eq name: c_local lname, the name of a row of the coefficient tables, in 16
*   columns: a stage row's block (Block_lab) or a detail row's term (Vlab) as shown, the overall
*   equation's rows with descriptive labels. Stored coefficient names do not change.
program define Row_lab
    args lname eq nm
    if inlist("`eq'", "overall", "") {
        local mtxt = cond("`e(ratio)'" != "", "ratio", "mean")
        local t = abbrev("`nm'", 16)
        if "`nm'" == "focal_group" local t "Focal `mtxt'"
        else if "`nm'" == "raked_focal" local t = cond("`e(ratio)'" != "", "Raked ratio", "Raked focal mean")
        else if "`nm'" == "ref_group" local t "Reference `mtxt'"
        else if "`nm'" == "difference" local t "Difference"
        else if "`nm'" == "explained" local t "Explained"
        else if "`nm'" == "unexplained" local t "Unexplained"
    }
    else if inlist("`eq'", "detail", "rev_detail") Vlab t "`nm'" 16
    else Block_lab t `nm' 16
    c_local `lname' `"`macval(t)'"'
end

* Eq_lab lname eq: the descriptive left heading of a coefficient equation.
program define Eq_lab
    args lname eq
    if substr("`eq'", 1, 4) == "rev_" local eq = substr("`eq'", 5, .)
    local t "By stage"
    if inlist("`eq'", "overall", "") local t "Overall"
    else if "`eq'" == "detail" local t "By term"
    else if "`e(stagetype)'" == "individually" local t "By block"
    c_local `lname' "`t'"
end

program define Display_table_note
    // The same identities hold for group means and ratios of weighted totals.
    di as txt "Difference = focal - reference. Explained = focal - raked focal."
end

* Ci_label lname type: c_local lname, the heading of the coefficient table's interval column
* for an interval type of bayesboot() (normal, t, percentile, bc, bca); Ci_code: its code in
* the bootstrap interval table, right-aligned in 6 columns
program define Ci_label
    args lname ty
    c_local `lname' = cond("`ty'" == "t", "t-based", cond("`ty'" == "percentile", "Percentile", ///
        cond("`ty'" == "bc", "Bias-corrected", cond("`ty'" == "bca", "BCa", "Normal-based"))))
end

program define Ci_code
    args lname ty
    c_local `lname' = cond("`ty'" == "t", "   (T)", cond("`ty'" == "percentile", "   (P)", ///
        cond("`ty'" == "bc", "  (BC)", cond("`ty'" == "bca", " (BCa)", "   (N)"))))
end

* Ci_key lname "types": c_local lname, the key of the codes of the interval table, the types
* shown (t with its degrees of freedom), ties adjusting the bias-corrected ones
program define Ci_key
    args lname types
    local key
    foreach ty of local types {
        if "`ty'" == "normal" local k "N: normal"
        else if "`ty'" == "t" local k "T: t with `=e(df_bs)' degrees of freedom"
        else if "`ty'" == "percentile" local k "P: percentile"
        else if "`ty'" == "bc" local k "BC: bias-corrected"
        else local k "BCa: bias-corrected and accelerated"
        local key = "`key'" + cond("`key'" != "", ", ", "") + "`k'"
    }
    if e(ties) == 1 & (`: list posof "bc" in types' | `: list posof "bca" in types') {
        local key "`key' (adjusted for ties)"
    }
    c_local `lname' "Key: `key'"
end

program define Table_headings
    // column headings follow the running Stata version's own coefficient tables (the estimate's
    // column titled Estimate, in the width of their Coefficient / Coef.); below a spanning
    // heading (the variance type) they are lowercase from Stata 17 on
    if c(stata_version) >= 17 {
        c_local h_coef  "   Estimate"
        c_local h_coefl "   Estimate"
        c_local h_se    "Std. err."
        c_local h_sel   "std. err."
        c_local h_ci    "conf. interval"
    }
    else {
        c_local h_coef  "  Estimate"
        c_local h_coefl "  Estimate"
        c_local h_se    "Std. Err."
        c_local h_sel   "Std. Err."
        c_local h_ci    "Conf. Interval"
    }
end

program define Display_stages
    // outcome(#): e(bystage#) for outcome # > 1 of outcomes(); reverse: the reverse order,
    // e(bystage_rev[#]), whose stage s adds block K + 1 - s
    syntax [, cf(string) pf(string) sf(string) outcome(integer 1) REVerse pair(integer 0) ]
    // chain, pair(#): the stage table of pair #, e(bystage[#]) and its equations
    // (p#:, p#_...), as an outcome's
    // individually: the table of the block contributions instead
    if "`e(stagetype)'" == "individually" {
        Display_contrib, cf(`cf') outcome(`outcome')
        exit
    }
    local cw = fmtwidth("`cf'")
    tempname S b
    local ns = e(N_stages)
    local mname = "bystage" + cond("`reverse'" != "", "_rev", "") + cond(`outcome' > 1, "`outcome'", "")
    local blocks `e(blocks)'
    if "`reverse'" != "" {
        local blocks
        forvalues s = `ns'(-1)1 {
            local blocks `blocks' `: word `s' of `e(blocks)''
        }
    }
    local ttl = "Decomposition by raking stage (sequentially" + cond("`e(fixothermargins)'" != "", ", fixothermargins", "") + ")"
    if "`reverse'" != "" local ttl "`ttl', reverse order"
    if "`e(outcomes)'" != "" {
        Vlab otl "`=word("`e(outcomes)'", `outcome')'" 16
        local ttl "`ttl': `otl'"
    }
    if `pair' > 0 {
        local mname = "bystage" + cond(`pair' > 1, "`pair'", "")
        Pair_title ttl `"`ttl'"' `pair'
    }
    // a title beyond 78 columns (fixothermargins, reverse order and an outcome) wraps, as
    // Pair_title wraps a pair's
    if `pair' == 0 & udstrlen("`ttl'") > 78 local ttl "{p 0 4 2 78}`ttl'{p_end}"
    capture confirm matrix e(`mname')
    if !_rc matrix `S' = e(`mname')
    else {
        // not posted (a replication prefix may drop it): the stage columns from e(b)
        local pfx
        if "`e(outcomes)'" != "" local pfx = word("`e(outcomes)'", `outcome') + "_"
        if `pair' > 0 local pfx p`pair'_
        if "`reverse'" != "" local pfx `pfx'rev_
        local eq0 = cond("`e(outcomes)'" != "", word("`e(outcomes)'", `outcome'), "overall")
        if `pair' > 0 local eq0 p`pair'
        matrix `b' = e(b)
        local g1 = `b'[1, colnumb(`b', "`eq0':focal_group")]
        local dd = `b'[1, colnumb(`b', "`eq0':difference")]
        matrix `S' = J(`ns', 27, .)
        forvalues s = 1/`ns' {
            local bn : word `s' of `blocks'
            matrix `S'[`s', 2] = `b'[1, colnumb(`b', "`pfx'cum_expl:`bn'")]
            matrix `S'[`s', 3] = `b'[1, colnumb(`b', "`pfx'cum_unexpl:`bn'")]
            matrix `S'[`s', 5] = `b'[1, colnumb(`b', "`pfx'marg_expl:`bn'")]
            matrix `S'[`s', 1] = `g1' - `S'[`s', 2]
            matrix `S'[`s', 4] = cond(`dd' != 0, 100 * `S'[`s', 2] / `dd', .)
            matrix `S'[`s', 27] = cond(`dd' != 0, 100 * `S'[`s', 5] / `dd', .)
        }
    }
    di as txt _n "`ttl'"
    di as txt "{hline 78}"
    // values end in columns 30, 41, 53, 60 (%), 71 and 78 (%)
    di as txt _col(33) "{hline 8} Cumulative {hline 8}" _col(62) "{hline 4} Marginal {hline 3}"
    local rlab = cond("`e(ratio)'" != "", "Raked ratio", "Raked mean")
    di as txt "Stage Block" _col(20) "`rlab'" _col(33) "Explained" _col(43) "Unexplained" ///
        _col(60) "%" _col(63) "Explained" _col(78) "%"
    di as txt "{hline 78}"
    // values below the rounding threshold are displayed as 0, and no percentages when the
    // difference itself is (Display_eps)
    Display_eps eps `=cond(`pair' > 0, `pair', `outcome')'
    local eq0 = cond("`e(outcomes)'" != "", word("`e(outcomes)'", `outcome'), "overall")
    if `pair' > 0 local eq0 p`pair'
    matrix `b' = e(b)
    local dz = abs(`b'[1, colnumb(`b', "`eq0':difference")]) < `eps'
    // detail (forward order): the terms of each block of two or more, indented under it, from
    // e(bystage_detail[#]) or, when it is not posted, from the detail elements of e(b); dcol#:
    // the column of term # in e(b), dnm#: its name there
    local det = ("`e(detail)'" != "" & "`reverse'" == "")
    local dmulti 0
    if `det' {
        tempname D
        local dmat = "bystage_detail" + cond(`outcome' > 1, "`outcome'", "")
        capture confirm matrix e(`dmat')
        local hasD = !_rc
        if `hasD' matrix `D' = e(`dmat')
        local deq = cond("`e(outcomes)'" != "", word("`e(outcomes)'", `outcome') + "_", "") + "detail"
        local bnames : colfullnames `b'
        local nd 0
        local j 0
        foreach fn of local bnames {
            local ++j
            gettoken eq nm : fn, parse(":")
            if "`eq'" != "`deq'" continue
            local ++nd
            local dcol`nd' `j'
            local dnm`nd' = substr("`nm'", 2, .)
        }
        local g1 = `b'[1, colnumb(`b', "`eq0':focal_group")]
        local dd = `b'[1, colnumb(`b', "`eq0':difference")]
    }
    local tt 0
    forvalues s = 1/`ns' {
        local nm : word `s' of `blocks'
        foreach c in 1 2 3 5 {
            local v`c' = `S'[`s', `c']
            if abs(`v`c'') < `eps' local v`c' 0
        }
        // percentages keep one decimal below 1,000 and stay within their 6 columns above
        // that (Display_pct6); one that rounds to zero is shown as 0.0, not -0.0
        local p4 = cond(`dz', ., `S'[`s', 4])
        local p27 = cond(`dz', ., `S'[`s', 27])
        if abs(`p4') < .05 local p4 0
        if abs(`p27') < .05 local p27 0
        Display_pct6 t4 `p4'
        Display_pct6 t27 `p27'
        Block_lab nml `nm' 12
        di as res %5.0f `s' as txt " " as res %-12s `"`nml'"' ///
            _col(`=31 - `cw'') `cf' `v1' _col(`=42 - `cw'') `cf' `v2' ///
            _col(`=54 - `cw'') `cf' `v3' _col(55) %6s "`t4'" ///
            _col(`=72 - `cw'') `cf' `v5' _col(73) %6s "`t27'"
        if !`det' continue
        // the block's terms: raked_focal at the term stage, the explained part accumulated within the
        // block (the term stage's minus the previous block stage's), the term's increment and its
        // percent of the difference
        local nbt : word count `e(block`s')'
        local cprev = cond(`s' > 1, `S'[`=`s' - 1', 2], 0)
        local cw0 0
        forvalues u = `=`tt' + 1'/`=`tt' + `nbt'' {
            if `nbt' < 2 continue
            local dmulti 1
            local inc = `b'[1, `dcol`u'']
            if `hasD' {
                local r1 = `D'[`u', 1]
                local cwb = `D'[`u', 2] - `cprev'
                local q = `D'[`u', 27]
            }
            else {
                local cwb = `cw0' + `inc'
                local r1 = `g1' - (`cprev' + `cwb')
                local q = cond(`dd' != 0, 100 * `inc' / `dd', .)
            }
            local cw0 `cwb'
            foreach x in r1 cwb inc {
                if abs(``x'') < `eps' local `x' 0
            }
            local q = cond(`dz', ., `q')
            if abs(`q') < .05 local q 0
            Display_pct6 tq `q'
            Vlab nml "`dnm`u''" 12
            di as res _col(9) %-12s `"`nml'"' ///
                _col(`=31 - `cw'') `cf' `r1' _col(`=42 - `cw'') `cf' `cwb' ///
                _col(`=72 - `cw'') `cf' `inc' _col(73) %6s "`tq'"
        }
        local tt = `tt' + `nbt'
    }
    di as txt "{hline 78}"
    // the footnote once: under the first outcome's forward table, and one line under the first
    // reverse table (chain: under the last pair's)
    if `outcome' > 1 exit
    if `pair' > 0 & `pair' < e(N_pairs) exit
    if "`reverse'" != "" {
        di as txt "{p 0 4 2 78}Reverse order: blocks added from last to first.{p_end}"
        exit
    }
    // (the variant, fixothermargins or not, is said in the header and in the table's title)
    di as txt "{p 0 4 2 78}%: percent of the difference.{p_end}"
    if `dmulti' {
        di as txt "{p 0 4 2 78}Indented: the block's terms, added in turn (Cumulative: within " ///
            "the block).{p_end}"
    }
end

program define Display_contrib
    // individually: for each block, its raking's raked_focal, explained part (the block's
    // contribution) and percent of the difference; then the sum of the blocks, the raking of all
    // blocks and the remainder. From e(bystage[#]), or rebuilt from e(b) when it is not posted
    syntax [, cf(string) outcome(integer 1) ]
    local cw = fmtwidth("`cf'")
    tempname S b
    local ns = e(N_stages)
    local nb = `ns' - 1
    local blocks `e(blocks)'
    local mname = "bystage" + cond(`outcome' > 1, "`outcome'", "")
    capture confirm matrix e(`mname')
    if !_rc matrix `S' = e(`mname')
    else {
        local pfx
        if "`e(outcomes)'" != "" local pfx = word("`e(outcomes)'", `outcome') + "_"
        local eq0 = cond("`e(outcomes)'" != "", word("`e(outcomes)'", `outcome'), "overall")
        matrix `b' = e(b)
        local g1 = `b'[1, colnumb(`b', "`eq0':focal_group")]
        local dd = `b'[1, colnumb(`b', "`eq0':difference")]
        matrix `S' = J(`ns', 27, .)
        forvalues s = 1/`ns' {
            local bn = cond(`s' < `ns', word("`blocks'", `s'), "_remainder")
            matrix `S'[`s', 5] = `b'[1, colnumb(`b', "`pfx'marg_expl:`bn'")]
            if `s' < `ns' matrix `S'[`s', 2] = `S'[`s', 5]
            else matrix `S'[`s', 2] = `b'[1, colnumb(`b', "`eq0':explained")]
            matrix `S'[`s', 1] = `g1' - `S'[`s', 2]
            matrix `S'[`s', 4] = cond(`dd' != 0, 100 * `S'[`s', 2] / `dd', .)
            matrix `S'[`s', 27] = cond(`dd' != 0, 100 * `S'[`s', 5] / `dd', .)
        }
    }
    local ttl = "Contribution of each block, raked alone (individually" + cond("`e(fixothermargins)'" != "", ", fixothermargins", "") + ")"
    if "`e(outcomes)'" != "" {
        Vlab otl "`=word("`e(outcomes)'", `outcome')'" 16
        local ttl "`ttl': `otl'"
    }
    // a title beyond 78 columns (fixothermargins and an outcome) wraps, as Pair_title wraps a pair's
    if udstrlen("`ttl'") > 78 local ttl "{p 0 4 2 78}`ttl'{p_end}"
    di as txt _n "`ttl'"
    di as txt "{hline 78}"
    // values end in columns 30 (raked_focal), 45 (explained) and 53 (%)
    local rlab = cond("`e(ratio)'" != "", "Raked ratio", "Raked mean")
    di as txt "Block" _col(20) "`rlab'" _col(37) "Explained" _col(53) "%"
    di as txt "{hline 78}"
    // values below the rounding threshold are displayed as 0, and no percentages when the
    // difference itself is (Display_eps)
    Display_eps eps `outcome'
    local eq0 = cond("`e(outcomes)'" != "", word("`e(outcomes)'", `outcome'), "overall")
    matrix `b' = e(b)
    local dz = abs(`b'[1, colnumb(`b', "`eq0':difference")]) < `eps'
    local sum 0
    local psum = cond(`dz', ., 0)
    forvalues s = 1/`ns' {
        // percentages as in the stage table: one decimal below 1,000, never -0.0
        local p = cond(`dz', ., cond(`s' < `ns', `S'[`s', 27], `S'[`s', 4]))
        if abs(`p') < .05 local p 0
        local fp`s' = cond(missing(`p') | abs(`p') < 1000, "%6.1f", cond(abs(`p') < 1e5, "%6.0f", "%6.0g"))
        local p`s' `p'
        if `s' < `ns' {
            local sum = `sum' + `S'[`s', 5]
            if !`dz' local psum = `psum' + `S'[`s', 27]
        }
        foreach c in 1 2 5 {
            local v`c'_`s' = `S'[`s', `c']
            if abs(`v`c'_`s'') < `eps' local v`c'_`s' 0
        }
    }
    if abs(`sum') < `eps' local sum 0
    forvalues s = 1/`nb' {
        local nm : word `s' of `blocks'
        Block_lab nml `nm' 21
        di as res %-21s `"`nml'"' _col(`=31 - `cw'') `cf' `v1_`s'' ///
            _col(`=46 - `cw'') `cf' `v5_`s'' _col(48) `fp`s'' `p`s''
    }
    di as txt "{hline 78}"
    local pr = cond(`dz', ., `S'[`ns', 27])
    foreach q in psum pr {
        if abs(``q'') < .05 local `q' 0
        local f`q' = cond(missing(``q'') | abs(``q'') < 1000, "%6.1f", cond(abs(``q'') < 1e5, "%6.0f", "%6.0g"))
    }
    di as txt "Sum of the blocks" as res _col(`=46 - `cw'') `cf' `sum' _col(48) `fpsum' `psum'
    di as txt "All blocks (overall)" as res _col(`=31 - `cw'') `cf' `v1_`ns'' ///
        _col(`=46 - `cw'') `cf' `v2_`ns'' _col(48) `fp`ns'' `p`ns''
    di as txt "Remainder" as res _col(`=46 - `cw'') `cf' `v5_`ns'' _col(48) `fpr' `pr'
    di as txt "{hline 78}"
    // the footnote once, under the first outcome's table
    if `outcome' > 1 exit
    local how = cond("`e(fixothermargins)'" != "", "are held at their focal group distributions", ///
        "are unconstrained")
    di as txt "{p 0 4 2 78}%: percent of the difference. Other blocks `how'. Remainder: " ///
        "explained part of all blocks less the sum of the blocks' parts.{p_end}"
end

* The diagnostic matrices supplement e(bystage) without changing its columns.
program define Display_diag_sets
    local sets
    capture confirm matrix e(bystage)
    if !_rc local sets base
    if "`e(chain)'" != "" {
        forvalues j = 2/`=e(N_pairs)' {
            capture confirm matrix e(bystage`j')
            if !_rc local sets `sets' `j'
        }
    }
    else {
        foreach q in rev detail {
            capture confirm matrix e(bystage_`q')
            if !_rc local sets `sets' _`q'
        }
    }
    c_local d_sets `sets'
end

program define Display_diag_get
    syntax, SET(string) S(name) C(name)
    local sx = cond("`set'" == "base", "", "`set'")
    matrix `s' = e(bystage`sx')
    capture matrix `c' = e(raking_checks`sx')
    local saved = (_rc == 0)
    if !`saved' {
        // Older saved estimates have no cell-level checks; do not invent their counts.
        matrix `c' = J(rowsof(`s'), 5, .)
        forvalues j = 1/`=rowsof(`s')' {
            matrix `c'[`j', 1] = `s'[`j', 6]
            if `s'[`j', 7] > 0 & `s'[`j', 6] == 0 matrix `c'[`j', 1] = .
            matrix `c'[`j', 2] = `s'[`j', 7]
        }
    }
    local ctx
    local caption
    local lp raking_miss`sx'
    if "`set'" == "_rev" {
        local ctx reverse order
        local caption Reverse order
    }
    else if "`set'" == "_detail" {
        local ctx detail
        local caption Term-by-term raking
    }
    else if "`e(chain)'" != "" {
        local j = cond("`set'" == "base", 1, real("`set'"))
        local ctx p`j'
        local caption Pair p`j'
        if `j' > 1 local lp raking_miss`j'_
    }
    else {
        capture confirm matrix e(bystage_rev)
        if !_rc {
            local ctx forward order
            local caption Forward order
        }
    }
    c_local d_sx `sx'
    c_local d_context `ctx'
    c_local d_caption `caption'
    c_local d_labelprefix `lp'
    c_local d_saved `saved'
end

program define Display_diag_scan
    Display_diag_sets
    local issue 0
    tempname S C
    foreach ds of local d_sets {
        Display_diag_get, set(`ds') s(`S') c(`C')
        forvalues j = 1/`=rowsof(`S')' {
            if `C'[`j', 1] != 1 | `C'[`j', 2] > 0 local issue 1
        }
    }
    c_local diag_issue `issue'
    c_local diag_available = ("`d_sets'" != "")
end

program define Display_raking_status
    Display_diag_sets
    if "`d_sets'" == "" exit
    tempname S C
    local baditer 0
    local badtarget 0
    local unknown 0
    local iw
    local tw
    local multi 0
    foreach ds of local d_sets {
        Display_diag_get, set(`ds') s(`S') c(`C')
        local ns = rowsof(`S')
        if `ns' > 1 | "`ds'" != "base" local multi 1
        local il
        local tl
        forvalues j = 1/`ns' {
            if `C'[`j', 1] < 1 {
                local baditer 1
                local il `il' `j'
            }
            if `C'[`j', 1] >= . local unknown 1
            if `C'[`j', 2] > 0 {
                local badtarget 1
                local tl `tl' `j'
            }
        }
        foreach k in i t {
            local loc
            if "``k'l'" != "" {
                local n : word count ``k'l'
                if `ns' > 1 {
                    local nums : subinstr local `k'l " " ", ", all
                    local nums = strtrim("`nums'")
                    if substr("`nums'", 1, 1) == "," local nums = strtrim(substr("`nums'", 2, .))
                    local loc = cond(`n' == 1, "stage ", "stages ") + "`nums'"
                }
                if "`d_context'" != "" local loc = "`d_context'" + cond("`loc'" != "", ": `loc'", "")
                if "`loc'" != "" local `k'w = "``k'w'" + cond("``k'w'" != "", "; ", "") + "`loc'"
            }
        }
    }
    local all = cond("`e(chain)'" != "", "all pairs and stages", "all stages")
    if !`baditer' & !`badtarget' & !`unknown' {
        local pfx = cond(`multi', "Raking (`all')", "Raking")
        di as txt "{p 0 4 2 78}`pfx': weight iterations converged; all targets met.{p_end}"
        exit
    }
    local it = cond(`baditer', "Convergence failure", cond(`unknown', "Iteration convergence unavailable", cond(`multi', "all weight iterations converged", "weight iterations converged")))
    if `baditer' & "`iw'" != "" local it "`it' (`iw')"
    local tt = cond(`badtarget', "targets *not* met", "all targets met")
    if `badtarget' & "`tw'" != "" local tt "`tt' (`tw')"
    di as txt "{p 0 4 2 78}Raking: `it'; `tt'.{p_end}"
end

program define Display_raking_footer
    syntax, SHOWN(integer) ISSUE(integer)
    Display_diag_sets
    if "`d_sets'" == "" exit
    if !`shown' | `issue' {
        di ""
        Display_raking_status
    }
    if `shown' exit
    tempname S C
    local lo .
    local hi 0
    local thin
    local nthin 0
    local multi 0
    foreach ds of local d_sets {
        Display_diag_get, set(`ds') s(`S') c(`C')
        local ns = rowsof(`S')
        if `ns' > 1 | "`ds'" != "base" local multi 1
        forvalues j = 1/`ns' {
            local lo = min(`lo', `S'[`j', 16])
            local hi = max(`hi', `S'[`j', 17])
        }
        if "`ds'" == "base" | "`e(chain)'" != "" {
            if `S'[`ns', 26] < . local nthin = `nthin' + `S'[`ns', 26]
            local n = strtrim(string(`S'[`ns', 26], "%12.0fc"))
            if "`e(chain)'" != "" {
                local thin = "`thin'" + cond("`thin'" != "", "; ", "") + "`d_context' = `n'"
            }
            else local thin `n'
        }
    }
    local lx = strtrim(string(`lo', "%9.4g"))
    local hx = strtrim(string(`hi', "%9.4g"))
    local all = cond("`e(chain)'" != "", "all pairs and stages", "all stages")
    local span = cond(`multi', " (`all')", "")
    di as txt "{p 0 4 2 78}Raking factors`span': `lx'-`hx'.{p_end}"
    // (no line when there is none)
    if `nthin' > 0 di as txt "{p 0 4 2 78}Thin target levels (<`e(mincell)' focal obs): `thin'.{p_end}"
end

program define Display_diag
    Display_diag_sets
    if "`d_sets'" == "" {
        di as txt _n "(raking diagnostics are unavailable for these saved results)"
        exit
    }
    Display_diag_scan
    local anyissue `diag_issue'
    tempname S C M
    di as txt _n "Raking diagnostics"
    foreach ds of local d_sets {
        Display_diag_get, set(`ds') s(`S') c(`C')
        local ns = rowsof(`S')
        local staged = (`ns' > 1)
        local dx = cond(`staged', 0, -5)
        local issue 0
        forvalues j = 1/`ns' {
            if `C'[`j', 1] != 1 | `C'[`j', 2] > 0 local issue 1
        }
        if "`d_caption'" != "" di as txt _n "`d_caption'"
        di as txt "{hline 78}"
        di as txt _col(`=9+`dx'') "-- Raking factor --" _col(`=35+`dx'') "Kish effective N" _col(`=55+`dx'') "ESS"
        local head = cond(`staged', "Stage", "")
        di as txt "`head'" _col(`=12+`dx'') "Min" _col(`=21+`dx'') "Max" _col(`=30+`dx'') "CV" ///
            _col(`=38+`dx'') "Base" _col(`=47+`dx'') "Raked" _col(`=55+`dx'') "ratio" ///
            _col(`=61+`dx'') "Zero wt" _col(`=70+`dx'') "At bound"
        di as txt "{hline 78}"
        forvalues j = 1/`ns' {
            local cv = `S'[`j', 19]
            if abs(`cv') < 1e-13 local cv 0
            local nf = cond(`S'[`j', 23] < 999999.5, "%10.1fc", "%10.0fc")
            local de = `S'[`j', 25]
            local df = cond(`de' < 100, "%7.3f", "%7.0g")
            local sl = cond(`staged', string(`j', "%5.0f"), "")
            di as res "`sl'" _col(`=6+`dx'') %9.4g `S'[`j', 16] %9.4g `S'[`j', 17] %8.4g `cv' ///
                `nf' `S'[`j', 23] `nf' `S'[`j', 24] `df' `de' %9.0fc `S'[`j', 21] %9.0fc `S'[`j', 22]
        }
        di as txt "{hline 78}"
        di as txt "{p 0 4 2 78}Raking factor = raked / base weight. ESS ratio = base / raked effective N.{p_end}"
        di as txt "{p 0 4 2 78}Zero wt: zero raked weights. At bound: distinct weights at a trimming bound.{p_end}"
        local nc = `S'[`ns', 26]
        di as txt "{p 0 4 2 78}Thin target levels (<`e(mincell)' focal obs): `=strtrim(string(`nc', "%12.0fc"))'.{p_end}"
        if !`issue' continue
        di as txt _n "Raking failure counts"
        di as txt "{hline 78}"
        di as txt "`head'" _col(`=8+`dx'') "Converged" _col(`=21+`dx'') "Targets met" ///
            _col(`=43+`dx'') "# Missed" _col(`=60+`dx'') "# Missed"
        di as txt _col(`=42+`dx'') "Variables" _col(`=60+`dx'') "Targets"
        di as txt "{hline 78}"
        forvalues j = 1/`ns' {
            local cv = cond(`C'[`j', 1] >= ., ".", cond(`C'[`j', 1] == 1, "yes", "no"))
            local tm = cond(`C'[`j', 2] == 0, "yes", "no")
            local sl = cond(`staged', string(`j', "%5.0f"), "")
            di as res "`sl'" _col(`=8+`dx'') "`cv'" _col(`=21+`dx'') "`tm'" ///
                _col(`=42+`dx'') %9.0f `C'[`j', 2] _col(`=58+`dx'') %9.0f `C'[`j', 3]
        }
        di as txt "{hline 78}"
        di as txt "{p 0 4 2 78}Converged: weight iterations. Targets met: within ctrltolerance().{p_end}"
        di as txt "{p 0 4 2 78}Missed variables: raking variables with one or more missed level targets.{p_end}"
        di as txt "{p 0 4 2 78}Missed targets: individual level totals outside ctrltolerance().{p_end}"
        if !`d_saved' di as txt "{p 0 4 2 78}Cell-level checks are unavailable for these saved results; re-estimate to display them.{p_end}"
        capture matrix `M' = e(raking_misses`d_sx')
        if _rc continue
        di as txt _n "Worst missed target for each affected variable"
        di as txt "{hline 78}"
        di as txt "`head'" _col(`=7+`dx'') "Variable" _col(`=20+`dx'') "Level" ///
            _col(`=44+`dx'') "Target" _col(`=54+`dx'') "Achieved" _col(`=66+`dx'') "Discrepancy"
        di as txt "{hline 78}"
        forvalues j = 1/`=rowsof(`M')' {
            mata: st_local("mv", st_global("e(`d_labelprefix'`j'_var)"))
            mata: st_local("ml", st_global("e(`d_labelprefix'`j'_level)"))
            Vlab dv "`mv'" 12
            mata: st_local("dl", seldecomp_smcl(seldecomp_vltrim(st_local("ml"), 18)))
            local sl = cond(`staged', string(`M'[`j', 1], "%5.0f"), "")
            di as res "`sl'" _col(`=7+`dx'') `"`macval(dv)'"' _col(`=20+`dx'') `"`macval(dl)'"' ///
                _col(`=39+`dx'') %10.4g `M'[`j', 4] _col(`=51+`dx'') %11.4g `M'[`j', 5] ///
                _col(`=65+`dx'') %11.4g `M'[`j', 6]
        }
        di as txt "{hline 78}"
        local ctl 1e-6
        if regexm("`e(rakeopts)'", "ctrltolerance\(([^)]+)\)") local ctl = real(regexs(1))
        di as txt "{p 0 4 2 78}Discrepancy = |s*achieved - s*target| / (|s*target| + 1), where s scales base weights to the raking sample size.{p_end}"
        di as txt "{p 0 4 2 78}Targets are met when discrepancy <= `=strtrim(string(`ctl', "%9.0g"))'.{p_end}"
    }
    if !`anyissue' {
        di ""
        Display_raking_status
    }
end

program define Eq_title
    // c_local mname: the title shown after the bar on the row of equation eq's name in the
    // coefficient tables, with a leading space (empty for the overall equation): what its
    // rows hold, one per stage or term
    args mname eq
    local rev = (substr("`eq'", 1, 4) == "rev_")
    if `rev' local eq = substr("`eq'", 5, .)
    local ind = ("`e(stagetype)'" == "individually")
    local mtxt = cond("`e(ratio)'" != "", "ratio", "mean")
    local t
    if "`eq'" == "marg_expl" {
        if `ind' local t "Contribution of each block raked alone, and the remainder"
        else local t "Marginal explained"
    }
    else if "`eq'" == "cum_expl" local t "Cumulative explained"
    else if "`eq'" == "cum_unexpl" local t "Cumulative unexplained"
    else if "`eq'" == "raked_focal" {
        if `ind' local t "Raked focal group `mtxt', each block raked alone"
        else local t "Raked focal group `mtxt'"
    }
    else if "`eq'" == "detail" local t "Marginal explained, by term"
    if `rev' & "`t'" != "" local t "`t', reverse order"
    if "`t'" != "" local t " `t'"
    c_local `mname' `"`t'"'
end

program define Display_bootci
    // the intervals (Display_ci) are at cilevel: the displayed level, or the estimation level
    // for results without e(bsdraws); outcome(#): the #th outcome of outcomes()
    syntax [, Level(cilevel) cf(string) pf(string) sf(string) ci(name) cilevel(string) ///
              outcome(integer 1) citypes(string) ]
    local cw = fmtwidth("`cf'")
    tempname b
    matrix `b' = e(b)
    local bsnames `e(bsnames)'
    local kb : word count `bsnames'
    Vlab depvar "`e(depvar)'" 16
    // outcomes(): the outcome's block of e(bsnames), equations shown without the prefix
    local nout = max(1, wordcount("`e(outcomes)'"))
    local j0 = (`outcome' - 1) * `kb' / `nout' + 1
    local j1 = `outcome' * `kb' / `nout'
    if `nout' > 1 {
        local yo : word `outcome' of `e(outcomes)'
        Vlab depvar "`yo'" 16
    }
    Display_eps eps `outcome'
    local lv = string(`cilevel', "%12.0g")
    di as txt _n "Bayesian bootstrap confidence intervals" _col(49) "Replications" _col(67) "= " ///
        as res %10.0fc e(N_reps)
    Table_headings
    di as txt "{hline 17}{c TT}{hline 60}"
    di as txt %16s "" " {c |}" _col(42) "Bootstrap"
    di as txt `"{ralign 16:`macval(depvar)'} {c |}"' _col(19) "`h_coefl'" _col(37) "Bias" _col(42) "`h_sel'" ///
        _col(54) "[`lv'% `h_ci']"
    di as txt "{hline 17}{c +}{hline 60}"
    // columns end at 29 (coefficient), 40 (bias), 50 (std. err.), 61 and 72 (interval)
    local c1 = 30 - `cw'
    local c2 = 41 - `cw'
    local c3 = 51 - `cw'
    local c4 = 62 - `cw'
    local c5 = 73 - `cw'
    local preveq
    forvalues j = `j0'/`j1' {
        local fn : word `j' of `bsnames'
        gettoken eq nm : fn, parse(":")
        gettoken colon nm : nm, parse(":")
        if `nout' > 1 {
            if "`eq'" == "`yo'" local eq overall
            else local eq = substr("`eq'", strlen("`yo'") + 2, .)
        }
        if "`eq'" != "`preveq'" {
            if "`preveq'" != "" di as txt "{hline 17}{c +}{hline 60}"
            Eq_title ttl `eq'
            Eq_lab eql `eq'
            di as res %-17s "`eql'" as txt "{c |}`ttl'"
            local preveq `eq'
        }
        local pj = colnumb(`b', "`fn'")
        // values below the rounding threshold are displayed as 0 (Display_eps)
        local x0 = `b'[1, `pj']
        local x1 = el(e(bias), 1, `j')
        local x2 = el(e(se), 1, `j')
        local k 2
        foreach M of local citypes {
            local x`++k' = el(`ci'_ci_`M', 1, `j')
            local x`++k' = el(`ci'_ci_`M', 2, `j')
        }
        forvalues i = 0/`k' {
            if abs(`x`i'') < `eps' local x`i' 0
        }
        // a row per interval type, the first beside the estimate, bias and std. err.
        local r 0
        foreach M of local citypes {
            local ++r
            local i1 = 1 + 2 * `r'
            local i2 = 2 + 2 * `r'
            Ci_code cc `M'
            if `r' == 1 {
                Row_lab nml `eq' `nm'
                di as txt `"{ralign 16:`macval(nml)'} {c |}"' as res _col(`c1') `cf' `x0' ///
                    _col(`c2') `cf' `x1' _col(`c3') `cf' `x2' ///
                    _col(`c4') `cf' `x`i1'' _col(`c5') `cf' `x`i2'' as txt "`cc'"
            }
            else di as txt %16s "" " {c |}" as res _col(`c4') `cf' `x`i1'' _col(`c5') `cf' `x`i2'' as txt "`cc'"
        }
    }
    di as txt "{hline 17}{c BT}{hline 60}"
    Ci_key key `"`citypes'"'
    di as txt "{p 0 4 2 78}`key'{p_end}"
    // outcomes(): said once, under the first outcome's table
    if `cilevel' != `level' & `outcome' == 1 {
        di as txt "{p 0 1 2 78}(intervals at the estimation level, `lv'%: level() cannot change them " ///
            "for these results){p_end}"
    }
end

* ============================================================================
* Display_chain: chain, the coefficient table in one piece (with bootci the
*   bootstrap interval table): for each pair a row with its name and the values of by() it
*   compares, then its equations indented, without the prefix (overall and with a staged mode
*   marg_expl; with inference also cum_expl, cum_unexpl and raked_focal); then the total alike.
*   After the coefficient table, the explained shares (Display_chain_share) and the footnotes
*   of Display_table
* ============================================================================
program define Display_chain
    syntax [, Level(cilevel) cf(string) pf(string) sf(string) ci(name) cilevel(string) ///
        atest(integer 1) bootci citypes(string) ]
    capture confirm matrix e(V)
    local hasV = (_rc == 0)
    Display_vce
    tempname b V
    matrix `b' = e(b)
    if `hasV' matrix `V' = e(V)
    Vlab depvar "`=word("`e(depvar)'", 1)'" 16
    local J = e(N_pairs)
    // t with the residual (design) df when the prefix set one (svy, jackknife), else z
    local df = e(df_r)
    local st = cond(`df' < ., "t", "z")
    if `df' < . local crit = invttail(`df', (100 - `level') / 200)
    else local crit = invnormal(1 - (100 - `level') / 200)
    // as Display_table: seldecomp's own bootstrap shows the first of the interval types shown
    local ownbs = ("`e(bsnames)'" != "")
    local mt normal
    local ccrit `crit'
    if `ownbs' {
        local mt : word 1 of `citypes'
        if "`mt'" == "normal" local ccrit = invnormal(1 - (100 - `level') / 200)
        else if "`mt'" == "t" local ccrit = invttail(e(df_bs), (100 - `level') / 200)
        else {
            tempname CM
            matrix `CM' = `ci'_ci_`mt'
        }
    }
    local cw = fmtwidth("`cf'")
    local pw = fmtwidth("`pf'")
    local sw = fmtwidth("`sf'")
    Table_headings
    if "`bootci'" != "" {
        local cols `e(bsnames)'
        local lv = string(`cilevel', "%12.0g")
        di as txt _n "Bayesian bootstrap confidence intervals" _col(49) "Replications" _col(67) "= " ///
            as res %10.0fc e(N_reps)
        di as txt "{hline 17}{c TT}{hline 60}"
        di as txt %16s "" " {c |}" _col(42) "Bootstrap"
        di as txt `"{ralign 16:`macval(depvar)'} {c |}"' _col(19) "`h_coefl'" _col(37) "Bias" _col(42) "`h_sel'" ///
            _col(54) "[`lv'% `h_ci']"
    }
    else {
        local cols : colfullnames `b'
        di ""
        di as txt "{hline 17}{c TT}{hline 60}"
        if `hasV' {
            // as Display_table: the variance type, all but its last word on a line above
            local nw : word count `vcehead'
            local vlast : word `nw' of `vcehead'
            local vfirst
            forvalues i = 1/`=`nw' - 1' {
                local vfirst `vfirst' `: word `i' of `vcehead''
            }
            local cilab "[`level'% `h_ci']"
            if !inlist("`mt'", "normal", "t") local cilab "[`cilevel'% `h_ci']"
            local cicol = 79 - udstrlen("`cilab'")
            if "`vfirst'" != "" {
                di as txt "{space 17}{c |}" _col(`=32 + floor((9 - udstrlen("`vfirst'")) / 2)') "`vfirst'"
            }
            local vcol = 32 + floor((9 - udstrlen("`vlast'")) / 2)
            if `vcebs' {
                Ci_label cihd `mt'
                di as txt "{space 17}{c |}" _col(`vcol') "`vlast'" _col(`=69 - floor(udstrlen("`cihd'") / 2)') "`cihd'"
                local hc `"`h_coefl'"'
            }
            else {
                di as txt "{space 17}{c |}" _col(`vcol') "`vlast'"
                local hc `"`h_coef'"'
            }
            di as txt `"{ralign 16:`macval(depvar)'} {c |}"' _col(19) "`hc'" _col(32) "`h_sel'" _col(49) "`st'" ///
                _col(51) "P>|`st'|" _col(`cicol') "`cilab'"
        }
        else di as txt `"{ralign 16:`macval(depvar)'} {c |}"' _col(19) "`h_coef'"
    }
    di as txt "{hline 17}{c +}{hline 60}"
    // the values of by() (user text), for the pairs' rows
    mata: _sd_lab = strtrim(ustrsplit(st_global("e(by_labels)"), "[|]"))
    local prevpan
    local preveq
    local jb 0
    foreach fn of local cols {
        local ++jb
        gettoken eq nm : fn, parse(":")
        gettoken colon nm : nm, parse(":")
        if !regexm("`eq'", "^(p[0-9]+|total)") continue
        local pan = regexs(1)
        local sub overall
        if regexm("`eq'", "^(p[0-9]+|total)_(.+)$") local sub = regexs(2)
        // the cumulative parts and raked means by stage are in the stage tables:
        // with inference only
        if inlist("`sub'", "cum_expl", "cum_unexpl", "raked_focal") & !`hasV' continue
        if "`pan'" != "`prevpan'" {
            if "`prevpan'" != "" di as txt "{hline 17}{c +}{hline 60}"
            local pi = cond("`pan'" == "total", 0, real(substr("`pan'", 2, .)))
            local k1 = cond(`pi' == 0, 1, `pi')
            local k2 = cond(`pi' == 0, `J' + 1, `pi' + 1)
            mata: st_local("ptx", seldecomp_vlab(st_global("e(by)")) + " " + _sd_lab[`k1'] + " to " + _sd_lab[`k2'] + ///
                (`pi' == 0 ? ", the sum of the pairs" : ""))
            mata: st_local("ptx", udstrlen(st_local("ptx")) > 59 ? usubstr(st_local("ptx"), 1, 56) + "..." : st_local("ptx"))
            local ptx : subinstr local ptx "{" "`=char(1)'", all
            local ptx : subinstr local ptx "}" "{c )-}", all
            local ptx : subinstr local ptx "`=char(1)'" "{c -(}", all
            di as res %-17s "`pan'" as txt "{c |} " `"`macval(ptx)'"'
            // values below the rounding threshold of the pair (the total: of every pair) read 0
            Display_eps eps `pi'
            local prevpan `pan'
            local preveq
        }
        if "`sub'" != "`preveq'" {
            Eq_title ttl `sub'
            Eq_lab eql `sub'
            di as res "  " %-15s "`eql'" as txt "{c |}`ttl'"
            local preveq `sub'
        }
        local j = colnumb(`b', "`fn'")
        local coef = `b'[1, `j']
        if abs(`coef') < `eps' local coef 0
        if "`bootci'" != "" {
            // as Display_bootci: bias, std. err., a row per interval type shown
            local x0 `coef'
            local x1 = el(e(bias), 1, `jb')
            local x2 = el(e(se), 1, `jb')
            local k 2
            foreach M of local citypes {
                local x`++k' = el(`ci'_ci_`M', 1, `jb')
                local x`++k' = el(`ci'_ci_`M', 2, `jb')
            }
            forvalues i = 1/`k' {
                if abs(`x`i'') < `eps' local x`i' 0
            }
            local c1 = 30 - `cw'
            local c2 = 41 - `cw'
            local c3 = 51 - `cw'
            local c4 = 62 - `cw'
            local c5 = 73 - `cw'
            local r 0
            foreach M of local citypes {
                local ++r
                local i1 = 1 + 2 * `r'
                local i2 = 2 + 2 * `r'
                Ci_code cc `M'
                if `r' == 1 {
                    Row_lab nml `sub' `nm'
                    di as txt `"{ralign 16:`macval(nml)'} {c |}"' as res _col(`c1') `cf' `x0' ///
                        _col(`c2') `cf' `x1' _col(`c3') `cf' `x2' ///
                        _col(`c4') `cf' `x`i1'' _col(`c5') `cf' `x`i2'' as txt "`cc'"
                }
                else di as txt %16s "" " {c |}" as res _col(`c4') `cf' `x`i1'' _col(`c5') `cf' `x`i2'' as txt "`cc'"
            }
            continue
        }
        Row_lab nml `sub' `nm'
        di as txt `"{ralign 16:`macval(nml)'} {c |}"' as res _col(`=30 - `cw'') `cf' `coef' _c
        if `hasV' {
            // a std. err. below the threshold is 0: the row shows its coefficient only
            local se = sqrt(`V'[`j', `j'])
            if `se' < `eps' local se 0
            if `se' < . & `se' > 0 {
                local z = `coef' / `se'
                if `df' < . local p = 2 * ttail(`df', abs(`z'))
                else local p = 2 * normal(-abs(`z'))
                if inlist("`mt'", "normal", "t") {
                    local ll = `coef' - `ccrit' * `se'
                    local ul = `coef' + `ccrit' * `se'
                }
                else {
                    local cc = colnumb(`CM', "`fn'")
                    local ll = `CM'[1, `cc']
                    local ul = `CM'[2, `cc']
                }
                if abs(`ll') < `eps' local ll 0
                if abs(`ul') < `eps' local ul 0
                di as res _col(`=41 - `cw'') `cf' `se' _col(`=50 - `sw'') `sf' `z' ///
                    _col(`=56 - `pw'') `pf' `p' _col(`=67 - `cw'') `cf' `ll' _col(`=79 - `cw'') `cf' `ul'
            }
            else di ""
        }
        else di ""
    }
    di as txt "{hline 17}{c BT}{hline 60}"
    if "`bootci'" == "" Display_table_note
    mata: mata drop _sd_lab
    if "`bootci'" == "" & `ownbs' & `hasV' & ("`mt'" == "t" | `df' < .) {
        local dfw = cond("`e(clustvar)'" != "", "clusters", "observations")
        di as txt "{p 0 1 2 78}(t: `=e(df_bs)' degrees of freedom, the number of `dfw' less 1){p_end}"
    }
    if "`bootci'" != "" {
        Ci_key key `"`citypes'"'
        di as txt "{p 0 4 2 78}`key'{p_end}"
        if `cilevel' != `level' {
            di as txt "{p 0 1 2 78}(intervals at the estimation level, `lv'%: level() cannot change them " ///
                "for these results){p_end}"
        }
        exit
    }
    Display_chain_share, ci(`ci') cilevel(`cilevel')
    if `hasV' {
        // the footnotes of Display_table: failed replicates at replay (at estimation one message,
        // Chain's note or Misreps_warn)
        local nf = cond(e(N_misreps) < ., e(N_misreps), 0)
        local ni = cond(e(N_iterlimit) < ., e(N_iterlimit), 0)
        local nr = e(N_reps) + `nf'
        if !`atest' & "`e(bsnames)'" != "" & (`nf' > 0 | `ni' > 0) {
            if `nf' > 0 {
                di as txt "{p 0 1 2 78}(`nf' of `nr' replications failed and " cond(`nf' == 1, "was", "were") ///
                    " dropped; bootstrap results may be incorrect){p_end}"
            }
            else {
                di as txt "{p 0 1 2 78}(`ni' of `nr' replications stopped before the weights converged){p_end}"
            }
        }
        else if !`atest' & "`e(bsnames)'" == "" & `nf' > 0 {
            local vce `e(vce)'
            local what = cond("`vce'" == "brr", "BRR", cond("`vce'" == "sdr", "SDR", "`vce'"))
            local dr = cond("`e(prefix)'" == "svy", "", " and " + cond(`nf' == 1, "was", "were") + " dropped")
            di as txt "{p 0 1 2 78}(`nf' of `nr' `what' replicates failed`dr'){p_end}"
        }
        // force: the replications kept although they missed raking targets, at replay
        if !`atest' & e(N_forced) > 0 & e(N_forced) < . {
            di as txt "{p 0 1 2 78}(force kept `=e(N_forced)' of `nr' replications that missed raking targets){p_end}"
        }
    }
end

* Display_chain_share: chain, the explained share of the difference of each pair and of the
*   total, one row each: the share, and with bayesboot() its std. err. and percentile CI (at
*   cilevel, from ci's matrix); n/a when the difference is 0 (Display_eps)
program define Display_chain_share
    syntax [, ci(string) cilevel(string) ]
    tempname b BP
    matrix `b' = e(b)
    matrix `BP' = e(by_pair)
    local J = e(N_pairs)
    local own = ("`e(bsnames)'" != "")
    local haspci 0
    if `own' {
        capture confirm matrix `ci'_pct_explained_ci
        local haspci = (_rc == 0)
    }
    local cse = colnumb(`BP', "pct_explained_se")
    di as txt _n "Explained share of the difference (%)"
    di as txt "{hline 17}{c TT}{hline 60}"
    if `own' {
        di as txt %16s "" " {c |}" _col(25) "Share" _col(32) "Std. err." _col(47) "[`cilevel'% percentile CI]"
    }
    else di as txt %16s "" " {c |}" _col(25) "Share"
    di as txt "{hline 17}{c +}{hline 60}"
    forvalues j = 1/`=`J' + 1' {
        local eq = cond(`j' <= `J', "p`j'", "total")
        Display_eps eps `=cond(`j' <= `J', `j', 0)'
        local jd = colnumb(`b', "`eq':difference")
        local dd = `b'[1, `jd']
        local sh = cond(`j' <= `J', `BP'[`j', 7], e(pct_explained))
        if abs(`dd') < `eps' local sh .
        local se .
        if `own' local se = cond(`j' <= `J', `BP'[`j', `cse'], e(pct_explained_se))
        if missing(`sh') {
            di as txt %16s "`eq'" " {c |}" _col(21) "n/a (the difference is 0)"
            continue
        }
        foreach q in sh se {
            if abs(``q'') < .05 local `q' 0
        }
        di as txt %16s "`eq'" " {c |}" as res _col(21) %9.1f `sh' _c
        if `own' {
            di as res _col(32) %9.1f `se' _c
            if `haspci' {
                local pl = el(`ci'_pct_explained_ci, `j', 1)
                local pu = el(`ci'_pct_explained_ci, `j', 2)
                foreach q in pl pu {
                    if abs(``q'') < .05 local `q' 0
                }
                di as res _col(46) %9.1f `pl' _col(57) %9.1f `pu'
            }
            else di ""
        }
        else di ""
    }
    di as txt "{hline 17}{c BT}{hline 60}"
end


* Pair_title <lmacname> "<title>" <#>: chain, the title of pair #'s stage
*   table or full diagnostics, "<title>: p#, <by> <value> to <value>"; the by() name left out
*   when the title would pass 78 columns, and a wrapped paragraph if it still does
program define Pair_title
    args lname base j
    local v1 = word("`e(by_levels)'", `j')
    local v2 = word("`e(by_levels)'", `j' + 1)
    Vlab byl "`e(by)'"
    local t "`base': p`j', `byl' `v1' to `v2'"
    if udstrlen("`t'") > 78 local t "`base': p`j', `v1' to `v2'"
    if udstrlen("`t'") > 78 local t "{p 0 4 2 78}`t'{p_end}"
    c_local `lname' "`t'"
end

* ============================================================================
* Display_estat: estat bootstrap after bayesboot(), reached from
*   seldecomp_estat.ado through the hidden replay option estatboot(): every element of
*   e(bsnames) in the layout of Stata's own estat bootstrap (_bs_display), with the intervals
*   of Display_ci at the level of the last display (e(level_replay), else e(level)) and the
*   formats and rounding rule of the citable table, so that both show the same digits; r() as
*   a replay at that level leaves it
* ============================================================================
program define Display_estat, rclass
    syntax [, all NORmal Percentile BC BCA noLegend Verbose ]
    mata: seldecomp_vlset(0)
    if "`all'" != "" {
        local normal normal
        local percentile percentile
        local bc bc
        if e(bca) == 1 local bca bca
    }
    if "`normal'`percentile'`bc'`bca'" == "" local bc bc
    local kinds
    foreach k in normal percentile bc bca {
        if "``k''" != "" local kinds `kinds' `k'
    }
    // the level as text: 99.9 has no exact binary form, and level() takes two decimals
    local lv = string(e(level), "%12.0g")
    if e(level_replay) < . local lv = string(e(level_replay), "%12.0g")
    tempname ci b
    Display_ci, level(`lv') prefix(`ci')
    Display_fmts
    local cw = fmtwidth("`cf'")
    matrix `b' = e(b)
    local bsnames `e(bsnames)'
    local kb : word count `bsnames'
    local nout = max(1, wordcount("`e(outcomes)'"))
    di ""
    di as txt "`e(title)'" _c
    di as txt _col(49) "Number of obs" _col(67) "= " as res %10.0fc e(N)
    di as txt _col(49) "Replications" _col(67) "= " as res %10.0f e(N_reps)
    if "`legend'" == "" {
        di ""
        if "`verbose'" != "" {
            // the command line as typed: SMCL braces shown as text, no macro expansion
            mata: st_local("cl", st_global("e(cmdline)"))
            local cl : subinstr local cl "{" "`=char(1)'", all
            local cl : subinstr local cl "}" "{c )-}", all
            local cl : subinstr local cl "`=char(1)'" "{c -(}", all
            di as txt "{p2colset 7 16 20 2}{...}"
            di as txt `"{p2col :Command:}{res:`macval(cl)'}{p_end}"'
            di as txt "{p2colreset}{...}"
        }
    }
    di ""
    if "`e(clustvar)'" != "" & e(N_clust) < . {
        di as txt "{ralign 78:(Replications based on {res:`=strtrim(string(e(N_clust), "%12.0fc"))'} " ///
            "clusters in `e(clustvar)')}"
    }
    // column titles, as _bs_display's (its interval heading for levels with decimals)
    Table_headings
    local depvar
    if wordcount("`e(depvar)'") == 1 local depvar = abbrev("`e(depvar)'", 16)
    local lvs = string(round(`cilevel', .01))
    if strlen("`lvs'") > 2 {
        local sp = cond(strlen("`lvs'") == 4, " ", "")
        local cilab = "   [`sp'`lvs'% " + subinstr("`h_ci'", "nterval", "nt.", 1) + "]"
    }
    else local cilab "  [`lvs'% `h_ci']"
    di as txt "{hline 17}{c TT}{hline 60}"
    di as txt _col(18) "{c |}" _col(42) "Bootstrap"
    di as txt `"{ralign 16:`macval(depvar)'} {c |}"' _col(`=30 - strlen("`h_coefl'")') "`h_coefl'" _col(37) "Bias" ///
        _col(42) "`h_sel'" _col(51) "`cilab'"
    // values end at 29 (coefficient), 40 (bias), 50 (std. err.), 61 and 72 (interval)
    local c1 = 30 - `cw'
    local c2 = 41 - `cw'
    local c3 = 51 - `cw'
    local c4 = 62 - `cw'
    local c5 = 73 - `cw'
    local nper = `kb' / `nout'
    local preveq
    local prevpan
    local prevo -1
    forvalues j = 1/`kb' {
        local fn : word `j' of `bsnames'
        gettoken eq nm : fn, parse(":")
        gettoken colon nm : nm, parse(":")
        // values below the rounding threshold of their outcome are displayed as 0 (Display_eps;
        // chain: of their pair, the total's of every pair, as in the citable table)
        local o = ceil(`j' / `nper')
        if "`e(chain)'" != "" {
            local o 0
            if regexm("`eq'", "^p([0-9]+)") local o = regexs(1)
        }
        if `o' != `prevo' {
            Display_eps eps `o'
            local prevo `o'
        }
        // Keep outcome and chain-pair panels while using the same descriptive equations
        // and row labels as the main and citable displays.
        local deq `eq'
        local pan
        if "`e(chain)'" != "" {
            if regexm("`eq'", "^(p[0-9]+|total)") local pan = regexs(1)
            local deq overall
            if regexm("`eq'", "^(p[0-9]+|total)_(.+)$") local deq = regexs(2)
        }
        else if `nout' > 1 {
            local pan : word `o' of `e(outcomes)'
            if "`eq'" == "`pan'" local deq overall
            else local deq = substr("`eq'", strlen("`pan'") + 2, .)
        }
        if "`eq'" != "`preveq'" {
            di as txt "{hline 17}{c +}{hline 60}"
            if "`pan'" != "`prevpan'" {
                if "`e(chain)'" != "" di as res %-17s "`pan'" as txt "{c |}"
                else di as res %-17s "Outcome" as txt "{c |} `pan'"
                local prevpan `pan'
            }
            Eq_lab eql `deq'
            Eq_title ttl `deq'
            di as res %-17s "`eql'" as txt "{c |}`ttl'"
            local preveq `eq'
        }
        local x0 = `b'[1, colnumb(`b', "`fn'")]
        local x1 = el(e(bias), 1, `j')
        local x2 = el(e(se), 1, `j')
        forvalues i = 0/2 {
            if abs(`x`i'') < `eps' local x`i' 0
        }
        local first 1
        foreach k of local kinds {
            local key = cond("`k'" == "normal", "N", cond("`k'" == "percentile", "P", ///
                cond("`k'" == "bc", "BC", "BCa")))
            local ll = `ci'_ci_`k'[1, `j']
            local ul = `ci'_ci_`k'[2, `j']
            if abs(`ll') < `eps' local ll 0
            if abs(`ul') < `eps' local ul 0
            if `first' {
                Row_lab nml `deq' `nm'
                di as txt `"{ralign 16:`macval(nml)'} {c |}"' as res _col(`c1') `cf' `x0' ///
                    _col(`c2') `cf' `x1' _col(`c3') `cf' `x2' _col(`c4') `cf' `ll' ///
                    _col(`c5') `cf' `ul' as txt _col(73) %6s "(`key')"
            }
            else {
                di as txt %16s "" " {c |}" as res _col(`c4') `cf' `ll' _col(`c5') `cf' `ul' ///
                    as txt _col(73) %6s "(`key')"
            }
            local first 0
        }
    }
    di as txt "{hline 17}{c BT}{hline 60}"
    if "`legend'" == "" Display_table_note
    // the key, as _bs_display's (aligned on the colons)
    local ties
    if e(ties) == 1 local ties " (adjusted for ties) "
    local cp = cond("`bca'" != "", 9, cond("`bc'" != "", 8, 7))
    di as txt "Key:" _c
    if "`normal'" != "" di as txt _col(`=`cp' - 1') "N: Normal"
    if "`percentile'" != "" di as txt _col(`=`cp' - 1') "P: Percentile"
    if "`bc'" != "" di as txt _col(`=`cp' - 2') "BC: Bias-corrected" "`ties'"
    if "`bca'" != "" di as txt _col(`=`cp' - 3') "BCa: Bias-corrected and accelerated" "`ties'"
    if e(N_misreps) > 0 & e(N_misreps) < . {
        local nf = e(N_misreps)
        local nr = e(N_reps) + `nf'
        di as txt "{p 0 1 2 78}(`nf' of `nr' replications failed and " cond(`nf' == 1, "was", "were") ///
            " dropped; bootstrap results may be incorrect){p_end}"
    }
    Display_jkfail
    // r(): as a replay at this level leaves it
    foreach m in ci_normal ci_percentile ci_bc ci_bca pct_explained_ci {
        capture confirm matrix `ci'_`m'
        if _rc continue
        return matrix `m' = `ci'_`m'
    }
    return scalar level = `cilevel'
end

* Note_stages <lmacname> <number of stages> <stage list>: the stages of a note, as
*   "stage 2", "stages 1-3" or "stages 1-2, 4" (empty for a run with one stage)
program define Note_stages
    gettoken lname 0 : 0
    gettoken ns 0 : 0
    local n : word count `0'
    local txt
    if `ns' > 1 & `n' > 0 {
        local a : word 1 of `0'
        local b `a'
        forvalues i = 2/`n' {
            local s : word `i' of `0'
            if `s' == `b' + 1 local b `s'
            else {
                local txt = "`txt'" + cond(`b' > `a', "`a'-`b'", "`a'") + ", "
                local a `s'
                local b `s'
            }
        }
        local txt = cond(`n' > 1, "stages ", "stage ") + "`txt'" + cond(`b' > `a', "`a'-`b'", "`a'")
    }
    c_local `lname' `"`txt'"'
end

program define Display_covdetail
    // for each level of each raking term: the Focal Group and Reference Group observations, the
    // weighted shares of the Focal Group, the Reference Group and the Focal Group under the final
    // raked weights, and the gap between the groups' shares in the order of the difference;
    // values end in columns 36, 44 (observations), 53, 61, 69 (shares) and 78 (gap)
    tempname C
    matrix `C' = e(covariates)
    local rows = rowsof(`C')
    local rnames : rowfullnames `C'
    di as txt _n "Raking variables: weighted shares by level"
    di as txt "{hline 78}"
    di as txt _col(30) "{hline 5} Obs. {hline 4}" _col(47) "{hline 9} Weighted share {hline 7}"
    di as txt "Level" _col(32) "Focal" _col(41) "Ref." _col(49) "Focal" _col(58) "Ref." ///
        _col(65) "Raked" _col(76) "Gap"
    di as txt "{hline 78}"
    local prevt
    forvalues r = 1/`rows' {
        local rn : word `r' of `rnames'
        gettoken t l : rn, parse(":")
        gettoken colon l : l, parse(":")
        if "`t'" != "`prevt'" {
            Vlab tl "`e(term`t')'"
            di as res `"`tl'"'
            local prevt `t'
            // split e(levels#) on the " | " separator in Mata (a label may contain quotes); the
            // cells of an interaction are shown by their components' value labels (or values),
            // joined by " # ", as the term line names the variables
            local tvn = subinstr("`e(term`t')'", "#", " ", .)
            mata: _sd_lab = strtrim(ustrsplit(st_global("e(levels`t')"), "[|]"))
            if `: word count `tvn'' > 1 {
                mata: for (_sd_i = 1; _sd_i <= cols(_sd_lab); _sd_i++) _sd_lab[_sd_i] = ///
                    seldecomp_ixlab(_sd_lab[_sd_i], tokens("`tvn'"))
            }
            local nlab = 0
            mata: st_local("nlab", strofreal(cols(_sd_lab)))
            forvalues ll = 1/`nlab' {
                mata: st_local("lab`ll'", _sd_lab[`ll'])
                local w`ll' = udstrlen(`"`macval(lab`ll')'"')
                local lab`ll' : subinstr local lab`ll' "{" "`=char(1)'", all
                local lab`ll' : subinstr local lab`ll' "}" "{c )-}", all
                local lab`ll' : subinstr local lab`ll' "`=char(1)'" "{c -(}", all
            }
            mata: mata drop _sd_lab
            capture mata: mata drop _sd_i
        }
        // the gap: Focal Group minus Reference Group share, never -0
        local gap = `C'[`r', 1] - `C'[`r', 2]
        if abs(`gap') < .00005 local gap 0
        if `w`l'' > 26 {
            // a long label gets its own line (wrapped) rather than being cut by abbrev()
            di as txt `"{p 2 4 2 78}`macval(lab`l')'{p_end}"'
        }
        else di as txt "  " `"`macval(lab`l')'"' _c
        di as res _col(30) %7.0fc `C'[`r', 4] _col(38) %7.0fc `C'[`r', 5] ///
            _col(47) %7.4f `C'[`r', 1] _col(55) %7.4f `C'[`r', 2] _col(63) %7.4f `C'[`r', 3] ///
            _col(71) %8.4f `gap'
    }
    di as txt "{hline 78}"
    di as txt "{p 0 4 2 78}Focal: focal group. Ref.: reference group. Raked: focal group share under " ///
        "the final raked weights. Gap: Focal minus Reference share (the order of the difference).{p_end}"
end


* ============================================================================
* Make_plot: the plot types of makeplot(): stages (explained and unexplained by stage),
* marginal (the marginal explained part of each stage), waterfall (the steps in units of the
* difference) and waterfalllevels (the steps in levels of the mean)
* ============================================================================
* Plot_spec plotspec: parse makeplot([plottype] [, options]) and refuse its own conflicts (Estimate
* calls it before any work, Make_plot before it draws). The plot type, the first word: stages,
* marginal, waterfall or waterfalllevels (empty: waterfall, or marginal with individually); the
* options after the comma: swapaxis, ci, reverse, chainstep(#), outcome(), format(), blocklabels() and
* twoway(), whose contents go to twoway. Its parts are returned in the caller's locals ps_<name>
program define Plot_spec
    local types stages marginal waterfall waterfalllevels
    local type
    local opts
    mata: st_local("lead", substr(strltrim(st_local("0")), 1, 1))
    if `"`macval(lead)'"' == "," {
        mata: st_local("opts", substr(strltrim(st_local("0")), 2, .))
    }
    else if `"`macval(lead)'"' != "" {
        gettoken type rest : 0, parse(", ") bind quotes
        mata: st_local("isty", strofreal(anyof(tokens(st_local("types")), st_local("type"))))
        if !`isty' {
            mata: st_local("isw", strofreal(ustrregexm(st_local("type"), "^[A-Za-z_][A-Za-z0-9_]*$")))
            local tw
            if `isw' local tw "`type' is not a plot type; "
            di as err "{p 0 4 2}makeplot(): `tw'type the plot type first (stages, marginal, waterfall, or waterfalllevels), then the options after a comma{p_end}"
            exit 198
        }
        mata: st_local("lead", substr(strltrim(st_local("rest")), 1, 1))
        if `"`macval(lead)'"' == "," {
            mata: st_local("opts", substr(strltrim(st_local("rest")), 2, .))
        }
        else if `"`macval(lead)'"' != "" {
            di as err "{p 0 4 2}makeplot(): put a comma between the plot type and the options, as in makeplot(`type', ci){p_end}"
            exit 198
        }
    }
    local 0 `", `macval(opts)'"'
    capture syntax [, SWAPaxis CI REVerse CHAINstep(string) OUTcome(string) FORmat(string) ///
        BLOCKlabels(string asis) TWoway(string asis) NORMal TINTerval PERCentile BC BCA * ]
    if _rc | `"`macval(options)'"' != "" {
        local bad
        if !_rc {
            gettoken bad : options, bind
            mata: st_local("bad", ustrregexra(st_local("bad"), "[(].*$", ""))
            mata: st_local("bad", ustrregexm(st_local("bad"), "^[A-Za-z_][A-Za-z0-9_]*$") ? st_local("bad") + " " : "")
        }
        di as err "{p 0 4 2}makeplot(): option `bad'not allowed; the options are swapaxis, ci, normal, tinterval, percentile, bc, bca, reverse, chainstep(), outcome(), format(), blocklabels(), and twoway(){p_end}"
        exit 198
    }
    local citype `normal' `tinterval' `percentile' `bc' `bca'
    local citype : list retokenize citype
    if `: word count `citype'' > 1 {
        di as err "makeplot(): specify only one interval type: normal, tinterval, percentile, bc, or bca"
        exit 198
    }
    if "`citype'" != "" & "`ci'" == "" {
        di as err "makeplot(): specify ci with an interval type"
        exit 198
    }
    if "`citype'" == "tinterval" local citype t
    // twoway(): by() would draw a graph per group, which the plots are not built for
    local 0 `", `macval(twoway)'"'
    capture syntax [, BY(string asis) * ]
    if !_rc & `"`macval(by)'"' != "" {
        di as err "makeplot(): by() is not allowed inside twoway()"
        exit 198
    }
    // chainstep(#): a positive integer (checked against the chain's steps by the callers)
    if `"`macval(chainstep)'"' != "" {
        capture confirm integer number `chainstep'
        if _rc | real(`"`chainstep'"') < 1 {
            di as err "makeplot(chainstep()) must be a positive integer"
            exit 198
        }
    }
    // outcome(): one name, the outcome drawn (checked against the outcomes by the callers); a list
    // would pass that check and fail as an equation name after the estimation
    if `"`macval(outcome)'"' != "" {
        mata: st_local("onw", strofreal(length(tokens(st_local("outcome")))))
        if `onw' > 1 {
            mata: st_local("on1", tokens(st_local("outcome"))[1])
            di as err `"makeplot(outcome()): `macval(outcome)' names `onw' outcomes; name one"'
            exit 198
        }
    }
    // format(): a numeric display format (%w.de, %w.df, %w.dg, with c) for the values labelled
    if `"`macval(format)'"' != "" {
        mata: st_local("fmtok", strofreal(regexm(st_local("format"), "^%-?0?[0-9]+[.,][0-9]+[efg]c?$")))
        if `fmtok' {
            capture confirm numeric format `format'
            local fmtok = !_rc
        }
        if !`fmtok' {
            mata: st_local("fmtx", subinstr(subinstr(st_local("format"), "{", "{c -(}"), "}", "{c )-}"))
            di as err `"makeplot(format()): `macval(fmtx)' is not a numeric display format"'
            exit 198
        }
    }
    // blocklabels(name "text" ...): pairs of a block name and its text (checked against the
    // block names by Plot_blabels)
    if `"`macval(blocklabels)'"' != "" {
        mata: st_local("bln", strofreal(length(tokens(st_local("blocklabels")))))
        if mod(`bln', 2) {
            di as err `"{p 0 4 2}makeplot(blocklabels()): type pairs of a block name and its text, as in blocklabels(demo "Demographics"){p_end}"'
            exit 198
        }
    }
    foreach m in type swapaxis ci citype reverse chainstep outcome format blocklabels twoway {
        c_local ps_`m' `"`macval(`m')'"'
    }
end

* Plot_blabels, blocks(namelist) [labels(string asis)]: the texts of blocklabels() for the blocks,
* in c_local bl_# (# the block's position) and their lengths in bll_#; the block names when a block
* has no text. Refuses (198) a name that is not a block, a block given twice and an empty text. The
* texts are shown as typed: $, `, { and } are written as {c S|}, {c 96}, {c -(} and {c )-} (twoway
* expands its text options again and reads SMCL tags in them), and every double quote as {c 34}
* (Mata seldecomp_ptxt, which every other text of the graphs goes through as well)
program define Plot_blabels
    syntax , blocks(string) [ labels(string asis) ]
    local nb : word count `blocks'
    forvalues k = 1/`nb' {
        local bl_`k' : word `k' of `blocks'
        local bll_`k' = strlen("`bl_`k''")
        local set_`k' 0
    }
    if `"`macval(labels)'"' != "" {
        mata: st_local("nt", strofreal(length(tokens(st_local("labels")))))
        forvalues j = 1(2)`nt' {
            mata: st_local("nm", tokens(st_local("labels"))[`j'])
            mata: st_local("nmx", seldecomp_ptxt(st_local("nm")))
            mata: st_local("k", strofreal(sum((1..`nb') :* (tokens(st_local("blocks")) :== st_local("nm")))))
            if `k' == 0 {
                di as err `"{p 0 4 2}makeplot(blocklabels()): `macval(nmx)' is not a block; the blocks are `blocks'{p_end}"'
                exit 198
            }
            if `set_`k'' {
                di as err "makeplot(blocklabels()): block `nm' is given twice"
                exit 198
            }
            local set_`k' 1
            mata: st_local("bl_`k'", tokens(st_local("labels"))[`j' + 1])
            // an empty text: twoway would show the row's axis position instead
            if `"`macval(bl_`k')'"' == "" {
                di as err `"{p 0 4 2}makeplot(blocklabels()): the text of `nm' is empty; for a blank label, type " "{p_end}"'
                exit 198
            }
            mata: st_local("bll_`k'", strofreal(ustrlen(st_local("bl_`k'"))))
            // a text that starts with a double quote is shown as typed too: twoway would read it
            // as a list of lines ("Degree" holders as two), so the quote is written {c 34}
            mata: st_local("bl_`k'", seldecomp_ptxt(st_local("bl_`k'")))
        }
    }
    forvalues k = 1/`nb' {
        c_local bl_`k' `"`macval(bl_`k')'"'
        c_local bll_`k' `bll_`k''
    }
end

* Plot_try, twoway_options: the twoway options of a plotspec tried on a graph of two plots, not
* drawn (a plot option such as msymbol() is refused there as in the graphs of makeplot), before
* the estimation; saving() is left out, and so are the user's nodraw and name() (twoway
* refuses them given twice), and the graph is dropped
program define Plot_try
    local 0 `", `macval(0)'"'
    syntax [, SAVing(string asis) NODRAW NAME(string asis) * ]
    if `"`macval(options)'"' == "" exit
    local G _seldecomp_try
    capture noisily twoway (scatteri 1 1) (scatteri 2 2), `macval(options)' nodraw name(`G', replace)
    local rc = _rc
    capture graph drop `G'
    if `rc' {
        di as err "makeplot(twoway()): twoway returned error `rc'"
        exit `rc'
    }
end

// Plot_tick lname lo hi: in local lname, the step of the value axis's ticks over lo to hi: a
// round step (1, 2, 2.5 or 5 times a power of 10) that gives about five of them
program define Plot_tick
    args lname lo hi
    local raw = (`hi' - `lo') / 5
    local mag = 10^floor(log10(`raw'))
    local f = `raw' / `mag'
    c_local `lname' = `mag' * cond(`f' <= 1, 1, cond(`f' <= 2, 2, cond(`f' <= 2.5, 2.5, cond(`f' <= 5, 5, 10))))
end

* Plot_fit lname "twoway options" "title" ...: c_local lname: the first of the titles that fits
* the width of the graph on one line, else the last. title() fills the width at about 60
* characters at the default size of the s1 and s2 schemes (5.5 by 4 inches), so 56 are taken to
* fit there; the text size follows the smaller dimension of the graph, so the room grows as
* xsize / min(xsize, ysize): with the st schemes' 7.5 by 4.5 (the session's scheme, or that of
* scheme() among the twoway options), and with xsize() and ysize() among them
program define Plot_fit
    gettoken lname 0 : 0
    gettoken opts 0 : 0
    local o0 : copy local 0
    local 0 `", `macval(opts)'"'
    capture syntax [, XSIZe(real -1) YSIZe(real -1) SCHeme(string) * ]
    if _rc {
        local xsize -1
        local ysize -1
        local scheme
    }
    if `"`scheme'"' == "" local scheme = c(scheme)
    local xs 5.5
    local ys 4
    if substr(strtrim(`"`scheme'"'), 1, 2) == "st" {
        local xs 7.5
        local ys 4.5
    }
    if `xsize' > 0 local xs `xsize'
    if `ysize' > 0 local ys `ysize'
    local w = 56 * (`xs' / min(`xs', `ys')) / (5.5 / 4)
    local 0 : copy local o0
    local t
    while `"`0'"' != "" {
        gettoken t 0 : 0
        if ustrlen(`"`t'"') <= `w' continue, break
    }
    c_local `lname' `"`t'"'
end

* Plot_ci: limits in e(b) column order, without changing estimation results or rerunning inference.
* A waterfall shifts the interval for the signed bar length to the bar's endpoint; these
* are not confidence intervals for a random cumulative endpoint. Levels plots reverse the
* step's sign and hence exchange/negate its asymmetric limits.
program define Plot_ci
    syntax, Level(cilevel) OUT(name) [ TYPE(string) ]
    capture confirm matrix e(V)
    if _rc {
        di as err "makeplot(ci): these results have no confidence intervals"
        exit 198
    }
    local ownbs = ("`e(bsnames)'" != "")
    if "`type'" == "" {
        if `ownbs' local type : word 1 of `e(citypes)'
        if "`type'" == "" {
            local type = cond((e(df_r) > 0 & e(df_r) < .) | ///
                (e(N_clust) > 0 & e(N_clust) < .), "t", "normal")
        }
    }
    tempname B V C P
    matrix `B' = e(b)
    matrix `V' = e(V)
    local names : colfullnames `B'
    local k = colsof(`B')
    matrix `out' = J(2, `k', .)
    matrix colnames `out' = `names'
    matrix rownames `out' = ll ul
    if inlist("`type'", "normal", "t") {
        local crit = invnormal(1 - (100 - `level') / 200)
        if "`type'" == "t" {
            local df = cond(`ownbs', e(df_bs), e(df_r))
            if !(`df' > 0 & `df' < .) {
                local units = cond(e(N_clust) > 0 & e(N_clust) < ., e(N_clust), e(N))
                local strata = cond(e(N_strata) > 0 & e(N_strata) < ., e(N_strata), 1)
                local df = `units' - `strata'
            }
            if !(`df' > 0 & `df' < .) {
                di as err "makeplot(tinterval): these results have no usable degrees of freedom"
                exit 198
            }
            local crit = invttail(`df', (100 - `level') / 200)
        }
        forvalues j = 1/`k' {
            if `V'[`j', `j'] >= 0 & `V'[`j', `j'] < . {
                matrix `out'[1, `j'] = `B'[1, `j'] - `crit' * sqrt(`V'[`j', `j'])
                matrix `out'[2, `j'] = `B'[1, `j'] + `crit' * sqrt(`V'[`j', `j'])
            }
        }
    }
    else {
        capture confirm matrix e(ci_`type')
        if _rc {
            di as err "makeplot(`type'): these results have no `type' intervals; request them when estimating inference"
            exit 198
        }
        if `ownbs' {
            if `level' != e(level) {
                capture confirm matrix e(bsdraws)
                if _rc {
                    di as err "makeplot(`type'): saved bootstrap draws are required for level(`level')"
                    exit 198
                }
            }
            Display_ci, level(`level') prefix(`P')
            matrix `C' = `P'_ci_`type'
            foreach m in ci_normal ci_t ci_percentile ci_bc ci_bca pct_explained_ci {
                capture matrix drop `P'_`m'
            }
        }
        else {
            if `level' != e(level) {
                di as err "makeplot(`type'): these results store intervals only at level(`=e(level)')"
                exit 198
            }
            matrix `C' = e(ci_`type')
        }
        forvalues j = 1/`k' {
            local c = colnumb(`C', "`: word `j' of `names''")
            if `c' < . {
                matrix `out'[1, `j'] = `C'[1, `c']
                matrix `out'[2, `j'] = `C'[2, `c']
            }
        }
    }
    local lab = cond("`type'" == "normal", "normal", cond("`type'" == "t", "t", ///
        cond("`type'" == "bc", "bias-corrected", cond("`type'" == "bca", "BCa", "percentile"))))
    c_local plotcilabel "`lab'"
end

program define Make_plot
    // The first word selects varlabels, the second is level(#), then the plotspec.
    gettoken vl 0 : 0
    gettoken plotlevel 0 : 0, bind
    local plotspec : copy local 0
    local 0 `", `plotlevel'"'
    syntax [, Level(cilevel) ]
    mata: seldecomp_vlset(`vl')
    Plot_spec `macval(plotspec)'
    foreach m in type swapaxis ci citype reverse chainstep outcome format blocklabels twoway {
        local `m' `"`macval(ps_`m')'"'
    }
    // the twoway options of twoway(); the graph is named seldecomp, replaced by the next one,
    // unless they name it
    local gname name(seldecomp, replace)
    local 0 `", `macval(twoway)'"'
    capture syntax [, NAME(string asis) * ]
    if !_rc & `"`macval(name)'"' != "" local gname
    local options : copy local twoway
    tempname S V b CI
    local ns = e(N_stages)
    local D = e(difference)
    // reverse: the stages of the reverse order (its blocks from the last to the first, its
    // equations rev_)
    local rev = ("`reverse'" != "")
    local rp
    local blocks `e(blocks)'
    if `rev' {
        if "`e(reverse)'" == "" {
            di as err "makeplot(reverse) requires results estimated with sequentially reverse"
            exit 198
        }
        local rp rev_
        local blocks
        forvalues s = `ns'(-1)1 {
            local blocks `blocks' `: word `s' of `e(blocks)''
        }
    }
    // blocklabels(): the text of each block's row (bl_#, its length bll_#)
    local nb0 : word count `blocks'
    Plot_blabels, blocks(`blocks') labels(`macval(blocklabels)')
    // varlabels: a block named after its terms and without a text of blocklabels() is shown by
    // their variable labels
    if `vl' {
        local blks0 `e(blocks)'
        forvalues k = 1/`nb0' {
            local bn : word `k' of `blocks'
            if `"`macval(bl_`k')'"' != "`bn'" continue
            local bk : list posof "`bn'" in blks0
            if `bk' == 0 continue
            local terms `"`e(block`bk')'"'
            if "`bn'" != "`terms'" & "`bn'" != subinstr("`terms'", "#", "_", .) continue
            mata: st_local("bl_`k'", seldecomp_vlab(st_local("terms")))
            mata: st_local("bll_`k'", strofreal(ustrlen(st_local("bl_`k'"))))
            mata: st_local("bl_`k'", seldecomp_ptxt(st_local("bl_`k'")))
        }
    }
    local bl_`=`nb0' + 1' Remainder
    local bll_`=`nb0' + 1' 9
    // format(): the values labelled; by default 3 significant digits
    local fmt = cond(`"`format'"' != "", `"`format'"', "%9.3g")
    // individually: the marginal plot alone (the rakings are not cumulative stages), with a row
    // for the remainder (marg_expl:_remainder)
    local variant = ("`e(stagetype)'" == "individually")
    if `variant' {
        if !inlist("`type'", "", "marginal") {
            di as err "makeplot(`type') not allowed with individually; use marginal"
            exit 198
        }
        local type marginal
        local blocks `blocks' _remainder
    }
    // ci: require fitted inference; Plot_ci below supplies the requested interval family.
    if "`ci'" != "" {
        capture confirm matrix e(V)
        if _rc {
            di as err "makeplot(ci): these results have no confidence intervals"
            exit 198
        }
    }
    matrix `b' = e(b)
    // outcomes(): the first outcome unless outcome() names another; its overall equation is
    // named after it (eq0) and its stage equations are prefixed by its name (pfx)
    local dv `e(depvar)'
    local eq0 overall
    local pfx
    if "`e(outcomes)'" != "" {
        local outs `e(outcomes)'
        if "`outcome'" == "" {
            local outcome : word 1 of `outs'
            di as txt "{p 0 1 2 78}(makeplot: outcome `outcome', the first of the outcomes; " ///
                "makeplot(, outcome(name)) draws another){p_end}"
        }
        else if !`: list outcome in outs' {
            di as err "{p 0 4 2}makeplot(outcome()): `outcome' is not one of the outcomes (`outs'){p_end}"
            exit 198
        }
        local eq0 `outcome'
        local pfx `outcome'_
        local dv `outcome'
        local D = _b[`outcome':difference]
    }
    else if "`outcome'" != "" & "`outcome'" != "`e(depvar)'" {
        di as err "makeplot(outcome()): `outcome' is not the outcome (`e(depvar)')"
        exit 198
    }
    // chain: chainstep(#) draws step # (the pair of the #th and the next value of
    // by()), whose equations are read as an outcome's (p#:, p#_...); without it the whole chain,
    // from the first value to the last (total:, total_marg_expl:), which has no group means and
    // no stages: the waterfall (the default) and the marginal plot
    local ch = ("`e(chain)'" != "")
    local pi 0
    if "`chainstep'" != "" & !`ch' {
        di as err "makeplot(chainstep()) requires results from a chain run"
        exit 198
    }
    if `ch' {
        if "`chainstep'" != "" {
            if real("`chainstep'") > e(N_pairs) {
                di as err "makeplot(chainstep()) must be 1 to `=e(N_pairs)'"
                exit 198
            }
            local pi = real("`chainstep'")
            local eq0 p`pi'
            local pfx p`pi'_
        }
        else {
            if inlist("`type'", "stages", "waterfalllevels") {
                di as err "makeplot(`type') requires chainstep(#) with chain results"
                exit 198
            }
            local eq0 total
            local pfx total_
        }
        local D = _b[`eq0':difference]
    }
    if "`type'" == "" local type waterfall
    // the orientation: marginal draws horizontal bars, the others vertical ones; swapaxis
    // swaps them
    local horizontal = cond(("`type'" == "marginal") != ("`swapaxis'" != ""), "horizontal", "")
    local hasVm 0
    if "`ci'" != "" {
        matrix `V' = e(V)
        local hasVm 1
    }
    local cinote
    if `hasVm' {
        Plot_ci, level(`level') type(`citype') out(`CI')
        local cinote note("`level'% `plotcilabel' intervals", size(vsmall))
    }
    // user text (labels, values) is copied without macro expansion and shown as typed: {, }, $,
    // ` and every " are written as SMCL codes (seldecomp_ptxt), since twoway expands its text
    // options again and reads SMCL tags in them (a label "Single {bf}x" would turn the text
    // bold from x on); the lengths below are those of the text as typed
    // the outcome: the variable label of depvar (the outcome drawn) as it stands, else its
    // name; ratio(): that of the denominator too
    local ymeas mean
    local rv `e(ratio)'
    if "`rv'" != "" local ymeas ratio
    foreach q in dv rv {
        if "``q''" == "" continue
        local `q'l ``q''
        capture confirm variable ``q'', exact
        if !_rc {
            mata: st_local("`q'l", st_varlabel(st_local("`q'")))
            mata: st_local("`q'l", st_local("`q'l") == "" ? st_local("`q'") : st_local("`q'l"))
        }
    }
    // the groups (gr#, as typed): their value labels (e(group_#_label) is "value (label)"), else
    // by = value (in parentheses in the subtitle and title, sr#), or a string value as it stands
    forvalues g = 1/2 {
        mata: st_local("gv", st_global("e(" + (`g' == 1 ? "focal" : "ref") + "_group)"))
        mata: st_local("gl", st_global("e(" + (`g' == 1 ? "focal" : "ref") + "_group_label)"))
        // chain: the pair's values of by(), or for the total the first and the last
        if `ch' {
            local kk = cond(`g' == 1, max(`pi', 1), cond(`pi' > 0, `pi' + 1, e(N_pairs) + 1))
            mata: st_local("gv", tokens(st_global("e(by_levels)"))[`kk'])
            mata: st_local("gl", strtrim(ustrsplit(st_global("e(by_labels)"), "[|]")[`kk']))
        }
        mata: st_local("gr`g'", st_local("gl") != st_local("gv") ? usubstr(st_local("gl"), ///
            ustrlen(st_local("gv")) + 3, ustrlen(st_local("gl")) - ustrlen(st_local("gv")) - 3) : ///
            (strtoreal(st_local("gv")) < . ? seldecomp_vlab(st_global("e(by)")) + " = " + st_local("gv") : st_local("gv")))
        mata: st_local("sr`g'", st_local("gl") == st_local("gv") & strtoreal(st_local("gv")) < . ? ///
            "(" + st_local("gr`g'") + ")" : st_local("gr`g'"))
    }
    // chain, the whole chain: the subtitle says that the difference is the pairs' sum; reverse:
    // that the stages are those of the reverse order
    local ssfx
    if `ch' & `pi' == 0 local ssfx " (sum of `=e(N_pairs)' pairs)"
    if `rev' local ssfx "`ssfx', reverse order"
    // the subtitle names the groups in the order of the difference, "A - B": on one line up to 72
    // characters, since s1mono and s2color fill the width at about 80; otherwise on two lines, "A -"
    // and "B", each group shortened to 70 characters with "..." (a longer line is cut at both
    // edges)
    mata: st_local("n", strofreal(ustrlen(st_local("sr1") + " - " + st_local("sr2") + st_local("ssfx"))))
    if `n' <= 72 {
        mata: st_local("sub", seldecomp_ptxt(st_local("sr1") + " - " + st_local("sr2") + st_local("ssfx")))
        local sub `"`"`macval(sub)'"'"'
    }
    else {
        mata: st_local("s1", seldecomp_ptxt(seldecomp_pfit(st_local("sr1"), 70) + " -"))
        mata: st_local("s2", seldecomp_ptxt(seldecomp_pfit(st_local("sr2"), 70) + st_local("ssfx")))
        local sub `"`"`macval(s1)'"' `"`macval(s2)'"'"'
    }
    // the value axis: "Hourly wage: difference in means" (the label is not put mid-phrase),
    // on two lines beyond 40 characters (seldecomp_pvt); joined in Mata, since the labels
    // are still as typed here and a one-line if is parsed twice, which dropped a `name' or a
    // $name of the labels with ratio()
    mata: st_local("yt0", st_local("dvl") + (st_local("rv") != "" ? " / " + st_local("rvl") : ""))
    mata: st_local("ytq", seldecomp_pvt(st_local("yt0"), "difference in `ymeas's"))
    // value labels below the rounding threshold of the outcome read 0, and there are no shares
    // when the difference itself does (Display_eps; the bars are drawn from e(b) as it stands)
    local oi 1
    if "`e(outcomes)'" != "" local oi : list posof "`outcome'" in outs
    if `ch' local oi `pi'
    Display_eps eps `oi'
    local Dz = (abs(`D') < `eps')
    local Dtxt = strtrim(string(cond(`Dz', 0, `D'), "`fmt'"))
    if regexm("`Dtxt'", "^-[0.,]+$") local Dtxt = substr("`Dtxt'", 2, .)

    if "`type'" == "waterfalllevels" {
        // the waterfall in levels: the Focal Group mean, a step for each block
        // (and with individually the remainder) to the raked Focal Group mean, then the
        // unexplained step to the Reference Group mean (ratios with ratio()); all from e(b). A step is
        // the change in the mean, minus its value in e(b) (the difference is the Focal Group - Reference Group),
        // and the step's share of the difference is its value over the difference. Spikes where
        // e(V) has the variance of a mean (around it) or of a step (around the step's end, as in
        // the waterfall)
        local nb : word count `blocks'
        local nr = `nb' + 4
        local names `eq0':focal_group
        foreach bn of local blocks {
            local names `names' `pfx'`rp'marg_expl:`bn'
        }
        local names `names' `eq0':raked_focal `eq0':unexplained `eq0':ref_group
        preserve
        qui drop _all
        qui set obs `nr'
        qui gen double pos = _n
        qui gen double val = .
        qui gen double se = .
        qui gen double cil = .
        qui gen double ciu = .
        // kind: 1 a mean, 3 block, 5 remainder, 4 unexplained
        qui gen byte kind = cond(inlist(_n, 1, `nb' + 2, `nr'), 1, cond(_n == `nb' + 3, 4, 3))
        if `variant' qui replace kind = 5 in `=`nb' + 1'
        // the row names are angled in the vertical view beyond seven bars (four blocks, as in the
        // waterfall in units of the difference), when a block name is longer than 10 characters,
        // or when a group name is longer than 7: the rows of the raked Focal Group mean and of the
        // Reference Group mean sit beside "Unexplained" (11), and 10 beside it collide at seven bars
        // under s1mono and s2color; the raked Focal Group row is named on one line when angled,
        // otherwise on two; a block's row is named by its text of blocklabels(), else by its name.
        // The names of the groups in the rows (the title names them in full): angled, the first
        // row's name, that of the Focal Group mean, is cut at the left edge from about 21 characters
        // with one block, 17 with three, 15 with eight to twelve (under s1mono and s2color), and
        // twoway does not angle a name of two lines, so the group names of the three rows are
        // shortened with "..." to 18 characters up to five rows (one block), 15 up to eight (four
        // blocks), 13 up to twelve, else 12 (the remainder row of individually counts as a
        // block; measured with value-axis labels of one character); in the horizontal view long
        // names narrow the plot until the labels of the bars are cut at the right edge, so a
        // name longer than 20 characters is put on two lines, and shortened beyond
        // (seldecomp_pwrap), up to nine rows: beyond, the two lines run into the neighbouring
        // rows' names (touching at ten rows, over them from eleven, under every scheme), so
        // the name is one line cut to 20; the raked row repeats the Focal Group name with ", raked"
        local len 0
        forvalues k = 1/`nb' {
            local len = max(`len', `bll_`k'')
        }
        mata: st_local("glen", strofreal(max((ustrlen(st_local("gr1")), ustrlen(st_local("gr2"))))))
        local ang = cond(`nr' > 7 | `len' > 10 | `glen' > 7, "angle(30)", "")
        if "`horizontal'" != "" & `nr' <= 9 {
            mata: st_local("rw1", seldecomp_pwrap(st_local("gr1"), 20, ""))
            mata: st_local("rw2", seldecomp_pwrap(st_local("gr2"), 20, ""))
            mata: st_local("rwr", seldecomp_pwrap(st_local("gr1"), 20, ", raked"))
        }
        else {
            local wl = cond("`horizontal'" != "", 20, cond(`nr' <= 5, 18, cond(`nr' <= 8, 15, cond(`nr' <= 12, 13, 12))))
            mata: st_local("rw1", seldecomp_ptxt(seldecomp_pfit(st_local("gr1"), `wl')))
            mata: st_local("rw2", seldecomp_ptxt(seldecomp_pfit(st_local("gr2"), `wl')))
            mata: st_local("rwr", seldecomp_ptxt(seldecomp_pfit(st_local("gr1"), `wl') + ", raked"))
        }
        mata: st_local("rw0", seldecomp_ptxt(st_local("gr1")))
        local tl
        forvalues r = 1/`nr' {
            local nm : word `r' of `names'
            local c = colnumb(`b', "`nm'")
            if `c' == . {
                di as err "makeplot: `nm' not found in e(b)"
                exit 111
            }
            qui replace val = `b'[1, `c'] in `r'
            if `hasVm' {
                local c = colnumb(`CI', "`nm'")
                if `c' < . {
                    qui replace cil = `CI'[1, `c'] in `r'
                    qui replace ciu = `CI'[2, `c'] in `r'
                    if `CI'[1, `c'] < . & `CI'[2, `c'] < . {
                        qui replace se = sqrt(`V'[`c', `c']) in `r'
                    }
                }
            }
            // the rows of the means are named by the groups (texts in compound quotes)
            if `r' == 1 local tl `"`macval(tl)' `r' `"`macval(rw1)'"'"'
            else if `r' <= `nb' + 1 {
                local tl `"`macval(tl)' `r' `"`macval(bl_`=`r' - 1')'"'"'
            }
            else if `r' == `nb' + 2 {
                if "`horizontal'" == "" & "`ang'" == "" {
                    local tl `"`macval(tl)' `r' `"`"`macval(rw0)'"' "raked""'"'
                }
                else local tl `"`macval(tl)' `r' `"`macval(rwr)'"'"'
            }
            else if `r' == `nb' + 3 local tl `"`macval(tl)' `r' "Unexpl.""'
            else local tl `"`macval(tl)' `r' `"`macval(rw2)'"'"'
        }
        // the steps (the change in the mean) and the bar ends: a step from the running level
        // (the Focal Group mean plus the steps before it), the unexplained step from the raked mean
        // to the Reference Group mean; the means' bars from the base, set below with the axis
        qui gen double stp = -val if kind != 1
        qui gen double w1 = val[1] in 1
        qui replace w1 = w1[_n - 1] + stp in 2/`=`nb' + 1'
        qui gen double w0 = w1[_n - 1] in 2/`=`nb' + 1'
        qui replace w1 = val in `=`nb' + 2'
        qui replace w0 = val[`nb' + 2] in `=`nb' + 3'
        qui replace w1 = val[`nr'] in `=`nb' + 3'
        qui replace w1 = val in `nr'
        // spikes around the bar ends, with ci (the means' around them, the steps' around the
        // step's end); labels: values below the rounding threshold read 0
        qui gen double ll = cond(kind == 1, cil, w1 + val - ciu)
        qui gen double ul = cond(kind == 1, ciu, w1 + val - cil)
        qui gen double vz = cond(abs(stp) < `eps', 0, stp) if kind != 1
        qui replace vz = cond(abs(val) < `eps', 0, val) if kind == 1
        // the means: 4 significant digits, and at least the decimals that show the smallest gap
        // between them (the explained and unexplained parts and the difference, rounding error
        // left out) to one significant digit, so that means that differ read differently (not
        // beyond 10 significant digits: the step labels then show the gap alone). Zeros in the
        // decimals beyond those of the gap are dropped, as %g does, and so is a leading 0;
        // thousands are separated when a mean reaches 10,000 (e-notation from 1e15, and where a
        // mean would need more than 20 decimals: below 1e-17, or a gap below 1e-20). The three
        // then share the most decimals any of them keeps, like a table column (3.800 beside
        // 3.815, not 3.80), but a mean gets back only the zeros it dropped: no more decimals than
        // its own 4 significant digits or the gap's (250.3 beside .0213, not 250.3000); a mean
        // that reads 0 stays 0
        local mx = max(abs(val[1]), abs(val[`nb' + 2]), abs(val[`nr']))
        local q1 = abs(val[`nb' + 2] - val[1])
        local q2 = abs(val[`nr'] - val[`nb' + 2])
        local q3 = abs(val[`nr'] - val[1])
        local gp = .
        forvalues j = 1/3 {
            if `q`j'' >= `eps' local gp = min(`gp', `q`j'')
        }
        local dg 0
        if `gp' < . {
            local dg = max(0, -floor(log10(`gp')))
            if floor(log10(`mx')) + 1 + `dg' > 10 local dg 0
        }
        local sep = cond(`mx' >= 9999.5, "c", "")
        qui gen lab = ""
        foreach r in 1 `=`nb' + 2' `nr' {
            local d = max(0, 3 - floor(log10(abs(vz[`r']))), `dg')
            local dm`r' `d'
            if vz[`r'] == 0 local s 0
            else if `mx' >= 1e15 | `d' > 20 local s = strtrim(string(vz[`r'], "%9.4g"))
            else {
                local s = strtrim(string(vz[`r'], "%24.`d'f`sep'"))
                local k = `d' - `dg'
                while `k' > 0 & substr("`s'", -1, 1) == "0" {
                    local s = substr("`s'", 1, strlen("`s'") - 1)
                    local --k
                }
                if substr("`s'", -1, 1) == "." local s = substr("`s'", 1, strlen("`s'") - 1)
                if substr("`s'", 1, 2) == "0." local s = substr("`s'", 2, .)
                else if substr("`s'", 1, 3) == "-0." local s = "-" + substr("`s'", 3, .)
            }
            qui replace lab = "`s'" in `r'
        }
        local dcol 0
        if `mx' < 1e15 {
            foreach r in 1 `=`nb' + 2' `nr' {
                local k = strpos(lab[`r'], ".")
                if `k' & !strpos(lab[`r'], "e") local dcol = max(`dcol', strlen(lab[`r']) - `k')
            }
            foreach r in 1 `=`nb' + 2' `nr' {
                local k = strpos(lab[`r'], ".")
                local k = cond(`k', strlen(lab[`r']) - `k', 0)
                local t = min(`dcol', `dm`r'')
                if vz[`r'] != 0 & `k' < `t' & !strpos(lab[`r'], "e") {
                    local s = strtrim(string(vz[`r'], "%24.`t'f`sep'"))
                    if substr("`s'", 1, 2) == "0." local s = substr("`s'", 2, .)
                    else if substr("`s'", 1, 3) == "-0." local s = "-" + substr("`s'", 3, .)
                    qui replace lab = "`s'" in `r'
                }
            }
        }
        // format(): the means in that format instead
        if `"`format'"' != "" {
            qui replace lab = strtrim(string(vz, "`fmt'")) if kind == 1
            qui replace lab = substr(lab, 2, .) if kind == 1 & regexm(lab, "^-[0.,]+$")
        }
        // the steps: 3 significant digits (or format()), a value that reads 0 without a sign; the
        // share of a step below the rounding threshold is 0, as in the waterfall in units of the
        // difference
        qui replace lab = strtrim(string(vz, "`fmt'")) if kind != 1
        qui replace lab = substr(lab, 2, .) if kind != 1 & regexm(lab, "^-[0.,]+$")
        qui replace lab = "+" + lab if kind != 1 & vz > 0 & !regexm(lab, "^[0.,]+$")
        if !`Dz' {
            qui replace lab = lab + " (" + string(cond(vz == 0 | abs(100 * val / `D') < .5, 0, ///
                round(100 * val / `D')), "%9.0f") + "%)" if kind != 1
        }
        // the value axis: the levels, the spikes and the labels, with a margin, and the numbers
        // of a value-axis option of the user's (the range() of yscale(); the values of ylabel(),
        // ytick(), ymlabel() and ymtick(), their texts left out; x... with horizontal), which
        // twoway would add to the axis. When 0 is among the levels, spikes and those numbers the
        // means' bars start at 0; otherwise at the end of the axis nearest 0 (which the ticks may
        // put at 0)
        local lo = .
        local hi = .
        foreach v in w0 w1 ll ul {
            summarize `v', meanonly
            local lo = min(`lo', r(min))
            local hi = max(`hi', r(max))
        }
        local rg = `hi' - `lo'
        if `rg' == 0 local rg = max(abs(`hi'), 1)
        local mlen 0
        forvalues r = 1/`nr' {
            local mlen = max(`mlen', strlen(lab[`r']))
        }
        // (horizontal: long row names narrow the plot, so that a character of a label spans more
        // of the axis; the room per character grows as the plot's width in characters of the
        // names, about 83 less the longest name, rn, falls below that beside short names, 70; the
        // labels of the bars are otherwise cut at the right edge beside group names of 26 characters)
        local hf 1
        if "`horizontal'" != "" {
            mata: st_local("rn", strofreal(max((11, min((ustrlen(st_local("gr1")), 20)) + 7, ///
                min((ustrlen(st_local("gr2")), 20))))))
            forvalues k = 1/`nb' {
                local rn = max(`rn', `bll_`k'')
            }
            local hf = max(1, 70 / max(83 - `rn', 20))
        }
        local pad = `rg' * cond("`horizontal'" != "", .02 + .018 * `mlen' * `hf', .08)
        local lo0 = `lo'
        local hi0 = `hi'
        local lo = `lo' - `pad'
        local hi = `hi' + `pad'
        if `lo0' >= 0 local lo = max(`lo', 0)
        if `hi0' <= 0 local hi = min(`hi', 0)
        local va = cond("`horizontal'" != "", "x", "y")
        local VA = upper("`va'")
        local u_opts : copy local options
        local u_rest : copy local options
        local u_num
        while `"`macval(u_rest)'"' != "" {
            local 0 `", `macval(u_rest)'"'
            capture syntax [, `VA'SCale(string asis) `VA'LABel(string asis) `VA'TIck(string asis) ///
                `VA'MLABel(string asis) `VA'MTIck(string asis) * ]
            if _rc continue, break
            local u_rest : copy local options
            local u_hit 0
            if `"`macval(`va'scale)'"' != "" {
                local u_hit 1
                local 0 `", `macval(`va'scale)'"'
                capture syntax [, Range(numlist) * ]
                if !_rc local u_num `u_num' `range'
            }
            foreach o in label tick mlabel mtick {
                if `"`macval(`va'`o')'"' == "" continue
                local u_hit 1
                local u_s : copy local `va'`o'
                gettoken u_rule : u_s, parse(",") quotes bind
                if `"`macval(u_rule)'"' == "," local u_rule
                local u_nl
                while `"`macval(u_rule)'"' != "" {
                    gettoken u_t u_rule : u_rule, quotes bind
                    if !inlist(substr(`"`macval(u_t)'"', 1, 1), char(34), char(96)) {
                        local u_nl `u_nl' `u_t'
                    }
                }
                capture numlist "`u_nl'"
                if !_rc local u_num `u_num' `r(numlist)'
            }
            if !`u_hit' continue, break
        }
        local options : copy local u_opts
        foreach u of local u_num {
            local lo0 = min(`lo0', `u')
            local hi0 = max(`hi0', `u')
            local lo = min(`lo', `u')
            local hi = max(`hi', `u')
        }
        if `hi' <= `lo' local hi = `lo' + 1
        // ticks at a round step (1, 2, 2.5 or 5 times a power of 10), about five of them
        Plot_tick tk `lo' `hi'
        local t0 = floor(`lo' / `tk' + 1e-9) * `tk'
        local t1 = ceil(`hi' / `tk' - 1e-9) * `tk'
        if `lo0' <= 0 & `hi0' >= 0 local base 0
        else if `lo0' > 0 local base = `t0'
        else local base = `t1'
        // the side of a label: that of the step, or of the mean from the base
        qui replace w0 = `base' if kind == 1
        qui gen double dir = cond(kind == 1, w1 - `base', stp)
        qui gen double anc = cond(se < ., cond(dir >= 0, max(w1, ul), min(w1, ll)), w1)
        local p1 = cond("`horizontal'" != "", 3, 12)
        local p2 = cond("`horizontal'" != "", 9, 6)
        qui gen byte lp = cond(dir >= 0, `p1', `p2')
        // (horizontal: a label beyond a bar's end near an end of the axis can run past it, under
        // s1mono and s2color over the row names at the left, since the ticks leave it whatever
        // room the rounding gives while its width in axis units grows with the axis's span; the
        // axis is extended by whole ticks until every label fits, a character taken as
        // .95 / (83 - rn) of the plot's width, measured under s1mono, plus one for the gap; an
        // end at 0 stays there, as the levels are then all on one side of it, so a label may
        // still run past it; the base and the means' starts follow the new end (dir is read for
        // its sign only, which does not change))
        if "`horizontal'" != "" {
            local cw = .95 / max(83 - `rn', 20)
            tempvar need
            qui gen double `need' = .
            forvalues it = 1/8 {
                qui replace `need' = (strlen(lab) + 1) * `cw' * (`t1' - `t0')
                qui count if lp == `p2' & anc - `need' < `t0'
                local xl = r(N) > 0 & `t0' != 0
                qui count if lp == `p1' & anc + `need' > `t1'
                local xr = r(N) > 0 & `t1' != 0
                if !`xl' & !`xr' continue, break
                if `xl' local t0 = cond(`t0' > 0, max(`t0' - `tk', 0), `t0' - `tk')
                if `xr' local t1 = cond(`t1' < 0, min(`t1' + `tk', 0), `t1' + `tk')
                // beyond eight ticks the step follows the wider span, the ends rounded to it
                // (an end at 0 stays): with the step of the span before the extension, an axis
                // extended many times (beside 50-character blocklabels() texts a character
                // spans about 3% of the axis) would carry overlapping tick labels; a single
                // extension keeps its ticks
                if round((`t1' - `t0') / `tk') + 1 > 8 {
                    Plot_tick tk `t0' `t1'
                    local t0 = floor(`t0' / `tk' + 1e-9) * `tk'
                    local t1 = ceil(`t1' / `tk' - 1e-9) * `tk'
                }
            }
            if `lo0' > 0 local base = `t0'
            else if `hi0' < 0 local base = `t1'
            qui replace w0 = `base' if kind == 1
        }
        qui gen double cy = w1 if _n < `nr'
        qui gen double c0 = pos + .3
        qui gen double c1 = pos + .7
        local ci
        qui count if se < .
        if r(N) {
            local ci (rcap ll ul pos if se < ., `=cond("`horizontal'" != "", "horizontal", "")' ///
                lcolor(foreground) lwidth(thin))
        }
        // the labels' shares: a share is the step over the difference, as
        // in the stage table, which is the step's part of the way from the Focal Group mean to the
        // Reference Group mean in either order of the difference; so on the way down a fall has a
        // positive share and a rise a negative one
        local bars
        // the blocks (one plot), the remainder, the unexplained step, the three means
        foreach kd in 3 5 4 1 {
            if `kd' == 5 & !`variant' continue
            if `kd' == 3 local sty pstyle(p1bar)
            else if `kd' == 5 local sty pstyle(p1bar) fintensity(35)
            else if `kd' == 4 local sty fcolor(gs12) lcolor(gs8) lwidth(thin)
            else local sty fcolor(none) lcolor(foreground) lwidth(medthin)
            local bars `bars' (rbar w0 w1 pos if kind == `kd', barwidth(.6) `sty' ///
                `=cond("`horizontal'" != "", "horizontal", "")')
        }
        mata: st_local("ytq", seldecomp_pvt(st_local("yt0"), "`ymeas'"))
        if "`horizontal'" != "" {
            local con (pcspike c0 cy c1 cy if cy < ., lcolor(gs8) lpattern(shortdash) lwidth(thin))
            local labs (scatter pos anc, msymbol(none) mlabel(lab) mlabvposition(lp) mlabsize(small) ///
                mlabcolor(foreground))
            local axes xline(`base', lcolor(gs8) lwidth(thin) lpattern(solid)) ///
                xscale(range(`t0' `t1')) xlabel(`t0'(`tk')`t1') ///
                ylabel(`macval(tl)', noticks nogrid angle(0)) yscale(reverse range(.5 `=`nr' + .5'))
            local tts xtitle(`macval(ytq)') ytitle("")
        }
        else {
            local con (pcspike cy c0 cy c1 if cy < ., lcolor(gs8) lpattern(shortdash) lwidth(thin))
            // smaller labels beyond seven bars, as the row names are angled
            local msz = cond(`nr' > 7, "vsmall", "small")
            local labs (scatter anc pos, msymbol(none) mlabel(lab) mlabvposition(lp) mlabsize(`msz') ///
                mlabcolor(foreground))
            local axes yline(`base', lcolor(gs8) lwidth(thin) lpattern(solid)) ///
                yscale(range(`t0' `t1')) ylabel(`t0'(`tk')`t1', angle(0)) ///
                xlabel(`macval(tl)', noticks nogrid `ang') xscale(range(.5 `=`nr' + .5'))
            local tts ytitle(`macval(ytq)') xtitle("")
        }
        // the title names the groups: on one line up to 56 characters (as typed), with
        // ", stage by stage" when that fits too, since s1mono and s2color fill the width at
        // about 60; beyond, on two lines, "From the <focal group> mean" and "to the <reference
        // group> mean", each group shortened to 42 characters with "..." (a line of 60
        // characters is cut at both edges)
        mata: st_local("tlen", strofreal(ustrlen("From the " + st_local("sr1") + " `ymeas' to the " + ///
            st_local("sr2") + " `ymeas'")))
        if `tlen' <= 56 {
            mata: st_local("ttl1", seldecomp_ptxt("From the " + st_local("sr1") + " `ymeas'"))
            mata: st_local("ttl2", seldecomp_ptxt("to the " + st_local("sr2") + " `ymeas'"))
        }
        else {
            mata: st_local("ttl1", seldecomp_ptxt("From the " + seldecomp_pfit(st_local("sr1"), 42) + " `ymeas'"))
            mata: st_local("ttl2", seldecomp_ptxt("to the " + seldecomp_pfit(st_local("sr2"), 42) + " `ymeas'"))
        }
        if `tlen' + 16 <= 56 {
            local ttl `"`"`macval(ttl1)' `macval(ttl2)', stage by stage"'"'
        }
        else if `tlen' <= 56 {
            local ttl `"`"`macval(ttl1)' `macval(ttl2)'"'"'
        }
        else {
            local ttl `"`"`macval(ttl1)'"' `"`macval(ttl2)'"'"'
        }
        // (the title names the groups, so a subtitle only says that the stages are the reverse order's)
        local rsub
        if `rev' local rsub subtitle("Reverse order", span)
        capture noisily twoway `bars' `con' `ci' `labs', `macval(axes)' `macval(tts)' ///
            title(`macval(ttl)', span) `rsub' legend(off) `cinote' `gname' `macval(options)'
    }
    else if "`type'" == "waterfall" {
        // the waterfall, in units of the difference: one bar per block (marg_expl:<block>,
        // with individually each block raked alone and then the remainder), each starting
        // where the one before ends (the running total of the steps from 0), then the explained
        // part from 0, the unexplained part from the explained part to the difference, and the
        // total difference from 0; all from e(b). Spikes where e(V) has the variance of a bar's
        // value (its length), around the bar's end; labels beyond the end on the side of the step
        local nb : word count `blocks'
        local nr = `nb' + 3
        local names
        foreach bn of local blocks {
            local names `names' `pfx'`rp'marg_expl:`bn'
        }
        local names `names' `eq0':explained `eq0':unexplained `eq0':difference
        preserve
        qui drop _all
        qui set obs `nr'
        qui gen double pos = _n
        qui gen double val = .
        qui gen double se = .
        qui gen double cil = .
        qui gen double ciu = .
        // kind: 3 block, 5 remainder, 2 explained, 4 unexplained, 1 total difference
        qui gen byte kind = cond(_n <= `nb', 3, cond(_n == `nb' + 1, 2, cond(_n == `nb' + 2, 4, 1)))
        if `variant' qui replace kind = 5 in `nb'
        // the row names are angled in the vertical view beyond six bars or when a block's name (its
        // text of blocklabels()) is longer than 10 characters
        local len 0
        forvalues k = 1/`nb' {
            local len = max(`len', `bll_`k'')
        }
        local ang = cond("`horizontal'" == "" & (`nr' > 6 | `len' > 10), "angle(30)", "")
        local tl
        forvalues r = 1/`nr' {
            local nm : word `r' of `names'
            local c = colnumb(`b', "`nm'")
            if `c' == . {
                di as err "makeplot: `nm' not found in e(b)"
                exit 111
            }
            qui replace val = `b'[1, `c'] in `r'
            if `hasVm' {
                local c = colnumb(`CI', "`nm'")
                if `c' < . {
                    qui replace cil = `CI'[1, `c'] in `r'
                    qui replace ciu = `CI'[2, `c'] in `r'
                    if `CI'[1, `c'] < . & `CI'[2, `c'] < . {
                        qui replace se = sqrt(`V'[`c', `c']) in `r'
                    }
                }
            }
            // (brace blocks: the texts of blocklabels() may hold quotes)
            local rl `"`macval(bl_`r')'"'
            if `r' == `nb' + 1 local rl Expl.
            else if `r' == `nb' + 2 local rl Unexpl.
            else if `r' == `nr' local rl Total diff
            local tl `"`macval(tl)' `r' `"`macval(rl)'"'"'
        }
        // the bar ends: the running total of the steps, then E, then from E to D, then D
        qui gen double w1 = sum(val) if _n <= `nb'
        qui gen double w0 = cond(_n == 1, 0, w1[_n - 1]) if _n <= `nb'
        qui replace w0 = 0 in `=`nb' + 1'
        qui replace w1 = val in `=`nb' + 1'
        qui replace w0 = val[`nb' + 1] in `=`nb' + 2'
        qui replace w1 = val[`nr'] in `=`nb' + 2'
        qui replace w0 = 0 in `nr'
        qui replace w1 = val in `nr'
        // spikes, with ci: around each bar's end (the explained part and the total difference
        // start at 0). The step of each bar (its signed length, the value labelled) sets the side
        // of its label
        qui gen double ll = w1 - val + cil
        qui gen double ul = w1 - val + ciu
        qui gen double anc = cond(se < ., cond(val >= 0, max(w1, ul), min(w1, ll)), w1)
        local p1 = cond("`horizontal'" != "", 3, 12)
        local p2 = cond("`horizontal'" != "", 9, 6)
        qui gen byte lp = cond(val >= 0, `p1', `p2')
        qui gen double vz = cond(abs(val) < `eps', 0, val)
        // the values: 3 significant digits, or format(); a value that reads 0 without a sign
        qui gen lab = strtrim(string(vz, "`fmt'"))
        qui replace lab = substr(lab, 2, .) if regexm(lab, "^-[0.,]+$")
        qui replace lab = "+" + lab if inlist(kind, 3, 4, 5) & vz > 0 & !regexm(lab, "^[0.,]+$")
        if !`Dz' {
            qui replace lab = lab + " (" + string(cond(abs(100 * vz / `D') < .5, 0, ///
                round(100 * vz / `D')), "%9.0f") + "%)" if kind != 1
        }
        // connectors: from each bar's end to the next bar, at the level of that end
        qui gen double cy = w1 if _n < `nr'
        qui gen double c0 = pos + .3
        qui gen double c1 = pos + .7
        local ci
        qui count if se < .
        if r(N) {
            local ci (rcap ll ul pos if se < ., `=cond("`horizontal'" != "", "horizontal", "")' ///
                lcolor(foreground) lwidth(thin))
        }
        // room for the labels beyond the bar ends
        local lo = min(0, `D')
        local hi = max(0, `D')
        foreach v in w0 w1 anc {
            summarize `v', meanonly
            local lo = min(`lo', r(min))
            local hi = max(`hi', r(max))
        }
        local mlen 0
        forvalues r = 1/`nr' {
            local mlen = max(`mlen', strlen(lab[`r']))
        }
        local pad = (`hi' - `lo') * cond("`horizontal'" != "", .02 + .018 * `mlen', .08)
        qui count if val < 0
        if r(N) local lo = `lo' - `pad'
        qui count if val >= 0
        if r(N) local hi = `hi' + `pad'
        local bars
        // the blocks and the explained part (one plot), the remainder, the unexplained part, the
        // total difference
        foreach kd in 3 5 4 1 {
            if `kd' == 5 & !`variant' continue
            local kif kind == `kd'
            if `kd' == 3 local kif inlist(kind, 2, 3)
            if `kd' == 3 local sty pstyle(p1bar)
            else if `kd' == 5 local sty pstyle(p1bar) fintensity(35)
            else if `kd' == 4 local sty fcolor(gs12) lcolor(gs8) lwidth(thin)
            else local sty fcolor(none) lcolor(foreground) lpattern(dash)
            local bars `bars' (rbar w0 w1 pos if `kif', barwidth(.6) `sty' ///
                `=cond("`horizontal'" != "", "horizontal", "")')
        }
        if "`horizontal'" != "" {
            local con (pcspike c0 cy c1 cy if cy < ., lcolor(gs8) lpattern(shortdash) lwidth(thin))
            local labs (scatter pos anc, msymbol(none) mlabel(lab) mlabvposition(lp) mlabsize(small) ///
                mlabcolor(foreground))
            local axes xline(0, lcolor(gs8) lwidth(thin) lpattern(solid)) xscale(range(`lo' `hi')) ///
                ylabel(`macval(tl)', noticks nogrid angle(0)) yscale(reverse range(.5 `=`nr' + .5'))
            local tts xtitle(`macval(ytq)') ytitle("")
        }
        else {
            local con (pcspike cy c0 cy c1 if cy < ., lcolor(gs8) lpattern(shortdash) lwidth(thin))
            // beyond six bars the labels of neighbouring bars at similar heights would touch
            local msz = cond(`nr' > 6, "vsmall", "small")
            local labs (scatter anc pos, msymbol(none) mlabel(lab) mlabvposition(lp) mlabsize(`msz') ///
                mlabcolor(foreground))
            local axes yline(0, lcolor(gs8) lwidth(thin) lpattern(solid)) yscale(range(`lo' `hi')) ///
                ylabel(, angle(0)) xlabel(`macval(tl)', noticks nogrid `ang') xscale(range(.5 `=`nr' + .5'))
            local tts ytitle(`macval(ytq)') xtitle("")
        }
        capture noisily twoway `bars' `con' `ci' `labs', `macval(axes)' `macval(tts)' ///
            title("Decomposition of the difference by raking stage", span) ///
            subtitle(`macval(sub)', span) legend(off) `cinote' `gname' `macval(options)'
    }
    else if "`type'" == "marginal" {
        // rows from the top (from the left with swapaxis): the total difference, the explained
        // part, the marginal explained part of each stage (marg_expl:<block>, the change in the
        // explained part when the block is added; with individually each block raked alone, then
        // the remainder) and the unexplained part, all from e(b); with ci their CI spikes
        local nb : word count `blocks'
        local nr = `nb' + 3
        local names `eq0':difference `eq0':explained
        foreach bn of local blocks {
            local names `names' `pfx'`rp'marg_expl:`bn'
        }
        local names `names' `eq0':unexplained
        preserve
        qui drop _all
        qui set obs `nr'
        qui gen double pos = .
        qui gen double val = .
        qui gen double se = .
        qui gen double cil = .
        qui gen double ciu = .
        qui gen byte kind = cond(_n == 1, 1, cond(_n == 2, 2, cond(_n < `nr', 3, 4)))
        local rem
        if `variant' {
            // the remainder: a lighter bar in the explained colour
            qui replace kind = 5 in `=`nr' - 1'
            local rem (rbar zero val pos if kind == 5, `horizontal' barwidth(.6) pstyle(p1bar) fintensity(35))
        }
        local tl
        local len 0
        forvalues r = 1/`nr' {
            local nm : word `r' of `names'
            local c = colnumb(`b', "`nm'")
            if `c' == . {
                di as err "makeplot: `nm' not found in e(b)"
                exit 111
            }
            qui replace val = `b'[1, `c'] in `r'
            if `hasVm' {
                local c = colnumb(`CI', "`nm'")
                if `c' < . {
                    qui replace cil = `CI'[1, `c'] in `r'
                    qui replace ciu = `CI'[2, `c'] in `r'
                    if `CI'[1, `c'] < . & `CI'[2, `c'] < . {
                        qui replace se = sqrt(`V'[`c', `c']) in `r'
                    }
                }
            }
            if `r' == 1 local p 1
            else if `r' < `nr' local p = `r' + .2
            else local p = `r' + .4
            qui replace pos = `p' in `r'
            // a block's row: its text of blocklabels(), else its name (brace blocks: quotes)
            if `r' == 1 local rl "Total diff"
            else if `r' == 2 local rl "Expl."
            else if `r' < `nr' {
                local rl `"`macval(bl_`=`r' - 2')'"'
                local len = max(`len', `bll_`=`r' - 2'')
            }
            else local rl "Unexpl."
            local tl `"`macval(tl)' `p' `"`macval(rl)'"'"'
        }
        qui gen double zero = 0
        qui gen double ll = cil
        qui gen double ul = ciu
        // labels: the value (3 significant digits, or format(); a value that reads 0 without a
        // sign) and its share of the total difference (rounded, never -0), beyond the end of the
        // bar or of its CI spike
        qui gen double anc = cond(se < ., cond(val >= 0, max(val, ul), min(val, ll)), val)
        local p1 = cond("`horizontal'" != "", 3, 12)
        local p2 = cond("`horizontal'" != "", 9, 6)
        qui gen byte lp = cond(val >= 0, `p1', `p2')
        qui gen double vz = cond(abs(val) < `eps', 0, val)
        qui gen double pct = 100 * vz / `D'
        qui gen lab = strtrim(string(vz, "`fmt'"))
        qui replace lab = substr(lab, 2, .) if regexm(lab, "^-[0.,]+$")
        if !`Dz' {
            qui replace lab = lab + " (" + string(cond(abs(pct) < .5, 0, round(pct)), "%9.0f") + "%)" if kind > 1
        }
        local ci
        qui count if se < .
        if r(N) local ci (rcap ll ul pos if se < ., `horizontal' lcolor(foreground) lwidth(thin))
        // room for the labels beyond the bar ends
        summarize anc, meanonly
        local lo = min(0, r(min))
        local hi = max(0, r(max))
        local mlen 0
        forvalues r = 1/`nr' {
            local mlen = max(`mlen', strlen(lab[`r']))
        }
        // (horizontal: long row names narrow the plot, so that a character of a label spans more
        // of the axis: the room grows as in the waterfall in levels, by 70 / (83 - rn), rn the
        // longest row name, at least that of "Total diff")
        local hf 1
        if "`horizontal'" != "" {
            local rn = max(10, `len')
            local hf = max(1, 70 / max(83 - `rn', 20))
        }
        local pad = (`hi' - `lo') * cond("`horizontal'" != "", (.02 + .018 * `mlen') * `hf', .08)
        qui count if val < 0
        if r(N) local lo = `lo' - `pad'
        qui count if val >= 0
        if r(N) local hi = `hi' + `pad'
        if "`horizontal'" != "" {
            local labs (scatter pos anc, msymbol(none) mlabel(lab) mlabvposition(lp) mlabsize(small) ///
                mlabcolor(foreground))
            local axes xline(0, lcolor(gs8) lwidth(thin) lpattern(solid)) xscale(range(`lo' `hi')) ///
                yscale(reverse range(.5 `=`nr' + .9')) ylabel(`macval(tl)', noticks nogrid angle(0))
            local tts ytitle("") xtitle(`macval(ytq)')
        }
        else {
            // as in the waterfall: the row names angled beyond six bars or when a block's name (its
            // text of blocklabels()) is longer than 10 characters, and smaller labels beyond six bars
            local ang = cond(`nr' > 6 | `len' > 10, "angle(30)", "")
            local msz = cond(`nr' > 6, "vsmall", "small")
            local labs (scatter anc pos, msymbol(none) mlabel(lab) mlabvposition(lp) mlabsize(`msz') ///
                mlabcolor(foreground))
            local axes yline(0, lcolor(gs8) lwidth(thin) lpattern(solid)) yscale(range(`lo' `hi')) ///
                ylabel(, angle(0)) xlabel(`macval(tl)', noticks nogrid `ang') xscale(range(.5 `=`nr' + .9'))
            local tts ytitle(`macval(ytq)') xtitle("")
        }
        // the title: the longest of four that fits the width of the graph (Plot_fit)
        Plot_fit ttl `"`macval(options)'"' "Marginal explanatory contribution at each raking stage" ///
            "Marginal explanatory contribution by raking stage" ///
            "Marginal explanatory contribution at each stage" "Marginal explanatory contribution by stage"
        capture noisily twoway ///
            (rbar zero val pos if kind == 1, `horizontal' barwidth(.6) fcolor(none) lcolor(foreground) lpattern(dash)) ///
            (rbar zero val pos if inlist(kind, 2, 3), `horizontal' barwidth(.6) pstyle(p1bar)) ///
            (rbar zero val pos if kind == 4, `horizontal' barwidth(.6) fcolor(gs12) lcolor(gs8) lwidth(thin)) ///
            `rem' `ci' `labs', `macval(axes)' `macval(tts)' ///
            title("`ttl'", span) subtitle(`macval(sub)', span) legend(off) `cinote' `gname' `macval(options)'
    }
    else {
        capture confirm matrix e(bystage)
        if _rc | "`pfx'`rp'" != "" {
            // after a replication prefix only e(b) survives, and an outcome of outcomes() and the
            // reverse order have their own equations: the cumulative values from e(b)
            matrix `S' = J(`ns', 3, .)
            forvalues s = 1/`ns' {
                local bn : word `s' of `blocks'
                matrix `S'[`s', 2] = `b'[1, colnumb(`b', "`pfx'`rp'cum_expl:`bn'")]
                matrix `S'[`s', 3] = `b'[1, colnumb(`b', "`pfx'`rp'cum_unexpl:`bn'")]
            }
        }
        else matrix `S' = e(bystage)
        // Every stage uses its own interval, including asymmetric and zero-width limits.
        local hasV `hasVm'
        preserve
        qui drop _all
        qui set obs `ns'
        qui gen int stage = _n
        qui gen double xe = stage - .2
        qui gen double xu = stage + .2
        qui gen double zero = 0
        qui gen double cum = .
        qui gen double unx = .
        qui gen double ell = .
        qui gen double eul = .
        qui gen double ull = .
        qui gen double uul = .
        local tl
        local len 0
        forvalues s = 1/`ns' {
            qui replace cum = `S'[`s', 2] in `s'
            qui replace unx = `S'[`s', 3] in `s'
            local bn : word `s' of `blocks'
            if `hasV' {
                local pe = colnumb(`CI', "`pfx'`rp'cum_expl:`bn'")
                local pu = colnumb(`CI', "`pfx'`rp'cum_unexpl:`bn'")
                if `pe' < . {
                    qui replace ell = `CI'[1, `pe'] in `s'
                    qui replace eul = `CI'[2, `pe'] in `s'
                }
                if `pu' < . {
                    qui replace ull = `CI'[1, `pu'] in `s'
                    qui replace uul = `CI'[2, `pu'] in `s'
                }
            }
            // the block's text of blocklabels(), else its name (a brace block: quotes)
            if `s' == 1 {
                local tl `"`macval(tl)' `s' `"`macval(bl_`s')'"'"'
            }
            else {
                local tl `"`macval(tl)' `s' `"+ `macval(bl_`s')'"'"'
            }
            local len = max(`len', `bll_`s'')
        }
        // bar labels: share of the total difference (rounded, never -0); none when it is 0 (or
        // displayed as 0, Display_eps)
        foreach v in cum unx {
            qui gen double p_`v' = 100 * `v' / `D'
            qui gen l_`v' = cond(!`Dz', string(cond(abs(p_`v') < .5, 0, round(p_`v')), "%9.0f") + "%", "")
        }
        // the labels hang off the bar end, or off the outer end of the CI cap when there is one;
        // above positive and below negative bars (right and left with horizontal)
        local hor
        local pos1 12
        local pos2 6
        if "`horizontal'" != "" {
            local hor horizontal
            local pos1 3
            local pos2 9
        }
        qui gen double ae = cum
        qui gen double au = unx
        qui gen byte pe = cond(cum >= 0, `pos1', `pos2')
        qui gen byte pu = cond(unx >= 0, `pos1', `pos2')
        local lg 1 "Expl." 2 "Unexpl." 3 "Total diff (`Dtxt')"
        local ci
        if `hasV' {
            qui replace ae = cond(cum >= 0, max(cum, eul), min(cum, ell)) if ell < . & eul < .
            qui replace au = cond(unx >= 0, max(unx, uul), min(unx, ull)) if ull < . & uul < .
            local ci (rcap ell eul xe, `hor' lcolor(foreground) lwidth(thin)) ///
                (rcap ull uul xu, `hor' lcolor(foreground) lwidth(thin))
        }
        // room for the labels beyond the bar ends
        local lo = min(0, `D')
        local hi = max(0, `D')
        foreach v in ae au {
            summarize `v', meanonly
            local lo = min(`lo', r(min))
            local hi = max(`hi', r(max))
        }
        local pad = (`hi' - `lo') * cond("`horizontal'" != "", .1, .08) * (!`Dz')
        qui count if cum < 0 | unx < 0
        if r(N) local lo = `lo' - `pad'
        qui count if cum >= 0 | unx >= 0
        if r(N) local hi = `hi' + `pad'
        // a label that the dashed line of the total would strike through (a share of 97% sits
        // on it) is moved beyond the line: a label runs from its anchor away from its bar, and
        // when the line lies within the label's extent (5% of the axis in the vertical view; in the
        // horizontal view the label's length, as the margin reserves it) the anchor moves to the line
        if !`Dz' {
            local rg = `hi' - `lo'
            foreach v in cum unx {
                local a = cond("`v'" == "cum", "ae", "au")
                local ext = cond("`horizontal'" != "", "(.02 + .018 * strlen(l_`v'))", ".05")
                qui replace `a' = `D' if cond(`v' >= 0, 1, -1) * (`D' - `a') > 0 & ///
                    cond(`v' >= 0, 1, -1) * (`D' - `a') < `ext' * `rg'
            }
        }
        local tot (function y = `D', range(.5 `=`ns' + .5') n(2) `hor' lpattern(dash) lcolor(foreground))
        if "`horizontal'" != "" {
            local labs (scatter xe ae, msymbol(none) mlabel(l_cum) mlabvposition(pe) mlabsize(small) mlabcolor(foreground)) ///
                (scatter xu au, msymbol(none) mlabel(l_unx) mlabvposition(pu) mlabsize(small) mlabcolor(gs6))
            local axes xline(0, lcolor(gs8) lwidth(thin) lpattern(solid)) xscale(range(`lo' `hi')) ///
                ylabel(`macval(tl)', noticks nogrid angle(0)) yscale(reverse range(.5 `=`ns' + .5'))
            local tts xtitle(`macval(ytq)') ytitle("")
        }
        else {
            local labs (scatter ae xe, msymbol(none) mlabel(l_cum) mlabvposition(pe) mlabsize(small) mlabcolor(foreground)) ///
                (scatter au xu, msymbol(none) mlabel(l_unx) mlabvposition(pu) mlabsize(small) mlabcolor(gs6))
            // the block names are angled beyond four stages or when a name is long
            local ang = cond(`ns' > 4 | `len' > 10, "angle(30)", "")
            local axes yline(0, lcolor(gs8) lwidth(thin) lpattern(solid)) yscale(range(`lo' `hi')) ///
                ylabel(, angle(0)) xlabel(`macval(tl)', noticks nogrid `ang') xscale(range(.5 `=`ns' + .5'))
            local tts xtitle("") ytitle(`macval(ytq)')
        }
        capture noisily twoway (rbar zero cum xe, `hor' barwidth(.34)) ///
               (rbar zero unx xu, `hor' barwidth(.34) fcolor(gs12) lcolor(gs8) lwidth(thin)) ///
               `tot' `ci' `labs', `macval(axes)' `macval(tts)' ///
            title("Explained and unexplained difference by raking stage", span) ///
            subtitle(`macval(sub)', span) legend(order(`lg') rows(1) position(6) size(small)) ///
            `cinote' `gname' `macval(options)'
    }
    if _rc {
        local prc = _rc
        di as err "{p 0 4 2}makeplot: twoway returned error `prc'; the results are posted; redraw the graph with seldecomp, makeplot(){p_end}"
        exit `prc'
    }
    restore
end


* ============================================================================
* Mata: jackknife acceleration for BCa
* ============================================================================
version 16
mata:
real rowvector seldecomp_accel(real matrix JK)
{
    real rowvector a
    real colvector v, d
    real scalar j, den
    a = J(1, cols(JK), .)
    for (j = 1; j <= cols(JK); j++) {
        v = JK[., j]
        v = select(v, v :< .)
        if (rows(v) < 2) continue
        d = mean(v) :- v
        den = 6 * (sum(d :^ 2))^1.5
        a[j] = (den > 0 ? sum(d :^ 3) / den : 0)
    }
    return(a)
}

// seldecomp_ixlab: for covariatedetail (Display_covdetail), the label of a cell of an
// interaction of the variables v, "v1=a (A) v2=b ...", shortened to "A # b ...": each
// component's value label, else its value (a string value or a c. bin as it stands); the
// label as it stands when it does not have that form
string scalar seldecomp_ixlab(string scalar lab, string rowvector v)
{
    real scalar j, p
    string scalar rest, out, piece
    string rowvector pc
    if (substr(lab, 1, strlen(v[1]) + 1) != v[1] + "=") return(lab)
    rest = substr(lab, strlen(v[1]) + 2, .)
    pc = J(1, cols(v), "")
    for (j = 2; j <= cols(v); j++) {
        p = strpos(rest, " " + v[j] + "=")
        if (p == 0) return(lab)
        pc[j - 1] = substr(rest, 1, p - 1)
        rest = substr(rest, p + strlen(v[j]) + 2, .)
    }
    pc[cols(v)] = rest
    out = ""
    for (j = 1; j <= cols(v); j++) {
        piece = pc[j]
        if (regexm(piece, "^-?[0-9]+ [(](.*)[)]$")) piece = regexs(1)
        out = out + (j > 1 ? " # " : "") + piece
    }
    return(out)
}

// varlabels: the switch of the display under way (an external, set by Display, Display_estat and
// Make_plot at their start), and the variable labels shown for variable names when it is on
void seldecomp_vlset(real scalar on)
{
    pointer(real scalar) scalar p
    if ((p = findexternal("_seldecomp_vl")) == NULL) p = crexternal("_seldecomp_vl")
    *p = on
}

real scalar seldecomp_vlon()
{
    pointer(real scalar) scalar p
    p = findexternal("_seldecomp_vl")
    return(p == NULL ? 0 : *p)
}

// seldecomp_vl1: the variable label of the variable named v (exactly; not abbreviated), else v
string scalar seldecomp_vl1(string scalar v)
{
    real scalar k
    string scalar l
    if (!st_isname(v)) return(v)
    k = _st_varindex(v, 0)
    if (k >= .) return(v)
    l = st_varlabel(k)
    return(strtrim(l) == "" ? v : l)
}

// seldecomp_vlab: with varlabels on, a list of variable names (or terms a#b) with each name
// replaced by its variable label, the parts of a term joined by " # " and the terms by ", ";
// otherwise the list as it is. The text as typed: smcl escapes its braces for the output
string scalar seldecomp_vlab(string scalar s)
{
    string rowvector t, q
    string scalar out, one
    real scalar i, j
    if (!seldecomp_vlon()) return(s)
    t = tokens(s)
    out = ""
    for (i = 1; i <= cols(t); i++) {
        q = tokens(t[i], "#")
        one = ""
        for (j = 1; j <= cols(q); j++) {
            if (q[j] == "#") one = one + " # "
            else one = one + seldecomp_vl1(q[j])
        }
        out = out + (i > 1 ? ", " : "") + one
    }
    return(out)
}

// seldecomp_smcl: literal display text, safe inside quoted strings and against a later
// macro expansion. Escape user braces before introducing the character codes.
string scalar seldecomp_smcl(string scalar s)
{
    s = subinstr(subinstr(s, "{", char(1)), "}", char(2))
    s = subinstr(subinstr(s, char(1), "{c -(}"), char(2), "{c )-}")
    s = subinstr(subinstr(s, char(36), "{c S|}"), char(96), "{c 96}")
    return(subinstr(s, char(34), "{c 34}"))
}

// seldecomp_vltrim: a text shortened to at most n characters, ending with "~" when it is cut
string scalar seldecomp_vltrim(string scalar s, real scalar n)
{
    return(ustrlen(s) <= n ? s : usubstr(s, 1, n - 1) + "~")
}

// seldecomp_vlnote: with varlabels on, a note with the variable names of its level listings
// (a name followed by " = " or "#", or after "#") replaced by their labels. The caller
// escapes the completed note once, so braces in a label remain literal braces.
string scalar seldecomp_vlnote(string scalar s)
{
    string scalar out, w, rest
    real scalar p
    if (!seldecomp_vlon()) return(s)
    out = ""
    rest = s
    while (ustrregexm(rest, "[A-Za-z_][A-Za-z0-9_]*")) {
        w = ustrregexs(0)
        p = strpos(rest, w)
        out = out + substr(rest, 1, p - 1)
        rest = substr(rest, p + strlen(w), .)
        if (substr(rest, 1, 3) == " = " | substr(rest, 1, 1) == "#" | substr(out, -1, 1) == "#") {
            out = out + seldecomp_vl1(w)
        }
        else out = out + w
    }
    return(out + rest)
}

// seldecomp_ptxt: for makeplot, a text (a label, a blocklabels() text) as twoway shows it as
// typed: {, }, $ and ` written as {c -(}, {c )-}, {c S|} and {c 96} (twoway expands its text
// options again and reads SMCL tags in them), and every double quote as {c 34} (leaving a
// closing quote unescaped can make twoway read subtitle options as additional text lines)
string scalar seldecomp_ptxt(string scalar s)
{
    return(seldecomp_smcl(s))
}

// seldecomp_pfit: a text (as typed) shortened to at most n characters, ending with "..." when
// it is cut: after the last whole word that fits when that keeps more than half of the room,
// else within the word; spaces, commas, semicolons and colons before the "..." dropped
string scalar seldecomp_pfit(string scalar s, real scalar n)
{
    real scalar p
    string scalar t
    if (ustrlen(s) <= n) return(s)
    t = usubstr(s, 1, n - 3)
    if (usubstr(s, n - 2, 1) != " ") {
        p = ustrrpos(t, " ")
        if (p > (n - 3) / 2) t = usubstr(t, 1, p - 1)
    }
    return(ustrregexra(t, "[ ,;:]+$", "") + "...")
}

// seldecomp_pvt: for makeplot, the title of the value axis as a title option takes it: "a: b" (a,
// the variable labels, as typed) on one line up to 40 characters; otherwise the two lines "a:"
// (shortened to 40 characters with "...") and "b". A vertical title is centred on the plot, which
// a two-line title or subtitle and angled names shorten: under s1mono it holds about 50
// characters, and runs into the title beyond about 40
string scalar seldecomp_pvt(string scalar a, string scalar b)
{
    string scalar lq, rq
    lq = char(96) + char(34)
    rq = char(34) + char(39)
    if (ustrlen(a + ": " + b) <= 40) return(lq + seldecomp_ptxt(a + ": " + b) + rq)
    return(lq + seldecomp_ptxt(seldecomp_pfit(a, 40) + ":") + rq + " " + lq + b + rq)
}

// seldecomp_pwrap: for makeplot's horizontal view, the text of a row named by a group (as typed)
// as twoway takes it, followed by tail (", raked" or ""): one line when the group is at most n
// characters; otherwise two lines split at a space, both at most n: after a comma when that fits
// ("Married," "spouse present"), else where the longer line is shortest; when no split fits, the
// words that fit in n on the first line and the rest on the second, shortened to n characters
// with "..." (seldecomp_pfit; a first word longer than n shortened likewise, and a group without
// spaces is one shortened line); tail after the last line; each line escaped (seldecomp_ptxt), two
// lines as `"a"' `"b"'
string scalar seldecomp_pwrap(string scalar s, real scalar n, string scalar tail)
{
    real scalar i, m, k, c, best, bk, bc
    string scalar a, b, lq, rq
    if (ustrlen(s) <= n) return(seldecomp_ptxt(s + tail))
    s = ustrtrim(s)
    m = ustrlen(s)
    if (m <= n) return(seldecomp_ptxt(s + tail))
    if (!ustrpos(s, " ")) return(seldecomp_ptxt(seldecomp_pfit(s, n) + tail))
    best = 0
    bk = .
    bc = 0
    for (i = 2; i < m; i++) {
        if (usubstr(s, i, 1) != " ") continue
        a = ustrrtrim(usubstr(s, 1, i - 1))
        b = ustrltrim(usubstr(s, i + 1, .))
        k = max((ustrlen(a), ustrlen(b)))
        c = (usubstr(a, ustrlen(a), 1) == ",")
        if (k <= n & (c > bc | (c == bc & k < bk))) {
            best = i
            bk = k
            bc = c
        }
    }
    if (best == 0) {
        for (i = min((n + 1, m)); i >= 2; i--) {
            if (usubstr(s, i, 1) == " ") break
        }
        if (i >= 2) best = i
    }
    if (best) {
        a = ustrrtrim(usubstr(s, 1, best - 1))
        b = ustrltrim(usubstr(s, best + 1, .))
    }
    else {
        i = ustrpos(s, " ")
        a = seldecomp_pfit(usubstr(s, 1, i - 1), n)
        b = ustrltrim(usubstr(s, i + 1, .))
    }
    b = seldecomp_pfit(b, n)
    lq = char(96) + char(34)
    rq = char(34) + char(39)
    return(lq + seldecomp_ptxt(a) + rq + " " + lq + seldecomp_ptxt(b + tail) + rq)
}

// seldecomp_group: a term's levels, as egen group() makes them: the distinct rows of the
// components (compv) on tousev, sorted in ascending order of the components (lexicographic),
// are numbered 1, 2, ...; each observation's number goes to codev (an existing variable),
// the distinct rows to the matrix uname
void seldecomp_group(string scalar compv, string scalar tousev, string scalar codev,
    string scalar uname)
{
    real matrix X, Xs
    real colvector o, nw, g
    real scalar n
    X = st_data(., tokens(compv), tousev)
    n = rows(X)
    o = order(X, 1..cols(X))
    Xs = X[o, .]
    nw = 1
    if (n > 1) nw = 1 \ (rowsum(Xs[|2, 1 \ n, .|] :!= Xs[|1, 1 \ n - 1, .|]) :> 0)
    g = J(n, 1, .)
    g[o] = runningsum(nw)
    st_store(., codev, tousev, g)
    st_matrix(uname, Xs[selectindex(nw), .])
}

// seldecomp_trim: the trimming rule: each weight w bounded to its upper bound, the smaller of
// w0 * hirel and hiabs, and to its lower bound, the larger of w0 * lorel and loabs (tb =
// hirel, lorel, hiabs, loabs; -1 an absent bound), so that every bound binds and the stricter
// of two bounds on one side wins; w0 the base weights. The bounds never contradict here:
// seldecomp_rake refuses bounds whose lower exceeds the upper for an observation (tclash).
// Nothing is renormalised. This clipping rule is not a bounded-KL optimization:
// intermediate clipping can change the final solution even when no final bound binds.
void seldecomp_trim(real colvector w, real colvector w0, real rowvector tb)
{
    real colvector lo, hi
    hi = J(rows(w), 1, .)
    lo = J(rows(w), 1, .)
    if (tb[1] != -1) hi = w0 * tb[1]
    if (tb[3] != -1) hi = rowmin((hi, J(rows(w), 1, tb[3])))
    if (tb[2] != -1) lo = w0 * tb[2]
    if (tb[4] != -1) lo = rowmax((lo, J(rows(w), 1, tb[4])))
    w = rowmin((w, hi))
    w = rowmax((w, lo))
}

// seldecomp_tclash: the observations whose bounds contradict, a lower bound above an upper
// one: (the number with loabs above w0 * hirel, the number with w0 * lorel above hiabs);
// trimhiabs() above trimloabs() and trimhirel() above trimlorel() are checked when parsed
real rowvector seldecomp_tclash(real colvector w0, real rowvector tb)
{
    real scalar c1, c2
    c1 = c2 = 0
    if (tb[4] != -1 & tb[1] != -1) c1 = sum(tb[4] :> w0 * tb[1])
    if (tb[2] != -1 & tb[3] != -1) c2 = sum(w0 * tb[2] :> tb[3])
    return((c1, c2))
}

// seldecomp_rake: raking by iterative proportional fitting as Kolenikov's ipfraking 1.3.74 does it
// (its PropAdjust, TrimWeights, CheckConvergence, CheckResults and ReportTrims: the reference
// algorithm the comments below refer to). A cycle visits the control matrices in
// order, and for each level multiplies the weights of its observations by target / current
// total (the totals of a margin are taken before its adjustment; the levels partition the
// sample). With a trimming bound (tfreq > 0) the weights are bounded (seldecomp_trim) after
// each margin (tfreq 3, often), after each cycle (2, sometimes) or once after the last cycle
// (1, once). After each cycle: converged when no weight changed by more than tol (relative,
// from the cycle's start to its end, trimming included); stop then, or after iterate cycles.
// Without clipping, ordinary IPF has the minimum-KL solution and its final fit is
// invariant to margin order up to tolerance. Clipping does not preserve those
// guarantees; convergence measures weight stabilization, not bounded optimality.
// Then, per margin, the
// largest reldif() between achieved and target totals; the margins above ctrltol are counted
// in badcontrols, and the worst level is the last one within epsfloat of the running maximum
// (the reference's rule); and per bound the weights at it within a relative 1e-6 (the
// reference's counts, on the final weights whatever tfreq). Arguments: wsv the base weights, tousev the
// raking sample, ctotal the control matrices (colnames = levels, rownames = the level
// variable), tb the bounds (hirel, lorel, hiabs, loabs; -1 absent), genv an existing double
// variable (the only data written, at the end), resname a 1 x 17 result (converged,
// badcontrols, maxreldif, maxctrl, worst margin, worst level, cycles, the weights at
// trimhiabs, trimloabs, trimhirel, trimlorel, the two counts of contradicting bounds, and the
// levels of target 0 missed: their observations kept weight, missed level targets,
// checked level targets, distinct weights at any final bound); local worstvar.
// Final-bound counts do not count weights clipped only in intermediate cycles.
// Optional missname records one row per affected
// margin: margin index, worst level, target, achieved total, discrepancy (scaled units).
void seldecomp_rake(string scalar wsv, string scalar tousev, string scalar ctotal,
    real scalar iterate, real scalar tol, real scalar ctrltol,
    real rowvector tb, real scalar tfreq, string scalar genv, real scalar trace,
    string scalar resname, | string scalar missname)
{
    string rowvector mats, cv
    string matrix cs
    real colvector w, w0, wprev, d, code, pos, f, atbound, hit
    real rowvector T, tg, lv, nl, tot, rd, ntr, tcl, achieved
    real matrix misses
    pointer(real colvector) rowvector ix, li
    pointer() scalar pb
    real scalar K, k, l, p, nlev, it, ncyc, conv, obj, badc, bad0, mrdmax, mr, wk, wl, epsf
    real scalar brk, badtargets, worst

    if (args() < 12) missname = ""
    if (missname != "") misses = J(0, 5, .)
    // test hook: with the external seldecomp_brk = n, the call after n
    // further calls stops after its first cycle with return code 1, as a Break does
    brk = 0
    if ((pb = findexternal("seldecomp_brk")) != NULL) {
        if (*pb <= 0) {
            brk = 1
            rmexternal("seldecomp_brk")
        }
        else *pb = *pb - 1
    }
    mats = tokens(ctotal)
    K = cols(mats)
    nlev = 0
    for (k = 1; k <= K; k++) nlev = nlev + cols(st_matrix(mats[k]))
    // one index vector per level of every margin, in the order of the matrices and columns,
    // and per margin the column of each observation's level (li)
    ix = J(1, nlev, NULL)
    li = J(1, K, NULL)
    tg = lv = J(1, nlev, .)
    nl = J(1, K, .)
    cv = J(1, K, "")
    p = 0
    for (k = 1; k <= K; k++) {
        T = st_matrix(mats[k])
        cs = st_matrixcolstripe(mats[k])
        cv[k] = st_matrixrowstripe(mats[k])[1, 2]
        code = st_data(., cv[k], tousev)
        nl[k] = cols(T)
        pos = J(max((max(code), max(strtoreal(cs[., 2])))), 1, .)
        for (l = 1; l <= nl[k]; l++) {
            p++
            lv[p] = strtoreal(cs[l, 2])
            tg[p] = T[l]
            ix[p] = &(selectindex(code :== lv[p]))
            pos[lv[p]] = l
        }
        li[k] = &(pos[code])
    }
    w = st_data(., wsv, tousev)
    if (tfreq) w0 = w
    // bounds that contradict for an observation: no raking, their counts in the result
    tcl = (tfreq ? seldecomp_tclash(w0, tb) : (0, 0))
    if (tcl[1] + tcl[2] > 0) {
        st_matrix(resname, (0, ., ., ., ., ., 0, J(1, 4, 0), tcl, 0, ., nlev, 0))
        st_local("worstvar", "")
        return
    }
    if (trace) {
        printf("{txt}Raking: %g margins, %g levels, %g observations\n",
            K, nlev, rows(w))
        if (tfreq) {
            printf("{txt}  trimming after each %s\n",
                (tfreq == 3 ? "margin" : (tfreq == 2 ? "cycle" : "raking (once)")))
        }
    }
    tot = J(1, nlev, .)
    conv = 0
    obj = .
    for (it = 1; it <= iterate; it++) {
        wprev = w
        p = 0
        for (k = 1; k <= K; k++) {
            // the factor of each level (1 for a level whose weights are all 0, which is
            // left alone), then every weight times its level's factor
            f = J(nl[k], 1, 1)
            for (l = 1; l <= nl[k]; l++) {
                tot[p + l] = quadsum(w[*ix[p + l]])
                if (tot[p + l] != 0) f[l] = tg[p + l] / tot[p + l]
            }
            w = w :* f[*li[k], 1]
            if (tfreq == 3) seldecomp_trim(w, w0, tb)
            p = p + nl[k]
        }
        if (tfreq == 2) seldecomp_trim(w, w0, tb)
        if (brk) exit(1)
        d = abs(w - wprev) :/ wprev
        conv = !any((d :> tol) :& (d :< .))
        d = select(d, d :< .)
        obj = (rows(d) ? max(d) : .)
        if (trace) {
            printf("{txt}  cycle %5.0f: largest relative change of a weight {res}%10.3e\n", it, obj)
        }
        if (conv) break
    }
    ncyc = (it > iterate ? iterate : it)
    if (tfreq == 1) seldecomp_trim(w, w0, tb)
    badc = bad0 = badtargets = 0
    mrdmax = 0
    wk = wl = .
    epsf = c("epsfloat")
    p = 0
    for (k = 1; k <= K; k++) {
        rd = J(1, nl[k], .)
        achieved = J(1, nl[k], .)
        for (l = 1; l <= nl[k]; l++) {
            achieved[l] = quadsum(w[*ix[p + l]])
            rd[l] = reldif(achieved[l], tg[p + l])
        }
        mr = max(rd)
        badtargets = badtargets + sum(rd :> ctrltol)
        if (mr > ctrltol) {
            badc++
            if (missname != "") {
                // Exact maximum within this term, last level in an exact tie. Unlike the
                // historical global worst-label tie rule below, this always names a miss.
                worst = 1
                for (l = 1; l <= nl[k]; l++) {
                    if (rd[l] == mr) worst = l
                }
                misses = misses \ (k, lv[p + worst], tg[p + worst], achieved[worst], rd[worst])
            }
        }
        // a level of target 0 (absent from the Reference Group) whose observations keep weight
        for (l = 1; l <= nl[k]; l++) {
            if (tg[p + l] == 0 & rd[l] > ctrltol) bad0++
        }
        mrdmax = max((mrdmax, mr))
        for (l = 1; l <= nl[k]; l++) {
            if (reldif(rd[l], mrdmax) < epsf) {
                wk = k
                wl = lv[p + l]
            }
        }
        p = p + nl[k]
    }
    // Per-bound trace counts retain the reference's tolerance. The public count
    // takes the union: each weight at one or more final bounds is counted once.
    ntr = J(1, 4, 0)
    atbound = J(rows(w), 1, 0)
    if (tb[3] != -1) {
        hit = (w :/ tb[3]) :> 1 - 1e-6
        ntr[1] = sum(hit)
        atbound = atbound :| hit
    }
    if (tb[4] != -1) {
        hit = (w :/ tb[4]) :< 1 + 1e-6
        ntr[2] = sum(hit)
        atbound = atbound :| hit
    }
    if (tb[1] != -1) {
        hit = (w :/ w0) :> (1 - 1e-6) * tb[1]
        ntr[3] = sum(hit)
        atbound = atbound :| hit
    }
    if (tb[2] != -1) {
        hit = (w :/ w0) :< (1 + 1e-6) * tb[2]
        ntr[4] = sum(hit)
        atbound = atbound :| hit
    }
    if (trace) {
        printf("{txt}  %s after %g cycles; margins missed by more than ctrltolerance(): {res}%g\n",
            (conv ? "converged" : "not converged"), ncyc, badc)
        printf("{txt}  largest relative difference from a control total {res}%10.3e\n", mrdmax)
        if (tfreq) {
            printf("{txt}  weights at a trimming bound: hiabs {res}%g{txt}, loabs {res}%g{txt}, hirel {res}%g{txt}, lorel {res}%g\n",
                ntr[1], ntr[2], ntr[3], ntr[4])
        }
    }
    st_store(., genv, tousev, w)
    st_matrix(resname, (conv, badc, obj, mrdmax, wk, wl, ncyc, ntr, tcl, bad0, badtargets, nlev, sum(atbound)))
    if (missname != "" & badc > 0) st_matrix(missname, misses)
    st_local("worstvar", (wk < . ? cv[wk] : ""))
}

// seldecomp_msd: mean, min, max, sd (n - 1), cv (sd / mean when the mean is positive)
real rowvector seldecomp_msd(real colvector x)
{
    real scalar m, sd
    m = quadsum(x) / rows(x)
    sd = (rows(x) > 1 ? sqrt(quadsum((x :- m) :^ 2) / (rows(x) - 1)) : .)
    return((m, min(x), max(x), sd, (m > 0 ? sd / m : .)))
}

// seldecomp_diag: Engine's weight diagnostics on the raking sample (rssv): mean, minimum,
// maximum, standard deviation (n - 1 denominator, as summarize) and coefficient of variation
// of the raked weights and of the factors raked / base, then the Kish effective sizes
// Wrss^2 / sum(w^2) and sum(rw)^2 / sum(rw^2); a 1 x 12 row in rname
void seldecomp_diag(string scalar rwv, string scalar wv, string scalar rssv, real scalar Wrss,
    string scalar rname)
{
    real colvector rw, w, f
    real rowvector a, b
    real scalar rwsq
    rw = st_data(., rwv, rssv)
    w = st_data(., wv, rssv)
    f = rw :/ w
    a = seldecomp_msd(rw)
    b = seldecomp_msd(f)
    rwsq = quadsum(rw :^ 2)
    st_matrix(rname, (a, b, Wrss^2 / quadsum(w :^ 2), (rwsq > 0 ? quadsum(rw)^2 / rwsq : .)))
}

// seldecomp_bw: Draw_weights' replicate weights on tousev, bw = w * g * den / num, num the
// total of w * g (denv: the total of w)
void seldecomp_bw(string scalar wv, string scalar gv, string scalar denv, string scalar tousev,
    string scalar bwv)
{
    real colvector wg
    wg = st_data(., wv, tousev) :* st_data(., gv, tousev)
    st_store(., bwv, tousev, wg :* st_data(., denv, tousev) :/ quadsum(wg))
}

// seldecomp_g2shares: Engine's Reference Group shares, once per call. P (terms x levels) holds the
// base-weight total of each level over Reference Group (selv: touse, Reference Group, positive weight)
// divided by W2, and 0 for a level without Reference Group observations
void seldecomp_g2shares(string scalar tvlist, string scalar nlevlist, string scalar wv,
    string scalar selv, real scalar W2, string scalar pname)
{
    string rowvector tv
    real rowvector nlev
    real colvector w, x, m
    real matrix P
    real scalar t, l
    tv = tokens(tvlist)
    nlev = strtoreal(tokens(nlevlist))
    w = st_data(., wv, selv)
    P = J(cols(tv), max(nlev), 0)
    for (t = 1; t <= cols(tv); t++) {
        x = st_data(., tv[t], selv)
        for (l = 1; l <= nlev[t]; l++) {
            m = (x :== l)
            if (any(m)) P[t, l] = quadsum(select(w, m)) / W2
        }
    }
    st_matrix(pname, P)
}

// seldecomp_ctrlmat: control matrix name (colnames = the levels in lev, coleq = the ones
// variable, rowname = the level variable) with targets W1 x share / (sum of the shares),
// the shares summed in level order. The coleq is cosmetic: seldecomp_rake reads the column
// names (the level codes) and the row name (the level variable) only
void seldecomp_ctrlmat(string scalar name, real rowvector sh, real rowvector lev,
    real scalar W1, string scalar onev, string scalar tvn)
{
    real rowvector M
    real scalar i, psum
    psum = 0
    for (i = 1; i <= cols(sh); i++) psum = psum + sh[i]
    M = J(1, cols(sh), 0)
    for (i = 1; i <= cols(sh); i++) M[i] = W1 * sh[i] / psum
    st_matrix(name, M)
    st_matrixcolstripe(name, (J(cols(lev), 1, onev), strofreal(lev')))
    st_matrixrowstripe(name, ("", tvn))
}

// seldecomp_stage1: one stage of Engine, the terms raked to the Reference Group (st). On the Focal Group
// observations with positive weight (rsv): the zero flag (zerov: 1 where a level of a stage
// term has no Reference Group weight) and the raking sample (rssv = rsv & !zerov); per stage term,
// the levels present in the raking sample, the thin target levels (fewer than mincell observations),
// the support failures (a level with Reference Group weight but none left in the raking sample: rows
// term, level, whether Focal Group has it at all, in matrix fname) and, while no failure has
// been found, the control matrix (seldecomp_ctrlmat; the names in mnames, one per term).
// incl: the observations with a level absent from the Reference Group stay in the raking sample
// (a lower trimming bound holds them), their levels' targets 0; the zero flag is then 0.
// Locals: zerolist (the cells of those observations), Nrake, Nin (the others), nsup (their
// number), nzero (the number set aside: 0 with incl), nthin, thinlist, nfail, mats.
void seldecomp_stage1(string scalar tvlist, string scalar nlevlist, string scalar stlist,
    string scalar rsv, string scalar zerov, string scalar rssv, string scalar pname,
    real scalar W1, real scalar mincell, string scalar onev, string scalar mnlist,
    string scalar fname, real scalar incl)
{
    string rowvector tv, mn
    string scalar zl, thl, mats
    real rowvector nlev, st, pres
    real colvector z, keep, x, m
    real matrix P, X, F
    real scalar j, t, l, c, nthin, nfail, nz
    tv = tokens(tvlist)
    mn = tokens(mnlist)
    nlev = strtoreal(tokens(nlevlist))
    st = strtoreal(tokens(stlist))
    P = st_matrix(pname)
    X = st_data(., tv[st], rsv)
    z = J(rows(X), 1, 0)
    zl = ""
    for (j = 1; j <= cols(st); j++) {
        t = st[j]
        for (l = 1; l <= nlev[t]; l++) {
            if (P[t, l] == 0) {
                m = (X[., j] :== l)
                c = sum(m)
                if (c > 0) {
                    z = z :| m
                    zl = zl + " " + strofreal(t) + ":" + strofreal(l) + ":" + strofreal(c, "%18.0g")
                }
            }
        }
    }
    nz = sum(z)
    keep = (incl ? J(rows(X), 1, 1) : !z)
    st_store(., zerov, J(st_nobs(), 1, 0))
    st_store(., zerov, rsv, (incl ? J(rows(X), 1, 0) : z))
    st_store(., rssv, J(st_nobs(), 1, 0))
    st_store(., rssv, rsv, keep)
    st_local("zerolist", strtrim(zl))
    st_local("Nrake", strofreal(sum(keep), "%18.0g"))
    st_local("Nin", strofreal(rows(X) - nz, "%18.0g"))
    st_local("nsup", strofreal(nz, "%18.0g"))
    st_local("nzero", strofreal(incl ? 0 : nz, "%18.0g"))
    nthin = nfail = 0
    thl = mats = ""
    F = J(0, 3, .)
    if (rows(X) - nz > 0) {
        for (j = 1; j <= cols(st); j++) {
            t = st[j]
            x = X[., j]
            pres = J(1, 0, .)
            for (l = 1; l <= nlev[t]; l++) {
                m = (x :== l)
                c = sum(keep :& m)
                if (c > 0) {
                    pres = pres, l
                    // a level absent from the Reference Group (incl) has its own note
                    if (c < mincell & P[t, l] > 0) {
                        nthin++
                        thl = thl + " " + strofreal(t) + ":" + strofreal(l) + ":" + strofreal(c, "%18.0g")
                    }
                }
                else if (P[t, l] > 0) {
                    nfail++
                    F = F \ (t, l, any(m))
                }
            }
            if (nfail == 0 & cols(pres) > 1) {
                seldecomp_ctrlmat(mn[t], P[t, pres], pres, W1, onev, tv[t])
                mats = mats + " " + mn[t]
            }
        }
    }
    st_local("nthin", strofreal(nthin))
    st_local("thinlist", strtrim(thl))
    st_local("nfail", strofreal(nfail))
    st_local("mats", strtrim(mats))
    if (nfail > 0) st_matrix(fname, F)
}

// seldecomp_stage2: with fixothermargins, the terms not yet added (ot) held at the Focal Group's own
// shares on the raking sample (rssv; total base weight Wrss): their control matrices are
// appended to the local mats
void seldecomp_stage2(string scalar tvlist, string scalar nlevlist, string scalar otlist,
    string scalar wv, string scalar rssv, real scalar Wrss, real scalar W1,
    string scalar onev, string scalar mnlist)
{
    string rowvector tv, mn
    string scalar mats
    real rowvector nlev, ot, pres, sh
    real colvector w, x, m
    real matrix X
    real scalar j, t, l
    tv = tokens(tvlist)
    mn = tokens(mnlist)
    nlev = strtoreal(tokens(nlevlist))
    ot = strtoreal(tokens(otlist))
    w = st_data(., wv, rssv)
    X = st_data(., tv[ot], rssv)
    mats = st_local("mats")
    for (j = 1; j <= cols(ot); j++) {
        t = ot[j]
        x = X[., j]
        pres = sh = J(1, 0, .)
        for (l = 1; l <= nlev[t]; l++) {
            m = (x :== l)
            if (any(m)) {
                pres = pres, l
                sh = sh, quadsum(select(w, m)) / Wrss
            }
        }
        if (cols(pres) > 1) {
            seldecomp_ctrlmat(mn[t], sh, pres, W1, onev, tv[t])
            mats = mats + " " + mn[t]
        }
    }
    st_local("mats", strtrim(mats))
}

// seldecomp_screen: insupport. On the estimation sample (tousev; every weight positive), the
// Focal Group observations (g1v == 1) with a level of any term (codes 1..nlev in tvlist) that has
// no Reference Group observation are flagged in outv (an existing variable). Locals: nexcl (their
// number), inslist ("t:l:count" for each such level found in the Focal Group), insfail ("t:l" for
// each level found in both groups whose Focal Group observations are all flagged)
void seldecomp_screen(string scalar tvlist, string scalar nlevlist, string scalar tousev,
    string scalar g1v, string scalar outv)
{
    string rowvector tv
    string scalar zl, fl
    real rowvector nlev
    real matrix X, IN2
    real colvector g, out, x, m
    real scalar t, l, c
    tv = tokens(tvlist)
    nlev = strtoreal(tokens(nlevlist))
    X = st_data(., tv, tousev)
    g = st_data(., g1v, tousev)
    out = J(rows(X), 1, 0)
    IN2 = J(max(nlev), cols(tv), 0)
    zl = fl = ""
    for (t = 1; t <= cols(tv); t++) {
        x = X[., t]
        for (l = 1; l <= nlev[t]; l++) {
            IN2[l, t] = any((x :== l) :& (g :== 0))
            if (IN2[l, t]) continue
            m = (x :== l) :& (g :== 1)
            c = sum(m)
            if (c > 0) {
                out = out :| m
                zl = zl + " " + strofreal(t) + ":" + strofreal(l) + ":" + strofreal(c, "%18.0g")
            }
        }
    }
    for (t = 1; t <= cols(tv); t++) {
        x = X[., t]
        for (l = 1; l <= nlev[t]; l++) {
            if (!IN2[l, t]) continue
            m = (x :== l) :& (g :== 1)
            if (any(m) & !any(m :& !out)) fl = fl + " " + strofreal(t) + ":" + strofreal(l)
        }
    }
    st_store(., outv, tousev, out)
    st_local("nexcl", strofreal(sum(out), "%18.0g"))
    st_local("inslist", strtrim(zl))
    st_local("insfail", strtrim(fl))
}

// seldecomp_commaby: Comma_check. In the command line s, the text before the first comma
// outside parentheses and double quotes; when it holds the word by( outside parentheses (by()
// typed before the comma), locals cbfound = 1, cbpre (the text before by) and cbpost (by()
// and what follows it, then the text after the comma); cbfound = 0 otherwise
void seldecomp_commaby(string scalar s)
{
    real scalar i, n, d, q, c, b
    string scalar ch, pre
    n = strlen(s)
    d = q = 0
    c = n + 1
    b = 0
    for (i = 1; i <= n; i++) {
        ch = substr(s, i, 1)
        if (ch == char(34)) q = !q
        else if (!q & ch == "(") {
            if (d == 0 & !b & i >= 3) {
                if (substr(s, i - 2, 2) == "by") {
                    if (i == 3) b = 1
                    else if (substr(s, i - 3, 1) == " ") b = i - 2
                }
            }
            d++
        }
        else if (!q & ch == ")") d--
        else if (!q & d == 0 & ch == ",") {
            c = i
            break
        }
    }
    st_local("cbfound", strofreal(b > 0))
    if (!b) return
    pre = substr(s, 1, c - 1)
    st_local("cbpre", substr(pre, 1, b - 1))
    st_local("cbpost", substr(pre, b, .) + (c <= n ? " " + substr(s, c + 1, .) : ""))
}

// seldecomp_prefail: the pre-flight check of an svy prefix's replicates,
// in its full-sample call. On the estimation sample (tousev; codes 1..nlev in tvlist,
// Focal Group where g1v == 1), a replicate's sample is the observations with a positive weight in
// it: replicate-weight variable r of rwlist, or, with psuv, all but PSU j (the delete-one-PSU
// jackknife; only the PSUs that hold all the Focal Group or all the Reference Group observations of a
// level can change anything). As Engine's last stage (every term raked to the Reference Group), a Focal Group
// observation with a level that has no Reference Group observation in the replicate is set aside, and
// a level with Reference Group observations but no Focal Group observation left fails the replicate
// (support failure), as does a replicate without Reference Group. In a replicate that does not fail,
// a level with Reference Group weight in the full sample and none in the replicate, while Focal Group has
// it there, has its Focal Group observations set aside (another estimand, silently). Locals:
// pf_R (replicates: the variables, or the PSUs), pf_nfail and pf_nvan (replicates failing,
// with a vanished level), pf_ng2 (replicates without Reference Group), pf_fail and pf_van ("t:l:n:o",
// term t of tvlist, level l, number of replicates, and with psuv a dataset observation of
// the PSU deleted in the first such replicate, else 0). chmode (optional, chain's pairs):
// a replicate without Focal Group fails without adding level items (its value of by()
// has no weight), and pf_g1ids and pf_g2ids list the replicates without Focal Group and Reference Group
void seldecomp_prefail(string scalar tvlist, string scalar nlevlist, string scalar tousev,
    string scalar g1v, string scalar rwlist, string scalar psuv, | real scalar chmode)
{
    real matrix X, U, P2, IN1, IN2, FR, F, V, FO, VO
    real rowvector nlev, k1
    real colvector g, idx, psu, cand, pos, z, keep, s, rw, in2, lv, fid, vid, g1id, g2id
    string rowvector rv
    real scalar T, L, R, nrep, r, t, l, i, q, n, nf, nv, ng2, g2, ob, nog1
    string scalar fl, vl

    if (args() < 7) chmode = 0
    X = st_data(., tokens(tvlist), tousev)
    g = st_data(., g1v, tousev)
    idx = selectindex(st_data(., tousev) :!= 0)
    nlev = strtoreal(tokens(nlevlist))
    T = cols(X)
    L = max(nlev)
    P2 = J(T, L, 0)
    for (t = 1; t <= T; t++) {
        s = select(X[., t], g :== 0)
        if (rows(s)) P2[t, s'] = J(1, rows(s), 1)
    }
    if (psuv != "") {
        psu = st_data(., psuv, tousev)
        R = rows(uniqrows(select(psu, psu :< .)))
        cand = J(0, 1, .)
        for (t = 1; t <= T; t++) {
            for (q = 0; q <= 1; q++) {
                U = uniqrows(select((X[., t], psu), g :== q))
                n = rows(U)
                for (i = 1; i <= n; i++) {
                    if ((i == 1 ? 1 : U[i - 1, 1] != U[i, 1]) & (i == n ? 1 : U[i + 1, 1] != U[i, 1])) {
                        if (U[i, 2] < .) cand = cand \ U[i, 2]
                    }
                }
            }
        }
        if (rows(cand)) cand = uniqrows(cand)
        nrep = rows(cand)
    }
    else {
        rv = tokens(rwlist)
        R = nrep = cols(rv)
    }
    F = V = FO = VO = J(T, L, 0)
    nf = nv = ng2 = 0
    fid = vid = g1id = g2id = J(0, 1, .)
    for (r = 1; r <= nrep; r++) {
        if (psuv != "") {
            pos = (psu :!= cand[r])
            ob = idx[selectindex(psu :== cand[r])[1]]
        }
        else {
            rw = st_data(., rv[r], tousev)
            pos = (rw :> 0) :& (rw :< .)
            ob = 0
        }
        IN1 = IN2 = J(T, L, 0)
        z = J(rows(X), 1, 0)
        for (t = 1; t <= T; t++) {
            s = select(X[., t], pos :& g :== 0)
            if (rows(s)) IN2[t, s'] = J(1, rows(s), 1)
            s = select(X[., t], pos :& g :== 1)
            if (rows(s)) IN1[t, s'] = J(1, rows(s), 1)
            in2 = IN2[t, .]'
            z = z :| !in2[X[., t], 1]
        }
        g2 = any(pos :& g :== 0)
        nog1 = (chmode ? !any(pos :& g :== 1) : 0)
        keep = pos :& (g :== 1) :& !z
        FR = J(T, L, 0)
        for (t = 1; t <= T; t++) {
            k1 = J(1, L, 0)
            s = select(X[., t], keep)
            if (rows(s)) k1[s'] = J(1, rows(s), 1)
            FR[t, .] = IN2[t, .] :& !k1
        }
        if (any(FR) | !g2 | nog1) {
            nf++
            if (!g2) ng2++
            if (!g2) g2id = g2id \ (psuv != "" ? cand[r] : r)
            if (nog1) g1id = g1id \ (psuv != "" ? cand[r] : r)
            if (!nog1) {
                FO = FO + ob :* (FR :& (F :== 0))
                F = F + FR
            }
            fid = fid \ (psuv != "" ? cand[r] : r)
        }
        else {
            FR = P2 :& !IN2 :& IN1
            if (any(FR)) {
                nv++
                VO = VO + ob :* (FR :& (V :== 0))
                V = V + FR
                vid = vid \ (psuv != "" ? cand[r] : r)
            }
        }
    }
    fl = vl = ""
    for (t = 1; t <= T; t++) {
        for (l = 1; l <= nlev[t]; l++) {
            if (F[t, l]) fl = fl + " " + strofreal(t) + ":" + strofreal(l) + ":" + strofreal(F[t, l]) + ":" + strofreal(FO[t, l], "%18.0g")
            if (V[t, l]) vl = vl + " " + strofreal(t) + ":" + strofreal(l) + ":" + strofreal(V[t, l]) + ":" + strofreal(VO[t, l], "%18.0g")
        }
    }
    st_local("pf_R", strofreal(R))
    st_local("pf_nfail", strofreal(nf))
    st_local("pf_nvan", strofreal(nv))
    st_local("pf_ng2", strofreal(ng2))
    st_local("pf_fail", strtrim(fl))
    st_local("pf_van", strtrim(vl))
    // the replicates themselves (their number, or with psuv the PSU deleted), for chain, whose
    // pairs are checked one at a time and a replicate failing in two of them counts once
    st_local("pf_fids", rows(fid) ? invtokens(strofreal(fid', "%18.0g")) : "")
    st_local("pf_vids", rows(vid) ? invtokens(strofreal(vid', "%18.0g")) : "")
    if (chmode) {
        st_local("pf_g1ids", rows(g1id) ? invtokens(strofreal(g1id', "%18.0g")) : "")
        st_local("pf_g2ids", rows(g2id) ? invtokens(strofreal(g2id', "%18.0g")) : "")
    }
}

// seldecomp_hadroute: the pre-flight check under svy brr with brrweight() variables svyset.
// Local pf_had is 1 when the prefix is using hadamard() instead of them, which its
// full-sample call cannot be told but can recognise: before the call the prefix has sorted the
// data by the svyset strata and PSU variables (strv, suv; the observation when suv is empty
// or _n), so that the caller's sort order (sortv, read at entry) is exactly those, and has
// made two working variables, over all observations: the stratum's index in that order and
// the PSU's number within the stratum. Both are looked for among the temporary variables
// (names beginning with __0); when the sort order or either variable is missing, 0, and the
// svyset weights are checked as usual
void seldecomp_hadroute(string scalar strv, string scalar suv, string scalar sortv)
{
    real colvector sid, pid, sch, pch, v
    real scalar n, i, k, fs, fp
    string scalar nm

    st_local("pf_had", "0")
    if (strv == "") return
    if (suv == "_n") suv = ""
    if (strtrim(stritrim(sortv)) != strtrim(strv + " " + suv)) return
    n = st_nobs()
    if (n < 2) return
    if (st_isstrvar(strv)) sch = (st_sdata(., strv)[2::n] :!= st_sdata(., strv)[1::(n - 1)])
    else sch = (st_data(., strv)[2::n] :!= st_data(., strv)[1::(n - 1)])
    if (suv == "") pch = J(n - 1, 1, 1)
    else if (st_isstrvar(suv)) pch = (st_sdata(., suv)[2::n] :!= st_sdata(., suv)[1::(n - 1)])
    else pch = (st_data(., suv)[2::n] :!= st_data(., suv)[1::(n - 1)])
    sid = 1 \ (1 :+ runningsum(sch))
    pid = J(n, 1, 1)
    for (i = 2; i <= n; i++) pid[i] = (sch[i - 1] ? 1 : pid[i - 1] + pch[i - 1])
    fs = fp = 0
    for (k = 1; k <= st_nvar(); k++) {
        nm = st_varname(k)
        if (substr(nm, 1, 3) != "__0" | st_isstrvar(k)) continue
        v = st_data(., k)
        if (!fs & v == sid) fs = 1
        else if (!fp & v == pid) fp = 1
        if (fs & fp) break
    }
    st_local("pf_had", strofreal(fs & fp))
}

// seldecomp_isopt: Comma_check. 1 when the word w is one of the options in bare (as -syntax-
// declares them: the leading lower-case letters and the capitals are the shortest
// abbreviation, as in DIAGnostics or SEQuentially), abbreviated to no fewer letters; else 0
real scalar seldecomp_isopt(string scalar w, string scalar bare)
{
    string rowvector o
    real scalar i, m
    o = tokens(bare)
    for (i = 1; i <= cols(o); i++) {
        m = (ustrregexm(o[i], "^[a-z]*[A-Z]+") ? strlen(ustrregexs(0)) : strlen(o[i]))
        if (strlen(w) >= m & strlen(w) <= strlen(o[i]) & substr(strlower(o[i]), 1, strlen(w)) == w) return(1)
    }
    return(0)
}

// seldecomp_other: the support check with the roles of the groups swapped, for the message of
// a support failure. On the estimation sample (tousev; every weight positive; codes 1..nlev
// in tvlist), the Reference Group observations (g1v == 0) with a level of any term that Focal Group lacks
// would get zero raked weight; raking Reference Group fails when a level of a term found in the Focal Group
// has no other Reference Group observation. All terms at once: a stage rakes a subset of them, which
// sets fewer observations aside. Local otherok: 1 when raking Reference Group would not fail
void seldecomp_other(string scalar tvlist, string scalar nlevlist, string scalar tousev,
    string scalar g1v)
{
    real matrix X
    real rowvector nlev
    real colvector g, keep, in1, in2, x
    real scalar t, ok
    X = st_data(., tokens(tvlist), tousev)
    g = st_data(., g1v, tousev)
    nlev = strtoreal(tokens(nlevlist))
    keep = (g :== 0)
    for (t = 1; t <= cols(X); t++) {
        x = X[., t]
        in1 = J(nlev[t], 1, 0)
        in1[select(x, g :== 1), 1] = J(sum(g :== 1), 1, 1)
        keep = keep :& in1[x, 1]
    }
    ok = any(keep)
    for (t = 1; t <= cols(X) & ok; t++) {
        x = X[., t]
        in1 = in2 = J(nlev[t], 1, 0)
        in1[select(x, g :== 1), 1] = J(sum(g :== 1), 1, 1)
        in2[select(x, keep), 1] = J(sum(keep), 1, 1)
        ok = !any(in1 :& !in2)
    }
    st_local("otherok", strofreal(ok))
}

// seldecomp_recode: insupport. One term's codes (codev) renumbered 1, 2, ... over the levels
// still present on tousev, in their order (so as seldecomp_group numbers them on that sample),
// missing elsewhere; the level matrix lvname keeps the rows of those levels, and the local
// rkeep lists their former numbers
void seldecomp_recode(string scalar codev, string scalar tousev, string scalar lvname)
{
    real colvector x, s, map
    real matrix LV
    LV = st_matrix(lvname)
    x = st_data(., codev, tousev)
    s = uniqrows(x)
    map = J(rows(LV), 1, .)
    map[s] = 1::rows(s)
    st_store(., codev, J(st_nobs(), 1, .))
    st_store(., codev, tousev, map[x, 1])
    st_matrix(lvname, LV[s, .])
    st_local("rkeep", invtokens(strofreal(s')))
}

// seldecomp_failcache: what failed the replicates of a replication prefix, kept in the external
// _seldecomp_failcache (a k x 3 string matrix: kind, text, the number of replicate calls) for
// the one message after the prefix's loop (Misreps_warn). kind: "support" (text: the first
// level without Focal Group weight, SMCL braces escaped), "iterate" (a margin missed, the raking
// stopped at the iterate() limit; text: the stage), "targets" (a margin missed, the raking
// converged), "group" (by() took one value; text: its name), "chain" (chain: a value of by()
// of the full-sample call missing; text: "<by> = <value>"), "other" (text: the error
// message). "reset" (every call but a replicate): drop the record; "add" (a replicate call
// that fails): count the cause; "get": locals fc_n (causes), fc_kind#, fc_text#, fc_count# in
// the order first seen
void seldecomp_failcache(string scalar mode, string scalar kind, string scalar text)
{
    pointer(transmorphic) scalar p
    string matrix F
    real scalar i

    p = findexternal("_seldecomp_failcache")
    if (mode == "reset") {
        if (p != NULL) rmexternal("_seldecomp_failcache")
        return
    }
    if (mode == "add") {
        if (p == NULL) {
            p = crexternal("_seldecomp_failcache")
            *p = J(0, 3, "")
        }
        F = *p
        for (i = 1; i <= rows(F); i++) {
            if (F[i, 1] == kind & F[i, 2] == text) {
                F[i, 3] = strofreal(strtoreal(F[i, 3]) + 1)
                *p = F
                return
            }
        }
        *p = F \ (kind, text, "1")
        return
    }
    F = (p == NULL ? J(0, 3, "") : *p)
    st_local("fc_n", strofreal(rows(F)))
    for (i = 1; i <= rows(F); i++) {
        st_local("fc_kind" + strofreal(i), F[i, 1])
        st_local("fc_text" + strofreal(i), F[i, 2])
        st_local("fc_count" + strofreal(i), F[i, 3])
    }
}

// seldecomp_gencache: the generate() variables made by a full-sample call under a replication
// prefix, recorded with the cut points (in _seldecomp_cutcache, same passes) under the key in
// the caller's local cckey, with the estimation sample's signature (seldecomp_gcsig). "check":
// the local gcok is 1 when every variable in vars was made by an earlier full-sample call of
// the current pass with the same specification and sample (a later one may then replace them;
// the prefixes rewrite if/in and the weight as temporary variables between these calls, so
// the sample is compared, not its text); "put": record the variables in vars
void seldecomp_gencache(string scalar mode, string scalar vars)
{
    pointer(transmorphic) scalar p
    string rowvector v
    string scalar key, tab
    real scalar i
    real rowvector sig
    transmorphic s

    tab = char(9)
    key = st_local("cckey") + tab + "generate()" + tab
    p = findexternal("_seldecomp_cutcache")
    v = tokens(vars)
    sig = seldecomp_gcsig(st_local("touse"), st_local("w"), st_local("depvar"))
    if (mode == "check") {
        st_local("gcok", "0")
        if (p == NULL) return
        if (asarray(*p, tab + "phase") == 1) return
        for (i = 1; i <= cols(v); i++) {
            if (!asarray_contains(*p, key + v[i])) return
            s = asarray(*p, key + v[i])
            if (eltype(s) != "real" | cols(s) != cols(sig)) return
            if (mreldif(s, sig) > 1e-12) return
        }
        st_local("gcok", "1")
        return
    }
    // put: after a replicate, a new pass (as in seldecomp_cutcache)
    if (p != NULL) {
        if (asarray(*p, tab + "phase") == 1) {
            rmexternal("_seldecomp_cutcache")
            p = NULL
        }
    }
    if (p == NULL) {
        p = crexternal("_seldecomp_cutcache")
        *p = asarray_create()
        asarray(*p, tab + "phase", 0)
    }
    for (i = 1; i <= cols(v); i++) asarray(*p, key + v[i], sig)
}

// seldecomp_gcsig: the signature of the estimation sample for seldecomp_gencache: its size and
// the sums of the weights and of the weighted dependent variable, which do not depend on the
// sort order (compared with a tolerance, as the order of summation may differ)
real rowvector seldecomp_gcsig(string scalar touse, string scalar w, string scalar y)
{
    real colvector wv

    wv = st_data(., w, touse)
    return((rows(wv), quadsum(wv), quadsum(wv :* st_data(., y, touse))))
}

// seldecomp_cutcache: the contvars() cut points under a replication prefix, kept in the
// external _seldecomp_cutcache (an asarray) under the key in the caller's local cckey (the
// specification). "put" (the prefix's full-sample call): record e(cuts#) for the variables in
// vars (e(cutvars)); "get" (a replicate call): the cut points of the variables in vars go to
// the matrices named in mats, and the local ccerr is "" (found), "none" or "ambiguous";
// "loop" (a replicate call without them); "reset" (a call without a prefix): drop the record.
// A pass is the full-sample calls before the first replicate: the first full-sample call after
// a replicate starts a new record. Within a pass, two full-sample calls with the same key and
// different cut points (a program calling seldecomp on two samples) make the key ambiguous,
// and its replicates then fail rather than use the wrong cut points
void seldecomp_cutcache(string scalar mode, string scalar vars, string scalar mats)
{
    pointer(transmorphic) scalar p
    string rowvector v, m
    string scalar key, tab
    real scalar i, st
    real rowvector c

    tab = char(9)
    key = st_local("cckey")
    p = findexternal("_seldecomp_cutcache")
    if (mode == "reset") {
        if (p != NULL) rmexternal("_seldecomp_cutcache")
        return
    }
    if (mode == "loop" | mode == "get") {
        st_local("ccerr", "none")
        if (p == NULL) return
        asarray(*p, tab + "phase", 1)
        if (mode == "loop" | !asarray_contains(*p, key)) return
        if (asarray(*p, key) == 2) {
            st_local("ccerr", "ambiguous")
            return
        }
        v = tokens(vars)
        m = tokens(mats)
        for (i = 1; i <= cols(v); i++) {
            if (!asarray_contains(*p, key + tab + v[i])) return
        }
        for (i = 1; i <= cols(v); i++) st_matrix(m[i], asarray(*p, key + tab + v[i]))
        st_local("ccerr", "")
        return
    }
    // put: after a replicate, a new pass
    if (p != NULL) {
        if (asarray(*p, tab + "phase") == 1) {
            rmexternal("_seldecomp_cutcache")
            p = NULL
        }
    }
    if (p == NULL) {
        p = crexternal("_seldecomp_cutcache")
        *p = asarray_create()
        asarray(*p, tab + "phase", 0)
    }
    v = tokens(vars)
    if (asarray_contains(*p, key)) {
        st = asarray(*p, key)
        for (i = 1; i <= cols(v); i++) {
            c = st_matrix("e(cuts" + strofreal(i) + ")")
            if (!asarray_contains(*p, key + tab + v[i])) st = 2
            else if (asarray(*p, key + tab + v[i]) != c) st = 2
        }
        asarray(*p, key, st)
        return
    }
    for (i = 1; i <= cols(v); i++) {
        asarray(*p, key + tab + v[i], st_matrix("e(cuts" + strofreal(i) + ")"))
    }
    asarray(*p, key, 1)
}

// seldecomp_chainlev: chain under a replication prefix: the values of by() of
// the prefix's full-sample call, kept in _seldecomp_cutcache (the passes of seldecomp_cutcache)
// under the key in the caller's local cckey. "put": record the caller's local glev (a further
// full-sample call with the same key and other values adds its set); "check" (a replicate
// call): the local clok is 1 when glev equals a recorded set or nothing is recorded, else 0,
// and clmiss then lists the values of the first recorded set that glev lacks
void seldecomp_chainlev(string scalar mode)
{
    pointer(transmorphic) scalar p
    string scalar key, tab, g, miss
    string colvector S
    string rowvector a, b
    real scalar i

    tab = char(9)
    key = st_local("cckey") + tab + "chain levels"
    g = st_local("glev")
    p = findexternal("_seldecomp_cutcache")
    if (mode == "check") {
        st_local("clok", "1")
        st_local("clmiss", "")
        if (p == NULL) return
        asarray(*p, tab + "phase", 1)
        if (!asarray_contains(*p, key)) return
        S = asarray(*p, key)
        if (anyof(S, g)) return
        st_local("clok", "0")
        a = tokens(S[1])
        b = tokens(g)
        miss = ""
        for (i = 1; i <= cols(a); i++) {
            if (!anyof(b, a[i])) miss = miss + " " + a[i]
        }
        st_local("clmiss", strtrim(miss))
        return
    }
    // put: after a replicate, a new pass
    if (p != NULL) {
        if (asarray(*p, tab + "phase") == 1) {
            rmexternal("_seldecomp_cutcache")
            p = NULL
        }
    }
    if (p == NULL) {
        p = crexternal("_seldecomp_cutcache")
        *p = asarray_create()
        asarray(*p, tab + "phase", 0)
    }
    if (!asarray_contains(*p, key)) asarray(*p, key, g)
    else {
        S = asarray(*p, key)
        if (!anyof(S, g)) asarray(*p, key, S \ g)
    }
}

// seldecomp_chainnotes: the notes of chain's pairs (Chain). "reset"; "add" (after pair j's
// call): its e(note#), a note with the lines under it (they start with four blanks) as one
// item; "get": every item once, in the order first met, its first line headed by the pairs
// that have it ("p1, p3: "), in the locals cnote1, cnote2, ... (cnote_n of them)
void seldecomp_chainnotes(string scalar mode, real scalar j)
{
    pointer(transmorphic) scalar p
    string matrix N
    string colvector u
    string rowvector it
    string scalar s, cur, pl
    real scalar n, m, i, k, q

    p = findexternal("_seldecomp_chainnotes")
    if (mode == "reset") {
        if (p != NULL) rmexternal("_seldecomp_chainnotes")
        return
    }
    if (p == NULL) {
        p = crexternal("_seldecomp_chainnotes")
        *p = J(0, 2, "")
    }
    if (mode == "add") {
        n = st_numscalar("e(N_notes)")
        cur = ""
        for (m = 1; m <= n; m++) {
            s = st_global("e(note" + strofreal(m) + ")")
            if (substr(s, 1, 4) == "    " & cur != "") cur = cur + char(10) + s
            else {
                if (cur != "") *p = *p \ (strofreal(j), cur)
                cur = s
            }
        }
        if (cur != "") *p = *p \ (strofreal(j), cur)
        return
    }
    N = *p
    u = J(0, 1, "")
    for (i = 1; i <= rows(N); i++) {
        if (!anyof(u, N[i, 2])) u = u \ N[i, 2]
    }
    k = 0
    for (i = 1; i <= rows(u); i++) {
        pl = invtokens(("p" :+ select(N[., 1], N[., 2] :== u[i]))', ", ")
        it = ustrsplit(u[i], "\n")
        for (q = 1; q <= cols(it); q++) {
            k++
            st_local("cnote" + strofreal(k), q == 1 ? pl + ": " + it[q] : it[q])
        }
    }
    st_local("cnote_n", strofreal(k))
    rmexternal("_seldecomp_chainnotes")
}

// seldecomp_chaincnt: chain's notes of the observations left out (the probe): the count of
// the observations marked by markv at each value of the by() variable byv, as "year = 1999:
// 4; 2001: 5", or "in 12 values of year" beyond six values
string scalar seldecomp_chaincnt(string scalar byv, string scalar markv)
{
    real colvector x, u
    real scalar i
    string scalar s

    x = st_data(., byv, markv)
    u = uniqrows(x)
    if (rows(u) > 6) return("in " + strofreal(rows(u)) + " values of " + byv)
    s = byv + " = "
    for (i = 1; i <= rows(u); i++) {
        s = s + (i > 1 ? "; " : "") + strtrim(strofreal(u[i], "%18.0g")) + ": " +
            strtrim(strofreal(sum(x :== u[i]), "%12.0fc"))
    }
    return(s)
}
end

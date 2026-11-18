* run_tests.do -- runs the complete seldecomp test suite
*
* Usage, from Stata with the working directory set to <repo>/tests:
*       . do run_tests.do
* (tests/README.md, "Running").
*
* seldecomp.ado is looked up in the directory named by the global
* SELDECOMP_ROOT if that global is set, and otherwise in the parent of the
* current working directory.  So either run this do-file from the tests
* directory, or set the global first:
*       . global SELDECOMP_ROOT "/path/to/seldecomp"
*       . cd /path/to/seldecomp/tests
* (the data and expected-value files are always read relative to the tests
* directory, so the working directory must be the tests directory either way).
*
* Requirements: Stata 16 or later.  oaxaca and ebalance (SSC) are needed only by
* test_oaxaca.do and test_ebalance.do; each of those two files is reported as SKIPPED
* when its command is not installed.
*
* Every test file is run with -capture noisily- so that all files run even if
* one fails; the summary at the end lists the failures and the do-file then
* exits with a nonzero return code.  A run in which a file was skipped ends with
* SELDECOMP TESTS PASSED (... skipped) rather than ALL SELDECOMP TESTS PASSED, so a
* full pass is only the latter.

clear all
set more off
set varabbrev off
version 16

* ---- locate seldecomp.ado -----------------------------------------------------
if `"$SELDECOMP_ROOT"' != "" adopath ++ `"$SELDECOMP_ROOT"'
else adopath ++ `"`c(pwd)'/.."'
capture which seldecomp
if _rc {
    di as err "seldecomp.ado not found."
    di as err "Run this do-file from the tests directory (cd <repo>/tests), or set"
    di as err `"the global SELDECOMP_ROOT to the repository root: global SELDECOMP_ROOT "<repo>""'
    exit 601
}
capture confirm file "data/seldecomp_test.csv"
if _rc {
    di as err "data/seldecomp_test.csv not found: run this do-file from the tests directory"
    exit 601
}

* ---- dependencies -------------------------------------------------------------
capture which oaxaca
local have_oaxaca = (_rc == 0)
if !`have_oaxaca' {
    di as txt "note: oaxaca (Jann 2008) is not installed; test_oaxaca.do will be skipped."
    di as txt "      Install it with:    ssc install oaxaca"
}
capture which ebalance
local have_ebalance = (_rc == 0)
if !`have_ebalance' {
    di as txt "note: ebalance (Hainmueller and Xu 2013) is not installed; test_ebalance.do will be skipped."
    di as txt "      Install it with:    ssc install ebalance"
}

* ---- run every test file, collecting failures ---------------------------------
local files test_reference.do test_analytic.do test_oracle.do test_oaxaca.do test_options.do ///
    test_bootstrap.do test_svy.do test_plot.do test_display.do test_invariance.do        ///
    test_ebalance.do test_features.do test_variants.do test_cutpoints.do test_detail.do  ///
    test_chain.do test_edge_cases.do test_bootstrap_cleanup.do test_single_level.do test_label_characters.do ///
    test_diagnostic_counts.do test_table_labels.do test_diagnostics_display.do ///
    test_diagnostic_bounds.do test_failure_messages.do test_plot_intervals.do
local nfiles : word count `files'
local failed
local skipped
local npass 0
foreach f of local files {
    if "`f'" == "test_oaxaca.do" & !`have_oaxaca' {
        di as txt _n "SKIPPED: `f' (oaxaca not installed)"
        local skipped `skipped' `f'
        continue
    }
    if "`f'" == "test_ebalance.do" & !`have_ebalance' {
        di as txt _n "SKIPPED: `f' (ebalance not installed)"
        local skipped `skipped' `f'
        continue
    }
    di as txt _n "{hline 78}"
    di as txt "seldecomp tests: running `f'"
    di as txt "{hline 78}"
    local test_linesize = c(linesize)
    set linesize 78
    capture noisily do `f'
    local rc = _rc
    set linesize `test_linesize'
    if `rc' {
        di as err _n "FAILED: `f' (return code `rc')"
        local failed `failed' `f'
    }
    else local ++npass
}

* ---- summary ------------------------------------------------------------------
di as txt _n "{hline 78}"
local nskip : word count `skipped'
if `nskip' di as txt "skipped: `skipped'"
if "`failed'" != "" {
    di as err "seldecomp tests: `npass' of `nfiles' files passed (`nskip' skipped)"
    di as err "FAILED FILES: `failed'"
    di as txt "{hline 78}"
    exit 9
}
if `nskip' di as res "SELDECOMP TESTS PASSED (`npass' of `nfiles' files; `nskip' skipped)"
else di as res "ALL SELDECOMP TESTS PASSED (`npass' files)"
di as txt "{hline 78}"

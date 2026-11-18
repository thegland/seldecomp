* display_cases.do -- the runs whose output is kept in reference/transcripts/<case>.txt
*
* Called by test_display.do (section 5):   do display_cases.do "<directory>"
* Each case writes one text log, <directory>/<case>.txt: the command, as a paragraph that
* starts with ". ", then everything seldecomp prints. test_display.do compares these files
* line by line with reference/transcripts/, ignoring only trailing blanks; with the global
* SELDECOMP_REGEN_TRANSCRIPTS set to 1 it writes them there instead (see tests/README.md).
* The cases are showcase runs plus the display options (sequentially, fixothermargins,
* bayesboot()). Each transcript follows the order of the display: the header (and its
* notes), then the raking diagnostics, then the
* covariate detail when asked for, then the decomposition (its title line, the stage table,
* the main table, the CI table with citable), each section after a rule of "=" (Display_sep).
* Nothing here depends on the Stata version (no graphs, no prefix output, no file paths).

version 16
args outdir
if `"`outdir'"' == "" {
    di as err "usage: do display_cases.do <directory>"
    exit 198
}
local ls0 = c(linesize)
local cf0 = c(cformat)
local pf0 = c(pformat)
local sf0 = c(sformat)
set linesize 80
set cformat
set pformat
set sformat

capture program drop Sd_case
program define Sd_case
    // Sd_case <directory> <case> <command>: <command> and its output in <directory>/<case>.txt
    gettoken dir 0 : 0
    gettoken name 0 : 0
    capture log close sdcase
    quietly log using `"`dir'/`name'.txt"', text replace name(sdcase) nomsg
    di as txt `"{p 0 2 2 78}.`macval(0)'{p_end}"'     // `0' starts with a blank
    capture noisily `0'
    local rc = _rc
    quietly log close sdcase
    if `rc' exit `rc'
end

* 1. help example 1
sysuse nlsw88, clear
Sd_case `"`outdir'"' 01_default seldecomp wage i.collgrad i.south c.ttl_exp, by(married) contvars(ttl_exp 4)

* 2. weights, covariatedetail and diagnostics
gen double wt = 1 + (age - 34) / 20
Sd_case `"`outdir'"' 02_weighted seldecomp wage i.collgrad i.south [pw=wt], by(married) covariatedetail diagnostics
drop wt

* 3. help example 2: sequentially fixothermargins (the Pylypchuk-Selden procedure) with named
* blocks and a binned variable
Sd_case `"`outdir'"' 03_sequentially_fixed seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) sequentially fixothermargins

* 4. sequentially fixothermargins with every stage bootstrapped
Sd_case `"`outdir'"' 04_sequentially_fixed_bs seldecomp wage (demo: race collgrad) (geo: south smsa), by(married) sequentially fixothermargins bayesboot(reps(100) rseed(1) nodots)

* 5. its replay at another level with the CI table
Sd_case `"`outdir'"' 05_replay seldecomp, level(90) citable

* 6. clusters and BCa (50 replicates), the dots included
Sd_case `"`outdir'"' 06_bca_cluster seldecomp wage i.collgrad i.south, by(married) bayesboot(reps(50) rseed(1) cluster(industry) bca)

* 7. sequentially alone (the other blocks free) with covariatedetail
Sd_case `"`outdir'"' 07_sequentially seldecomp wage (demo: race collgrad) (geo: south smsa) c.age, by(married) contvars(age 5) sequentially covariatedetail

* 8. redisplay after the jackknife prefix: t with the prefix's df
quietly jackknife, cluster(industry) nodots: seldecomp wage i.collgrad i.south if industry < ., by(married)
Sd_case `"`outdir'"' 08_jackknife_replay seldecomp

* 9. notes: a level absent from the Reference Group (change of support) and thin target levels
gen byte race2 = race
replace race2 = 9 if married == 0 & _n <= 30
Sd_case `"`outdir'"' 09_notes seldecomp wage (a: collgrad) (b: race2) (c: south), by(married) sequentially fixothermargins mincell(300)
drop race2

* 10. a share that rounds to zero is 0.0%, never -0.0%
clear
set obs 2000
gen byte grp = _n > 1000
gen byte x = cond(grp, _n <= 1498, _n <= 502)
gen byte z = mod(_n, 2)
gen double y = 10 * grp + x
Sd_case `"`outdir'"' 10_zero_share seldecomp y (a: x) (b: z), by(grp) sequentially fixothermargins

set linesize `ls0'
set cformat `cf0'
set pformat `pf0'
set sformat `sf0'

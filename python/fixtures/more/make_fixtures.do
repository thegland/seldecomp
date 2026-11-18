* make_fixtures.do -- runs seldecomp (Stata) on the cases of cases.txt and writes, for each, the
* results that the Python tests compare with (sd_dump.do of the folder above writes the JSON).
*   <name> @<dataset> : <seldecomp command>     on data/<dataset>.dta; dataset t is the test
*                                               data with the variables of sd_data.do
*   + <name> : <command>                        run on the results of the case before it,
*                                               when that case ran without error
* Run in Stata from a copy of the fixtures folder, in its subfolder more, with the globals
*     SD_ADO   the folder that holds seldecomp.ado
*     SD_DATA  the path of tests/data/seldecomp_test.csv
*     SD_OUT   the folder to write to
* Then pack.py <SD_OUT> gathers the files into stata.json.gz and weights.npz.

clear all
set more off
set varabbrev off
version 16
adopath ++ "$SD_ADO"
do ../sd_dump.do
do ../sd_data.do

tempname ch
file open `ch' using cases.txt, read text
file read `ch' line
local baseok 0
while r(eof) == 0 {
    local line = strtrim(`"`macval(line)'"')
    if substr(`"`macval(line)'"', 1, 1) == "+" {
        if `baseok' {
            local line = strtrim(substr(`"`macval(line)'"', 2, .))
            sd_follow `macval(line)'
        }
    }
    else if `"`macval(line)'"' != "" & substr(`"`macval(line)'"', 1, 1) != "*" {
        gettoken head rest : line, parse(":")
        gettoken name ds : head, parse("@")
        local name = strtrim("`name'")
        local ds = strtrim(subinstr("`ds'", "@", "", 1))
        if "`ds'" == "t" sd_data
        else use "data/`ds'.dta", clear
        ereturn clear
        sd_case `name' `macval(rest)'
        local baseok = ("`e(cmd)'" == "seldecomp")
        // generated weights (generate(rk)) and replicate weights (bayesboot(generate(bw)))
        capture unab gl : rk*
        if !_rc {
            format rk* %24.17g
            export delimited id `gl' using "$SD_OUT/`name'_weights.csv", replace datafmt
        }
        capture unab bl : bw_*
        if !_rc & `baseok' {
            format bw_* %24.17g
            export delimited id `bl' if e(sample) using "$SD_OUT/`name'_bsweights.csv", replace datafmt
        }
    }
    file read `ch' line
}
file close `ch'

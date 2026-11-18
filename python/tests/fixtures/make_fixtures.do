* make_fixtures.do -- runs seldecomp (Stata) on the test data and writes the results that
* the Python tests compare with: one JSON file per case, and CSV files of the generated
* weights and of the bootstrap replicate weights.
*
* Run in Stata from a copy of this folder, with the globals
*     SD_ADO   the folder that holds seldecomp.ado
*     SD_DATA  the path of tests/data/seldecomp_test.csv
*     SD_OUT   the folder to write to
* The cases are listed in cases.txt, one per line:  <name> : <seldecomp command>
* A line that starts with + is run on the results of the case before it (a replay, estat
* bootstrap):  + <name> : <command>

clear all
set more off
set varabbrev off
version 16
adopath ++ "$SD_ADO"
do sd_dump.do
do sd_data.do

tempname ch
file open `ch' using cases.txt, read text
file read `ch' line
while r(eof) == 0 {
    local line = strtrim(`"`macval(line)'"')
    if substr(`"`macval(line)'"', 1, 1) == "+" {
        // a follow-up on the results of the case before it
        local line = strtrim(substr(`"`macval(line)'"', 2, .))
        sd_follow `macval(line)'
    }
    else if `"`macval(line)'"' != "" & substr(`"`macval(line)'"', 1, 1) != "*" {
        gettoken name rest : line, parse(":")
        local name = strtrim("`name'")
        sd_data
        sd_case `macval(line)'
        // generated weights (generate(rk)) and replicate weights (bayesboot(generate(bw)))
        capture unab gl : rk*
        if !_rc {
            format rk* %24.17g
            export delimited id `gl' using "$SD_OUT/`name'_weights.csv", replace datafmt
        }
        capture unab bl : bw_*
        if !_rc {
            format bw_* %24.17g
            export delimited id `bl' if e(sample) using "$SD_OUT/`name'_bsweights.csv", replace datafmt
        }
    }
    file read `ch' line
}
file close `ch'

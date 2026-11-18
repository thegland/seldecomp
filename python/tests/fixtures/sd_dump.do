* sd_dump.do -- programs that write the results of a seldecomp command to a JSON file, for
* the tests of the Python version. Numbers are written as strings with 17 significant
* digits (a missing value as null), so that they are read back exactly.
*
*   sd_case <name> : <command>       runs the command, writes <name>.json to $SD_OUT: the
*                                    return code, the printed output, and e() when it ran
*   sd_follow <name> : <command>     the same for a command run on the results in memory
*                                    (a replay, estat bootstrap): return code and output
*   sd_num <expression>              c_local sd_n: the JSON text of a number

capture program drop sd_num
program define sd_num
    args x
    if missing(`x') c_local sd_n null
    else c_local sd_n `""`=strtrim(string(`x', "%24.17g"))'""'
end

capture program drop sd_str
program define sd_str
    // c_local sd_s: text as a JSON string (backslash and double quote escaped)
    args s
    mata: st_local("s", subinstr(subinstr(st_local("s"), "\", "\\"), char(34), "\" + char(34)))
    // an apostrophe or a backquote next to the quotes of the macro would end the text early:
    // they are written as JSON escapes
    mata: st_local("s", subinstr(subinstr(st_local("s"), char(39), "\u0027"), char(96), "\u0060"))
    c_local sd_s `""`macval(s)'""'
end

capture program drop sd_matrix
program define sd_matrix
    // sd_matrix <handle> <matrix name> : {"rows": [...], "cols": [...], "data": [[...]]}
    args fh M
    local rn : rowfullnames `M'
    local cn : colfullnames `M'
    file write `fh' `"{"rows": ["'
    local i 0
    foreach r of local rn {
        local ++i
        file write `fh' `"`=cond(`i' > 1, ", ", "")'"`r'""'
    }
    file write `fh' `"], "cols": ["'
    local i 0
    foreach c of local cn {
        local ++i
        file write `fh' `"`=cond(`i' > 1, ", ", "")'"`c'""'
    }
    file write `fh' `"], "data": ["'
    forvalues i = 1/`=rowsof(`M')' {
        file write `fh' `"`=cond(`i' > 1, ", ", "")'["'
        forvalues j = 1/`=colsof(`M')' {
            if missing(`M'[`i', `j']) local sd_n null
            else local sd_n `""`=strtrim(string(`M'[`i', `j'], "%24.17g"))'""'
            file write `fh' `"`=cond(`j' > 1, ", ", "")'`sd_n'"'
        }
        file write `fh' "]"
    }
    file write `fh' "]}"
end

capture program drop sd_case
program define sd_case
    gettoken name 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    local name = strtrim("`name'")
    local cmd = strtrim(`"`0'"')
    tempname fh M lh
    tempfile lg
    // the command's output, as printed at 80 columns
    local ls = c(linesize)
    set linesize 80
    quietly log using `"`lg'"', text replace name(sdcase)
    capture noisily `cmd'
    local rc = _rc
    quietly log close sdcase
    set linesize `ls'
    file open `fh' using `"$SD_OUT/`name'.json"', write replace text
    sd_str `"`macval(cmd)'"'
    file write `fh' `"{"case": "`name'", "cmd": `macval(sd_s)', "rc": `rc'"'
    file write `fh' `", "output": ["'
    file open `lh' using `"`lg'"', read text
    file read `lh' line
    local i 0
    while r(eof) == 0 {
        local ++i
        sd_str `"`macval(line)'"'
        file write `fh' `"`=cond(`i' > 1, ", ", "")'`macval(sd_s)'"'
        file read `lh' line
    }
    file close `lh'
    file write `fh' "]"
    if `rc' == 0 {
        // scalars
        file write `fh' `", "scalars": {"'
        local i 0
        foreach s in `: e(scalars)' {
            local ++i
            sd_num e(`s')
            file write `fh' `"`=cond(`i' > 1, ", ", "")'"`s'": `sd_n'"'
        }
        file write `fh' "}"
        // macros
        file write `fh' `", "macros": {"'
        local i 0
        foreach m in `: e(macros)' {
            if inlist("`m'", "rngstate", "cmdline", "bsfile") continue
            local ++i
            mata: st_local("v", st_global("e(`m')"))
            sd_str `"`macval(v)'"'
            file write `fh' `"`=cond(`i' > 1, ", ", "")'"`m'": `macval(sd_s)'"'
        }
        file write `fh' "}"
        // matrices
        file write `fh' `", "matrices": {"'
        local i 0
        foreach m in `: e(matrices)' {
            local ++i
            matrix `M' = e(`m')
            file write `fh' `"`=cond(`i' > 1, ", ", "")'"`m'": "'
            sd_matrix `fh' `M'
        }
        file write `fh' "}"
    }
    file write `fh' "}" _n
    file close `fh'
end

capture program drop sd_follow
program define sd_follow
    // sd_follow <name> : <command>   a command run on the results in memory (a replay,
    // estat bootstrap): writes its return code and printed output to <name>.json
    gettoken name 0 : 0, parse(":")
    gettoken colon 0 : 0, parse(":")
    local name = strtrim("`name'")
    local cmd = strtrim(`"`0'"')
    tempname fh lh
    tempfile lg
    local ls = c(linesize)
    set linesize 80
    quietly log using `"`lg'"', text replace name(sdcase)
    capture noisily `cmd'
    local rc = _rc
    quietly log close sdcase
    set linesize `ls'
    file open `fh' using `"$SD_OUT/`name'.json"', write replace text
    sd_str `"`macval(cmd)'"'
    file write `fh' `"{"case": "`name'", "cmd": `macval(sd_s)', "rc": `rc', "output": ["'
    file open `lh' using `"`lg'"', read text
    file read `lh' line
    local i 0
    while r(eof) == 0 {
        local ++i
        sd_str `"`macval(line)'"'
        file write `fh' `"`=cond(`i' > 1, ", ", "")'`macval(sd_s)'"'
        file read `lh' line
    }
    file close `lh'
    file write `fh' "]}" _n
    file close `fh'
end

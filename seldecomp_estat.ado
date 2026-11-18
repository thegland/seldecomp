*! version 1.0.0  18nov2026
*! seldecomp_estat: the estat subcommands after seldecomp (e(estat_cmd)); estat bootstrap
*! displays the results of seldecomp's bayesboot() in the layout of Stata's estat bootstrap

program define seldecomp_estat, rclass
    version 16
    if "`e(cmd)'" != "seldecomp" error 301
    gettoken key rest : 0, parse(", ")
    local lkey = strlen(`"`key'"')
    if `"`key'"' != "" & `"`key'"' == substr("bootstrap", 1, max(4, `lkey')) {
        if "`e(bsnames)'" == "" {
            // after the bootstrap: or bayesboot: prefix, Stata's own table
            if inlist("`e(prefix)'", "bootstrap", "bayesboot") {
                estat_default `0'
                return add
                exit
            }
            di as err "estat bootstrap requires results from bayesboot() or the bootstrap prefix"
            exit 321
        }
        local 0 `"`rest'"'
        syntax [, ALL NORmal Percentile BC BCA noLegend Verbose Level(string) ]
        if `"`level'"' != "" {
            di as err "{p 0 0 2}option level() not allowed: estat bootstrap uses the level of " ///
                "the last display; type seldecomp, level(#) first{p_end}"
            exit 198
        }
        if "`bca'" != "" & e(bca) != 1 {
            di as err "option bca requires results from bayesboot() with bca"
            exit 198
        }
        if "`all'`normal'`percentile'`bc'`bca'" == "" local bc bc
        seldecomp, estatboot(`all' `normal' `percentile' `bc' `bca' `legend' `verbose')
        return add
        exit
    }
    // estat vce (e(V)), and a missing subcommand: Stata's own handling
    if `"`key'"' == "vce" | `"`key'"' == "" | `"`key'"' == "," {
        estat_default `0'
        return add
        exit
    }
    di as err `"{bf:estat} {bf:`key'} not valid"'
    exit 321
end

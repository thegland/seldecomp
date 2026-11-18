* sd_data.do -- program sd_data: the test data ($SD_DATA, the path of
* tests/data/seldecomp_test.csv) with the derived variables that the fixture cases use.
* The Python tests build the same variables (tests/common.py).

capture program drop sd_data
program define sd_data
    // the test data and the derived variables that the cases use (the Python tests build the
    // same variables: tests/common.py)
    import delimited using "$SD_DATA", clear varnames(1) case(preserve) stringcols(8) asdouble
    gen byte g3 = mod(id, 3)
    gen double hrs = 20 + mod(id, 25)
    gen double earn = y * hrs
    gen byte xs = x1
    replace xs = 9 if group == 0 & mod(id, 40) == 0
    gen byte xr = x1
    replace xr = 8 if group == 1 & mod(id, 50) == 0
    gen double wz = wt
    replace wz = 0 if mod(id, 37) == 0
    replace wz = . if mod(id, 91) == 0
    gen double ym = y
    replace ym = . if mod(id, 53) == 0
    gen byte ta = x2
    gen byte tb = 1 - x2
    gen byte xc = x1
    replace xc = 7 if g3 == 0 & mod(id, 30) == 0
    gen str5 gs = cond(group == 0, "alpha", "beta")
    gen byte xj = x1
    replace xj = 5 if id == 9 | (group == 1 & mod(id, 20) == 0)
    gen double hm = hrs
    replace hm = . if mod(id, 61) == 0
    gen double cm = cluster
    replace cm = . if mod(id, 71) == 0
    gen double xt = round(x3)
    label define grp 0 "Group zero" 1 "Group one"
    label values group grp
    label variable y "Outcome"
    label variable x1 "First factor"
    label variable x2 "Second factor"
end

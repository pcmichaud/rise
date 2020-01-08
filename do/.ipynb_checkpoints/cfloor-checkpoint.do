* computing cmin time series from Moffitt (2002)
* http://www.econ2.jhu.edu/people/moffitt/ben_doc.pdf

* selection of variable, the real modified benefit is in 1996 dollars
* same var as in Scholz et al. (2006)
import delimited "$src/raw/other/ben_data_moffitt.csv", clear

* first collapse by year using pop as weight
replace cmin_real = . if cmin_real == -1
collapse cmin_real [fw=pop], by(year)

* adjust to 2005 dollars
replace cmin_real = cmin_real * (194.82/156.9)

* since benefit is monthly, make annual
replace cmin_real = cmin_real * 12
gen id = 1
reshape wide cmin_real, i(id) j(year)
forvalues j = 1961/1998 {
	local j1 = `j'-1
	replace cmin_real`j' = cmin_real`j1' if cmin_real`j'==. & cmin_real`j1'!=.
}

forvalues j = 1999/2060 {
	gen cmin_real`j' = cmin_real1998
} 
reshape long cmin_real, i(id) j(year)
drop if year<1965
gen age = year - 1940

sum cmin_real if year<=2005
egen cmin_mean = mean(cmin_real)
#d ;
twoway (line cmin_real year if year<=2005, lpattern(solid)) 
(line cmin_mean year if year<=2005, lpattern(dash)),
	xtitle("year") ytitle("ressource floor")
	$fig
	legend(label(1 "real ressource floor") label(2 "average"));
#d cr
graph export figures/xfloor.eps, as(eps) replace
list age cmin_real
keep cmin_real
outsheet using "params/input/cmin.csv", replace comma nonames nolabel

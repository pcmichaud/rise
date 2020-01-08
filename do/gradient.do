	* Statistics for Table 9
	capture log close
    set more off
    
    use $src/clean/psid_extract.dta, clear
	* sample cuts
	keep if rage>=35&rage<=84
	keep if inlist(year,2005)
    replace ass = 0 if ass<0
	drop if rshlt3==.
	gen poor = rshlt3==1 if rshlt3!=.
	egen qw = cut(ass), at(0 100e3 350e3 500e3 1e6)
	tab qw
	drop if missing(hhidpn)
    tabstat poor [aw=wgid_cs] if rage>=55&rage<=65, by(qw)
	

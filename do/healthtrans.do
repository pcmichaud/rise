*************************************************
* Health Profiles 				    			*
*************************************************


use $src/clean/psid_extract.dta, clear
* sample cuts
keep if rage>=35&rage<=84
egen rage5 = cut(rage), at(35(5)85) 
drop if rshlt3==.
drop if missing(hhidpn)
char rbyr_c[omit] 3
keep if year<=1997
tsset hhidpn year
gen frshlt3 = f.rshlt3
label def hlth 1 "poor" 2 "good" 3 "vgood"
label values frshlt3 hlth
gen poor = frshlt3==1 if frshlt3!=.
gen good = frshlt3==2 if frshlt3!=.
gen vgood = frshlt3==3 if frshlt3!=.

tab rshlt3 frshlt3, row nofreq

foreach num of numlist 1 3 {
	graph bar (mean) poor good vgood if rshlt3==`num', $fig over(rage5) legend(rows(1) label(1 "poor") label(2 "good") label(3 "very good")) stack ytitle("fraction health(t+1)") 
	graph export "figures/healthtrans_psid_`num'.eps", as(eps) replace
}


use $src/temp/meps-nhis_long_extract.dta, clear
* sample cuts
keep if age>=35&age<=84
egen age5 = cut(age), at(35(5)85) 
drop if hlth==.
rename dupersid hhidpn
destring hhidpn, replace
drop if missing(hhidpn)
tsset hhidpn t
replace hlth = hlth-1
gen fhlth = f.hlth
label def hlth 1 "poor" 2 "good" 3 "vgood"
label values fhlth hlth
gen poor = fhlth==1 if fhlth!=.
gen good = fhlth==2 if fhlth!=.
gen vgood = fhlth==3 if fhlth!=.

tab hlth fhlth, row nofreq

foreach num of numlist 1 3 {
	graph bar (mean) poor good vgood if hlth==`num', $fig over(age5) legend(rows(1) label(1 "poor") label(2 "good") label(3 "very good")) stack ytitle("fraction health(t+1)") 
	graph export "figures/healthtrans_meps_`num'.eps", as(eps) replace
}

exit

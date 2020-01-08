*************************************************
* Moments used in Estimation (and others) 	    *
*************************************************


***********************************
* Step 1 PSID
***********************************

preserve 
	use $src/clean/psid_extract.dta, clear
	* sample cuts
	keep if rage>=35&rage<=84
	replace ass = 0 if ass<0
	drop if rshlt3==.
	drop if missing(hhidpn)

	char rbyr_c[omit] 3
	char rage[omit] 35
	* wealth
	gen pass = .
	xi: xtreg ass i.rage , fe
	gen insamp = e(sample)
	matrix beta = e(b)'
	predict u if insamp, ue
	qui sum u [aw=wgid_cs] if rbyr_c==3&insamp
	local mu = r(mean)
	* compute adjusted prediction
	matrix bcons = beta["_cons",1]
	qui sum u [aw=wgid_cs] if insamp
	qui replace pass = bcons[1,1] + (u - r(mean)) + `mu'  if rage==35&insamp

	forvalues t = 36(1)84 {
			matrix bage = beta["_Irage_`t'",1]
			qui replace pass = (u - r(mean)) +  bage[1,1] + bcons[1,1] + `mu'  if rage==`t'&insamp
	}
	drop insamp u
	
	* pass is the adjusted wealth variable

	* work
	replace rwork = . if rage>69
	gen pwork = .
	foreach h of numlist 1/3 {
		qui xi: reg rwork i.rage i.rbyr_c [aw=wgid_cs] if rshlt3==`h'
		gen insamp = e(sample)
		predict u if rshlt3==`h'&insamp, res
		matrix beta = e(b)'
	    matrix bcons = beta["_cons",1]
	    replace pwork = u + bcons[1,1] if rage==35&rshlt3==`h'&insamp
		forvalues t = 36(1)69 {
			matrix bage = beta["_Irage_`t'",1]
			qui replace pwork = u + bage[1,1] + bcons[1,1] if rage==`t'&rshlt3==`h'&insamp
		}	
		drop u	
		drop insamp
	}
	* pwork has the adjusted data for work

	* very good health
	gen vgood = rshlt3==3 
	reg vgood i.rage ib3.rbyr_c [aw=wgid_cs]
	gen insamp = e(sample)
	predict u if insamp, res
	matrix beta = e(b)'
	gen pvgood = .
    matrix bcons = beta["_cons",1] 
    replace pvgood = bcons[1,1] + u if rage==35&insamp
	forvalues t = 36(1)84 {
		matrix bage = beta["`t'.rage",1]
		qui replace pvgood = bage[1,1] + bcons[1,1] + u if rage==`t'&insamp
	}	
	drop insamp u
	* pvgood has the adjusted data for vgood

	* poor health
	gen vpoor = rshlt3==1 
	reg vpoor i.rage ib3.rbyr_c [aw=wgid_cs]
	gen insamp = e(sample)
	predict u if insamp, res
	matrix beta = e(b)'
	gen ppoor = .
    matrix bcons = beta["_cons",1] 
    replace ppoor = bcons[1,1] + u if rage==35&insamp
	forvalues t = 36(1)84 {
		matrix bage = beta["`t'.rage",1]
		qui replace ppoor = bage[1,1] + bcons[1,1] + u if rage==`t'&insamp
	}	
	drop insamp u

	*keep if inlist(rage,35,40,45,50,55,60,65,70,75,80,84)
	tsset hhidpn year
	xtdes
	sum ppoor pvgood pwork pass
	gen psid = 1
	keep hhidpn rage ppoor pvgood pwork pass rshlt3 wgid_cs psid

	save $src/clean/psid_adjusted.dta, replace

restore


preserve
	set more off
	use $src/temp/meps-nhis_long_extract.dta, clear
	destring dupersid, replace
	rename dupersid hhidpn
	drop if missing(hhidpn)
	sort hhidpn year
	gen rbyr = year - age
	keep if age>=35
	#d ;
	recode rbyr (min/1925=1) (1925/1935=2) 
	(1936/1945=3) (1946/1955=4) (1956/1965=5) (1966/max=6), gen(rbyr_c);
	#d cr
	char rbyr_c[omit] 3
	char age[omit] 35
	tab rbyr_c , gen(rbyr_c)
	drop if hlth==.
	reg mepsexp i.age ib3.rbyr_c
	gen insamp = e(sample)
	predict u if insamp, res
	matrix beta = e(b)'			
	gen pmed = .
	matrix bcons = beta["_cons",1] 
	replace pmed = u + bcons[1,1] if age==35&insamp
	forvalues t = 36(1)84 {
		matrix bage = beta["`t'.age",1]
		qui replace pmed = u + bage[1,1] + bcons[1,1] if age==`t'&insamp
	}
	drop insamp u
	rename age rage
	rename hlth rshlt3
	replace rshlt3 = rshlt3 - 1
	rename perwt wgid_cs
	gen meps = 1
	keep hhidpn rage pmed rshlt3 wgid_cs meps
	save $src/clean/meps_adjusted.dta, replace
restore


use $src/clean/psid_adjusted.dta, clear
append using $src/clean/meps_adjusted.dta
replace meps = 0 if meps==.
replace psid = 0 if psid==.
gen age = rage
merge n:1 age using $src/other/mortality-ssa-cohort.dta
keep if _merge==3
keep if age>=25&age<=84
gen pdead = (uniform()<pmx)
drop pmx age
egen id = group(hhidpn psid meps)
drop _merge
sort id rage
by id: gen t = _n
tsset id t
xtsum id
global nobs = r(n)
save $src/clean/data_adjusted.dta, replace


use $src/clean/data_adjusted.dta, clear

local nages = 84 - 35 + 1
local nwork = 69 - 35 + 1
local i = 1
* work in poor health
forvalues j = 1/`nwork' {
	local a = 35 + `j' - 1
	gen mom_`i' = pwork  * (rage==`a'&rshlt3==1&psid==1)
	gen in_`i' = (rage==`a'&rshlt3==1&psid==1&mom_`i'!=.)
	replace mom_`i' = 0 if in_`i'==0
	local i = `i' + 1
}

* work in good health
forvalues j = 1/`nwork' {
	local a = 35 + `j' - 1
	gen mom_`i' = pwork * (rage==`a'&rshlt3==2&psid==1) 
	gen in_`i' = (rage==`a'&rshlt3==2&psid==1&mom_`i'!=.)
	replace mom_`i' = 0 if in_`i'==0 
	local i = `i' + 1
}

* work in very good health
forvalues j = 1/`nwork' {
	local a = 35 + `j' - 1
	gen mom_`i' = pwork * (rage==`a'&rshlt3==3&psid==1) 
	gen in_`i' = (rage==`a'&rshlt3==3&psid==1&mom_`i'!=.)
	replace mom_`i' = 0 if in_`i'==0
	local i = `i' + 1
}

* medical spending
forvalues j = 1/`nages' {
	local a = 35 + `j' - 1
	gen mom_`i' = pmed * (rage==`a' & meps==1) 
	gen in_`i' = (rage==`a' & meps==1 & mom_`i'!=.)
	replace mom_`i' = 0 if in_`i'==0
	local i = `i' + 1
}

* mortality
forvalues j = 1/`nages' {
	local a = 35 + `j' - 1
	gen mom_`i' = pdead * (rage==`a') 
	gen in_`i' = (rage==`a' & mom_`i'!=.)
	replace mom_`i' = 0 if in_`i'==0
	local i = `i' + 1
}

* poor
forvalues j = 1/`nages' {
	local a = 35 + `j' - 1
	gen mom_`i' = ppoor * (rage==`a') 
	gen in_`i' = (rage==`a' & psid==1 & mom_`i'!=.)
	replace mom_`i' = 0 if in_`i'==0
	local i = `i' + 1
}


* vgood
forvalues j = 1/`nages' {
	local a = 35 + `j' - 1
	gen mom_`i' = pvgood * (rage==`a') 
	gen in_`i' = (rage==`a' & psid==1 & mom_`i'!=.)
	replace mom_`i' = 0 if in_`i'==0	
	local i = `i' + 1
}

* wealth
forvalues j = 1/`nages' {
	local a = 35 + `j' - 1
	gen mom_`i' = pass * (rage==`a' & psid==1) 
	gen in_`i' = (rage==`a' & psid==1 & mom_`i'!=.)	
	replace mom_`i' = 0 if in_`i'==0
	local i = `i' + 1
}

keep mom_* in_*
order mom_* in_*

*forvalues i = 1/300 {
*	sum mom_`i' if in_`i'==1
*}



preserve
	keep mom_*
	sum mom_*
	outsheet using params/input/data_adjusted.csv, replace comma nonames nolabel
restore

preserve
	keep in_*
	sum in_*
	outsheet using params/input/select_adjusted.csv, replace comma nonames nolabel
restore




global ntot = _N
file open obs using "params/input/nobs.dat", write text replace
file write obs ("$nobs")  _n
file write obs ("$ntot")  _n
file close obs

use $src/clean/data_adjusted.dta, clear

capture program drop onemom
program onemom, eclass 
  tempname mom
  	local nages = 84 - 35 + 1
  	local nwork = 69 - 35 + 1

  	matrix mom = J(1,5*`nages' + 3*`nwork',0) 
  	local i = 1
	tabstat pwork if rshlt3==1&rage<=69, by(rage) save
	forvalues j = 1/`nwork' {
		matrix mom[1,`i'] = r(Stat`j')
		local i = `i' + 1
	}
	tabstat pwork if rshlt3==2, by(rage) save 
	forvalues j = 1/`nwork' {
		matrix mom[1,`i'] = r(Stat`j')
		local i = `i' + 1
	}
	tabstat pwork if rshlt3==3, by(rage) save 
	forvalues j = 1/`nwork' {
		matrix mom[1,`i'] = r(Stat`j')
		local i = `i' + 1
	}
	tabstat pmed , by(rage) save 
	forvalues j = 1/`nages' {
		matrix mom[1,`i'] = r(Stat`j')
		local i = `i' + 1
	}
	tabstat pdead, by(rage) save 
	forvalues j = 1/`nages' {
		matrix mom[1,`i'] = r(Stat`j')
		local i = `i' + 1
	}
	tabstat ppoor, by(rage) save 
	forvalues j = 1/`nages' {
		matrix mom[1,`i'] = r(Stat`j')
		local i = `i' + 1
	}
	tabstat pvgood, by(rage) save
	forvalues j = 1/`nages' {
		matrix mom[1,`i'] = r(Stat`j')
		local i = `i' + 1
	} 
	tabstat pass, by(rage) save
	forvalues j = 1/`nages' {
		matrix mom[1,`i'] = r(Stat`j')
		local i = `i' + 1
	}
	ereturn post mom 
end

tsset, clear
bootstrap _b, reps(100) seed(10101) nodots nowarn cluster(id): onemom

matrix result = (e(b)',e(se)')
matrix colnames result = moments se_moments
svmat result
keep if result1!=.
keep result1 result2
outsheet using params/input/moments.csv, replace comma nolabel nonames


exit

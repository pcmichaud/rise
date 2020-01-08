*****************************************************
* Extracting from MEPS for analysis
*****************************************************

global mepsout "$src/temp"

**We merge MEPS and MEPS_NHIS
**Using raw data - our construction of MEPS and MEPS_NHIS
**do files to construct MEPS9608 are in mesps_selected9608.do
**do files to construct MEPS_NHIS_all.dta are in transfermepshis.do

use "$mepsout/MEPS9608.dta", clear

* count number of interviews (drop if less than 2, balanced panel)
bysort dupersid: gen count = _N
keep if count<=2
sort dupersid yr
bysort dupersid: gen t = _n

global hlth "hearte diabe hibpe stroke lunge adl1p iadl1"
* select variables and reshape wide
#d ;
keep dupersid year t count totinc hlth age obese educ male mepsexp mepsslf perwt inscov panel $hlth; 
#d cr	

save $mepsout/meps-nhis_long_extract.dta, replace

#d ;
reshape wide year totinc hlth educ age obese male mepsexp mepsslf perwt panel inscov $hlth, i(dupersid) j(t);
#d cr

* keep only wave 1 characteristics 
drop male2 obese2 educ2 age2 perwt2 inscov2 hearte2 diabe2 hibpe2 stroke2 lunge2 adl1p2 iadl12
rename male1 male
rename age1 age
rename obese1 obese
rename totinc1 totinc
rename educ1 educ
rename year1 year
rename panel1 panel
rename perwt1 perwt
rename hlth1 hlth
rename hlth2 fhlth
rename totinc2 ftotinc
rename mepsexp2 fmepsexp
rename mepsexp1 mepsexp
rename mepsslf2 fmepsslf
rename mepsslf1 mepsslf
rename hearte1 hearte
rename diabe1 diabe
rename hibpe1 hibpe
rename stroke1 stroke
rename lunge1 lunge
rename adl1p1 adl1p
rename iadl11 iadl1

rename inscov1 inscov

* procedure to obtain mortality and smoking data from NHIS
sort dupersid year
capture drop _merge
di _N
des
tab year
tab panel
merge  1:1 dupersid panel  using "$mepsout/his_meps_all.dta"
keep if _merge==3
tab year smokev

forvalues t = 1995/2008 {
	merge n:1 publicid using $src/nhis/nhis_`t'_pmort.dta, update gen(deathmatch`t')
	tab year deathmatch`t'
	drop if deathmatch`t'==2
	drop deathmatch`t'
}
di _N


* health
* 1 = dead, 2 = poor-fair, 3 = good, 4 = vgood-excellent
gen died = dodyear==(year+1) if dodyear!=. 
replace died = . if dodyear<=year

replace fhlth = 1 if died==1
label def hlth 1 "dead" 2 "poor-fair" 3 "good" 4 "vgood-excellent"
label values hlth hlth
label values fhlth hlth

* merge life-table (age,sex,year from HMD)
gen rage = age
* gen birth year 
gen byr = year - age
	#d ;
	recode byr (min/1920=1) (1921/1930=2) (1931/1935=3) 
	(1936/1945=4) (1946/1955=5) (1956/1965=6) (1966/max=7), gen(byr_c);
	#d cr
	char byr_c[omit] 4
merge m:1 male year rage using $src/other/mortality.dta, gen(lifetable)
* drop if ages not present in MEPS-NHIS
drop if lifetable==2

* destring id
destring dupersid, replace

* sample selection and cleanup
#d ;
keep dupersid fhlth age hlth fmepsexp smoken smokev obese totinc mepsexp
	mx male obese byr byr_c educ perwt year  mepsslf inscov $hlth;
order dupersid fhlth age hlth fmepsexp smoken smokev obese totinc mepsexp
	mx male obese byr byr_c educ perwt year mepsslf inscov $hlth;
#d cr
save $src/clean/meps-nhis_wide_extract.dta, replace

exit



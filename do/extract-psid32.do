*************************************************
* Constructing PSID Panel                       *
* Inputs
* PSID Cornell Equivalence Files  1980-2005     *
* wealth files (1984,1989,1994,1999,2001(2)2005 *
* Individual files (1986, 1994-2005)            *
* Outputs
* clean/psid_extract.dta							*
*************************************************

global fig "graphregion(color(white)) bgcolor(white)"
**************************
* append equivalence files
* from PSID directory
**************************
use $src/psid/pequiv_1980.dta, clear
rensfix _1980
gen year = 1980
save $src/temp/psid_panel.dta, replace

foreach s of numlist 1981/1997 1999 2001 2003 2005{
	use $src/psid/pequiv_`s'.dta, clear
	rensfix _`s'
	gen year = `s'
	save $src/temp/psid_panel_`s', replace
	use $src/temp/psid_panel.dta, clear
	append using $src/temp/psid_panel_`s'.dta
	save $src/temp/psid_panel.dta, replace
	erase $src/temp/psid_panel_`s'.dta
}

save $src/temp/psid_panel.dta, replace

**************************
*  select sample
**************************
use $src/temp/psid_panel.dta, clear
* generate death indicator (not used)
gen death = x11103==2 if x11103!=0
* Selecting sample
* only household heads
keep if d11105==1
* age range
keep if d11101>21&d11101<100
* drop oversample
keep if x11104ll==11
save $src/temp/psid_panel.dta, replace


*****************************
* merge wealth information 
* from wealth files
*****************************

* merge back wealth data to psid_panel
use $src/temp/psid_panel.dta, clear
sort x11102 year
foreach s of numlist 1984 1989 1994 1999 2001 2003 2005 {
	capture drop _merge
	merge m:1 x11102 year using $src/psid/wlth_`s'.dta, update
	drop _merge
	sort x11102 year
}

save $src/temp/psid_panel.dta, replace

****************************
* merge CPI information (BLS)
****************************
use $src/temp/psid_panel.dta, clear
sort year
merge year using $src/other/cpi.dta, nokeep
drop c19* c20*
save $src/temp/psid_panel.dta, replace

*****************************
* rename variables
*****************************

use $src/temp/psid_panel.dta, clear
rename x11101ll hhidpn
rename x11102 hhid
rename d11103 hrace
rename d11106 hhsize
rename d11107 hhchild
rename l11101 hgstate
global ids "hhidpn hhid year hgstate hrace hhsize hhchild"

rename d11101 rage
rename d11102ll rsex
rename d11112ll rrace
rename d11104  rmstat
rename d11108 reduc
rename d11109 redyrs
global demo "rage rsex rrace rmstat reduc redyrs"

rename w11102 wghh
rename w11101 wgid_cs
rename w11103 wgid_long
global weights "wghh wgid_*"

rename i11103 hiearn
rename i11104 hicap
rename i11106 hipriv
rename i11107 hipub
rename i11108 hissben
rename i11117 hipen
rename i11118 hiothr
rename i11110 riearn
global inc "hi* riearn"

rename h11103 hnkid1
rename h11104 hnkid4
rename h11105 hnkid7
rename h11106 hnkid12
rename h11110 hnothr
rename h11112 hspouse
global hhcomp "hnkid* hnothr hspouse"

rename e11101 rhours
rename e11102 rwork
rename e11103 rpart
rename e11105 rjcocc
rename e11106 rjcind
global work "rhours rwork rpart rjcocc rjcind"

rename m11101 rhstay
rename m11102 rhdays
rename m11104 rhvigact
rename m11105 rstroke
rename m11106 rhibpe
rename m11107 rdiabe
rename m11108 rcance
rename m11109 rpsyce
rename m11110 rarthe
rename m11111 rheare
rename m11112 radl_breathe
rename m11113 radl_climb
rename m11114 radl_bathing
rename m11115 radl_dress
rename m11116 radl_bed
rename m11117 riadl_shop
rename m11118 riadl_walk
rename m11119 riadl_hwork
rename m11120 riadl_kneel
rename m11121 riadl_vig
rename m11122 rhgt
rename m11123 rwgt
rename m11124 rdis
rename m11126 rshlt
global medc "rhstay rhdays"
global risk "rhvigact rhgt rwgt"
global cond "rstroke rhibpe rdiabe rcance rpsyce rarthe rheare death "
global adls "radl_*"
global iadls "riadl_*"
global shlt "rdis rshlt" 

rename S16 hatotn
rename S17 hatota
global wlth "hatotn hatota"

***************************
* Keep variables of interest
***************************

keep $ids $demo $weights $inc $hhcomp $work $medc $risk $cond $adls $iadls $shlt $wlth cpi
order $ids $demo $weights $inc $hhcomp $work $medc $risk $cond $adls $iadls $shlt $wlth cpi
sort hhidpn year
by hhidpn: gen time = _n
by hhidpn: gen period = _N
by hhidpn: gen fyear = year[1]
by hhidpn: gen lyear = year[_N]
tsset hhidpn time

************************************
* convert to real 2005 dollars
************************************
foreach var of varlist $inc $wlth {
	replace `var' = `var'/cpi*195.3
} 
gen rbyr = year - rage
save $src/temp/psid_panel.dta, replace


******************************
* getting labor force status
* from raw individual file (not necessariy to run)
******************************

global empq "ER30353 ER30382 ER30411 ER30441 ER30474 ER30509 ER30545 ER30580 ER30616 ER30653 ER30699 ER30744 ER30816 ER33111 ER33211 ER33311 ER33411 ER33512 ER33612 ER33712 ER33813"
use ER30001 ER30002  $empq using $src/psid/indfile_fam.dta, clear
gen hhidpn = ER30001*1000 + ER30002
sort hhidpn
rename ER30353 status1981
rename ER30382 status1982
rename ER30411 status1983
rename ER30441 status1984
rename ER30474 status1985
rename ER30509 status1986
rename ER30545 status1987
rename ER30580 status1988
rename ER30616 status1989
rename ER30653 status1990
rename ER30699 status1991
rename ER30744 status1992
rename ER30816 status1993
rename ER33111 status1994
rename ER33211 status1995
rename ER33311 status1996
rename ER33411 status1997
rename ER33512 status1999
rename ER33612 status2001
rename ER33712 status2003
rename ER33813 status2005
reshape long status, i(hhidpn) 
rename _j year
gen unemp = status==3 if status!=9&status!=0
gen disab = status==5 if status!=9&status!=0
sort hhidpn year
save $src/temp/status.dta, replace

*******************************************
* update status for 1994 to 2005 
* using files created using individual data
*******************************************

use $src/temp/psid_panel.dta, clear
sort hhidpn year
merge m:1 hhidpn year using $src/temp/status.dta
drop if _merge==2
drop _merge
save $src/temp/psid_panel.dta, replace
erase $src/temp/status.dta

***********************************************
* clean variables and create final PSID dataset
***********************************************

use $src/temp/psid_panel.dta, clear

* ids
global ids "hhidpn hhid year wgid_cs"

* equivalence scale
egen totkids  = rowtotal(hnkid*)
gen hnadults = hhsize - totkids
gen eqscale = (hnadults + 0.7*totkids)^0.7
recode hhsize (6/max=6)
recode hspouse (2=.), gen(spouse)
label var spouse "spouse in household"

* Education 
gen rskill=1 if (reduc==3)
replace rskill=0 if (reduc==1|reduc==2)
label def skill 0 "no college" 1 "college"
label values rskill skill
* cohort groups (3 is reference for all analysis)
#d ;
recode rbyr (min/1925=1) (1926/1935=2) (1935/1945=3) (1946/1955=4) (1956/1965=5) (1966/max=6), gen(rbyr_c);
#d cr

char rbyr_c[omit] 3
label var rbyr_c "birth year group (1935-1945=3, ref)"

recode rsex (1=1) (2=0), gen(male)
label var male "resp. is male"
global group "rage rskill rbyr_c eqscale hhsize hnadults eqscale male spouse"

* Assets 
rename hatota ass
replace ass = . if ass>2e6
replace ass = ass

* Earnings
replace riearn = riearn
gen logriearn = log(riearn)  if riearn!=0&missing(rhours)==0
* drop the obs with earnings > 350k < 5k
replace logriearn = . if riearn>3.5e5&missing(riearn)==0
replace logriearn = . if riearn<5e3
* drop those who work too many or too few hours
replace logriearn = . if rhours>6000
replace logriearn = . if rhours<100
* wage
gen rwage = riearn / rhours if rhours>0
gen logrwage = log(rwage)
global earn "riearn logriearn rwage logrwage"

* other household income
replace hiearn = hiearn
replace hipub = hipub
replace hissben = hissben
replace hipen = hipen
replace hiothr = hiothr
gen otherinc = (hiearn -riearn) 
label var otherinc "other household income (spouse)"
global other "otherinc hipen hissben"

* Health status and mortality
replace rshlt = 6 - rshlt
gen rshlt3=.
replace rshlt3=1 if (rshlt==1|rshlt==2)
replace rshlt3=2 if (rshlt==3)
replace rshlt3=3 if (rshlt==4|rshlt==5)

label define rshlt3f 1 "bad" 2 "good" 3 "excellent"
label value rshlt3 rshlt3f
label variable rshlt3 "health status (3 cat)"

gen rshlt2=.
replace rshlt2=1 if (rshlt3==2|rshlt3==3)
replace rshlt2=0 if (rshlt3==1)
label def rshlt2 0 "fair/poor" 1 "good or better"
label values rshlt2 rshlt2

gen bmi = rwgt/(rhgt)^2
label var bmi "bmi of respondent"
gen obese = bmi>=30 if missing(bmi)==0
label var obese "respondent is obese (bmi>30)"
global health "rshlt* death bmi obese"

* claiming
sort hhidpn year
bysort hhidpn: gen n = _n
tsset hhidpn n
gen ssclaim  = f.hissben>0 if hissben==0&missing(f.hissben)==0
replace ssclaim = . if rage<=61
replace ssclaim = . if rage>=70

* generate not participating (unemployed or disabled), not used
gen nlfp = 0 if status==1|status==2
replace nlfp = 1 if status==3

gen nlfp2 = 0 if status==1|status==2
replace nlfp2 = 1 if status==3

global lfp "rwork ssclaim rpart"

order $ids $group ass $earn $other $health $lfp
keep $ids $group ass $earn $other $health $lfp
label var ass "net household assets"
label var year "year of interview"
label var death "dead next wave"
label var eqscale "equivalence scale"
label var rskill "education level, 2 cat"
label var logriearn "log earnings"
label var rshlt2 "health status (2 cat)"
label var ssclaim "claiming status (2 cat)"

* drop if age is missing
drop if rage==.

* keep males
keep if male==1

* drop if health not available
drop if rshlt==.

* merge life table information (HMD, mortality.org: michaud@rand.org, pass = usual)
merge m:1 male year rage using $src/other/mortality.dta
drop if _merge==2
label var mx "one-year mortality rate from HMD, age,gender,year"

* save final dataset
label data "Fonseca, Michaud, Kapteyn and Galama (2017) - PSID extract, ver32"
save $src/clean/psid_extract.dta, replace

exit 





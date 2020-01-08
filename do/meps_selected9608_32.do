**************************
* Consolidated MEPS data

* Period 1996-2008
* http://meps.ahrq.gov/mepsweb/
**********************************

global mepsdir  "$src/meps"
global mepsout "$src/temp"

*LIST OF FILE NUMBERS FOR ANNUAL CONSOLIDATED FILE 1996 - 2008
*global csdlist  12 20 28 38 50 60	70	79	89	97	105	113	121
global csdlist   50 60	70	79	89	97	105	113	121
*LIST OF FILE NUMBERS FOR MEDICAL CONDITIONS FILE 2000 = 2008
*Note: I changed the name hc006r to h6 for the year 1996
*global condlist 6 18 27 37 52	61	69	78	87	96	104	112	120
global condlist 52	61	69	78	87	96	104	112	120
*LIST OF TWO-DIGIT YEARS - NEEDED TO RENAME VARIABLES
*global yrlist 96 97 98 99 00 01 02 03 04 05 06 07 08 
global yrlist 00 01 02 03 04 05 06 07 08 

*******************************************
* Selected cost and disease variables from MEPS

*******************************************

***1996
fdause "$mepsdir/h12.ssp"

	ren inscope1 inscop31
	ren inscope2 inscop42
	ren marry96x marry53x

	* Rename variables
	ren inscop96 inscopend
	ren age96x age
	*ren perwt96f perwt
	ren wtdper96 perwt
	gen male = sex == 1 if sex < .
	*ren famwt96f famwt 
	ren famsze96 famsze 
	ren educyr96 educyear
	ren region96 region53
	
*** MEDICAL Expenditures;
	foreach item in exp slf mcr mcd prv va ofd wcp opr opu osr{
		ren tot`item'96 meps`item'
	}
	gen yr = 1996

*** Utilization variables ; 
	foreach item in obtotv ipdis ipngtd {  
		ren `item'96 `item' 
	}
	
	ren obtotv doctim 
	ren ipdis  hsptim
	ren ipngtd hspnit

**Health
		ren iadlhlp1 iadlhp53 
		ren adlhelp1 adlhlp53
		ren rtehlth1 rthlth53
		ren mnthlth1 mnhlth53	
	
***Total income
        ren ttlpnx  totinc	
	
***Insurance variables;
*	foreach item in triev mcrev mcdev opaev opbev unins inscov {  
*		ren `item'96 `item' 
*	}
*perwt not found for 1996
*diabdx53 hibpdx53 chddx53 angidx53 midx53 ohrtdx53 strkdx53 emphdx53 asthdx53 bmindx53 adsmok42
	gen panel96 = 1
	#d;
	keep dupersid perwt yr male educyear age racex hispanx marry53x
	meps* doctim hsptim hspnit  
	duid famidyr famrfpyr famszeyr 
	iadlhp53 adlhlp53 
	region53 
	rthlth53 mnhlth53 panel96
	totinc
	empst* selfcm* numemp* insc* held* offer* ;
	* just the last years indcat* occc*; 
	#d cr
	
    save "$mepsout/m1996.dta", replace	
	
	
*****1997

fdause "$mepsdir/h20.ssp"


	* Rename variables
	ren inscop97 inscopend
	ren age97x age
	*ren perwt97f perwt
	ren wtdper97 perwt
	gen male = sex == 1 if sex < .
	*ren famwt97f famwt 
	ren famsze97 famsze 
	ren educyr97 educyear
*** MEDICAL Expenditures;
	foreach item in exp slf mcr mcd prv va ofd wcp opr opu osr{
		ren tot`item'97 meps`item'
	}
	gen yr = 1997

*** Utilization variables ; 
	foreach item in obtotv ipdis ipngtd {  
		ren `item'97 `item' 
	}
	
	ren obtotv doctim 
	ren ipdis  hsptim
	ren ipngtd hspnit

	
***Total income
        ren ttlp97x  totinc	
		
***Insurance variables;
*	foreach item in triev mcrev mcdev opaev opbev unins inscov {  
*		ren `item'97 `item' 
*	}
*perwt  no found in 1997
*	diabdx53 hibpdx53 chddx53 angidx53 midx53 ohrtdx53 strkdx53 emphdx53 asthdx53 bmindx53 adsmok42

	#d;
	keep dupersid perwt yr male educyear age racex hispanx marry53x
	meps* doctim hsptim hspnit 
	duid famidyr famrfpyr famszeyr 
	iadlhp53 adlhlp53 
	region53 
	rthlth53 mnhlth53 panel97
	totinc
	empst* selfcm* numemp* insc* held* offer* ;
	* just the last years indcat* occc*; 
	#d cr
	
    save "$mepsout/m1997.dta", replace		
	
*****1998
fdause "$mepsdir/h28.ssp


	* Rename variables
	ren inscop98 inscopend
	ren age98x age
	*ren perwt98f perwt
	ren wtdper98 perwt
	gen male = sex == 1 if sex < .
	*ren famwt98f famwt 
	ren famsze98 famsze 
	ren educyr98 educyear
*** MEDICAL Expenditures;
	foreach item in exp slf mcr mcd prv va ofd wcp opr opu osr{
		ren tot`item'98 meps`item'
	}
	gen yr = 1998

*** Utilization variables ; 
	foreach item in obtotv ipdis ipngtd {  
		ren `item'98 `item' 
	}
	
	ren obtotv doctim 
	ren ipdis  hsptim
	ren ipngtd hspnit
	
***Total income
        ren ttlp98x  totinc		
	
***Insurance variables;
*	foreach item in triev mcrev mcdev opaev opbev unins inscov {  
*		ren `item'98 `item' 
*	}
*perwt  no found in 1998
*	diabdx53 hibpdx53 chddx53 angidx53 midx53 ohrtdx53 strkdx53 emphdx53 asthdx53 bmindx53 adsmok42

	#d;
	keep dupersid perwt yr male educyear age racex hispanx marry53x
	meps* doctim hsptim hspnit 
	duid famidyr famrfpyr famszeyr panel98
	iadlhp53 adlhlp53 
	region53 
	rthlth53 mnhlth53
	totinc
	empst* selfcm* numemp* insc* held* offer* ;
	* just the last years indcat* occc*; 
	#d cr
	
    save "$mepsout/m1998.dta", replace		

*****1999

fdause "$mepsdir/h38.ssp"


	* Rename variables
	ren inscop99 inscopend
	ren age99x age
	ren perwt99f perwt
	gen male = sex == 1 if sex < .
	ren famwt99f famwt 
	ren famsze99 famsze 

*** MEDICAL Expenditures;
	foreach item in exp slf mcr mcd prv va ofd wcp opr opu osr{
		ren tot`item'99 meps`item'
	}
	gen yr = 1999

*** Utilization variables ; 
	foreach item in obtotv ipdis ipngtd {  
		ren `item'99 `item' 
	}
	
	ren obtotv doctim 
	ren ipdis  hsptim
	ren ipngtd hspnit
	
***Total income
        ren ttlp99x  totinc	
		
***Insurance variables;
*	foreach item in triev mcrev mcdev opaev opbev unins inscov {  
*		ren `item'99 `item' 
*	}
*perwt  no found in 1998
*	diabdx53 hibpdx53 chddx53 angidx53 midx53 ohrtdx53 strkdx53 emphdx53 asthdx53 bmindx53 adsmok42

	#d;
	keep dupersid perwt yr male educyear age racex hispanx marry53x
	meps* doctim hsptim hspnit 
	duid famidyr famrfpyr famszeyr panel99
	iadlhp53 adlhlp53 
	region53 
	rthlth53 mnhlth53
	totinc
	empst* selfcm* numemp* insc* held* offer* ;
	* just the last years indcat* occc*; 
	#d cr
	
    save "$mepsout/m1999.dta", replace		

	
	
	
*******2000-2008
forvalues frsyr = 2000/2008  {
	
	local fnumid = `frsyr' - 2000 + 1
	local v =  word("$yrlist", `fnumid')
	local fnum = word("$csdlist", `fnumid')
	local i = `frsyr' - 2000
	local pan = "0`i'"
	drop _all
	fdause "$mepsdir/h`fnum'.ssp"
	if `frsyr' == 1996 {
		ren inscope1 inscop31
		ren inscope2 inscop42
		ren marry96x marry53x
		ren iadlhlp1 iadlhp53 
		ren adlhelp1 adlhlp53
	}
	
	if `frsyr' < 1999{
		ren wtdper`v' perwt`v'f
		ren educyr`v' educyear
	}
	
	if `frsyr' == 2000 {
		gen bmindx53 = (weight53 * 0.4545)/(hghtft53*0.3048+hghtin53*0.0254)^2
		replace bmindx53 = . if weight53 <0 |hghtft53 < 0 | hghtin53  < 0 
	}
	
	if `frsyr' > 2004 {
		ren educyr educyear
	}
	
	if `frsyr' > 2006 {
		foreach oldvar in diabdx hibpdx chddx angidx midx ohrtdx strkdx emphdx asthdx  { 
			ren `oldvar' `oldvar'53
		}
	}
	* Only keep those with positive weights and in scope
	* keep if perwt`v'f > 0 & perwt`v'f < . & insc1231 == 1
	* keep if perwt`v'f > 0 & perwt`v'f < .
	
	* Rename variables
	ren inscop`v' inscopend
	ren age`v'x age
	ren perwt`v'f perwt
	gen male = sex == 1 if sex < .
	*ren famwt`v'f famwt 
	ren famsze`v' famsze
	ren ttlp`v'  totinc
	
*** MEDICAL Expenditures;
	foreach item in exp slf mcr mcd prv va ofd wcp opr opu osr{
		ren tot`item'`v' meps`item'
	}
	gen yr = `frsyr'

*** Utilization variables ; 
	foreach item in obtotv ipdis ipngtd {  
		ren `item'`v' `item' 
	}
	
	ren obtotv doctim 
	ren ipdis  hsptim
	ren ipngtd hspnit
	
***Insurance variables;
*	foreach item in triev mcrev mcdev opaev opbev unins inscov {  
*		ren `item'`v' `item' 
*	}

	
	
*no famwt 96 	
	if (`frsyr' >= 2005) {
		gen panel`pan' = panel
	}	
	#d;
	keep dupersid yr male educyear age racex hispanx marry53x
	meps* doctim hsptim hspnit perwt 
	duid famidyr famrfpyr famszeyr 
	iadlhp53 adlhlp53 diabdx53 hibpdx53 chddx53 angidx53 midx53 ohrtdx53 strkdx53 emphdx53 asthdx53 bmindx53 adsmok42
	region53 
	rthlth53 mnhlth53 panel`pan'
	totinc
	empst* selfcm* numemp* insc* held* offer* ;
	* just the last years indcat* occc*; 
	#d cr
	save "$mepsout/m`frsyr'.dta", replace
}

* Combine mutiple years of cost data

drop _all
set obs 1
gen t = 1
forvalues i = 2008(-1)1996{
	append using "$mepsout/m`i'.dta"
	erase "$mepsout/m`i'.dta"
}
drop if t == 1
drop t
save "$mepsout/csd", replace


*****************************************************
****** From 1996 to 2008
******SELF REPORTED CONDITIONS
*****************************************************


clear
set more off
clear matrix
set seed 2344234
set trace off



*LIST OF FILE NUMBERS FOR ANNUAL CONSOLIDATED FILE 1996 - 2008
global csdlist  12 20 28 38 50 60	70	79	89	97	105	113	121
*global csdlist   50 60	70	79	89	97	105	113	121
*LIST OF FILE NUMBERS FOR MEDICAL CONDITIONS FILE 2000 = 2008
*Note: I changed the name hc006r to h6 for the year 1996
global condlist 6 18 27 37 52	61	69	78	87	96	104	112	120
*global condlist 52	61	69	78	87	96	104	112	120
*LIST OF TWO-DIGIT YEARS - NEEDED TO RENAME VARIABLES
global yrlist 96 97 98 99 00 01 02 03 04 05 06 07 08 
*global yrlist 00 01 02 03 04 05 06 07 08 




*******************************************
* Selected self-reported conditions from conditions file
*******************************************
#d;
	use "$mepsout/csd", clear ;
	sort dupersid yr, stable ; 
	tempfile old ; 
	save `old', replace ; 

set trace off ;
forvalues i = 1996/2008 { ;
	dis "year `i'" ;
	local fnumid = `i' - 1996 + 1 ;
	local fnum = word("$condlist", `fnumid');
	
	fdause "$mepsdir/h`fnum'.ssp", clear ; 
	keep dupersid cccodex ;
	* Destring cccodex ; 
	destring cccodex, gen(ncccodex)  ; 
	gen cancrecr = 1 if inrange(ncccodex, 11, 21) | inrange(ncccodex, 24,45) ; 
	gen heartecr = 1 if inrange(ncccodex, 96,97) | inrange(ncccodex, 100,108) ; 
	* gen lungecr  = 1 if inrange(ncccodex, 127,127) | inrange(ncccodex, 129,134) ; 
	gen lungecr  = 1 if inlist(ncccodex, 127,129,130,131,132) ; 	
	gen diabecr  = 1 if inlist(ncccodex, 49,50) ;
	gen hibpecr  = 1 if inlist(ncccodex, 98,99) ; 
	gen strokecr = 1 if inlist(ncccodex, 109,110,112,113) ; 
	gen kdny1cr = 1 if inlist(ncccodex,156,158,161) ;
	gen kdny2cr = 1 if inlist(ncccodex,156,158,160,161);
	
	* One observation per person ; 
	sort dupersid, stable; 
	foreach v in diabe hibpe stroke hearte lunge cancre kdny1 kdny2 { ; 
		by dupersid: egen cum = total(`v'cr == 1) ; 
		by dupersid: replace `v' = cum  >= 1 ; 
		drop cum ; 
	} ; 
	by dupersid: keep if _n == 1 ; 
	
	gen yr = `i' ; 
	merge 1:1 dupersid yr using `old' ; 
	tab _merge ; 
	qui count if _merge == 1 ; 
	if r(N) > 0 { ; 
		dis "NO matched file in CSD: " r(N); 
	}; 
	drop if _merge == 1 ; 
	drop _merge; 
	keep if yr == `i' ;
	
	foreach v in diabe hibpe stroke hearte lunge cancre kdny1 kdny2 { ; 
		replace `v'cr = 0 if `v'cr!= 1 ; 
	};
	
	drop *codex ;
	sort dupersid yr, stable ; 
	save "$mepsout/c`i'", replace ; 

} ; 
#d cr

drop _all
set obs 1
gen t = 1
forvalues i = 2008(-1)1996{
	append using "$mepsout/c`i'.dta"
	erase "$mepsout/c`i'.dta"
}
drop if t == 1
drop t
	
	label var diabecr "Diabetes CCC code 049/050" 
	label var hibpecr "hypertension CCC code 098/099" 
	label var strokecr "Stroke CCC code 109/110/112/113"  
	label var lungecr "Lung disease CCC code 127/129-134"  
	label var heartecr "Heart disease CCC code 96/97/100 to 108"  
	label var cancrecr "Cancer (except skin) CCC code 11-21/24-45"  
	label var kdny1cr "CKD CCC code 156/158/161" 	
	label var kdny2cr "CKD CCC code 156/158/160/161" 
	 		
		
	******************************  
	* Recode variables 	        
	******************************  
	
	* Recode health conditions in MEPS
	gen hearte = 0 
	label var hearte "CHD/ANGINA/MI/OTHER heart problems"
	foreach v in chd angi mi ohrt { 
		replace hearte = 1 if `v'dx53 == 1
	}
	foreach v in chd angi mi ohrt { 
		replace hearte = . if `v'dx53 < 0 & hearte == 0
	}
	
	gen diabe = diabdx53 == 1 if inlist(diabdx53,1,2)
	label var diabe "Ever diagnosed with diabetes"
	gen hibpe = hibpdx53 == 1 if inlist(hibpdx53,1,2)
	label var hibpe "Ever diagnosed with high blood pressure"
	gen stroke = strkdx53 == 1 if inlist(strkdx53,1,2)
	label var stroke "Ever diagnosed with stroke"
	gen lunge   = emphdx53 == 1 if inlist(emphdx53,1,2)
	label var lunge "Ever diagnosed with emphysema"
	
	gen adl1p = adlhlp53 == 1 if inlist(adlhlp53,1,2)
	label var adl1p "One or more ADLs"

	gen iadl1 = iadlhp53 == 1 if inlist(iadlhp53,1,2)
	replace iadl1 = 0 if adl1p == 1 
	label var iadl1 "IADL only"
		
	************
	* Recode demographic variables 
	************
	
	gen age5659 = inrange(age,56,59) if age < .
	label var age5659 "Aged 56-59"
	gen age6064 = inrange(age,60,64) if age < . 
	label var age6064 "Aged 60-64"
	
	ren hispanx hispan
	replace hispan = 0 if hispan == 2
	gen black = racex == 2 & hispan == 0 if inrange(racex,1,6)
	label var black "R is non-hispanic black"
	
	*** Note: in year 2000 and 2001, coding of racex is 1 to 5, and black is 4
	replace black = racex == 4 & hispan == 0 if inrange(racex,1,5) & inlist(yr,2000,2001)
	
	#d;
	recode educyear (0/11 = 1 "1 less than HS") (12 = 2 "2 HS grad")
	(13/15 = 3 "3 Some college") (16/17 = 4 "4 College grad") (nonmissing = .),gen(educ) ; 
	label var educ "Education recoded";
	#d cr
	
	gen hsdrop = educ == 1 if educyear < . 
	label var hsdrop "Less than HS"
	gen somecol = educ == 3 if educyear < . 
	label var somecol "Some college"
	gen collgrad = educ == 4 if educyear < . 
	label var collgrad "College grad"
	
	gen regnth = region53 == 1 if region53>0 & region53<.
	label var regnth "Census region: Northeast"
	
	gen regmid = region53 == 2 if region53>0 & region53<. 
	label var regmid "Census region: Midwest"
	
	gen regwst = region53 == 4 if region53>0 & region53<. 
	label var regwst "Census region: West"
	
	gen widowed = inlist(marry53,2,8) if marry53>0&marry53<.&marry53!=6
	label var widowed "Marital status:widowed"
	
	gen single = inlist(marry53,3,4,5,9,10) if marry53>0&marry53<.&marry53!=6
	label var single "Marital status: single"
	
	************
	* Recode smoking and obesity 
	************

	
	gen bmi = bmindx53 if bmindx53>0&bmindx53<50
	label var bmi "BMI if 0-50"
	*** obesity
	gen obese = 0*bmi
	replace obese = 1 if bmi >=30 & bmi < .
	label var obese "whether obese (bmi>=30)"
	
	*** overweight
	gen overwt = 0*bmi
	replace overwt = 1 if bmi >= 25 & bmi < 30
	label var overwt "whether over weight (25<=bmi<30)"
	
	*** normal weight
	gen normalwt = 0*bmi
	replace normalwt = 1 if bmi >= 18.5 & bmi < 25
	label var normalwt "whether normal weight (20<=bmi<25)"
	
	*** underweight
	generate underwt = 0*bmi
	replace underwt = 1 if bmi < 18.5 & bmi > 0
	label var underwt "whether under-weight (bmi< 20)"
	
	*** exclusive weight status
	gen wtstate = 0*bmi
	replace wtstate = 1 if underwt  == 1
	replace wtstate = 2 if normalwt == 1
	replace wtstate = 3 if overwt == 1
	replace wtstate = 4 if obese == 1
	label var wtstate "bmi status"
	
	*** panel id
	gen panel = 1 if panel96==1|panel97==1
	replace panel = 2 if panel97==2|panel98==2
	replace panel = 3 if panel98==3|panel99==3
	replace panel = 4 if panel99==4|panel00==4
	replace panel = 5 if panel00==5|panel01==5
	replace panel = 6 if panel01==6|panel02==6
	replace panel = 7 if panel02==7|panel03==7
	replace panel = 8 if panel03==8|panel04==8
	replace panel = 9 if panel04==9|panel05==9
	replace panel = 10 if panel05==10|panel06==10
	replace panel = 11 if panel06==11|panel07==11
	replace panel = 12 if panel07==12|panel08==12
	replace panel = 13 if panel08==13
	
	************
	* Save the file 
	************	
	
	************
	* adjust all monetary amouts
	************
	gen year = yr
	recode totinc  (500000/max = .) (min/0 = 0)
	
	gen cpi = 156.9 if year==1996
	replace cpi = 160.5 if year==1997
	replace cpi = 163.0 if year==1998
	replace cpi = 163.0 if year==1998
	replace cpi = 166.6 if year==1999
	replace cpi = 172.2 if year==2000
	replace cpi = 177.1 if year==2001
	replace cpi = 179.9 if year==2002
	replace cpi = 184 if year==2003
	replace cpi = 188.9 if year==2004
	replace cpi = 195.3 if year==2005
	replace cpi = 201.6 if year==2006
	replace cpi = 207.3 if year==2007
	replace cpi = 215.3 if year==2008
	global cpi2005 = 195.3

	* NHE (corrected to substract NH spending, see historical.ipynb)
	gen nhe = 3130.258303 if year==1996
	replace nhe = 3266.978102 if year==1997
	replace nhe = 3417.046931 if year==1998
	replace nhe = 3602.609319 if year==1999
	replace nhe = 3817.290780 if year==2000
	replace nhe = 4107.445614 if year==2001
	replace nhe = 4434.229965 if year==2002
	replace nhe = 4748.579310 if year==2003
	replace nhe = 5059.372014 if year==2004
	replace nhe = 5370.477966 if year==2005
	replace nhe = 5665.624161 if year==2006
	replace nhe = 5958.471761 if year==2007
	replace nhe = 6184.875000 if year==2008
	
	* correct NHE to exclude NH spending, 6.3% share from NHE, adjust to 2005 dollars
	gen nhe_real = nhe * ($cpi2005/cpi)
	* adjust other amounts to 2005 prices
	gen mepsexp_real = mepsexp * $cpi2005/cpi
	gen totinc_real = totinc * $cpi2005/cpi
	foreach x in slf mcr mcd prv va ofd wcp opr opu osr {
		gen meps`x'_real = meps`x' * $cpi2005/cpi
	}

	gen mepsexp_nhe = .
	forvalues j = 1996/2008 {
		sum mepsexp_real [aw=perwt] if year==`j'
		global totmeps = r(mean)
		sum nhe_real [aw=perwt] if year==`j'
		global totnhe = r(mean)
		replace mepsexp_real = mepsexp_real * ($totnhe/$totmeps) if year==`j'
		replace mepsslf_real = mepsslf_real * ($totnhe/$totmeps) if year==`j'
	}	
	
	drop mepsexp mepsslf
	rename mepsexp_real mepsexp
	label var mepsexp "MEPS total spending corrected to match NHE, $2005"
	rename mepsslf_real mepsslf
	label var mepsslf "MEPS oop spending corrected to match NHE, $2005"
	
	* insurance coverage
	gen inscov = .
	replace inscov = inscov96 if year==1996
	replace inscov = inscov97 if year==1997
	replace inscov = inscov98 if year==1998
	replace inscov = inscov99 if year==1999
	replace inscov = inscov00 if year==2000
	replace inscov = inscov01 if year==2001
	replace inscov = inscov02 if year==2002
	replace inscov = inscov03 if year==2003
	replace inscov = inscov04 if year==2004
	replace inscov = inscov05 if year==2005
	replace inscov = inscov06 if year==2006
	replace inscov = inscov07 if year==2007
	replace inscov = inscov08 if year==2008

	* keep males
	keep if male==1
	* age range
	keep if age>=25&age<=84
	* adjust total mepsexp for outliers	
	recode mepsexp (250000/max=.)
	recode mepsslf (250000/max=.)
	
	* recode health 
	* 1 = dead, 2 = poor-fair, 3 = good, 4 = vgood-excellent
	recode rthlth53 (1/2=4) (3=3) (4/5=2) (else=.), gen(hlth)

	label data "MEPS1996-2008" 
	saveold "$mepsout/MEPS9608", replace 
	erase "$mepsout/csd.dta"
		

	
	
	

*****************************************************
* Estimation of production function using MEPS-NHIS *
*****************************************************
graph drop _all
* extract from MEPS-NHIS data
use $src/clean/meps-nhis_wide_extract.dta, clear

* first stage variable (do log(1+x) to deal with zeros)
gen logmed  = log(1+fmepsexp)
gen llogmed = log(1+mepsexp)
gen loginc  = log(1+totinc)
gen loginc2 = loginc*loginc
* instrument set and regressors
global x "age ib2.hlth obese smokev"
global z "loginc"
* first stage
xi: reg logmed $x $z, robust
test $z
predict ulogmed, residuals
gen ulogmed2 = ulogmed^2
gen ulogmed3 = ulogmed^3
* second-stage (mlogit with control function), can do robust by adding those ulogmed2 ulogmed3
gen logmed2 = logmed*logmed
mlogit fhlth ulogmed logmed logmed2 $x if fhlth!=1, base(2)

matrix b = e(b)
matrix bhealth = J(10,3,.)
local i = 1
forvalues h = 1/3 {
	forvalues x=1/10 {
		matrix bhealth[`x',`h'] = b[1,`i']
		local i = `i'+1
	}
}
matrix bhealth = bhealth[2..10,1..3]
matrix list bhealth

margins , dydx(logmed) predict(outcome(2))
margins , dydx(logmed) predict(outcome(3))
margins , dydx(logmed) predict(outcome(4))


* predictions for transitions at different rates of spending
preserve
gen hlth3 = hlth==3 if hlth!=.
gen hlth4 = hlth==4 if hlth!=.
mlogit fhlth ulogmed logmed logmed2 age hlth3 hlth4 obese smokev if fhlth!=1, base(2)
replace hlth3 = 0
replace hlth4 = 1
replace ulogmed = 0
predict p1_f p2_f p3_f, pr
replace logmed = 0
replace logmed2 = 0
predict p1 p2 p3, pr
tabstat p1_f p2_f p3_f p1 p2 p3 if e(sample)==1, by(age)
restore

* computing corrected standard errors using bootstrap
gen hlth3 = hlth==3 if hlth!=.
gen hlth4 = hlth==4 if hlth!=.

global x "age hlth3 hlth4 obese smokev"
global z "loginc"
capture program drop controlfunc
program controlfunc, eclass
  version 10.1
  tempname b
  reg logmed $x $z, robust
  capture drop ulogmed
  predict ulogmed, residuals
  mlogit fhlth ulogmed logmed logmed2  $x if fhlth!=1, base(2)
  matrix `b' = e(b)
  ereturn post `b' 
end
preserve
	bootstrap _b, reps(100) seed(10101) nodots nowarn: controlfunc
restore

* other robustness checks (for editor)
tab educ, gen(educ_)
global x "age hlth3 hlth4 obese smokev educ_2 educ_3 educ_4"
global z "loginc"
capture program drop controlfunc
program controlfunc, eclass
  version 10.1
  tempname b
  reg logmed $x $z, robust
  capture drop ulogmed
  predict ulogmed, residuals
  mlogit fhlth ulogmed logmed logmed2  $x if fhlth!=1, base(2)
  matrix `b' = e(b)
  ereturn post `b' 
end
preserve
	bootstrap _b, reps(100) seed(10101) nodots nowarn: controlfunc
restore



**** Model for mortality *********
gen dead = fhlth==1 if missing(fhlth)==0
global xd "age ib2.hlth"
cloglog dead $xd 
matrix bdead = e(b)'
predict pdead, pr

global x "age ib2.hlth obese smokev"
* effect medical spending on mortality (Figure 2 in paper)
preserve
	local j = 1
	set more off
	cloglog dead $xd
	margins hlth if age>=65, post
	matrix pd = e(b)'
	matrix list pd
	matrix result = J(13,3,0)
	mlogit fhlth ulogmed logmed logmed2  $x if fhlth!=1, base(2)
		foreach med of numlist 0 10 100 250 500 750 1000 2500 5000 7500 10000 15000 20000 {
			local logmed = log(1+`med')
			local logmed2 = log(1+`med')^2
			
			matrix ph = J(3,3,.)
			di `med'			
			foreach h of numlist 2/4 {
				qui margins hlth if age>=65, at(logmed=`logmed' logmed2=`logmed2') predict(outcome(`h')) 
				local h1 = `h'-1
				matrix ph[1,`h1'] = r(b)'
			}		
			matrix list  ph
			matrix y = ph*pd
			matrix result[`j',1]  = y'
			local j = `j'+1
		}
	matrix list result
	gen med = .
	replace med = 0 in 1
	replace med = 10 in 2
	replace med = 100 in 3
	replace med = 250 in 4
	replace med = 500 in 5
	replace med = 750 in 6
	replace med = 1000 in 7
	replace med = 2500 in 8
	replace med = 5000 in 9
	replace med = 7500 in 10
	replace med = 10000 in 11
	replace med = 15000 in 12
	replace med = 20000 in 13
	svmat result
	keep result* med

	* marginal cost
	forvalues i = 1/3 {
		gen mc`i' = .
		forvalues j = 1/12 {
			replace mc`i' = -(med[`j'+1] - med[`j'])/(result`i'[`j'+1] - result`i'[`j'])*result`i'[`j']*0.1 if _n==`j'
		}
	}

	#d ;
		twoway 
			(line result1 med if med>=0)
			(line result2 med if med>=0)
			(line result3 med if med>=0)
			,
			xtitle("medical expenditures")
			ytitle("")
			title("mean mortality rate")
			subtitle("by current health status and medical expenditures")
			$fig
			ylabel(,labsize(small))
			xlabel(1000 2500 5000 7500 10000 15000 20000, labsize(small))
			legend(label(1 "poor-fair") label(2 "good") label(3 "vgood-excellent")  
			size(vsmall) rows(1));
			
	#d cr
	graph export figures/mortality-effect.eps, as(eps) replace
	forvalues j = 1/3 {
		replace result`j' = 1-result`j'
	}
	#d ;
		twoway 
			(line mc1 result1 if med>=0)
			(line mc2 result2 if med>=0)
			(line mc3 result3 if med>=0)
			,
			xtitle("survival probability")
			ytitle("")
			title("marginal cost of reducing mortality by 10%")
			subtitle("by current health status and survival risk")
			$fig
			ylabel(,labsize(small))
			xlabel(, labsize(small))
			legend(label(1 "poor-fair") label(2 "good") label(3 "vgood-excellent")  
			size(vsmall) rows(1));
			
	#d cr
	graph export figures/mcost-effect.eps, as(eps) replace

restore



* save productivity parameters
preserve
	svmat bhealth
	keep bhealth*
	keep if _n<=2
	outsheet using params/input/production-theta.csv, replace comma nonames nolabel	
restore	


* save parameters other than productivity
preserve
	svmat bhealth
	keep bhealth*
	keep if bhealth1!=.
	drop if _n<=2
	outsheet using params/input/production-health.csv, replace comma nonames nolabel	
restore	

* save death coefficients
preserve
	svmat bdead
	keep bdead*
	keep if bdead1!=.
	outsheet using params/input/production-death.csv, replace comma nonames nolabel	
restore	

merge n:1 age using $src/other/mortality-ssa-cohort.dta
drop if _merge==2
egen mpdead = mean(pdead), by(age)
gen pi = (pmx/mpdead)
tabstat pmx mpdead pi, by(age)

preserve
* estimate a log linear model with spline at 60, good fit for smooth and extrapo
keep pi age dead mx pmx
collapse pi dead mx pmx, by(age)
gen logpi = log(pi)
gen age60m = min(age,60)
gen age60p = max(age-60,0)
reg logpi age60m age60p
matrix pimat = e(b)'
predict plogpi
gen ppi = exp(plogpi)
tabstat ppi, by(age)
gen checkpmx = dead*ppi
* graph of mortality compared to period life-table in MEPS (Figure B.1)
#d ;
twoway (line mx age if age<84, sort lpattern(solid)) (line pmx age if age<84,sort lpattern(dash))
(line dead age if age<84,sort lpattern(dot)) 
	, 
	xtitle("age of respondent")
	ytitle("") $fig
	xlabel(25(5)85, labsize(small))
	title("one-year mortality rate")
	legend(label(1 "period mx") label(2 "cohort mx") label(3 "MEPS-NHIS mx"));
#d cr
graph export figures/mortality-comparison.eps, as(eps) replace

* intrapolation of adjustment factor
twoway (line ppi age, sort) (line pi age, sort)
	svmat pimat
	keep pimat
	outsheet using params/input/adjustdead.csv, replace comma nonames nolabel
restore

* data on obesity and smoking
preserve
graph drop _all
	gen age2 = age*age
	char byr_c[omit] 4
	xi: reg obese age age2
	matrix b = e(b)
	matrix bcons = b[1,"_cons"]
	matrix bage = b[1,"age"]
	matrix bage2 = b[1,"age2"]
	gen pobese = bage[1,1]*age + bage2[1,1]*age2 + bcons[1,1] 
	sum pobese if age==75
	replace pobese = r(mean) if age>=75
	xi: reg smokev age age2 
	matrix b = e(b)
	matrix bcons = b[1,"_cons"]
	matrix bage = b[1,"age"]
	matrix bage2 = b[1,"age2"]
	gen psmoke = bage[1,1]*age + bage2[1,1]*age2 + bcons[1,1] 
	sum psmoke if age==75
	replace psmoke = r(mean) if age>=75
	collapse psmoke pobese smokev obese, by(age) 	
	#d ;
	twoway (scatter obese age ) (line pobese age ), name(obese) $fig;	
	graph export figures/risk-obese.eps, as(eps) replace;
	twoway (scatter smokev age ) (line psmoke age), name(smoke) $fig;
	graph export figures/risk-smoke.eps, as(eps) replace;
	#d cr		
	local samp = 120-25+1 
	set obs `samp'
	replace age = age[_n-1]+1 if age==.
	replace psmoke = psmoke[_n-1] if psmoke==.
	replace pobese = pobese[_n-1] if pobese==.
	keep psmoke pobese
	order psmoke pobese
	outsheet using params/input/risk.csv, replace comma nonames nolabel
restore

exit




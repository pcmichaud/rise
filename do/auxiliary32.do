***********************************************************
* Estimation of Auxiliary Processes
* creates parameter files:
* - params/input/earnings.csv
* - params/input/otherinc.csv
* - params/input/initsample.csv
***********************************************************

***************************************************
* Step 1: Earnings process PSID                   *
***************************************************

* minimum distance function for (b = (rho, sige))		
mata 
  void mdfunc(todo, b, mCovData, func, g, H)
  {
  	// map parameters
  	T = cols(mCovData)
  	rho= tanh(b[1,1])
  	sige = exp(b[1,2])
	sigv = exp(b[1,3])
  	// Get "True" Covariance Matrix
  	mCovTrue = J(T,T,.)
  	for (r = 1; r<=T; ++r) {
  		for (c=1; c<=T; ++c) {
  			mCovTrue[r,c] = (r==c)*sigv + (rho^abs(r-c))*sige/(1-rho^2)

  		}
  	} 
	d = vech(mCovData :- mCovTrue) 	  	
  	func = d'*d
  }
end

preserve
		use $src/clean/psid_extract.dta, clear
		* sample selection
		*keep if rwork==1
		keep if logriearn!=.
		keep if rage>=25&rage<=65
		keep if year<=1997
		gen rage2 = rage^2
		sort hhidpn year
		* declare as panel
		bysort hhidpn: gen t = _n
		bysort hhidpn: gen T = _N 
		bysort hhidpn: gen last = t==T
		tsset hhidpn t
		xtdes
		tab T
		table rage rshlt3, content(mean riearn)
		keep if T>=3
		* deterministic component of earnings	
		xtreg logriearn rage rage2, robust fe
		matrix b = e(b)
		predict ue, ue	
		predict e, e
		* get average fixed effect for those born 1935-45
		predict u, u
		qui sum u if rbyr_c==3
		local mu = r(mean)
		matrix b[1,3] = b[1,3] + `mu'		
		* smearing for retransformation
			gen plogriearn = b[1,3] + b[1,1]*rage + b[1,2]*rage2
			gen priearn = exp(plogriearn) 
			sum riearn [aw=wgid_cs]
			local mriearn = r(mean)
			sum priearn [aw=wgid_cs] 
			local mpriearn = r(mean)
			local smear = log(`mriearn'/`mpriearn')
			matrix b[1,3] = b[1,3] + `smear'
			replace priearn = exp(`smear')*priearn
			matrix list b	
			sum priearn riearn
		* collect residuals	for minimum distance (estimation of covariance)
		sum ue
		di _N
		forvalues t = 1/12 {
			by hhidpn: gen ue_`t' = ue[`t'] 
		}
		keep if last==1
		di _N
		corr ue_1-ue_12, cov
		matrix mCov = r(C)
		* optimization for estimating covariance
		mata {
			// covariance matrix from data
			mCovData = st_matrix("mCov")
			D = optimize_init()
			optimize_init_which(D , "min")
			optimize_init_evaluator(D, &mdfunc())
			optimize_init_evaluatortype(D, "d0")
			optimize_init_technique(D, "bfgs")
			optimize_init_argument(D, 1, mCovData)	
			// Starting values
			b0 = (0.9,log(0.03),log(0.08))
			optimize_init_params(D,b0)
			// Optimize
			b = optimize(D)
			Vb = optimize_result_V(D)
			//Retrieve parameters
			b[1,1] = tanh(b[1,1])
			b[1,2] = exp(b[1,2])
			b[1,3] = exp(b[1,3])
			st_matrix("cb",b)

		}
	
		* prepare for posting
		matrix b = (b,cb)'
		matrix rownames b = age age2 constant rho sige sigv
	restore

di "parameters of earnings process"

matrix list b


* prepare graph of average earnings profiles with 95% band
preserve
	graph drop _all
	use $src/clean/psid_extract.dta, clear
	keep if logriearn!=.
	keep if rage>=25&rage<=65
	keep if year<=1997
	gen priearn = exp(b[1,1]*rage + b[2,1]*rage^2 + b[3,1]) 
	scalar varu = b[5,1]/(1-b[4,1]^2)
	di "stationary variance is ", varu
	gen up_priearn = priearn*exp(sqrt(varu)*2) 
	gen low_priearn = priearn*exp(-sqrt(varu)*2) 
	sum priearn low_* up_*
	
	collapse (mean) priearn up_priearn low_priearn, by(rage)
	replace priearn = priearn * 1e-3
	replace up_priearn = up_priearn * 1e-3
	replace low_priearn = low_priearn * 1e-3
	keep if rage<=70
	#d ;
	twoway  (line up_priearn rage, lpattern(dash))
	 	   (line priearn rage, lpattern(solid))
	 	   (line low_priearn rage, lpattern(dash)),
	 	   legend(label(1 "+2se") label(2 "mean") label(3 "-2ses") rows(1) size(small))
	 	   xtitle("age of household head") $fig
	 	   ytitle("earnings (000)")
	 	   xlabel(25(5)70)
	 	   ylabel(0(25)200, labsize(small))
	 	   name(earnings)
	 	   ;
	#d cr

	graph export figures/earnings.eps, as(eps) replace
restore

* transfer parameters to ascii to be read in fortran
svmat b
keep b1
keep if b1!=.
outsheet using params/input/earnings.csv, replace comma nonames nolabel


**************************************
* Step 2: Spouse earnings + pensions *
* regression by age, education	     *
* and earnings       (PSID)          *
**************************************
preserve
	graph drop _all
	use $src/clean/psid_extract.dta, clear
	keep if rage>=25 & rage<=84
	gen rage2 = rage*rage
	* own income
	gen own = hissben + riearn	
	replace own = . if own>=250000
	replace otherinc = . if otherinc>=250000
	gen prbyr_c = rbyr_c
	* instrument
	ivregress 2sls otherinc (own = rskill) rage rage2 ib3.prbyr_c
	replace prbyr_c = 3
	matrix beta = e(b)
	matrix beta = (beta[1,1],beta[1,2], beta[1,3], beta[1,10])'	
		gen potherinc = beta[4,1] + beta[1,1]*own + beta[2,1]*rage + beta[3,1]*rage2
		collapse (mean) potherinc (p5) potherinc_5=potherinc (p95) potherinc_95=potherinc, by(rage)
		replace  potherinc_5 = max(potherinc_5,0)
		replace potherinc_95 = potherinc_95 * 1.e-3
		replace potherinc = potherinc * 1.e-3
		replace potherinc_5 = potherinc_5 * 1.e-3
		#d ;
		twoway  (line potherinc_95 rage if rage<=85, lpattern(dash))
	 	   (line potherinc rage if rage<=85, lpattern(solid))
	 	   (line potherinc_5 rage if rage<=85, lpattern(dash)),
	 	   legend(label(1 "95th pct") label(2 "mean") label(3 "5 pct") rows(1) size(small))
	 	   xtitle("age of household head") $fig
	 	   xlabel(25(5)85)
	 	   ytitle("other income (000)") 
	 	   ylabel(0(10)50, labsize(small))
	 	   name(otherinc)
	 	   ;
		#d cr
	graph export figures/otherincome.eps, as(eps) replace
	svmat beta
	keep if beta1!=.
	keep beta1
	list
	outsheet using params/input/otherinc.csv, replace comma nonames nolabel
restore

**************************************
* Step 3: fraction with spouse *
* regression by age, education	     *
* and earnings       (PSID)          *
**************************************
preserve
	graph drop _all
	use $src/clean/psid_extract.dta, clear
	keep if rage>=25 & rage<=100
	* only use data prior to age 55, will impute from reg to older ages
	gen rage2 = rage*rage
	gen rage3 = rage*rage*rage
	gen rage4 = rage*rage*rage*rage
	gen rage5 = rage4*rage
	gen prbyr_c = rbyr_c
	* regression
	regress spouse rage rage2 rage3 rage4 rage5 ib3.prbyr_c
	replace prbyr_c = 3
	matrix s = e(b)
	matrix beta = (s[1,1],s[1,2], s[1,3], s[1,4], s[1,5], s[1,12])
		#d ;
		gen pspouse = beta[1,6] + beta[1,1]*rage + beta[1,2]*rage2 + beta[1,3]*rage3 
			+ beta[1,4]*rage4 + beta[1,5]*rage5;
		#d cr
		collapse (mean) pspouse, by(rage)
		#d ;
		twoway  
	 	   (line pspouse rage if rage<=100, lpattern(solid)),
	 	   legend(label(1 "mean") rows(1) size(small))
	 	   xtitle("age of household head")
	 	   xlabel(25(5)85) $fig
	 	   title("fraction with spouse") 
	 	   ylabel(0(0.1)1, labsize(small))
	 	   name(spouse)
	 	   ;
		#d cr
	graph export figures/spouse.eps, as(eps) replace
	matrix beta = beta'
	svmat beta
	keep if beta1!=.
	keep beta1
	list
	outsheet using params/input/spouse.csv, replace comma nonames nolabel
restore

**************************************
* Step 3: Insurance status 
* in HRS (for cohorts borh 1935-1945)
* assume coverage constant with age
**************************************
preserve
	use $src/hrs/rndhrsg_all.dta, clear
	*borh 1935-1945
	keep if rabyear>=1935&rabyear<=1945
	* 50 to 55 years old
	keep if ragey_b>=50&ragey_b<=55
	* keep males
	keep if ragender==1
	
	* get insurance status	
	egen insured = anymatch(rcovr rcovs), values(1)
	replace insured = 2 if insured==1&rcovrt==1
	label def insured 0 "no coverage" 1 "tied coverage" 2 "retiree coverage"
	label values insured insured
	replace insured = . if rhigov==1
	
	* get DB pension
	egen dbpen = anymatch(rptyp*), values(1)
	
	gen rskill=1 if (raeduc==4|raeduc==5)
	replace rskill=0 if (raeduc==1|raeduc==2|raeduc==3)

	* how many have db (for referee)
	tabstat dbpen [aw=rwtresp], by(insured)
	
	* obtaining the fractions of each type
	tab insured [aw=rwtresp], nofreq matcell(data)
	matrix list data
	scalar total = data[1,1]+data[2,1]+data[3,1]
	matrix data[1,1] = data[1,1]/total
	matrix data[2,1] = data[2,1]/total
	matrix data[3,1] = data[3,1]/total
		
	matrix data[2,1] = data[2,1] + data[1,1]
	matrix data[3,1] = 1
	
	* data for using in initial conditions
	matrix rownames data = no tied retiree
	matrix list data
restore


**********************************************
* Step 5: Initial Distribution of 
* (health,educ,assets,earn,
* aime,insurance), PSID
**********************************************
preserve
	use $src/clean/psid_extract.dta, clear
	* use age 24 to 28 to boost initial sample
	keep if rage>=24&rage<=28
	* variables
	keep logriearn ass rshlt3
	* rescale rshlt
	replace rshlt3 = rshlt3 + 1	
	* make sure none missing
	reg logriearn ass rshlt3
	gen miss = e(sample)==0
	drop if miss==1
	* dummy education var 
	gen educ = 1	
	
	* assets (only positive)
	replace ass = 0 if ass<=0
	
	tabstat ass, by(rshlt3) statistics(p10 p25 p50 p75 p90 mean)
	
	* draw insurance status using HRS computations
	gen u = uniform()
	gen insured = .
	replace insured = 1 if u<=data[1,1]&insured==.
	replace insured = 2 if u>data[1,1]&u<=data[2,1]&insured==.
	replace insured = 3 if u>data[2,1]&u<=data[3,1]&insured==.
	drop u miss
	tab insured	
	tab insured if insured!=1
	* replicate twice and sample with replacement 10000 for initial sample	
	expand 10
	bsample 10000
	* order for reading in fortran
	tabstat logriearn, by(rshlt3)
	order logriearn educ insured ass rshlt3
	* save file
	save $src/clean/initsample.dta, replace
	* transfer dataset to ascii (to be read in Fortran)
	outsheet using params/input/initsample.csv, replace comma nonames nolabel
restore

**********************************
* Step 6: Insurance co-payments
* from MEPS, Table 1 of paper 
**********************************
preserve
	use $src/clean/meps-nhis_wide_extract.dta, clear
	keep if age>=25&age<=84
	keep if male==1
	replace mepsexp = . if mepsexp>100000
	* how much paid out-of-pocket
	gen copay = mepsslf/mepsexp
	replace copay = 1 if copay>1&copay!=.	
	label def ins 1 "private" 2 "public" 3 "none"
	label values inscov ins
	matrix results = J(4,1,.)
	sum copay [aw=perwt] if age<65&inscov==1, d
	matrix results[1,1] = r(p50)
	sum copay [aw=perwt] if age<65&inscov==2, d
	matrix results[2,1] = r(p50)
	sum copay [aw=perwt] if age<65&inscov==3, d
	matrix results[3,1] = r(p50)
	sum copay [aw=perwt] if age>=65&inscov==2, d
	matrix results[4,1] = r(p50)
	svmat results
	keep results1
	keep if results1!=.
	list
	outsheet using params/input/coinsurance.csv, replace comma nonames nolabel	
restore

exit







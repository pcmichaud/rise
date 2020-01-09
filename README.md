# Accounting for the Rise of Health Spending and Longevity
Fonseca, Michaud, Galama and Kapteyn 
[Journal of European Economic Association](https://academic.oup.com/jeea)

December 2019

## Software needed 

Four pieces of software are needed to reproduce the results in the paper:
* Fortran compiler: we use [gfortran](https://gcc.gnu.org/wiki/GFortran)
* MPI: we use [openmpi](https://www.open-mpi.org/) 
* Python: we use [Anaconda Python distribution](https://www.anaconda.com/distribution/) 
* Stata: [Stata MP v15](https://www.stata.com/) 
  
Other versions and compilers are possible but users will need to make appropriate changes to programs and makefile. The main fortran executables were executed on an HPCC running the scheduler [SLURM](https://slurm.schedmd.com/documentation.html). All other programs can run on a standard machine. Paths need to be set at various places clearly identified at the beginning of each program. 

## Data Access

The following raw micro datasets are used for this project. 

* The Panel Study of Income Dynamics: We use the files provided by [The Cross National Equivalent File project](https://cnef.ehe.osu.edu/). 
* The Medical Expenditure Panel Study: We use the files from this [website](https://www.meps.ahrq.gov/mepsweb/). 
* The National Health Interview Survey Mortality Follow-Up Study: We use the files available at this [website](https://www.cdc.gov/nchs/data-linkage/mortality.htm). 
* The Health and Retirement Study: We have used version G of [RANDHRS](https://www.rand.org/well-being/social-and-behavioral-policy/centers/aging/dataprod/hrs-data.html)
* The National Health Expenditure Accounts: [website](https://www.cms.gov/Research-Statistics-Data-and-Systems/Statistics-Trends-and-Reports/NationalHealthExpendData/index)
* Other small datasets and parameter files described in paper

They are included in the data folder which will appear once unziped (data.zip). But users are responsible for obtaining registration to use these files with the relevant organizations (for example HRS) and the linked mortality data for MEPS-NHIS.

## Code and Data Archive 

This archive has the following directories: 

* [do](do/): contains do (stata program files) to generate all the inputs and some of the figures and tables. The main do-file is master-dataprep.do which runs all the stata scripts. 

* [py](py/): contains all python programs and notebooks that generate also inputs as well as tables and figures in paper. 
  
* data: once unzipped, described above. contains all data used in paper. the folder clean will contain stata datasets created by the programs. 

* [params](params/): contains parameter files that are produced from do and py programs as well as settings to run the model. 
  
* [figures](figures/): empty but upon execution contain the outputs from do and py programs. 

* [src](src/): contains the source fortran files to run the model

* [runtime](runtime/): contains the executables, once make has been run. 
  
* [libs](libs/): contains two libraries used by the fortran programs: [dcdflib](https://person.hst.aau.dk/magnus/pkgsrc-kolga/math/dcdflib.f/), which is a statistical library and [newuoa](https://en.wikipedia.org/wiki/NEWUOA) which contains the optimization routines. 

* [output.zip](output.zip) contains all the simulated datasets and decision rules in the paper. Once unzipped, creates output directory. 

## Compilation of Fortran Code 

To build executables, there is a [Makefile](Makefile) that is also supplied in the root directory of this archive. Please check carefully to adapt to the particular environment in which we run calculations. Running make at the terminal will create the executables estimate and generate. 

## Estimation 

We run estimation in cohort mode. 

* Set [params/commons/solve.info](params/common/solve.info) to nz = 1 and zmin, zmax = 1940. 
*  Set an params/common/init_params_scnname.info where scnname is the scenario name you want to have. Check also that the [params/common/moments.info](params/common/estimation.info) file uses the moments you want in estimation. 
*  There are settings for estimation that can be set in [params/common/estimation.info](params/common/estimation.info)
* executing mpiexec ./estimate scnname will execute your scenario. 
* output will be directed to output/estimation

Depending on starting values you use, this may take a number of days/weeks to run even if parallel computing is used. 

## Simulating from a Scenario 

* You can simulate in cohort or population mode. This is controlled by nz, zmin and zmax in [params/common/solve.info](params/common/solve.info) 
* Set your scenario in params/scenarios/scenario_name.info. Multiple flags are possible. 
* make sure your estimated params have been copied to params/common/esti_params.info. These will be used for simulation. 
* run mpiexec ./generate scnname. Output is found in output/simulation.  

## Tables and Figures 

The following list contains the source program used to create each table and figures. 

* Figure 1: [py/fig1-Historical.ipynb](do/fig1-Historical.ipynb)
* Figure 2: [do/production32.do](do/production32.do)
* Figure 3: [py/fig3-fit.py](py/fig3-fit.py)
* Figure 4: [py/fig4-production.py](py/production.py)

* Table 1: [do/auxiliary32.do](do/auxiliary32.do)
* Table 2: [do/production32.do](do/production32.do)
* Table 3: [do/production32.do](do/production32.do)
* Table 4: [py/tab4-estimates.py](py/tab4-estimates.py)
* Table 5: [py/tab5-price-elasticity.py](py/tab5-price-elasticity.py)
* Table 6: [py/tab6-income-elasticity.py](py/tab6-income-elasticity.py)
* Table 7: [py/tab7-growth-pop.py](py/tab7-growth-pop.py) 
* Table 8: [py/tab8-welfare.py](py/tab8-welfare.py) 
* Table 9: [py/tab9-gradient.py](py/tab9-gradient.py)

Appendix

* Figure A.1: [py/figA1-rules.py](py/figA1-rules.py)
* Figure A.2: [py/fig3-fit.py](py/fig3-fit.py)
* Figure B.2: [do/production32.do](do/production32.do) 
* Figure C.1: [do/auxiliary32.do](do/auxiliary32.do) 
* Figure C.2: [do/auxiliary32.do](do/auxiliary32.do)
* Figure C.2: [do/healthtrans.do](do/healthtrans.do)


* Table B.1: [do/stats-psid32.do](do/stats-psid32.do)
* Table B.2: [do/stats-meps32.do](do/stats-meps32.do) 
* Table C.2: [do/auxiliary32.do](do/auxiliary32.do) 
* Table C.3: [do/auxiliary32.do](do/auxiliary32.do) 

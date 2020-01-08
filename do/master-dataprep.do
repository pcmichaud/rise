***********************************************************
* Master program to prepare data, moments and parameters
***********************************************************

clear all
capture log close
set more off

* setting home directory of project 
cd ~/cedia/Projets/jeea-fonseca

* setup options and data
do do/setup.do

* construct PSID extract
do do/extract-psid32.do

* construct MEPS extract (matched with NHIS for mortality follow-up)
do do/meps_selected9608_32.do
do do/getsmokenhis-32.do
do do/transfermepsnhis32.do
do do/extract-meps32.do

* obtain estimates of production function
do do/production32.do

* obtain auxiliary parameters and initial sample
do do/auxiliary32.do

* obtain estimates for AIME
do do/aime32.do

* obtain moments for plotting
do do/moments32.do

* stats from PSID and MEPS for appendix
do do/stats-meps32.do
do do/stats-psid32.do

* graph of health transitions for appendix
do do/healthtrans.do

* PSID gradient 
do do/gradient.do

exit

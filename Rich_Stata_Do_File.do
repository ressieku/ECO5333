********************************************************************************************
**********This dofile simulates a hypothetical trade integration between the US and the EU**
**************using the GEPPML procedure developed by Anderson, Larch and Yotov*************
********************************************************************************************

* Clear memory and set parameters
clear all
clear matrix
set matsize 8000
set maxvar 30000
set more off
pause on

* Close and create log
capture log close
log using agoa_log, text replace

* Open and use data
use "\Users\richmondessieku\Desktop\FALL2023\Terem Paper\data.csv"

* Generate output of exporting countries
bys exporter: egen y = sum(trade)
* Generate aggregate expenditure of importing countries
bys importer: egen e = sum(trade)

* Select a reference country for normalization process. In our case, Germany
replace importer = "ZZZ" if importer == "DEU"
replace exporter = "ZZZ" if exporter == "DEU"
g e_deu_bln = e if importer == "ZZZ"
egen e_deu = mean(e_deu_bln)

* Assign noc to the number of countries in our sample
qui tab exporter
global noc = r(r)
global noc_1 = $noc - 1

qui egen pairid = group(importer exporter)

* Generate the fixed effects
qui tab exporter, g(exp_fe)
qui tab importer, g(imp_fe) 

**************************************
*Step 1a: Estimate the Baseline Model*
**************************************
ppml trade ln_dist contig comlang colony brdr eu nafta exp_fe* imp_fe1-imp_fe$noc_1, cluster(pairid) nocons iter(30)

* Predict trade flows in the baseline scenario
predict tradehat_bln, mu
* Generate residuals
g resid = trade/tradehat_bln

* Based on the estimated exporter and importer fixed effects, create the actual set of fixed effects
forvalues k = 1 (1) $noc_1 {
qui replace exp_fe`k' = exp_fe`k' * exp(_b[exp_fe`k'])
qui replace imp_fe`k' = imp_fe`k' * exp(_b[imp_fe`k'])
}
* Create the exporter and importer fixed effects for the country of reference (Germany)
qui replace exp_fe$noc = exp_fe$noc * exp(_b[exp_fe$noc])
qui replace imp_fe$noc = imp_fe$noc * exp(0)
* Creat two variables holding each of all estimated exporter and importer fixed effects in the baseline
egen exp_pi_0 = rowtotal(exp_fe1-exp_fe$noc)
egen exp_chi_0 = rowtotal(imp_fe1-imp_fe$noc)

* Generate the baseline bilateral trade costs
g tij_bln = exp(_b[ln_dist]*ln_dist+_b[contig]*contig+_b[comlang]*comlang+_b[colony]*colony+_b[brdr]*brdr+_b[eu]*eu+_b[nafta]*nafta)

*****************************************
*Step 1b: Construct the Baseline Indexes*
*****************************************
scalar sigma = 7
* Generate the Outward Multilateral Resistance (OMR) index
g omr_bln = y*e_deu / exp_pi_0
* Generate the Inward Multilateral Resistance (IMR) index
g imr_bln = e / (exp_chi_0 * e_deu)
* Generate the real GDP index
g rgdp_bln_temp = y/((imr_bln)^(1/(1-sigma))) if importer == exporter
bys exporter: egen rgdp_bln = sum(rgdp_bln_temp)
* Generate the total exports index
g exports_bln = tradehat_bln if importer != exporter
bys exporter: egen tot_exports_bln = sum(exports_bln)

*********************************************
* Step 2: Define the Counterfactual Scenario*
*********************************************
* Introduce a trade integration between the US and the EU
replace eu = 1 if exporter == "USA" & importer == "AUT"
replace eu = 1 if importer == "USA" & exporter == "AUT"
replace eu = 1 if exporter == "USA" & importer == "BEL"
replace eu = 1 if importer == "USA" & exporter == "BEL"
replace eu = 1 if exporter == "USA" & importer == "BGR"
replace eu = 1 if importer == "USA" & exporter == "BGR"
replace eu = 1 if exporter == "USA" & importer == "CYP"
replace eu = 1 if importer == "USA" & exporter == "CYP"
replace eu = 1 if exporter == "USA" & importer == "CZE"
replace eu = 1 if importer == "USA" & exporter == "CZE"
replace eu = 1 if exporter == "USA" & importer == "DNK"
replace eu = 1 if importer == "USA" & exporter == "DNK"
replace eu = 1 if exporter == "USA" & importer == "EST"
replace eu = 1 if importer == "USA" & exporter == "EST"
replace eu = 1 if exporter == "USA" & importer == "FIN"
replace eu = 1 if importer == "USA" & exporter == "FIN"
replace eu = 1 if exporter == "USA" & importer == "FRA"
replace eu = 1 if importer == "USA" & exporter == "FRA"
replace eu = 1 if exporter == "USA" & importer == "IRL"
replace eu = 1 if importer == "USA" & exporter == "IRL"
replace eu = 1 if exporter == "USA" & importer == "ITA"
replace eu = 1 if importer == "USA" & exporter == "ITA"
replace eu = 1 if exporter == "USA" & importer == "LTU"
replace eu = 1 if importer == "USA" & exporter == "LTU"
replace eu = 1 if exporter == "USA" & importer == "LUX"
replace eu = 1 if importer == "USA" & exporter == "LUX"
replace eu = 1 if exporter == "USA" & importer == "LVA"
replace eu = 1 if importer == "USA" & exporter == "LVA"
replace eu = 1 if exporter == "USA" & importer == "MLT"
replace eu = 1 if importer == "USA" & exporter == "MLT"
replace eu = 1 if exporter == "USA" & importer == "NLD"
replace eu = 1 if importer == "USA" & exporter == "NLD"
replace eu = 1 if exporter == "USA" & importer == "POL"
replace eu = 1 if importer == "USA" & exporter == "POL"
replace eu = 1 if exporter == "USA" & importer == "PRT"
replace eu = 1 if importer == "USA" & exporter == "PRT"
replace eu = 1 if exporter == "USA" & importer == "ESP"
replace eu = 1 if importer == "USA" & exporter == "ESP"
replace eu = 1 if exporter == "USA" & importer == "SWE"
replace eu = 1 if importer == "USA" & exporter == "SWE"
replace eu = 1 if exporter == "USA" & importer == "HRV"
replace eu = 1 if importer == "USA" & exporter == "HRV"
replace eu = 1 if exporter == "USA" & importer == "GRC"
replace eu = 1 if importer == "USA" & exporter == "GRC"
replace eu = 1 if exporter == "USA" & importer == "HUN"
replace eu = 1 if importer == "USA" & exporter == "HUN"
replace eu = 1 if exporter == "USA" & importer == "ROM"
replace eu = 1 if importer == "USA" & exporter == "ROM"
replace eu = 1 if exporter == "USA" & importer == "SVK"
replace eu = 1 if importer == "USA" & exporter == "SVK"
replace eu = 1 if exporter == "USA" & importer == "SVN"
replace eu = 1 if importer == "USA" & exporter == "SVN"

************************************************************************************************
* Step 3a: Estimate the Conditional Gravity Model and Construct the General Equilibrium Indexes*
************************************************************************************************
* Generate the counterfactual bilateral trade cost
g tij_cfl = exp(_b[ln_dist]*ln_dist+_b[contig]*contig+_b[comlang]*comlang+_b[colony]*colony+_b[brdr]*brdr+_b[eu]*eu+_b[nafta]*nafta)
g ln_tij_cfl = ln(tij_cfl)

* Drop fixed effects from the baseline scenario
drop exp_fe* imp_fe*
* Regenerate new fixed effects for the counterfactual model
qui tab exporter, g(exp_fe)
qui tab importer, g(imp_fe)

* Estimate the constrained counterfactual gravity model with the PPML estimator offsetting trade costs
ppml trade exp_fe* imp_fe1-imp_fe$noc_1, cluster(pairid) nocons offset(ln_tij_cfl) iter(30)
* Predict trade in the counterfactual scenario
predict tradehat_cfl_cdl, mu

* Based on the estimated exporter and importer fixed effects, create the actual set of fixed effects in the counterfactual scenario
forvalues j = 1 (1) $noc_1 {
qui replace exp_fe`j' = exp_fe`j' * exp(_b[exp_fe`j'])
qui replace imp_fe`j' = imp_fe`j' * exp(_b[imp_fe`j'])
}
* Create the exporter and importer fixed effects for the country of reference (South Africa) in the counterfactual scenario
qui replace exp_fe$noc = exp_fe$noc * exp(_b[exp_fe$noc])
qui replace imp_fe$noc = imp_fe$noc * exp(0)
* Creat two variables holding each of all estimated exporter and importer fixed effects in the counterfactual scenario
egen exp_pi_1 = rowtotal(exp_fe1-exp_fe$noc)
egen exp_chi_1 = rowtotal(imp_fe1-imp_fe$noc)

* Construct the indexes of interest from the conditional counterfactual scenario
g omr_cfl_cdl = y * e_deu / exp_pi_1
g imr_cfl_cdl = e / (exp_chi_1 * e_deu)
g rgdp_cfl_cdl_temp = y/((imr_cfl_cdl)^(1/(1-sigma))) if importer == exporter
bys exporter: egen rgdp_cfl_cdl = sum(rgdp_cfl_cdl_temp)
g exports_cfl_cdl = tradehat_cfl_cdl if importer != exporter
bys exporter: egen tot_exports_cfl_cdl = sum(exports_cfl_cdl)

****************************************************************************
*Step 4a: Construct the Conditional General Equilibrium Indexes of Interest*
****************************************************************************
preserve
* Obtain percentage changess in the OMR, real GDP and total exports indexes
collapse(mean) omr_bln omr_cfl_cdl rgdp_bln rgdp_cfl_cdl tot_exports_bln tot_exports_cfl_cdl, by(exporter)
g omr_index_cdl = (omr_cfl_cdl^(1/(1 - sigma)) - omr_bln^(1/(1 - sigma)))/omr_bln^(1/(1 - sigma))*100
g welfare_index_cdl = (rgdp_cfl_cdl - rgdp_bln)/rgdp_bln*100
g tot_exports_index_cdl = ((tot_exports_cfl_cdl - tot_exports_bln)/tot_exports_bln)*100
rename exporter country
replace country = "DEU" if country == "ZZZ"
keep country omr_index_cdl welfare_index_cdl tot_exports_index_cdl
save exp_indexes_cdl.dta, replace
restore
preserve
* Obtain percentage changes in the IMR index
collapse(mean) imr_bln imr_cfl_cdl, by(importer)
g imr_index_cdl = (imr_cfl_cdl^(1/(1 - sigma)) - imr_bln^(1/(1 - sigma)))/imr_bln^(1/(1 - sigma))*100
rename importer country
replace country = "DEU" if country == "ZZZ"
keep country imr_index_cdl
save imp_index_cdl.dta, replace

* Merge indexes of interest from the exporting side to the IMR index of interest
joinby country using "exp_indexes_cdl.dta"
save indexes_cdl.dta, replace
restore

***************************************************************************************************
* Step 3b: Estimate the Full Endowment Gravity Model and Construct the General Equilibrium Indexes*
***************************************************************************************************
* Define variables for the loop
local i = 3
local diff_exp_pi_sd = 1
local diff_exp_pi_max = 1
g double trade_1_pred = tradehat_cfl_cdl
g double y_bln = y
g double e_bln = e
g double phi = e/y if exporter == importer
g double temp = exp_pi_0 if exporter == importer
bys importer: egen double exp_pi_0_imp = mean(temp)
drop temp*
g double e_temp_1 = phi * y if exporter == importer
bys importer: egen double e_1 = mean(e_temp_1)
drop *temp*
g double e_deu01_1 = e_1 if importer == "ZZZ"
egen double e_deu1_1 = mean(e_deu01_1)
egen double e_deu_1 = mean(e_deu01_1)
g double temp = exp_pi_1 if exporter == importer
bys importer: egen double exp_pi_1_imp = mean(temp)
drop *temp*
g double p_full_exp_0 = 0
g double p_full_exp_1 = (exp_pi_1/exp_pi_0)^(1/(1-sigma))
g double p_full_imp_1 = (exp_pi_1_imp/exp_pi_0_imp)^(1/(1-sigma))
g double imr_full_1 = e_1/(exp_chi_1 * e_deu_1)
g double imr_full_ch_1 = 1
g double omr_full_1 = y * e_deu_1/(exp_pi_1)
g double omr_full_ch_1 = 1

* Start loop with convergence criteria that the mean and the standard devition of the difference in each of the factory-gate prices
* between two subsequent iterations is smaller than the pre-defined tolerance criterium of 0.05

while (`diff_exp_pi_sd' > 0.05) | (`diff_exp_pi_max') > 0.05 {

local i_1 = `i' - 1
local i_2 = `i' - 2
local i_3 = `i' - 3

* Allow for endogenous output/income, expenditures and trade
g double trade_`i_1' = trade_`i_2'_pred * p_full_exp_`i_2' * p_full_imp_`i_2'/(omr_full_ch_`i_2'*imr_full_ch_`i_2')

* Drop fixed effects and regenerate the new fixed effects for the full endowment scenario
drop exp_fe* imp_fe*
qui tab exporter, g(exp_fe)
qui tab importer, g(imp_fe)

* Estimate the structural gravity model and predict trade
capture glm trade_`i_1' exp_fe* imp_fe*, family(poisson) offset(ln_tij_cfl) nocons irls iter(30) 
predict trade_`i_1'_pred, mu

* Obtain estimates of the fixed effects from the full endowment scenario
forvalues m = 1 (1) $noc {
qui replace exp_fe`m' = exp_fe`m' * exp(_b[exp_fe`m'])
qui replace imp_fe`m' = imp_fe`m' * exp(_b[imp_fe`m'])
}
egen double exp_pi_`i_1' = rowtotal(exp_fe1-exp_fe$noc) 
egen double exp_chi_`i_1' = rowtotal(imp_fe1-imp_fe$noc)

* Update output
bys exporter: egen double y_`i_1' = total(trade_`i_1'_pred)

* Update expenditure of reference country
bys importer: egen double e_check_`i_1' = total(trade_`i_1'_pred)
g double e_deu0_`i_1' = e_check_`i_1' if importer == "ZZZ"
egen double e_deu_`i_1' = mean(e_deu0_`i_1')
g double temp = exp_pi_`i_1' if exporter == importer
bys importer: egen double exp_pi_`i_1'_imp = mean(temp)
drop temp*

* Update factory-gate prices and Multilateral Resistance terms
g double p_full_exp_`i_1' = ((exp_pi_`i_1'/exp_pi_`i_2')/(e_deu_`i_1'/e_deu_`i_2'))^(1/(1-sigma))
g double p_full_imp_`i_1' = ((exp_pi_`i_1'_imp/exp_pi_`i_2')/(e_deu_`i_1'/e_deu_`i_2'))^(1/(1-sigma))
g double omr_full_`i_1' = y_`i_1'/exp_pi_`i_1'
g double omr_full_ch_`i_1' = omr_full_`i_1'/omr_full_`i_2'

* Update expenditure
g double e_temp_`i_1' = phi * y_`i_1' if exporter == importer
bys importer: egen double e_`i_1' = mean(e_temp_`i_1')
g double imr_full_`i_1' = e_`i_1'/(exp_chi_`i_1' * e_deu_`i_1')
g double imr_full_ch_`i_1' = imr_full_`i_1'/imr_full_`i_2'

* Iterate until the change in factory-gate prices converges to zero or less than the tolerance criterium of 0.005
g double diff_p_full_exp_`i_1' = p_full_exp_`i_2' - p_full_exp_`i_3'
su diff_p_full_exp_`i_1'
local diff_exp_pi_sd = r(sd)
local diff_exp_pi_max = abs(r(max))
local i = `i' + 1
}
*********************************************End of Iterative Procedure*********************************************
drop y e e_deu_bln
rename e_deu e_deu_bln
bys exporter: egen double y = sum(trade_`i_2'_pred)
g double e_temp = phi * y if exporter == importer
bys importer: egen double e = mean(e_temp)
g double e_deu0 = e if importer == "ZZZ"
egen double e_deu = mean(e_deu0)

* Construct the indexes of interest from the full endowment counterfactual scenario
g double p_full = ((exp_pi_1/exp_pi_0)/(e_deu/e_deu_bln))^(1/(1-sigma))
g double y_full = p_full * y_bln
g double omr_full = y_full * e_deu/(exp_pi_1)
g double imr_full = e/(exp_chi_1 * e_deu)
g double rgdp_full_temp = y_full/imr_full^(1/(1-sigma)) if exporter == importer
bys exporter: egen double rgdp_full = sum(rgdp_full_temp)
g double e_full_temp = phi * y_full if exporter == importer
bys importer: egen double e_full = mean(e_full_temp)
g double trade_full = (y_full * e_full * tij_cfl)/(imr_full * omr_full)
g double exports_full_temp = trade_full if exporter != importer
bys exporter: egen double tot_exports_full = sum(exports_full_temp)

* Creat log of aggregate output for the graphs
g ln_y_bln = ln(y_bln)

*******************************************************************************
* Step 4b: Construct the Full Endowment General Equilibrium Indexes of Interst*
*******************************************************************************
preserve
* Obtain percentage changess in the OMR, real GDP total exports and factory-gate price indexes
collapse(mean) omr_bln omr_full rgdp_bln rgdp_full tot_exports_bln tot_exports_full p_full ln_y_bln, by(exporter)
g omr_index_full = (omr_full^(1/(1 - sigma)) - omr_bln^(1/(1 - sigma)))/omr_bln^(1/(1 - sigma))*100
g welfare_index_full = ((rgdp_full - rgdp_bln)/rgdp_bln)*100
g tot_exports_index_full = ((tot_exports_full - tot_exports_bln)/tot_exports_bln)*100
g fgp_index_full = (p_full - 1)*100
rename exporter country
replace country = "DEU" if country == "ZZZ"
keep country omr_index_full welfare_index_full tot_exports_index_full fgp_index_full ln_y_bln
save exp_indexes_full.dta, replace
restore
preserve
* Obtain percentage changes in the IMR index
collapse(mean) imr_bln imr_full, by(importer)
g imr_index_full = (imr_full^(1/(1 - sigma)) - imr_bln^(1/(1 - sigma)))/imr_bln^(1/(1 - sigma))*100
rename importer country
replace country = "DEU" if country == "ZZZ"
keep country imr_index_full
save imp_index_full.dta, replace

* Merge indexes from the exporting side to the IMR index of interest
joinby country using "exp_indexes_full.dta"
save indexes_full.dta, replace
* Merge full endowment gravity effects with conditional gravity effects in Step 4a
joinby country using "indexes_cdl.dta"
save ttip_geppml_indexes.dta, replace

* Export results to MS Excel
export excel using "ttip_geppml_indexes.xls ", firstrow(variables) replace

* Open geppml indexes and plot graphs
use "C:\Users\Awuse Nicholas\Desktop\Borger\ttip_geppml_indexes.dta"

* Keeping all member countries (EU and the US)
preserve
keep country ln_y_bln imr_index_full omr_index_full welfare_index_full tot_exports_index_full fgp_index_full
keep if fgp_index_full>0 | country == "USA"

* Scatter plot of output and the full endowment indexes for member countries
twoway scatter ln_y_bln imr_index_full, mlabel(country) mlabsize(small)
twoway scatter ln_y_bln fgp_index_full, mlabel(country) mlabsize(small)
twoway scatter ln_y_bln welfare_index_full , mlabel(country) mlabsize(small)
twoway scatter ln_y_bln tot_exports_index_full , mlabel(country) mlabsize(small)
restore

* Keeping all non-member countries in the data
preserve
keep country ln_y_bln imr_index_full omr_index_full welfare_index_full tot_exports_index_full fgp_index_full
keep if fgp_index_full<0
drop if country == "USA"

* Scatter plot of output and the full endowment indexes for non-member countries
twoway scatter ln_y_bln imr_index_full, mlabel(country) mlabsize(small)
twoway scatter ln_y_bln fgp_index_full, mlabel(country) mlabsize(small)
twoway scatter ln_y_bln welfare_index_full , mlabel(country) mlabsize(small)
twoway scatter ln_y_bln tot_exports_index_full , mlabel(country) mlabsize(small)
restore

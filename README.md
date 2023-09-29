## Investigating inequalities in HIV testing in sub-Saharan Africa: insights from a spatial analysis of 25 countries

##### Pearl Anne Ante-Testard, Gabriel Carrasco-Escobar, Tarik Benmarhnia, Laura Temime, Kévin Jean

##### Corresponding author: Pearl Anne Ante-Testard (pearlannemante@gmail.com/ pearl.ante@ucsf.edu)

## Systems Requirement
All analyses were running using R software version 4.2.1 on Mac OSX Ventura using the RStudio IDE (https://www.rstudio.com).

> sessionInfo()

R version 4.2.1 (2022-06-23)

Platform: x86_64-apple-darwin17.0 (64-bit)

Running under: macOS Ventura 13.4.1

## Installation Guide
You can download and install R from CRAN: https://cran.r-project.org

You can download and install RStudio from their website: https://www.rstudio.com

The installation time should be < 10 minutes total on a typical desktop computer.

## Main R files

##### analysis.Rproj
The R project for this study.
##### analysis_wealth.R
This calculates and maps the relative index of inequality and slope index of inequality at the cluster level.
##### analysis_wealth_reg.R
This calculates and maps the relative index of inequality and slope index of inequality at the regional/province level.
##### analysis_wealth_national.R
This calculates and maps the relative index of inequality and slope index of inequality at the national level.
##### globalG.R
Assessed the global G* statistic before the hotspot/coldspot spatial analysis.
##### hivtest_ineq_fun_2.R
Function needed in analysis_wealth.R to calculate the RII and SII at the national and province level.
##### hivtest_ineq_linear_clust_fun.R
Function needed in analysis_wealth.R to calculate the RII and SII at the cluster level.
##### spatialanalysis2.R
This conducts spatial analysis of inequalities at the cluster level using the Getis-Ord G* statistics.
##### ggscatter_hivprev_test.R
This plots the scatter plot between HIV prevalence and HIV testing at the national, province and cluster level.
##### plot_reg_hivprev_hivtest.R 
Maps the HIV prevalence and HIV testing at the province level for all countries and sexes.
##### plots_reg_sii_rii.R 
Maps the SII and RII at the province level for all countries and sexes.
##### suppl_distribution.R
Plots the sample size and % of HIV testing distribution at the cluster level. 


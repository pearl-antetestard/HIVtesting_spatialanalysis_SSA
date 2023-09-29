

################
## This script calculates inequality indicators (RII and SII) and maps them at the national level
## Author: Pearl Ante-Testard
## email: pearlannemante@gmail.com/ pearl.ante@ucsf.edu
###############

# cleans environment
rm(list = ls())

# packages
library(tidyverse)
library(sf)
library(MASS)
library(dplyr)
library(weights)
library(scales)
library(here)
library(cowplot)
library(spData)
library(colorspace)
library(Rcpp)
library(Amelia)
library(MatchIt)
library(VGAM)
library(Zelig)
library(msm)
library(geepack)
#remotes::install_version("Zelig", "5.1.7")
#install.packages("hadley/scales")
#library(rineq)

# reading csv files
gender = "female"
gender = "male"
gender = "all"


dat <- readRDS(here::here("mergedata_allcountries.rds")) 

  
if(gender=="female"){
  dat <- dat %>%
    dplyr::filter(gender=="female") %>%
    dplyr::filter(country!="BF") %>%
    mutate(w = sampleweight,
           cluster_number = cluster) %>%
    mutate(hivstatpos = case_when(hivstatneg == 1 ~ 0,
                                  hivstatneg == 0 ~ 1))
} else if(gender=="male"){
  dat <- dat %>%
    dplyr::filter(gender=="male") %>%
    dplyr::filter(country!="BF") %>%
    mutate(w = sampleweight,
           cluster_number = cluster) %>%
    mutate(hivstatpos = case_when(hivstatneg == 1 ~ 0,
                                  hivstatneg == 0 ~ 1))
} else if(gender=="all"){
  dat <- dat %>%
    dplyr::filter(country!="BF") %>%
    mutate(w = sampleweight,
           cluster_number = cluster) %>%
    mutate(hivstatpos = case_when(hivstatneg == 1 ~ 0,
                                  hivstatneg == 0 ~ 1))
}

# calling function to calculate RII and SII
source(here::here("hivtest_ineq_fun_2.R"))

################

# recoding variables
dat$w2=dat$w/1000000
dat$wealthindex=fct_relevel(dat$wealthindex,c("poorest","poorer","middle","richer","richest"))

#prev_val=aggregate(dat$HIVtest_12num,by=list(dat$REGNAME),mean,na.rm=T)
#dat$prev_reg = prev_val$x[match(dat$REGNAME, prev_val[,1] )]
#dat$prev_reg=dat$prev_reg*100 

################

# preparing data
if(gender=="female"){
  
  dat.n.f <- dat %>%
    group_by(country) %>%
    nest() %>%
    #mutate(sample.size.N = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12,HIVtest_12num,wealth_rank,wealthindex,wealthscore,w2) %>%
    mutate(sample.size.N = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12,HIVtest_12num,hivstatpos,
                                                              wealth_rank,wealthindex,wealthscore,w2,region.x,
                                                              cluster_number
                                                              ) %>%
                                 dplyr::count(name = "N") %>% unlist()),
           sample.size.m = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12num) %>% 
                                 filter(!is.na(HIVtest_12num)) %>%
                                 dplyr::tally(HIVtest_12num,name = "N_m",sort = T) %>% unlist()),
           #sample.size.p = map(.x = data, .f = ~dplyr::select(.x, hivstatpos) %>% 
           # filter(!is.na(hivstatpos)) %>%
           #dplyr::tally(hivstatpos,name = "N_p",sort = T) %>% unlist()),
           #dat_s = map(.x = data, .f = ~dplyr::select(.x,HIVtest_12,HIVtest_12num,wealth_rank,wealthindex,wealthscore,w2) %>% 
           dat_s = map(.x = data, .f = ~dplyr::select(.x,HIVtest_12,HIVtest_12num,hivstatpos,
                                                      wealth_rank,wealthindex,wealthscore,w2,region.x,
                                                      cluster_number
                                                      ) %>% 
                         filter(!is.na(HIVtest_12)),
                         filter(!is.na(wealthindex))),
                         #filter(complete.cases(.))),
           sample.size = map(.x = dat_s, .f = ~dplyr::count(.x) %>% unlist()),
           prev = map(.x = dat_s, .f = ~wpct(.x$HIVtest_12, .x$w2)[2] %>% unlist()), #weighted frequency table divided by its sum
           #prev = map(.x = dat_s, .f = ~mean(.x$HIVtest_12, na.rm=T) %>% unlist()),
           HIVprev = map(.x = dat_s, .f = ~wpct(.x$hivstatpos, .x$w2)[2] %>% unlist()),
           wealth_cats = map(.x = dat_s, .f = ~distinct(.x, wealthindex) %>% nrow())) %>%
    # FILTERS =======================
  dplyr::filter(prev > 0 & prev < 1) %>%
    #dplyr::filter(HIVprev > 0 & HIVprev < 1) #%>%
    #dplyr::filter(sample.size.m >=15) %>%
    dplyr::filter(sample.size>=10) %>%
    dplyr::filter(wealth_cats>1)
  
  #sample.size.N = sample size per cluster
  #sample.size.m = number of people tested for HIV per cluster
  
} else if (gender=="male") {
  
  dat.n.m <- dat %>%
    group_by(country) %>%
    nest() %>%
    #mutate(sample.size.N = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12,HIVtest_12num,wealth_rank,wealthindex,wealthscore,w2) %>%
    mutate(sample.size.N = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12,HIVtest_12num,hivstatpos,wealth_rank,wealthindex,wealthscore,w2,region.x,cluster_number) %>%
                                 dplyr::count(name = "N") %>% unlist()),
           sample.size.m = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12num) %>% 
                                 filter(!is.na(HIVtest_12num)) %>%
                                 dplyr::tally(HIVtest_12num,name = "N_m",sort = T) %>% unlist()),
           #sample.size.p = map(.x = data, .f = ~dplyr::select(.x, hivstatpos) %>% 
           # filter(!is.na(hivstatpos)) %>%
           #dplyr::tally(hivstatpos,name = "N_p",sort = T) %>% unlist()),
           #dat_s = map(.x = data, .f = ~dplyr::select(.x,HIVtest_12,HIVtest_12num,wealth_rank,wealthindex,wealthscore,w2) %>% 
           dat_s = map(.x = data, .f = ~dplyr::select(.x,HIVtest_12,HIVtest_12num,hivstatpos,wealth_rank,wealthindex,wealthscore,w2,region.x,cluster_number) %>% 
                         filter(!is.na(HIVtest_12)),
                         filter(!is.na(wealthindex))),
                         #filter(complete.cases(.))),
           sample.size = map(.x = dat_s, .f = ~dplyr::count(.x) %>% unlist()),
           prev = map(.x = dat_s, .f = ~wpct(.x$HIVtest_12, .x$w2)[2] %>% unlist()), #weighted frequency table divided by its sum
           #prev = map(.x = dat_s, .f = ~mean(.x$HIVtest_12, na.rm=T) %>% unlist()),
           HIVprev = map(.x = dat_s, .f = ~wpct(.x$hivstatpos, .x$w2)[2] %>% unlist()),
           wealth_cats = map(.x = dat_s, .f = ~distinct(.x, wealthindex) %>% nrow())) %>%
    # FILTERS =======================
  dplyr::filter(prev > 0 & prev < 1) %>%
    #dplyr::filter(HIVprev > 0 & HIVprev < 1) #%>%
    #dplyr::filter(sample.size.m >=15) %>%
    dplyr::filter(sample.size>=10) %>%
    dplyr::filter(wealth_cats>1)
  
  #sample.size.N = sample size per cluster
  #sample.size.m = number of people tested for HIV per cluster
  
} else if (gender=="all") {
  
  dat.n.all <- dat %>%
    group_by(country) %>%
    nest() %>%
    #mutate(sample.size.N = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12,HIVtest_12num,wealth_rank,wealthindex,wealthscore,w2) %>%
    mutate(sample.size.N = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12,HIVtest_12num,hivstatpos,wealth_rank,wealthindex,wealthscore,w2,region.x,cluster_number) %>%
                                 dplyr::count(name = "N") %>% unlist()),
           sample.size.m = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12num) %>% 
                                 filter(!is.na(HIVtest_12num)) %>%
                                 dplyr::tally(HIVtest_12num,name = "N_m",sort = T) %>% unlist()),
           #sample.size.p = map(.x = data, .f = ~dplyr::select(.x, hivstatpos) %>% 
           # filter(!is.na(hivstatpos)) %>%
           #dplyr::tally(hivstatpos,name = "N_p",sort = T) %>% unlist()),
           #dat_s = map(.x = data, .f = ~dplyr::select(.x,HIVtest_12,HIVtest_12num,wealth_rank,wealthindex,wealthscore,w2) %>% 
           dat_s = map(.x = data, .f = ~dplyr::select(.x,HIVtest_12,HIVtest_12num,hivstatpos,wealth_rank,wealthindex,wealthscore,w2,region.x,cluster_number) %>% 
                         filter(!is.na(HIVtest_12)),
                         filter(!is.na(wealthindex))),
                         #filter(complete.cases(.))),
           sample.size = map(.x = dat_s, .f = ~dplyr::count(.x) %>% unlist()),
           prev = map(.x = dat_s, .f = ~wpct(.x$HIVtest_12, .x$w2)[2] %>% unlist()), #weighted frequency table divided by its sum
           #prev = map(.x = dat_s, .f = ~mean(.x$HIVtest_12, na.rm=T) %>% unlist()),
           HIVprev = map(.x = dat_s, .f = ~wpct(.x$hivstatpos, .x$w2)[2] %>% unlist()),
           wealth_cats = map(.x = dat_s, .f = ~distinct(.x, wealthindex) %>% nrow())) %>%
    # FILTERS =======================
  dplyr::filter(prev > 0 & prev < 1) %>%
    #dplyr::filter(HIVprev > 0 & HIVprev < 1) #%>%
    #dplyr::filter(sample.size.m >=15) %>%
    dplyr::filter(sample.size>=10) %>%
    dplyr::filter(wealth_cats>1)
  
  #sample.size.N = sample size per cluster
  #sample.size.m = number of people tested for HIV per cluster
  
  
}


################

# calculating RII and SII
if(gender=="female"){
  
  dat_ci.f <- dat.n.f  %>%
    mutate(
      h_calc = map(.x = dat_s, .f = ~h_ineq(dat = .x, var_soc = wealthindex, var_outcome = HIVtest_12num
                                            ,clust = cluster_number
      )) 
    )
  
  dat_ci_f <- dat_ci.f %>%
    dplyr::select(country,dat_s,sample.size.N, sample.size.m, sample.size, prev, HIVprev, 
                  wealth_cats, h_calc) %>%
    unnest(cols = c(sample.size.N, sample.size.m,sample.size, prev, HIVprev, 
                    wealth_cats,h_calc)) %>%
    ungroup() %>%
    dplyr::filter(rii > 0.1 & rii < 300)
  
  
} else if(gender=="male"){
  
  dat_ci.m <- dat.n.m %>%
    mutate(
      h_calc = map(.x = dat_s, .f = ~h_ineq(dat = .x, var_soc = wealthindex, var_outcome = HIVtest_12num
                                            ,clust = cluster_number
      )) 
    )
  
  dat_ci_m <- dat_ci.m %>%
    dplyr::select(country,dat_s,sample.size.N, sample.size.m, sample.size, prev, HIVprev, wealth_cats, h_calc) %>%
    unnest(cols = c(sample.size.N, sample.size.m,sample.size, prev, HIVprev, wealth_cats,h_calc)) %>%
    ungroup() %>%
    dplyr::filter(rii > 0.1 & rii < 300)
  
} else if(gender=="all"){
  
  
  dat_ci.all <- dat.n.all %>%
    mutate(
      h_calc = map(.x = dat_s, .f = ~h_ineq(dat = .x, var_soc = wealthindex, var_outcome = HIVtest_12num
                                            ,clust = cluster_number
      )) 
    )
  
  dat_ci_all <- dat.n.all %>%
    dplyr::select(country,dat_s,sample.size.N, sample.size.m, sample.size, prev, HIVprev, wealth_cats, h_calc) %>%
    unnest(cols = c(sample.size.N, sample.size.m,sample.size, prev, HIVprev, wealth_cats,h_calc)) %>%
    ungroup() %>%
    dplyr::filter(rii > 0.1 & rii < 300)
} 


################

saveRDS(dat_ci_f, paste0(here::here("dat_sii-rii_national_map_latest_v3"),gender,".rds"))
write_csv(dat_ci_f, paste0(here::here("dat_sii-rii_national_map_latest_v3"),gender,".csv"))

saveRDS(dat_ci_m, paste0(here::here("dat_sii-rii_national_map_latest_v3"),gender,".rds"))
write_csv(dat_ci_m, paste0(here::here("dat_sii-rii_national_map_latest_v3"),gender,".csv"))
#saveRDS(dat_ci_map, paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/dat_sii-rii_linear_map2_HIVprev",gender,".rds"))

################

saveRDS(dat_ci_f, paste0(here::here("dat_sii-rii_national_map_latest_v4"),gender,".rds"))
write_csv(dat_ci_f, paste0(here::here("dat_sii-rii_national_map_latest_v4"),gender,".csv"))

saveRDS(dat_ci_m, paste0(here::here("dat_sii-rii_national_map_latest_v4"),gender,".rds"))
write_csv(dat_ci_m, paste0(here::here("dat_sii-rii_national_map_latest_v4"),gender,".csv"))


################
## This script calculates inequality indicators (RII and SII) and maps them
## Author: Pearl Ante-Testard
## email: pearlannemante@gmail.com/ pearl.ante@ucsf.edu
###############

# cleans environment
rm(list = ls())

# packages
library(codebook)
library(tidyverse)
library(sf)
library(MASS)
library(dplyr)
library(weights)
library(scales)
library(here)
library(geepack)
library(cowplot)
library(sf)
library(mapview)
library(viridis)
library(spData)
library(colorspace)
#install.packages("hadley/scales")
library(scales)
#library(rineq)

# reading csv files
gender = "female"
gender = "male"
gender = "all"

if(gender=="female"){
  dat <- read.csv(here("dat_mal_ineq_v3.csv")) %>%
    dplyr::filter(gender=="female") %>%
    dplyr::filter(country!="BF")
} else if(gender=="male"){
  dat <- read.csv(here("dat_mal_ineq_v3.csv")) %>%
    dplyr::filter(gender=="male")%>%
    dplyr::filter(country!="BF")
} else if(gender=="all"){
  dat <- read.csv(here("dat_mal_ineq_v3.csv")) #%>%
    #dplyr::filter(gender=="male")
}

# calling function to calculate RII and SII
source(here("hivtest_ineq_fun_2.R"))
  
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
    group_by(country,REGNAME) %>%
    nest() %>%
    mutate(sample.size.N = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12,HIVtest_12num,wealth_rank,wealthindex,wealthscore,w2,cluster_number) %>%
    #mutate(sample.size.N = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12,HIVtest_12num,hivstatpos,wealth_rank,wealthindex,wealthscore,w2,cluster_number) %>%
                                        dplyr::count(name = "N") %>% unlist()),
           sample.size.m = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12) %>% 
                                 filter(!is.na(HIVtest_12)) %>%
                                 dplyr::tally(HIVtest_12,name = "N_m",sort = T) %>% unlist()),
           #sample.size.p = map(.x = data, .f = ~dplyr::select(.x, hivstatpos) %>% 
                               # filter(!is.na(hivstatpos)) %>%
                                #dplyr::tally(hivstatpos,name = "N_p",sort = T) %>% unlist()),
           dat_s = map(.x = data, .f = ~dplyr::select(.x,HIVtest_12,HIVtest_12num,wealth_rank,wealthindex,wealthscore,w2,cluster_number) %>% 
           #dat_s = map(.x = data, .f = ~dplyr::select(.x,HIVtest_12,HIVtest_12num,hivstatpos,wealth_rank,wealthindex,wealthscore,w2,cluster_number) %>% 
                         filter(complete.cases(.))),
           sample.size = map(.x = dat_s, .f = ~dplyr::count(.x) %>% unlist()),
           prev = map(.x = dat_s, .f = ~wpct(.x$HIVtest_12, .x$w2)[2] %>% unlist()), #weighted frequency table divided by its sum
           #prev = map(.x = dat_s, .f = ~mean(.x$HIVtest_12, na.rm=T) %>% unlist()),
           #HIVprev = map(.x = dat_s, .f = ~wpct(.x$hivstatpos, .x$w2)[2] %>% unlist()),
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
    group_by(country,REGNAME) %>%
    nest() %>%
    mutate(sample.size.N = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12,HIVtest_12num,wealth_rank,wealthindex,wealthscore,w2,cluster_number) %>%
                                 #mutate(sample.size.N = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12,HIVtest_12num,hivstatpos,wealth_rank,wealthindex,wealthscore,w2,cluster_number) %>%
                                 dplyr::count(name = "N") %>% unlist()),
           sample.size.m = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12) %>% 
                                 filter(!is.na(HIVtest_12)) %>%
                                 dplyr::tally(HIVtest_12,name = "N_m",sort = T) %>% unlist()),
           #sample.size.p = map(.x = data, .f = ~dplyr::select(.x, hivstatpos) %>% 
                                 #filter(!is.na(hivstatpos)) %>%
                                 #dplyr::tally(hivstatpos,name = "N_p",sort = T) %>% unlist()),
           dat_s = map(.x = data, .f = ~dplyr::select(.x,HIVtest_12,HIVtest_12num,wealth_rank,wealthindex,wealthscore,w2,cluster_number) %>% 
                         #dat_s = map(.x = data, .f = ~dplyr::select(.x,HIVtest_12,HIVtest_12num,hivstatpos,wealth_rank,wealthindex,wealthscore,w2,cluster_number) %>% 
                         filter(complete.cases(.))),
           sample.size = map(.x = dat_s, .f = ~dplyr::count(.x) %>% unlist()),
           prev = map(.x = dat_s, .f = ~wpct(.x$HIVtest_12, .x$w2)[2] %>% unlist()), #weighted frequency table divided by its sum
           #prev = map(.x = dat_s, .f = ~mean(.x$HIVtest_12, na.rm=T) %>% unlist()),
           #HIVprev = map(.x = dat_s, .f = ~wpct(.x$hivstatpos, .x$w2)[2] %>% unlist()),
           wealth_cats = map(.x = dat_s, .f = ~distinct(.x, wealthindex) %>% nrow())) %>%
    # FILTERS =======================
    dplyr::filter(prev > 0 & prev < 1) %>%
    #dplyr::filter(HIVprev > 0 & HIVprev < 1) #%>%
    #dplyr::filter(sample.size.m >=15) %>%
    dplyr::filter(sample.size>=10) %>%
    dplyr::filter(wealth_cats>1)
  
  #sample.size.N = sample size per cluster
  #sample.size.m = number of people tested for HIV per cluster
  #sample.size.N = sample size per cluster
  #sample.size.m = number of people tested for HIV per cluster
} else if (gender=="all") {
  
  dat.n.all <- dat %>%
    group_by(country, REGNAME) %>%
    nest() %>%
    mutate(sample.size.N = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12,HIVtest_12num,wealth_rank,wealthindex,wealthscore,w2,cluster_number) %>%
                                 dplyr::count(name = "N") %>% unlist()),
           sample.size.m = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12) %>% 
                                 filter(!is.na(HIVtest_12)) %>%
                                 dplyr::tally(HIVtest_12,name = "N_m",sort = T) %>% unlist()),
           #sample.size.p = map(.x = data, .f = ~dplyr::select(.x, hivstatpos) %>% 
                                 #filter(!is.na(hivstatpos)) %>%
                                 #dplyr::tally(hivstatpos,name = "N_p",sort = T) %>% unlist()),
           dat_s = map(.x = data, .f = ~dplyr::select(.x, HIVtest_12,HIVtest_12num,wealth_rank,wealthindex,wealthscore,w2,cluster_number2) %>% 
                         filter(complete.cases(.))),
           sample.size = map(.x = dat_s, .f = ~dplyr::count(.x) %>% unlist()),
           prev = map(.x = dat_s, .f = ~wpct(.x$HIVtest_12, .x$w2)[2] %>% unlist()),
           #HIVprev = map(.x = dat_s, .f = ~wpct(.x$hivstatpos, .x$w2)[2] %>% unlist()),
           wealth_cats = map(.x = dat_s, .f = ~distinct(.x, wealthindex) %>% nrow())) %>%
    # FILTERS =======================
    dplyr::filter(prev > 0 & prev < 1) %>%
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
      h_calc = map(.x = dat_s, .f = ~h_ineq(dat = .x, var_soc = wealthindex, var_outcome = HIVtest_12num,
                                            clust = cluster_number))
    )
  
    dat_ci_f <- dat_ci.f %>%
    dplyr::select(country,REGNAME,dat_s,sample.size.N, sample.size.m, sample.size, prev,wealth_cats,h_calc) %>%
    unnest(cols = c(sample.size.N, sample.size.m,sample.size, prev,wealth_cats,h_calc)) %>%
    ungroup() %>%
    dplyr::filter(rii > 0.1 & rii < 300)

  
} else if(gender=="male"){

  dat_ci.m <- dat.n.m %>%
    mutate(
      h_calc = map(.x = dat_s, .f = ~h_ineq(dat = .x, var_soc = wealthindex, var_outcome = HIVtest_12,
                                            clust = cluster_number))
    )
  
  dat_ci_m <- dat_ci.m %>%
    dplyr::select(country,REGNAME,dat_s,sample.size.N, sample.size.m, sample.size, prev,wealth_cats,h_calc) %>%
    unnest(cols = c(sample.size.N, sample.size.m,sample.size, prev,wealth_cats,h_calc)) %>%
    ungroup() %>%
    dplyr::filter(rii > 0.1 & rii < 300)
  
} else if(gender=="all"){

  
  dat_ci.all <- dat.n.all  %>%
    #ungroup()
    mutate(
      h_calc = map(.x = dat_s, .f = ~h_ineq(dat = .x, var_soc = wealth_rank, var_outcome = HIVtest_12)))
  
  
  
  dat_ci_all <- dat_ci.all %>%
    dplyr::select(country, REGNAME, cluster_number, dat_s,sample.size.N, sample.size.m, sample.size, prev,wealth_cats,h_calc) %>%
    unnest(cols = c(sample.size.N, sample.size.m,sample.size, prev,HIVprev,wealth_cats,h_calc)) %>%
    ungroup() 
} 

################
  
# distributions of wealth 

dat_ci = dat_ci_f
dat_ci = dat_ci_m
#dat_ci = dat_ci_all

a <- hist_ineq(dat_ci, sample.size, "Sample size per cluster")
b <- hist_ineq(dat_ci, prev, "Proportion of HIV Testing")
c <- bi_hist_ineq(dat_ci, var_x = prev, var_y = sample.size, lab_x = "Proportion of HIV Testing", lab_y = "Sample size per cluster")

(dist <- plot_grid(a,b,c, labels = c("A)", "B)", "C)"), nrow = 1))

ggsave(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/res/",gender,"_psu2",".jpg"),
       units = "cm", width = 30, height = 17, dpi = 600)

################

# joining DHS dataset and spatial data

if(gender=="female"){
  
  dat_ci_map <- dat_ci_f %>%
    inner_join(read.csv("/Users/pearlanneante-testard/Dropbox/PhD/mapping/dat_gps_flat.csv"), by = c("country","cluster_number")) %>%
    st_as_sf(coords = c("LONGNUM", "LATNUM"), crs = 4326) %>%
    dplyr::mutate(lat = sf::st_coordinates(.)[,2],
                  long = sf::st_coordinates(.)[,1])
} else if(gender=="male"){
  
  dat_ci_map <- dat_ci_m %>%
    inner_join(read.csv("/Users/pearlanneante-testard/Dropbox/PhD/mapping/dat_gps_flat.csv"), by = c("country","cluster_number")) %>%
    st_as_sf(coords = c("LONGNUM", "LATNUM"), crs = 4326) %>%
    dplyr::mutate(lat = sf::st_coordinates(.)[,2],
                  long = sf::st_coordinates(.)[,1])
} else if(gender=="all"){
  
  dat_ci_map <- dat_ci_all %>%
    inner_join(read.csv("/Users/pearlanneante-testard/Dropbox/PhD/mapping/dat_gps_flat.csv"), by = c("country","cluster_number")) %>%
    st_as_sf(coords = c("LONGNUM", "LATNUM"), crs = 4326) %>%
    dplyr::mutate(lat = sf::st_coordinates(.)[,2],
                  long = sf::st_coordinates(.)[,1])
}


saveRDS(dat_ci_map, paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/dat_sii-rii_linear_map_latest2",gender,".rds"))
#saveRDS(dat_ci_map, paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/dat_sii-rii_linear_map2_HIVprev",gender,".rds"))
 
  
if(gender=="female"){
  
  #dat_adm <- st_read("/Users/pearlanneante-testard/Dropbox/PhD/mapping/SHP/union/DHS_adm.shp") %>%
  dat_adm <- st_read(here("union", "DHS_adm.shp")) %>%
    inner_join(dat_ci_f, by = c("country","REGNAME")) 
} else if(gender=="male"){
  
  dat_adm <- st_read("/Users/pearlanneante-testard/Dropbox/PhD/mapping/SHP/union/DHS_adm.shp") %>%
    inner_join(dat_ci_m, by = c("country", "REGNAME")) 
  
} else if(gender=="all"){
  
  dat_adm <- st_read("/Users/pearlanneante-testard/Dropbox/PhD/mapping/SHP/union/DHS_adm.shp") %>%
    inner_join(dat_ci_all, by = c("country", "REGNAME")) 
  
}

saveRDS(dat_adm, paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/dat_sii-rii_linear_adm_ggscatter",gender,".rds"))


### SII and RII by region
#dat_ci_adm <- dat_adm %>%
  #group_by(country, REGNAME) %>%
  #mutate(N = sample.size.N,
            #N_m = sample.size.m,
            #n = sample.size)
            #HIVprev_val = median(HIVprev, na.rm = T),
            #ci_val = median(c_index, na.rm = T),
            #sd_ci = sd(c_index, na.rm = T),
            #ci_val_er = median(ci_val_er, na.rm = T),
            #ci_val_rel = median(ci_val_rel, na.rm = T))
            #sii_val = median(sii, na.rm = T),
            #sd_sii = sd(sii, na.rm = T),
            #sii_low = median(sii_low, na.rm = T),
            #sii_up = median(sii_up, na.rm = T),
            #rii_val = median(rii, na.rm = T))
            #sd_rii = sd(rii, na.rm = T)) 
            #rii_low = median(rii_low, na.rm = T),
            #rii_up = median(rii_up, na.rm = T)) 

saveRDS(dat_ci_adm, paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/",gender,"_dat_ci_adm_v3.rds"))
#write.csv(dat_ci_adm, paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/",gender,"_countrylevel.csv"))
#d=read.csv("/Users/pearlanneante-testard/Dropbox/PhD/mapping/dat_mal_ineq_v3.csv")

################
  
# mapping SII

sPDF <- world %>%
  filter(continent == "Africa")

#sPDF$country <- car::recode(sPDF$iso_a2, "'AO'='AO';'BI'='BU';'CD'='CD';'CI'='CI';'CM'='CM';
               #            'ET'='ET';'GA'='GA';'GH'='GH';'GN'='GN';'LR'='LB';'LS'='LS';
                #           'ML'='ML';'MW'='MW';'MZ'='MZ';'NA'='NM';'RW'='RW';'SL'='SL';
                 #          'SN'='SN';'TD'='TD';'TG'='TG';'TZ'='TZ';'UG'='UG';'ZA'='ZA';
                   #        'ZM'='ZM';'ZW'='ZW'")

# female

# country level
fs_val=ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_adm %>%
          filter(!is.na(sii)), #%>%
          #mutate(log_sii = log(sii)) %>%
          #filter(log_sii != Inf & log_sii != -Inf), 
          aes(fill = sii),
          size = 0) +
  geom_sf(data = sPDF, fill = NA) +
  scale_fill_continuous_diverging(palette = "Blue-Red 3") +
  geom_sf(data = sPDF, fill = NA) +
  labs(fill = paste0("SII, ",gender)) +
  theme_bw()

# region level
fsii=ggplot() +
  geom_sf(data = map, fill = "grey") + 
  geom_sf(data = dat_ci_map, #%>%
          #filter(!is.na(sii)) %>%
          #mutate(log_sii = log(sii)), #%>%
          #filter(log_sii != Inf & log_sii != -Inf), 
          aes(col = sii),
          size = 0) +
  geom_sf(data = map, fill = NA) +
  scale_colour_continuous_diverging(palette = "Blue-Red 3") +
  geom_sf(data = map, fill = NA) +
  labs(col = paste0("PSU SII, ",gender)) +
  theme_bw()

# male

# country level
ms_val=ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_adm, # %>%
          #filter(!is.na(sii)), #%>%
          #mutate(log_sii = log(sii)) %>%
          #filter(log_sii != Inf & log_sii != -Inf), 
          aes(fill = sii),
          size = 0) +
  geom_sf(data = sPDF, fill = NA) +
  scale_fill_continuous_diverging(palette = "Blue-Red 3") +
  geom_sf(data = sPDF, fill = NA) +
  labs(fill = paste0("SII, ",gender)) +
  theme_bw()

# region level
msii=ggplot() +
  geom_sf(data = map, fill = "grey") + 
  geom_sf(data = dat_ci_map, #%>%
          #filter(!is.na(sii)) %>%
          #mutate(log_sii = log(sii)), #%>%
          #filter(log_sii != Inf & log_sii != -Inf), 
          aes(col = sii),
          size = 0) +
  geom_sf(data = map, fill = NA) +
  scale_colour_continuous_diverging(palette = "Blue-Red 3") +
  geom_sf(data = map, fill = NA) +
  labs(col = paste0("PSU SII, ",gender)) 

# mapping RII

source("/Users/pearlanneante-testard/Dropbox/Pearl - PhD project/mapping/log_noNaNs.R")

# female

# country level
frii_val = ggplot() +
  geom_sf(data = sPDF, fill = "grey") +
  geom_sf(data = dat_adm, #%>%
            #filter(!is.na(rii)) %>%
            #mutate(log_rii = transform_to_log_scale(rii)) %>%
            #filter(log_rii != Inf & log_rii != -Inf), 
          aes(fill = rii),
          size = 0) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_fill_continuous_diverging(palette = "Blue-Red 3") +
  geom_sf(data = sPDF, fill = NA) + 
  labs(fill = paste0("RII, ",gender)) +
  theme_bw()

#breaks = 10**(1:10)
#scale_y_log10(breaks = breaks, labels = comma(breaks))

# region level
frii=ggplot() +
  geom_sf(data = map, fill = "grey") + 
  geom_sf(data = dat_ci_map %>%
            filter(!is.na(rii)) %>%
            mutate(log_rii = log(rii)) %>%
            filter(log_rii != Inf & log_rii != -Inf), 
          aes(col = log_rii),
          size = 0) +
  geom_sf(data = map, fill = NA) +
  scale_colour_continuous_diverging(palette = "Blue-Red 3") +
  geom_sf(data = map, fill = NA) +
  labs(col = paste0("PSU log RII,\n",gender)) +
  theme_bw()

# male

# country level
mrii_val = ggplot() +
  geom_sf(data = sPDF, fill = "grey") +
  #geom_sf(data = dat_ci_adm, aes(fill = rii_val),
  #   size = 0) +
  geom_sf(data = dat_adm %>%
            filter(!is.na(rii)) %>%
            mutate(log_rii = transform_to_log_scale(rii)) %>%
            filter(log_rii != Inf & log_rii != -Inf), 
          aes(fill = log_rii),
          size = 0) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_fill_continuous_diverging(palette = "Blue-Red 3") +
  geom_sf(data = sPDF, fill = NA) + 
  labs(fill = paste0("log RII,",gender)) +
  theme_bw()

# region level
mrii=ggplot() +
  geom_sf(data = sPDF, fill = "grey") +
  #geom_sf(data = dat_ci_adm, aes(fill = rii_val),
  #   size = 0) +
  geom_sf(data = dat_adm, #%>%
          #filter(!is.na(rii)) %>%
          #mutate(log_rii = transform_to_log_scale(rii)) %>%
          #filter(log_rii != Inf & log_rii != -Inf), 
          aes(fill = rii),
          size = 0) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_fill_continuous_diverging(palette = "Blue-Red 3",trans = log10_trans()
                                  ) +
  geom_sf(data = sPDF, fill = NA) + 
  labs(fill = paste0("RII, ",gender)) +
  theme_bw()


###########
# combining by gender, region


ind="rii"
lev="reg"

fem_reg_sii_rii <- plot_grid(fs_val, frii_val,ncol = 2,
                             labels = c("A", "B"))

male_reg_sii_rii <- plot_grid(ms_val, mrii,ncol = 2,
                              labels = c("C", "D"))

fm_reg_sii_rii <- plot_grid(fem_reg_sii_rii, male_reg_sii_rii, nrow=2)

ggsave(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/res/fm_reg_sii_rii5_geeglm.jpg"),
       units = "cm", width = 40, height = 20, dpi = 650, fm_reg_sii_rii)


# combining by gender, PSU

fem_psu_sii_rii <- plot_grid(fsii, frii,ncol = 2,
                             labels = c("A", "B"))

male_psu_sii_rii <- plot_grid(msii, mrii,ncol = 2,
                              labels = c("C", "D"))

fm_psu_sii_rii <- plot_grid(fem_psu_sii_rii, male_psu_sii_rii, nrow=2)

ggsave(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/res/fm_psu_sii_rii.jpg"),
       units = "cm", width = 40, height = 20, dpi = 650, fm_psu_sii_rii)


###########

# mapping HIV testing uptake

fprev_val=ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_adm %>%
            mutate(prev_perc = prev*100),
          aes(fill = prev_perc),
          size = 0) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_fill_continuous_diverging(palette = "Vik", rev = T,
                                  breaks=c(0,20,40,60,80)) +
  #scale_fill_binned(type = "viridis") +
  labs(fill = paste0("Recent HIV testing\nuptake (%), ",gender)) +
  theme_bw()
fprev_val

mprev_val=ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_adm %>%
            mutate(prev_perc = prev*100),
          aes(fill = prev_perc),
          size = 0) +
  #scale_fill_gradientn(
  #colors=c("violet","blue"))+
  #rescale(mini_map:maxi_map))+
  #values=rescale(c(mini,med,maxi)))+ 
  geom_sf(data = sPDF, fill = NA) + 
  scale_fill_continuous_diverging(palette = "Vik", rev = T)+
                                  #breaks=c(0,20,40,60,80))+
  labs(fill = paste0("Recent HIV testing\nuptake (%), ",gender)) +
  theme_bw()
mprev_val

pr10 <- plot_grid(fprev_val, mprev_val, ncol = 2,
                  labels = c("C", "D"))

ind="percHIVtest"
ggsave(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/res/",ind,"3.jpg"),
       units = "cm", width = 35, height = 17, dpi = 650, pr10)



# mapping HIV prevalence

fhivprev_val=ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_adm %>%
            #filter(!is.na(HIVprev)) %>%
            mutate(HIVprev_perc = HIVprev*100),
            aes(fill = HIVprev_perc),
            size = 0) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_fill_continuous_diverging(palette = "Green-Orange") +
  labs(fill = paste0("Weighted HIV \nprevalence (%), ",gender)) +
  theme_bw()
fhivprev_val

mhivprev_val=ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_adm %>%
            #filter(!is.na(HIVprev)) %>%
            mutate(HIVprev_perc = HIVprev*100),
          aes(fill = HIVprev_perc),
          size = 0) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_fill_continuous_diverging(palette = "Green-Orange") +
  labs(fill = paste0("Weighted HIV \nprevalence (%), ",gender)) +
  theme_bw()
mhivprev_val

###########

# joining HIV prevalence for female and male
pr11 <- plot_grid(fhivprev_val, mhivprev_val, ncol = 2,
                  labels = c("A", "B"))

# saving 
ind="hivprev"
ggsave(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/res/",ind,"3.jpg"),
       units = "cm", width = 35, height = 17, dpi = 650, pr11)


# joining HIV testing uptake for female and male
pr12 <- plot_grid(pr11, pr10, nrow = 2)

ind="hivprev-hivtest"
ggsave(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/res/",ind,"map.jpg"),
       units = "cm", width = 35, height = 17, dpi = 650, pr12)

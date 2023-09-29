
# mapping HIV prevalence

rm(list = ls())

library(sp)
library(spdep)
library(dplyr)
library(ggplot2)
library(colorspace)
library(scales)
library(cowplot)

gender = "female"
gender = "male"

if(gender=="female"){
  dat_ci_f=readRDS(paste0(here::here("dat_sii-rii_region_map_latest_v4"),gender,".rds")) 
  dat_ci_f <- dat_ci_f %>%
    st_as_sf(coords = c("long", "lat"),crs = 4326) %>%
    dplyr::mutate(lat = sf::st_coordinates(.)[,2],
                  long = sf::st_coordinates(.)[,1]) %>%
    mutate(region = toupper(region))
  
  dat_adm <- st_read(here::here("union", "DHS_adm.shp"))
  dat_adm <- st_make_valid(dat_adm) %>%
    mutate(region = REGNAME) 
  dat_ci_adm <- st_join(dat_adm, dat_ci_f) 
} else if(gender=="male"){
  dat_ci_m=readRDS(paste0(here::here("dat_sii-rii_region_map_latest_v4"),gender,".rds"))
  dat_ci_m <- dat_ci_m %>%
    st_as_sf(coords = c("long", "lat"),crs = 4326) %>%
    dplyr::mutate(lat = sf::st_coordinates(.)[,2],
                  long = sf::st_coordinates(.)[,1]) %>%
    mutate(region = toupper(region))
  
  dat_adm <- st_read(here::here("union", "DHS_adm.shp"))
  dat_adm <- st_make_valid(dat_adm) %>%
    mutate(region = REGNAME) 
  dat_ci_adm <- st_join(dat_adm, dat_ci_m) 
} 

dat_ci_adm <-dat_ci_adm[!duplicated(dat_ci_adm$REG_ID),]

################ female

sPDF <- world %>%
  filter(continent == "Africa")


fhivprev_val=ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_ci_adm %>%
            #filter(!is.na(HIVprev)) %>%
            mutate(HIVprev_perc = HIVprev*100),
          aes(fill = HIVprev_perc),
          size = 0) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_fill_continuous_diverging(palette = "Green-Orange") +
  labs(fill = paste0("Weighted HIV \nprevalence (%), ",gender), tag = "A") +
  theme_bw()
fhivprev_val <- fhivprev_val + theme(legend.text = element_text(size=8), 
                            plot.tag = element_text(face="bold")) 


fhivtest_val=ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_ci_adm %>%
            #filter(!is.na(HIVprev)) %>%
            mutate(HIVtest_perc = prev*100),
          aes(fill = HIVtest_perc),
          size = 0) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_fill_continuous_diverging(palette = "Vik") +
  labs(fill = paste0("Recent HIV testing\nuptake (%), ",gender), tag = "B") +
  theme_bw()
fhivtest_val <- fhivtest_val + theme(legend.text = element_text(size=8), 
                                     plot.tag = element_text(face="bold")) 


########### male

mhivprev_val=ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_ci_adm %>%
            #filter(!is.na(HIVprev)) %>%
            mutate(HIVprev_perc = HIVprev*100),
          aes(fill = HIVprev_perc),
          size = 0) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_fill_continuous_diverging(palette = "Green-Orange") +
  labs(fill = paste0("Weighted HIV \nprevalence (%), ",gender), tag = "C") +
  theme_bw()
mhivprev_val <- mhivprev_val + theme(legend.text = element_text(size=8), 
                                     plot.tag = element_text(face="bold"))

mhivtest_val=ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_ci_adm %>%
            #filter(!is.na(HIVprev)) %>%
            mutate(HIVtest_perc = prev*100),
          aes(fill = HIVtest_perc),
          size = 0) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_fill_continuous_diverging(palette = "Vik") +
  labs(fill = paste0("Recent HIV testing\nuptake (%), ",gender), tag = "D") +
  theme_bw()
mhivtest_val <- mhivtest_val + theme(legend.text = element_text(size=8), 
                                     plot.tag = element_text(face="bold")) 

###########
# combining by gender, region

freg_hivprev_hivtest <- plot_grid(fhivprev_val, fhivtest_val,ncol = 2)

mreg_hivprev_hivtest <- plot_grid(mhivprev_val, mhivtest_val,ncol = 2)

fm_reg_hivprev_hivtest <- plot_grid(freg_hivprev_hivtest, mreg_hivprev_hivtest, nrow=2)

ggsave(here::here(paste0("fm_reg_hivprev_hivtest.pdf")),
       units = "cm", width = 30, height = 20, dpi = 1200, fm_reg_hivprev_hivtest)

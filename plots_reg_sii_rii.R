
rm(list = ls())

library(sp)
library(spdep)
library(dplyr)
library(ggplot2)
library(colorspace)
library(scales)
library(cowplot)

gender = "female"
dat_ci_f=readRDS(paste0(here::here("dat_sii-rii_region_map_latest_v4"),gender,".rds")) 
dat_ci_f <- dat_ci_f %>%
  st_as_sf(coords = c("long", "lat"),crs = 4326) %>%
  dplyr::mutate(lat = sf::st_coordinates(.)[,2],
                long = sf::st_coordinates(.)[,1]) %>%
  mutate(region = toupper(region))

gender = "male"
dat_ci_m=readRDS(paste0(here::here("dat_sii-rii_region_map_latest_v4"),gender,".rds"))
dat_ci_m <- dat_ci_m %>%
  st_as_sf(coords = c("long", "lat"),crs = 4326) %>%
  dplyr::mutate(lat = sf::st_coordinates(.)[,2],
                long = sf::st_coordinates(.)[,1]) %>%
  mutate(region = toupper(region))


if(gender=="female"){
  dat_adm <- st_read(here::here("union", "DHS_adm.shp"))
  dat_adm <- st_make_valid(dat_adm) %>%
    mutate(region = REGNAME) 
  dat_ci_adm <- st_join(dat_adm, dat_ci_f) 
} else if(gender=="male"){
  dat_adm <- st_read(here::here("union", "DHS_adm.shp"))
  dat_adm <- st_make_valid(dat_adm) %>%
    mutate(region = REGNAME) 
  dat_ci_adm <- st_join(dat_adm, dat_ci_m) 
} else if(gender=="all"){
  dat_adm <- st_read(here::here("union", "DHS_adm.shp"))
  dat_adm <- st_make_valid(dat_adm) %>%
    mutate(region = REGNAME) 
  dat_ci_adm <- st_join(dat_adm, dat_ci_all) 
  
  
}

dat_ci_adm <-dat_ci_adm[!duplicated(dat_ci_adm$REG_ID),]

################

# mapping SII

sPDF <- world %>%
  filter(continent == "Africa")

# female

fsii=ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_ci_adm, #%>%
            #filter(!is.na(sii)), #%>%
          aes(fill = sii),
          size = 0) +
  geom_sf(data = sPDF, fill = NA) +
  scale_fill_continuous_diverging(palette = "Blue-Red 3") +
  geom_sf(data = sPDF, fill = NA) +
  labs(fill = paste0("SII, ",gender), tag = "A") +
  theme_bw()
fsii <- fsii + theme(legend.text = element_text(size=8), 
                     plot.tag = element_text(face="bold")) 


frii=ggplot() +
  geom_sf(data = sPDF, fill = "grey") +
  geom_sf(data = dat_ci_adm, 
          aes(fill = rii),
          size = 0) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_fill_continuous_diverging(palette = "Blue-Red 3",trans = log10_trans()
  ) +
  geom_sf(data = sPDF, fill = NA) + 
  labs(fill = paste0("RII, ",gender), tag = "B") +
  theme_bw()
frii <- frii + theme(legend.text = element_text(size=8), 
                     plot.tag = element_text(face="bold")) 

# male 

msii=ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_ci_adm, #%>%
            #filter(!is.na(sii)), #%>%
          aes(fill = sii),
          size = 0) +
  geom_sf(data = sPDF, fill = NA) +
  scale_fill_continuous_diverging(palette = "Blue-Red 3") +
  geom_sf(data = sPDF, fill = NA) +
  labs(fill = paste0("SII, ",gender), tag = "C") +
  theme_bw()
msii <- msii + theme(legend.text = element_text(size=8), 
                     plot.tag = element_text(face="bold")) 

mrii=ggplot() +
  geom_sf(data = sPDF, fill = "grey") +
  geom_sf(data = dat_ci_adm, 
          aes(fill = rii),
          size = 0) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_fill_continuous_diverging(palette = "Blue-Red 3",trans = log10_trans()
  ) +
  geom_sf(data = sPDF, fill = NA) + 
  labs(fill = paste0("RII, ",gender), tag = "D") +
  theme_bw()
mrii <- mrii + theme(legend.text = element_text(size=8), 
                     plot.tag = element_text(face="bold")) 


###########
# combining by gender, region

freg_sii_rii <- plot_grid(fsii, frii,ncol = 2)

mreg_sii_rii <- plot_grid(msii, mrii,ncol = 2)

fm_reg_sii_rii <- plot_grid(freg_sii_rii, mreg_sii_rii, nrow=2)


ggsave(here::here(paste0("fm_reg_sii_rii.pdf")),
       units = "cm", width = 30, height = 20, dpi = 1200, fm_reg_sii_rii)


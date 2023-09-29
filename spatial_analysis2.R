
rm(list = ls())
library(tidyverse)
library(Rcpp)
library(spdep)
library(raster)
library(tmap)
library(sf)
library(leaflet)
library(dplyr)
library(mapview)
library(spData)
library(cowplot)
library(colorspace)

gender = "female"
gender = "male"
gender = "all"

iso = "CM"
iso = "CD"
iso = "CI"
iso = "AO"
iso = "BI"
iso = "ET"
iso = "GA"
iso = "GH"
iso = "GN"
iso = "LB"
iso = "LS"
iso = "ML"
iso = "MW"
iso = "MZ"
iso = "NM"
iso = "SL"
iso = "SN"
iso = "RW"
iso = "TD"
iso = "TG"
iso = "TZ"
iso = "UG"
iso = "ZA"
iso = "ZM"
iso = "ZW"

# region
dat_ci_map <- readRDS(paste0(here::here("dat_sii-rii_linear_map_latest_v4"),gender,".rds")) %>%
  #readRDS(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/dat_sii-rii_linear_map_latest2",gender,".rds")) %>%
  #dat_ci_map <- readRDS(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/dat_sii-rii_linear_map3_HIVprev",gender,".rds")) %>%
  dplyr::filter(sample.size>=10) %>%
  #dplyr::filter(sample.size>=20) %>%
  #dplyr::filter(sample.size>=30) %>%
  dplyr::filter(rii !=Inf & rii>0) %>%
  #dplyr::filter(country %in% c("AO","BU","CM","ET","GN","MW","MZ","SL","SN","ZA","ZM","ZW")) 
  dplyr::filter(country %in% c("CD","CI","GA","GH","LB","LS","ML","NA","RW","TD","TG","TZ","UG")) 
  #dplyr::filter(country %in% c(iso)) 
#dat_ci_map <- subset(dat_ci_map, gender=="female", select=c(colnames(dat_ci_map)))

#source("/Users/pearlanneante-testard/Dropbox/Pearl - PhD project/mapping/mal_ineq_fun.R") 

breaks <- c(-Inf, -2.58, -1.96, -1.65, 1.65, 1.96, 2.58, Inf)
labels <- c("Cold spot: 99% confidence", "Cold spot: 95% confidence", "Cold spot: 90% confidence", "PSU Not significant","Hot spot: 90% confidence", "Hot spot: 95% confidence", "Hot spot: 99% confidence")


### Wealth

#dat_ci_sp_f <- as(dat_ci_map_f, "Spatial")
dat_ci_map_2 <-dat_ci_map[!duplicated(dat_ci_map$ID),]
xy <- dat_ci_map_2[,c("long","lat")]
dat_ci_sp <- SpatialPointsDataFrame(coords = xy, data = dat_ci_map_2,
                                    proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
dat_ci_sp <- as(st_as_sf(dat_ci_sp), "Spatial")
#dat_ci_sp_2 <-dat_ci_sp[!duplicated(dat_ci_sp$ID),]
coords<-coordinates(dat_ci_sp)
IDs<-row.names(as.data.frame(dat_ci_sp))

system.time(Neigh_nb<-knn2nb(knearneigh(coords, k=1, longlat = TRUE), row.names=IDs)) #knn2nb function converts a knn object returned by knearneigh into a neighbours list of class nb with a list of integer vectors containing neighbour region number ids.
dsts<-unlist(nbdists(Neigh_nb,coords)) #returns the Euclidean distances along the links in a list of the same form as the
max_1nn<-max(dsts) #max distance
Neigh_kd1<-dnearneigh(coords,d1=0, d2=max_1nn, row.names=IDs) #identifies neighbors of region points by Euclidean distance between lower and upper
self<-nb2listw(Neigh_kd1, style="W", zero.policy = T) #supple

w.getis <- dat_ci_sp %>% as.data.frame() %>%
  #st_set_geometry(NULL) %>%
  nest(data = everything()) %>%
  crossing(var1 = c("sii", "rii")) %>%
  mutate(data_sub = map2(.x = data, .y = var1, .f = ~dplyr::select(.x, .y, country, cluster_number, 
                                                                   prev, lat, long) %>% 
                           rename(output = 1)),
         LISA = map(.x = data_sub, .f = ~localG(.x$output, self)),
         clust_LISA = map(.x = LISA, .f = ~cut(.x, include.lowest = TRUE, 
                                               breaks = breaks, labels = labels))) %>%
  dplyr::select(data_sub, var1, LISA, clust_LISA) %>%
  unnest(cols = c(data_sub, LISA, clust_LISA))

#saveRDS(w.getis, paste0(here::here("dat_localG_cluster"),gender,".rds")) 
#write_csv(w.getis, paste0(here::here("dat_localG_cluster"),gender,".csv")) 

dat_w_map_sii <- w.getis %>%
  dplyr::filter(var1 =="sii") %>%
  st_as_sf(coords = c("long", "lat"),crs = 4326) %>%
  dplyr::filter(!is.na(LISA)) %>%
  dplyr::mutate(lat = sf::st_coordinates(.)[,2],
                long = sf::st_coordinates(.)[,1])

sPDF <- world %>%
  filter(continent == "Africa")# %>%
  #filter(iso_a2 == iso)


f3_a <- ggplot() +
  geom_sf(data = sPDF, fill = "darkgrey") + 
  geom_sf(data = dat_w_map_sii,
          aes(col = clust_LISA), size = 0.7, alpha =1) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_color_discrete_diverging(palette = "Green-Brown",
                                 limits = labels) +
  #labs(color = paste0("LISA - SII, "," ",gender), tag = "A") +
  #labs(color = paste0(iso,", LISA - SII, "," ",gender)) +
  #labs(color = paste0("LISA - SII (PSU>=20), "," ",gender), tag = "A") +
  #labs(color = paste0("LISA - SII (PSU>=30), "," ",gender), tag = "A") +
  #labs(color = paste0("LISA - SII (2015-19),"," ",gender), tag = "A") +
  labs(color = paste0("LISA - SII (2011-14),"," ",gender), tag = "A") +
  #labs(color = paste0("LISA - SII (HIV prev),"," ",gender)) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme_minimal()
f3_a <- f3_a + theme(legend.text = element_text(size=8), 
                     plot.tag = element_text(face="bold")) 

m3_a <- ggplot() +
  geom_sf(data = sPDF, fill = "darkgrey") + 
  geom_sf(data = dat_w_map_sii,
          aes(col = clust_LISA), size = 0.7, alpha =1) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_color_discrete_diverging(palette = "Green-Brown",
                                 limits = labels) +
  #labs(color = paste0("LISA - SII, "," ",gender), tag = "C") +
  #labs(color = paste0("LISA - SII (2015-19),"," ",gender), tag = "C") +
  labs(color = paste0("LISA - SII (2011-14),"," ",gender), tag = "C") +
  #labs(color = paste0(iso,", LISA - SII, "," ",gender)) +
  #labs(color = paste0("LISA - SII (PSU>=20), "," ",gender), tag = "C") +
  #labs(color = paste0("LISA - SII (PSU>=30), "," ",gender), tag = "C") +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme_minimal()
m3_a <- m3_a + theme(legend.text = element_text(size=8), 
                     plot.tag = element_text(face="bold")) 

#######
#######
#######

dat_w_map_rii <- w.getis %>%
  filter(var1 =="rii") %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  dplyr::filter(!is.na(LISA)) %>%
  dplyr::mutate(lat = sf::st_coordinates(.)[,2],
                long = sf::st_coordinates(.)[,1])

#pal <- colorFactor("RdBu", dat_w_map_sii$clust_LISA, reverse = T)

f3_b <- ggplot() +
  geom_sf(data = sPDF, fill = "darkgrey") + 
  geom_sf(data = dat_w_map_rii,
          aes(col = clust_LISA), size = .7, alpha =1) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_color_discrete_diverging(palette = "Green-Brown",
                                 limits = labels) +
  #labs(color = paste0("LISA - RII, "," ",gender), tag = "B") +
  #labs(color = paste0("LISA - RII (2015-19),"," ",gender), tag = "B") +
  labs(color = paste0("LISA - RII (2011-14),"," ",gender), tag = "B") +
  #labs(color = paste0(iso,", LISA - RII, "," ",gender)) +
  #labs(color = paste0("LISA - RII (PSU>=20), "," ",gender), tag = "B") +
  #labs(color = paste0("LISA - RII (PSU>=30), "," ",gender), tag = "B") +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme_minimal()
f3_b <- f3_b + theme(legend.text = element_text(size=8), 
                     plot.tag = element_text(face="bold")) 

m3_b <- ggplot() +
  geom_sf(data = sPDF, fill = "darkgrey") + 
  geom_sf(data = dat_w_map_rii,
          aes(col = clust_LISA), size = .7, alpha =1) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_color_discrete_diverging(palette = "Green-Brown", 
                                 limits = labels) +
  #labs(color = paste0("LISA - RII, "," ",gender), tag = "D") +
  #labs(color = paste0("LISA - RII (2015-19),"," ",gender), tag = "D") +
  labs(color = paste0("LISA - RII (2011-14),"," ",gender), tag = "D") +
  #labs(color = paste0(iso,", LISA - RII, "," ",gender)) +
  #labs(color = paste0("LISA - RII (PSU>=20), "," ",gender), tag = "D") +
  #labs(color = paste0("LISA - RII (PSU>=30), "," ",gender), tag = "D") +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme_minimal()
m3_b <- m3_b + theme(legend.text = element_text(size=8), 
                     plot.tag = element_text(face="bold")) 


###### by gender
spf <- plot_grid(f3_a, f3_b, ncol = 2,  align = "hv", rel_widths = c(1, 1))

spm <- plot_grid(m3_a, m3_b, ncol = 2, align = "hv", rel_widths = c(1, 1))

sp_fm <- plot_grid(spf, spm, align = "hv", nrow = 2)

#library(egg)      
#ggarrange(spf,
        #  spm,
        #  ncol = 1)

ind="sii-rii"
#ggsave(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/res/",ind,"spatial_latest2.jpg"),
       #units = "cm", width = 40, height = 20, dpi = 650, sp_fm)

ggsave(here::here(paste0(ind,"spatial_latest4.pdf")),
       units = "cm", width = 30, height = 20, dpi = 1200, sp_fm)

ggsave(here::here(paste0(ind,"spatial_latest4_atleast20.pdf")),
       units = "cm", width = 30, height = 20, dpi = 1200, sp_fm)

ggsave(here::here(paste0(ind,"spatial_latest4_atleast30.pdf")),
       units = "cm", width = 30, height = 20, dpi = 1200, sp_fm)

ggsave(here::here(paste0(ind,"spatial_latest4_2015_2019.pdf")),
       units = "cm", width = 30, height = 20, dpi = 1200, sp_fm)

ggsave(here::here(paste0(ind,"spatial_latest4_2011_2014.pdf")),
       units = "cm", width = 30, height = 20, dpi = 1200, sp_fm)

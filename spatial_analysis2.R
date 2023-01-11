
################
## This script conducts Getis-Ord G* statistics and outputs maps
## Author: Pearl Ante-Testard
## email: pearlannemante@gmail.com/ pearl.ante@ucsf.edu
###############

# cleans environment
rm(list = ls())


# packages
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
library(spData)

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
#dat_ci_map <- readRDS(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/dat_sii-rii_linear_map_latest2",gender,".rds")) %>%
  #dat_ci_map <- readRDS(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/dat_sii-rii_linear_map3_HIVprev",gender,".rds")) %>%
dat_ci_map <- readRDS(file = here(paste0("dat_sii-rii_linear_map_latest2", gender, ".rds"))) %>%  
  dplyr::filter(sample.size>=10) %>%
  dplyr::filter(rii !=Inf & rii>0) %>%
  #dplyr::filter(country %in% c("AO","BI","CM","ET","GN","LB","ML","MW","MZ","SL","SN","ZA","ZM","ZW")) 
  dplyr::filter(country %in% c("CD","CI","GA","GH","LS","NM","RW","TD","TG","TZ","UG")) 
  #dplyr::filter(country %in% c(iso)) 
#dat_ci_map <- subset(dat_ci_map, gender=="female", select=c(colnames(dat_ci_map)))

#source("/Users/pearlanneante-testard/Dropbox/Pearl - PhD project/mapping/mal_ineq_fun.R") 

breaks <- c(-Inf, -2.58, -1.96, -1.65, 1.65, 1.96, 2.58, Inf)
labels <- c("Cold spot: 99% confidence", "Cold spot: 95% confidence", "Cold spot: 90% confidence", "PSU Not significant","Hot spot: 90% confidence", "Hot spot: 95% confidence", "Hot spot: 99% confidence")


### Wealth

#dat_ci_sp_f <- as(dat_ci_map_f, "Spatial")
dat_ci_sp <- as(dat_ci_map, "Spatial")

coords<-coordinates(dat_ci_sp)  # retrieving spatial coordinates
IDs<-row.names(as.data.frame(coords)) # row names
Neigh_nb<-knn2nb(knearneigh(coords, k=2, longlat = TRUE), row.names=IDs) #knn2nb function converts a knn object returned by 
# knearneigh into a neighbours list of class nb with a list of integer vectors containing neighbour region number ids.
# knearneigh returns a matrix with indices of points belonging to the set of the k nearest neighbor of each other
# k = number of nearest neighbors returned
dsts<-unlist(nbdists(Neigh_nb,coords)) #returns the Euclidean distances along the links in a list of the same form as the
#neighbors list
max_1nn<-max(dsts) # max distance
Neigh_kd1<-dnearneigh(coords,d1=0, d2=max_1nn, row.names=IDs) # identifies neighbors of region points by Euclidean distance 
# between lower and upper bounds, longlat=TRUE, by Great Circle distance in kilometersNeighbour list object:
# Number of regions: 5080 
# Number of nonzero links: 1736086 
# Percentage nonzero weights: 6.727347 
# Average number of links: 341.7492 

self<-nb2listw(Neigh_kd1, style="W", zero.policy = T) #supplements a neighbours list with spatial weights for the chosen coding scheme
#Characteristics of weights list object:
# Neighbour list object:
# Number of regions: 5080 
# Number of nonzero links: 1736086 
# Percentage nonzero weights: 6.727347 
# Average number of links: 341.7492 
# Weights style: W 
# Weights constants summary:
# n       nn   S0       S1       S2
# W 5080 25806400 5080 81.52199 20459.81

w.getis <- dat_ci_map %>%
  st_set_geometry(NULL) %>%
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

###
# A LISA (Anselin, 1995) is seen as having two important characteristics. First, it provides a statistic for each location with an assessment of significance. 
# Second, it establishes a proportional relationship between the sum of the local statistics and a corresponding global statistic.
# ref: https://geodacenter.github.io/workbook/6a_local_auto/lab6a.html#lisa-principle


dat_w_map_sii <- w.getis %>%
  filter(var1 =="sii") %>%
  st_as_sf(coords = c("long", "lat"),crs = 4326) %>%
  filter(!is.na(LISA)) %>%
  dplyr::mutate(lat = sf::st_coordinates(.)[,2],
                long = sf::st_coordinates(.)[,1])

#pal <- colorFactor("RdBu", dat_w_map_sii$clust_LISA, reverse = T)
#sPDF <- world %>%
  #filter(continent == "Africa")

#library(rgdal)
#my_spdf <- readOGR( 
  #dsn= "/Users/pearlanneante-testard/Dropbox/PhD/mapping/world_shape_file/TM_WORLD_BORDERS_SIMPL-0.3.shp") 
#layer=("TM_WORLD_BORDERS_SIMPL-0.3",verbose=FALSE)
# Select Africa only
#africa <- my_spdf[my_spdf@data$REGION==2 , ]
#map <- st_as_sf(africa)

sPDF <- world %>%
  filter(continent == "Africa") #%>%
  #filter(iso_a2 == iso)

### Plots SII

f3_a <- ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_w_map_sii,
          aes(col = clust_LISA), size = 0.7, alpha =1) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_color_discrete_diverging(palette = "Green-Brown",
                                 limits = labels) +
  #labs(color = paste0(iso,", LISA - SII, "," ",gender)) +
  labs(color = paste0("LISA - SII, "," ",gender)) +
  #labs(color = paste0("LISA - SII (PSU>=20), "," ",gender)) +
  #labs(color = paste0("LISA - SII (PSU>=30), "," ",gender)) +
  #labs(color = paste0("LISA - SII (2015-19),"," ",gender)) +
  #labs(color = paste0("LISA - SII (2011-14),"," ",gender)) +
  #labs(color = paste0("LISA - SII (HIV prev),"," ",gender)) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme_bw()

m3_a <- ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_w_map_sii,
          aes(col = clust_LISA), size = 0.7, alpha =1) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_color_discrete_diverging(palette = "Green-Brown",
                                 limits = labels) +
  #labs(color = paste0("LISA - SII (2015-19),"," ",gender)) +
  #labs(color = paste0("LISA - SII (2011-14),"," ",gender)) +
  #labs(color = paste0(iso,", LISA - SII, "," ",gender)) +
  labs(color = paste0("LISA - SII, "," ",gender)) +
  #labs(color = paste0("LISA - SII (PSU>=20), "," ",gender)) +
  #labs(color = paste0("LISA - SII (PSU>=30), "," ",gender)) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme_bw()

# Merging plots
sp_w_f <- plot_grid(f3_a, m3_a, ncol = 2, 
                    labels = c("A","B"))

ind="sii"
ggsave(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/res/",ind,"spatialnoBF_2.jpg"),
       units = "cm", width = 35, height = 17, dpi = 650, sp_w_f)


### Plots RII

dat_w_map_rii <- w.getis %>%
  filter(var1 =="rii") %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  filter(!is.na(LISA)) %>%
  dplyr::mutate(lat = sf::st_coordinates(.)[,2],
                long = sf::st_coordinates(.)[,1])

#pal <- colorFactor("RdBu", dat_w_map_sii$clust_LISA, reverse = T)

f3_b <- ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_w_map_rii,
          aes(col = clust_LISA), size = .7, alpha =1) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_color_discrete_diverging(palette = "Green-Brown",
                                 limits = labels) +
  #labs(color = paste0("LISA - RII (2015-19),"," ",gender)) +
  #labs(color = paste0("LISA - RII (2011-14),"," ",gender)) +
  #labs(color = paste0(iso,", LISA - RII, "," ",gender)) +
  labs(color = paste0("LISA - RII, "," ",gender)) +
  #labs(color = paste0("LISA - RII (PSU>=20), "," ",gender)) +
  #labs(color = paste0("LISA - RII (PSU>=30), "," ",gender)) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme_bw()

m3_b <- ggplot() +
  geom_sf(data = sPDF, fill = "grey") + 
  geom_sf(data = dat_w_map_rii,
          aes(col = clust_LISA), size = .7, alpha =1) +
  geom_sf(data = sPDF, fill = NA) + 
  scale_color_discrete_diverging(palette = "Green-Brown", 
                                 limits = labels) +
  #labs(color = paste0("LISA - RII (2015-19),"," ",gender)) +
  #labs(color = paste0("LISA - RII (2011-14),"," ",gender)) +
  #labs(color = paste0(iso,", LISA - RII, "," ",gender)) +
  labs(color = paste0("LISA - RII, "," ",gender)) +
  #labs(color = paste0("LISA - RII (PSU>=20), "," ",gender)) +
  #labs(color = paste0("LISA - RII (PSU>=30), "," ",gender)) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme_bw()


# Merging plots
sp_w <- plot_grid(f3_b, m3_b, ncol = 2, 
                    labels = c("C","D"))
ind="rii"

ggsave(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/res/",ind,"spatial_2.jpg"),
       units = "cm", width = 35, height = 17, dpi = 650, sp_w)


### Merging SII and RII plots

sp_w_mf <- plot_grid(sp_w_f, sp_w, nrow = 2)


### Saves plots

ind="sii-rii"
#ggsave(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/res/",ind,"spatial_2011-14.jpg"),
      # units = "cm", width = 35, height = 17, dpi = 650, sp_w_mf)

ggsave(file = here("hotspot_mapping.pdf"),
       units = "cm", width = 35, height = 17, dpi = 650, sp_w_mf)

#ggsave(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/res/",ind,"spatial_HIVprev.jpg"),
      # units = "cm", width = 35, height = 17, dpi = 650, sp_w_mf)

###### by gender
spf <- plot_grid(f3_a, f3_b, ncol = 2,  align = "hv",
                  labels = c("A","B"), rel_widths = c(1, 1))

spm <- plot_grid(m3_a, m3_b, ncol = 2, align = "hv",
                 labels = c("C","D"), rel_widths = c(1, 1))

sp_fm <- plot_grid(spf, spm, align = "hv", nrow = 2, rel_widths = c(1, 1))

#library(egg)      
#ggarrange(spf,
        #  spm,
        #  ncol = 1)

ind="sii-rii"
ggsave(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/res/",ind,"spatial_latest2.jpg"),
       units = "cm", width = 40, height = 20, dpi = 650, sp_fm)

#ggsave(file = here("hotspot_mapping_bygender_2011-14.pdf"),
       #units = "cm", width = 40, height = 20, dpi = 650, sp_fm)

#ggsave(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/res/",ind,iso,"spatial_latest3.jpg"),
      # units = "cm", width = 40, height = 20, dpi = 650, sp_fm)

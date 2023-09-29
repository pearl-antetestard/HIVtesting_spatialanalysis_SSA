
rm(list = ls())

library( RColorBrewer)
library(ggpubr)
library(cowplot)
library(tidyverse)
library(sf)

gender="female"
#dat_ci_map_f=readRDS(paste0(here::here("dat_sii-rii_linear_map_latest_v4"),gender,".rds")) #cluster
#dat_ci_map_f=readRDS(paste0(here::here("dat_sii-rii_region_map_latest_v4"),gender,".rds"))  #region
dat_ci_map_f=readRDS(paste0(here::here("dat_sii-rii_national_map_latest_v4"),gender,".rds")) #national
dat_ci_map_f$Gender = "Female"

# cluster
dat_ci_map_f_2 <-dat_ci_map_f[!duplicated(dat_ci_map_f$ID),]
dat_ci_map_f_2 <- dat_ci_map_f_2 %>%
  dplyr::select(country,region,cluster_number,ID,year,lat,long,dat_s,sample.size.N, sample.size.m, sample.size, 
                prev, HIVprev, wealth_cats, sii,rii,Gender) %>%
  unnest(cols = c(country,region,cluster_number,ID,year,lat,long,dat_s,sample.size.N, sample.size.m, sample.size, 
                  prev, HIVprev, wealth_cats, sii,rii,Gender)) %>%
  ungroup() 

# regional
dat_ci_map_f <- dat_ci_map_f %>%
  st_as_sf(coords = c("long", "lat"),crs = 4326) %>%
  dplyr::mutate(lat = sf::st_coordinates(.)[,2],
                long = sf::st_coordinates(.)[,1]) %>%
  mutate(region = toupper(region))

if(gender=="female"){
  dat_adm <- st_read(here::here("union", "DHS_adm.shp"))
  dat_adm <- st_make_valid(dat_adm) %>%
    mutate(region = REGNAME) 
  dat_ci_adm <- st_join(dat_adm, dat_ci_map_f) 
  
  dat_ci_map_f_2 <-dat_ci_adm[!duplicated(dat_ci_adm$REG_ID),] %>% #region
    filter(country.x != "BF")
  
} else if(gender=="male"){
  dat_adm <- st_read(here::here("union", "DHS_adm.shp"))
  dat_adm <- st_make_valid(dat_adm) %>%
    mutate(region = REGNAME) 
  dat_ci_adm <- st_join(dat_adm, dat_ci_map_f) 
  
  dat_ci_map_m_2 <-dat_ci_adm[!duplicated(dat_ci_adm$REG_ID),] %>% #region
    filter(country.x != "BF")
} 

# national
dat_ci_map_f_2 <- dat_ci_map_f %>%
  dplyr::select(country,dat_s,sample.size.N, sample.size.m, sample.size, 
                prev, HIVprev, wealth_cats, sii,rii,Gender) %>%
  unnest(cols = c(country,dat_s,sample.size.N, sample.size.m, sample.size, 
                  prev, HIVprev, wealth_cats, sii,rii,Gender)) %>%
  ungroup() 

gender="male"
#dat_ci_map_m=readRDS(paste0(here::here("dat_sii-rii_linear_map_latest_v4"),gender,".rds")) #cluster
#dat_ci_map_m=readRDS(paste0(here::here("dat_sii-rii_region_map_latest_v4"),gender,".rds")) #region
dat_ci_map_m=readRDS(paste0(here::here("dat_sii-rii_national_map_latest_v4"),gender,".rds")) #national
dat_ci_map_m$Gender = "Male"

# cluster
dat_ci_map_m_2 <-dat_ci_map_m[!duplicated(dat_ci_map_m$ID),]
dat_ci_map_m_2 <- dat_ci_map_m_2 %>%
  dplyr::select(country,region,cluster_number,ID,year,lat,long,dat_s,sample.size.N, sample.size.m, sample.size, 
                prev, HIVprev, wealth_cats, sii,rii,Gender) %>%
  unnest(cols = c(country,region,cluster_number,ID,year,lat,long,dat_s,sample.size.N, sample.size.m, sample.size, 
                  prev, HIVprev, wealth_cats, sii,rii,Gender)) %>%
  ungroup() 

# regional
dat_ci_map_m <- dat_ci_map_m %>%
  st_as_sf(coords = c("long", "lat"),crs = 4326) %>%
  dplyr::mutate(lat = sf::st_coordinates(.)[,2],
                long = sf::st_coordinates(.)[,1]) %>%
  mutate(region = toupper(region))

if(gender=="female"){
  dat_adm <- st_read(here::here("union", "DHS_adm.shp"))
  dat_adm <- st_make_valid(dat_adm) %>%
    mutate(region = REGNAME) 
  dat_ci_adm <- st_join(dat_adm, dat_ci_map_f) 
  
  dat_ci_map_f_2 <-dat_ci_adm[!duplicated(dat_ci_adm$REG_ID),] %>% #region
    filter(country.x != "BF")
  
} else if(gender=="male"){
  dat_adm <- st_read(here::here("union", "DHS_adm.shp"))
  dat_adm <- st_make_valid(dat_adm) %>%
    mutate(region = REGNAME) 
  dat_ci_adm <- st_join(dat_adm, dat_ci_map_m) 
  
  dat_ci_map_m_2 <-dat_ci_adm[!duplicated(dat_ci_adm$REG_ID),] %>% #region
    filter(country.x != "BF")
} 


# national
dat_ci_map_m_2 <- dat_ci_map_m %>%
  dplyr::select(country,dat_s,sample.size.N, sample.size.m, sample.size, 
                prev, HIVprev, wealth_cats, sii,rii,Gender) %>%
  unnest(cols = c(country,dat_s,sample.size.N, sample.size.m, sample.size, 
                  prev, HIVprev, wealth_cats, sii,rii,Gender)) %>%
  ungroup() 


dat_ci_map = rbind(dat_ci_map_f_2,dat_ci_map_m_2) %>%
             mutate(HIVtestprev = prev*100,
                   HIVprev = HIVprev*100) %>%
             filter(!is.na(HIVtestprev))

#write_csv2(dat_ci_map, "/Users/pearlante/Dropbox/PhD/mapping/dat_hivtest_prev_national.csv")
#write_csv2(dat_ci_map, "/Users/pearlante/Dropbox/PhD/mapping/dat_hivtest_prev.csv")

#colourC = 27L
#getPalette = colorRampPalette(brewer.pal(8, "Set2"))

options(scipen=999)

# Custom formatting function
format_pval <- function(pval){
  pval <- scales::pvalue(pval, accuracy= 0.001, add_p = TRUE)
  gsub(pattern = "(<)", replacement = " \\1 ", x = pval)
}

s=ggscatter(dat_ci_map,
          x= "HIVprev", y= "HIVtestprev",
          add = "reg.line",                         # Add regression line
          #conf.int = TRUE,                          # Add confidence interval
          col="Gender",repel = F, fullrange = F, size=2,
          palette = c("orange", "darkblue"),
          #facet.by = c("country","wealthindex"), # cluster and regional
          #facet.by = "country",
          #facet.by = "country.x", # region
          legend = "right", 
          xlab = "National-level HIV prevalence (%, weighted)",
          ylab = "National-level self-reported HIV testing uptake in the previous 12 months (%, weighted)",
          #xlab = "Province-level HIV prevalence (%, weighted)",
          #ylab = "Province-level self-reported HIV testing uptake in the previous 12 months (%, weighted)",
          #xlab = "PSU-level HIV prevalence (%, weighted)",
          #ylab = "PSU-level self-reported HIV testing uptake in the previous 12 months (%, weighted)",
          #cor.method = "pearson", 
          #add.params = list(linetype = "Gender"), 
          cor.method = "spearman"
          #font.label = c(7, "plain"),
          ) +#scales="free_x")+
          #scale_colour_manual(values = getPalette(colourC))+
  #stat_cor(aes(col = Gender,label = ..p.label.., ),p.accuracy = 0.001,digits=3,size = 5)
  stat_cor(aes(col = Gender,label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001,r.accuracy = 0.001,digits=3,
           size = 5, # country
           #size = 2.5, # region, cluster
           #size = 3, # cluster x wealth
           label.x.npc = 0.5, label.y.npc = 0.40 # country
           #label.x.npc = 0.35, label.y.npc = 0.2 # region
           #label.x.npc = 0.47, label.y.npc = 0.18 # cluster
           #label.x.npc = 0.53, label.y.npc = 0.4 # cluster x wealth
           ) 

sl=s+ #theme(legend.position = "none")+
  ggtitle("Correlation between national-level HIV prevalence and self-reported recent uptake of HIV testing")
  #ggtitle("Correlation between province-level HIV prevalence and self-reported recent uptake of HIV testing")
  #ggtitle("Correlation between PSU-level HIV prevalence and self-reported recent uptake of HIV testing")
  #ggtitle("Correlation between PSU-level HIV prevalence and self-reported recent uptake of HIV testing by wealth")
sl
#msl=s+theme(legend.position = "none")+
  #gtitle("Relationship between PSU HIV prevalence and PSU self-reported recent (< 12 months) uptake of HIV testing among men")
#fm_sp <-plot_grid(fsl, sl, ncol = 1, 
                  #labels = c("Female","Male"),rel_widths = c(1, 1))

ggsave(here::here("ggscatter_HIVprev_test_spearman_cluster_wealth.pdf"),
          units = "cm", width = 40, height = 60, dpi = 1200, sl)
ggsave(here::here("ggscatter_HIVprev_test_spearman_cluster.pdf"),
       units = "cm", width = 30, height = 20, dpi = 1200, sl)

ggsave(here::here("ggscatter_HIVprev_test_spearman_region_wealth.pdf"),
       units = "cm", width = 40, height = 60, dpi = 1200, sl)
ggsave(here::here("ggscatter_HIVprev_test_spearman_region.pdf"),
       units = "cm", width = 30, height = 20, dpi = 1200, sl)


ggsave(here::here("ggscatter_HIVprev_test_spearman_national.pdf"),
       units = "cm",width = 30, height = 20, dpi = 1200, sl)
#ggsave(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/res/ggscatter_HIVprev_test",gender,".jpg"),
       #units = "cm", width = 30, height = 19, dpi = 600, msl)

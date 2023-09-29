
rm(list = ls())

library(sf)
library(spdep)
library(tmap)

gender = "female"
gender = "male"
gender = "all"

options(scipen=999)

# region
dat_ci_map <- readRDS(paste0(here::here("dat_sii-rii_linear_map_latest_v4"),gender,".rds")) %>%
  #readRDS(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/dat_sii-rii_linear_map_latest2",gender,".rds")) %>%
  #dat_ci_map <- readRDS(paste0("/Users/pearlanneante-testard/Dropbox/PhD/mapping/dat_sii-rii_linear_map3_HIVprev",gender,".rds")) %>%
  dplyr::filter(sample.size>=10) %>%
  dplyr::filter(rii !=Inf & rii>0) #%>%
  #dplyr::filter(country %in% c("AO","BI","CM","ET","GN","LB","ML","MW","MZ","SL","SN","ZA","ZM","ZW")) 
  #dplyr::filter(country %in% c("CD","CI","GA","GH","LS","NM","RW","TD","TG","TZ","UG")) 
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

###hypothesis: the sii values are randomly distributed across regions following a completely random process”

system.time(Neigh_nb<-knn2nb(knearneigh(coords, k=1, longlat = TRUE), row.names=IDs)) #knn2nb function converts a knn object returned by knearneigh into a neighbours list of class nb with a list of integer vectors containing neighbour region number ids.
dsts<-unlist(nbdists(Neigh_nb,coords)) #returns the Euclidean distances along the links in a list of the same form as the
max_1nn<-max(dsts) #max distance
Neigh_kd1<-dnearneigh(coords,d1=0, d2=max_1nn, row.names=IDs) #identifies neighbors of region points by Euclidean distance between lower and upper
self<-nb2listw(Neigh_kd1, style="W", zero.policy = T) #supple

#self$weights[1]

k=10

if(k==1){
  
  Neigh_nb<-knn2nb(knearneigh(coords, k=1, longlat = TRUE), row.names=IDs) #knn2nb function converts a knn object returned by knearneigh into a neighbours list of class nb with a list of integer vectors containing neighbour region number ids.
  dsts<-unlist(nbdists(Neigh_nb,coords)) #returns the Euclidean distances along the links in a list of the same form as the
  max_1nn<-max(dsts) #max distance
  Neigh_kd1<-dnearneigh(coords,d1=0, d2=max_1nn, row.names=IDs) #identifies neighbors of region points by Euclidean distance between lower and upper
  self<-nb2listw(Neigh_kd1, style="W", zero.policy = T) #supple
  
  
  k1=globalG.test(as.numeric(abs(dat_ci_sp$sii)),self,alternative="two.sided",zero.policy=TRUE)
  
  #k1=moran.test(as.numeric(dat_ci_sp_2$sii),self, alternative="two.sided")
  k_1=k1$estimate[1]
  
  #k1_rii=moran.test(as.numeric(dat_ci_sp$rii),self, alternative="two.sided")
  k1_rii=globalG.test(as.numeric(dat_ci_sp$rii),self,alternative="two.sided",zero.policy=TRUE)
  k_1_rii=k1_rii$estimate[1]
  
}else if(k==2) {
  
  Neigh_nb<-knn2nb(knearneigh(coords, k=2, longlat = TRUE), row.names=IDs) #knn2nb function converts a knn object returned by knearneigh into a neighbours list of class nb with a list of integer vectors containing neighbour region number ids.
  dsts<-unlist(nbdists(Neigh_nb,coords)) #returns the Euclidean distances along the links in a list of the same form as the
  max_1nn<-max(dsts) #max distance
  Neigh_kd1<-dnearneigh(coords,d1=0, d2=max_1nn, row.names=IDs) #identifies neighbors of region points by Euclidean distance between lower and upper
  self<-nb2listw(Neigh_kd1, style="W", zero.policy = T) #supple
  
  k2=globalG.test(as.numeric(abs(dat_ci_sp$sii)),self,alternative="two.sided",zero.policy=TRUE)
  k_2=k2$estimate[1]
  
  k2_rii=globalG.test(as.numeric(dat_ci_sp$rii),self,alternative="two.sided",zero.policy=TRUE)
  k_2_rii=k2_rii$estimate[1]
  
}else if(k==3) {
  
  Neigh_nb<-knn2nb(knearneigh(coords, k=3, longlat = TRUE), row.names=IDs) #knn2nb function converts a knn object returned by knearneigh into a neighbours list of class nb with a list of integer vectors containing neighbour region number ids.
  dsts<-unlist(nbdists(Neigh_nb,coords)) #returns the Euclidean distances along the links in a list of the same form as the
  max_1nn<-max(dsts) #max distance
  Neigh_kd1<-dnearneigh(coords,d1=0, d2=max_1nn, row.names=IDs) #identifies neighbors of region points by Euclidean distance between lower and upper
  self<-nb2listw(Neigh_kd1, style="W", zero.policy = T) #supple
  
  k3=globalG.test(as.numeric(abs(dat_ci_sp$sii)),self,alternative="two.sided",zero.policy=TRUE)
  k_3=k3$estimate[1]
  
  k3_rii=globalG.test(as.numeric(dat_ci_sp$rii),self,alternative="two.sided",zero.policy=TRUE)
  k_3_rii=k3_rii$estimate[1]
  
}else if(k==4) {
  
  Neigh_nb<-knn2nb(knearneigh(coords, k=4, longlat = TRUE), row.names=IDs) #knn2nb function converts a knn object returned by knearneigh into a neighbours list of class nb with a list of integer vectors containing neighbour region number ids.
  dsts<-unlist(nbdists(Neigh_nb,coords)) #returns the Euclidean distances along the links in a list of the same form as the
  max_1nn<-max(dsts) #max distance
  Neigh_kd1<-dnearneigh(coords,d1=0, d2=max_1nn, row.names=IDs) #identifies neighbors of region points by Euclidean distance between lower and upper
  self<-nb2listw(Neigh_kd1, style="W", zero.policy = T) #supple
  
  k4=globalG.test(as.numeric(abs(dat_ci_sp$sii)),self,alternative="two.sided",zero.policy=TRUE)
  k_4=k4$estimate[1]
  
  k4_rii=globalG.test(as.numeric(dat_ci_sp$rii),self,alternative="two.sided",zero.policy=TRUE)
  k_4_rii=k4_rii$estimate[1]
  
}else if(k==5) {
  
  Neigh_nb<-knn2nb(knearneigh(coords, k=5, longlat = TRUE), row.names=IDs) #knn2nb function converts a knn object returned by knearneigh into a neighbours list of class nb with a list of integer vectors containing neighbour region number ids.
  dsts<-unlist(nbdists(Neigh_nb,coords)) #returns the Euclidean distances along the links in a list of the same form as the
  max_1nn<-max(dsts) #max distance
  Neigh_kd1<-dnearneigh(coords,d1=0, d2=max_1nn, row.names=IDs) #identifies neighbors of region points by Euclidean distance between lower and upper
  self<-nb2listw(Neigh_kd1, style="W", zero.policy = T) #supple
  
  k5=globalG.test(as.numeric(abs(dat_ci_sp$sii)),self,alternative="two.sided",zero.policy=TRUE)
  k_5=k5$estimate[1]
  
  k5_rii=globalG.test(as.numeric(dat_ci_sp$rii),self,alternative="two.sided",zero.policy=TRUE)
  k_5_rii=k5_rii$estimate[1]
  
}else if(k==6) {
  
  Neigh_nb<-knn2nb(knearneigh(coords, k=6, longlat = TRUE), row.names=IDs) #knn2nb function converts a knn object returned by knearneigh into a neighbours list of class nb with a list of integer vectors containing neighbour region number ids.
  dsts<-unlist(nbdists(Neigh_nb,coords)) #returns the Euclidean distances along the links in a list of the same form as the
  max_1nn<-max(dsts) #max distance
  Neigh_kd1<-dnearneigh(coords,d1=0, d2=max_1nn, row.names=IDs) #identifies neighbors of region points by Euclidean distance between lower and upper
  self<-nb2listw(Neigh_kd1, style="W", zero.policy = T) #supple
  
  k6=globalG.test(as.numeric(abs(dat_ci_sp$sii)),self,alternative="two.sided",zero.policy=TRUE)
  k_6=k6$estimate[1]
  
  k6_rii=globalG.test(as.numeric(dat_ci_sp$rii),self,alternative="two.sided",zero.policy=TRUE)
  k_6_rii=k6_rii$estimate[1]
  
}else if(k==7) {
  
  Neigh_nb<-knn2nb(knearneigh(coords, k=7, longlat = TRUE), row.names=IDs) #knn2nb function converts a knn object returned by knearneigh into a neighbours list of class nb with a list of integer vectors containing neighbour region number ids.
  dsts<-unlist(nbdists(Neigh_nb,coords)) #returns the Euclidean distances along the links in a list of the same form as the
  max_1nn<-max(dsts) #max distance
  Neigh_kd1<-dnearneigh(coords,d1=0, d2=max_1nn, row.names=IDs) #identifies neighbors of region points by Euclidean distance between lower and upper
  self<-nb2listw(Neigh_kd1, style="W", zero.policy = T) #supple
  
  k7=globalG.test(as.numeric(abs(dat_ci_sp$sii)),self,alternative="two.sided",zero.policy=TRUE)
  k_7=k7$estimate[1]
  
  k7_rii=globalG.test(as.numeric(dat_ci_sp$rii),self,alternative="two.sided",zero.policy=TRUE)
  k_7_rii=k7_rii$estimate[1]
  
}else if(k==8) {
  
  Neigh_nb<-knn2nb(knearneigh(coords, k=8, longlat = TRUE), row.names=IDs) #knn2nb function converts a knn object returned by knearneigh into a neighbours list of class nb with a list of integer vectors containing neighbour region number ids.
  dsts<-unlist(nbdists(Neigh_nb,coords)) #returns the Euclidean distances along the links in a list of the same form as the
  max_1nn<-max(dsts) #max distance
  Neigh_kd1<-dnearneigh(coords,d1=0, d2=max_1nn, row.names=IDs) #identifies neighbors of region points by Euclidean distance between lower and upper
  self<-nb2listw(Neigh_kd1, style="W", zero.policy = T) #supple
  
  k8=globalG.test(as.numeric(abs(dat_ci_sp$sii)),self,alternative="two.sided",zero.policy=TRUE)
  k_8=k8$estimate[1]
  
  k8_rii=globalG.test(as.numeric(dat_ci_sp$rii),self,alternative="two.sided",zero.policy=TRUE)
  k_8_rii=k8_rii$estimate[1]
  
}else if(k==9) {
  
  Neigh_nb<-knn2nb(knearneigh(coords, k=9, longlat = TRUE), row.names=IDs) #knn2nb function converts a knn object returned by knearneigh into a neighbours list of class nb with a list of integer vectors containing neighbour region number ids.
  dsts<-unlist(nbdists(Neigh_nb,coords)) #returns the Euclidean distances along the links in a list of the same form as the
  max_1nn<-max(dsts) #max distance
  Neigh_kd1<-dnearneigh(coords,d1=0, d2=max_1nn, row.names=IDs) #identifies neighbors of region points by Euclidean distance between lower and upper
  self<-nb2listw(Neigh_kd1, style="W", zero.policy = T) #supple
  
  k9=globalG.test(as.numeric(abs(dat_ci_sp$sii)),self,alternative="two.sided",zero.policy=TRUE)
  k_9=k9$estimate[1]
  
  k9_rii=globalG.test(as.numeric(dat_ci_sp$rii),self,alternative="two.sided",zero.policy=TRUE)
  k_9_rii=k9_rii$estimate[1]
  
}else if(k==10) {
  
  Neigh_nb<-knn2nb(knearneigh(coords, k=10, longlat = TRUE), row.names=IDs) #knn2nb function converts a knn object returned by knearneigh into a neighbours list of class nb with a list of integer vectors containing neighbour region number ids.
  dsts<-unlist(nbdists(Neigh_nb,coords)) #returns the Euclidean distances along the links in a list of the same form as the
  max_1nn<-max(dsts) #max distance
  Neigh_kd1<-dnearneigh(coords,d1=0, d2=max_1nn, row.names=IDs) #identifies neighbors of region points by Euclidean distance between lower and upper
  self<-nb2listw(Neigh_kd1, style="W", zero.policy = T) #supple
  
  k10=globalG.test(as.numeric(abs(dat_ci_sp$sii)),self,alternative="two.sided",zero.policy=TRUE)
  k_10=k10$estimate[1]
  
  k10_rii=globalG.test(as.numeric(dat_ci_sp$rii),self,alternative="two.sided",zero.policy=TRUE)
  k_10_rii=k10_rii$estimate[1]
}

k_nearest_neighbor = c(1:10)
GlobalG_statistic = c(k_1, k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9, k_10)
GlobalG_statistic_rii = c(k_1_rii, k_2_rii, k_3_rii, k_4_rii, k_5_rii, k_6_rii, k_7_rii, 
                          k_8_rii, k_9_rii, k_10_rii)

ind="SII"
pdf(file=here::here(paste0("globalG",ind,gender,".pdf")))
plot(k_nearest_neighbor,GlobalG_statistic, xlab="Number of nearest neighbors",
     ylab="Global G statistic", main=paste0("Global G statistic per number of nearest neighbor",
                                            ", ",ind,", ",gender))
dev.off()

ind="RII"
pdf(file=here::here(paste0("globalG",ind,gender,".pdf")))
plot(k_nearest_neighbor,GlobalG_statistic_rii, xlab="Number of nearest neighbors",
       ylab="Global G statistic", main=paste0("Global G statistic per number of nearest neighbor",
                                               ", ",ind,", ",gender))
dev.off()


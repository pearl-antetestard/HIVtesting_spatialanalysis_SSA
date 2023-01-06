
# h_ineq ============================================================

## function to calculate RII and SII
## when computing country- and province-level RII-SII, select modified Poisson regression.
## when computing cluster-level RII-SII, select linear regression.

h_ineq <- function(dat, var_soc, var_outcome, clust, full = F, scale = 1) {

  soc <- enquo(var_soc)
  health <- enquo(var_outcome)
  cluster <- enquo(clust)
  
  # Master Table =============================================
  c_index.d <- dat %>% 
    ungroup() %>%
    mutate(group := !!soc) %>%
    group_by(group) %>%
    summarise(indiv=n(), prop_w=indiv/nrow(dat), cases=sum(!!health), prop_cases=cases/indiv) %>%
    mutate(cum_w=cumsum(indiv), pt=cum_w/sum(woman), cum_c=cumsum(cases), Lt=cum_c/sum(cases),
           con_i = (pt*lead(Lt, default=0))-(lead(pt, default=0)*Lt))
  
  # CI: FROM LOW TO HIGH SES =============================================
  c_index <- sum(c_index.d$con_i)
  
  # SII: FROM LOW TO HIGH SES =============================================
  ##### for HIV testing
  c_index.d2 <- dat %>% 
    ungroup() %>%
    mutate(group = !!soc) %>%
    group_by(group) %>%
    #mutate(group=factor(group),
           #group=fct_relevel(group,c("poorest","poorer","middle","richer","richest"))) %>% 
    #mutate(group=fct_relevel(group,c("poorest","poorer","middle","richer","richest")))  %>%
    summarise(indiv=n(), prop_w=indiv/nrow(dat), 
              testprevperSEP=mean(!!health,na.rm=T), #%>%
                                  test=(!!health),
                                  clusid=(!!cluster)) %>%
              #RANK = rank(!!soc)/indiv) %>%
    #dplyr::filter(testprevperSEP>0.000) %>%
    #arrange(group) %>%
    ungroup() %>%
    arrange(group) %>%
    mutate(testperSEP = testprevperSEP*scale,
      #logtest = log(testperSEP),
      cum_w=cumsum(as.numeric(indiv)), pt=cum_w/sum(indiv),#)
      #RANK = rank(group)/indiv)
      RANK = (pt + lag(pt, default = 0))/2)
     
  
  # SII: FROM HIGH TO LOW SES =============================================
  #### for HIV prevalence
  #c_index.d2 <- dat %>% 
    #ungroup() %>%
    #mutate(group = !!soc) %>%
    #group_by(group) %>%
    #mutate(group=factor(group)) %>% 
    #mutate(group=fct_relevel(group,c("poorest","poorer","middle","richer","richest")))  %>%
    #mutate(group=fct_relevel(group,c("richest","richer","middle","poorer","poorest")))  %>% #HIVprev
    #summarise(indiv=n(), prop_w=indiv/nrow(dat), #cases=sum(hivstatpos), prop_cases=cases/indiv) %>%
              #Tested=mean(!!health,na.rm=T)) %>%
    #arrange(desc(group)) %>%
    #mutate(Mort = Mort*scale, #NEW
           #tested=mean(HIVtest_12,na.rm=T),
           #cum_w=cumsum(indiv), pt=cum_w/sum(indiv), cum_c=cumsum(Mort), Lt=cum_c/sum(Mort),
           #con_i = (pt*lead(Lt, default=0))-(lead(pt, default=0)*Lt),
           #RANK = (pt + lag(pt, default = 0))/2, #NEW
           #cum_prop_cases = cum_c/cum_w #NEW
    #)
  
  library(geepack)
  lm1 <- geeglm(test~RANK, data = c_index.d2, family = "poisson"(link="log"),
                corstr = "exchangeable", id=clusid) #modified poisson for national and regional levels
  #library(sandwich) # to get robust estimators
  #library(lmtest) # to test coefficients
  #coeftest(lm1, vcov = sandwich)
  
  #lm1 <- lm(testperSEP~RANK, data = c_index.d2) 
  #https://data.library.virginia.edu/interpreting-log-transformations-in-a-linear-model/
  #the relationship of RANK and logtest become multiplicative since we are using the log of y
  
  #summary(lm1)
  
  #(SII <- coef(lm1)[[2]])
  (SII <- exp(coef(summary(lm1))[,"Estimate"][1]+coef(summary(lm1))[,"Estimate"][2]) - exp(coef(summary(lm1))[,"Estimate"][1]))
  
  #library(Zelig)
  #library(msm)
  
  #b <- coef(lm1)
  #v = vcov_gee(lm1)
  #SII_sd = deltamethod(~exp(x1+x2)-exp(x1),b,v)
  #(SII_low<-SII-1.96*SII_sd)
  #(SII_up<-SII+1.96*SII_sd)
  
  # RII: FROM LOW TO HIGH SES ============================================
  #p.low <- predict(lm1, data.frame(RANK=0), type="response")
  #p.low <- abs(p.low)
  #p.low=rescale(p.low, to=c(0,max(p.low)))
  #p.high <- predict(lm1, data.frame(RANK=1), type="response")
  #p.high <- abs(p.high)
  #p.high=rescale(p.high, to=c(0,max(p.high)))
  
  #if(p.low<0){
    #p.low = 0.001
  #}
  #if(p.high<0){
    #p.high = 0.001
  #}
  
  #(RII <- p.high/p.low)
  #(RII <- rescale(RII, to=c(0,max(RII))))
  
  #RII_sd = sd(RII)
  #(RII_low<-RII-1.96*RII_sd)
  #(RII_up<-RII+1.96*RII_sd)
  
  #b <- coef(lm1)
  #v = vcov(lm1)
  #form <- sprintf("~p.high/p.low")
  #RII_sd = deltamethod(as.formula(form),b,v)
  #(RII_low<-RII-1.96*RII_sd)
  #(RII_up<-RII+1.96*RII_sd)
  
  
  (RII <- exp(coef(summary(lm1))[,"Estimate"][2]))
  
  #beta_1 = coef(summary(lm1))[,"Estimate"][2]
  #se_b1 = coef(summary(lm1))[2,"Std.err"]
  #(RII_low<-exp(beta_1-1.96*se_b1))
  #(RII_up<-exp(beta_1+1.96*se_b1))
  
  # RII: FROM HIGH TO LOW SES =============================================
  #p.low <- predict(lm1, data.frame(RANK=1))
  #p.high <- predict(lm1, data.frame(RANK=0))
  
  #(RII <- p.low/p.high)
  
  # output =============================================
  dat_out <- data.frame(c_index = c_index, sii = SII, #sii_low = SII_low, sii_up = SII_up,
                        rii = RII)#, rii_low = RII_low, rii_up = RII_up)

  out <- dat_out

  if (isTRUE(full)) { out <- list(dat_out, c_index.d, c_index.d2)}

  return(out)
}



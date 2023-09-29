
rm(list = ls())

library(cowplot)
library(ggthemes)
library(ggExtra)

hist_ineq <- function(dat, var, lab) {
  
  var1 <- enquo(var)
  
  dat %>%
    ggplot(aes(x = !!var1)) +
    geom_histogram(fill = "white", col= "black") +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = lab) +
    theme_few()
}


bi_hist_ineq <- function(dat, var_x, var_y, lab_x, lab_y, trim.y.l = 0, trim.y.u = 1) {
  
  var1 <- enquo(var_x)
  var2 <- enquo(var_y)
  
  p <- dat %>%
    filter(quantile(!!var2, trim.y.l, na.rm = T)<=!!var2) %>%
    filter(quantile(!!var2, trim.y.u, nna.rm = T)>=!!var2) %>%
    ggplot(aes(x = !!var1, y = !!var2)) +
    geom_point(alpha = .2) +
    labs(x = lab_x, y = lab_y) +
    theme_few()
  
  ggExtra::ggMarginal(p, type = "histogram")
  
}

######################
gender = "female"

dat_ci_f=readRDS(paste0(here::here("dat_sii-rii_linear_map_latest_v4"),gender,".rds")) 

a <- hist_ineq(dat_ci_f, sample.size, "Sample size per cluster")
b <- hist_ineq(dat_ci_f, prev*100, "Recent HIV testing (%)")
#c <- bi_hist_ineq(dat_ci_f, var_x = prev*100, var_y = sample.size, lab_x = "Recent HIV testing (%)", lab_y = "Sample size per cluster")

(dist <- plot_grid(a,b,#c, 
                   labels = c("A)", "B)"#, "C)"
                                     ), nrow = 1))

ggsave(here::here(paste0(gender,"distribution_suppl.pdf")),
       units = "cm", width = 30, height = 20, dpi = 1200, dist)


######################

gender = "male"

dat_ci_m=readRDS(paste0(here::here("dat_sii-rii_linear_map_latest_v4"),gender,".rds")) 

d <- hist_ineq(dat_ci_m, sample.size, "Sample size per cluster")
e <- hist_ineq(dat_ci_m, prev*100, "Recent HIV testing (%)")
#f <- bi_hist_ineq(dat_ci_m, var_x = prev*100, var_y = sample.size, lab_x = "Recent HIV testing (%)", lab_y = "Sample size per cluster")

(dist <- plot_grid(d,e,#f, 
                   labels = c("A)", "B)", #"C)"
                              ), 
                   nrow = 1))

ggsave(here::here(paste0(gender,"distribution_suppl.pdf")),
       units = "cm", width = 30, height = 20, dpi = 1200, dist)

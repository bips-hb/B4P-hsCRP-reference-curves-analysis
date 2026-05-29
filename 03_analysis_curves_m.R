## ****************************************************************************
##
## Project:       Biomarkers4Pediatrics
## Program name:  03_analysis_curves_m.R
## Author:        Timm Intemann
## R-Version:     4.5.2
## version/ Date: V 1.0/ 2026-05-13
##
## Purpose:       Calculate reference curves for CRP (males)
##
## Program specification: 
## - Prerequisites: none
## - Read:  "./30_Analysis/CRP/data/ana_data_20260424.rds"

## ***************************************************************************

## ***************************************************************************
# 0 Packages, sources, directories ------
## ***************************************************************************

require(sas7bdat)
require(foreign)
require(plyr)
require(dplyr)
require(tidyverse)
require(gamlss)
require(weights)
require(Hmisc)
require(colorspace)
require(lme4)

setwd("E:/Projects/B4P")
source("./30_Analysis/CRP/codes/functions.R")

sessionInfo()

## ****************************************************************************
#1 Load analysis data and data management ----------------------------------
## ****************************************************************************

ana_data_CRP <- readRDS("./30_Analysis/CRP/data/ana_data_20260424.rds")

# select males
ana_data_CRP_m <- subset(ana_data_CRP, sex == 1)

# gamlss analysis data  for two part models
g1_data_CRP_m <- ana_data_CRP_m[, c("age", "crp_le02", "study_acronym", "W")]
g1_data_CRP_m$study_acronym <- as.factor(g1_data_CRP_m$study_acronym)

g2_data_CRP_m <- subset(ana_data_CRP_m,
                        crp_le02 == 0, 
                        select = c("age",
                                   "crp",
                                   "study_acronym",
                                   "W")
                        )

g2_data_CRP_m$shift_crp <- g2_data_CRP_m$crp - 0.02

g2_data_CRP_m$study_acronym <- as.factor(g2_data_CRP_m$study_acronym)

studies <- unique(ana_data_CRP_m$study_acronym)

#set colors
colors <- qualitative_hcl(length(studies), palette = "Dark 3")

## ****************************************************************************
# 2 Analysis: Two-part model  -------------------------------------------------
## ****************************************************************************

#general settings and prediction data
con <- gamlss.control(n.cyc = 500, trace = TRUE)
i.con <- glim.control(bf.cyc = 300)
k <- log(length(g1_data_CRP_m$age))

age_seq <- 10:225/10
newdata <- data.frame(age = age_seq)

## 2.1 First part: mixed log reg ---------------------

mod_logRA_m <- gamlss(crp_le02 ~ pb(age, max.df = 8, method = "GAIC", k=k) + random(study_acronym),
                       data = g1_data_CRP_m, 
                       family = BI(), weights = W, control = con)

### 2.1.1 Plot 1st part of model ---------------
svg("./30_Analysis/CRP/results/Fig_S5A_fitted_log_reg_1st_Part_m.svg")

# Empty plot initializing
plot(NA, xlim = c(1, 24), ylim = c(0, 1),
     xlab = "Age (yrs)", ylab = "Probability (hs-CRP <= 0.2 mg/L)",
     main = "", xaxt = "n", yaxt = "n", type = "n")
title("A   Males", adj = 0)
add_gridlines_simple(x = 0:30, y = 0:10/10)

xg <- 0:12*2 + 1
yg <- seq(0, 1, by = 0.1)
axis(1, at = xg, labels = xg)
axis(2, at = yg, labels = yg, las = 1)

legend("topleft",
       inset = 0.02,
       legend = c(studies, "Overall"),
       col = c(colors, "black"),
       lwd = c(rep(2, length(colors)), 3),
       lty = 1,
       bty = "o", box.lwd = 1.2, bg = "white",
       cex = 0.75, y.intersp = 0.8, ncol = 2
       )

# Plot global curve
newdata <- data.frame(age = age_seq,
                        study_acronym = "")
mu_terms_1_part_global <- predict(mod_logRA_m, newdata = newdata, type = "terms")
pred_prob_global <- plogis(mu_terms_1_part_global[,1] + 
                             attributes(mu_terms_1_part_global)$constant)
lines(age_seq, pred_prob_global, col = "black", lwd = 3)

# Plots studies separately
for (i in seq_along(studies)) {
  newdata_i <- data.frame(age = age_seq,
                          study_acronym = studies[i])
  pred_prob <- predict(mod_logRA_m, newdata = newdata_i, type = "response")
  lines(age_seq, pred_prob, col = colors[i], lwd = 2, lty = 2)
}
dev.off()

## 2.2 Second part ----------------------------

#set new parameters
con <- gamlss.control(n.cyc = 5000, trace = TRUE)
cent <- c(5, 50, 95)
k2 <- log(length(g2_data_CRP_m$age))

###2.2.1 Candidate distributions ------------

start_time <- Sys.time()
mod_GG_m <- gamlss(crp ~ 1,
                   data = g2_data_CRP_m, 
                   family = GG(), weights = W, control = con) 
end_time <- Sys.time()
time_mod_GG_m <- end_time - start_time

start_time <- Sys.time()
mod_GG_m_shift <- gamlss(shift_crp ~ 1,
                         data = g2_data_CRP_m, 
                         family = GG(), weights = W, control = con) 
end_time <- Sys.time()
time_mod_GG_m_shift <- end_time - start_time

start_time <- Sys.time()
mod_GB2_m <- gamlss(crp ~ 1,
                    data = g2_data_CRP_m, 
                    family = GB2(), weights = W, control = con) 
end_time <- Sys.time()
time_mod_GB2_m <- end_time - start_time

start_time <- Sys.time()
mod_GB2_m_shift <- gamlss(shift_crp ~ 1,
                          data = g2_data_CRP_m, 
                          family = GB2(), weights = W, control = con) 
end_time <- Sys.time()
time_mod_GB2_m_shift <- end_time - start_time

start_time <- Sys.time()
mod_GIG_m_shift <- gamlss(shift_crp ~ 1,
                          data = g2_data_CRP_m, 
                          family = GIG(), weights = W, control = con) 
end_time <- Sys.time()
time_mod_GIG_m_shift <- end_time - start_time

start_time <- Sys.time()
mod_GB1_m_shift <- gamlss(shift_crp ~ 1,
                          data = g2_data_CRP_m, 
                          family = GB1(), weights = W, control = con) 
end_time <- Sys.time()
time_mod_GB1_m_shift <- end_time - start_time

####2.2.1.1 Diagnostics (no age) -------------
BIC_CRP_noage_m <- AIC(mod_GG_m,
                       mod_GG_m_shift,
                       mod_GB2_m,
                       mod_GB2_m_shift,
                       mod_GIG_m_shift,
                       mod_GB1_m_shift, 
                       k = k2)

write.csv(BIC_CRP_noage_m, "./30_Analysis/CRP/results/Tab_S8A_BIC_CRP_2nd_Part_noage_m.csv")

#### 2.2.1.2 Histogram and densities (no age) -------------

colors_dist <- colors[c(1, 9)]
svg("./30_Analysis/CRP/results/Fig_S6A_hist_dens_noage_m.svg")
#histogram
wtd.hist(g2_data_CRP_m$shift_crp, 
         weight = g2_data_CRP_m$W, 
         breaks = 100, 
         xlab = "hs-CRP- shifted (mg/L)", 
         main = " ",
         freq = FALSE,
         xaxt = "n",
         ylim = c(0, 16))
 xg_mg_dL <- 0:10/10
 xg_mg_L <- xg_mg_dL*10
 axis(1, at = xg_mg_dL, labels = xg_mg_L)

title("A   Males", adj = 0)

abline(v = wtd.quantile(g2_data_CRP_m$shift_crp, 
                        weight = g2_data_CRP_m$W,
                        probs = 0.5), 
       col = "black",
       lty = 1,
       lwd = 2,
       cex = 3)
abline(v = wtd.quantile(g2_data_CRP_m$shift_crp, 
                        weight = g2_data_CRP_m$W,
                        probs = 0.9), 
       col = "black",
       lty = 1,
       lwd = 3,
       cex = 3)

#densities
distributions <- list(dGG, 
                      dGB1) 
models <- list(mod_GG_m_shift, mod_GB1_m_shift) 

quantile_funcs <- list(qGG, qGB1)

x_vals <- 1:100/100

for (i in seq_along(distributions)) {
  dist_fun <- distributions[[i]]
  model <- models[[i]]
  
  # Arg for density and quantile fct
  args <- list()
  if (!is.null(model$mu.fv))    args$mu    <- model$mu.fv[1]
  if (!is.null(model$sigma.fv)) args$sigma <- model$sigma.fv[1]
  if (!is.null(model$nu.fv))    args$nu    <- model$nu.fv[1]
  if (!is.null(model$tau.fv))   args$tau   <- model$tau.fv[1]
  
  # Density plot
  y_vals <- do.call(dist_fun, c(list(x = x_vals), args))
  lines(x_vals, y_vals, col = colors_dist[i], lwd = 3)
  
  # Calculate model quantiles
  qfun <- quantile_funcs[[i]]
  q50 <- do.call(qfun, c(list(p = 0.5), args))
  q90 <- do.call(qfun, c(list(p = 0.9), args))
  
  # Abline in each color
  abline(v = q50, col = colors_dist[i], lty = 2, lwd = 3)
  abline(v = q90, col = colors_dist[i], lty = 4, lwd = 3)
}

legend_labels <- c("GG - shifted", "GB1 - shifted") 
legend("topright", legend = legend_labels, col = colors_dist, lty = 2, lwd = 3, cex = 1)
dev.off()

### 2.2.2 Age-specific modelling -----------------------
# using the best distribution from before: GG, GB1 (both shifted)

#complete complex models
start_time <- Sys.time()
mod_age_GB1_m_all <- gamlss(shift_crp ~ pb(age, max.df = 8, method = "GAIC", k=k2) + random(study_acronym),
                            sigma.fo = ~ pb(age, max.df = 4, method = "GAIC", k=k2),
                            nu.fo = ~ pb(age, max.df = 4, method = "GAIC", k=k2),
                            tau.fo = ~ pb(age, max.df = 4, method = "GAIC", k=k2),
                            data = g2_data_CRP_m, 
                            family = GB1(), weights = W, control = con) 
end_time <- Sys.time()
time_mod_age_GB1_m_all <- end_time - start_time

start_time <- Sys.time()
mod_age_GB1_m_all2 <- gamlss(shift_crp ~ pb(age, max.df = 8, method = "GAIC", k=k2) + random(study_acronym),
                            sigma.fo = ~ pb(age, max.df = 6, method = "GAIC", k=k2),
                            nu.fo = ~ pb(age, max.df = 4, method = "GAIC", k=k2),
                            tau.fo = ~ pb(age, max.df = 4, method = "GAIC", k=k2),
                            data = g2_data_CRP_m, 
                            family = GB1(), weights = W, control = con) 
end_time <- Sys.time()
time_mod_age_GB1_m_all2 <- end_time - start_time

start_time <- Sys.time()
mod_age_GG_m_pre_all <- gamlss(shift_crp ~ pb(age, max.df = 8, method = "GAIC", k=k2) + random(study_acronym),
                           sigma.fo = ~ pb(age, max.df = 4, method = "GAIC", k=k2),
                           nu.fo = ~ 1,
                           data = g2_data_CRP_m, 
                           family = GG(), weights = W, control = con) 
end_time <- Sys.time()
time_mod_age_GG_m_pre_all <- end_time - start_time

#### 2.2.2.1 Diagnostics ------------------------
svg("./30_Analysis/CRP/results/Fig_S7A1_WP_final_2nd_part_m.svg")
wgt_wp(mod_age_GG_m_pre_all$residuals,
       mod_age_GG_m_pre_all$weights, 
       mod_name_title = paste("(shifted GG with", round(mod_age_GG_m_pre_all$df.fit), "degrees of freedom)"))
dev.off()
svg("./30_Analysis/CRP/results/Fig_S7A2_WP_final_2nd_part_m.svg")
wgt_wp(mod_age_GB1_m_all$residuals,
       mod_age_GB1_m_all$weights, 
       mod_name_title = paste("(shifted GB1 with", round(mod_age_GB1_m_all$df.fit), "degrees of freedom)"))
dev.off()
svg("./30_Analysis/CRP/results/Fig_S7A3_WP_final_2nd_part_m.svg")
wgt_wp(mod_age_GB1_m_all2$residuals,
       mod_age_GB1_m_all2$weights, 
       mod_name_title = paste("(shifted GB1 with", round(mod_age_GB1_m_all2$df.fit), "degrees of freedom)"))
dev.off()

BIC_summary_CRP_best_m <- AIC(mod_age_GB1_m_all, 
                              mod_age_GB1_m_all2, 
                              mod_age_GG_m_pre_all, k = k2)

write.csv(BIC_summary_CRP_best_m, "./30_Analysis/CRP/results/Tab_S9A_BIC_part2_m.csv")

#### 2.2.2.2 Curve plotting -----------------

newdata <- data.frame(age = 10:225/10)
newdata$study_acronym <- "ABCD" #as placeholder
cent <- c(5, 50, 95)

svg("./30_Analysis/CRP/results/Fig_S8A_2nd_mod-comp_m.svg")
plot(c(-100, 10000), rep(20, times = 2), 
     ylim = c(0, 1), xlim = range(newdata$age) * c(1, 1.06), 
     col = grey(0.9),
     xlab = "Age (yrs)", ylab = "hs-CRP (shifted) (mg/L)",
     main = " ", , xaxt = "n", yaxt = "n", type = "n")
title("A  Males", adj = 0)

add_gridlines(x = 0:30, y = 0:10/10,
              P = 23.5, cent = cent, y_P = c(0, 0.05, 0.5))

xg <- 0:12*2 + 1
yg_mg_dl <- seq(0, 1, by = 0.1)
yg_mg_l <- yg_mg_dl*10
axis(1, at = xg, labels = xg)
axis(2, at = yg_mg_dl, labels = yg_mg_l, las = 1)

add_plot_global_curves_GG(mod = mod_age_GG_m_pre_all, newdata, cent, 
                          mu_term_name = 'pb(age, max.df = 8, method = "GAIC", k = k2)',
                          col_curve = "darkgreen")

add_plot_global_curves_GB1(mod = mod_age_GB1_m_all, newdata, cent, 
                       mu_term_name = 'pb(age, max.df = 8, method = "GAIC", k = k2)',
                       col_curve = "darkblue")

add_plot_global_curves_GB1(mod = mod_age_GB1_m_all2, newdata, cent, 
                       mu_term_name = 'pb(age, max.df = 8, method = "GAIC", k = k2)',
                       col_curve = "black")
legend("topleft", 
       inset = 0.02,
       legend = c(paste("Shifted GG with", round(mod_age_GG_m_pre_all$df.fit), "degrees of freedom"), 
                  paste("Shifted GB1 with", round(mod_age_GB1_m_all$df.fit), "degrees of freedom"), 
                  paste("Shifted GB1 with", round(mod_age_GB1_m_all2$df.fit), "degrees of freedom")
                  ), 
       col = c("darkgreen", "darkblue", "black"), 
       lty = 1, 
       lwd = 2, title = "Models",
       box.lwd = 1.2, 
       bg = "white",
       cex = 0.75, 
       y.intersp = 0.8, 
       ncol = 1
)
dev.off()# -> final model: mod_age_GB1_m_all2

## save final model 
save(mod_age_GB1_m_all2, file = "./30_Analysis/CRP/data/mod_age_GB1_m_all2.RData")
#load("./30_Analysis/CRP/data/mod_age_GB1_m_all2.RData")

#Plot final 2nd part of model including study specific curves
newdata <- data.frame(age = 10:225/10)
newdata$study_acronym <- " "

svg("./30_Analysis/CRP/results/Fig_S9A_2nd_Part_final_m.svg")
plot(c(-100, 10000), rep(20, times = 2), 
     ylim = c(0, 1), xlim = range(newdata$age) * c(1, 1.06), 
     col = grey(0.9),
     xlab = "Age (yrs)", ylab = "hs-CRP (shifted) (mg/L)",
     main = "", xaxt = "n", yaxt = "n", type = "n")

title("A  Males", adj = 0)

add_gridlines(x = 0:30, y = 0:10/10,
              P = 23.5, cent = cent, y_P = c(0, 0.05, 0.5))

xg <- 0:12*2 + 1
yg_mg_dl <- seq(0, 1, by = 0.1)
yg_mg_l <- yg_mg_dl*10
axis(1, at = xg, labels = xg)
axis(2, at = yg_mg_dl, labels = yg_mg_l, las = 1)

add_plot_all_study_curves_GB1(mod = mod_age_GB1_m_all2, 
                              newdata = newdata, 
                              cent = cent,
                              list_studies = studies,
                              lwd = 2, 
                              color_vector = colors)

add_plot_global_curves_GB1(mod = mod_age_GB1_m_all2, 
                           newdata = newdata, 
                           cent = cent, 
                           mu_term_name = 'pb(age, max.df = 8, method = "GAIC", k = k2)',
                           col_curve = "black")

legend("topleft",
       inset = 0.02,
       legend = c(studies, "Overall"),
       col = c(colors, "black"),
       lwd = c(rep(2, length(colors)), 3),
       lty = 1,
       bty = "o", box.lwd = 1.2, bg = "white",
       cex = 0.75, y.intersp = 0.8, ncol = 2
       )

dev.off()

## ****************************************************************************
# 3 Combine 1st and 2nd part ------------------------
## ****************************************************************************

## 3.1 Compile data for plotting ---------------------------

All_perc_comb_m <- calc_all_perc_hurdle(age_seq = age_seq, 
                                        mod1 = mod_logRA_m, 
                                        mod2 = mod_age_GB1_m_all2,
                                        studies = studies)
#Select columns
Tab_S7A <- All_perc_comb_m$global[,c("age",
                                     "P_1",
                                     "Perc_comb_50",
                                     "Perc_comb_75",
                                     "Perc_comb_85",
                                     "Perc_comb_90",
                                     "Perc_comb_95",
                                     "Perc_comb_97",
                                     "Perc_comb_cutoff")]

#Unit transformation
Tab_S7A <-  Tab_S7A * rep(c(1, 100, rep(10, times = 7)), each = nrow(Tab_S7A))
#Rounding
Tab_S7A <- round(Tab_S7A, 2)
#Name columns
names(Tab_S7A) <- c("Age (yrs)",
                    "Prob. (hs-CRP <= 0.2 mg/L)",
                    "P50", "P75", "P85", "P90", "P95", "P97", "Cut-off")
#Write csv
write.csv(Tab_S7A, "./30_Analysis/CRP/results/Tab_S7A_All_perc_comb_m.csv", row.names = FALSE)

##3.2 Plotting  ----------------------

###3.2.1 Cutoffs ------------------

pdf("./30_Analysis/CRP/results/Fig_3A_CRP_curves_cutoff_m.pdf")
par(mar = c(5, 4, 4, 5) + 0.1)
plot(c(-100, 10000), rep(20, times = 2), 
     ylim = c(0, 1), xlim = c(0.5, 23.5) , 
     col = grey(0.9),
     xlab = "Age (yrs)", ylab = "hs-CRP (mg/L)", 
     main = "", xaxt = "n", yaxt = "n", type = "n")

title("A  Males", adj = 0)

xg <- 0:12*2 + 1
yg_mg_dl <- seq(0, 1, by = 0.1)
yg_mg_l <- yg_mg_dl*10
axis(1, at = xg, labels = xg)
axis(2, at = yg_mg_dl, labels = yg_mg_l, las = 1)

# add right axis in nmol/L
conv <- 95.2381   # mg/L -> nmol/L
yticks <- pretty(seq(0, 1, length.out = 5))   # ticks in mg/L scale
axis(side = 4,
     at = yticks,
     labels = round(yticks * conv, 1), las = 1)        # convert to nmol/L
mtext("hs-CRP (nmol/L)", side = 4, line = 3)

add_gridlines(x = 0:30, y = 0:10/10,
              P = 23, cent = round(All_perc_comb_m$global$Perc_comb_18[1]*100, 1), 
              y_P = All_perc_comb_m$global$Perc_comb_cutoff[All_perc_comb_m$global$age==22.5] - 0.04)

add_gridlines_simple(x = c(18), y = 0.3, col = "darkgrey", lwd = 3)

lines(All_perc_comb_m$global$age, 
      All_perc_comb_m$global$Perc_comb_cutoff, 
      lwd = 3, lty = 1, col = "black")

for (i in 1:(length(studies))){
  lines(All_perc_comb_m[[i]]$age, 
        All_perc_comb_m[[i]]$Perc_comb_cutoff, 
        lwd = 2, lty = 2, col = colors[i])
}

lines(All_perc_comb_m$global$age, 
      All_perc_comb_m$global$Perc_comb_cutoff, 
      lwd = 3, lty = 1, col = "black")

legend("topleft",
       inset = 0.02,
       legend = c(studies, "Overall"),
       col = c(colors, "black"),
       lwd = c(rep(2, length(colors)), 3),
       lty = 1,
       bty = "o", box.lwd = 1.2, bg = "white",
       cex = 0.75, y.intersp = 0.8, ncol = 2
) 
dev.off()

###3.2.2 All global curves ----------------------------
pdf("./30_Analysis/CRP/results/Fig_2A_curves_global_total_m.pdf")
par(mar = c(5, 4, 4, 5) + 0.1)
plot(c(-100, 10000), rep(20, times = 2), 
     ylim = c(0, 1), xlim = c(0.5, 23.5) , 
     col = grey(0.9),
     xlab = "Age (yrs)", ylab = "hs-CRP (mg/L)", 
     main = "", xaxt = "n", yaxt = "n", type = "n")

title("A  Males", adj = 0)

add_gridlines(x = 0:30, y = 0:10/10,
              P = 23.5, cent = c(50, 75, 85, 90, 95, 97), 
              y_P = c(All_perc_comb_m$global$Perc_comb_50[All_perc_comb_m$global$age==22.5],
                      All_perc_comb_m$global$Perc_comb_75[All_perc_comb_m$global$age==22.5],
                      All_perc_comb_m$global$Perc_comb_85[All_perc_comb_m$global$age==22.5],
                      All_perc_comb_m$global$Perc_comb_90[All_perc_comb_m$global$age==22.5],
                      All_perc_comb_m$global$Perc_comb_95[All_perc_comb_m$global$age==22.5],
                      All_perc_comb_m$global$Perc_comb_97[All_perc_comb_m$global$age==22.5])
              )

xg <- 0:12*2 + 1
yg_mg_dl <- seq(0, 1, by = 0.1)
yg_mg_l <- yg_mg_dl*10
axis(1, at = xg, labels = xg)
axis(2, at = yg_mg_dl, labels = yg_mg_l, las = 1)

# add right axis in nmol/l 
conv <- 95.2381   # mg/L -> nmol/L
yticks <- pretty(seq(0, 1, length.out = 5))   # ticks in mg/L scale
axis(side = 4,
     at = yticks,
     labels = round(yticks * conv, 1), las = 1)        # convert to nmol/L
mtext("hs-CRP (nmol/L)", side = 4, line = 3)

lines(All_perc_comb_m$global$age, 
      All_perc_comb_m$global$Perc_comb_50, 
      lwd = 3, lty = 1, col = "black")
lines(All_perc_comb_m$global$age, 
      All_perc_comb_m$global$Perc_comb_75, 
      lwd = 3, lty = 1, col = "black")
lines(All_perc_comb_m$global$age, 
      All_perc_comb_m$global$Perc_comb_85, 
      lwd = 3, lty = 1, col = "black")
lines(All_perc_comb_m$global$age,
      All_perc_comb_m$global$Perc_comb_90,
      lwd = 3, lty = 1, col = "black")
lines(All_perc_comb_m$global$age, 
      All_perc_comb_m$global$Perc_comb_95, 
      lwd = 3, lty = 1, col = "black")
lines(All_perc_comb_m$global$age,
      All_perc_comb_m$global$Perc_comb_97, 
      lwd = 3, lty = 1, col = "black")

dev.off()

###3.2.3 P50 and P90 as study-specific curves ---------------------

# all global curves
svg("./30_Analysis/CRP/results/Fig_S10A_CRP_curves_studyspec_total_m.svg")
par(mar = c(5, 4, 4, 5) + 0.1)
plot(c(-100, 10000), rep(20, times = 2), 
     ylim = c(0, 1), xlim = c(0.5, 23.5) , 
     col = grey(0.9),
     xlab = "Age (yrs)", ylab = "hs-CRP (mg/L)",
     main = "", xaxt = "n", yaxt = "n", type = "n")

title("A   Males", adj = 0)
add_gridlines(x = 0:30, y = 0:10/10,
              P = 26, cent = cent, y_P = c(0.0, 0.1, 0.65))
xg <- 0:12*2 + 1
yg_mg_dl <- seq(0, 1, by = 0.1)
yg_mg_l <- yg_mg_dl*10
axis(1, at = xg, labels = xg)
axis(2, at = yg_mg_dl, labels = yg_mg_l, las = 1)

# add right axis in nmol/l
conv <- 95.2381   # mg/L -> nmol/L
yticks <- pretty(seq(0, 1, length.out = 5))   # ticks in mg/L scale
axis(side = 4,
     at = yticks,
     labels = round(yticks * conv, 1), las = 1)        # convert to nmol/L
mtext("hs-CRP (nmol/L)", side = 4, line = 3)

#global part
lines(All_perc_comb_m$global$age, 
      All_perc_comb_m$global$Perc_comb_50, 
      lwd = 3, lty = 1, col = "black")
lines(All_perc_comb_m$global$age,
      All_perc_comb_m$global$Perc_comb_90,
      lwd = 3, lty = 1, col = "black")

text(23.5, All_perc_comb_m$global$Perc_comb_50[All_perc_comb_m$global$age==22.5], "P50")
text(23.5, All_perc_comb_m$global$Perc_comb_90[All_perc_comb_m$global$age==22.5], "P90")

for (i in seq_along(studies)) {
  nm <- studies[i]
  df <- All_perc_comb_m[[i]]
  sel <- df$age < 22.5
  
  lines(df$age[sel], df$Perc_comb_50[sel], lwd = 2, lty = 2, col = colors[i])
  lines(df$age[sel], df$Perc_comb_90[sel], lwd = 2, lty = 2, col = colors[i])
}

legend("topleft",
       inset = 0.02,
       legend = c(studies, "Overall"),
       col = c(colors, "black"),
       lwd = c(rep(2, length(colors)), 3),
       lty = 1,
       bty = "o", box.lwd = 1.2, bg = "white",
       cex = 0.75, y.intersp = 0.8, ncol = 2
)

dev.off()

## 3.3 Compilation for data tables ---------------------------

All_perc_comb_long_m <- calc_all_perc_hurdle(age_seq = 100:2250/100, 
                                             mod1 = mod_logRA_m, 
                                             mod2 = mod_age_GB1_m_all2,
                                             studies = studies)

write.csv(All_perc_comb_long_m, "./30_Analysis/CRP/data/All_perc_comb_long_m.csv", row.names = FALSE)

#4 Association of repeated CRP values in the reference system ------------------

## 4.1 Select repeated measurements data -------
rep_ana_data <- ana_data_CRP_m[ana_data_CRP_m$age<=22.5, 
                               c("b4p_id", "age", "crp", "study_acronym", "crp_le02")] %>%
  group_by(b4p_id) %>%
  filter(n() > 1) %>%
  ungroup()

## 4.2 Calculate individual percentile ranks -------
rep_ana_data <- calc_individual_P(df = rep_ana_data, 
                                  mod1 = mod_logRA_m, 
                                  mod2 = mod_age_GB1_m_all2)

## 4.3 Calculate ICCs -------------

#fit linear mixed models
mod_P <- lmer(P_crp ~ 1 + (1 | b4p_id), data = rep_ana_data, REML = TRUE)
mod_raw <- lmer(crp ~ 1 + (1 | b4p_id), data = rep_ana_data, REML = TRUE)

#extract variances
vc_P <- as.data.frame(VarCorr(mod_P))
vc_raw <- as.data.frame(VarCorr(mod_raw))

sigma_person_P <- vc_P$vcov[vc_P$grp == "b4p_id"]
sigma_resid_P  <- vc_P$vcov[vc_P$grp == "Residual"]

sigma_person_raw <- vc_raw$vcov[vc_raw$grp == "b4p_id"]
sigma_resid_raw  <- vc_raw$vcov[vc_raw$grp == "Residual"]

ICC_P <- sigma_person_P / (sigma_person_P + sigma_resid_P) 
ICC_raw <- sigma_person_raw / (sigma_person_raw + sigma_resid_raw) 

ICC_P
ICC_raw
dim(rep_ana_data)

## 4.4. ORs for having an elevated CRP values (ie >=0.3) when 18 -----------

#data preparation: wide format with indicators for cutoffs
rep_ana_data_over18 <- subset(rep_ana_data, age>=18, select = c("b4p_id", "crp"))
rep_ana_data_under18 <- subset(rep_ana_data, age<18, select = c("b4p_id", "crp", "P_crp"))

rep_ana_data_over18$crp_0.3_over18 <- rep_ana_data_over18$crp>=0.3
rep_ana_data_under18$crp_high_under18 <- rep_ana_data_under18$P_crp>=0.935

rep_ana_data_under18_max <- rep_ana_data_under18 %>%
  group_by(b4p_id) %>%
  summarise(max_crp = max(crp_high_under18, na.rm = TRUE)
            )

ana_data_over_under <- merge(rep_ana_data_under18_max, 
                             rep_ana_data_over18, 
                             by = "b4p_id",
                             all = FALSE)

# fit log model
mod_log <- glm(crp_0.3_over18 ~ max_crp,
               data = ana_data_over_under, family = "binomial")     

#extract results
exp(coefficients(mod_log)[2]) #OR
exp(confint.default(mod_log)[2, ]) #corresponding 95%-CI

length(mod_log$residuals) #sample size
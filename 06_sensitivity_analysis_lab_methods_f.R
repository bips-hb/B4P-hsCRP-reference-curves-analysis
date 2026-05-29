## ****************************************************************************
##
## Project:       Biomarkers4Pediatrics
## Program name:  06_sensitivity_analysis_lab_methods_f.R
## Author:        Jiayi Zeng
## R-Version:     4.5.2
## version/ Date: V 1.0/ 2026-05-25
##
## Purpose:       Sensitivity analysis for different lab methods (females)
##
## Program specification: 
## - Prerequisites: none
## - Read:  *Data: "./data/ana_data_20260424.rds"
##          *Cutoff file: "./30_Analysis/CRP/data/All_perc_comb_long_f.csv"
##          *Model: "./30_Analysis/CRP/data/mod_age_GB1_f_all2.RData"
## ***************************************************************************

## ***************************************************************************
# 1 Packages, sources, directories ------
## ***************************************************************************

require(sas7bdat)
require(foreign)
require(plyr)
require(dplyr)
require(tidyverse)
require(gamlss)
require(weights)
require(Hmisc)
library(scico)
library(colorspace)

setwd("E:/Projects/B4P") 
source("./30_Analysis/CRP/codes/functions.R")
sessionInfo()

## ****************************************************************************
#2 Load analysis data and data management ----------------------------------
## ****************************************************************************
ana_data_CRP <- readRDS("./30_Analysis/CRP/data/ana_data_20260424.rds") #data
load("./30_Analysis/CRP/data/mod_age_GB1_f_all2.RData")   #model, object: mod_age_GB1_f_all2 

# select females
ana_data_CRP_f <- subset(ana_data_CRP, sex == 2)

# gamlss analysis data for two part models
g1_data_CRP_f <- ana_data_CRP_f[, c("age", "methods","crp_le02", "study_acronym", "W")]
g1_data_CRP_f$study_acronym <- as.factor(g1_data_CRP_f$study_acronym)

g2_data_CRP_f <- subset(ana_data_CRP_f,
                        crp_le02 == 0, 
                        select = c("age",
                                   "crp",
                                   "study_acronym",
                                   "methods",
                                   "W")
)

g1_data_CRP_f$methods <- droplevels(factor(g1_data_CRP_f$methods))
g2_data_CRP_f$methods <- droplevels(factor(g2_data_CRP_f$methods))

g2_data_CRP_f$shift_crp <- g2_data_CRP_f$crp - 0.02
g2_data_CRP_f$study_acronym <- as.factor(g2_data_CRP_f$study_acronym)

# methods reference setting
ref_meth <- names(sort(table(g1_data_CRP_f$methods), decreasing = TRUE))[1]
g1_data_CRP_f$methods <- relevel(g1_data_CRP_f$methods, ref = ref_meth)
g2_data_CRP_f$methods <- relevel(g2_data_CRP_f$methods, ref = ref_meth)

## ****************************************************************************
# 3 Analysis: Two-part model  -------------------------------------------------
## ****************************************************************************

## 3.1 First part: log reg ---------------------
#general settings
con <- gamlss.control(n.cyc = 500, trace = TRUE)
i.con <- glim.control(bf.cyc = 300)
age_seq <- 10:250/10
newdata <- data.frame(age = age_seq)

# define basis parameter
k <- log(length(g1_data_CRP_f$age))

### 3.1.1 Without random effect -------------- 
mod_log_f <- gamlss(crp_le02 ~ pb(age, max.df = 8, method = "GAIC", k=k),
                    data = g1_data_CRP_f, 
                    family = BI(), weights = W, control = con)

### 3.1.2 Adding study-specific random effect ------------------
mod_logRA_f <- gamlss(crp_le02 ~ pb(age, max.df = 8, method = "GAIC", k=k) + random(study_acronym),
                      data = g1_data_CRP_f,
                      family = BI(), weights = W, control = con)

### 3.1.3 Adding study random effect and methods------------

mod_logRA_meth_f <- gamlss(
  crp_le02 ~ pb(age, max.df = 8, method = "GAIC", k=k) + methods + random(study_acronym),
  data = g1_data_CRP_f, family = BI(), weights = W, control = con
)


## 3.2 Second part ----------------------------

# set parameters
con <- gamlss.control(n.cyc = 5000, trace = TRUE)
i.con <- glim.control(bf.cyc = 300)
cent <- c(5, 50, 95)
k2 <- log(length(g2_data_CRP_f$age))

### 3.2.1 Sensitivity model (+methods in mu) ------------

#fit the model 
mod_age_GB1_f_all2_m <- gamlss(
   shift_crp ~ pb(age, max.df = 8, method = "GAIC", k=k2) + methods + random(study_acronym),
   sigma.fo = ~ pb(age, max.df = 6, method = "GAIC", k=k2),
   nu.fo    = ~ pb(age, max.df = 4, method = "GAIC", k=k2),
   tau.fo   = ~ pb(age, max.df = 4, method = "GAIC", k=k2),
   data = g2_data_CRP_f, family = GB1(), weights = W, control = con
 )

#save the model (mod_age_GB1_f_all2_m)
save(mod_age_GB1_f_all2_m, file = "./30_Analysis/CRP/data/methods_added_mod_age_GB1_f_all2_meth.RData")
# load("./30_Analysis/CRP/data/methods_added_mod_age_GB1_f_all2_meth.RData")   

## ****************************************************************************
# 4 Combine 1st and 2nd Part ------------------------
## ****************************************************************************

#All_perc_comb_f.rds saved
df <- read.csv("./30_Analysis/CRP/data/All_perc_comb_long_f.csv")
df_global <- df %>%
  select(starts_with("global."))
names(df_global) <- sub("^global\\.", "", names(df_global))
All_perc_comb_f <- list()
All_perc_comb_f[["global"]] <- df_global
write_rds(All_perc_comb_f, "./30_Analysis/CRP/data/All_perc_comb_long_f.rds")

## 4.1 Methods added global data compiling ---------------
 
levels_methods <- levels(model.frame(mod_logRA_meth_f)$methods)

levels_study   <- levels(g1_data_CRP_f$study_acronym)

All_perc_comb_f_meth <- list()

for (m in levels_methods) {
  newdata <- data.frame(
    age = age_seq,
    study_acronym = factor(levels_study[1], levels = levels_study),
    methods = factor(m, levels = levels_methods)
  )

  ## Logistic part - use only fixed effects (age, method, intercept),
  ## exclude study random effects, then transform with plogis() for global curve

  tt <- predict(mod_logRA_meth_f, newdata = newdata, type = "terms")
  keep <- !grepl("(^re\\(|^random\\()|study_acronym", colnames(tt))
  eta1 <- rowSums(tt[, keep, drop = FALSE]) + attr(tt, "constant")
  P_1  <- plogis(eta1)

  P2_50 <- Calc_Perc_2(P_1, 0.5)
  P2_75 <- Calc_Perc_2(P_1, 0.75)
  P2_85 <- Calc_Perc_2(P_1, 0.85)
  P2_90 <- Calc_Perc_2(P_1, 0.9)
  P2_95 <- Calc_Perc_2(P_1, 0.95)
  P2_97 <- Calc_Perc_2(P_1, 0.97)

  ## mu of GB1: use only fixed effects (age, method, intercept),
  ## exclude study random effects, then transform with plogis() for global curve

  tt_mu <- predict(mod_age_GB1_f_all2_m, newdata = newdata, type = "terms", what = "mu")
  keep_mu <- !grepl("(^re\\(|^random\\()|study_acronym", colnames(tt_mu))
  mu_eta <- rowSums(tt_mu[, keep_mu, drop = FALSE]) + attr(tt_mu, "constant")
  mu <- plogis(mu_eta)

  ## sigma/nu/tau:no methods or random(study) ---
  sigma <- predict(mod_age_GB1_f_all2_m, newdata = newdata, type = "response", what = "sigma")
  nu    <- predict(mod_age_GB1_f_all2_m, newdata = newdata, type = "response", what = "nu")
  tau   <- predict(mod_age_GB1_f_all2_m, newdata = newdata, type = "response", what = "tau")

  ## combine and shift
  Perc_comb_50 <- qGB1(P2_50, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
  Perc_comb_75 <- qGB1(P2_75, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
  Perc_comb_85 <- qGB1(P2_85, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
  Perc_comb_90 <- qGB1(P2_90, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
  Perc_comb_95 <- qGB1(P2_95, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
  Perc_comb_97 <- qGB1(P2_97, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02

  All_perc_comb_f_meth[[m]] <- data.frame(
    age = age_seq,
    P_1 = P_1,
    P2_50 = P2_50, P2_75 = P2_75, P2_85 = P2_85, P2_90 = P2_90, P2_95 = P2_95, P2_97 = P2_97,
    mu = mu, sigma = sigma, nu = nu, tau = tau,
    Perc_comb_50 = Perc_comb_50,
    Perc_comb_75 = Perc_comb_75,
    Perc_comb_85 = Perc_comb_85,
    Perc_comb_90 = Perc_comb_90,
    Perc_comb_95 = Perc_comb_95,
    Perc_comb_97 = Perc_comb_97
  )
}

All_perc_comb_f_meth <- add_cutoffs(cutoff = 0.3, shift = 0.02, all_perc_list = All_perc_comb_f_meth)

write.csv(All_perc_comb_f_meth, "./30_Analysis/CRP/data/sensitivity_methods_added_All_perc_comb_f.csv")
write_rds(All_perc_comb_f_meth, "./30_Analysis/CRP/data/sensitivity_methods_added_All_perc_comb_f.rds")

## 4.2 Generate 3-methods and global curves-----------------

### 4.2.1 Data load----------
# All_perc_comb_f <- readRDS("./30_Analysis/CRP/data/All_perc_comb_long_f.rds") 
# All_perc_comb_f_meth <- readRDS("./30_Analysis/CRP/data/sensitivity_methods_added_All_perc_comb_f.rds")

orig_global <- All_perc_comb_f[["global"]]
methods_vec <- names(All_perc_comb_f_meth)

### 4.2.2 Detail setting--------

# colors
meth_cols <- qualitative_hcl(13, "Dark 3")
base_cols <- meth_cols[c(1, 5, 9)]
meth_cols <- setNames(base_cols[seq_along(methods_vec)], methods_vec)

# y axis
y_all <- c(orig_global$Perc_comb_50, orig_global$Perc_comb_90)
for (m in methods_vec) {
  dff <- All_perc_comb_f_meth[[m]]
  y_all <- c(y_all, dff$Perc_comb_50, dff$Perc_comb_90)
}
ylim_use <- c(0, 1)

# age until 22.5
sel_g <- orig_global$age < 22.5

### 4.2.3 SVG generation----------
svg("./30_Analysis/CRP/results/Fig_S11B_methods_f.svg",
    width = 7, height = 6)

par(mar = c(5, 5, 3, 5) + 0.1)

# 
plot(
  NA,
  xlim = c(0.5, 23.5),
  ylim = ylim_use,           # 0-1 mg/dL 
  xlab = "Age (yrs)",
  ylab = "",
  xaxt = "n",
  yaxt = "n",
  las  = 1,
  type = "n"
)
title("B  Females", adj = 0)
add_gridlines(x = 0:30, y = 0:10/10,
              P = 26, cent = cent, y_P = c(0.0, 0.1, 0.65))

## x: 1,3,5,...,23
axis(1, at = seq(1, 23, by = 2), labels = seq(1, 23, by = 2))

## 
axis(2,
     at     = seq(0, 1, 0.1),    # 0-1(mg/dL)
     labels = seq(0, 10, 1),     # 0-10(mg/L)
     las = 1)
mtext("hs-CRP (mg/L)", side = 2, line = 3.2)

## 
conv <- 95.2381
yticks <- seq(0, 1, 0.1)

raw_ticks <- c(0, 0.25, 0.50, 0.75, 1.00)
converted_ticks <- round(raw_ticks * conv, 1)
# 
axis(4,
     at     = raw_ticks,
     labels = round(converted_ticks, 1),
     las    = 1)
mtext("hs-CRP (nmol/L)", side = 4, line = 3)

## grid
abline(v = 0:30,                lty = 1, lwd = 0.5, col = "grey85")
abline(h = seq(0, 1, by = 0.1), lty = 1, lwd = 0.5, col = "grey85")

# global
lines(orig_global$age[sel_g], orig_global$Perc_comb_50[sel_g], lwd = 3, lty = 1, col = "black")
lines(orig_global$age[sel_g], orig_global$Perc_comb_90[sel_g], lwd = 3, lty = 1, col = "black")

# method specific percentile curves
for (m in methods_vec) {
  dff <- All_perc_comb_f_meth[[m]]
  sel <- dff$age < 22.5
  lines(dff$age[sel], dff$Perc_comb_50[sel], lwd = 2, lty = 2, col = meth_cols[m])
  lines(dff$age[sel], dff$Perc_comb_90[sel], lwd = 2, lty = 2, col = meth_cols[m])
}

## calculate sample size and number of waves per method 
df_counts_f <- ana_data_CRP_f %>%
  dplyr::filter(methods %in% methods_vec) %>%
  dplyr::mutate(
    wave_norm = ifelse(is.na(wave), "__MISSING__", as.character(wave)),
    key_sw    = paste(study_acronym, wave_norm, sep = "|")
  ) %>%
  dplyr::group_by(methods) %>%
  dplyr::summarise(
    n     = dplyr::n(),
    waves = dplyr::n_distinct(key_sw),
    .groups = "drop"
  )

# join method vector with summary statistics
summ_f <- data.frame(methods = as.character(methods_vec)) %>%
  left_join(df_counts_f, by = "methods")

#prepare legends including summary statistics
legend_methods <- sprintf("%s (n = %s, waves = %d)",
                          summ_f$methods, format(summ_f$n, big.mark=","), summ_f$waves)

##  label in the figure
x_last <- max(orig_global$age[sel_g])
y50    <- tail(orig_global$Perc_comb_50[sel_g], 1)
y90    <- tail(orig_global$Perc_comb_90[sel_g], 1)
text(x_last, y50, "P50", pos = 4, cex = 0.8, xpd = TRUE, offset = 0.3)
text(x_last, y90, "P90", pos = 4, cex = 0.8, xpd = TRUE, offset = 0.3)

##  legend
legend("topleft",
       legend = legend_methods,
       col    = unname(meth_cols[methods_vec]),
       lty    = 2, lwd = 2,
       bg     = "white", bty = "o", box.lwd = 1.2,
       cex    = 0.8, inset = 0.02)
dev.off()

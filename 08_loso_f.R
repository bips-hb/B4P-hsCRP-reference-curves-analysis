## ****************************************************************************
##
## Project:       Biomarkers4Pediatrics
## Program name:  08_loso_f.R
## Author:        Jiayi Zeng
## R-Version:     4.5.2
## version/ Date: V 1.0/ 2026-05-25
##
## Purpose:       leave one study out validation (females)
##
## Program specification: 
## - Prerequisites: none
## - Read:  *Data: "./data/ana_data_20260424.rds"
##          *Model: "./30_Analysis/CRP/data/mod_age_GB1_f_all2.RData"
##
## ***************************************************************************

## ***************************************************************************
# 0 Packages, sources, directories ------
## ***************************************************************************

library(colorspace)
require(sas7bdat)
require(foreign)
require(plyr)
require(dplyr)
require(tidyverse)
require(gamlss)
require(weights)
require(Hmisc)

setwd("E:/Projects/B4P") 
source("./30_Analysis/CRP/codes/functions.R")
sessionInfo()

## ***************************************************************************
# 1 Data load-------------
## ***************************************************************************
ana_data_CRP <- readRDS("./30_Analysis/CRP/data/ana_data_20260424.rds")
load("./30_Analysis/CRP/data/mod_age_GB1_f_all2.RData")   # model, object: mod_age_GB1_f_all2

# select females
ana_data_CRP_f <- subset(ana_data_CRP, sex == 2)

# gamlss analysis data  for two part models
g1_data_CRP_f <- ana_data_CRP_f[, c("age", "crp_le02", "study_acronym", "W")]
g1_data_CRP_f$study_acronym <- as.factor(g1_data_CRP_f$study_acronym)

g2_data_CRP_f <- subset(ana_data_CRP_f,
                        crp_le02 == 0, 
                        select = c("age",
                                   "crp",
                                   "study_acronym",
                                   "W")
)

g2_data_CRP_f$shift_crp <- g2_data_CRP_f$crp - 0.02

g2_data_CRP_f$study_acronym <- as.factor(g2_data_CRP_f$study_acronym)

## ****************************************************************************
# 2 Preliminary analysis -------------------------------
## ****************************************************************************

## 2.1 Define basis parameter----
con <- gamlss.control(n.cyc = 500, trace = TRUE)
i.con <- glim.control(bf.cyc = 300)
age_seq <- 10:250/10
newdata <- data.frame(age = age_seq)
k <- log(length(g1_data_CRP_f$age))

## 2.2 First part-------
###  adding study-specific random effect 
mod_logRA_f <- gamlss(crp_le02 ~ pb(age, max.df = 8, method = "GAIC", k=k) + random(study_acronym),
                      data = g1_data_CRP_f, 
                      family = BI(), weights = W, control = con)

### first part of LOSO 
list_mods_1stpart <- LOSO_CV_mods_only(mod = mod_logRA_f,
                                       analysis_data = g1_data_CRP_f,
                                       ana_data_factor = g1_data_CRP_f$study_acronym)

save(list_mods_1stpart, file = "./30_Analysis/CRP/data/LOSO_crp_part1_f.RData")
# load("./30_Analysis/CRP/data/LOSO_crp_part1_f.RData")   

##  2.3 Second part ----------------------------

#set parameters
con <- gamlss.control(n.cyc = 5000, trace = TRUE)
i.con <- glim.control(bf.cyc = 300)
cent <- c(5, 50, 95)

k2 <- log(length(g2_data_CRP_f$age))

#2nd part of LOSO
list_mods_2ndpart <- LOSO_CV_mods_only(mod = mod_age_GB1_f_all2,
                                       analysis_data =g2_data_CRP_f ,
                                       ana_data_factor = g2_data_CRP_f$study_acronym)


save(list_mods_2ndpart, file = "./30_Analysis/CRP/data/LOSO_crp_part2_f.RData")
# load("./30_Analysis/CRP/data/LOSO_crp_part2_f.RData")   

## ****************************************************************************
# 3 Compile and plotting-----
## ****************************************************************************

## 3.1 Settings ----
shift_val <- 0.02
term_name_part1 <- 'pb(age, max.df = 8, method = "GAIC", k = k)'
term_name_part2 <- 'pb(age, max.df = 8, method = "GAIC", k = k2)'

# align folds across part1 & part2
left1 <- vapply(list_mods_1stpart, find_leftout, "", data = g1_data_CRP_f)
left2 <- vapply(list_mods_2ndpart, find_leftout, "", data = g2_data_CRP_f)
idx_ok <- which(!is.na(left1) & !is.na(left2) & left1 == left2)
stopifnot(length(idx_ok) > 0)
folds <- left1[idx_ok]

## 3.2  Plotting: P95 & P50, age-terms only ------------

svg("./30_Analysis/CRP/results/Fig_S12B_CRP_LOSO_f.svg")
op_mar <- par(mar = c(5, 4, 4, 5) + 0.1)  # 

# 
x_end <- 22.5

# empty frame
plot(NA,
     xlim = c(1, 24),     
     ylim = c(0, 1),
     xlab = "Age (yrs)", ylab = "hs-CRP (mg/L)",
     xaxt = "n", yaxt = "n", type = "n")
title("B  Females", adj = 0)
add_gridlines(x = 0:30, y = 0:10/10,
              P = 26, cent = cent, y_P = c(0.0, 0.1, 0.65))
# grid 
xg <- seq(0, 30, by = 1)  
yg <- seq(0, 1, by = 0.1)
axis(1, at = xg, labels = xg)
axis(2, at = yg, labels = yg * 10, las = 1)
abline(v = xg, col = "grey85", lwd = 0.5)
abline(h = yg, col = "grey85", lwd = 0.5)

# (mg/dL to nmol/L)
conv <- 95.2381                # 
yticks <- yg                   # 

raw_ticks <- c(0, 0.25, 0.50, 0.75, 1.00)
converted_ticks <- round(raw_ticks * conv, 1)
# 
axis(4,
     at     = raw_ticks,
     labels = round(converted_ticks, 1),
     las    = 1)
mtext("hs-CRP (nmol/L)", side = 4, line = 3)

### 3.2.1 Overall black line------------
ndf <- data.frame(
  age = age_seq,
  study_acronym = factor(levels(g1_data_CRP_f$study_acronym)[1],
                         levels = levels(g1_data_CRP_f$study_acronym))
  )

t1f <- predict(mod_logRA_f, newdata = ndf, type = "terms")
P1f <- plogis(t1f[, term_name_part1] + attr(t1f, "constant"))
t2f_mu <- predict(mod_age_GB1_f_all2, newdata = ndf, type = "terms", what = "mu")

muf    <- plogis(t2f_mu[, term_name_part2] + attr(t2f_mu, "constant")) #without random effect
sigmaf <- as.numeric(predict(mod_age_GB1_f_all2, newdata = ndf, type = "response", what = "sigma"))
nuf    <- as.numeric(predict(mod_age_GB1_f_all2, newdata = ndf, type = "response", what = "nu"))
tauf   <- as.numeric(predict(mod_age_GB1_f_all2, newdata = ndf, type = "response", what = "tau"))

P2f_50 <- Calc_Perc_2(P1f, 0.50)
Qf50 <- gamlss.dist::qGB1(P2f_50, muf, sigmaf, nuf, tauf) + 0.02
P2f_95 <- Calc_Perc_2(P1f, 0.95)
Qf95 <- gamlss.dist::qGB1(P2f_95, muf, sigmaf, nuf, tauf) + 0.02

# 
label_x <- x_end * 0.985

### 3.2.2 LOSO colorful lines till 22.5 ----------------

cols <- qualitative_hcl(13, "Dark 3")

downgraded <- character()

for (j in seq_along(idx_ok)) {
  i <- idx_ok[j]
  left_out <- folds[j]
  
  # part 1(terms)
  data_LOSO <- droplevels(subset(g1_data_CRP_f, study_acronym != left_out))
  nd1 <- data.frame(
    age = age_seq,
    study_acronym = factor(levels(data_LOSO$study_acronym)[1],
                           levels = levels(data_LOSO$study_acronym))
  )
  t1  <- predict(list_mods_1stpart[[i]], newdata = nd1, type = "terms")
  P1  <- plogis(t1[, term_name_part1] + attr(t1, "constant"))
  
  # part 2(mu: terms;others: response)
  data_LOSO <- droplevels(subset(g2_data_CRP_f, study_acronym != left_out))
  nd2 <- data.frame(
    age = age_seq,
    study_acronym = factor(levels(data_LOSO$study_acronym)[1],
                           levels = levels(data_LOSO$study_acronym))
  )
  tmu <- predict(list_mods_2ndpart[[i]], newdata = nd2, type = "terms", what = "mu")
  mu  <- plogis(tmu[, term_name_part2] + attr(tmu, "constant"))
  sigma <- as.numeric(predict(list_mods_2ndpart[[i]], newdata = nd2, type = "response", what = "sigma"))
  nu    <- as.numeric(predict(list_mods_2ndpart[[i]], newdata = nd2, type = "response", what = "nu"))
  tau   <- as.numeric(predict(list_mods_2ndpart[[i]], newdata = nd2, type = "response", what = "tau"))
  
  P2_50 <- Calc_Perc_2(P1, 0.50)
  Q50 <- gamlss.dist::qGB1(P2_50, mu, sigma, nu, tau) + 0.02
  P2_95 <- Calc_Perc_2(P1, 0.95)
  Q95 <- gamlss.dist::qGB1(P2_95, mu, sigma, nu, tau) + 0.02
  
  # 
  idx <- nd1$age <= x_end
  lines(nd1$age[idx], Q95[idx], col = cols[j], lwd = 2, lty = 2)
  lines(nd1$age[idx], Q50[idx], col = cols[j], lwd = 2, lty = 2)
}

# black lines
idx_f <- ndf$age <= x_end
lines(ndf$age[idx_f], Qf95[idx_f] , col = "black", lwd = 3)
lines(ndf$age[idx_f], Qf50[idx_f] , col = "black", lwd = 3)

# label of percentile
text(label_x, tail(Qf95[idx_f] , 1), "P95", pos = 4, cex = 0.85, xpd = TRUE, offset = 0.6)
text(label_x, tail(Qf50[idx_f] , 1), "P50", pos = 4, cex = 0.85, xpd = TRUE, offset = 0.6)

# legend
op_xpd <- par(xpd = NA)
legend("topleft",
       inset = 0.02,
       legend = c(folds, "Overall"),
       col    = c(cols, "black"),
       lwd    = c(rep(2, length(cols)), 3),
       lty    = 1,
       bty    = "o", box.lwd = 1.2, bg = "white",
       cex    = 0.75, y.intersp = 0.8, ncol = 2)
par(op_xpd)
par(op_mar)
dev.off()

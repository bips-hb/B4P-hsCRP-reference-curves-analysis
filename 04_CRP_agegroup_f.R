## ****************************************************************************
##
## Project:       Biomarkers4Pediatrics
## Program name:  04_CRP_agegroup_f.R
## Author:        Jiayi Zeng
## R-Version:     4.5.2
## version/ Date: V 1.0/ 2026-04-22
##
## Purpose:       Prevalence calculation of elevated CRP by applying 3 mg/L and age-/sex-specific cutoffs at age 1-22.5 (females)
##
## Program specification: 
## - Prerequisites: none
## - Read:  
##  "./20_Data/DataFiles/alspac_V_20260128.csv", 
##  "./20_Data/DataFiles/chns_V_20250807.csv", 
##  "./20_Data/DataFiles/degs1_V_20260424.csv",
##  "./20_Data/DataFiles/fuprecol_V_20250807.csv", 
##  "./20_Data/DataFiles/helena_V_20250807.csv", 
##  "./20_Data/DataFiles/kiggs_V_20250807.csv", 
##  "./20_Data/DataFiles/lsac_V_20250807.csv", 
##  "./20_Data/DataFiles/idefifam_V_20260129.csv", 
##  "./20_Data/DataFiles/sutas_V_20250807.csv", 
##  "./20_Data/DataFiles/usnhanes_V_20250807.csv", 
##  "./20_Data/DataFiles/knhanes_V_20250605.csv", 
##  "./20_Data/DataFiles/champs_V_20250807.csv", 
##  "./20_Data/DataFiles/mthatha_V_20250807.csv",
##  "./20_Data/DataFiles/ensanut_ecu_V_20250814.csv", 
##  "./20_Data/DataFiles/ensanut_V_20250829.csv"
## *Cutoff file
##"./30_Analysis/CRP/data/All_perc_comb_long_f.csv" [from 02_...]
##
## ***************************************************************************

## ***************************************************************************
# 1 Packages, sources, directories ------
## ***************************************************************************
library(readr)
library(haven)
library(dplyr)
library(ggplot2)
library(gamlss)
library(grDevices)
library(scales)
library(purrr)
library(writexl)
library(tidyverse)
library(readxl)
library(survey)

sessionInfo()
setwd("E:/Projects/B4P") 
source("./30_Analysis/CRP/codes/functions.R")

## ***************************************************************************
# 2 Data preparation------------------------------
## ***************************************************************************

selected_files <- c("./20_Data/DataFiles/alspac_V_20260128.csv", 
                    "./20_Data/DataFiles/chns_V_20250807.csv", 
                    "./20_Data/DataFiles/degs1_V_20260424.csv",
                    "./20_Data/DataFiles/fuprecol_V_20250807.csv", 
                    "./20_Data/DataFiles/helena_V_20250807.csv", 
                    "./20_Data/DataFiles/kiggs_V_20250807.csv", 
                    "./20_Data/DataFiles/lsac_V_20250807.csv", 
                    "./20_Data/DataFiles/idefifam_V_20260129.csv", 
                    "./20_Data/DataFiles/sutas_V_20250807.csv", 
                    "./20_Data/DataFiles/usnhanes_V_20250807.csv", 
                    "./20_Data/DataFiles/knhanes_V_20250605.csv", 
                    "./20_Data/DataFiles/champs_V_20250807.csv", 
                    "./20_Data/DataFiles/mthatha_V_20250807.csv",
                    "./20_Data/DataFiles/ensanut_ecu_V_20250814.csv", 
                    "./20_Data/DataFiles/ensanut_V_20250829.csv") 

list_data <- lapply(selected_files, function(f) {
  df <- read.csv(f)
  df$wave <- as.character(df$wave)
  df$study_id <- as.character(df$study_id)
  return(df)
})

full_data <- bind_rows(list_data)

## some basic checking for confirmation
unique(full_data$country)
unique(full_data$subregion)
unique(full_data$region)
unique(full_data$study_acronym)

## ensure agegroup2 exists (match function input)
full_data <- full_data %>%
  filter(age >= 1 & age <= 22.5)

full_data <- full_data %>%
  mutate(agegroup2 = case_when(
    age >= 1 & age < 5   ~ " 1-<5",
    age >= 5 & age < 10  ~ " 5-<10",
    age >= 10 & age < 15 ~ "10-<15",
    age >= 15 & age < 20 ~ "15-<20",
    age >= 20 & age <= 22.5 ~ "20-22.5"
  ))

## only keep wave 2009 of CHNS
full_data <- full_data %>%
  filter(study_acronym != "CHNS" | (study_acronym == "CHNS" & wave == 2009))

## CRP and age unification
### age
full_data <- full_data %>%
  filter(age >= 1) %>%
  mutate(age = round(age, 2))

### crp_si convert to mg/dL
full_data$crp_si_to_mgdL <- full_data$crp_si / 95.2381

### round crp_si_to_mgdL to 6 digits 
##for those with CRP available, use CRP as the value of crp_round and keep 6 digits;for those with only CRP_si available, use the calculated CRP_si/95.2381 as the value of crp_round and keep 6 digits
full_data <- full_data %>%
  mutate(
    crp_round = case_when(
      !is.na(crp) ~ round(crp, 6),
      !is.na(crp_si) & is.na(crp)~ round(crp_si_to_mgdL, 6),
      TRUE ~ NA_real_
    )
  )

### remove crp, crp_si_tomgdl, add final analysed variable crp
full_data <- full_data %>%
  select(-crp_si_to_mgdL, -crp, -X)%>% rename(crp = crp_round)

### select females
full_data <- subset(full_data, sex == 2)

## ****************************************************************************
#3 Merge the cutoff table with dataset  -------------------------------
## ****************************************************************************

## 3.1 Merging data preparation-------------
## load the dataset
cutoff <- read.csv("./30_Analysis/CRP/data/All_perc_comb_long_f.csv")[, c("global.age", "global.Perc_comb_90","global.Perc_comb_cutoff")]

##  merge 
full_data <- merge(full_data, cutoff, by.x = "age", by.y = "global.age", all.x = TRUE)

## generate crp_high variable
full_data$crp_high <- ifelse(full_data$crp >= full_data$global.Perc_comb_cutoff, 1, 0)
full_data$crp_p90 <- ifelse(full_data$crp >= full_data$global.Perc_comb_90, 1, 0)

## generate crp_high_0.3 variable
full_data$crp_0.3 <- ifelse(full_data$crp >= 0.3, 1, 0)

## total population with crp
ana_data_f <- full_data %>%
  filter(
    !is.na(crp))

## 3.2 Table of the percentage of elevated CRP in age groups------
by_age <- ana_data_f %>%
  group_by(agegroup2) %>%
  summarise(
    n = n(),
    prop_p89.9 = percent(mean(crp_high, na.rm = TRUE), accuracy = 0.01),
    prop_03    = percent(mean(crp_0.3,   na.rm = TRUE), accuracy = 0.01)
  )

overall <- ana_data_f %>%
  summarise(
    agegroup2  = "Overall",
    n          = n(),
    prop_p89.9 = percent(mean(crp_high, na.rm = TRUE), accuracy = 0.01),
    prop_03    = percent(mean(crp_0.3,   na.rm = TRUE), accuracy = 0.01)
  )

results <- bind_rows(by_age, overall)
write_xlsx(results, "./30_Analysis/CRP/results/Table2_cutoff_f.xlsx")
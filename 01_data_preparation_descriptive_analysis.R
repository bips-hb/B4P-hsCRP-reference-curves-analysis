## ****************************************************************************
##
## Project:       Biomarkers4Pediatrics
## Program name:  01_data_preparation_descriptive_analysis.R
## Author:        Jiayi Zeng
## R-Version:     4.5.2
## version/ Date: V 1.0/ 2026-04-22
##
## Purpose:        Data preparation for CRP reference curves + descriptive analysis
##
## Program specification: 
## - Prerequisites: none
##
## ***************************************************************************

## ***************************************************************************
# 1 Packages, sources, directories ------
## ***************************************************************************
library(readr)
library(patchwork)
library(haven)
library(dplyr)
library(ggplot2)
library(gamlss)
library(scales)
library(purrr)
library(abind)
library(writexl)
library(tidyverse)
library(readxl)
library(survey)
library(Hmisc)

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

##some basic checking for confirmation
unique(full_data$country)
unique(full_data$subregion)
unique(full_data$region)
unique(full_data$study_acronym)
unique(full_data$wave)

## 2.1 Study name changes--------------
full_data$study_acronym[full_data$study_acronym == "SUTAS"] <- "ZUTAS"
full_data$study_acronym[full_data$study_acronym == "SA_Hyper"] <- "MSALMSA"

## 2.2 Create agegroup2  (match function input)----
full_data <- full_data %>%
  mutate(agegroup2 = case_when(
    age >= 0 & age < 5   ~ "0-<5",
    age >= 5 & age < 10  ~ "5-<10",
    age >= 10 & age < 15 ~ "10-<15",
    age >= 15 & age < 20 ~ "15-<20",
    age >= 20 & age <= 25 ~ "20-25"
  ))

## 2.3 Set standard age group proportions--pop_agegroup2 for weighting------
pop_agegroup2 <- data.frame(
  agegroup2 = c("0-<5", "5-<10", "10-<15", "15-<20", "20-25"),
  Freq = c(0.2, 0.2, 0.2, 0.2, 0.2)
)
pop_agegroup2$agegroup2 <- as.factor(pop_agegroup2$agegroup2)

## 2.4 Only keep wave 2009 of CHNS-------
full_data <- full_data %>%
  filter(study_acronym != "CHNS" | (study_acronym == "CHNS" & wave == 2009))

## 2.5 Add variable comb_bmi_who-------------
full_data$comb_BMI_WHO <- ifelse(full_data$age <= 19,
                                 full_data$BMI_WHO_cat,
                                 full_data$BMI_cat)

## 2.6 Generate obesity/overweight variable ------
## obesity definition of who:0-<=19 BMI for age >2SD   | 19-25  BMI>= 30
## overweight definition of who: 0-<=19 BMI > 1SD|  19-25  BMI >=25
## generate obesity variable
full_data <- full_data %>%
  mutate(
    obesity = case_when(
      age <= 19 & !is.na(z_BMI_WHO) & z_BMI_WHO > 2 ~ TRUE,
      age > 19  & !is.na(bmi) & bmi >= 30 ~ TRUE,
      TRUE ~ FALSE
    )
  )

## generate overweight variable
full_data <- full_data %>%
  mutate(
    overweight = case_when(
      age <= 19 & !is.na(z_BMI_WHO) & z_BMI_WHO > 1 & z_BMI_WHO <= 2 ~ TRUE,
      age > 19  & !is.na(bmi) & bmi >= 25 & bmi < 30 ~ TRUE,
      TRUE ~ FALSE
    )
  )

##generate normal weight and underweight variables
full_data<- full_data %>%
  mutate(
  normal_weight = case_when(
  age <= 19 & !is.na(z_BMI_WHO) & z_BMI_WHO >= -2 & z_BMI_WHO <= 1 ~ TRUE,
  age > 19 & !is.na(bmi) & bmi >= 18.5 & bmi < 25 ~ TRUE,
  TRUE ~ FALSE
),
  underweight = case_when(
  age <= 19 & !is.na(z_BMI_WHO) & z_BMI_WHO < -2 ~ TRUE,
  age > 19 & !is.na(bmi) & bmi < 18.5 ~ TRUE,
  TRUE ~ FALSE
)
)

## 2.7 IL-6 and CRP both high population------
# generate il-6 high  variable
full_data <- full_data %>%
  mutate(
    il_6_high = ifelse(is.na(il_6_high), FALSE, il_6_high),
  )

# define the studies with IL-6 data
studies_with_il6 <- c("HELENA", "IDEFICS/I.FAMILY", "SUTAS", "ALSPAC")

# circulate each study
for (study in studies_with_il6) {
  # calculate 90p
  crp_90 <- quantile(full_data$crp_si[full_data$study_acronym == study], 0.9, na.rm = TRUE)
  il6_90 <- quantile(full_data$il_6[full_data$study_acronym == study], 0.9, na.rm = TRUE)
  
  # mark the reason of exclusion
  full_data <- full_data %>%
    mutate(
      il_6_high = ifelse(
        study_acronym == study & crp_si > crp_90 & il_6 > il6_90,
        TRUE,
        il_6_high
      )
    )
}

## ***************************************************************************
# 3 Reference population------
## ***************************************************************************

## 3.1 CRP variable preparation--------
## crp_si convert to mg/dL
full_data$crp_si_to_mgdL <- full_data$crp_si / 95.2381

## for those with CRP available, use CRP as the value of crp_round and keep 6 digits;for those with only CRP_si available, use the calculated CRP_si/95.2381 as the value of crp_round and keep 6 digits
full_data <- full_data %>%
  mutate(
    crp_round = case_when(
      !is.na(crp) ~ round(crp, 6),
      !is.na(crp_si) & is.na(crp)~ round(crp_si_to_mgdL, 6),
      TRUE ~ NA_real_
    )
  )

## remove crp, crp_si_tomgdL, add final analysed variable crp
full_data <- full_data %>%
  select(-crp_si_to_mgdL, -crp, -X) %>%
  dplyr::rename(crp = crp_round)
full_data$crp_le02 <- full_data$crp<=0.02 
full_data %>% distinct(b4p_id) %>% nrow() ## 404852 observations with 310835 n

## 3.2 Exclusion criteria: acute disease, crp over 1 mg/dL, obesity, il-6 high and crp high---------

##   (109558 n with 131981 observations) with crp

## 123857 observations  n = 103861 CRP and bmi/age/sex // 
ref_data <- full_data %>%
  filter(
    !is.na(crp)
  )
ref_data %>% distinct(b4p_id) %>% nrow()

ref_data <- ref_data %>%
  filter(
    !is.na(age),
    !is.na(sex),
    !is.na(bmi)
  )
ref_data %>% distinct(b4p_id) %>% nrow()

## 119900 with crp <1 / n = 100710 / descreased 3957 n=3151
ref_data %>%
  filter(crp>= 1) %>%
  count(study_acronym, sort = TRUE)

ref_data <- ref_data %>%
  filter(
    crp < 1
  )
ref_data %>% distinct(b4p_id) %>% nrow()

##  with 105452 non obesity / n = 88301 / descreased 14448 n=12409 
ref_data %>%
  filter(obesity== TRUE) %>%
  count(study_acronym, sort = TRUE)

ref_data <- ref_data %>%
  filter(
    !obesity
  )
ref_data %>% distinct(b4p_id) %>% nrow()

## with 103739 non acute disease / n = 87731 / descreased 1713 n=570
ref_data %>%
  filter(acute_disease== TRUE) %>%
  count(study_acronym, sort = TRUE)

ref_data <- ref_data %>%
  filter(
    (is.na(acute_disease) | acute_disease == FALSE)
  )
ref_data %>% distinct(b4p_id) %>% nrow()

## with 101962 no current medication / n = 86960 / descreased 1777 n=771
ref_data %>%
  filter(current_medic== TRUE) %>%
  count(study_acronym, sort = TRUE)

ref_data <- ref_data %>%
  filter(
    (is.na(current_medic) | current_medic == FALSE)
  )
ref_data %>% distinct(b4p_id) %>% nrow()

## with 100079 non simutaneously increasing il-6 / n = 86420 / descreased 1883 n=540
ref_data %>%
  filter(il_6_high== TRUE) %>%
  count(study_acronym, sort = TRUE)

ref_data <- ref_data %>%
  filter(
    (is.na(il_6_high) | il_6_high == FALSE)
  )
ref_data %>% distinct(b4p_id) %>% nrow()

## 3.3 Sex-specific outlier detection (after the reference population definition)---------
## result: 0
data_male <- ref_data[ref_data$sex == 1, ]
data_female   <- ref_data[ref_data$sex == 2, ]

excl_f <- flag_mod_based_outlier(data_female$crp, data_female$age, data_female)
excl_m <- flag_mod_based_outlier(data_male$crp, data_male$age, data_male)

table(excl_m)
table(excl_f)

## ***************************************************************************
# 4 Limit of detection (LOD) -----
## ***************************************************************************
##check for every wave of each study about the LOD percentage

## 4.1 Add wave-specific LOD ------
study_wave_lod <- ref_data %>%
  distinct(study_acronym, wave) 

study_wave_lod <- study_wave_lod %>%
  mutate(
    lod_value = case_when(
      study_acronym == "ALSPAC"  ~ 0.010,
      study_acronym == "CHNS"  ~ 0.010,
      study_acronym == "DEGS" ~ 0.015,
      study_acronym == "HELENA" ~ 0.0007,
      study_acronym == "IDEFICS/I.FAMILY" & wave =="t0" ~ 0.016,
      study_acronym == "IDEFICS/I.FAMILY" & wave =="t1" ~ 0.016,
      study_acronym == "IDEFICS/I.FAMILY" & wave =="t3" ~ 0.00001,
      study_acronym == "KIGGS" & wave =="2" ~ 0.010,
      study_acronym == "KIGGS" & wave =="0" ~ 0.020,
      study_acronym == "KNHANES" & wave =="15" ~ 0.010,
      study_acronym == "KNHANES" & wave =="16" ~ 0.010,
      study_acronym == "KNHANES" & wave =="17" ~ 0.015,
      study_acronym == "KNHANES" & wave =="18" ~ 0.015,
      study_acronym == "KNHANES" & wave =="22" ~ 0.020,
      study_acronym == "KNHANES" & wave =="23" ~ 0.020,
      study_acronym == "LSAC"  ~ 0.001,
      study_acronym == "ZUTAS"  ~ 0.002,
      study_acronym == "USNHANES" & wave %in% c("A", "B", "C", "D", "E", "F") ~ 0.020,
      study_acronym == "USNHANES" & wave %in% c("P", "L") ~ 0.015,
      study_acronym == "USNHANES" & wave =="I" ~ 0.011,
      study_acronym == "CHAMPS" ~ 0.015,
      study_acronym == "MSALMSA" ~ 0.015,
      study_acronym == "ENSANUT" & wave =="2006" ~ 0.017,
      TRUE ~ NA_real_  # others tagged as NA
    )
  )

ref_data <- ref_data %>%
  left_join(study_wave_lod, by = c("study_acronym", "wave"))

ref_data %>%
  filter(is.na(lod_value)) %>%
  count(study_acronym, wave)

sum(!is.na(ref_data$lod_value))

## 4.2 LOD percentage study-specific------
lod_summary_table_new <- ref_data %>%
  filter(!is.na(lod_value)) %>%    # only those with LOD info
  group_by(study_acronym, wave) %>%
  summarise(
    lod_value = first(lod_value),
    n_total = n(),
    n_leq_lod = sum(crp <= lod_value),
    pct_leq_lod = mean(crp <= lod_value) * 100,
    .groups = "drop"
  ) %>%
  rename(
    `n(<=LOD)` = n_leq_lod,
    `%(<=LOD)` = pct_leq_lod
      )

## ***************************************************************************
# 5 Measurement methods------
## ***************************************************************************

measurement_methods <- ref_data %>%
  distinct(study_acronym, wave) 

measurement_methods <- measurement_methods %>%
  mutate(
    methods = case_when(
      study_acronym == "ALSPAC"  ~ "Immunoturbidimetry",
      study_acronym == "CHNS"  ~ "Immunoturbidimetry",
      study_acronym == "DEGS" ~ "Nephelometry",
      study_acronym == "HELENA" ~ "Immunoturbidimetry",
      study_acronym == "IDEFICS/I.FAMILY" & wave =="t0" ~ "Nephelometry",
      study_acronym == "IDEFICS/I.FAMILY" & wave =="t1" ~ "Nephelometry",
      study_acronym == "IDEFICS/I.FAMILY" & wave =="t3" ~ "Chemiluminescence",
      study_acronym == "KIGGS"  ~ "Immunoturbidimetry",
      study_acronym == "KNHANES" ~ "Immunoturbidimetry",
      study_acronym == "LSAC"  ~ "Immunoturbidimetry",
      study_acronym == "ZUTAS"  ~ "Chemiluminescence",
      study_acronym == "USNHANES" & wave %in% c("A", "B", "C", "D", "E", "F") ~ "Nephelometry",
      study_acronym == "USNHANES" & wave %in% c("I", "P", "L") ~ "Immunoturbidimetry",
      study_acronym == "FUPRECOL" ~ "Immunoturbidimetry",
      study_acronym == "CHAMPS" ~ "Immunoturbidimetry",
      study_acronym == "MSALMSA" ~ "Immunoturbidimetry",
      study_acronym == "ENSANUT_ECU" ~ "Nephelometry",
      study_acronym == "ENSANUT" & wave == "2018"~ "Chemiluminescence",
      study_acronym == "ENSANUT"& wave %in% c("2006", "2021", "2022")~ "Nephelometry",
      TRUE ~ NA_character_  # others tagged as NA
    )
  )
ref_data <- ref_data %>%
  left_join(measurement_methods, by = c("study_acronym", "wave"))

## ***************************************************************************
#6 LOD study-specific assigning--------------
## ***************************************************************************

## 6.1 replace those lower than LOD and 0 values with the study specific LOD---------

ref_data <- ref_data %>%
  mutate(crp = if_else(!is.na(crp) & !is.na(lod_value) & crp <= lod_value,
                       lod_value, crp))

## 6.2 Two zero values with unknown LOD from ENSANUT 2018-- observation = 100077 n=86418--------
ref_data %>%
  filter(ref_data$crp == 0) %>%
  count(study_acronym, wave, sort = TRUE)

ref_data %>%
  filter(crp == 0, study_acronym == "ENSANUT") %>%
  select(b4p_id, age, crp) 

ref_data<- ref_data[ref_data$crp > 0, ]

ref_data %>% distinct(b4p_id) %>% nrow()

## ***************************************************************************
#7 Weighting ------
## ***************************************************************************

## 7.1 Delete ecu and fuprecol  ------------
## observations decreased by 11198 = 100077 - 88879 
## n decreased by 11198 = 86418-75220
ref_data_without_ecu_fu <- subset(ref_data, !(study_acronym %in% c("ENSANUT_ECU", "FUPRECOL")))
ref_data_without_ecu_fu %>% distinct(b4p_id) %>% nrow()

## 7.2 Delete two participants under age 1 ------------
ref_data_without_ecu_fu %>%
  filter(ref_data_without_ecu_fu$age< 1) %>%
  count(study_acronym, sort = TRUE)

ref_data_without_ecu_fu <- ref_data_without_ecu_fu[ref_data_without_ecu_fu$age >= 1, ]

## final dataset observation=88877 n=75218
df_main <- ref_data_without_ecu_fu 
ref_data_without_ecu_fu %>% distinct(b4p_id) %>% nrow()

## 7.3 Load population reference data (3 levels)----
pop_country <- read.csv("./20_Data/DataFiles/additional_data/add_data_hdi_V_20250602.csv")
pop_country <- pop_country[pop_country$year == 2008, c("country", "pop_total")]
names(pop_country) <- c("country", "pop")

pop_region <- read.csv("./20_Data/DataFiles/additional_data/add_data_ncdrisc_region_V_20250324.csv")
pop_region <- pop_region[, c("region", "pop_total")]
names(pop_region) <- c("region", "pop")

pop_subregion <- read.csv("./20_Data/DataFiles/additional_data/add_data_ncdrisc_subregion_V_20250324.csv")
pop_subregion <- pop_subregion[, c("subregion", "pop_total")]
names(pop_subregion) <- c("subregion", "pop")

# change it to variable + Freq 
pop_country   <- pop_country   |> dplyr::transmute(country,  Freq = pop)
pop_region    <- pop_region    |> dplyr::transmute(region,   Freq = pop)
pop_subregion <- pop_subregion |> dplyr::transmute(subregion,Freq = pop)

## Manually choose one of the three for sensitivity comparison
pop_subreg_var <- pop_subregion   # or pop_region / pop_country

## automatically filter the existing sub region data in ref_data
pop_subreg_var <- pop_subreg_var %>%
  dplyr::filter(subregion %in% unique(df_main$subregion))

##
df_main$subregion <- factor(as.character(df_main$subregion),
                            levels = sort(unique(pop_subreg_var$subregion)))

## 7.4 Run weighting-----
##new weighting dataset
weighted_data <- survey_weighting_control_opt(
  df = df_main,
  pop_subreg_var = pop_subreg_var,
  control = list(maxit = 40)  
)

# Optional: quick check
summary(weighted_data$W)

# weighting comparison tables
target_age <- pop_agegroup2 %>%
  mutate(target_pct = Freq / sum(Freq) * 100) %>%
  select(agegroup2, target_pct)

target_subregion <- pop_subreg_var %>%
  mutate(target_pct = Freq / sum(Freq) * 100) %>%
  select(subregion, target_pct)

age_table <- weighted_data %>%
  group_by(agegroup2) %>%
  summarise(
    unweighted_n = n(),
    weighted_n = sum(W, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    unweighted_pct = unweighted_n / sum(unweighted_n) * 100,
    weighted_pct = weighted_n / sum(weighted_n) * 100
  ) %>%
  left_join(target_age, by = "agegroup2")

subregion_table <- weighted_data %>%
  mutate(subregion = as.integer(as.character(subregion))) %>% 
  group_by(subregion) %>%
  summarise(
    unweighted_n = n(),
    weighted_n = sum(W, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    unweighted_pct = unweighted_n / sum(unweighted_n) * 100,
    weighted_pct = weighted_n / sum(weighted_n) * 100
  ) %>%
  left_join(target_subregion, by = "subregion")

age_table_out <- age_table %>%
  transmute(
    dimension = "Age group (5-year)",
    category  = as.character(agegroup2),
    unweighted_pct,
    weighted_pct,
    target_pct
  )

subregion_table_out <- subregion_table %>%
  transmute(
    dimension = "Subregion",
    category  = as.character(subregion),
    unweighted_pct,
    weighted_pct,
    target_pct
  )

weighting_diagnostic_table <- bind_rows(
  age_table_out,
  subregion_table_out
)

## save document
saveRDS(weighted_data, file = "./30_Analysis/CRP/data/ana_data_20260424.rds") 

lod_summary_table_new<-lod_summary_table_new %>%
  mutate(lod_mg_L=lod_value*10)

write.csv(weighting_diagnostic_table,
          "./30_Analysis/CRP/results/Tab_S4_weighting_diagnostics.csv",
          row.names = FALSE)

write.csv(lod_summary_table_new,
          "./30_Analysis/CRP/results/Tab_S3_wave_lod.csv",
          row.names = FALSE)

## ***************************************************************************
#8 Descriptive tables-----------------------------------
## ***************************************************************************

## 8.1 Tab_1_b4p_sample_size_comparison-------------
obs_count <- tapply(full_data$b4p_id, full_data$study_acronym, length) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("study_acronym") %>%
  rename(observation_N = ".")

unique_count <- tapply(ref_data$b4p_id, ref_data$study_acronym, function(x) length(unique(x))) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("study_acronym") %>%
  rename(unique_participants = ".")

merged_table <- full_join(obs_count, unique_count, by = "study_acronym")

write.csv(merged_table, "./30_Analysis/CRP/results/Tab_1_b4p_sample_size_comparison.csv", row.names = FALSE)

#change sex from number to factors
ref_data$sex <- factor(ref_data$sex, levels = c(1, 2), labels = c("Male", "Female"))

## 8.2 Heat map----------------------
# agegroup2 factor order
weighted_data$agegroup2 <- factor(
  weighted_data$agegroup2,
  levels = c("0-<5", "5-<10", "10-<15", "15-<20", "20-25")
)

###8.2.1 Fig_S2_age_study_heatmap----------------
p1 <- age_country_plot(weighted_data, " ")
svg("./30_Analysis/CRP/results/Fig_S2_age_study_heatmap.svg", width = 8, height = 6)
print(p1)
dev.off()

###8.2.2 Fig_S3_age_country_heatmap----------------------

p2 <- age_study_plot(weighted_data, "")

svg("./30_Analysis/CRP/results/Fig_S3_age_country_heatmap.svg", width = 8, height = 6)
print(p2)
dev.off()

## 8.3 CRP distribution--------------- 

###8.3.1 Fig_S4_hist_each_study--------------
# base R, Dec 1st
svg("./30_Analysis/CRP/results/Fig_S4_hist_each_study.svg", width = 12, height = 8)

par(
  mfrow = c(3, 5),
  mar = c(3, 3, 2, 1),
  oma = c(0, 0, 3, 1),
  mgp = c(1.5, 0.5, 0),
  lwd = 1.3,
  cex.axis = 1.2,
  cex.lab  = 1.1,
  font.lab = 1.1,
  font.axis = 1.2
)

for (z in unique(ref_data$study_acronym)) {
  
  x <- ref_data$crp[ref_data$study_acronym == z]   # mg/dL
  
  hist(x,
       main = paste("Study:", z),
       cex.main = 1.3,
       xlab = "CRP (mg/L)",    # mg/L label
       col = "skyblue", border = "white",
       freq = FALSE,
       breaks = seq(0, 1, by = 0.02),   
       xlim   = c(0, 1),
       ylim   = c(0, 30),
       xaxt = "n")   # 
  
  # 
  axis(1,
       at = seq(0, 1, by = 0.2),     
       labels = seq(0, 10, by = 2))  
  
  #0.02 mg/dL = 0.2 mg/L
  abline(v = 0.02, col = "red", lty = 2, lwd = 1.2)
}

dev.off()

###8.3.2 Tab_S6_studywise_crp_summary-------------------
crp_summary <- weighted_data %>%
  group_by(study_acronym, sex) %>%
  summarise(
    Min    = min(crp, na.rm = TRUE) * 10,
    Q1     = quantile(crp, 0.25, na.rm = TRUE) * 10,
    Median = median(crp, na.rm = TRUE) * 10,
    Mean   = mean(crp, na.rm = TRUE) * 10,
    Q3     = quantile(crp, 0.75, na.rm = TRUE) * 10,
    Max    = max(crp, na.rm = TRUE) * 10,
    N      = n()
  )

write_xlsx(crp_summary, "./30_Analysis/CRP/results/Tab_S6_studywise_crp_summary.xlsx")

##8.4 Tab_S5_characteristics--------------
###8.4.1 Finalize data for descriptive tables -----------
dat0 <- weighted_data %>%
  filter(!is.na(b4p_id), !is.na(W)) %>%
  mutate(
    sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    region = factor(region, levels = 1:7,
                    labels = c(
                      "East Asia and Pacific", "South Asia",
                      "Central & Eastern Europe and Central Asia",
                      "North Africa and Middle East", "Sub-Saharan Africa",
                      "Latin America and Caribbean", "High-income regions"
                    )) %>%droplevels(),
    subregion = factor(
            subregion,
            levels = c(
              11, 12, 13, 21,
              31, 32, 33, 41,
              51, 52, 53, 54,
              61, 62, 63, 64, 65,
              71, 72, 73, 74
            ),
            labels = c(
              "Southeast Asia",
              "East Asia",
              "Oceania",
              "South Asia",
              "Central Asia",
              "Central Europe",
              "Eastern Europe",
              "North Africa and Middle East",
              "Central Africa",
              "East Africa",
              "Southern Africa",
              "West Africa",
              "Andean Latin America",
              "Central Latin America",
              "Southern Latin America",
              "Tropical Latin America",
              "Caribbean",
              "Asia-Pacific (high-income)",
              "Australasia",
              "Western Europe",
              "North America (high-income)"
            )
          ) %>%droplevels(),
    country = factor(country),
    study_acronym = factor(study_acronym),
    agegroup2 = factor(agegroup2,
                       levels = c("0-<5", "5-<10", "10-<15", "15-<20", "20-25"))
  )


### 8.4.2 Combine Male / Female / All for table S5 ----
tab_m <- one_col(filter(dat0, sex=="Male"))   %>% rename(Male_unw=Unweighted, Male_w=Weighted)
tab_f <- one_col(filter(dat0, sex=="Female")) %>% rename(Female_unw=Unweighted, Female_w=Weighted)
tab_a <- one_col(dat0)                        %>% rename(All_unw=Unweighted,   All_w=Weighted)

table1 <- tab_m %>%
  full_join(tab_f, by="Characteristics") %>%
  full_join(tab_a, by="Characteristics")

### 8.4.3 Export ----
write.csv(table1, "./30_Analysis/CRP/results/Tab_S5_characteristics.csv", row.names = FALSE)
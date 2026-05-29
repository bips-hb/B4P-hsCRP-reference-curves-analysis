## ****************************************************************************
##
## Project:       Biomarkers4Pediatrics
## Program name:  functions.R
## Authors:       Timm Intemann, Jiayi Zeng
## R-Version:     4.5.2
## version/ Date: V 1.0/ 2026-04-28
##
## Purpose:       Store functions for analyzing reference curves
##
## ***************************************************************************

## ****************************************************************************
# 1 Exclusion: Plausis and outlier  ---------------
## ****************************************************************************

## 1.1 Model based Outlier   --------------------
# Description: Function to flag model-based outliers
# Arguments: 
## y: outcome variable
## x: predictor variable (age)
## dataset: dataset used for model
## lower: lower residual threshold
## upper: upper residual threshold

flag_mod_based_outlier <- function(y, x, dataset, lower = -3.5, upper = 10){
  
  mod <- lm(y ~ bs(x, 4), data = dataset)

x_mod_out_f <- ((mod$residuals/summary(mod)$sigma) > upper) + 
  ((mod$residuals/summary(mod)$sigma) < (lower))
return(x_mod_out_f)
}

## ****************************************************************************
# 2 Descriptive analysis  -------------
## ****************************************************************************

## 2.1 Age and country rays --------
# Description: Function to generate plot describing sample size by age group and country
# Arguments: 
## dataset: dataset used for function
## title_text: title text

age_country_plot <- function(dataset, title_text = ""){
n_country_agegroup <- dataset %>% group_by(country, agegroup2) %>% dplyr::summarise( 
  N_obs = sum(!is.na(age)))

plot <- ggplot(n_country_agegroup, aes(x = agegroup2, y = country, fill = N_obs)) +
  geom_tile(color = "white", size = 0.5) + # Create tiles with borders
  geom_text(aes(label = N_obs), color = "black", size = 4) + # Add text for sample size
  scale_fill_gradient(
    low = "white", # Brightest color
    high = "steelblue" # Darker, but not too dark
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, hjust = 0),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    legend.position = "none", # Remove legend
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = title_text,
    x = "Age group (yrs)",
    y = "Country"
  )
return(plot)
}

## 2.2 Age and study rays --------
# Description: Function to generate plot describing sample size by age group and study
# Arguments: 
## dataset: dataset used for function
## title_text: title text

age_study_plot <- function(dataset, title_text = " "){
  n_study_agegroup <- dataset %>% group_by(study_acronym, agegroup2) %>% dplyr::summarise( 
    N_obs = sum(!is.na(age)))
  
  # Create the plot
  plot <- ggplot(n_study_agegroup, aes(x = agegroup2, y = study_acronym, fill = N_obs)) +
    geom_tile(color = "white", size = 0.5) + # Create tiles with borders
    geom_text(aes(label = N_obs), color = "black", size = 4) + # Add text for sample size
    scale_fill_gradient(
      low = "white", 
      high = "steelblue" 
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 10, hjust = 0),
      axis.text.y = element_text(size = 12),
      panel.grid = element_blank(),
      legend.position = "none",
      strip.text = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = title_text,
      x = "Age group (yrs)",
      y = "Study"
    )
  return(plot)
}

## *********************************************************************************************
# 3 Weighting ----------------------------------
## *********************************************************************************************
# Description: Function to generate survey weights by participant, age group and subregion
# Arguments
## df: dataset used for weighting
## pop_subreg_var: predefined subregion population margins
## control: control options for ranking procedures

survey_weighting_control_opt <- function(df, pop_subreg_var, control = list(partial = TRUE)){
  
  #preparation
  n_participant <- df %>% group_by(b4p_id)  %>% dplyr::summarise( 
    N_obs = sum(!is.na(age)), .groups="keep")
  
  df <- merge(df, n_participant, by = "b4p_id", all.x = TRUE)
  df$W_participant <- 1/df$N_obs
  
  #as factors
  df$subregion <- as.factor(df$subregion)
  df$agegroup2 <- as.factor(df$agegroup2)
  
  #starting design - each participants weights same
  design <- svydesign(ids = ~1, data = df, weights = ~W_participant)
  
  # Rake to the predefined marginals
  design_raked <- rake(
    design,
    sample.margins = list(~agegroup2, ~subregion),
    population.margins = list(pop_agegroup2, pop_subreg_var),
    control = control
  )
  
  # Normalize weights
  w <- weights(design_raked)
  w_norm <- w * (length(w) / sum(w))
  design_raked_norm <- svydesign(ids = ~1, data = df, weights = ~w_norm)
  
  #Trim weights
  design_raked_norm_trim <- trimWeights(design_raked_norm, lower = 0.1, upper = 10, strict = TRUE)
  
  #save
  df$W <- weights(design_raked_norm_trim)
  return(df)
}


## *****************************************************************************
# 4 Analysis / plotting -------------------------------------------------------
## *****************************************************************************

# Description: Function to add global percentile curves from a GAMLSS with a Generalized Gamma distribution (GG) 
# Arguments: 
## mod: fitted gamlss model
## newdata: dataset used for prediction
## cent: percentiles to be plotted
## mu_term_name: name of the mu model term
## sigma_term_name: name of the sigma model term
## color_curve: curve color
## lwd: line width

add_plot_global_curves_GG <- function(mod, newdata, cent, 
                                      mu_term_name = "age", sigma_term_name = "age",
                                      col_curve = "darkgrey", lwd = 3){
  
  mu_terms <- predict(mod, newdata = newdata, type = "terms", what ="mu")
  mu <- exp(attributes(mu_terms)$constant + mu_terms[, mu_term_name])
  
  if (length(mod$sigma.coefficients)<=2){
    sigma <- predict(mod, newdata = newdata, type = "response", what ="sigma")
  }
  if (length(mod$sigma.coefficients)>=3){
    sigma_terms <- predict(mod, newdata = newdata, type = "terms", what ="sigma")
    sigma <- exp(attributes(sigma_terms)$constant + sigma_terms[, sigma_term_name])
  }
  
  nu <- predict(mod, newdata = newdata, type = "response", what ="nu")
  
  quants <- qGG(p = rep(cent/100, each = length(mu)), 
                mu = rep(mu, times = length(cent)), 
                sigma = rep(sigma, times = length(cent)), 
                nu = rep(nu, times = length(cent))
  )
  quants <- matrix(quants, ncol = length(cent))
  
  for( i in 1:length(cent)){
    lines(newdata$age, quants[, i], lwd = lwd, col = col_curve) 
  }
}


# Description: Function to add study specific percentile curves from a GAMLSS with a generalized Beta type 1 distribution (GB1) 
# Arguments: 
## mod: fitted gamlss model
## newdata: dataset used for prediction
## cent: percentiles to be plotted
## list_studies: studies to be plotted
## lwd: line width
## color_vectors: colors for study-specific curves

add_plot_all_study_curves_GB1 <- function(mod, 
                                          newdata, 
                                          cent, 
                                          list_studies, 
                                          lwd = 2, 
                                          color_vector){
  color_vector <- color_vector
  col <- color_vector[1]
  col_num <- 1
  for (j in list_studies){
    newdata$study_acronym <- j
    col <- color_vector[col_num]
    add_plot_study_curve_GB1(mod = mod, newdata = newdata, 
                         cent = cent, col_curve = col, lwd = lwd)
    col_num <- col_num + 1
  }

}
 

# Description: Function to add global percentile curves from a GAMLSS with a generalized Beta type 1 distribution (GB1) 
# Arguments: 
## mod: fitted gamlss model
## newdata: dataset used for prediction
## cent: percentiles to be plotted
## mu_term_name: name of the mu model term
## sigma_term_name: name of the sigma model term
## color_curve: curve color
## lwd: line width

add_plot_global_curves_GB1 <- function(mod, 
                                       newdata, 
                                       cent, 
                                       mu_term_name = "age", 
                                       sigma_term_name = "age",
                                       col_curve = "darkgrey", 
                                       lwd = 3){
  
  mu_terms <- predict(mod, newdata = newdata, type = "terms", what ="mu")
  mu <-plogis(attributes(mu_terms)$constant + mu_terms[, mu_term_name])
  
  if (length(mod$sigma.coefficients)<=2){
    sigma <- predict(mod, newdata = newdata, type = "response", what ="sigma")
  }
  if (length(mod$sigma.coefficients)>=3){
    sigma_terms <- predict(mod, newdata = newdata, type = "terms", what ="sigma")
    sigma <- plogis(attributes(sigma_terms)$constant + sigma_terms[, sigma_term_name])
  }
  
  nu <- predict(mod, newdata = newdata, type = "response", what ="nu")
  tau <- predict(mod, newdata = newdata, type = "response", what ="tau")
  
  quants <- qGB1(p = rep(cent/100, each = length(mu)), 
                 mu = rep(mu, times = length(cent)), 
                 sigma = rep(sigma, times = length(cent)), 
                 nu = rep(nu, times = length(cent)), 
                 tau = rep(tau, times = length(cent))
  )
  quants <- matrix(quants, ncol = length(cent))
  
  for( i in 1:length(cent)){
    lines(newdata$age, quants[, i], lwd = lwd, col = col_curve) 
  }
}

# Description: Function to add a study specific percentile curve from a GAMLSS with a generalized Beta type 1 distribution (GB1) 
# Arguments: 
## mod: fitted gamlss model
## newdata: dataset used for prediction
## cent: percentiles to be plotted
## lwd: line width
## color_curve: curve color

add_plot_study_curve_GB1 <- function(mod, newdata, cent, 
                                     col_curve = "darkgrey", lwd = 2){
  
  mu <- predict(mod, newdata = newdata, type = "response", what ="mu")
  sigma <- predict(mod, newdata = newdata, type = "response", what ="sigma")
  nu <- predict(mod, newdata = newdata, type = "response", what ="nu")
  tau <- predict(mod, newdata = newdata, type = "response", what ="tau")
  
  quants <- qGB1(p = rep(cent/100, each = length(mu)), 
                 mu = rep(mu, times = length(cent)), 
                 sigma = rep(sigma, times = length(cent)), 
                 nu = rep(nu, times = length(cent)), 
                 tau = rep(tau, times = length(cent))
  )
  quants <- matrix(quants, ncol = length(cent))
  
  for( i in 1:length(cent)){
    lines(newdata$age, quants[, i], lwd = lwd, lty = 2, col = col_curve) 
  }
}


# Description: Function to create a weighted worm plot
# Arguments: 
## res: model residuals
## wgt: weights
## mod_name_title: model names shown in the plot title

wgt_wp <- function(res, wgt, mod_name_title = ""){
  wgt_residuals <- data.frame(res = res, 
                            wgt = wgt)

wgt_residuals <- wgt_residuals[order(wgt_residuals$res), ]
wgt_residuals$prob <-  wgt_residuals$wgt/(sum(wgt_residuals$wgt) + 1)
wgt_residuals$cum <-  cumsum(wgt_residuals$prob)
wgt_residuals$qNOcum <- qnorm(wgt_residuals$cum)
wgt_residuals$diff <- wgt_residuals$res - wgt_residuals$qNOcum

plot(wgt_residuals$qNOcum, wgt_residuals$diff, 
     ylim = c(-0.15, 0.15), xlim = c(-4, 4),
     main = paste("Weighted worm plot", mod_name_title),
     xlab = "Unit normal weighted quantile",
     ylab = "Deviation")

lines(c(-100, 100), c(0, 0),  col = "red", lty = 2)
lines(c(-100, 100), c(-0.1, -0.1), col = "grey")
lines(c(-100, 100), c(0.1, 0.1), col = "grey")

lines(c(0, 0),c(-100, 100),  col = "red", lty = 2)
lines(c(1, 1),c(-100, 100),  col = "grey")
lines(c(2, 2),c(-100, 100),  col = "grey")
lines(c(3, 3),c(-100, 100),  col = "grey")
lines(c(4, 4),c(-100, 100),  col = "grey")

lines(c(-1, -1),c(-100, 100),  col = "grey")
lines(c(-2, -2),c(-100, 100),  col = "grey")
lines(c(-3, -3),c(-100, 100),  col = "grey")
lines(c(-4, -4),c(-100, 100),  col = "grey")

model <- lm(wgt_residuals$diff ~ poly(wgt_residuals$qNOcum, 3, raw = TRUE), weights = wgt_residuals$wgt )
lines(wgt_residuals$qNOcum, 
      predict(model, newdata = wgt_residuals), col = "red", lwd = 2)

a <- model$coef

names(a) <- c("Intercept", "x^1", "x^2", "x^3")

# format coefficients
a_fmt <- formatC(a, format = "f", digits = 3)
# build expression text
coef_text <- paste(names(a), "=", a_fmt, collapse = "   ")
# add as margin text under x-axis
mtext(coef_text, side = 1, line = 4, cex = 0.9)

return(a)
}


# Description: Function to add gridlines and percentile labels
# Arguments: 
## x: x axis positions for vertical gridlines
## y: y axis positions for horizontal gridlines
## p: x axis positions for percentile labels
## cent: percentiles to be labeled
## y_P: y axis positions for percentile labels

add_gridlines <- function(x = c(0, 5, 10, 15, 20 ,25), 
                          y = 0:10*10, 
                          P = 26, 
                          cent = cent, 
                          y_P){
  for (i in 1:length(x)){
    lines(c(x[i], x[i]), c(-100, 1000), col = "lightgrey", lwd = 0.5)
  }
  for (i in 1:length(y)){
    lines(c(-100, 1000), c(y[i], y[i]),  col = "lightgrey", lwd = 0.5)
  }
  for (i in 1:length(cent)){
    text(P, y_P[i] , paste("P", round(cent[i], 1)))
  }
}

# Description: Function to add gridlines
# Arguments: 
## x: x axis positions for vertical gridlines
## y: y axis positions for horizonal gridlines
## col: color
## lwd: line width

add_gridlines_simple <- function(x = c(0, 5, 10, 15, 20 ,25), 
                                 y = 0:10*10, 
                                 col = "lightgrey", 
                                 lwd = 0.5){
  for (i in 1:length(x)){
    lines(c(x[i], x[i]), c(-100, 1000), col = col, lwd = lwd)
  }
  for (i in 1:length(y)){
    lines(c(-100, 1000), c(y[i], y[i]),  col = col, lwd = lwd)
  }
}

# Description: Function to calculate weighted and unweighted mean and SD
# Arguments: 
## df: dataset used for calculation
## var: variable name

mean_sd_cols <- function(df, var){
  x <- df[[var]]
  w <- df$W
  mu_unw <- mean(x, na.rm=TRUE)
  sd_unw <- sd(x, na.rm=TRUE)
  mu_w   <- weighted.mean(x, w, na.rm=TRUE)
  v_w <- sum(w * (x - mu_w)^2, na.rm=TRUE) / sum(w, na.rm=TRUE)
  sd_w <- sqrt(v_w)
  c(sprintf("%.2f (%.2f)", mu_unw, sd_unw),
    sprintf("%.2f (%.2f)", mu_w, sd_w))
}

# Description: Function to calculate weighted and unweighted proportions
# Arguments: 
## df: dataset used for calculation
## flag_var: flag variable name

prop_cols <- function(df, flag_var){
  x <- df[[flag_var]]
  w <- df$W
  n_unw <- sum(as.numeric(x) > 0, na.rm = TRUE)
  p_unw <- mean(as.numeric(x) > 0, na.rm = TRUE) * 100
  p_w   <- sum(w[as.numeric(x) > 0], na.rm = TRUE) / sum(w, na.rm = TRUE) * 100
  c(sprintf("%d (%.1f%%)", n_unw, p_unw),   # Unweighted column: n (unweighted %)
    sprintf("%d (%.1f%%)", n_unw, p_w))    # Weighted  column: n (weighted %)
}

# Description: Function to calculate weighted and unweighted crp summaries
# Arguments: 
## df: dataset used for calculation

crp_cols <- function(df){
  x <- df$crp*10
  w <- df$W
  med_unw <- median(x, na.rm=TRUE)
  q_unw   <- quantile(x, c(.25,.75), na.rm=TRUE, names=FALSE)
  q_w     <- Hmisc::wtd.quantile(x, weights=w, probs=c(.5,.25,.75), na.rm=TRUE)
  c(sprintf("%.2f (%.2f-%.2f)", med_unw, q_unw[1], q_unw[2]),
    sprintf("%.2f (%.2f-%.2f)", q_w[1], q_w[2], q_w[3]))
}

# Description: Function to calculate weighted and unweighted proportions of crp below 3 mg/L
# Arguments: 
## df: dataset used for calculation

crp_below3_cols <- function(df){
  x <- df$crp * 10   # convert to mg/L
  w <- df$W
  
  n_unw <- sum(x < 3, na.rm = TRUE)
  p_unw <- mean(x < 3, na.rm = TRUE) * 100
  p_w   <- sum(w[x < 3], na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE) * 100
  
  c(
    sprintf("%d (%.1f%%)", n_unw, p_unw),
    sprintf("%d (%.1f%%)", n_unw, p_w)
  )
}

# Description: Function to create a summary column
# Arguments: 
## df: dataset used for summary

one_col <- function(df){
  out <- list()
  out[["N"]] <- c(format(nrow(df), big.mark=","), "")
  out[["Age (years)"]]               <- mean_sd_cols(df, "age")
  out[["Weight (kg)"]]               <- mean_sd_cols(df, "weight")
  out[["Height (cm)"]]               <- mean_sd_cols(df, "height")
  out[["BMI (kg/m2)"]]               <- mean_sd_cols(df, "bmi")
  out[["BMI-for-age z-score (WHO)"]] <- mean_sd_cols(df, "z_BMI_WHO")
  out[["Overweight (%)"]]            <- prop_cols(df, "overweight")
  out[["Normal weight (%)"]]         <- prop_cols(df, "normal_weight")
  out[["Underweight (%)"]]           <- prop_cols(df, "underweight")
  out[["CRP (mg/L)"]]            <- crp_cols(df)
  out[["CRP < 3 mg/L (%)"]]      <- crp_below3_cols(df)
  
  # Age group 
  for (lev in levels(df$agegroup2)) {
    n_unw <- sum(df$agegroup2 == lev, na.rm = TRUE)
    p_unw <- mean(df$agegroup2 == lev, na.rm = TRUE) * 100
    p_w   <- sum(df$W[df$agegroup2 == lev], na.rm = TRUE) / sum(df$W, na.rm = TRUE) * 100
    out[[paste0("Age group (years) - ", lev)]] <- c(
      sprintf("%d (%.1f%%)", n_unw, p_unw),  # Unweighted: n (unweighted %)
      sprintf("%d (%.1f%%)", n_unw, p_w)     # Weighted:   n (weighted %)
    )
  }
  
  # ---- Region / Subregion / Country / Study 
  out[["Region"]] <- c("", "")   
  for (lev in levels(dat0$region)) {
    n_unw <- sum(df$region == lev, na.rm = TRUE)
    p_unw <- mean(df$region == lev, na.rm = TRUE) * 100
    p_w   <- sum(df$W[df$region == lev], na.rm = TRUE) / sum(df$W, na.rm = TRUE) * 100
    out[[lev]] <- c(
      sprintf("%d (%.1f%%)", n_unw, p_unw),
      sprintf("%d (%.1f%%)", n_unw, p_w)
    )
  }
  
  # ---- Subregion 
  out[["Subregion"]] <- c("", "")
  for (lev in levels(dat0$subregion)) {
    n_unw <- sum(df$subregion == lev, na.rm = TRUE)
    p_unw <- mean(df$subregion == lev, na.rm = TRUE) * 100
    p_w   <- sum(df$W[df$subregion == lev], na.rm = TRUE) / sum(df$W, na.rm = TRUE) * 100
    out[[lev]] <- c(
      sprintf("%d (%.1f%%)", n_unw, p_unw),
      sprintf("%d (%.1f%%)", n_unw, p_w)
    )
  }

  # ---- Country 
  out[["Country"]] <- c("", "")
  for (lev in levels(dat0$country)) {
    n_unw <- sum(df$country == lev, na.rm = TRUE)
    p_unw <- mean(df$country == lev, na.rm = TRUE) * 100
    p_w   <- sum(df$W[df$country == lev], na.rm = TRUE) / sum(df$W, na.rm = TRUE) * 100
    out[[lev]] <- c(
      sprintf("%d (%.1f%%)", n_unw, p_unw),
      sprintf("%d (%.1f%%)", n_unw, p_w)
    )
  }
  
  # ---- Study
  out[["Study"]] <- c("", "")
  for (lev in levels(dat0$study_acronym)) {
    n_unw <- sum(df$study_acronym == lev, na.rm = TRUE)
    p_unw <- mean(df$study_acronym == lev, na.rm = TRUE) * 100
    p_w   <- sum(df$W[df$study_acronym == lev], na.rm = TRUE) / sum(df$W, na.rm = TRUE) * 100
    out[[lev]] <- c(
      sprintf("%d (%.1f%%)", n_unw, p_unw),
      sprintf("%d (%.1f%%)", n_unw, p_w)
    )
  }
  
  tibble(
    Characteristics = names(out),
    Unweighted = sapply(out, `[`, 1),
    Weighted   = sapply(out, `[`, 2)
  )
}


# Description: Function to find out which study is left out
# Arguments: 
## mod: fitted LOSO model
## data: full dataset
## fac: grouping variables for leave-one-out

find_leftout <- function(mod, data, fac = "study_acronym") {
  nfit <- try(nobs(mod), silent = TRUE)
  if (inherits(nfit, "try-error") || is.null(nfit)){
    nfit <- length(fitted(mod))
  } 
    
  levs <- levels(data[[fac]])
  for (s in levs) {
    n_sub <- nrow(subset(data, data[[fac]] != s))
    if (n_sub == nfit){
      return(s)
    } 
  }
  NA_character_
}

## *****************************************************************************
# 5 Leave-one-study-out (LOSO) cross validation -------------
## *****************************************************************************

# Description: Function to fit loso models
# Arguments: 
## mod: fitted model to be updated
## analysis_data: dataset used for model fitting
## ana_data_factor: grouping variables for leave-one-out 

LOSO_CV_mods_only <- function(mod,    
                              analysis_data, 
                              ana_data_factor){
  
  mod_list <- vector(mode = "list", length = length(unique(ana_data_factor)))
  
  i <- 1
  
  for (fact in unique(ana_data_factor)){
    print(fact)
    
    data_LOSO <- subset(analysis_data, ana_data_factor != fact)
    data_one_study <- subset(analysis_data, ana_data_factor == fact)
    
    update_mod <- update(mod, data = data_LOSO)
    mod_list[[i]] <- update_mod
    i <- i +1
    
  }
  
  return(mod_list)
}


## *****************************************************************************
# 6 Combining 1st and 2nd part model (hurdle model) --------------------------
## *****************************************************************************

# Description: Function to calculate the percentile rank for 2nd part of the hurdle model
# Arguments: 
## p_1: percentile rank of the 1st part model
## perc_comb: combined percentile rank

Calc_Perc_2 <- function(P_1, Perc_comb){
  Perc_2 <- (Perc_comb - P_1) / (1 - P_1)
  Perc_2 <- pmax(Perc_2, 0.0000001)
  return(Perc_2)
}

# Description: Function to add the cutoff based percentile
# Arguments: 
## cutoff: crp cutoff value
## shift: shift value (to shift GB1 distribution)
## all_perc_list: list of percentile prediction results (distribution parameters)

add_cutoffs <- function(cutoff = 0.3, shift = 0.02, all_perc_list){
  cutoff_shifted <- cutoff - shift
  
  for (i in 1:length(all_perc_list)){
    
    #Find the combined percentile rank at age 18
    P2_18 <- pGB1(cutoff_shifted, 
                  mu = all_perc_list[[i]]$mu[all_perc_list[[i]]$age==18],
                  sigma = all_perc_list[[i]]$sigma[all_perc_list[[i]]$age==18],
                  nu = all_perc_list[[i]]$nu[all_perc_list[[i]]$age==18],
                  tau = all_perc_list[[i]]$tau[all_perc_list[[i]]$age==18])
    
    P1_18 <- all_perc_list[[i]]$P_1[all_perc_list[[i]]$age==18] # percentage part 1 mod
    Perc_comb_18 <- P1_18 + ((1-P1_18)*P2_18)# required combined percentage
    
    # Transfer the percentile cut-off to the whole age range
    P2_cutoff <- Calc_Perc_2(all_perc_list[[i]]$P_1, Perc_comb_18) #calc age depending part 2 cutoff
    all_perc_list[[i]]$Perc_comb_18 <- Perc_comb_18 #save fix cut-off in table
    
    #Final cutoff as concrete quantile of the 2nd part model (shifted) 
    all_perc_list[[i]]$Perc_comb_cutoff <- qGB1(P2_cutoff,
                                                mu = all_perc_list[[i]]$mu, 
                                                sigma = all_perc_list[[i]]$sigma,
                                                nu = all_perc_list[[i]]$nu, 
                                                tau = all_perc_list[[i]]$tau) + shift
  }
  return(all_perc_list)
}


# Description: Function to combine part 1 and part 2 models for final curves
# Arguments: 
## age_seq: age consequences
## mod1: the fitted first part of the model
## mod2: the fitted second part of the model
## studies: the studies used for prediction

calc_all_perc_hurdle <- function(age_seq = age_seq, 
                                 mod1 = mod_logRA_f, 
                                 mod2 = mod_age_GB1_f_all2,
                                 studies = studies){
  
  All_perc_comb <- list()
  for (i in seq_along(studies)) {
    newdata <- data.frame(age = age_seq,
                          study_acronym = studies[i]
    )  
    # Preds calc
    P_1 <- predict(mod1, newdata = newdata, type = "response")
    P2_50 <- Calc_Perc_2(P_1, 0.5)
    P2_75 <- Calc_Perc_2(P_1, 0.75)
    P2_85 <- Calc_Perc_2(P_1, 0.85)
    P2_90 <- Calc_Perc_2(P_1, 0.9)
    P2_95 <- Calc_Perc_2(P_1, 0.95)
    P2_97 <- Calc_Perc_2(P_1, 0.97)
    
    #Age-specific parameters
    mu <- predict(mod2, newdata = newdata, type = "response", what ="mu")
    sigma <- predict(mod2, newdata = newdata, type = "response", what ="sigma")
    nu <- predict(mod2, newdata = newdata, type = "response", what ="nu")
    tau <- predict(mod2, newdata = newdata, type = "response", what ="tau")
    
    #Percentile rank of shifted distribution
    Perc_comb_50 <- qGB1(P2_50, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
    Perc_comb_75 <- qGB1(P2_75, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
    Perc_comb_85 <- qGB1(P2_85, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
    Perc_comb_90 <- qGB1(P2_90, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
    Perc_comb_95 <- qGB1(P2_95, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
    Perc_comb_97 <- qGB1(P2_97, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
    
    # Dataframe combine
    out <- cbind(newdata, 
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
    
    #Write in list
    All_perc_comb[[ as.character(studies[i]) ]] <- out
  }
  
  ## "global"
  mu_terms <- predict(mod1, newdata = newdata, type = "terms")
  P_1 <- plogis(mu_terms[, 'pb(age, max.df = 8, method = "GAIC", k = k)'] + attributes(mu_terms)$constant)
  
  P2_50 <- Calc_Perc_2(P_1, 0.5)
  P2_75 <- Calc_Perc_2(P_1, 0.75)
  P2_85 <- Calc_Perc_2(P_1, 0.85) 
  P2_90 <- Calc_Perc_2(P_1, 0.9)
  P2_95 <- Calc_Perc_2(P_1, 0.95)
  P2_97 <- Calc_Perc_2(P_1, 0.97)
  
  mu_terms_GB1 <- predict(mod2, newdata = newdata, type = "terms", what ="mu")
  mu <- plogis(mu_terms_GB1[, 'pb(age, max.df = 8, method = "GAIC", k = k2)'] + 
                 attributes(mu_terms_GB1)$constant)
  sigma <- predict(mod2, newdata = newdata, type = "response", what ="sigma")
  nu <- predict(mod2, newdata = newdata, type = "response", what ="nu")
  tau <- predict(mod2, newdata = newdata, type = "response", what ="tau")
  
  #Percentile of shifted distribution
  Perc_comb_50 <- qGB1(P2_50, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
  Perc_comb_75 <- qGB1(P2_75, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
  Perc_comb_85 <- qGB1(P2_85, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
  Perc_comb_90 <- qGB1(P2_90, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
  Perc_comb_95 <- qGB1(P2_95, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
  Perc_comb_97 <- qGB1(P2_97, mu = mu, sigma = sigma, nu = nu, tau = tau) + 0.02
  
  out <- data.frame(age = age_seq,
                    study_acronym = "",
                    P_1 = P_1,
                    P2_50 = P2_50, P2_75 = P2_75, P2_85 = P2_85, P2_90 = P2_90, P2_95 = P2_95, P2_97 = P2_97,
                    mu = mu, sigma = sigma, nu = nu, tau = tau,
                    Perc_comb_50 = Perc_comb_50,
                    Perc_comb_75 = Perc_comb_75,
                    Perc_comb_85 = Perc_comb_85,
                    Perc_comb_90 = Perc_comb_90,
                    Perc_comb_95 = Perc_comb_95,
                    Perc_comb_97 = Perc_comb_97)
  #Write in list
  All_perc_comb[[ "global" ]] <- out
  
  #calculate cutoff Percentiles
  All_perc_comb <- add_cutoffs(cutoff = 0.3, shift = 0.02, all_perc_list = All_perc_comb)
}


# Description: Function to calculate the individual's crp percentile rank
# Arguments: 
## df: dataset used for calculation
## mod1: fitted first part model
## mod2: fitted second part model

calc_individual_P <- function(df = rep_ana_data, mod1, mod2){
  #shifting crp
  df$shift_crp <- df$crp - 0.02
  
  #creating empty new data only including age
  newdata <- data.frame(age = df$age, study_acronym = "")
  
  # predict first part model parameters 
  mu_terms <- predict(mod1, newdata = newdata, type = "terms")
  P_1 <- plogis(mu_terms[, 'pb(age, max.df = 8, method = "GAIC", k = k)'] + attributes(mu_terms)$constant)
  
  # predict 2nd part model parameters
  mu_terms_GB1 <- predict(mod2, newdata = newdata, type = "terms", what ="mu")
  mu <- plogis(mu_terms_GB1[, 'pb(age, max.df = 8, method = "GAIC", k = k2)'] + 
                 attributes(mu_terms_GB1)$constant)
  sigma <- predict(mod2, newdata = newdata, type = "response", what ="sigma")
  nu <- predict(mod2, newdata = newdata, type = "response", what ="nu")
  tau <- predict(mod2, newdata = newdata, type = "response", what ="tau")
  
  # combine and calculate F_total from P1 and P2 
  df$P_crp <- ifelse(df$crp <= 0.02, 
                     P_1/2, 
                     P_1 + ((1 - P_1) * pGB1(df$shift_crp, 
                                             mu = mu,
                                             sigma = sigma,
                                             nu = nu,
                                             tau = tau))
  )
  return(df)
}
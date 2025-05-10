################################################################################
# File Name: sensitivity_results                                               #
#                                                                              #
# Purpose:   Model Influenza dynamics from symptom dynamics in Virginia,       #
#            Maryland counties, and during COVID-19.                           #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load and transform data                                        # 
#            3. Create Virginia plot                                           #
#            4. Create Maryland counties plots                                 #
#            5. Create COVID-19 plot                                           #
#                                                                              #
# Project:   Syndromic Influenza                                               #
# Author:    Ronan Corgel                                                      #
################################################################################

####################
# 1. SET-UP SCRIPT #
####################

# Start with a clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(lubridate)
library(mvnfast)
library(mvtnorm)
library(splines)
library(MCMCpack)
library(mcmcplots)

# Set the seed
set.seed(123456)

# Set the directory
setwd('/Users/rcorgel/Documents/symptomatic-flu-proj/')

##############################
# 2. LOAD AND TRANSFORM DATA #
##############################

# Load the formatted Maryland data
flu_md <- readRDS('./tmp/flu_symp_maryland.rds')

# Convert all count data to log(count)
flu_md$flu_log <- log(flu_md$flu + 1)
flu_md$fever_log <- log(flu_md$fever + 1)
flu_md$cough_log <- log(flu_md$cough + 1)
flu_md$sore_throat_log <- log(flu_md$sore_throat + 1)
flu_md$fatigue_log <- log(flu_md$fatigue + 1)
flu_md$myalgia_log <- log(flu_md$myalgia + 1)
flu_md$hypoxemia_log <- log(flu_md$hypoxemia + 1)
flu_md$short_breath_log <- log(flu_md$short_breath + 1)
flu_md$bronchitis_log <- log(flu_md$bronchitis + 1)
flu_md$chest_pain_log <- log(flu_md$chest_pain + 1)
flu_md$nausea_vom_log <- log(flu_md$nausea_vom + 1)
flu_md$headache_log <- log(flu_md$headache + 1)
flu_md$diarrhea_log <- log(flu_md$diarrhea + 1)
flu_md$congestion_log <- log(flu_md$congestion + 1)
flu_md$sneezing_log <- log(flu_md$sneezing + 1)

###########################
# 2. CREATE VIRGINIA PLOT #
###########################

# Set data for the model
Y <- flu_md$flu_log
X <- model.matrix(Y ~ fever_log + sore_throat_log + cough_log +
                    fatigue_log + myalgia_log + hypoxemia_log +
                    short_breath_log + bronchitis_log + chest_pain_log +
                    nausea_vom_log + headache_log + diarrhea_log +
                    congestion_log + sneezing_log,
                  data = flu_md) # Parametric
Theta <- bs(flu_md$feelslike, df = 3) # Non-parametric
n <- length(Y)
j <- ncol(Theta)

# Set number of iterations
B <- 30000

# Store parameter values in matrices
beta <- matrix(0, nrow = B, ncol = ncol(X))
sigma_sq <- matrix(0, nrow = B, ncol = 1)
lambda <- matrix(0, nrow = B, ncol = ncol(Theta))
tau_sq <- matrix(0, nrow = B, ncol = 1)

# Set initial values (All 1's)
beta[1, ] <- rep(1, ncol(X))
sigma_sq[1, ] <- 1
lambda[1, ] <- rep(1, ncol(Theta))
tau_sq[1, ] <- 1

# Some matrix algebra outside of the loop
XtX_inv <- solve(t(X) %*% X)
ThtTh <- t(Theta) %*% Theta

# Run the Gibbs sampler
for (b in 2:B) { 
  
  # Sample betas
  beta[b,] <- c(mvtnorm::rmvnorm(1, XtX_inv%*%t(X)%*%(Y - Theta%*%lambda[b-1,]),
                                 sigma_sq[b-1,]*XtX_inv))
  
  # Sample sigma_sq
  sigma_sq[b, ] <- rinvgamma(1, n/2, sum((Y - X%*%beta[b-1,] - 
                                            Theta%*%lambda[b-1,])^2)/2)
  
  # Sample lambdas
  sigma_lambda <- solve((ThtTh/sigma_sq[b-1, ]) + diag(j) / tau_sq[b-1,])
  mean_lambda <- (t(Theta)%*%(Y-X%*%beta[b-1,]))/sigma_sq[b-1,]
  lambda[b,] <- c(mvtnorm::rmvnorm(1,  sigma_lambda%*%mean_lambda, sigma_lambda))
  
  # Sample tau_sq
  tau_sq[b, ] <- rgamma(1, j/2, sum(lambda[b-1,]^2)/2)
  
}

# Set empty vectors for point-wise WAIC calculations
lppdv <- vector(length = n)
pwaicv <- vector(length = n)

# Loop through all weeks
for (i in 1:n) {
  Yi <- Y[i]
  mu <- beta[c(15001:30000),] %*% X[i,] + lambda[c(15001:30000),] %*% Theta[i,]
  liki <- dnorm(Yi, mu, sqrt(sigma_sq[c(15001:30000),]))
  lppdv[i] <- 1/(B/2)*sum(liki)
  pwaicv[i]	<- (1/((B/2)-1))*sum((log(liki) - mean(log(liki)))^2)
}

# Calculate WAIC
lppd	<- sum(log(lppdv))
pwaic	<- sum(pwaicv)
WAIC_all_temp	<- -2*lppd + 2*pwaic

beta_pred <- beta
lambda_pred <- lambda

# Load the formatted Maryland data
flu_va <- readRDS('./tmp/flu_symp_virginia.rds')

# Convert all count data to log(count)
flu_va$flu_log <- log(flu_va$flu + 1)
flu_va$fever_log <- log(flu_va$fever + 1)
flu_va$cough_log <- log(flu_va$cough + 1)
flu_va$sore_throat_log <- log(flu_va$sore_throat + 1)
flu_va$fatigue_log <- log(flu_va$fatigue + 1)
flu_va$myalgia_log <- log(flu_va$myalgia + 1)
flu_va$hypoxemia_log <- log(flu_va$hypoxemia + 1)
flu_va$short_breath_log <- log(flu_va$short_breath + 1)
flu_va$bronchitis_log <- log(flu_va$bronchitis + 1)
flu_va$chest_pain_log <- log(flu_va$chest_pain + 1)
flu_va$nausea_vom_log <- log(flu_va$nausea_vom + 1)
flu_va$headache_log <- log(flu_va$headache + 1)
flu_va$diarrhea_log <- log(flu_va$diarrhea + 1)
flu_va$congestion_log <- log(flu_va$congestion + 1)
flu_va$sneezing_log <- log(flu_va$sneezing + 1)


pred_values <- NULL

Y <- flu_va$flu_log
X <- model.matrix(Y ~ fever_log + sore_throat_log + cough_log +
                    fatigue_log + myalgia_log + hypoxemia_log +
                    short_breath_log + bronchitis_log + chest_pain_log +
                    nausea_vom_log + headache_log + diarrhea_log +
                    congestion_log + sneezing_log,
                  data = flu_va) # Parametric
Theta <- bs(flu_va$feelslike, df = 3) # Non-parametric
n <- length(Y)
for (i in 1:n) {
  Yi <- Y[i]
  mu <- beta[c(15001:30000),] %*% X[i,] + lambda[c(15001:30000),] %*% Theta[i,]
  est <- t(apply(mu, 2, quantile, probs = c(0.5, 0.05, 0.95)))
  pred_values <- rbind(pred_values, est)
}

flu_va <- cbind(flu_va, as.data.frame(pred_values))
flu_va$median <- exp(flu_va$`50%`) - 1
flu_va$low <- exp(flu_va$`5%`) - 1
flu_va$high <- exp(flu_va$`95%`) - 1


va_plot <- ggplot(data = flu_va[2:(nrow(flu_va) - 1),]) + 
  geom_line(aes(x = week_date, y = flu), 
            color = 'black', linewidth = 1.2, alpha = 0.8, lty = 2) + 
  geom_line(aes(x = week_date, y = median), 
            color = '#C7A939', linewidth = 1.2, alpha = 0.7) +
  geom_ribbon(aes(x = week_date, ymin = low, ymax = high), fill = "#C7A939", alpha = 0.35) +
  ylab('Flu Count') + xlab('Week') + 
  ggtitle('Predicted vs. Observed Virginia Flu Dynamics') +
  theme_minimal() + 
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

flu_va$test <- ifelse(flu_va$flu < flu_va$high & flu_va$flu > flu_va$low, 1, 0)
sum(flu_va$test) / 210


library(splus2R)

flu_va$obs_peak <- peaks(flu_va$flu, span=21)
obs_peak <- flu_va |> filter(obs_peak == TRUE) |>
  dplyr::select(c(week_date, obs_peak))
flu_va$pred_peak <- peaks(flu_va$median, span=21)
pred_peak <- flu_va |> filter(pred_peak == TRUE) |>
  dplyr::select(c(week_date, pred_peak)) |>
  rename('week_date_2' = 'week_date')
peaks <- cbind(obs_peak, pred_peak)
peaks$diff_vals <- peaks$week_date - peaks$week_date_2

#####################################
# 3. CREATE MARYLAND COUNTIES PLOTS #
#####################################

# Load the formatted Maryland data
flu_md <- readRDS('./tmp/flu_symp_maryland_counties.rds')

# Convert all count data to log(count)
flu_md$flu_log <- log(flu_md$flu + 1)
flu_md$fever_log <- log(flu_md$fever + 1)
flu_md$cough_log <- log(flu_md$cough + 1)
flu_md$sore_throat_log <- log(flu_md$sore_throat + 1)
flu_md$fatigue_log <- log(flu_md$fatigue + 1)
flu_md$myalgia_log <- log(flu_md$myalgia + 1)
flu_md$hypoxemia_log <- log(flu_md$hypoxemia + 1)
flu_md$short_breath_log <- log(flu_md$short_breath + 1)
flu_md$bronchitis_log <- log(flu_md$bronchitis + 1)
flu_md$chest_pain_log <- log(flu_md$chest_pain + 1)
flu_md$nausea_vom_log <- log(flu_md$nausea_vom + 1)
flu_md$headache_log <- log(flu_md$headache + 1)
flu_md$diarrhea_log <- log(flu_md$diarrhea + 1)
flu_md$congestion_log <- log(flu_md$congestion + 1)
flu_md$sneezing_log <- log(flu_md$sneezing + 1)


county_list <- unique(flu_md$county_fips)

percent_correct <- NULL
peak_difference <- NULL

flu_md$year <- year(flu_md$week_date)
flu_md$season <- ifelse(flu_md$year == 2016 & flu_md$month > 8, '2016 - 2017', '')
flu_md$season <- ifelse(flu_md$year == 2017 & flu_md$month < 9, '2016 - 2017', flu_md$season)
flu_md$season <- ifelse(flu_md$year == 2017 & flu_md$month > 8, '2017 - 2018', flu_md$season)
flu_md$season <- ifelse(flu_md$year == 2018 & flu_md$month < 9, '2017 - 2018', flu_md$season)
flu_md$season <- ifelse(flu_md$year == 2018 & flu_md$month > 8, '2018 - 2019', flu_md$season)
flu_md$season <- ifelse(flu_md$year == 2019 & flu_md$month < 9, '2018 - 2019', flu_md$season)
flu_md$season <- ifelse(flu_md$year == 2019 & flu_md$month > 8, '2019 - 2020', flu_md$season)
flu_md$season <- ifelse(flu_md$year == 2020 & flu_md$month < 9, '2019 - 2020', flu_md$season)

for (i in county_list) {
  county_dat <- flu_md |> filter(county_fips == i)
  Y <- county_dat$flu_log
  X <- model.matrix(Y ~ fever_log + sore_throat_log + cough_log +
                      fatigue_log + myalgia_log + hypoxemia_log +
                      short_breath_log + bronchitis_log + chest_pain_log +
                      nausea_vom_log + headache_log + diarrhea_log +
                      congestion_log + sneezing_log,
                    data = county_dat) # Parametric
  Theta <- bs(county_dat$feelslike, df = 3) # Non-parametric
  n <- length(Y)
  pred_values <- NULL
  for (i in 1:n) {
    Yi <- Y[i]
    mu <- beta[c(15001:30000),] %*% X[i,] + lambda[c(15001:30000),] %*% Theta[i,]
    est <- t(apply(mu, 2, quantile, probs = c(0.5, 0.05, 0.95)))
    pred_values <- rbind(pred_values, est)
  }
  
  county_dat <- cbind(county_dat, as.data.frame(pred_values))
  county_dat$median <- exp(county_dat$`50%`) - 1
  county_dat$low <- exp(county_dat$`5%`) - 1
  county_dat$high <- exp(county_dat$`95%`) - 1
  
  county_dat$test <- ifelse(county_dat$flu < county_dat$high & county_dat$flu > county_dat$low, 1, 0)
  percent_correct <- c(percent_correct, (sum(county_dat$test) / 210))
  
  county_dat$obs_peak <- peaks(county_dat$flu, span=41)
  obs_peak <- county_dat |> filter(obs_peak == TRUE) |>
    dplyr::select(c(week_date, obs_peak, season))
  county_dat$pred_peak <- peaks(county_dat$median, span=41)
  pred_peak <- county_dat |> filter(pred_peak == TRUE) |>
    dplyr::select(c(week_date, pred_peak, season)) |>
    rename('week_date_2' = 'week_date') |>
    rename('county_fips_2' = 'county_fips')
  peaks <- left_join(obs_peak, pred_peak, by = c('season' = 'season'))
  peaks$diff_vals <- peaks$week_date_2 - peaks$week_date
  peaks <- peaks |> dplyr::select(county_fips, week_date, week_date_2, diff_vals)
  peak_difference <- rbind(peak_difference, peaks)
  
}

flu_md_county <- flu_md
flu_md <- flu_md |> filter(county_fips == '24027')

# Set data for the model
Y <- flu_md$flu_log
X <- model.matrix(Y ~ fever_log + sore_throat_log + cough_log +
                    fatigue_log + myalgia_log + hypoxemia_log +
                    short_breath_log + bronchitis_log + chest_pain_log +
                    nausea_vom_log + headache_log + diarrhea_log +
                    congestion_log + sneezing_log,
                  data = flu_md) # Parametric
Theta <- bs(flu_md$feelslike, df = 3) # Non-parametric
n <- length(Y)
j <- ncol(Theta)

# Set number of iterations
B <- 30000

# Store parameter values in matrices
beta <- matrix(0, nrow = B, ncol = ncol(X))
sigma_sq <- matrix(0, nrow = B, ncol = 1)
lambda <- matrix(0, nrow = B, ncol = ncol(Theta))
tau_sq <- matrix(0, nrow = B, ncol = 1)

# Set initial values (All 1's)
beta[1, ] <- rep(1, ncol(X))
sigma_sq[1, ] <- 1
lambda[1, ] <- rep(1, ncol(Theta))
tau_sq[1, ] <- 1

# Some matrix algebra outside of the loop
XtX_inv <- solve(t(X) %*% X)
ThtTh <- t(Theta) %*% Theta

# Run the Gibbs sampler
for (b in 2:B) { 
  
  # Sample betas
  beta[b,] <- c(mvtnorm::rmvnorm(1, XtX_inv%*%t(X)%*%(Y - Theta%*%lambda[b-1,]),
                                 sigma_sq[b-1,]*XtX_inv))
  
  # Sample sigma_sq
  sigma_sq[b, ] <- rinvgamma(1, n/2, sum((Y - X%*%beta[b-1,] - 
                                            Theta%*%lambda[b-1,])^2)/2)
  
  # Sample lambdas
  sigma_lambda <- solve((ThtTh/sigma_sq[b-1, ]) + diag(j) / tau_sq[b-1,])
  mean_lambda <- (t(Theta)%*%(Y-X%*%beta[b-1,]))/sigma_sq[b-1,]
  lambda[b,] <- c(mvtnorm::rmvnorm(1,  sigma_lambda%*%mean_lambda, sigma_lambda))
  
  # Sample tau_sq
  tau_sq[b, ] <- rgamma(1, j/2, sum(lambda[b-1,]^2)/2)
  
}

# Set empty vectors for point-wise WAIC calculations
lppdv <- vector(length = n)
pwaicv <- vector(length = n)

# Loop through all weeks
for (i in 1:n) {
  Yi <- Y[i]
  mu <- beta[c(15001:30000),] %*% X[i,] + lambda[c(15001:30000),] %*% Theta[i,]
  liki <- dnorm(Yi, mu, sqrt(sigma_sq[c(15001:30000),]))
  lppdv[i] <- 1/(B/2)*sum(liki)
  pwaicv[i]	<- (1/((B/2)-1))*sum((log(liki) - mean(log(liki)))^2)
}

# Calculate WAIC
lppd	<- sum(log(lppdv))
pwaic	<- sum(pwaicv)
WAIC_county	<- -2*lppd + 2*pwaic

county_list <- unique(flu_md_county$county_fips)

percent_correct_c <- NULL
peak_difference_c <- NULL

flu_md$year <- year(flu_md$week_date)
flu_md$season <- ifelse(flu_md$year == 2016 & flu_md$month > 8, '2016 - 2017', '')
flu_md$season <- ifelse(flu_md$year == 2017 & flu_md$month < 9, '2016 - 2017', flu_md$season)
flu_md$season <- ifelse(flu_md$year == 2017 & flu_md$month > 8, '2017 - 2018', flu_md$season)
flu_md$season <- ifelse(flu_md$year == 2018 & flu_md$month < 9, '2017 - 2018', flu_md$season)
flu_md$season <- ifelse(flu_md$year == 2018 & flu_md$month > 8, '2018 - 2019', flu_md$season)
flu_md$season <- ifelse(flu_md$year == 2019 & flu_md$month < 9, '2018 - 2019', flu_md$season)
flu_md$season <- ifelse(flu_md$year == 2019 & flu_md$month > 8, '2019 - 2020', flu_md$season)
flu_md$season <- ifelse(flu_md$year == 2020 & flu_md$month < 9, '2019 - 2020', flu_md$season)

for (i in county_list) {
  county_dat <- flu_md_county |> filter(county_fips == i)
  Y <- county_dat$flu_log
  X <- model.matrix(Y ~ fever_log + sore_throat_log + cough_log +
                      fatigue_log + myalgia_log + hypoxemia_log +
                      short_breath_log + bronchitis_log + chest_pain_log +
                      nausea_vom_log + headache_log + diarrhea_log +
                      congestion_log + sneezing_log,
                    data = county_dat) # Parametric
  Theta <- bs(county_dat$feelslike, df = 3) # Non-parametric
  n <- length(Y)
  pred_values <- NULL
  for (i in 1:n) {
    Yi <- Y[i]
    mu <- beta[c(15001:30000),] %*% X[i,] + lambda[c(15001:30000),] %*% Theta[i,]
    est <- t(apply(mu, 2, quantile, probs = c(0.5, 0.05, 0.95)))
    pred_values <- rbind(pred_values, est)
  }
  
  county_dat <- cbind(county_dat, as.data.frame(pred_values))
  county_dat$median <- exp(county_dat$`50%`) - 1
  county_dat$low <- exp(county_dat$`5%`) - 1
  county_dat$high <- exp(county_dat$`95%`) - 1
  
  county_dat$test <- ifelse(county_dat$flu < county_dat$high & county_dat$flu > county_dat$low, 1, 0)
  percent_correct_c <- c(percent_correct_c, (sum(county_dat$test) / 210))
  
  county_dat$obs_peak <- peaks(county_dat$flu, span=41)
  obs_peak <- county_dat |> filter(obs_peak == TRUE) |>
    dplyr::select(c(week_date, obs_peak, season))
  county_dat$pred_peak <- peaks(county_dat$median, span=41)
  pred_peak <- county_dat |> filter(pred_peak == TRUE) |>
    dplyr::select(c(week_date, pred_peak, season)) |>
    rename('week_date_2' = 'week_date') |>
    rename('county_fips_2' = 'county_fips')
  peaks <- left_join(obs_peak, pred_peak, by = c('season' = 'season'))
  peaks$diff_vals <- peaks$week_date_2 - peaks$week_date
  peaks <- peaks |> dplyr::select(county_fips, week_date, week_date_2, diff_vals)
  peak_difference_c <- rbind(peak_difference_c, peaks)
  
}

peak_difference$Cat <- 'State-Trained'
peak_difference_c$Cat <- 'County-Trained'
peak_difference_chart <- rbind(peak_difference, peak_difference_c)

percent_correct_dat <- data.frame('State-Trained', percent_correct)
percent_correct_dat <- percent_correct_dat |>
  rename('Cat' = 'X.State.Trained.')
percent_correct_c_dat <- data.frame('County-Trained', percent_correct_c)
percent_correct_c_dat <- percent_correct_c_dat |>
  rename('Cat' = 'X.County.Trained.') |>
  rename('percent_correct' = 'percent_correct_c')
perc_correct <- rbind(percent_correct_dat, percent_correct_c_dat)

density <- ggplot(perc_correct , aes(percent_correct, fill = Cat, colour = Cat)) +
  geom_density(alpha = 0.4) + scale_fill_manual("Model", values=c('#c36272', '#347dab')) + 
  scale_color_manual('Model', values=c('#c36272', '#347dab')) +
  theme_minimal() + xlab('Percent Correct') + ylab('Density') +
  ggtitle('County-Level Percent Correct') +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

colors <- c('0-4' = '#c36272', '5-12' = '#C7A939', '13-17' = '#5a9374', 
            '18-49' = '#347dab', '50-64' = '#9086ba', '65+' = 'darkgrey')

box <- ggplot(peak_difference_chart, aes(x=Cat, y=diff_vals)) + 
  geom_boxplot(data = peak_difference_chart, aes(color = Cat, fill = Cat), alpha = 0.4) +
  scale_fill_manual("Model", values=c('#c36272', '#347dab')) + 
  scale_color_manual('Model', values=c('#c36272', '#347dab')) +
  theme_minimal() + xlab('Model') + ylab('Difference in Days') +
  ggtitle('County-Level Peak Time Difference') +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

###########################
# 4. CREATE COVID-19 PLOT #
###########################

# Load the formatted Maryland data
flu_md_2020 <- readRDS('./tmp/flu_symp_maryland_2020.rds')

# Convert all count data to log(count)
flu_md_2020$flu_log <- log(flu_md_2020$flu + 1)
flu_md_2020$fever_log <- log(flu_md_2020$fever + 1)
flu_md_2020$cough_log <- log(flu_md_2020$cough + 1)
flu_md_2020$sore_throat_log <- log(flu_md_2020$sore_throat + 1)
flu_md_2020$fatigue_log <- log(flu_md_2020$fatigue + 1)
flu_md_2020$myalgia_log <- log(flu_md_2020$myalgia + 1)
flu_md_2020$hypoxemia_log <- log(flu_md_2020$hypoxemia + 1)
flu_md_2020$short_breath_log <- log(flu_md_2020$short_breath + 1)
flu_md_2020$bronchitis_log <- log(flu_md_2020$bronchitis + 1)
flu_md_2020$chest_pain_log <- log(flu_md_2020$chest_pain + 1)
flu_md_2020$nausea_vom_log <- log(flu_md_2020$nausea_vom + 1)
flu_md_2020$headache_log <- log(flu_md_2020$headache + 1)
flu_md_2020$diarrhea_log <- log(flu_md_2020$diarrhea + 1)
flu_md_2020$congestion_log <- log(flu_md_2020$congestion + 1)
flu_md_2020$sneezing_log <- log(flu_md_2020$sneezing + 1)

pred_values <- NULL

Y <- flu_md_2020$flu_log
X <- model.matrix(Y ~ fever_log + sore_throat_log + cough_log +
                    fatigue_log + myalgia_log + hypoxemia_log +
                    short_breath_log + bronchitis_log + chest_pain_log +
                    nausea_vom_log + headache_log + diarrhea_log +
                    congestion_log + sneezing_log,
                  data = flu_md_2020) # Parametric
Theta <- bs(flu_md_2020$feelslike, df = 3) # Non-parametric
n <- length(Y)
for (i in 1:n) {
  Yi <- Y[i]
  mu <- beta_pred[c(15001:30000),] %*% X[i,] + lambda_pred[c(15001:30000),] %*% Theta[i,]
  est <- t(apply(mu, 2, quantile, probs = c(0.5, 0.05, 0.95)))
  pred_values <- rbind(pred_values, est)
}

flu_md_2020 <- cbind(flu_md_2020, as.data.frame(pred_values))
flu_md_2020$median <- exp(flu_md_2020$`50%`) - 1
flu_md_2020$low <- exp(flu_md_2020$`5%`) - 1
flu_md_2020$high <- exp(flu_md_2020$`95%`) - 1


covid <- ggplot(data = flu_md_2020[2:(nrow(flu_md_2020) - 1),]) + 
  geom_line(aes(x = week_date, y = flu), 
            color = 'black', linewidth = 1.2, alpha = 0.8, lty = 2) + 
  geom_line(aes(x = week_date, y = median), 
            color = '#5a9374', linewidth = 1.2, alpha = 0.7) +
  geom_ribbon(aes(x = week_date, ymin = low, ymax = high), fill = "#5a9374", alpha = 0.35) +
  ylab('Flu Count') + xlab('Week') + 
  ggtitle('Predicted vs. Observed Maryland Flu Dynamics') +
  theme_minimal() + 
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

flu_md_2020$test <- ifelse(flu_md_2020$flu < flu_md_2020$high & flu_md_2020$flu > flu_md_2020$low, 1, 0)
sum(flu_md_2020$test) / 210

flu_md_2020$obs_peak <- peaks(flu_md_2020$flu, span=41)
obs_peak <- flu_md_2020 |> filter(obs_peak == TRUE) |>
  dplyr::select(c(week_date, obs_peak))
flu_md_2020$pred_peak <- peaks(flu_md_2020$median, span=41)
pred_peak <- flu_md_2020 |> filter(pred_peak == TRUE) |>
  dplyr::select(c(week_date, pred_peak)) |>
  rename('week_date_2' = 'week_date')
peaks <- cbind(obs_peak, pred_peak)
peaks$diff_vals <- peaks$week_date - peaks$week_date_2

## COMBINE FIGURES FOR MANUSCRIPT ##
sub <- cowplot::plot_grid(density, box, nrow = 1,
                          labels = c('B', 'C'),
                          label_size = 26, hjust = 0)

figure_3 <- cowplot::plot_grid(va_plot, sub, covid,
                               nrow = 3,
                               labels = c('A', '', 'D'),
                               label_size = 26, hjust = 0) 

ggsave('./figs/figure_3.jpg', plot = figure_3, height = 12, width = 11)

################################################################################
################################################################################

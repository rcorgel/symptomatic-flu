################################################################################
# File Name: semi_parametric_model                                             #
#                                                                              #
# Purpose:   Model Influenza dynamics from symptom dynamics in Maryland.       #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load and transform data                                        # 
#            3. Model with ILI and dew point                                   #
#            4. Model with ILI and temperature                                 #
#            5. Model with all symptoms and dew point                          #
#            6. Model with all symptoms and temperature                        #
#            7. Perform model checks                                           #
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

##############################
# 3. ILI AND DEW POINT MODEL #
##############################

# Set data for the model
Y <- flu_md$flu_log
X <- model.matrix(Y ~ fever_log + sore_throat_log + cough_log,
                  data = flu_md) # Parametric
Theta <- bs(flu_md$dew, df = 3) # Non-parametric
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
WAIC_ili_dew	<- -2*lppd + 2*pwaic

################################
# 4. ILI AND TEMPERATURE MODEL #
################################

# Set data for the model
Y <- flu_md$flu_log
X <- model.matrix(Y ~ fever_log + sore_throat_log + cough_log,
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
WAIC_ili_temp	<- -2*lppd + 2*pwaic

##############################
# 5. ALL AND DEW POINT MODEL #
##############################

# Set data for the model
Y <- flu_md$flu_log
X <- model.matrix(Y ~ fever_log + sore_throat_log + cough_log +
                    fatigue_log + myalgia_log + hypoxemia_log +
                    short_breath_log + bronchitis_log + chest_pain_log +
                    nausea_vom_log + headache_log + diarrhea_log +
                    congestion_log + sneezing_log,
                  data = flu_md) # Parametric
Theta <- bs(flu_md$dew, df = 3) # Non-parametric
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
WAIC_all_dew	<- -2*lppd + 2*pwaic

################################
# 6. ALL AND TEMPERATURE MODEL #
################################

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

################################
# 7. CHECK THE DF OF THE MODEL #
################################

# Set empty vector to fill with WAIC values
WAIC <- data.frame(rep(NA, 10))

# Loop through degrees of freedom values
for (i in 3:10) {
  
  # Set data for the model
  Y <- flu_md$flu_log
  X <- model.matrix(Y ~ fever_log + sore_throat_log + cough_log +
                      fatigue_log + myalgia_log + hypoxemia_log +
                      short_breath_log + bronchitis_log + chest_pain_log +
                      nausea_vom_log + headache_log + diarrhea_log +
                      congestion_log + sneezing_log,
                    data = flu_md) # Parametric
  Theta <- bs(flu_md$feelslike, df = i) # Non-parametric
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
    lambda[b,] <- c(rmvn(1,  sigma_lambda%*%mean_lambda, sigma_lambda))
    
    # Sample tau_sq
    tau_sq[b, ] <- rgamma(1, j/2, sum(lambda[b-1,]^2)/2)
    
  }
  
  # Set empty vectors for point-wise WAIC calculations
  lppdv <- vector(length = n)
  pwaicv <- vector(length = n)
  
  # Loop through all weeks
  for (j in 1:n) {
    Yi <- Y[j]
    mu <- beta[c(15001:30000),] %*% X[j,] + lambda[c(15001:30000),] %*% Theta[i,]
    liki <- dnorm(Yi, mu, sqrt(sigma_sq[c(15001:30000),]))
    lppdv[j] <- 1/(B/2)*sum(liki)
    pwaicv[j]	<- (1/((B/2)-1))*sum((log(liki) - mean(log(liki)))^2)
  }
  
  # Calculate WAIC
  lppd	<- sum(log(lppdv))
  pwaic	<- sum(pwaicv)
  WAIC_val	<- -2*lppd + 2*pwaic
  WAIC[i, 1] <- WAIC_val
  
}

# Create a visual for the change in degree
WAIC[3, 1] <- WAIC_all_temp
WAIC <- WAIC |>
  mutate(Degree = row_number()) |>
  filter(!is.na(rep.NA..10.)) |>
  rename('WAIC' = 'rep.NA..10.') 

ggplot() + 
  geom_line(data = WAIC, 
            aes(x = Degree, y = WAIC), 
            color = '#9086ba', linewidth = 1.2, alpha = 0.8) + 
  ylab('WAIC') + xlab('Degrees of Freedom') + 
  ggtitle('Best Fit Model WAIC by Degrees of Freedom') +
  theme_minimal() + 
  scale_x_continuous(breaks = c(3, 4, 5, 6, 7, 8, 9, 10))

################################################################################
################################################################################

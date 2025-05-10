################################################################################
# File Name: best_fit_results                                                  #
#                                                                              #
# Purpose:   Model Influenza dynamics from symptom dynamics in Maryland.       #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load and transform data                                        # 
#            3. Create regression forest plot                                  #
#            4. Create Geweke diagnostic plot                                  #
#            5. Create Maryland estimate plot                                  #
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

#####################################
# 2. CREATE REGRESSION FORREST PLOT #
#####################################

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

coefs <- apply(beta[-c(1:(B/2)),], 2, quantile, probs = c(0.5, 0.05, 0.95))
colnames(coefs) <- colnames(X)
coef_dat <- data.frame(round(t(coefs), 4))
coef_dat$Variable <- c('Intercept', 'Fever', 'Sore Throat', 'Cough', 'Fatigue', 
                       'Myalgia', 'Hypoxemia', 'Short Breath', 'Bronchitis', 
                       'Chest Pain', 'Nausea', 'Headache', 'Diarrhea', 
                       'Congestion', 'Sneezing')
# Create forest plot
forrest <- ggplot(data = coef_dat[2:15,]) +
  geom_pointrange(aes(x=reorder(Variable, -X50.), y=X50., ymin=X5., ymax=X95.), 
                  size = 0.75, alpha = 1.0, color = '#5a9374') + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ggtitle('Coefficient Estimates') +
  xlab("") + ylab("Estimate (90% CI)") +
  theme_minimal() +  # use a white background
  theme(legend.text = element_text(size = 16),
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
        plot.title = element_text(size=20, hjust = 0))

####################################
# 3. CREATE GEWEKE DIAGNOSTIC PLOT #
####################################

# Create bar plot
geweke <- apply(beta[-c(1:(B/2)),c(2:15)], 2, geweke.diag)
gewekw_dat <- data.frame(sapply(geweke, "[[", "z"), 
                         c('Fever', 'Sore Throat', 'Cough', 'Fatigue', 
                           'Myalgia', 'Hypoxemia', 'Short Breath', 'Bronchitis', 
                           'Chest Pain', 'Nausea', 'Headache', 'Diarrhea', 
                           'Congestion', 'Sneezing'))
colnames(gewekw_dat) <- c('Geweke Value', 'Symptoms')
gewekw_dat$order <- coef_dat[2:15,]$X50.


geweke_plot <- ggplot(gewekw_dat, aes(x = reorder(Symptoms, order), y = `Geweke Value`)) +
  geom_bar(stat = "identity", fill = '#5a9374', width = 0.5,
           show.legend = FALSE) +
  ggtitle('Geweke Convergence Diagnostic') +
  xlab("") +
  geom_hline(yintercept=1.96, lty=2) +
  geom_hline(yintercept=-1.96, lty=2) +
  ylab("Geweke Value") +
  theme_minimal() +  # use a white background
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size=20),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0))

####################################
# 3. CREATE MARYLAND ESTIMATE PLOT #
####################################

pred_values <- NULL

for (i in 1:n) {
  Yi <- Y[i]
  mu <- beta[c(15001:30000),] %*% X[i,] + lambda[c(15001:30000),] %*% Theta[i,]
  est <- t(apply(mu, 2, quantile, probs = c(0.5, 0.05, 0.95)))
  pred_values <- rbind(pred_values, est)
}

flu_md <- cbind(flu_md, as.data.frame(pred_values))
flu_md$median <- exp(flu_md$`50%`) - 1
flu_md$low <- exp(flu_md$`5%`) - 1
flu_md$high <- exp(flu_md$`95%`) - 1

estimate <- ggplot(data = flu_md[2:(nrow(flu_md) - 1),]) + 
  geom_line(aes(x = week_date, y = flu), 
            color = 'black', linewidth = 1.2, alpha = 0.8, lty = 2) + 
  geom_line(aes(x = week_date, y = median), 
            color = '#5a9374', linewidth = 1.2, alpha = 0.7) +
  geom_ribbon(aes(x = week_date, ymin = low, ymax = high), fill = "#5a9374", alpha = 0.35) +
  ylab('Flu Count') + xlab('Week') + 
  ggtitle('Predicted vs. Observed Maryland Flu Dynamics') +
  theme_minimal() + 
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
        plot.title = element_text(size=20, hjust = 0))

flu_md$test <- ifelse(flu_md$flu < flu_md$high & flu_md$flu > flu_md$low, 1, 0)
sum(flu_md$test) / 210

library(splus2R)
flu_md$obs_peak <- peaks(flu_md$flu, span=41)
flu_md$pred_peak <- peaks(flu_md$median, span=41)

## COMBINE FIGURES FOR MANUSCRIPT ##
sub <- cowplot::plot_grid(forrest, geweke_plot, nrow = 1,
                          labels = c('A', 'B'),
                          label_size = 26, hjust = 0)

figure_2 <- cowplot::plot_grid(sub, estimate,
                               nrow = 2,
                               labels = c('', 'C'),
                               label_size = 26, hjust = 0) 

ggsave('./figs/figure_2.jpg', plot = figure_2, height = 12, width = 11)

################################################################################
################################################################################

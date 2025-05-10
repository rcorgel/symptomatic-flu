################################################################################
# File Name: load_symptom_data                                                 #
#                                                                              #
# Purpose:   Load and process Influenza medical claims data from 2016-2020.    #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load, append, check, and save data                             #
#            3. Impute and clean data                                          #
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
library(assertr)
library(ISOweek)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Documents/symptomatic-flu-proj/')

################
# 2. LOAD DATA #
################

# Influenza symptom data
flu_16_17 <- read.csv('raw/flu_2016-09_2017-08_p2.csv')
flu_17_18 <- read.csv('raw/flu_2017-09_2018-08_p2.csv')
flu_18_19 <- read.csv('raw/flu_2018-09_2019-08_p2.csv')
flu_19_20 <- read.csv('raw/flu_2019-09_2020-08_p2.csv')
# Don't include 2020-2023 due to low flu circulation & covid-19 co-circulation

# Append
flu <- rbind(flu_16_17, flu_17_18, flu_18_19, flu_19_20)
remove(flu_16_17, flu_17_18, flu_18_19, flu_19_20)

# Run sanity checks on the data
flu |>
  # There should be 48 months
  verify(length(unique(year_week)) == 210) |>
  # There should be more than 2500 county groups
  verify(length(unique(county_fips)) > 2500) |>
  # There are unknown age groups
  verify(length(unique(age_grp)) == 7) |>
  # There are unknown genders
  verify(length(unique(patient_gender_code)) == 3) |>
  # Check disease, symptom, and count vars are not missing
  assert(not_na, flu:patient_count)

# Save full data
saveRDS(flu, './tmp/flu_dat_raw.rds')

# Impute data
# Check the percent of rows that will be imputed (80.4%)
nrow(flu[flu$patient_count == "<=5",]) / nrow(flu)
# Impute rows where patient count is <=5 to a random number 1 to 5
flu$patient_count_imp <- flu$patient_count
flu$patient_count_imp[flu$patient_count_imp == '<=5'] <- 
  sample(1:5, length(flu$patient_count_imp[flu$patient_count_imp == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
flu$patient_count_imp <- as.numeric(flu$patient_count_imp)
# Check the percent of patients that were imputed (7.9%)
sum(flu[flu$patient_count_imp <= 5,]$patient_count_imp) / sum(flu$patient_count_imp)

# Create new variables
flu <- flu |>
  # Convert week date to date
  mutate(week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", year_week)) |>
  mutate(week_date = ISOweek2date(week_date)) |>
  # Create state fips
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  # Create symptom count
  mutate(symp_count = fever + myalgia + cough + 
           sore_throat + short_breath + hypoxemia + 
           chest_pain + bronchitis + nausea_vom + 
           diarrhea + fatigue + headache + congestion + sneezing)

# Check the percent of patients that will be dropped for:
# Missing county or not US state & DC (0.3%)
sum(flu[flu$state_fips %in% c("72","99",""),]$patient_count_imp) / sum(flu$patient_count_imp)
# Unknown or missing ages (0.8%)
sum(flu[flu$age_grp %in% c("U",""),]$patient_count_imp) / sum(flu$patient_count_imp)
# Unknown or missing genders (<0.1%)
sum(flu[flu$patient_gender_code %in% c("U",""),]$patient_count_imp) / sum(flu$patient_count_imp)
# Asymptomatic or non-coded individuals (84.2%)
sum(flu[flu$symp_count == 0,]$patient_count_imp) / sum(flu$patient_count_imp)

# Drop based on:
flu_filt <- flu |>
  # Missing county or not US state & DC
  filter(!(state_fips %in% c("72","99",""))) |>
  # Unknown or missing ages
  filter(!(age_grp %in% c("U",""))) |>
  # Unknown or missing genders
  filter(!(patient_gender_code %in% c("U","")))
# Check the number of individuals kept (98.2%)
sum(flu_filt$patient_count_imp) / sum(flu$patient_count_imp)

# Confirm drop was completed
flu_filt |>
  # There should be more than 2500 county groups
  verify(state_fips != "99" | state_fips != "72") |>
  # There are unknown age groups
  verify(age_grp != "U") |>
  # There are unknown genders
  verify(patient_gender_code != "U") |>
  # Check disease, symptom, and count vars are not missing
  assert(not_na, county_fips:symp_count)

# Save filtered data
saveRDS(flu_filt, './tmp/flu_dat_proc.rds')

################################################################################
################################################################################

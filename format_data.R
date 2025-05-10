################################################################################
# File Name: format_data                                                       #
#                                                                              #
# Purpose:   Format Influenza medical claims data from 2016-2020.              #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Limit data to symptomatic                                      #
#            3. Collapse data to county-level                                  #
#            4. Merge on weather data                                          #
#            5. Merge on all cause data                                        #
#            6. Restrict data to specific counties                             #      
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
library(ISOweek)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Documents/symptomatic-flu-proj/')

################################
# 2. LIMIT DATA TO SYMPTOMATIC #
################################

# Load data
flu <- readRDS('./tmp/flu_dat_proc.rds')

# Keep symptomatic individuals
flu_symp <- flu |>
  filter(symp_count > 0)
# Check the number of individuals kept (15.7%)
sum(flu_symp$patient_count_imp) / sum(flu$patient_count_imp)
remove(flu)

# Save symptomatic version
saveRDS(flu_symp, './tmp/flu_dat_proc_symp.rds')

####################################
# 3. COLLAPSE DATA TO COUNTY-LEVEL #
####################################

# Bring to county-week-age level
flu_symp_county_age <- flu_symp |> group_by(county_fips, week_date, age_grp, 
                                        state_fips, month) |>
  mutate(flu = sum(flu*patient_count_imp),
         fever = sum(fever*patient_count_imp),
         cough = sum(cough*patient_count_imp),
         sore_throat = sum(sore_throat*patient_count_imp),
         fatigue = sum(fatigue*patient_count_imp),
         myalgia = sum(myalgia*patient_count_imp),
         hypoxemia = sum(hypoxemia*patient_count_imp),
         short_breath = sum(short_breath*patient_count_imp),
         bronchitis = sum(bronchitis*patient_count_imp),
         chest_pain = sum(chest_pain*patient_count_imp),
         nausea_vom = sum(nausea_vom*patient_count_imp),
         headache = sum(headache*patient_count_imp),
         diarrhea = sum(diarrhea*patient_count_imp),
         congestion = sum(congestion*patient_count_imp),
         sneezing = sum(sneezing*patient_count_imp)) |>
  distinct(county_fips, week_date, age_grp, state_fips, month, flu, fever, 
           myalgia, hypoxemia, short_breath, cough, bronchitis, chest_pain, 
           nausea_vom, sore_throat, fatigue, diarrhea, headache, congestion, 
           sneezing, .keep_all = FALSE)

flu_symp_county_age$age_group <- factor(flu_symp_county_age$age_grp,
                           levels = c('0', '1', '2', '3', '4', '5'),
                           labels = c('0-4', '5-12', '13-17', '18-49', '50-64', '65+'))

# Save data
saveRDS(flu_symp_county_age, './tmp/flu_symp_county_age.rds')

# Bring to county-week level
flu_symp_county <- flu_symp_county_age |> group_by(county_fips, week_date, 
                                            state_fips, month) |>
  mutate(flu = sum(flu),
         fever = sum(fever),
         cough = sum(cough),
         sore_throat = sum(sore_throat),
         fatigue = sum(fatigue),
         myalgia = sum(myalgia),
         hypoxemia = sum(hypoxemia),
         short_breath = sum(short_breath),
         bronchitis = sum(bronchitis),
         chest_pain = sum(chest_pain),
         nausea_vom = sum(nausea_vom),
         headache = sum(headache),
         diarrhea = sum(diarrhea),
         congestion = sum(congestion),
         sneezing = sum(sneezing)) |>
  distinct(county_fips, week_date, state_fips, month, flu, fever, 
           myalgia, hypoxemia, short_breath, cough, bronchitis, chest_pain, 
           nausea_vom, sore_throat, fatigue, diarrhea, headache, congestion, 
           sneezing, .keep_all = FALSE)

# Save data
saveRDS(flu_symp_county, './tmp/flu_symp_county.rds')

# Bring to state-week level
flu_symp_state <- flu_symp_county_age |> group_by(state_fips, week_date, month) |>
  mutate(flu = sum(flu),
         fever = sum(fever),
         cough = sum(cough),
         sore_throat = sum(sore_throat),
         fatigue = sum(fatigue),
         myalgia = sum(myalgia),
         hypoxemia = sum(hypoxemia),
         short_breath = sum(short_breath),
         bronchitis = sum(bronchitis),
         chest_pain = sum(chest_pain),
         nausea_vom = sum(nausea_vom),
         headache = sum(headache),
         diarrhea = sum(diarrhea),
         congestion = sum(congestion),
         sneezing = sum(sneezing)) |>
  distinct(state_fips, week_date, month, flu, fever, 
           myalgia, hypoxemia, short_breath, cough, bronchitis, chest_pain, 
           nausea_vom, sore_throat, fatigue, diarrhea, headache, congestion, 
           sneezing, .keep_all = FALSE)

# Save data
saveRDS(flu_symp_state, './tmp/flu_symp_state.rds')

# Remove full flu data 
remove(flu_symp)

################################
# 4. MERGE ON ALL WEATHER DATA #
################################

# Re-load formatted data
flu_symp_county <- readRDS('./tmp/flu_symp_county.rds')
flu_symp_state <- readRDS('./tmp/flu_symp_state.rds')

# Load weather data
weather_week <- readRDS('./tmp/county_week_weather.rds') 
weather_state_week <- readRDS('./tmp/state_week_weather.rds') 

# Merge on the weather data
flu_symp_county <- left_join(flu_symp_county, weather_week, 
                             by = c('county_fips' = 'county_group',
                                    'week_date' = 'week_date'))
flu_symp_state <- left_join(flu_symp_state, weather_state_week, 
                             by = c('state_fips' = 'state_fips',
                                    'week_date' = 'week_date'))

##############################
# 5. MERGE ON ALL CAUSE DATA #
##############################

# Load all cause data
all_cause_week <- read.csv('raw/county_week_ac_v2.csv')

# Impute data
all_cause_week$all_cause_imp <- all_cause_week$all_cause
all_cause_week$all_cause_imp[all_cause_week$all_cause == '<=5'] <- 
  sample(1:5, length(all_cause_week$all_cause_imp[all_cause_week$all_cause == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
all_cause_week$all_cause_imp <- as.numeric(all_cause_week$all_cause_imp)

# Alter week date
all_cause_week <- all_cause_week |>
  # Convert week date to date
  mutate(week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", year_week)) |>
  mutate(week_date = ISOweek2date(week_date)) |>
  select(c(week_date, county_fips, all_cause_imp))

# Merge weekly all cause to disease data
flu_symp_county <- left_join(flu_symp_county, all_cause_week, 
                             by = c('week_date' = 'week_date',
                                    'county_fips' = 'county_fips'))

# Collapse all cause data to the state level
all_cause_state <- all_cause_week |>
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  group_by(week_date, state_fips) |> mutate(all_cause_state = sum(all_cause_imp)) |>
  distinct(week_date, state_fips, all_cause_state)

# Merge state all cause to disease data
flu_symp_state <- left_join(flu_symp_state, all_cause_state, 
                             by = c('week_date' = 'week_date',
                                    'state_fips' = 'state_fips'))

###########################################
# 5. RESTRICT DATA TO COUNTIES AND STATES #
###########################################

# Restrict to DC
flu_dc <- flu_symp_county |> filter(county_fips == '11001')

# Restrict to Maryland counties
flu_md_counties <- flu_symp_county |> filter(substr(county_fips, 1, 2) == '24')

# Restrict to Virginia counties
flu_va_counties <- flu_symp_county |> filter(substr(county_fips, 1, 2) == '51')

# Restrict to Maryland
flu_md <- flu_symp_state |> filter(state_fips == '24')

# Restrict to Virginia
flu_va <- flu_symp_state |> filter(state_fips == '51')

# Save data
saveRDS(flu_md, './tmp/flu_symp_maryland.rds')
saveRDS(flu_dc, './tmp/flu_symp_dc.rds')
saveRDS(flu_md_counties, './tmp/flu_symp_maryland_counties.rds')
saveRDS(flu_va, './tmp/flu_symp_virginia.rds')
saveRDS(flu_va_counties, './tmp/flu_symp_virginia_counties.rds')

################################################################################
################################################################################

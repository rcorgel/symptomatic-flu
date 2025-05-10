################################################################################
# File Name: load_weather_data                                                 #
#                                                                              #
# Purpose:   Load weather data from: DATA/visual_crossing_weather_bulk/.       #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load weather data                                              #
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
library(assertr)
library(ISOweek)
library(lubridate)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Documents/symptomatic-flu-proj/')

########################
# 2. LOAD WEATHER DATA #
########################

# Load weather data
temp_16 <- read_csv('raw/tot_2016_crosswalk_county.csv')
temp_17 <- read_csv('raw/tot_2017_crosswalk_county.csv')
temp_18 <- read_csv('raw/tot_2018_crosswalk_county.csv')
temp_19 <- read_csv('raw/tot_2019_crosswalk_county.csv')
temp_20 <- read_csv('raw/tot_2020_crosswalk_county.csv')
temp_21 <- read_csv('raw/tot_2021_crosswalk_county.csv')
temp_22 <- read_csv('raw/tot_2022_crosswalk_county.csv')
temp_23 <- read_csv('raw/tot_2023_crosswalk_county.csv')

# Combine to one dataset
temp_16 <- temp_16 |> dplyr::select(c(county_fips, datetime, dew, feelslike))
temp_17 <- temp_17 |> dplyr::select(c(county_fips, datetime, dew, feelslike))
temp_18 <- temp_18 |> dplyr::select(c(county_fips, datetime, dew, feelslike))
temp_19 <- temp_19 |> dplyr::select(c(county_fips, datetime, dew, feelslike))
temp_20 <- temp_20 |> dplyr::select(c(county_fips, datetime, dew, feelslike))
temp_21 <- temp_21 |> dplyr::select(c(county_fips, datetime, dew, feelslike))
temp_22 <- temp_22 |> dplyr::select(c(county_fips, datetime, dew, feelslike))
temp_23 <- temp_23 |> dplyr::select(c(county_fips, datetime, dew, feelslike))
weather <- rbind(temp_16, temp_17, temp_18, temp_19, temp_20, temp_21, temp_22, temp_23)
remove(temp_16, temp_17, temp_18, temp_19, temp_20, temp_21, temp_22, temp_23)

# Add leading 0's to county fips
weather <- weather |> 
  mutate(county_fips = ifelse(nchar(county_fips) == 4, 
                              paste0("0", county_fips), county_fips))

# Drop unknown counties and territories
weather$state_fips <- substr(weather$county_fips, 1, 2)
weather <- weather |> filter(county_fips != '99999') |>
  filter(state_fips != '78') |> filter(state_fips != '72') |>
  filter(state_fips != '60')

# Convert date variable from string to date
weather$date <- as.Date(weather$datetime, format = '%Y-%m-%d')

# Collapse data to county group to align with claims
# Load crosswalk
x_walk <- readRDS('./tmp/county_group_xwalk.rds')
x_walk <- x_walk |> rename('county_group' = 'county_fips')
# Join county to county group
weather <- left_join(weather, x_walk[, c(1, 2)],
                     by = c('county_fips' = 'fips'))
# Assert that all counties are merged on
weather |> assert(not_na, county_group)
# Collapse to county group
weather_group <- weather |> group_by(county_group, date) |>
  mutate(feelslike = mean(feelslike),
         dew = mean(dew)) |>
  distinct(county_group, date, feelslike, dew)

# Collapse to week level
weather_group$week_date_ISO <- as.character(ISOweek(weather_group$date))
weather_group$week_date <- sub("W", "", weather_group$week_date_ISO)
weather_group$week_date <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", weather_group$week_date)
weather_group$week_date <- ISOweek2date(weather_group$week_date)
weather_week <- weather_group |> group_by(county_group, week_date) |>
  mutate(feelslike = mean(feelslike),
         dew = mean(dew)) |>
  distinct(county_group, week_date, feelslike, dew)

# Collapse to state-week level
weather_group$state_fips <- substr(weather_group$county_group, 1, 2)
weather_state_week <- weather_group |> group_by(state_fips, week_date) |>
  mutate(feelslike = mean(feelslike),
         dew = mean(dew)) |>
  distinct(state_fips, week_date, feelslike, dew)

# Save weather data
saveRDS(weather_week, './tmp/county_week_weather.rds') 
saveRDS(weather_state_week, './tmp/state_week_weather.rds') 

################################################################################
################################################################################

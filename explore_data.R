################################################################################
# File Name: explore_data                                                      #
#                                                                              #
# Purpose:   Explore Influenza medical claims data from 2016-2020 for the      #
#            state of Maryland.                                                #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Explore the data                                               #    
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
library(sf)
library(ISOweek)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Documents/symptomatic-flu-proj/')

###################
# 2. EXPLORE DATA #
###################

# Load the formatted Maryland data
flu_md <- readRDS('./tmp/flu_symp_maryland.rds')

## PLOT FLU DYNAMICS ##

flu <- ggplot() + 
  geom_line(data = flu_md[2:(nrow(flu_md) - 1),], 
            aes(x = week_date, y = flu), 
            color = '#C7A939', linewidth = 1.2, alpha = 0.8) + 
  ylab('Flu Count') + xlab('Week') + 
  ggtitle('Maryland Flu Dynamics, 2016-2020') +
  theme_minimal()

## PLOT SYMPTOM DYNAMICS ##

colors <- c("Fever" = "#c36272", "Cough" = "#5a9374", "Sore Throat" = "#347dab")
symp <- ggplot()+ 
  geom_line(data = flu_md[2:(nrow(flu_md) - 1),], aes(x = week_date, y = sore_throat, 
                               color = "Sore Throat"), 
            linewidth = 1.2, alpha = 0.8) + 
  geom_line(data = flu_md[2:(nrow(flu_md) - 1),], aes(x = week_date, y = cough, 
                               color = "Cough"), 
            linewidth = 1.2, alpha = 0.8) +
  geom_line(data = flu_md[2:(nrow(flu_md) - 1),], aes(x = week_date, y = fever, 
                               color = "Fever"),
            linewidth = 1.2, alpha = 0.8) +
  ggtitle('Maryland Symptom Dynamics, 2016-2020') +
  theme_minimal() +
  labs(x = "Week", y = "Symptom Count", color = "Symptoms") +
  scale_color_manual(values = colors) + 
  theme(legend.position = 'bottom')

## PLOT WEATHER DYNAMICS ## 

weather_colors <- c("Temperature" = "#c36272", "Humidity" = "#9086ba")
ggplot()+ 
  geom_line(data = flu_md[2:(nrow(flu_md) - 1),], aes(x = week_date, y = feelslike, 
                                                      color = "Temperature"), 
            linewidth = 1.2, alpha = 0.8) + 
  geom_line(data = flu_md[2:(nrow(flu_md) - 1),], aes(x = week_date, y = dew, 
                                                      color = "Humidity"), 
            linewidth = 1.2, alpha = 0.8) +
  ggtitle('Maryland Weather Dynamics, 2016-2020') +
  theme_minimal() +
  labs(x = "Week", y = "Temperature", color = "Weather Variable") +
  scale_color_manual(values = weather_colors) + 
  theme(legend.position = 'bottom')

temp <- ggplot() + 
  geom_line(data = flu_md[2:(nrow(flu_md) - 1),], 
            aes(x = week_date, y = feelslike), 
            color = '#9086ba', linewidth = 1.2, alpha = 0.8) + 
  ylab('Temerature') + xlab('Week') + 
  ggtitle('Maryland Temperature Dynamics, 2016-2020') +
  theme_minimal()

## PLOT SPATIAL AND TEMPORAL REPRESENTITIVENESS ##

# Load the relevant all claims data and county-group xwalk
county_groups <- readRDS('./tmp/county_group_xwalk.rds')
all_cause_season <- read_csv('./raw/county_season_ac_v2.csv')
county_pop <- read_csv('./raw/co-est2019-alldata.csv')

# Create population county FIPS variable
county_pop$county_fips <- paste0(county_pop$STATE, county_pop$COUNTY, sep = '')
county_pop <- county_pop |> dplyr::select(c(county_fips, POPESTIMATE2019))

# Create weekly state percent variable
all_cause_week <- flu_md |>
  mutate(percent = all_cause_state / 6177224)
# Maryland 2020 population from: https://www.census.gov/library/stories/state-by-state/maryland.html

# Limit season data to Maryland and merge on county population
all_cause_season <- all_cause_season |> 
  filter(substr(county_fips, 1, 2) == '24') |>
  filter(season == '2019-2020')
all_cause_season <- left_join(all_cause_season, county_pop, 
                              by = c('county_fips' = 'county_fips'))
all_cause_season$percent <- all_cause_season$all_cause / 
  all_cause_season$POPESTIMATE2019

# Load maps
usa_albers_state <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds')) # convert to sf

# Collapse county map to county group
usa_albers_county_x_walk <- left_join(county_groups[, c(1,2)], usa_albers_county[, c(5, 6)], by = c('fips' = 'GEOID'))
usa_albers_county_group <- usa_albers_county_x_walk %>% 
  group_by(county_fips) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

# Limit map data to Maryland
usa_albers_county_group <- usa_albers_county_group |>
  filter(substr(county_fips, 1, 2) == '24')
usa_albers_state <- usa_albers_state |>
  filter(STATEFP == '24')

# Merge on map data
usa_albers_county_group <- left_join(usa_albers_county_group, all_cause_season,
                                     by = c('county_fips' = 'county_fips'))
# Create map
map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = percent, group = county_fips), color= 'black', linewidth = 0.10) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.3) +
  scale_fill_gradient('Percent    ', low = "white", high = 'black') + ggtitle('Percent of Population Represented, 2019-2020') + 
  theme_void() + theme(legend.position = 'right',
                       plot.title = element_text(size = 16, hjust = 0.5),
                       panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                       legend.title = element_text(size = 12),
                       legend.text = element_text(size = 12),
                       legend.key.height = unit(1.2, 'cm'),
                       legend.key.width = unit(0.5, "cm")) 

# Create line plot
# Filter data to appropriate time frame
all_cause_week_filt <- all_cause_week |>
  filter(week_date > '2016-08-29') |>
  filter(week_date < '2020-09-07')

ggplot() + 
  geom_line(data = all_cause_week_filt, 
            aes(x = week_date, y = percent), 
            color = 'black', linewidth = 1.2, alpha = 0.8) + 
  ylab('Percent Observed') + xlab('Week') + 
  ggtitle('Weekly Maryland Population Observed, 2016-2020') +
  theme_minimal()

################################################################################
################################################################################
  
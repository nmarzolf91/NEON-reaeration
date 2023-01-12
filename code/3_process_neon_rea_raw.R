## Header ----
## Script name: process_neon_rea_raw.R
##
## Purpose of script: process NEON reaeration data from raw download into usable chunks
##
## Author: Nick Marzolf
## Date Created: 2022-11-16
## Date Modified: 
## Email: nicholas.marzolf@duke.edu
##
## load packages:  
library(tidyverse)
library(dplyr)
library(ggplot2)
library(feather)
library(data.table)
library(glue)

# this function will use the core processing in KC_reaerationformat.R as a guide and template
# this function/script initiates much of the data collating and processing and merging


process_neon_rea_raw <- function(dir = 'data/raw/') {
  
  injection_summary <- summarise_neon_rea_fieldData()
  
  rea_dir <- glue::glue(dir, '/Reaeration field and lab collection/')
  dsc_dir <- glue::glue(dir, '/Discharge field collection/')
  
  for(i in 1:length(unique(injection_summary$eventID))) {
    
    # define the individual event
    id <- injection_summary$eventID[i]
    site <- substr(id, 1,4)
    date <- substr(id, 6, 15)
    
    # read in the data and filter to the specific event, based on the date
    ## Reaeration data
    event_rea_fieldData <- try(feather::read_feather(
      glue::glue(rea_dir, site, '/rea_fieldData.feather')) %>% 
        filter(lubridate::date(collectDate) == date)
    )
    if(inherits(event_fieldData, 'try-error')) {
      cat(glue::glue('No reaeration field data from eventID: {site}_{date}'))
    }
    
    ## plateau data
    event_plateauMeasurementFieldData <- try(feather::read_feather(
      glue::glue(dir, site, '/rea_plateauMeasurementFieldData.feather')) %>% 
        filter(lubridate::date(collectDate) == date)
    )
    if(inherits(event_plateauMeasurementFieldData, 'try-error')) {
      cat(glue::glue('No plateauMeasurementFieldData from eventID: {site}_{date}'))
    } 
    
    ## background conductivity
    event_backgroundFieldCondData <- try(feather::read_feather(
      glue::glue(dir, site, '/rea_backgroundFieldCondData.feather')) %>% 
        filter(lubridate::date(startDate) == date))
    if(inherits(event_backgroundFieldCondData, 'try-error')) {
      cat(glue::glue('No event_backgroundFieldCondData from eventID: {site}_{date}'))
    } 
    
    ## lab salt 
    event_externalLabDataSalt <- try(feather::read_feather(
      glue::glue(dir, site, '/rea_externalLabDataSalt.feather')) %>% 
        filter(lubridate::date(startDate) == date))
    if(inherits(event_externalLabDataSalt, 'try-error')) {
      cat(glue::glue('No externalLabDataSalt from eventID: {site}_{date}'))
    } 
    
    ## field salt
    event_backgroundFieldSaltData <- try(feather::read_feather(
      glue::glue(dir, site, '/rea_backgroundFieldSaltData.feather')) %>% 
        filter(lubridate::date(collectDate) == date))
    if(inherits(event_externalLabDataSalt, 'try-error')) {
      cat(glue::glue('No backgroundFieldSaltData from eventID: {site}_{date}'))
    } 
    
    ## Discharge field collection
    ### ADCP data
    event_dsc_fieldData <- try(feather::read_feather(
      glue::glue(dsc_dir, site, '/dsc_fieldDataADCP.feather')) %>% 
        filter(lubridate::date(collectDate) == date)
    )
    if(inherits(event_fieldData, 'try-error')) {
      cat(glue::glue('No dsc field data from eventID: {site}_{date}'))
    } 
    
    event_dsc_fieldDataADCP <- try(feather::read_feather(
      glue::glue(dsc_dir, site, '/dsc_fieldDataADCP.feather')) %>% 
        filter(lubridate::date(collectDate) == date)
    )
    if(inherits(event_fieldDataADCP, 'try-error')) {
      cat(glue::glue('No dsc ADCP field data from eventID: {site}_{date}'))
    } 
    
    event_dsc_individualFieldData <- try(feather::read_feather(
      glue::glue(dsc_dir, site, '/dsc_dsc_individualFieldData.feather')) %>% 
        filter(lubridate::date(collectDate) == date)
    )
    if(inherits(event_fieldDataADCP, 'try-error')) {
      cat(glue::glue('No dsc_individualFieldData from eventID: {site}_{date}'))
    } 
    
    
    
    
    
  } # end for loop    
} # end function

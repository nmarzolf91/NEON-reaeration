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
{library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(feather)
  library(data.table)
  library(glue)}

# this function will use the core processing in KC_reaerationformat.R as a guide and template
# this function/script initiates much of the data collating and processing and merging


process_neon_rea_raw <- function(dir = 'data/raw/') {
  
  injection_summary <- summarise_neon_rea_fieldData()
  
  rea_dir <- glue::glue(dir, 'Reaeration field and lab collection/')
  dsc_dir <- glue::glue(dir, 'Discharge field collection/')
  
  tracker <- matrix(nrow = length(unique(injection_summary$eventID)),
                    ncol = 12)
  
  colnames(tracker) <- c('eventID', 'Has bkgnd fieldcond', 'Has bkgnd fieldsalt', 'Has fielddata',
                         'Has plateau  meas', 'Has plateau sample', 'Has external lab salt',
                         'Has external lab gas', 'Has width', 'Has dsc field', 'Has dsc ind field',
                         'Has dsc ACDP')
  
  out_names <- data.frame(
    'eventID' ,
    'site' ,
    'date' ,
    'station' ,
    'station_dist_downstream' ,
    'slugTracerMass',
    'slugPourTime' ,
    'dripStartTime' ,
    'dripEndTime' ,
    'backgroundSalt' ,
    'plateauSalt' ,
    'meanPlatSalt' ,
    'inj_gas' ,
    'platGas' ,
    'meanPlatGas' ,
    'mean_width' ,
    'water_temp' ,
    'hoboID' ,
    'fieldDischarge' 
  )
  
  out <- data.frame(matrix(ncol = length(out_names),
                           nrow = length(unique(injection_summary$eventID))))
  
  for(i in 1:length(unique(injection_summary$eventID))) {
    
    # define the individual event
    id <- injection_summary$eventID[i]
    site <- substr(id, 1,4)
    date <- substr(id, 6, 15)
    
    tracker[i,1] <- id
    
    # read in the data and filter to the specific event, based on the date
    ## Reaeration data
    
    ### background conductivity
    event_backgroundFieldCondData <- try(feather::read_feather(
      glue::glue(rea_dir, site, '/rea_backgroundFieldCondData.feather')) %>% 
        filter(lubridate::date(startDate) == date))
    if(inherits(event_backgroundFieldCondData, 'try-error')) {
      cat(glue::glue('No event_backgroundFieldCondData from eventID: {site}_{date}'))
      tracker[i, 2] <- 'N'
    } else tracker[i, 2] <- 'Y'
    
    ### background field salt
    event_backgroundFieldSaltData <- try(feather::read_feather(
      glue::glue(rea_dir, site, '/rea_backgroundFieldSaltData.feather')) %>% 
        filter(lubridate::date(collectDate) == date))
    if(inherits(event_backgroundFieldSaltData, 'try-error')) {
      cat(glue::glue('No backgroundFieldSaltData from eventID: {site}_{date}'))
      tracker[i, 3] <- 'N'
    } else tracker[i, 3] <- 'Y'
    
    
    ### field Data
    event_rea_fieldData <- try(feather::read_feather(
      glue::glue(rea_dir, site, '/rea_fieldData.feather')) %>% 
        filter(lubridate::date(collectDate) == date)
    )
    if(inherits(event_rea_fieldData, 'try-error')) {
      cat(glue::glue('No reaeration field data from eventID: {site}_{date}'))
      tracker[i, 4] <- 'N'
    } else tracker[i, 4] <- 'Y'
    
    
    ### plateau field data
    event_plateauMeasurementFieldData <- try(feather::read_feather(
      glue::glue(rea_dir, site, '/rea_plateauMeasurementFieldData.feather')) %>% 
        filter(lubridate::date(collectDate) == date)
    )
    if(inherits(event_plateauMeasurementFieldData, 'try-error')) {
      cat(glue::glue('No plateauMeasurementFieldData from eventID: {site}_{date}'))
      tracker[i, 5] <- 'N'
    } else tracker[i, 5] <- 'Y'
    
    
    ### plateau sample field data
    event_plateauSampleFieldData <- try(feather::read_feather(
      glue::glue(rea_dir, site, '/rea_plateauSampleFieldData.feather')) %>% 
        filter(lubridate::date(startDate) == date)
    )
    if(inherits(event_plateauSampleFieldData, 'try-error')) {
      cat(glue::glue('No plateau sample field data from eventID: {site}_{date}'))
      tracker[i, 6] <- 'N'
    } else tracker[i, 6] <- 'Y'
    
    
    ### external lab salt 
    event_externalLabDataSalt <- try(feather::read_feather(
      glue::glue(rea_dir, site, '/rea_externalLabDataSalt.feather')) %>% 
        filter(lubridate::date(startDate) == date))
    if(inherits(event_externalLabDataSalt, 'try-error')) {
      cat(glue::glue('No externalLabDataSalt from eventID: {site}_{date}'))
      tracker[i, 7] <- 'N'
    } else tracker[i, 7] <- 'Y'
    
    
    ### external lab gas
    event_externalLabGas <- try(feather::read_feather(
      glue::glue(rea_dir, site, '/rea_externalLabDataGas.feather')) %>% 
        filter(lubridate::date(startDate) == date))
    if(inherits(event_externalLabGas, 'try-error')) {
      cat(glue::glue('No externalLabDataGas from eventID: {site}_{date}'))
      tracker[i, 8] <- 'N'
    } else tracker[i, 8] <- 'Y'
    
    
    ### width
    event_widthFieldData <- try(feather::read_feather(
      glue::glue(rea_dir, site, '/rea_widthFieldData.feather')) %>% 
        filter(lubridate::date(collectDate) == date))
    if(inherits(event_widthFieldData, 'try-error')) {
      cat(glue::glue('No width from eventID: {site}_{date}'))
      tracker[i, 9] <- 'N'
    } else tracker[i, 9] <- 'Y'
    
    
    ## Discharge field collection
    ### field data
    event_dsc_fieldData <- try(feather::read_feather(
      glue::glue(dsc_dir, site, '/dsc_fieldData.feather')) %>% 
        filter(lubridate::date(startDate) == date)
    )
    if(inherits(event_dsc_fieldData, 'try-error')) {
      cat(glue::glue('No dsc field data from eventID: {site}_{date}'))
      tracker[i, 10] <- 'N'
    } else tracker[i, 10] <- 'Y'
    
    ### ADCP data
    event_dsc_fieldDataADCP <- try(feather::read_feather(
      glue::glue(dsc_dir, site, '/dsc_fieldDataADCP.feather')) %>% 
        filter(lubridate::date(startDate) == date)
    )
    if(inherits(event_dsc_fieldDataADCP, 'try-error')) {
      cat(glue::glue('No dsc ADCP field data from eventID: {site}_{date}'))
      tracker[i, 11] <- 'N'
    } else tracker[i, 11] <- 'Y'
    
    ### individual data
    event_dsc_individualFieldData <- try(feather::read_feather(
      glue::glue(dsc_dir, site, '/dsc_individualFieldData.feather')) %>% 
        filter(lubridate::date(collectDate) == date)
    )
    if(inherits(event_dsc_individualFieldData, 'try-error')) {
      cat(glue::glue('No dsc_individualFieldData from eventID: {site}_{date}'))
      tracker[i, 12] <- 'N'
    } else tracker[i, 12] <- 'Y'
    
    
    if(tracker[i, 'Has bkgnd fieldsalt'] == 'N') {
      logger_data <- cbind(event_backgroundFieldCondData,
                           event_rea_fieldData)
    } else {
      logger_data <- cbind(event_backgroundFieldSaltData,
                           event_rea_fieldData) 
    }
    
    # start pulling data out
    ## define sample type in lab salt
    event_externalLabDataSalt$sampleType <- ifelse(event_externalLabDataSalt$saltSampleID %in% event_rea_fieldData$injectateSampleID,
                                                   'injectate',
                                                   ifelse(event_externalLabDataSalt$saltSampleID %in% event_backgroundFieldSaltData$saltBackgroundSampleID,
                                                          'background',
                                                          'plateau')
    )
    
    event_dsc_fieldData_calc <- try(stageQCurve::conv.calc.Q(stageData = event_dsc_fieldData,
                                                             dischargeData = event_dsc_individualFieldData)
                                    )
    mean_width <- readr::read_csv(glue::glue('data/derived/mean_width/{site}_mean_width.csv')) %>% 
      filter(eventID %in% id) %>% 
      pull(mean_width)
    
    stream_hydraulics <- event_dsc_individualFieldData %>% 
      arrange(stationNumber) %>%
      summarise(mean_depth_m = mean(waterDepth, na.rm = TRUE), 
                mean_velocity_ms = mean(averageVelocity, na.rm = TRUE)) %>% 
      mutate(mean_width_m = mean_width,
             field_Q = event_dsc_fieldData_calc$calcQ/1000)
      
      
    
    stations <- logger_data$namedLocation
    
    hobos <- logger_data %>% 
      select(namedLocation, stationToInjectionDistance)
    
    bkgnd_salt <- event_externalLabDataSalt %>% 
      dplyr::filter(namedLocation %in% stations,
                    sampleType == 'background')
    
    plat_salt <- event_externalLabDataSalt %>% 
      dplyr::filter(namedLocation %in% stations,
                    sampleType == 'plateau')
    
    mean_plat_salt <- plat_salt %>% 
      group_by(namedLocation) %>% 
      summarise(mean_salt = mean(finalConcentration, na.rm = TRUE))
    
    inj_gas <- logger_data$gasTracerType
    
    plat_gas <- event_externalLabGas %>% 
      filter(namedLocation %in% stations)
    
    mean_plat_gas <- plat_gas %>% 
      group_by(namedLocation) %>% 
      summarise(mean_gas = mean(gasTracerConcentration, na.rm = TRUE))
    
    out <- out %>% 
      add_row(
        eventID = id,
        siteID = site,
        date = date,
        
      )
  } # end for loop    
  
  readr::write_csv(data.frame(tracker),
                   'data/derived/data_availability.csv')
} # end function

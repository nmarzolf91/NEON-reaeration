estimate_hydraulics <- function(dir = 'data/raw/') {
  
  dsc_dir <- glue::glue(dir, 'Discharge field collection/')
  rea_dir <- glue::glue(dir, 'Reaeration field and lab collection/')
  
  dsc_sites <- list.files(dsc_dir)
  rea_sites <- list.files(rea_dir)
  
  hyd_sites <- intersect(dsc_sites, rea_sites)
  
  for(i in 1:length(hyd_sites)) {
    
    site <- hyd_sites[i]
    
    ## Discharge field collection
    ### field data
    event_dsc_fieldData <- try(feather::read_feather(
      glue::glue(dsc_dir, site, '/dsc_fieldData.feather')) %>% 
        mutate(eventID = paste(siteID, lubridate::date(collectDate), sep = '_'))
    )
    if(inherits(event_dsc_fieldData, 'try-error')) {
      cat(glue::glue('No dsc field data from eventID: {site}')) 
    }
    
    ### ADCP data
    # event_dsc_fieldDataADCP <- try(feather::read_feather(
    #   glue::glue(dsc_dir, site, '/dsc_fieldDataADCP.feather')) %>% 
    #     filter(lubridate::date(startDate) == date)
    # )
    # if(inherits(event_dsc_fieldDataADCP, 'try-error')) {
    #   cat(glue::glue('No dsc ADCP field data from eventID: {site}_{date}'))
    #   tracker[i, 11] <- 'N'
    # } else tracker[i, 11] <- 'Y'
    
    ### individual data
    event_dsc_individualFieldData <- try(feather::read_feather(
      glue::glue(dsc_dir, site, '/dsc_individualFieldData.feather')) %>% 
        mutate(eventID = paste(siteID, lubridate::date(collectDate), sep = '_'))
    )
    if(inherits(event_dsc_individualFieldData, 'try-error')) {
      cat(glue::glue('No dsc_individualFieldData from eventID: {site}'))
    }
    
    
    event_dsc_fieldData_calc <- try(stageQCurve::conv.calc.Q(stageData = event_dsc_fieldData,
                                                             dischargeData = event_dsc_individualFieldData) %>% 
                                      mutate(calcQ_m3s = calcQ/1000,
                                             eventID = paste(siteID, lubridate::date(collectDate), sep = '_')) %>% 
                                      select(eventID, calcQ_m3s)
                                    )
    if(inherits(event_dsc_fieldData_calc, 'try-error')){
      cat(glue('{site} has no field discharge data'))
      next
    }
    
    mean_width <- try(readr::read_csv(glue::glue('data/derived/mean_width/{site}_mean_width.csv')))
    if(inherits(mean_width, 'try-error')){
      cat(glue::glue('No widths from {site}'))
    }
    
    stream_hydraulics <- event_dsc_individualFieldData %>% 
      group_by(eventID) %>% 
      arrange(stationNumber) %>%
      summarise(mean_depth_m = mean(waterDepth, na.rm = TRUE), 
                mean_velocity_ms = mean(averageVelocity, na.rm = TRUE)) %>% 
      left_join(., mean_width %>% 
                  select(eventID, mean_width_m = mean_width), 
                by = 'eventID') %>% 
      left_join(., event_dsc_fieldData_calc, by = 'eventID') %>% 
      mutate(est_mean_width_m = calcQ_m3s/(mean_depth_m*mean_velocity_ms))
    
    write_dir <- 'data/derived/hydraulics/'
    
    if(!dir.exists(write_dir))
      dir.create(write_dir)
    
    readr::write_csv(stream_hydraulics,
                     glue::glue(write_dir, '{site}_hydraulics.csv'))
    
  } # end for loop
} # end function

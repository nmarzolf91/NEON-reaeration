## Header ----
## Script name: summarise_neon_rea_fieldData.R
##
## Purpose of script: append the fieldData file from all reaeration experiments
##
## Author: Nick Marzolf
## Date Created: 2022-11-10
## Date Modified: 
## Email: nicholas.marzolf@duke.edu
##
## load packages:  
{library(tidyverse)
  library(dplyr)
  library(ggplot2)
}

summarise_neon_rea_fieldData <- function(dir = 'data/raw/Reaeration field and lab collection/'){
  
  sites <- list.files(dir)
  
  list_fieldDat <- list()
  
  for(s in 1:length(sites)){
    
    site <- sites[s]
    
    site_fD <- read_feather(glue(dir, site, '/rea_fieldData.feather'))
    
    site_fD_sub <- site_fD %>% 
      mutate(eventID = paste(site, lubridate::date(collectDate),
                             sep = '_'), .before = uid) %>% 
      select(-uid, -siteID, -coordinateUncertainty, -elevationUncertainty,-injectateSampleID,
             -injectateSampleCode,
             -samplingProtocolVersion, -recordedBy,-collectedBy, -publicationDate)  

    list_fieldDat[[s]] <- site_fD_sub
    
  } # end for loop
  
  fieldDat_summary <- bind_rows(list_fieldDat)
  return(fieldDat_summary)
  
} # end function

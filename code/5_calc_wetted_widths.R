## Header ----
## Script name: calc_wetted_widths.R
##
## Purpose of script: calculate mean wetted width from each injection and experiment
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
  library(glue)
}

calc_wetted_widths <- function(raw_dir = 'data/raw/Reaeration field and lab collection/') {
  sites <- list.files(raw_dir)
  
  for(i in 1:length(sites)) {
    site <- sites[i]
    
    widths <- read_feather(glue(raw_dir,'{site}/rea_widthFieldData.feather')) 
    
    mean_widths <- widths %>% 
      mutate(eventID = paste(siteID, lubridate::date(collectDate), sep = '_')) %>% 
      arrange(collectDate, widthMeasurementNumber) %>% 
      group_by(eventID) %>% 
      summarise(n_meas = n(),
                mean_width = mean(wettedWidth, na.rm = TRUE),
                median_width = median(wettedWidth, na.rm = TRUE),
                min_width = min(wettedWidth, na.rm = TRUE),
                max_width = max(wettedWidth, na.rm = TRUE),
                sd_width = sd(wettedWidth, na.rm = TRUE),
                se_width = sd_width/sqrt(n_meas))
    
    save_dir <- 'data/derived/'
    
    if(!dir.exists(save_dir))
      dir.create(save_dir)
    
    write_dir <- glue(save_dir, 'mean_width')
    
    if(!dir.exists(write_dir))
      dir.create(write_dir)
    
    write_csv(mean_widths,
              glue::glue(write_dir, '/{site}_mean_width.csv'))
    
  } # end for loop
} # end function


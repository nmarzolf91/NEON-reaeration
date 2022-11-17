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

# this function will use the core processing in KC_reaerationformat.R as a guide and template
# this function/script initiates much of the data collating and processing and merging


process_neon_rea_raw <- function(raw_dir = 'data/raw/Reaeration field and lab collection/') {
  sites <- list.files(raw_dir)
  
  for(i in 1:length(sites)) {
    site <- sites[i]
site
   

    
  } # end for loop    
    # read in the files collected from each 
} # end function
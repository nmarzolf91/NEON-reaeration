# get reaeration data product datasets from NEON

# this function will create a file structure in a data folder
# the file structure will be used in subsequent functions to evaluate reaeration data


get_reaeration_data <- function() {
  # load packages
  library(tidyverse)
  library(glue)
  library(feather)
  library(neonUtilities)
  
  # load helper functions
  source('code/helpers.R')
  source('code/neon_helpers.R')
  
  # load NEON site ID from NEON
  neon_metadata_all <- read_csv('https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv') %>% 
    data.frame()
  
  # subset sites to wadable stream sites where reaeration injections occur
  neon_wadable <- neon_metadata_all %>% 
    filter(field_site_subtype %in% c('Wadeable Stream'))
  
  # and pull the codes
  neon_wadable_codes <- neon_lotic[,'field_site_id']

  # define the product code
  product_code <- "DP1.20190.001"

  # get the data; this will take a while!
  # this function creates the file structure under a 'data' folder
  # each raw product from NEON will be housed in 'data/raw' folder
  # in 'data/raw', data is stored by site, with each site in a new folder
  # data is stored as feather files for the sake of efficiency 
  get_neon_product(product_codes = product_code)
}

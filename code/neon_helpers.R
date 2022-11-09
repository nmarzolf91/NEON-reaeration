# source('functions/get_avail_neon_product_sets.R')
# source('functions/serialize_list_to_dir.R')
# source('functions/read_neon_feathers.R')

get_neon_discharge <- function(gdrive = "MacroSheds",
                               gdrive_filename = 'simulated_discharge_2022-01-16',
                               file_ext = '.zip',
                               site_folder = 'Q_by_site') {
  
  # first ----
  # download THE MOST RECENT NEON discharge simulations from Google Drive
  # for me: NCSU GDrive/Shared Drives/macroSHEDS/data/Q_simulations
  # save to data folder within RProject
  
  # connect to MS gdrive and locate desired file
  ms_gdrive = googledrive::drive_find(gdrive)
  ## ms_gdrive = googledrive::shared_drive_find(gdrive)
  gdrive_file = paste0(gdrive_filename, file_ext)
  ms_sim_q = drive_find(gdrive_file)
  
  # download into temp folder, unzip
  temp <- tempfile(fileext = file_ext)
  drive_download(ms_sim_q$id, temp, overwrite = TRUE)
  out <- unzip(temp, exdir = tempdir())
  
  print('ATTENTION: printing README:')
  print(readLines(out[1]))
  
  # read data in ----
  this_dir <- paste0(tempdir(), '/', gdrive_filename, '/', site_folder, '/')
  print(paste('loading in site data from', this_dir))
  
  # the map function allows us to read in multiple files at once
  site_discharge <- purrr::map_dfr(
    list.files(this_dir,
               full.names = TRUE), readr::read_csv)
  
  
  return(site_discharge)
  ## NOTE: dont know what this is
  ## legit_q_bool = ! is.na(q_obs)
  ## q_sim[legit_q_bool] = q_obs[legit_q_bool]
}

# NEON webpage for water quality
# https://data.neonscience.org/data-products/DP1.20288.001

get_neon_product <- function(product_codes = 'DP0.20288.001',
                             dest_fp = NULL, site_filter = NULL) {
  
  # checking for user defined filepath
  if(is.null(dest_fp)) {
    dest_fp <- file.path(getwd(), 'data', 'raw')
  }
  
  if(!dir.exists(dest_fp)){
    # create direcotry if doesnt exists
    print(paste('Directory does not exist, creating directory here:', dest_fp))
    dir.create(dest_fp)
  }
  
  for(product_code in product_codes) {
    # get NEON product name, used for filepath
    product_name <- neonUtilities::getProductInfo(product_code)$productName
    data_fp <- file.path(dest_fp, product_name)
    products_n = length(product_codes)
    
    if(!dir.exists(data_fp)){
      # create direcotry if doesnt exists
      print(paste('Directory does not exist, creating directory here:', data_fp))
      dir.create(data_fp)
    }
    
    writeLines(paste('retrieving NEON data for', products_n, 'data products',
                     'and saving results at:\n', data_fp))
    
    # call to NEON for avilable data of code type
    avail_sets <- get_avail_neon_product_sets(product_code)
    
    # sites form this record
    avail_sites <- unique(avail_sets$site_name)
    
    if(!is.null(site_filter)) {
      print('values suppleid to site_filter argument, filtering queryable sites to those in site_filter values')
      avail_sites <- avail_sites[avail_sites %in% site_filter]
    }
    
    avail_sites_n <- length(avail_sites)
    if(avail_sites_n == 0) {
      stop('no available sites')
    }
    
    
    # download product for each site
    writeLines(paste('\nquerying for NEON product code:', product_code,
                     '\n  total sites to query:', avail_sites_n))
    for(j in 1:length(avail_sites)){
      writeLines(paste('  querying site', j, 'of', avail_sites_n))
      
      # Defines what site to retrieve
      site_name <- avail_sites[j]
      writeLines(paste('  site name:', site_name))
      
      # Download data product for the site
      tryCatch(
        expr = {
          data_pile <- neonUtilities::loadByProduct(dpID = product_code,
                                                    site = site_name,
                                                    package='basic',
                                                    check.size=FALSE)
        },
        error = function(e) {
          writeLines(paste('ERROR at', site_name, 'trying again with expanded package'))
          data_pile <- try(neonUtilities::loadByProduct(dpID = product_code,
                                                        site = site_name,
                                                        package='expanded',
                                                        check.size=FALSE))
          
        }
      )
      
      if(inherits(data_pile, 'try-error')) {
        writeLines(paste('no data downloaded skipping', site_name))
        next
      } else {
        # Create the file path for where the files will be saved. This can be changed
        # if you  want to save the files in a different location
        raw_data_dest <- file.path(data_fp, site_name)
        
        # Save the list of all dataframes for the site as individual feather files
        # in the file path raw_data_dest
        serialize_list_to_dir(data_pile, raw_data_dest)
      }
    }
  }
}


sample_neon_product <- function(product_codes = 'DP0.20288.001', product_name = 'neon_data',
                                dest_fp = NULL, site_filter = NULL) {
  # checking for user defined filepath
  if(is.null(dest_fp)) {
    dest_fp <- file.path(getwd(), 'data', 'raw')
  }
  
  data_fp <- file.path(dest_fp, product_name)
  products_n = length(product_codes)
  
  writeLines(paste('retrieving NEON data for', products_n, 'data products',
                   'and saving results at:\n', data_fp))
  
  for(product_code in first(product_codes)) {
    # call to NEON for avilable data of code type
    avail_sets <- get_avail_neon_product_sets(product_code)[1,]
    # sites form this record
    avail_sites <- unique(avail_sets$site_name)
    
    if(!is.null(site_filter)) {
      print('values suppleid to site_filter argument, filtering queryable sites to those in site_filter values')
      avail_sites <- avail_sites[avail_sites$site_name %in% site_filter,]
    }
    
    avail_sites_n <- length(avail_sites)
    if(avail_sites_n == 0) {
      stop('no available sites')
    }
    
    
    # download product for each site
    writeLines(paste('\nquerying for NEON product code:', product_code,
                     '\n  total sites to query:', avail_sites_n))
    for(j in 1:length(avail_sites)){
      writeLines(paste('  querying site', j, 'of', avail_sites_n))
      
      # Defines what site to retrieve
      site_name <- avail_sites[j]
      
      # Download data product for the site
      data_pile <- neonUtilities::loadByProduct(dpID = product_code,
                                                site = site_name,
                                                package='basic',
                                                check.size=FALSE)
      
      # Create the file path for where the files will be saved. This can be changed
      # if you  want to save the files in a different location
      ## raw_data_dest <- file.path(data_fp, site_name)
      # Save the list of all dataframes for the site as individual feather files
      # in the file path raw_data_dest
      ## serialize_list_to_dir(data_pile, raw_data_dest)
      
      return(data_pile)
    }
  }
}

get_neon_eval_q <- function(dest_fp = 'data/raw/qaqc/', dest_fn = 'neon_q_eval.csv') {
  destfile <- paste0(dest_fp, dest_fn)
  
  
  if(!dir.exists(dest_fp)){
    # create direcotry if doesnt exists
    print(paste('Directory does not exist, creating dir at:', dest_fp))
    dir.create(dest_fp)
  }
  
  # read in the discharge evaluation from Rhea et al. (in prep)
  download.file(
    url = 'https://www.hydroshare.org/resource/03c52d47d66e40f4854da8397c7d9668/data/contents/neon_q_eval.csv',
    destfile = destfile
  )
  
  q_eval <- read.csv(destfile)
  
  # filter the sites with bad discharge data
  good_q <- q_eval %>%
    filter(!final_qual %in% c('drift', 'regression_flag', 'Tier3'))
  return(good_q)
}


get_neon_site_data <- function(arg = 'details') {
  # create the site codes
  us_states <- USAboundaries::us_states()
  
  site_data <- macrosheds::ms_download_site_data() %>%
    filter(domain == 'neon',
           site_type == 'stream_gauge')
  
  site_geo <- site_data %>%
    sf::st_as_sf(coords = c('longitude', 'latitude'), crs = 4326)
  
  region_code = sf::st_join(site_geo, us_states) %>%
    dplyr::pull(state_abbr)
  
  sites <- site_data %>%
    dplyr::pull(site_code)
  
  # build a data frame with sites, site codes, and dates that will carry out of the loop
  # sp_code is the look-up site name in streampulse
  site_deets <- data.frame(
    site_code = sites,
    sp_code = paste0(region_code, '_', sites, 'temp'),
    start_date = NA,
    end_date = NA
  )
  if(arg == 'details') {
    return(site_deets)
  } else {
    return(site_data)
  }
}

# testing
## stream_width_sample <- sample_neon_product(product_codes = stream_width_pc)

## some NEON products are incomplete, and require 'check-and-fill' workflows to
## provide extant data where possible, two such are
##  pressure: provide NOAA gauge pressure where not in NEON
##  Q: provide USGS dataRetrieval Q for TOMB specifically
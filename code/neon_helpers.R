
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

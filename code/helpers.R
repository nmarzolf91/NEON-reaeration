## Title: Helper Functions
## Purpose: general functions, and staging ground for functions that may end up in a
##          more specific script at some point
## Author: Weston Slaughter
## Email: wslaughter@berkeley.edu
##

# TODO: get 'theme_nick' from NM, and define here manually, then import into rest of scripts

# functions from pre-existing scripts
# all scripts seem to start with these, no matter what
library(tidyverse)
library(dplyr)
library(ggplot2)

# source folder: functions

# function for converting from radians to degrees
# rad2deg <- function(x){
#   return(x * 180 / pi)
#   }
# 
# plot_metab_inputs <- function(df) {
#   df %>%
#     pivot_longer(cols = DO.obs:light,
#                  names_to = 'variable',
#                  values_to = 'value') %>%
#     ggplot(., aes(x = solar.time,
#                   y = value))+
#     geom_point()+
#     facet_wrap(variable ~ .,
#                scales = 'free')
# }
# full_days <- function(df) {
#   df_sum <- df %>%
#     dplyr::group_by(date = lubridate::date(solar.time)) %>%
#     dplyr::summarise(obs = dplyr::n(),
#               min = min(solar.time),
#               max = max(solar.time))
#   #return(df_sum)
# 
#   bad_days <- df_sum %>%
#     dplyr::filter(obs == 97) %>%
#     dplyr::pull(max)
#   #return(bad_days)
# 
#   # short_days <- df_sum %>%
#   #   dplyr::filter(obs < 96) %>%
#   #   dplyr::pull(max) %>%
#   #   date()
# 
#   df_out <- df %>%
#     dplyr::filter(!solar.time %in% bad_days) %>%
#     dplyr::mutate(date = date(solar.time)) %>%
#     #dplyr::filter(!date %in% short_days) %>%
#     dplyr::select(-date)
#   return(df_out)
# }

# below libs for serialize_list_to_dir func
library(neonUtilities)
library(glue)
library(httr)
library(jsonlite)
library(stringr)
library(readr)
library(feather)
library(tidyverse)


# serialize_list_to_dir: take a list of data tables downloaded from
# NEON and save that list as individual feather files for each dataframe in
# the list. Feather files are much like CSVs and can be loaded with
# read_feather(file_path) or saved with write_feather(dataframe, file_path)

serialize_list_to_dir <- function(l, dest){

  #l must be a named list
  #dest is the path to a directory that will be created if it doesn't exist

  #list element classes currently handled: data.frame, character

  elemclasses = lapply(l, class)

  handled = lapply(elemclasses,
                   function(x) any(c('character', 'data.frame') %in% x))

  if(! all(unlist(handled))){
    stop('Unhandled class encountered')
  }

  dir.create(dest, showWarnings=FALSE, recursive=TRUE)

  for(i in 1:length(l)){

    if('data.frame' %in% elemclasses[[i]]){

      fpath = paste0(dest, '/', names(l)[i], '.feather')
      write_feather(l[[i]], fpath)

    } else if('character' %in% x){

      fpath = paste0(dest, '/', names(l)[i], '.txt')
      readr::write_file(l[[i]], fpath)
    }
  }

}
# below libs for neon_feathers func
library(neonUtilities)
library(glue)
library(httr)
library(jsonlite)
library(stringr)
library(readr)
library(feather)
library(tidyverse)

# read_neon_feather will either read files into a list for one site or will
# read all sites and combine them into a list. The structure of this list mimics
# how a list downloaded from neon would look
#    file_path = where the data is stored. ether data/raw/product_name/site or
#                data/raw/product_name
#    by_site = TRUE or FALSE, TRUE to load one site or FALSE to load
#              all sites. Make sure to include the site you want to load in
#              file_path if by_site = TRUE

read_neon_feathers <- function(file_path, by_site){

  if(by_site == TRUE){

    neon_files <- list.files(file_path, full.names = TRUE)

    file_names <- list.files(file_path)

    neon_list <- lapply(neon_files, read_feather)

    names(neon_list) <- file_names

    return(neon_list)

  } else{

    neon_files <- list.files(file_path, full.names = TRUE, recursive = TRUE)
    file_names <- str_split_fixed(neon_files, '.feather', n = Inf)[,1]
    file_names <- str_split_fixed(file_names, '[/]', n = Inf)
    inv_file_names <- file_names[,dim(file_names)[2]]
    site_name <- file_names[,dim(file_names)[2]-1]

    inv_file_names <- unique(inv_file_names)
    site_names <- unique(site_name)

    #categoricalCodes
    #readme
    #validation
    #variables

    data_files <-  inv_file_names[!grepl('categoricalCodes|readme|validation|variables', inv_file_names)]

    final_list <- list()
    for(i in 1:length(data_files)){

      all_files <- tibble()
      for(s in 1:length(site_names)){

        path <- glue('{f}/{s}/{i}.feather',
                     f = file_path,
                     i = data_files[i],
                     s = site_names[s])

        one_site <- read_feather(path)

        all_files <- rbind(all_files, one_site)
      }

      all_site_files <- list(all_files)

      names(all_site_files) <- data_files[i]

      final_list <- c(final_list, all_site_files)
    }

    meta_data_files <-  inv_file_names[grepl('categoricalCodes|readme|validation|variables', inv_file_names)]

    for(i in 1:length(meta_data_files)){

      path <- glue('{f}/{s}/{i}.feather',
                   f = file_path,
                   i = meta_data_files[i],
                   s = site_names[1])

      meta_file <- read_feather(path)

      meta_list <- list(meta_file)
      names(meta_list) <- meta_data_files[i]

      final_list <- c(final_list, meta_list)
    }
    return(final_list)
  }

}

# libraries for get_avail_neon_product_sets
library(neonUtilities)
library(glue)
library(httr)
library(jsonlite)
library(stringr)
library(readr)
library(feather)
library(tidyverse)

# get_avail_neon_product_sets: retrieve all the sites and components available
# from the NEON website for the given data product. Neon organizes their data into
# year-month components for most data products, where every year-month has a data
# file (or many) with all corresponding meta data files for that month-year.

get_avail_neon_product_sets <- function(prodcode_full){

  #returns: tibble with url, site_name, component columns

  avail_sets = tibble()

  req = httr::GET(paste0("http://data.neonscience.org/api/v0/products/",
                         prodcode_full))

  txt = httr::content(req, as="text")

  neondata = jsonlite::fromJSON(txt, simplifyDataFrame=TRUE, flatten=TRUE)

  urls = unlist(neondata$data$siteCodes$availableDataUrls)

  avail_sets = stringr::str_match(urls,
                                  '(?:.*)/([A-Z]{4})/([0-9]{4}-[0-9]{2})') %>%
    as_tibble(.name_repair='unique') %>%
    rename(url=`...1`, site_name=`...2`, component=`...3`)

  return(avail_sets)
}

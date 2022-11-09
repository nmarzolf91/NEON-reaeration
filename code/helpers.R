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

# FindandCollect_airpres = function(lat, long, start_datetime, end_datetime) {
#   #get df of all available air pressure stations
#   tf = tempfile()
#   download.file("ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-history.txt", tf, mode="wb")
#   noaa.sites <- read.fwf(tf, skip = 22, header = F,
#                          # widths = c(6,-1,5,-1,30, 5, 3, 6, 8, 9, 8, 9, 8), comment.char = "",
#                          widths = c(6,-1,5,-45, 8, 9,-8, 9, 8), comment.char = "",
#                          col.names = c("USAF", "WBAN", "LAT", "LON", "BEGIN", "END"),
#                          # col.names = c("USAF", "WBAN", "STATION NAME", "CTRY", "ST", "CALL", "LAT", "LON", "ELEV(M)", "BEGIN", "END"),
#                          flush = TRUE, colClasses=c('USAF'='character', 'WBAN'='character'))
#   noaa.sites <- na.omit(noaa.sites)
#   #narrow them down to those within 5 lats/longs
#   noaa.sites <- noaa.sites %>%
#     mutate(LAT = as.numeric(as.character(LAT))) %>%
#     mutate(LON = as.numeric(as.character(LON))) %>%
#     filter(LAT < (lat + 5) & LAT > (lat - 5) & LON < (long + 5) & LON > (long - 5))
#   #filter by coverage, order by distance
#   pt1 <- cbind(rep(long, length.out = length(noaa.sites$LAT)),
#                rep(lat, length.out = length(noaa.sites$LAT)))
#   pt2 <- cbind(noaa.sites$LON, noaa.sites$LAT)
#   dist <- diag(geosphere::distm(pt1, pt2, fun=geosphere::distHaversine))/1000
#   noaa.sites$dist <- dist
#   tmp <- which((as.numeric(substr(noaa.sites$END,1,4)) >=
#                   as.numeric(substr(end_datetime, 1, 4))) &
#                  as.numeric(substr(noaa.sites$BEGIN,1,4)) <=
#                  as.numeric(substr(start_datetime, 1, 4)))
#   noaa.sites <- noaa.sites[tmp,]
#   noaa.sites <- noaa.sites[with(noaa.sites, order(dist)),]
#   yrs <- seq(as.numeric(substr(start_datetime, 1, 4)),as.numeric(substr(end_datetime, 1, 4)), by = 1)
#   for (i in 1:length(noaa.sites$dist)) {
#     k <- i
#     available <- vector(mode = 'logical', length = length(yrs))
#     USAF <- as.character(noaa.sites$USAF[i])
#     if(nchar(as.character(noaa.sites$WBAN[i])) == 5){
#       WBAN <- as.character(noaa.sites$WBAN[i])
#     } else {
#       WBAN <- paste0(0,as.character(noaa.sites$WBAN[i]))
#     }
#     y <- as.data.frame(matrix(NA, nrow = 1, ncol = 12))
#     for(j in 1:length(yrs)){
#       tf = tempfile()
#       res = tryCatch(suppressWarnings(download.file(paste0("ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-lite/",
#                                                            yrs[j], "/", USAF, "-", WBAN, "-", yrs[j], ".gz"), tf, mode="wb")),
#                      error=function(e){
#                        # message('NCDC download failed; trying next closest station')
#                        return('download failed')
#                      })
#       if(exists('res') && res == 'download failed'){
#         break #try next station
#       }
#       x = read.table(tf)
#       x[x==-9999] = NA
#       if(length(which(!is.na(x$V7))) >= 0.9 * length(x$V7)) {
#         available[j] <- TRUE
#         y <- rbind(x,y)
#       }else {
#         break #too many NAs, move to next station
#       }
#     }
#     if(length(yrs) == length(which(available))){
#       break #got one
#     }
#   }
#   y <- y[!is.na(y$V1),]
#   colnames(y) = c("y","m","d","h","air_temp","dewtemp","air_kPa","winddir","sindspeed","skycover","precip1h","precip6h")
#   y$air_kPa = y$air_kPa/100
#   y$air_temp = y$air_temp/10
#   y$DateTime_UTC = readr::parse_datetime(paste0(y$y,"-",
#                                                 sprintf("%02d",y$m),"-",sprintf("%02d",y$d)," ",sprintf("%02d",y$h),
#                                                 ":00:00"), "%F %T")
#   y <- y[with(y, order(DateTime_UTC)),]
#   y = tibble::as_tibble(y) %>% select(DateTime_UTC,air_temp,air_kPa)
#   ss = tibble::tibble(DateTime_UTC=seq(y$DateTime_UTC[1],
#                                        y$DateTime_UTC[nrow(y)], by=900))
#   xx = left_join(ss, y, by = "DateTime_UTC")
#   xx = mutate(xx, air_temp=zoo::na.approx(air_temp),
#               air_kPa=zoo::na.approx(air_kPa))
#   daterng = c(start_datetime, end_datetime)
#   xtmp = xx %>%
#     filter(DateTime_UTC>=daterng[1] & DateTime_UTC<=daterng[2]) %>%
#     mutate(air_mb = air_kPa*10) # convert to mbar
#   # select(xtmp, DateTime_UTC, air_kPa, air_temp)
#   # print(noaa.sites[k,])
#   return(select(xtmp, DateTime_UTC, air_mb, air_temp))
# }
# 
# # source folder: folder: prep_neon_data
# # source script: prep_daily_Q.R
# prep_daily_discharge <- function(site_data, domain = 'neon') {
#     # pull in site data
#     site_code <- read_csv('data/site_data/site_data.csv') %>%
#     # ensure domain is neon
#       filter(domain == 'neon') %>%
#     # extract site codes
#       pull(site_code)
# 
#     # loop through data
#     for(i in 1:length(site_code)){
# 
#       site <- site_code[i]
# 
#       d <- try(read_csv(glue('data/discharge/Q_simulations/simulated_discharge_2022-01-25/Q_by_site/{site}.csv')) %>%
#                  filter(year(date) > 2014) %>%
#                  mutate(discharge_m3s_sim = discharge_Ls_sim/1000) %>%
#                  select(site_code, date, discharge_m3s_sim))
# 
#       if(inherits(d, 'try-error')) {
#         next
#       }
# 
#       write_csv(d,
#                 glue('data/sm_ready_dailyQ/{site}_daily_Q.csv'))
#     }
# }

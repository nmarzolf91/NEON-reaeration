# get reaeration data product datasets

library(tidyverse)
library(neonUtilities)

source('code/helpers.R')
source('code/neon_helpers.R')


site_data <- get_neon_site_data(arg = 'data')


product_code <- "DP1.20190.001"
# product_name <- getProductInfo(dp)$productName

get_neon_product(product_codes = product_code,
                 product_name = product_name)



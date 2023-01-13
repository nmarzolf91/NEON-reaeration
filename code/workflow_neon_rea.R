## Header ----
## Script name: NEON_rea_workflow.R
##
## Purpose of script: workflow for analyzing NEON reaeration product
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
##
## clear the environment if needed
rm(list = ls())
##
## set the ggplot theme
source("C:/Users/Nick Marzolf/Desktop/Research/R code/theme_nick.R")
theme_set(theme_nick())

# Step 0: get necessary functions
source('code/helpers.R')
source('code/neon_helpers.R')

# Step 1: get the data
## this is a wrapper function around neonUtilities::loadByProduct()
## all files are returned as feather files
get_reaeration_data()

# get the field discharge product
get_neon_product(product_codes = 'DP1.20048.001')

# Step 2: summarize injection experiments
## Isolate each experiment and maintain some of the meta-data from these experiments
injection_summary <- summarise_neon_rea_fieldData()

# Step 3: Print the time series of conductivity from each experiment


# Step X: calculate wetted widths for each injection
calc_wetted_widths()



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

# get the data
get_reaeration_data()

# summarize injection experiments
summarise_neon_rea_fieldData()

# calculate wetted widths for each injection
calc_wetted_widths()



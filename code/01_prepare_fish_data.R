library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggplot2)
library(tidync)
library(sdmTMB)
library(marmap)
library(sf)

setwd("~/Dropbox/GitHub/wsg-choke-species")
source("code/helper_funs.R")
source("code/util_funs.R")

#Note: current taxa options for calculating metabolic index: "teleostei", "elasmobranchii", "perciformes", "gadidae"
#Read in 
species <- read_excel("data/species_table.xlsx")
species$common_name <- tolower(species$common_name)
species$scientific_name <- tolower(species$scientific_name)

#Ones w/ IPHC data
sub_species <- c("sablefish", "pacific cod", "pacific halibut", "yelloweye rockfish", "longnose skate", "big skate", "spiny dogfish", "rougheye rockfish")
sub <-filter(species, common_name %in% sub_species)

spcs <- tolower(sub$common_name)
sci_names <- tolower(sub$scientific_name)
taxas <- tolower(sub$MI_Taxa)
file_names <- spcs
#Ones to include IPHC data
for(i in 1:length(spcs)){
  spc <- spcs[i]
  print(spc)
  sci_name <- sci_names[i]
  file_name <- file_names[i]
  taxa <- taxas[i]
  try(prepare_data(spc=spc, sci_name=sci_name, mi=T, iphc=T, file_name=file_name, taxa=taxa))
  
}

#To not include IPHC data
sub_species <- c("sablefish", "pacific cod", "pacific halibut", "yelloweye rockfish", "longnose skate", "big skate", "spiny dogfish", "rougheye rockfish")
sub <-filter(species, !(common_name %in% sub_species))

spcs <- tolower(sub$common_name)
sci_names <- tolower(sub$scientific_name)
taxas <- tolower(sub$MI_Taxa)
file_names <- spcs
#Ones to include IPHC data
for(i in 1:length(spcs)){
  spc <- spcs[i]
  print(spc)
  sci_name <- sci_names[i]
  file_name <- file_names[i]
  taxa <- taxas[i]
  try(prepare_data(spc=spc, sci_name=sci_name, mi=T, iphc=F, file_name=file_name, taxa=taxa))
  gc()
  
}

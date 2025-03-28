library(dplyr)
library(tidyr)
library(purrr)
library(readxl)

# Pick a taxonomic group to simulate pcrit based on each temperature
setwd("~/Dropbox/GitHub/wsg-choke-species")
source("code/helper_funs.R")

##Create sequence of metabolic index values for marginal effects, so same for all
#Load data
files <- list.files(path = "data/processed_data/fish2", pattern = ".rds", full.names=T)
dat <- map(files,readRDS)
dat <- bind_rows(dat)

#Remove NAs (and remove IPHC by removing cpue)
dat <- dat  %>%
  drop_na(depth,year, mi1,mi2,mi3, X, Y, catch_weight)
#Taxa lookup
taxa <- read_excel("data/species_table.xlsx")
taxa$MI_Taxa <- tolower(taxa$MI_Taxa)
#Taxa lookup
taxa <- read_excel("data/species_table.xlsx")
taxa$MI_Taxa <- tolower(taxa$MI_Taxa)
taxa$common_name <- tolower(taxa$common_name)

#Assign
taxas <- unique(taxa$MI_Taxa)
W <- 1

#For range
t.range <- seq(min(dat$temp), max(dat$temp), length.out = 100)

##Pull from Tim's code
setwd("~/Dropbox/metabolic_index-main")
source("code/fit_model_funs.R")

#Run for each species
for(i in 1:nrow(taxa)){
taxa.2.use <- taxa[i,]$MI_Taxa
species.2.use <- taxa[i,]$common_name
print(species.2.use)
pcrit_df <- lookup_taxa_t(taxa.2.use, t.range = t.range, w.2.use = W)

# make an x_y positions df based on inner 50% range
positions50 <- data.frame(
  temp = c(pcrit_df$Temp, rev(pcrit_df$Temp)),
  ys = c(pcrit_df$lower50s, rev(pcrit_df$upper50s)),
  species = species.2.use,
  mi_taxa = taxa.2.use
)

positions90 <- data.frame(
  temp = c(pcrit_df$Temp, rev(pcrit_df$Temp)),
  ys = c(pcrit_df$lower90s, rev(pcrit_df$upper90s)),
  species = species.2.use,
  mi_taxa = taxa.2.use
)

if(i==1){
  dats50 <- positions50
  dats90 <- positions90
} else {
  dats50 <- rbind(dats50, positions50)
  dats90 <-  rbind(dats90, positions90)
}
}

setwd("~/Dropbox/GitHub/wsg-choke-species")
saveRDS(dats50, "data/lab_ests50.rds")
saveRDS(dats90, "data/lab_ests90.rds")


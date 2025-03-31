library(readxl)
library(sdmTMB)
library(dplyr)
library(tidyr)
library(gt)
library(openxlsx2)

setwd("~/Dropbox/GitHub/wsg-choke-species")

#Load functions
source("code/helper_funs.R")

#Output folder
output_folder <- "region_comp"

#Models
models <- c("model7", "model8", "model9", "model10", "model11","model12", "model13", "model14","model15", "model16")
aic_table = as.data.frame(matrix(NA, 1, length(models)+5))
dat_names <- c("cc", "bc", "goa", "ebs", "coastwide")
species <- read_excel("data/species_table.xlsx")
species$common_name <- tolower(species$common_name)
species$scientific_name <- tolower(species$scientific_name)
species <- unique(species$common_name)

for(i in 1:length(species)) {
  this_species = species[i]
  for(h in 1:length(dat_names)){
    temp <- matrix(NA, 1, length(models)+5)
    this_dat <- dat_names[h]
    print(this_dat)
    dat <- try(readRDS(file = paste0("output/region_comp/", this_species, "_", this_dat, "_dat.rds")))
    if(class(dat)!="try-error"){
      for(j in 1:length(models)) {
        this_model <- models[j]
        fit <- try(readRDS(file = paste0("output/region_comp/", this_species,"_",this_dat, "_", this_model,".rds")))
        if(class(fit)!="try-error"){
          s <- try(sanity(fit, silent=TRUE))
          if(class(s)!="try-error"){
            if(s$hessian_ok + s$eigen_values_ok +s$nlminb_ok== 3){
              temp[1,j] <- AIC(fit)
            }
          } else{
            temp[1,j] <- NA
          }
        }
      }
      mins <- apply(temp, 1, min, na.rm=T)
      for(c in 1:length(models)){
        temp[1,c] <- abs(temp[1,c]-mins[1])
      }
      
      temp[1,length(models)+1] <- this_species
      temp[1,length(models)+2] <- this_dat
      temp[1,length(models)+3]  <- nrow(filter(dat, catch>0))
      temp[1,length(models)+4]  <- length(unique(dat$year))
      temp[1,length(models)+5] <- length(unique(dat$region))
      aic_table <- bind_rows(aic_table, as.data.frame(temp))
    }
  }
}

aic_table <- aic_table[2:nrow(aic_table),]
colnames(aic_table) <- c("model7", "model8", "model9", "model10", "model11","model12", "model13", "model14","model15", "model16","species", "data type", "N obs", "N years", "N regions")
#Round
aic_table <- aic_table %>%
  mutate(across(1:length(models), is.numeric))
aic_table <- aic_table %>%
  mutate(across(1:length(models), round, 5))


#Replace NA with --
aic_table[is.na(aic_table)] <- "--"
aic <- aic_table %>%
  gt() %>%
  tab_style(
    style = cell_fill(color = "gray"),
    locations = cells_body(columns = model7, rows = model7==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "indianred"),
    locations = cells_body(columns = model8, rows = model8==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "thistle2"),
    locations = cells_body(columns = model9, rows = model9==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "thistle"),
    locations = cells_body(columns = model10, rows = model10==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "thistle3"),
    locations = cells_body(columns = model11, rows = model11==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_body(columns = model12, rows = model12==0)
  ) %>%
      tab_style(
        style = cell_fill(color = "thistle2"),
        locations = cells_body(columns = model13, rows = model13==0)
      ) %>%
      tab_style(
        style = cell_fill(color = "thistle"),
        locations = cells_body(columns = model14, rows = model14==0)
      ) %>%
      tab_style(
        style = cell_fill(color = "thistle3"),
        locations = cells_body(columns = model15, rows = model15==0)
      ) %>%
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_body(columns = model16, rows = model16==0)
  )

#save
write_xlsx(aic_table, filename="output/", output_folder, "/aic_table.xlsx")
gtsave(aic, filename="output/", output_folder, "/aic_table.html")
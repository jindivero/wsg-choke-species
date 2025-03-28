library(readxl)
library(sdmTMB)
library(dplyr)
library(tidyr)
library(gt)
library(openxlsx2)

setwd("~/Dropbox/GitHub/wsg-choke-species")

#Load functions
source("code/helper_funs.R")

#Bottom trawl--quadratic depth
#Just bottom trawl data, and the main five models 
##Make AIC table
models <- c("m1", "m2", "m3", "m4", "m5")
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
        fit <- try(readRDS(file = paste0("output/region_comp/", this_species,"_",this_dat, "_model", j, ".rds")))
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
colnames(aic_table) <- c("model1", "model2", "model3", "model4", "model5", "species", "data type", "N obs", "N years", "N regions")
#Round
aic_table$model1 <- round(as.numeric(aic_table$model1), digits=1)
aic_table$model2 <- round(as.numeric(aic_table$model2), digits=1)
aic_table$model3 <- round(as.numeric(aic_table$model3), digits=1)
aic_table$model4 <- round(as.numeric(aic_table$model4), digits=1)
aic_table$model5 <- round(as.numeric(aic_table$model5), digits=1)

aic <- aic_table %>%
  gt() %>%
  tab_style(
    style = cell_fill(color = "gray"),
    locations = cells_body(columns = model1, rows = model1==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "darkgoldenrod1"),
    locations = cells_body(columns = model2, rows = model2==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_body(columns = model3, rows = model3==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue2"),
    locations = cells_body(columns = model4, rows = model4==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue3"),
    locations = cells_body(columns = model5, rows = model5==0)
  )

write_xlsx(aic_table, "output/region_comp/aic_table_region_comp_priors_goodonly.xlsx")
gtsave(aic, filename="output/region_comp/aic_table_region_comp_priors_goodonly.html")

##For all models
models <- c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12")
aic_table = as.data.frame(matrix(NA, 1, length(models)+5))
dat_names <- c("cc", "bc", "goa", "ebs", "coastwide")

for(i in 1:length(species)) {
  this_species = species[i]
  for(h in 1:length(dat_names)){
    temp <- matrix(NA, 1, length(models)+5)
    this_dat <- dat_names[h]
    print(this_dat)
    dat <- try(readRDS(file = paste0("output/region_comp/", this_species, "_", this_dat, "_dat.rds")))
    if(class(dat)!="try-error"){
      for(j in 1:length(models)) {
        fit <- try(readRDS(file = paste0("output/region_comp/", this_species,"_",this_dat, "_model", j, ".rds")))
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
colnames(aic_table) <- c("model1", "model2", "model3", "model4", "model5","model6", "model7", "model8","model9", "model10", "model11", "model12", "species", "data type", "N obs", "N years", "N regions")
#Round
aic_table$model1 <- round(as.numeric(aic_table$model1), digits=1)
aic_table$model2 <- round(as.numeric(aic_table$model2), digits=1)
aic_table$model3 <- round(as.numeric(aic_table$model3), digits=1)
aic_table$model4 <- round(as.numeric(aic_table$model4), digits=1)
aic_table$model5 <- round(as.numeric(aic_table$model5), digits=1)
aic_table$model6 <- round(as.numeric(aic_table$model6), digits=1)
aic_table$model7 <- round(as.numeric(aic_table$model7), digits=1)
aic_table$model8 <- round(as.numeric(aic_table$model8), digits=1)
aic_table$model9 <- round(as.numeric(aic_table$model9), digits=1)
aic_table$model10 <- round(as.numeric(aic_table$model10), digits=1)
aic_table$model11 <- round(as.numeric(aic_table$model11), digits=1)
aic_table$model12 <- round(as.numeric(aic_table$model12), digits=1)

#save
write_xlsx(aic_table, "output/region_comp/aic_table_region_comp_all.xlsx")

##The models with cubic depth--just bottom trawl
##Make AIC table
models <- c("m7", "m8", "m9", "m10", "m11")
aic_table = as.data.frame(matrix(NA, 1, length(models)+5))
dat_names <- c("cc", "bc", "goa", "ebs", "coastwide")
species <- read_excel("data/species_table.xlsx")
species$common_name <- tolower(species$common_name)
species$scientific_name <- tolower(species$scientific_name)
species <- unique(species$common_name)

for(i in 1:length(species)) {
  this_species =species[i]
  print(this_species)
  for(h in 1:length(dat_names)){
    temp <- matrix(NA, 1, length(models)+5)
    this_dat <- dat_names[h]
    print(this_dat)
    dat <- try(readRDS(file = paste0("output/region_comp/", this_species, "_", this_dat, "_dat.rds")))
    if(class(dat)!="try-error"){
      for(j in 1:length(models)) {
        fit <- try(readRDS(file = paste0("output/region_comp/", this_species,"_",this_dat, "_model", j, ".rds")))
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
colnames(aic_table) <- c("model7", "model8", "model9", "model10", "model11", "species", "data type", "N obs", "N years", "N regions")
#Round
aic_table$model7 <- round(as.numeric(aic_table$model7), digits=1)
aic_table$model8 <- round(as.numeric(aic_table$model8), digits=1)
aic_table$model9 <- round(as.numeric(aic_table$model9), digits=1)
aic_table$model10 <- round(as.numeric(aic_table$model10), digits=1)
aic_table$model11 <- round(as.numeric(aic_table$model11), digits=1)

#Replace NA with --
aic_table[is.na(aic_table)] <- "--"

aic <- aic_table %>%
  gt() %>%
  tab_style(
    style = cell_fill(color = "gray"),
    locations = cells_body(columns = model7, rows = model7==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "darkgoldenrod1"),
    locations = cells_body(columns = model8, rows = model8==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_body(columns = model9, rows = model9==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue3"),
    locations = cells_body(columns = model10, rows = model10==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "deepskyblue3"),
    locations = cells_body(columns = model11, rows = model11==0)
  )

##Just species with IPHC and cubic depth
models <- c("model7", "model8", "model9", "model10", "model11")
aic_table = as.data.frame(matrix(NA, 1, length(models)+5))
dat_names <- c("cc", "bc", "goa", "ebs", "coastwide", "cc _iphc", "bc _iphc", "goa _iphc", "ebs _iphc", "coastwide _iphc")
sub_species <- c("sablefish", "pacific cod", "pacific halibut", "yelloweye rockfish", "longnose skate", "big skate", "spiny dogfish", "rougheye rockfish")

for(i in 1:length(sub_species)) {
  this_species =sub_species[i]
  print(this_species)
  for(h in 1:length(dat_names)){
    temp <- matrix(NA, 1, length(models)+5)
    this_dat <- dat_names[h]
    print(this_dat)
    dat <- try(readRDS(file = paste0("output/region_comp/", this_species, "_", this_dat, "_dat.rds")))
    if(class(dat)!="try-error"){
      for(j in 1:length(models)) {
        model_name <- models[j]
        fit <- try(readRDS(file = paste0("output/region_comp/", this_species,"_",this_dat, "_", model_name, ".rds")))
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

#Round
aic_table <- aic_table[2:nrow(aic_table),]
colnames(aic_table) <- c("model7", "model8", "model9", "model10", "model11", "species", "data type", "N obs", "N years", "N regions")

aic_table$model7 <- round(as.numeric(aic_table$model7), digits=3)
aic_table$model8 <- round(as.numeric(aic_table$model8), digits=3)
aic_table$model9 <- round(as.numeric(aic_table$model9), digits=3)
aic_table$model10 <- round(as.numeric(aic_table$model10), digits=3)
aic_table$model11 <- round(as.numeric(aic_table$model11), digits=3)

#Replace NA with --
aic_table[is.na(aic_table)] <- "--"

aic <- aic_table %>%
  gt() %>%
  tab_style(
    style = cell_fill(color = "gray"),
    locations = cells_body(columns = model7, rows = model7==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "darkgoldenrod1"),
    locations = cells_body(columns = model8, rows = model8==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_body(columns = model9, rows = model9==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue3"),
    locations = cells_body(columns = model10, rows = model10==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "deepskyblue3"),
    locations = cells_body(columns = model11, rows = model11==0)
  )


write_xlsx(aic_table, "output/region_comp/aic_table_region_comp_priors_goodonly_iphc.xlsx")
gtsave(aic, filename="output/region_comp/aic_table_region_comp_priors_goodonly_iphc.html")

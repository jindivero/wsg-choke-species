remove.packages("sdmTMB")
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE,  ref="newbreakpt")
library(sdmTMB)
library(sp)
library(dplyr)
library(tidyr)
library(pals)
library(purrr)
library(gt)
library(openxlsx2)
library(stringr)

set.seed(9876)

setwd("~/Dropbox/GitHub/wsg-choke-species")

#Load functions
source("code/helper_funs.R")

#Load data
files <- list.files(path = "data/processed_data/fish2", pattern = ".rds", full.names=T)
dat <- map(files,readRDS)
dat <- bind_rows(dat)

#Combine NBS & EBS into a BS region?
combine_bs <- F
if(combine_bs){
  dat$region <- ifelse((dat$region=="ebs"|dat$region=="nbs"), "bs", dat$region)
}

#Remove NBS 
remove_nbs <- T
if(remove_nbs){
  dat <- filter(dat, region!="nbs")
}

#Remove Aleutian Islands
remove_ai <- T
if(remove_ai){
  dat <- filter(dat, region!="ai")
}

#remove iphc
remove_iphc <- F
if(remove_iphc){
  dat <- filter(dat, survey!="iphc")
}

#Remove other missing rows
dat <- dat  %>%
  drop_na(depth,year, mi1,mi2,mi3, X, Y)

#Combine catch data
dat$catch_weight_combined <- ifelse(is.na(dat$catch_weight), dat$cpue_weight, dat$catch_weight)
dat$catch_count_combined <- ifelse(is.na(dat$catch_numbers), dat$cpue_count, dat$catch_numbers)

#Remove weird depths
dat <- filter(dat, depth>0)

#Remove oxygen outliers
dat <- filter(dat, O2_umolkg<1500)

#Remove catch outliers?
remove_outlier <- T

#Years
dat <- filter(dat, year>=2003)

##Clean up survey names
dat <- dat %>%
  mutate(across('survey', str_replace, 'EBS', 'afsc_bsai'),across('survey', str_replace, 'GOA', 'afsc_goa') )

#Species to fit--just ones with IPHC data
species<- c("sablefish", "pacific cod", "pacific halibut", "yelloweye rockfish", "longnose skate", "big skate", "spiny dogfish", "rougheye rockfish")
#fit models
for (h in 1:length(species)) {
  this_species = species[h]
  print(this_species)
  dat.2.use <- dplyr::filter(dat, common_name == this_species)
  
  ##Use numbers for all except halibut use weight
  if(this_species=="pacific halibut"){
    dat.2.use$catch <- dat.2.use$catch_weight_combined
  } else {
    dat.2.use$catch <- dat.2.use$catch_count_combined
  }
  
  ##Remove missing catch
  dat.2.use <- drop_na(dat.2.use, catch)
  
  #Years for adding extra_time later
  min_year <- min(dat.2.use$year)
  max_year <- max(dat.2.use$year)
  
  ##Separate out regions for 
  #Make list to store filtered data
  subs <- list()
  regions <- unique(dat.2.use$region)
  for(i in 1:length(regions)){
    region.2.use <- regions[i]
    sub <- filter(dat.2.use, region==region.2.use)
    num <- filter(sub, catch>0)
    #Only include if there are more than 50 observations in the region
    if(nrow(num)>50){
      subs[[i]] <- sub
      names(subs)[[i]] <- region.2.use
    }
  }
  #Remove elements with no obs
  subs <- subs[(names(subs)=="cc"|names(subs)=="bc"|names(subs)=="ebs"|names(subs)=="goa")]
  
  #Add global dataset
  dat.global <- bind_rows(subs)
  subs[[length(subs)+1]] <- dat.global
  
  ##Add variables and clean up datasets
  dats <- list()
  for(j in 1:length(subs)){
    sub <- subs[[j]]
    #Quadratic and scaled temp
    sub$temp_scaled <- scale(sub$temperature_C)
    sub$temp_scaled2 <- sub$temp_scaled^2
    
    #Quadratic and scaled depth
    sub$log_depth_scaled <- scale(log(sub$depth))
    sub$log_depth_scaled2 <- sub$log_depth_scaled^2
    
    #Re-scale MIs & po2
    sub$mi1_s <- scale(sub$mi1)
    sub$mi2_s <- scale(sub$mi2)
    sub$mi3_s <- scale(sub$mi3)
    sub$po2_s <- scale(sub$po2)
    
    #Remove catch outliers
    if(remove_outlier==T){
      #Remove outliers = catch > 10 sd above the mean
      sub$catch_s <- scale(sub$catch)
      sub <- dplyr::filter(sub, catch_s <=20)
    }
    dats[[j]] <- sub
  }
  
  #Create objects to test
  #List data names
  dat_names <- unique(dat.global$region)
  dat_names[length(dat_names)+1] <- "coastwide"
  dat_names <- paste(dat_names, "_iphc")
  for(i in 1:length(dats)){
    sub <- dats[[i]]
    print(dat_names[i])
    #Make mesh
    bnd <- INLA::inla.nonconvex.hull(cbind(sub$X, sub$Y), 
                                     convex = -0.05)
    inla_mesh <- INLA::inla.mesh.2d(
      boundary = bnd,
      max.edge = c(150, 1000),
      offset = -0.1, # default -0.1
      cutoff = 50,
      min.angle = 5 # default 21
    )
    spde <- make_mesh(sub, c("X", "Y"), mesh = inla_mesh)
    
    priors <- sdmTMBpriors(
      matern_s = pc_matern(
        range_gt = 50, range_prob = 0.05, #A value one expects the range is greater than with 1 - range_prob probability.
        sigma_lt = 25, sigma_prob = 0.05 #A value one expects the marginal SD (sigma_O or sigma_E internally) is less than with 1 - sigma_prob probability.
      ),
      matern_st = pc_matern(
        range_gt = 50, range_prob = 0.05,
        sigma_lt = 25, sigma_prob = 0.05
      ),
      #  ar1_rho = normal(0.7,0.1),
      #tweedie_p = normal(1.5,0.2)
    )
    # refactor to avoid identifiability errors
    sub$survey <- as.factor(as.character(sub$survey))
    
    models <- c("m1", "m2", "m3", "m4", "m5", "m6")
    
    saveRDS(sub, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_dat.rds"))
    
    # Model 1: null
    print(paste(this_species))
    print("fitting m1")
      formula = "catch ~ -1 + survey+ log_depth_scaled+ log_depth_scaled2"
    start = Sys.time()
    m1 <- try(sdmTMB(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = "iid",
      control = sdmTMBcontrol(normalize = TRUE,
                              multiphase = TRUE, 
                              newton_loops = 2),
      extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m1, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model1.rds"))
    
    # Model 2: quadratic temp (uniform across regions)
    print("fitting m2")
    formula = "catch ~ -1 + survey +temp_scaled + temp_scaled2 + log_depth_scaled+ log_depth_scaled2"
    start = Sys.time()
    m2 <- try(sdmTMB(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = "iid",
      control = sdmTMBcontrol(normalize = TRUE,
                              multiphase = TRUE,
                              newton_loops = 3),
      extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m2, file = paste0("output/region_comp/", this_species,  "_", dat_names[i],"_model2.rds"))
    
    # Model 4: breakpoint MI low
    print("fitting m3")
    formula = "catch ~ -1 + survey +breakpt(mi1_s)+ log_depth_scaled+ log_depth_scaled2"
    start = Sys.time()
    m3 <-try(sdmTMB(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = "iid",
      control = sdmTMBcontrol(normalize = TRUE,
                              multiphase = TRUE,
                              newton_loops = 2),
      extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m3, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model3.rds"))
    
    # Model 5: Breakpoint(mi median)
    print("fitting m4")
    start = Sys.time()
      formula = "catch ~ -1 + survey +breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2"
    m4 <-try(sdmTMB(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = "iid",
      control = sdmTMBcontrol(normalize = TRUE,
                              multiphase = TRUE,
                              newton_loops = 2),
      extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m4, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model4.rds"))
    
    # Model 5: breakpoint(mi high)
    print("fitting m5")
    start = Sys.time()
      formula = "catch ~ -1 +survey +breakpt(mi3_s)+ log_depth_scaled+ log_depth_scaled2"
    m5 <- try(sdmTMB(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = "iid",
      control = sdmTMBcontrol(normalize = TRUE,
                              multiphase = TRUE,
                              newton_loops = 2),
      extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m5, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model5.rds"))
  }
  print("fitting m6")
  start = Sys.time()
  formula = "catch ~ -1 +survey+temp_scaled + temp_scaled2 +breakpt(po2_s)+ log_depth_scaled+ log_depth_scaled2"
  m6 <- try(sdmTMB(
    formula = as.formula(formula),
    mesh = spde,
    time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = "iid",
    control = sdmTMBcontrol(normalize = TRUE,
                            multiphase = TRUE,
                            newton_loops = 2),
    extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  print( Sys.time() - start )
  saveRDS(m6, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model6.rds"))
  # Model 1: null
  print(paste(this_species))
  print("fitting m7")
  if(dat_names[i]!="coastwide"){
    formula =   "catch ~ 1 + log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
  } else {
    formula = "catch ~ -1 + region+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
  }
  
  start = Sys.time()
  m7 <- try(sdmTMB(
    formula = as.formula(formula),
    mesh = spde,
    time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = "iid",
    control = sdmTMBcontrol(normalize = TRUE,
                            multiphase = TRUE, 
                            newton_loops = 2),
    extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  print( Sys.time() - start )
  saveRDS(m7, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model7.rds"))
  
  # Model 2: quadratic temp (uniform across regions)
  print("fitting m8")
  if(dat_names[i]!="coastwide"){
    formula =   "catch ~ 1+temp_scaled + temp_scaled2 + log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
  } else {
    formula = "catch ~ -1 + region +temp_scaled + temp_scaled2 + log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
  }
  start = Sys.time()
  m8 <- try(sdmTMB(
    formula = as.formula(formula),
    mesh = spde,
    time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = "iid",
    control = sdmTMBcontrol(normalize = TRUE,
                            multiphase = TRUE,
                            newton_loops = 3),
    extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  print( Sys.time() - start )
  saveRDS(m8, file = paste0("output/region_comp/", this_species,  "_", dat_names[i],"_model8.rds"))
  
  # Model 9: breakpoint MI low
  print("fitting m9")
  if(dat_names[i]!="coastwide"){
    formula =   "catch ~ 1 +breakpt(mi1_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
  } else {
    formula = "catch ~ -1 + region +breakpt(mi1_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
  }
  start = Sys.time()
  m9 <-try(sdmTMB(
    formula = as.formula(formula),
    mesh = spde,
    time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = "iid",
    control = sdmTMBcontrol(normalize = TRUE,
                            multiphase = TRUE,
                            newton_loops = 2),
    extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  print( Sys.time() - start )
  saveRDS(m9, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model9.rds"))
  
  # Model 10: Breakpoint(mi median)
  print("fitting m10")
  start = Sys.time()
  if(dat_names[i]!="coastwide"){
    formula =   "catch ~ 1 +breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
  } else {
    formula = "catch ~ -1 + region +breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
  }
  m10 <-try(sdmTMB(
    formula = as.formula(formula),
    mesh = spde,
    time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = "iid",
    control = sdmTMBcontrol(normalize = TRUE,
                            multiphase = TRUE,
                            newton_loops = 2),
    extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  print( Sys.time() - start )
  saveRDS(m4, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model10.rds"))
  
  # Model 11: breakpoint(mi high)
  print("fitting m11")
  start = Sys.time()
  if(dat_names[i]!="coastwide"){
    formula =   "catch ~ 1 +breakpt(mi3_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
  } else {
    formula = "catch ~ -1 +region +breakpt(mi3_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
  }
  m11 <- try(sdmTMB(
    formula = as.formula(formula),
    mesh = spde,
    time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = "iid",
    control = sdmTMBcontrol(normalize = TRUE,
                            multiphase = TRUE,
                            newton_loops = 2),
    extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  print( Sys.time() - start )
  saveRDS(m11, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model11.rds"))
  # Model 11: breakpoint(mi high)
  print("fitting m12")
  start = Sys.time()
  if(dat_names[i]!="coastwide"){
    formula =   "catch ~ 1 +breakpt(po2_s)+temp_scaled + temp_scaled2+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
  } else {
    formula = "catch ~ -1 +region +breakpt(po2_s)+temp_scaled + temp_scaled2+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
  }
  m12 <- try(sdmTMB(
    formula = as.formula(formula),
    mesh = spde,
    time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = "iid",
    control = sdmTMBcontrol(normalize = TRUE,
                            multiphase = TRUE,
                            newton_loops = 2),
    extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  print( Sys.time() - start )
  saveRDS(m12, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model12.rds"))
}
  gc()

##Make AIC table
models <- c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12")
aic_table = as.data.frame(matrix(NA, 1, length(models)+5))
dat_names <- c("cc_iphc", "bc_iphc", "goa_iphc", "ebs_iphc", "coastwide_iphc")

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
    style = cell_fill(color = "lightblue3"),
    locations = cells_body(columns = model4, rows = model4==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "deepskyblue4"),
    locations = cells_body(columns = model5, rows = model5==0)
  )

write_xlsx(aic_table, "output/region_comp/aic_table_region_comp_priors_goodonly.xlsx")
gtsave(aic, filename="output/region_comp/aic_table_region_comp_priors_goodonly.html")

#Include all, including those that may have failed to converge
aic_table = as.data.frame(matrix(NA, 1, length(models)+5))

for(i in 1:length(species)) {
  this_species = species[i]
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
          temp[1,j] <- try(AIC(fit), silent=T)
        } else{
          temp[1,j] <- NA
        }
      }
      mins <- apply(temp, 1, min, na.rm=T)
      for(c in 1:length(models)){
        temp[1,c] <- try(abs(temp[1,c]-mins[1]))
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
write.csv(aic_table, "output/region_comp/aic_table_region_comp_priors_all.csv")
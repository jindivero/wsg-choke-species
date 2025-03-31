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

#Combine catch data
dat$catch_weight_combined <- ifelse(is.na(dat$catch_weight), dat$cpue_weight, dat$catch_weight)
dat$catch_count_combined <- ifelse(is.na(dat$catch_numbers), dat$cpue_count, dat$catch_numbers)

#Dataframe for doing depth
dat_depth <- dat

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
  drop_na(depth,year, survey, po2,mi1,mi2,mi3, X, Y)

#Remove weird depths
dat <- filter(dat, depth>0)
dat_depth <- filter(dat_depth, depth>0)

#Remove oxygen outliers
dat <- filter(dat, O2_umolkg<1500)

#Remove catch outliers?
remove_outlier <- T

#Years
dat <- filter(dat, year>=2003)

##Clean up survey names
dat <- dat %>%
  mutate(across('survey', str_replace, 'EBS', 'afsc_bsai'),across('survey', str_replace, 'GOA', 'afsc_goa') )

#Spatio-temporal variation--if T, uses IID; if F, adds an annual fixed year effect 
spatio_temp <- F

#Filter depths?
filter_depth <- T

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
    dat_depth$catch <- dat_depth$catch_weight_combined
    dat_depth$catch<- replace(dat_depth$catch, dat_depth$catch == "Inf", NA)
  } else {
    dat.2.use$catch <- dat.2.use$catch_count_combined
    dat_depth$catch <- dat_depth$catch_count_combined
    dat_depth$catch <- replace(dat_depth$catch, dat_depth$catch == "Inf", NA)
  }
  
  ##Remove missing catch
  dat.2.use <- drop_na(dat.2.use, catch)
  
  #Years for adding extra_time later
  min_year <- min(dat.2.use$year)
  max_year <- max(dat.2.use$year)
  
  #Get depth habitat
  if(filter_depth){
    dat_depth2 <- filter(dat_depth, common_name==this_species)
    dat_depth2 <- dat_depth2  %>% drop_na(depth, catch)
    dat_depth2$catch <- replace(dat_depth2$catch, dat_depth2$catch == "Inf", NA)
    
    # Sort by depth
    dat_depth2 <- dat_depth2[order(dat_depth2$depth), ]
    
    #Calculate the cumulative sum of catch by depth
    dat_depth2$cumsum_catch <- cumsum(dat_depth2$catch)
    
    #Calculate the proportional cumulative sum
    dat_depth2$prop_cumsum_catch <- dat_depth2$cumsum_catch / sum(dat_depth2$catch, na.rm=T)
    
    # Find the index of the closest value to 98% (0.97) in prop_cumsum_var1
    closest_index <- which.min(abs(dat_depth2$prop_cumsum_catch - 0.99))
    
    # Get the depth at 97% catch
    closest_value <- dat_depth2$depth[closest_index]
    
    #Print
    print(paste("Depth filtered to", closest_value, "m"))
    dat.2.use$species_depth_habitat <- closest_value
    
    dat.2.use <- filter(dat.2.use, depth<closest_value)
  }
  
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
    sub$log_depth_scaled3 <- sub$log_depth_scaled^3
    
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
    sub$region <- as.factor(as.character(sub$region))
    sub$year <- as.factor(as.character(sub$year))
    
  #  models <- c("m7", "m8", "m9","m10", "m11", "m12")
    
    saveRDS(sub, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_dat.rds"))
    
    # Model 7: null
    print(paste(this_species))
    print("fitting m7")
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1+survey+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula = "catch ~ -1 + survey+year+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      } else{
        st_type="iid"
        formula = "catch ~ -1 + survey+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      }
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
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(normalize = TRUE,
                              multiphase = TRUE, 
                              newton_loops = 2)
      # extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m7, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model7.rds"))
    
    # Model 8: quadratic temp (uniform across regions)
    print("fitting m8")
    if(dat_names[i]!="coastwide"){
      formula =  "catch ~ -1+survey+temp_scaled + temp_scaled2 + log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula ="catch ~ -1 + survey +year+temp_scaled + temp_scaled2 + log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      } else{
        st_type="iid"
        formula ="catch ~ -1 + survey +temp_scaled + temp_scaled2 + log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      }
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
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(normalize = TRUE,
                              multiphase = TRUE,
                              newton_loops = 3)
      # extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m8, file = paste0("output/region_comp/", this_species,  "_", dat_names[i],"_model8.rds"))
    
    # Model 9: breakpoint MI low
    print("fitting m9")
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1 +survey+breakpt(mi1_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula = "catch ~ -1 + survey +year+breakpt(mi1_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      } else{
        st_type="iid"
        formula = "catch ~ -1 + survey +breakpt(mi1_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      }
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
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(normalize = TRUE,
                              multiphase = TRUE,
                              newton_loops = 2)
      #extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m9, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model9.rds"))
    
    # Model 10: Breakpoint(mi median)
    print("fitting m10")
    start = Sys.time()
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1 +survey+breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula =  "catch ~ -1 + survey+year+breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      } else{
        st_type="iid"
        formula =  "catch ~ -1 + survey+breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      }
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
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(normalize = TRUE,
                              multiphase = TRUE,
                              newton_loops = 2)
      # extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m10, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model10.rds"))
    
    # Model 11: breakpoint(mi high)
    print("fitting m11")
    start = Sys.time()
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1 +survey+breakpt(mi3_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula =  "catch ~ -1 + survey +year+breakpt(mi3_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      } else{
        st_type="iid"
        formula =  "catch ~ -1 + survey +breakpt(mi3_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      }
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
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(normalize = TRUE,
                              multiphase = TRUE,
                              newton_loops = 2)
      #extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m11, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model11.rds"))
    # Model 12: breakpoint(o2)
    print("fitting m12")
    start = Sys.time()
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1+ survey+breakpt(po2_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      st_type="iid"
    } else {
      formula =  "catch ~ -1 + survey+year +breakpt(po2_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      if(!spatio_temp){
        st_type <- "off"
      } else{
        st_type="iid"
      }
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
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(normalize = TRUE,
                              multiphase = TRUE,
                              newton_loops = 2)
      # extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m12, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model12.rds"))
    # Model 13: breakpoint MI low+quad temp
    print("fitting m13")
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1 +survey+breakpt(mi1_s)+ +temp_scaled + temp_scaled2+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula = "catch ~ -1 + survey +year+breakpt(mi1_s)+ +temp_scaled + temp_scaled2+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      } else{
        st_type="iid"
        formula = "catch ~ -1 + survey +breakpt(mi1_s)+ +temp_scaled + temp_scaled2+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      }
    }
    start = Sys.time()
    m13 <-try(sdmTMB(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(normalize = TRUE,
                              multiphase = TRUE,
                              newton_loops = 2)
      #extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m13, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model13.rds"))
    
    # Model 14: Breakpoint(mi median)
    print("fitting m14")
    start = Sys.time()
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1+survey+breakpt(mi2_s)+temp_scaled + temp_scaled2+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      st_type="iid"
    } else {
      formula =  "catch ~ -1 + survey+year +breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      if(!spatio_temp){
        st_type <- "off"
      } else{
        st_type="iid"
      }
    }
    m14 <-try(sdmTMB(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(normalize = TRUE,
                              multiphase = TRUE,
                              newton_loops = 2)
      # extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m14, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model14.rds"))
    
    # Model 15: breakpoint(mi high)
    print("fitting m15")
    start = Sys.time()
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1+survey +breakpt(mi3_s)+ temp_scaled + temp_scaled2+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula =  "catch ~ -1 + survey +year+breakpt(po2_s)+temp_scaled + temp_scaled2+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      } else{
        st_type="iid"
        formula =  "catch ~ -1 + survey +breakpt(po2_s)+temp_scaled + temp_scaled2+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      }
    }
    m15 <- try(sdmTMB(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(normalize = TRUE,
                              multiphase = TRUE,
                              newton_loops = 2)
      #extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m15, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model15.rds"))
    # Model 16: breakpoint(o2)
    print("fitting m16")
    start = Sys.time()
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1+survey +breakpt(po2_s)+ temp_scaled + temp_scaled2+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula =   "catch ~ -1+survey +year+breakpt(po2_s)+ temp_scaled + temp_scaled2+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      } else{
        st_type="iid"
        formula =   "catch ~ -1+survey +breakpt(po2_s)+ temp_scaled + temp_scaled2+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      }
    }
    m16 <- try(sdmTMB(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(normalize = TRUE,
                              multiphase = TRUE,
                              newton_loops = 2)
      # extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m16, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_model16.rds"))
    gc()
  }
}

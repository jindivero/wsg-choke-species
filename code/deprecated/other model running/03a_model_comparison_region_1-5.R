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
remove_iphc <- T
if(remove_iphc){
  dat <- filter(dat, survey!="iphc")
}

#Remove other missing rows
dat <- dat  %>%
  drop_na(depth,year, mi1,mi2,mi3, X, Y, catch_weight)

#Remove weird depths
dat <- filter(dat, depth>0)

#Remove oxygen outliers
dat <- filter(dat, O2_umolkg<1500)

##Species to do--full list
species <- c(unique(dat$common_name))

#Remove catch outliers?
remove_outlier <- T

#Years
dat <- filter(dat, year>=2003)

#fit models
for (h in 1:length(species)) {
  this_species = species[h]
  print(this_species)
  dat.2.use <- dplyr::filter(dat, common_name == this_species)
  
  #Years for adding extra_time later
  min_year <- min(dat.2.use$year)
  max_year <- max(dat.2.use$year)
  
  #Rename column
  dat.2.use$catch <- dat.2.use$catch_weight
  
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
    
    #Re-scale MIs
    sub$mi1_s <- scale(sub$mi1)
    sub$mi2_s <- scale(sub$mi2)
    sub$mi3_s <- scale(sub$mi3)
    
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
    
    models <- c("m1", "m2", "m3", "m4", "m5")
    
    saveRDS(sub, file = paste0("output/region_comp/", this_species, "_", dat_names[i], "_dat.rds"))
    
    # Model 1: null
    print(paste(this_species))
    print("fitting m1")
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ 1 + log_depth_scaled+ log_depth_scaled2"
    } else {
      formula = "catch ~ -1 + region+ log_depth_scaled+ log_depth_scaled2"
    }

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
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ 1+temp_scaled + temp_scaled2 + log_depth_scaled+ log_depth_scaled2"
    } else {
      formula = "catch ~ -1 + region +temp_scaled + temp_scaled2 + log_depth_scaled+ log_depth_scaled2"
    }
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
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ 1 +breakpt(mi1_s)+ log_depth_scaled+ log_depth_scaled2"
    } else {
      formula = "catch ~ -1 + region +breakpt(mi1_s)+ log_depth_scaled+ log_depth_scaled2"
    }
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
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ 1 +breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2"
    } else {
      formula = "catch ~ -1 + region +breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2"
    }
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
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ 1 +breakpt(mi3_s)+ log_depth_scaled+ log_depth_scaled2"
    } else {
      formula = "catch ~ -1 +region +breakpt(mi3_s)+ log_depth_scaled+ log_depth_scaled2"
    }
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
  gc()
}
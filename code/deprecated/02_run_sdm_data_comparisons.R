remove.packages("sdmTMB")
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE,  ref="newbreakpt")
library(sdmTMB)
library(sp)
library(dplyr)
install.packages('pals')
library(pals)
library(purrr)

set.seed(9876)

setwd("~/Dropbox/GitHub/wsg-choke-species")

#Load functions
source("code/helper_funs.R")

#Load data
files <- list.files(path = "data/processed_data/fish", pattern = ".rds", full.names=T)
dat <- map(files,readRDS)
dat <- bind_rows(dat)

#Remove NAs
dat <- dat  %>%
  drop_na(depth,year, mi1,mi2,mi3, X, Y)

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

##Species to do--full list
species <- c(unique(dat$common_name))

##Species for IPHC vs bottom trawl comparison
species <- c("pacific cod", "sablefish", "yelloweye rockfish", "spiny dogfish", "pacific halibut", "rougheye rockfish")

#Remove catch outliers?
remove_outlier <- T
#fit models
for (h in 1:length(species)) {
  this_species = species[h]
  dat.2.use <- dplyr::filter(dat, common_name == this_species)
  
  #Years for adding extra_time later
  min_year <- min(dat.2.use$year)
  max_year <- max(dat.2.use$year)
  
  ##Separate out IPHC and bottom trawl data
  #Bottom trawl
  sub1 <- dat.2.use %>%
    drop_na(cpue_kg_km2)
  sub1$catch <- sub1$cpue_kg_km2
  sub1$type <- "bottom trawl"
  
  #IPHC
  sub2 <- dat.2.use %>%
    drop_na(cpue_count)
  if(this_species=="pacific halibut"){
    sub2$catch <- sub2$cpue_weight
    sub2$catch <- ifelse(sub2$catch=="Inf", NA, sub2$catch)
    sub2 <- sub2 %>%
      drop_na(catch)
  } else{
    sub2$catch <- sub2$cpue_count
  }
  sub2$type <- "iphc"
  
  #List datasets
  subs <- list(sub1,sub2)

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
  dat_names <- c("bottomtrawl", "iphc")
  
  for(i in 1:length(dats)){
  sub <- dats[[i]]
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
  
  # Model 1: null
  print(paste(this_species))
  print("fitting m1")
  m1 <- try(sdmTMB(
    formula = catch ~ -1 + region + year + log_depth_scaled + log_depth_scaled2,
    mesh = spde,
   # time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = "off",
   control = sdmTMBcontrol(normalize = TRUE,
                            multiphase = TRUE, 
                            newton_loops = 2),
    extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  
  saveRDS(m1, file = paste0("output/data_type/", this_species, "_", dat_names[i], "_model1.rds"))
  
  # Model 2: quadratic temp (uniform across regions)
  print("fitting m2")
  m2 <- try(sdmTMB(
    catch ~ -1 + region +year+temp_scaled + temp_scaled2 + log_depth_scaled+ log_depth_scaled2,
    mesh = spde,
   # time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = "off",
    control = sdmTMBcontrol(normalize = TRUE,
                            multiphase = TRUE,
                           newton_loops = 3),
   extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  saveRDS(m2, file = paste0("output/data_type/", this_species,  "_", dat_names[i],"_model2.rds"))

  # Model 4: breakpoint MI low
  print("fitting m3")
  m3 <- try(sdmTMB(
    catch ~ -1 + region + year+ breakpt(mi1_s)+ log_depth_scaled+ log_depth_scaled2,
    mesh = spde,
   # time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = "off",
    control = sdmTMBcontrol(normalize = TRUE,
                            multiphase = TRUE,
                           newton_loops = 2),
    extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  saveRDS(m3, file = paste0("output/data_type/", this_species, "_", dat_names[i], "_model3.rds"))
  
  # Model 5: Breakpoint(mi median)
  print("fitting m4")
  m4 <- try(sdmTMB(
    catch ~ -1 + region + year+breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2,
    mesh = spde,
  #  time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = "off",
    control = sdmTMBcontrol(normalize = TRUE,
                            multiphase = TRUE,
                           newton_loops = 2),
    extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  saveRDS(m4, file = paste0("output/data_type/", this_species, "_", dat_names[i], "_model4.rds"))
  
  # Model 5: breakpoint(mi high)
  print("fitting m6")
  m5 <- try(sdmTMB(
    catch ~ -1 + region + year+ breakpt(mi3_s)+ log_depth_scaled+ log_depth_scaled2,
    mesh = spde,
    #time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = "off",
    control = sdmTMBcontrol(normalize = TRUE,
                            multiphase = TRUE,
                            newton_loops = 2),
    extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  saveRDS(m5, file = paste0("output/data_type/", this_species, "_", dat_names[i], "_model5.rds"))
  
  }
}

##Make AIC table
aic_table = as.data.frame(matrix(NA, 1, length(models)+2))

for(i in 1:length(species)) {
temp <- matrix(NA, length(dat_names), length(models)+2)
this_species = species[i]
  
for(j in 1:length(models)) {
fit <- try(readRDS(file = paste0("output/data_type/", this_species, "_bottomtrawl_model", j, ".rds")))
fit2 <- try(readRDS(file = paste0("output/data_type/", this_species,"_iphc_model", j,".rds")))

if(class(fit)!="try-error"){
  s <- try(sanity(fit, silent=TRUE))
if(s$hessian_ok + s$eigen_values_ok + s$nlminb_ok == 3)
  temp[1,j] <- AIC(fit)
}

if(class(fit2)!="try-error"){
  s2 <- try(sanity(fit2, silent=TRUE))
if(s2$hessian_ok + s2$eigen_values_ok + s2$nlminb_ok == 3){
  temp[2,j] = AIC(fit2)
}
}
}
mins <- apply(temp, 1, min, na.rm=T)
for(c in 1:length(models)){
  temp[1,c] <- abs(temp[1,c]-mins[1])
  temp[2,c] <- abs(temp[2,c]-mins[2])
}

temp[1,length(models)+1] <- this_species
temp[2,length(models)+1] <- this_species
temp[1,length(models)+2] <- "bottom trawl"
temp[2,length(models)+2] <- "iphc"
aic_table <- bind_rows(aic_table, as.data.frame(temp))
}

aic_table <- aic_table[2:nrow(aic_table),]
colnames(aic_table) <- c("null", "quadratic temp", "quadratic temp: region", "breakpt(mi_low)", "breakpt(mi_mid)", "breakpt(mi_high)", "species", "data type")

write.csv(aic_table, "output/data_type/aic_table_datatype_priors.csv")

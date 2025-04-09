remove.packages("sdmTMB")
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE,  ref="newbreakpt")
library(sdmTMB)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

###wd
setwd("~/Dropbox/GitHub/wsg-choke-species")

#Load functions
source("code/helper_funs.R")

###Oxygen data
##Combine with in situ data
insitu <- readRDS("data/processed_data/o2/insitu_combined.rds")
#Just columns of interest
dat <- insitu[,c("survey", "depth", "year", "date", "latitude", "longitude", "X", "Y", "temperature_C", "do_mlpL", "salinity_psu", "sigma0_kgm3", "O2_umolkg", "event_id")]
#Drop if NA in O2
dat <- drop_na(dat, latitude, longitude, O2_umolkg, temperature_C, salinity_psu, year)
#Calculate pO2
#Calculate pO2 from umol kg
dat$po2 <- calc_po2_sat(salinity=dat$salinity_psu, temp=dat$temperature_C, depth=dat$depth, oxygen=dat$O2_umolkg, lat=dat$latitude, long=dat$longitude, umol_m3=T, ml_L=F)

##Add region
# load regional polygons
regions.hull <- readRDS("data/processed_data/regions_hull.rds")
#make dataframe an sf object
dat_df <-  st_as_sf(dat, coords = c("longitude", "latitude"), crs = st_crs(4326))
dat_df$latitude <- dat$latitude
dat_df$longitude <- dat$longitude
# cycle through all regions
region_list <- c("ai", "bc", "cc", "ebs", "goa")
dats <- list()
for (i in 1:length(region_list)) {
  region <- region_list[i]
  poly <- regions.hull[i,2]
  # pull out observations within each region
  region_dat  <- st_filter(dat_df, poly)
  region_dat$region <- paste(region_list[i])
  dats[[i]] <- as.data.frame(region_dat)
}
#Bind back together
dat <- bind_rows(dats)
dat <- unique(dat)

#Remove weird depths
dat <- filter(dat, depth>0)

#Remove oxygen outliers
dat <- filter(dat, O2_umolkg<1500)

#Remove ai from list
region_list <- region_list[-1]

#Log depth
dat$depth_ln <- log(dat$depth)

#Fit model just to year (each region sep)
year_models <- list()
for(i in 1:length(region_list)){
  dat.2.use <- filter(dat, region==region_list[i])
  dat.2.use$year <- as.factor(dat.2.use$year)
  #Create mesh
  ## Make Mesh and fit model ####
  spde <- make_mesh(data = dat.2.use,
                    xy_cols = c("X", "Y"),
                    cutoff = 45)
  #Fit model
    formula <- po2 ~ 1+s(depth_ln)+as.factor(year)
  
  m1 <- try(sdmTMB(
    formula = as.formula(formula),
    mesh = spde,
    time="year",
    data = dat.2.use,
    family = gaussian(),
    spatial = "on",
    spatiotemporal  = "iid"
  ))
  year_models[[i]] <- m1
}

#What are the annual effects?
for(i in 1:length(year_models)){
  m <- year_models[[i]]
  #Get the coefficients
  summary(m)
  #Create a dataframe of the coefficients
  if(i==1){
    coefs <- data.frame(coef(m))
    colnames(coefs) <- region_list[i]
    coefs$term <- row.names(coefs)
    coefs <- select(coefs, term, bc)
  }
  if(i>1){
    coefs.temp <- data.frame(coef(m))
    colnames(coefs.temp) <- region_list[i]
    coefs.temp$term <- row.names(coefs.temp)
    coefs <- left_join(coefs, coefs.temp, by="term")
  }
}

#What is the SD of spatial & spatiotmeporal variation
for(i in 1:length(year_models)){
  m <- year_models[[i]]
  #Get the coefficients
  a <- tidy(m, effects = "ran_pars")
  a <- a[,1:2]
  colnames(a) <- c("term", paste(region_list[i]))
  #Add to the year terms
  if(i==1){
  pars <- a
  } else {
  pars <- left_join(pars, a, by="term")
  }
}


###Predict to grid
##Load grids
load("data/nwfsc_grid.rda")
load("data/afsc_grid.rda")
ebs_grid <- filter(afsc_grid, survey=="Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey"|survey=="Eastern Bering Sea Slope Bottom Trawl Survey")       
goa_grid <- filter(afsc_grid, survey=="Gulf of Alaska Bottom Trawl Survey")
load("data/dfo_synoptic_grid.rda")

###Add depth to ebs and goa grids from NOAA bathymetry data
##Bathymetry data
# load bathymetry data
bathy_all <- readRDS("data/bathymetry_regions_from_grids.rds")
#Just goa and ebs
bathy_all <- bathy_all[4:5]

# fit models for each region
make_depth_model <- function(bathydat) {
  
  spde <- make_mesh(data = as.data.frame(bathydat), xy_cols = c("X", "Y"), n_knots = 300)
  
  depth_model <- sdmTMB(log(noaadepth) ~ 1,
                        data = as.data.frame(bathydat),
                        spatial = "on", 
                        mesh = spde,
                        family = gaussian()
  )
}
depth_models <- lapply(X = bathy_all,
                       FUN = make_depth_model)


# load all data
# cycle through the alaska grids
afsc_grids <- list(ebs_grid, goa_grid)

for (i in 1:length(afsc_grids)) {
  grid2use <- afsc_grids[[i]]
  #Get X and Y
  grid2use$longitude <- grid2use$lon
  grid2use$latitude <- grid2use$lat
  # UTM transformation
  dat_ll = grid2use
  sp::coordinates(dat_ll) <- c("lon", "lat")
  sp::proj4string(dat_ll) <- sp::CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = sp::spTransform(dat_ll, 
                            sp::CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # convert back from sp object to data frame
  dat3 = as.data.frame(dat_utm)
  dat3= dplyr::rename(dat3, X = coords.x1,
                       Y = coords.x2)
  dat3$lat <- dat3$latitude
  dat3$lon <- dat3$longitude
  model.2.use <- depth_models[[i]]
  
  # get predicted log(depth) for each observation, based on model fit to that region
  region_dat_predict <- predict(model.2.use, as.data.frame(dat3))
  # add the predicted log(depth) to the data frame
  region_dat_predict$depth_m <- exp(region_dat_predict$est)
  dat2 <- select(region_dat_predict, area, depth_m, survey, survey_domain_year, lon, lat)
  
  if(i==1){
  ebs_grid <- dat2
  } else {
  goa_grid <- dat2
  }
}

##Make a list of grids
grids <- list(dfo_synoptic_grid, nwfsc_grid, ebs_grid, goa_grid)

##Predict pO2 across grid
preds <- list()
for(i in 1:length(region_list)){
  grid <- grids[[i]]
  m <- year_models[[i]]
  region.2.use <- region_list[i]
  dat.2.use <- filter(dat, region==region.2.use)
  print(region_list[i])
  #Remove negative depths
  grid <- filter(grid, depth_m>0)
  #Make depth ln
  grid$depth_ln <- log(grid$depth_m)
  #Get X and Y
  grid$longitude <- grid$lon
  grid$latitute <- grid$lat
  # UTM transformation
  dat_ll = grid
  sp::coordinates(dat_ll) <- c("lon", "lat")
  sp::proj4string(dat_ll) <- sp::CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = sp::spTransform(dat_ll, 
                            sp::CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # convert back from sp object to data frame
  grid = as.data.frame(dat_utm)
  grid = dplyr::rename(grid, X = coords.x1,
                        Y = coords.x2)
  #Expand by years
  grid <- grid %>%
    expand_grid(year = unique(dat.2.use$year)) %>%
    mutate(year = as.factor(year))
  grid$year <- as.factor(grid$year)
  
  #Add region column
  grid$region <- region_list[i]
  #Predict pO2
  p <- predict(m, newdata = grid)
  #Add to list
  preds[[i]] <- p
}

#Bind together
preds <- bind_rows(preds)

#Map spatial variation in each region
###Map#########
map_data <- rnaturalearth::ne_countries(scale = "large",
                                        returnclass = "sf",
                                        continent = "North America")

us_coast_proj <- sf::st_transform(map_data, crs = 32610)

for(i in 1:length(region_list)){
  grid2use<- filter(preds, region==region_list[i])
  grid2use<- filter(grid2use, year=="2021")
  print(region_list[i])

ggplot(us_coast_proj) + geom_sf() +
    geom_point(grid2use, mapping=aes(x=X*1000, y=Y*1000,colour=omega_s), size=0.1)+
    xlim(min(grid2use$X)*1000, max(grid2use$X)*1000)+
    ylim(min(grid2use$Y)*1000, max(grid2use$Y)*1000)+
    theme_minimal(base_size=12)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(axis.text.x=element_blank())+
    scale_colour_viridis()

ggsave(
  paste("output/plots/o2_omega_",region_list[i], ".png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

}

##Spatiotemporal effects for each year
for(i in 1:length(region_list)){
  grid2use<- filter(preds, region==region_list[i])
  print(region_list[i])
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(grid2use, mapping=aes(x=X*1000, y=Y*1000,colour=epsilon_st), size=0.1)+
    xlim(min(grid2use$X)*1000, max(grid2use$X)*1000)+
    ylim(min(grid2use$Y)*1000, max(grid2use$Y)*1000)+
      facet_wrap("year", ncol=4)+
    theme_minimal(base_size=12)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(axis.text.x=element_blank())+
    scale_colour_viridis()
  
  ggsave(
    paste("output/plots/o2_epsilon_",region_list[i], ".png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height = 11,
    units = c("in"),
    dpi = 600,
    limitsize = TRUE, bg="white"
  )
  
}

#Calculate the coefficient of variation of epsilon_st for each latitude & longitude combination
epsilon_cv <- preds %>%
  group_by(X, Y, region) %>%
  summarise(epsilon_st_sd = sd(epsilon_st),epsilon_st_mean = mean(epsilon_st),epsilon_st_var = var(epsilon_st), epsilon_st_cv = (sd(epsilon_st)/mean(epsilon_st))) %>%
  ungroup()

#plot
for(i in 1:length(region_list)){
  grid2use<- filter(epsilon_cv, region==region_list[i])
  print(region_list[i])
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(grid2use, mapping=aes(x=X*1000, y=Y*1000,colour=epsilon_st_var), size=0.1)+
    xlim(min(grid2use$X)*1000, max(grid2use$X)*1000)+
    ylim(min(grid2use$Y)*1000, max(grid2use$Y)*1000)+
    theme_minimal(base_size=12)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(axis.text.x=element_blank())+
    scale_colour_viridis()
  
  ggsave(
    paste("output/plots/o2_epsilon_cv_",region_list[i], ".png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height = 5,
    units = c("in"),
    dpi = 600,
    limitsize = TRUE, bg="white"
  )
  
}

#Time series
ggplot(filter(preds,year=="2021"), aes(x=omega_s, y=est))+
  geom_point()+
  facet_wrap("region", scales="free")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=12))

  ggsave(
    paste("output/plots/data_fit_mapping/map_",this_species,"_", this_dat,".png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height = 5,
    units = c("in"),
    dpi = 600,
    limitsize = TRUE, bg="white"
  )
  
ggplot(preds, aes(x=epsilon_st, y=est))+
  geom_point()+
  facet_wrap("region")

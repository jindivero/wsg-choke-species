remove.packages("sdmTMB")
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE,  ref="newbreakpt")
library(sdmTMB)
library(sp)
library(dplyr)
library(tidyr)
install.packages('pals')
library(pals)
library(purrr)
library(ggplot2)
library(ggstance)
library(readxl)
library(viridis)
library(ggstance)
library(ggh4x)
library(ggnewscale)

set.seed(9876)

setwd("~/Dropbox/GitHub/wsg-choke-species")

#ggplot themes
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

#Load functions
source("code/helper_funs.R")

#IPHC or bottom trawl?
iphc <- F
both <- F

####Plot data fit in models
##Data names
if(!iphc){
  dat_names <- c("cc", "bc", "goa", "ebs", "coastwide")
  species <- read_excel("data/species_table.xlsx")
  species$common_name <- tolower(species$common_name)
  species$scientific_name <- tolower(species$scientific_name)
  species <- unique(species$common_name)
}
if(iphc){
  dat_names <- c("cc _iphc", "bc _iphc", "goa _iphc", "ebs _iphc", "coastwide _iphc")
  species <- c("sablefish", "pacific cod", "pacific halibut", "yelloweye rockfish", "longnose skate", "big skate", "spiny dogfish", "rougheye rockfish")
}
if(both){
  dat_names <- c("cc", "bc", "goa", "ebs", "coastwide", "cc _iphc", "bc _iphc", "goa _iphc", "ebs _iphc", "coastwide _iphc")
  species <- read_excel("data/species_table.xlsx")
  species$common_name <- tolower(species$common_name)
  species$scientific_name <- tolower(species$scientific_name)
  species <- unique(species$common_name)
}
#Output folder
output_folder <- "presence_absence"

#Metabolic index models to use
mi_models2use <- c("model13", "model14", "model15")

##Pull data for each species and plot
map_data <- rnaturalearth::ne_countries(scale = "large",
                                        returnclass = "sf",
                                        continent = "North America")

us_coast_proj <- sf::st_transform(map_data, crs = 32610)
for(i in 1:length(species)) {
  this_species = species[i]
  print(this_species)
  for(h in 1:length(dat_names)){
    this_dat <- dat_names[h]
    print(this_dat)
    dat2plot <- try(readRDS(file = paste0("output/",output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
    if(is.data.frame(dat2plot)){
      ggplot(us_coast_proj) + geom_sf() +
        geom_point(filter(dat2plot,catch_weight>0),mapping=aes(x=X*1000, y=Y*1000,colour=survey), size=0.1)+
        xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
        ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
        facet_wrap("year", ncol=5)+
        theme_minimal(base_size=12)+
        xlab("Longitude")+
        ylab("Latitude")+
        ggtitle(paste(this_species, this_dat, sep=" "))+
        theme(axis.text.x=element_blank())
      
      
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
    }
  }
}

##Create sequence of metabolic index values for marginal effects, so same for all

#Load AIC table for model output
aic <- as.data.frame(read_excel(paste0("output/",output_folder, "/aic_table.xlsx")))
#Make first 6 columns in AIC numeric
aic[,1:(ncol(aic)-5)] <- sapply(aic[,1:(ncol(aic)-5)], as.numeric)
#Filter aic to only if a 0 in any of columns named in mi_models2use
aic <- filter(aic, aic[,mi_models2use[1]]==0 | aic[,mi_models2use[2]]==0 | aic[,mi_models2use[3]]==0)

#Species to run
species <- unique(aic$species)

#Load data
files <- list.files(path = "data/processed_data/fish2", pattern = ".rds", full.names=T)
dat <- map(files,readRDS)
dat <- bind_rows(dat)

#filter
dat <- filter(dat, common_name %in% species)

#Remove any rows with necessary data missing
dat <- dat %>%
  drop_na(depth, mi1, temperature_C, salinity_psu, X, Y, year)

#Remove weird depths
dat <- filter(dat, depth>0)

#Remove oxygen outliers
dat <- filter(dat, O2_umolkg<1500)
#Sequence of MI values
#Remove NAs (and remove IPHC by removing cpue)
mi1_pred = seq(min(dat$mi1), max(dat$mi1), length.out = 300)
mi2_pred = seq(min(dat$mi2), max(dat$mi2), length.out = 300)
mi3_pred = seq(min(dat$mi3), max(dat$mi3), length.out = 300)

#Taxa lookup
taxa <- read_excel("data/species_table.xlsx")
taxa$MI_Taxa <- tolower(taxa$MI_Taxa)
taxa$common_name <- tolower(taxa$common_name)

##Add range from phylogenetic imputation
dats50 <- readRDS("data/lab_ests50.rds")
dats90 <- readRDS("data/lab_ests90.rds")
#Closest temp to 12
temps <- unique(dats50$temp)
temp_12 <- temps[which(abs(temps-12)==min(abs(temps-12)))]
dats50_12 <- filter(dats50, temp==temp_12)
dats90_12 <- filter(dats90, temp==temp_12)
dats50_12$type <- rep(c("min", "max"), times=(nrow(dats50_12)/2))
dats90_12$type <- rep(c("min", "max"), times=(nrow(dats90_12)/2))
#pivot wider
dats50_12 <- pivot_wider(dats50_12, id_cols=species,values_from=ys, names_from=type)
dats90_12 <- pivot_wider(dats90_12, id_cols=species,values_from=ys, names_from=type)


##Loop to pull each set of comparisons that converged, pull the best-fitting model if includes MI, predict to a marginal effect prediction grid, and make a dataframe with them
#Create dataframe
cond_effects_preds = as.data.frame(matrix(NA, 1, 18))
for(i in 1:length(species)) {
  this_species = species[i]
  print(this_species)
  this_aic <- as.data.frame(filter(aic, species==this_species))
  this_aic$data.type <- this_aic[,"data type"]
  dat_names <- unique(this_aic$data.type)
  for(h in 1:length(dat_names)){
    this_dat <- dat_names[h]
    print(this_dat)
    this_data <- as.data.frame(filter(this_aic, data.type==this_dat))
    #Identify if MI is best-fitting model and pull that model
    mi_models <- this_data[,c(mi_models2use)]
    #If any cells are equal to 0
    if(any(mi_models==0 & !is.na(mi_models))){
      #Extract the name of the column which equals 0 
      best_model <- as.data.frame(which(mi_models==0, arr.ind=TRUE))
      best_model <- best_model$col
      best_model <- colnames(mi_models)[best_model]
      #Pull the model output file
      fit <- try(readRDS(file = paste0("output/", output_folder,"/", this_species,"_",this_dat,"_", best_model, ".rds")))
      #Pull the data file
      this_datframe <- try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
      ##Predict to marginal effect grid
      #Scale sequence of MI by the mean & sd of the dataset
      #Replace dat with this_datframe
      mi1_s_pred = (mi1_pred - mean(this_datframe$mi1))/sd(this_datframe$mi1)
      mi2_s_pred = (mi2_pred - mean(this_datframe$mi2))/sd(this_datframe$mi2)
      mi3_s_pred = (mi3_pred - mean(this_datframe$mi3))/sd(this_datframe$mi3)
      mi_best <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") mi1_pred else if(best_model=="model4"|best_model=="model10"|best_model=="model14") mi2_pred else if(best_model=="model5"|best_model=="model11"|best_model=="model15") mi3_pred
      mi_best_s <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") mi1_s_pred else if(best_model=="model4"|best_model=="model10"|best_model=="model14") mi2_s_pred else if(best_model=="model5"|best_model=="model11"|best_model=="model15") mi3_pred
      pred_year <- unique(this_datframe$year)[1]
      log_depth_scaled_mean <- mean(this_datframe$log_depth_scaled)
      log_depth_scaled_mean2 <- log_depth_scaled_mean^2
      log_depth_scaled_mean3 <- log_depth_scaled_mean^3
      this_region <- unique(this_datframe$region)[1]
      nd_po2 <- data.frame(mi1_s = mi1_s_pred,
                           mi2_s =mi2_s_pred,
                           mi3_s =mi3_s_pred,
                           mi1=mi1_pred,
                           mi2=mi2_pred,
                           mi3=mi3_pred,
                           mi_best=mi_best,
                           temp_scaled = 0,
                           temp_scaled2 = 0,
                           log_depth_scaled = log_depth_scaled_mean,
                           log_depth_scaled2 = log_depth_scaled_mean2,
                           log_depth_scaled3 = log_depth_scaled_mean3,
                           year = pred_year,
                           survey="iphc",
                           region=this_region)
      #nd_po2 <- convert_class(nd_po2)
      p1 <- try(predict(fit, newdata = nd_po2, se_fit = TRUE, re_form = NA))
      if(is.data.frame(p1)){
        p1$model <- paste(best_model)
        p1$data <- paste(this_dat)
        p1$species <- paste(this_species)
        p1 <-  p1 %>%
          mutate(est_sc= exp(est)/max(exp(est), na.rm=T),
                 est_se_sc1 = (exp(est)-exp(est_se))/max(exp(est), na.rm=T),
                 est_se_sc2 = (exp(est)+exp(est_se))/max(exp(est), na.rm=T))
        # p1$est_se_sc1 <- ifelse((p1$est_se_sc1=="NaN"|p1$est_se_sc1=="Inf"), NA, p1$est_se_sc1)
        #  p1$est_se_sc2 <- ifelse((p1$est_se_sc2=="NaN"|p1$est_se_sc2=="Inf"), NA, p1$est_se_sc2)
        cond_effects_preds <- bind_rows(cond_effects_preds, p1)
      }
    }
  }
}

##Clean up dataframe
#Remove first row
cond_effects_preds <- cond_effects_preds[-1,]
#Remove first 18 columns
cond_effects_preds <- cond_effects_preds[,19:ncol(cond_effects_preds)]


if(!iphc){
  cond_effects_preds$data <- factor(cond_effects_preds$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
  labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
  names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")
  
  ##Plot conditional effects
  ggplot(cond_effects_preds, aes(mi_best, y=est_sc))+
    facet_wrap("species", ncol=4)+
    geom_line(aes(colour=data))+
    # geom_ribbon(aes(ymin = est_se_sc1, ymax = est_se_sc2, fill=data), alpha=0.4)+
    scale_x_continuous(limits=c(0,15, by=5))+
    labs(x = bquote('Metabolic Index'), y = bquote('Conditional Effect Population Density'~(kg~km^-2)))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    scale_colour_manual(labels=labs, values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"), drop=FALSE)+
    theme(legend.key.height = unit(2, "lines"))+
    theme(panel.spacing = unit(1, "lines"))
  
  ggsave(
    paste0("output/", output_folder, "/cond_effects_region_mi_scaled_no_se.png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height = 11,
    units = c("in"),
    bg="white",
    dpi = 600,
    limitsize = TRUE
  )
  
}

if(iphc){
  ##Plot conditional effects
  ggplot(cond_effects_preds, aes(mi_best, y=est_sc))+
    facet_wrap("species", ncol=4)+
    geom_line(aes(colour=data))+
    #geom_ribbon(aes(ymin = est_se_sc1, ymax = est_se_sc2, fill=data), alpha=0.4)+
    scale_x_continuous(limits=c(0,5, by=1))+
    labs(x = bquote('Metabolic Index'), y = bquote('Conditional Effect Population Density'~(kg~km^-2)))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    theme(legend.key.height = unit(2, "lines"))+
    theme(panel.spacing = unit(1, "lines"))
}
if(both){
  cond_effects_preds$data <- factor(cond_effects_preds$data, levels=c("cc", "bc", "goa", "ebs", "coastwide", "cc _iphc", "bc _iphc", "goa _iphc", "ebs _iphc", "coastwide _iphc"))
  labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide", "British Columbia IPHC", "California Current IPHC", "Gulf of Alaska IPHC", "Eastern Bering Sea IPHC", "Coastwide IPHC")
  names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide", "cc _iphc", "bc _iphc", "goa _iphc", "ebs _iphc", "coastwide _iphc")
  
  ##Plot conditional effects
  ggplot(cond_effects_preds, aes(mi_best, y=est_sc))+
    facet_wrap("species", ncol=4)+
    geom_line(aes(colour=data))+
    # geom_ribbon(aes(ymin = est_se_sc1, ymax = est_se_sc2, fill=data), alpha=0.4)+
    scale_x_continuous(limits=c(0,15, by=5))+
    labs(x = bquote('Metabolic Index'), y = bquote('Conditional Effect Population Density'~(kg~km^-2)))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    scale_colour_manual(labels=labs, values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6","#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6" ))+
    theme(legend.key.height = unit(2, "lines"))+
    theme(panel.spacing = unit(1, "lines"))
  
  ggsave(
    paste0("output/", output_folder, "/cond_effects_region_mi_scaled_no_se.png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height = 11,
    units = c("in"),
    bg="white",
    dpi = 600,
    limitsize = TRUE
  )
  
}


ggsave(
  paste0("output/", output_folder, "/cond_effects_region_mi_scaled_no_se.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 11,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

##Unscaled
ggplot(cond_effects_preds, aes(mi_best, y=exp(est)))+
  facet_wrap("species")+
  geom_line(aes(colour=data))+
  #geom_ribbon(aes(ymin = est_se_sc1, ymax = est_se_sc2, fill=data), alpha=0.4)+
  xlim(0,5)+
  labs(x = bquote('Metabolic Index'), y = bquote('Conditional Effect Population Density'~(kg~km^-2)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=12))+
  theme(legend.position="top")+
  theme(text=element_text(size=15))

###Alternative by manually plotting
cond_effects_preds = as.data.frame(matrix(NA, 1, 18))
for(i in 1:length(species)) {
  this_species = species[i]
  print(this_species)
  this_aic <- as.data.frame(filter(aic, species==this_species))
  dat_names <- unique(this_aic$data.type)
  for(h in 1:length(dat_names)){
    this_dat <- dat_names[h]
    this_data <- as.data.frame(filter(this_aic, data.type==this_dat))
    #Identify if MI is best-fitting model and pull that model
    mi_models <- this_data[,c(mi_models2use)]
    #If any cells are equal to 0
    if(any(mi_models==0 & !is.na(mi_models))){
      #Extract the name of the column which equals 0 
      best_model <- as.data.frame(which(mi_models==0, arr.ind=TRUE))
      best_model <- best_model$col
      best_model <- colnames(mi_models)[best_model]
      #Pull the model output file
      fit <- try(readRDS(file = paste0("output/", output_folder,"/", this_species,"_",this_dat,"_", best_model, ".rds")))
      #Pull the data file
      this_datframe <- try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
      ##Predict to marginal effect grid
      #Scale sequence of MI by the mean & sd of the dataset
      mi1_s_pred = (mi1_pred - mean(this_datframe$mi1))/sd(this_datframe$mi1)
      mi2_s_pred = (mi2_pred - mean(this_datframe$mi2))/sd(this_datframe$mi2)
      mi3_s_pred = (mi3_pred - mean(this_datframe$mi3))/sd(this_datframe$mi3)
      mi_best <- if(best_model=="model3") mi1_pred else if(best_model=="model4") mi2_pred else mi3_pred
      mi_best_s <- if(best_model=="model3") mi1_s_pred else if(best_model=="model4") mi2_s_pred else mi3_s_pred
      nd_po2 <- data.frame(mi_s=mi_best_s,
                           mi=mi_best)
      
      #Pull breakpoint and slope
      pars <- as.data.frame(tidy(fit, effects="fixed", conf.int=T))
      slope <- filter(pars, grepl("slope", term))
      thresh <- filter(pars,grepl("breakpt", term))
      #Function to calculate breakpoint effect
      breakpoint_calc <- function(x, b_slope, b_thresh){
        if (x < b_thresh) {
          pred = x  * exp(b_slope)
          #pred = (x-b_thresh) * b_slope
        } else {
          pred=b_thresh * exp(b_slope)
          #pred = 0
        }
        return(pred)
      }
      
      nd_po2$est <- sapply(nd_po2$mi_s,breakpoint_calc, slope$estimate,thresh$estimate)
      nd_po2$est2 <- sapply(nd_po2$mi_s,breakpoint_calc, (slope$estimate-slope$std.error),(thresh$estimate))
      nd_po2$est3 <- sapply(nd_po2$mi_s,breakpoint_calc, (slope$estimate+slope$std.error), (thresh$estimate))
      nd_po2$model <- paste(best_model)
      nd_po2$data <- paste(this_dat)
      nd_po2$species <- paste(this_species)
      cond_effects_preds <- bind_rows(cond_effects_preds, nd_po2)
    }
  }
}

##Clean data
#Remove first row
cond_effects_preds <- cond_effects_preds[-1,]
#Remove first 18 columns
cond_effects_preds <- cond_effects_preds[,19:ncol(cond_effects_preds)]
cond_effects_preds$data <- factor(cond_effects_preds$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")

ggplot(cond_effects_preds, aes(x=mi, y=exp(est)))+
  facet_wrap("species", scales="free_y")+
  geom_line(aes(colour=data))+
  #geom_ribbon(aes(ymin = exp(est2), ymax = exp(est3), fill=data), alpha=0.4)+
  scale_x_continuous(limits=c(0,15, by=2))+
  labs(x = bquote('Metabolic Index'), y = bquote('Conditional Effect Population Density'~(kg~km^-2)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=12))+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(labels=labs, values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
  theme(legend.key.height = unit(2, "lines"))+
  theme(panel.spacing = unit(1, "lines"))

#_-------------------------------------------------------------------------------------------
##Breakpoint estimates of MI into line plots
for(i in 1:length(species)) {
  this_species = species[i]
  print(this_species)
  this_aic <- as.data.frame(filter(aic, species==this_species))
  this_aic$data.type <- this_aic[,"data type"]
  dat_names <- unique(this_aic$data.type)
  for(h in 1:length(dat_names)){
    this_dat <- dat_names[h]
    print(this_dat)
    this_data <- as.data.frame(filter(this_aic, data.type==this_dat))
    #Identify if MI is best-fitting model and pull that model
    mi_models <- this_data[,c(mi_models2use)]
    #If any cells are equal to 0
    if(any(mi_models==0 & !is.na(mi_models))){
      #Extract the name of the column which equals 0 
      best_model <- as.data.frame(which(mi_models==0, arr.ind=TRUE))
      best_model <- best_model$col
      best_model <- colnames(mi_models)[best_model]
      #Pull the model output file
      fit <- try(readRDS(file = paste0("output/", output_folder,"/", this_species,"_",this_dat,"_", best_model, ".rds")))
      #Pull the data file
      this_datframe <- try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
      #Calculate mean of this_datframe$mi based on best model
      mean_mi <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") mean(this_datframe$mi1) else if(best_model=="model4"|best_model=="model10"|best_model=="model14") mean(this_datframe$mi2) else mean(this_datframe$mi3)
      sd_mi <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") sd(this_datframe$mi1) else if(best_model=="model4"|best_model=="model10"|best_model=="model14") sd(this_datframe$mi2) else sd(this_datframe$mi3)
      #Pull breakpoint and slope
      pars <- as.data.frame(tidy(fit, effects="fixed", conf.int=T))
      slope <- filter(pars, grepl("slope", term))
      thresh <- filter(pars,grepl("breakpt", term))
      thresh$est <- (thresh$estimate*sd_mi)+mean_mi
      thresh$low <- ((thresh$estimate-thresh$std.error)*sd_mi)+mean_mi 
      thresh$high <- ((thresh$estimate+thresh$std.error)*sd_mi)+mean_mi
      if(!exists("bp_est")){
        bp_est <- data.frame(species=this_species, data=this_dat, model=best_model, breakpt=thresh$est, breakpt_se1=thresh$low, breakpt_se2=thresh$high)
      } else {
        bp_est <- bind_rows(bp_est, data.frame(species=this_species, model=best_model, data=this_dat, breakpt=thresh$est, breakpt_se1=thresh$low, breakpt_se2=thresh$high))
      }
    }
  }
}

##Plot line plot of breakpoint estimates
if(!iphc){
  bp_est$data <- factor(bp_est$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
  labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
  names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")
  
  ggplot(bp_est, aes(y=species, x=breakpt, colour=data))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data), size=3, position=ggstance::position_dodgev(height=0.4))+
    #Can add back shape
    geom_linerange(aes(xmin = breakpt_se1, xmax = breakpt_se2, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    #xlim(-1,10)+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"), drop=FALSE)+
    guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(5, "lines"))+ #Make more space between species
    #scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    xlab("Metabolic Index Breakpoint Estimate")+
    ylab("Species")+
    geom_vline(xintercept=0, linetype="dashed")+
    geom_vline(xintercept=1, linetype="dashed")
}
if(both){
  bp_est$data <- factor(bp_est$data, levels=c("cc", "bc", "goa", "ebs", "coastwide", "cc _iphc", "bc _iphc", "goa _iphc", "ebs _iphc", "coastwide _iphc"))
  labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide", "British Columbia IPHC", "California Current IPHC", "Gulf of Alaska IPHC", "Eastern Bering Sea IPHC", "Coastwide IPHC")
  names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide", "cc _iphc", "bc _iphc", "goa _iphc", "ebs _iphc", "coastwide _iphc")
  #Add column for IPHC versus not
  bp_est$type <- ifelse(grepl("iphc", bp_est$data), "iphc included", "no iphc")
  
  ggplot(bp_est, aes(y=species, x=breakpt, colour=data))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data, shape=type), size=3, position=ggstance::position_dodgev(height=0.4))+
    #Can add back shape
    geom_linerange(aes(xmin = breakpt_se1, xmax = breakpt_se2, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    #xlim(-1,10)+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6", "#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(5, "lines"))+ #Make more space between species
    #scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    xlab("Metabolic Index Breakpoint Estimate")+
    ylab("Species")+
    geom_vline(xintercept=0, linetype="dashed")+
    geom_vline(xintercept=1, linetype="dashed")
}

ggsave(
  paste0("output/", output_folder, "/breakpt_est_all.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height =11,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

##Capped at sensible numbers
if(!iphc){
  ggplot(bp_est, aes(y=species, x=breakpt, colour=data))+ #can add back shape
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data), size=3, position=ggstance::position_dodgev(height=0.4))+
    geom_linerange(aes(xmin = breakpt_se1, xmax = breakpt_se2, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"), drop=FALSE)+
    guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(5, "lines"))+ #Make more space between species
    scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    xlab("Metabolic Index Breakpoint Estimate")+
    ylab("Species")+
    scale_x_continuous(limits=c(0,25))
  # coord_cartesian(xlim=c(-20, 30))+
  #  geom_vline(xintercept=0, linetype="dashed")+
  # geom_vline(xintercept=1, linetype="dashed")
}
if(iphc){
  ggplot(bp_est, aes(y=species, x=breakpt, colour=data, shape=model))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data, shape=model), size=3, position=ggstance::position_dodgev(height=0.4))+
    geom_linerange(aes(xmin = breakpt_se1, xmax = breakpt_se2, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    #guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(2, "lines"))+ #Make more space between species
    scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    xlab("Metabolic Index Breakpoint Estimate")+
    ylab("Species")+
    coord_cartesian(xlim=c(0, 5))+
    geom_vline(xintercept=0, linetype="dashed")+
    geom_vline(xintercept=1, linetype="dashed")
}
ggsave(
  paste0("output/", output_folder, "/breakpt_est_truncated.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height =11,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

##pO2 crit: Calculate pO2 at a reference temperature and body size
#Reference temp (12 deg C)
ref_temp <- 15
#Body size at 1kg
body_size <- 2

for(i in 1:nrow(bp_est)){
  temp <- bp_est[i,]
  #find taxa from species
  taxa.2.use <- taxa$MI_Taxa[taxa$common_name==temp$species]
  #Calc invtemp
  kelvin = 273.15
  boltz = 0.000086173324
  tref <- 15
  invtemp.2.use <- (1 / boltz)  * ( 1 / (tref + 273.15) - 1 / (tref + 273.15))
  #Model
  model.2.use <- temp$model
  #calculate pO2 at a reference temperature and body size
  temp$est_o2 <- calc_po2_crit(invtemp.2.use,taxa.2.use,temp$breakpt,body_size, temp$model, fancy=F)
  temp$est_o2_low <-calc_po2_crit(invtemp.2.use,taxa.2.use,temp$breakpt_se1,body_size, temp$model, fancy=F)
  temp$est_o2_high<-calc_po2_crit(invtemp.2.use,taxa.2.use,temp$breakpt_se2,body_size, temp$model, fancy=F)
  temp$est_o2_lab <- calc_po2_crit(invtemp.2.use,taxa.2.use,1,body_size, temp$model, fancy=T)
  if(i==1){
    bp_est2 <- temp
  } else {
    bp_est2 <- bind_rows(bp_est2, temp)
  }
}

#Plot
##Plot line plot of pO2-crit estimates
if(!iphc){
  bp_est2$data <- factor(bp_est2$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
  labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
  names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")
  
  ggplot(bp_est2, aes(y=species, x=est_o2, colour=data))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data), size=2, position=ggstance::position_dodgev(height=0.4))+
    #Can add shape back
    geom_linerange(aes(xmin = est_o2_low, xmax = est_o2_high, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    # xlim(0,50)+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"), drop=FALSE)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(5, "lines"))+ #Make more space between species
    # scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    xlab(bquote(pO[2]~"(kPa)")) +
    ylab("Species")
}
if(iphc){
  ggplot(bp_est2, aes(y=species, x=est_o2, colour=data))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data, shape=model), size=2, position=ggstance::position_dodgev(height=0.4))+
    geom_linerange(aes(xmin = est_o2_low, xmax = est_o2_high, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    # xlim(0,50)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    #guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(2, "lines"))+ #Make more space between species
    scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    xlab(bquote(pO[2]~"(kPa)")) +
    ylab("Species")
}
if(both){
  bp_est2$data <- factor(bp_est2$data, levels=c("cc", "bc", "goa", "ebs", "coastwide", "cc _iphc", "bc _iphc", "goa _iphc", "ebs _iphc", "coastwide _iphc"))
  labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide", "British Columbia IPHC", "California Current IPHC", "Gulf of Alaska IPHC", "Eastern Bering Sea IPHC", "Coastwide IPHC")
  names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide", "cc _iphc", "bc _iphc", "goa _iphc", "ebs _iphc", "coastwide _iphc")
  #Add column for IPHC versus not
  bp_est2$type <- ifelse(grepl("iphc", bp_est2$data), "iphc included", "no iphc")
  
  ggplot(bp_est2, aes(y=species, x=est_o2, colour=data))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data, shape=type), size=2, position=ggstance::position_dodgev(height=0.4))+
    #Can add shape back
    geom_linerange(aes(xmin = est_o2_low, xmax = est_o2_high, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    # xlim(0,50)+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6", "#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(5, "lines"))+ #Make more space between species
    # scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    xlab(bquote(pO[2]~"(kPa)")) +
    ylab("Species")
}

ggsave(
  paste0("output/", output_folder, "/breakpt_est_o2_all.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height =11,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

#Only show reasonable ones
if(!iphc){
  ggplot(filter(bp_est2, est_o2>0), aes(y=species, x=est_o2, colour=data))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data), size=2, position=ggstance::position_dodgev(height=0.4))+
    geom_linerange(aes(xmin = est_o2_low, xmax = est_o2_high, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    # xlim(0,50)+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(5, "lines"))+ #Make more space between species
    scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    xlab(bquote(pO[2]~"(kPa)")) +
    ylab("Species")+
    coord_cartesian(xlim=c(0, 75))+
    scale_x_continuous(expand=c(0,0))
}
if(iphc){
  ggplot(bp_est2, aes(y=species, x=est_o2, colour=data))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data, shape=model), size=2, position=ggstance::position_dodgev(height=0.4))+
    geom_linerange(aes(xmin = est_o2_low, xmax = est_o2_high, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    # xlim(0,50)+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(2, "lines"))+ #Make more space between species
    scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    xlab(bquote(pO[2]~"(kPa)")) +
    ylab("Species")+
    coord_cartesian(xlim=c(0, 75))+
    scale_x_continuous(expand=c(0,0))
}

ggsave(
  paste0("output/", output_folder, "/breakpt_est_o2_truncated.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height =11,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)


#Add laboratory estimates (just single points)
##Add laboratory lines
#Only show reasonable ones
if(!iphc){
  ggplot(filter(bp_est2, est_o2>0),aes(y=species, x=est_o2, colour=data))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data), size=2, position=ggstance::position_dodgev(height=0.4))+
    geom_point(mapping=aes(x=est_o2_lab), colour="black", shape="triangle")+
    geom_linerange(aes(xmin = est_o2_low, xmax = est_o2_high, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    # xlim(0,50)+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    # guides(colour=guide_legend(nrow=1,byrow=TRUE),shape=guide_legend(nrow=1,byrow=TRUE))+
    guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(5, "lines"))+ #Make more space between species
    xlab(bquote(pO[2]~"(kPa)")) +
    ylab("Species")+
    scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    coord_cartesian(xlim=c(0, 75))+
    scale_x_continuous(expand=c(0,0))
}
if(iphc){
  ggplot(bp_est2, aes(y=species, x=est_o2, colour=data))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data, shape=model), size=2, position=ggstance::position_dodgev(height=0.4))+
    geom_point(mapping=aes(x=est_o2_lab), colour="black", shape="triangle")+
    geom_linerange(aes(xmin = est_o2_low, xmax = est_o2_high, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    # xlim(0,50)+
    # scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    # guides(colour=guide_legend(nrow=1,byrow=TRUE),shape=guide_legend(nrow=1,byrow=TRUE))+
    # guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(2, "lines"))+ #Make more space between species
    xlab(bquote(pO[2]~"(kPa)")) +
    ylab("Species")+
    scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    coord_cartesian(xlim=c(0, 75))+
    scale_x_continuous(expand=c(0,0))
}
ggsave(
  paste0("output/", output_folder, "/breakpt_est_o2_lab_ests.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height =11,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

#Add to dataset
dats50_12$data <- "laboratory phylogenetic imputation"
dats90_12$data <- "laboratory phylogenetic imputation"
dats50_12$est_o2_low <- dats50_12$min
dats50_12$est_o2_high <- dats50_12$max
dats90_12$est_o2_low <- dats90_12$min
dats90_12$est_o2_high <- dats90_12$max
#Filter for species included
dats50_12 <- filter(dats50_12, species %in% unique(bp_est2$species))
dats90_12 <- filter(dats90_12, species %in% unique(bp_est2$species))

bp_est5 <- bind_rows(bp_est2, dats50_12, dats90_12)

##Add ribbons
if(!iphc){
  ggplot(filter(bp_est5,est_o2>0), aes(y=species, x=est_o2, colour=data))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data), size=2, position=ggstance::position_dodgev(height=0.4))+
    geom_point(mapping=aes(x=est_o2_lab), colour="black", shape="triangle")+
    geom_linerange(filter(bp_est5, data!="laboratory phylogenetic imputation"),mapping=aes(xmin = est_o2_low, xmax = est_o2_high, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    geom_linerange(filter(bp_est5, data=="laboratory phylogenetic imputation"),mapping=aes(xmin = est_o2_low, xmax = est_o2_high), colour="black", alpha=0.5, linetype="dashed")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    # xlim(0,50)+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6", "black"))+
    # guides(colour=guide_legend(nrow=1,byrow=TRUE),shape=guide_legend(nrow=1,byrow=TRUE))+
    guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(5, "lines"))+ #Make more space between species
    xlab(bquote(pO[2]~"(kPa)")) +
    ylab("Species")+
    scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    #coord_cartesian(xlim=c(0, 75))+
    scale_x_continuous(expand=c(0,0))
}
if(iphc){
  ggplot(bp_est5, aes(y=species, x=est_o2, colour=data))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data, shape=model), size=2, position=ggstance::position_dodgev(height=0.4))+
    geom_point(mapping=aes(x=est_o2_lab), colour="black", shape="triangle")+
    geom_linerange(filter(bp_est5, data!="laboratory phylogenetic imputation"),mapping=aes(xmin = est_o2_low, xmax = est_o2_high, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    geom_linerange(filter(bp_est5, data=="laboratory phylogenetic imputation"),mapping=aes(xmin = est_o2_low, xmax = est_o2_high), colour="black", alpha=0.5, linetype="dashed")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    # xlim(0,50)+
    # scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6", "black"))+
    # guides(colour=guide_legend(nrow=1,byrow=TRUE),shape=guide_legend(nrow=1,byrow=TRUE))+
    # guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(2, "lines"))+ #Make more space between species
    xlab(bquote(pO[2]~"(kPa)")) +
    ylab("Species")+
    scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    coord_cartesian(xlim=c(0, 75))+
    scale_x_continuous(expand=c(0,0))
}

ggsave(
  paste0("output/", output_folder, "/breakpt_est_o2_lab_ests_ribbon.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height =11,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

#Truncated
if(!iphc){
  ggplot(filter(bp_est5, est_o2>0), aes(y=species, x=est_o2, colour=data))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data), size=2, position=ggstance::position_dodgev(height=0.4))+
    geom_point(mapping=aes(x=est_o2_lab), colour="black", shape="triangle")+
    geom_linerange(filter(bp_est5, data!="laboratory phylogenetic imputation"),mapping=aes(xmin = est_o2_low, xmax = est_o2_high, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    geom_linerange(filter(bp_est5, data=="laboratory phylogenetic imputation"),mapping=aes(xmin = est_o2_low, xmax = est_o2_high), colour="black", alpha=0.5, linetype="dashed")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    # xlim(0,50)+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6", "black"))+
    # guides(colour=guide_legend(nrow=1,byrow=TRUE),shape=guide_legend(nrow=1,byrow=TRUE))+
    guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(5, "lines"))+ #Make more space between species
    xlab(bquote(pO[2]~"(kPa)")) +
    ylab("Species")+
    scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    coord_cartesian(xlim=c(0, 10))+
    scale_x_continuous(expand=c(0,0))
}
if(iphc){
  ggplot(bp_est5, aes(y=species, x=est_o2, colour=data))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data, shape=model), size=2, position=ggstance::position_dodgev(height=0.4))+
    geom_point(mapping=aes(x=est_o2_lab), colour="black", shape="triangle")+
    geom_linerange(filter(bp_est5, data!="laboratory phylogenetic imputation"),mapping=aes(xmin = est_o2_low, xmax = est_o2_high, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    geom_linerange(filter(bp_est5, data=="laboratory phylogenetic imputation"),mapping=aes(xmin = est_o2_low, xmax = est_o2_high), colour="black", alpha=0.5, linetype="dashed")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    # xlim(0,50)+
    #scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6", "black"))+
    # guides(colour=guide_legend(nrow=1,byrow=TRUE),shape=guide_legend(nrow=1,byrow=TRUE))+
    # guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(2, "lines"))+ #Make more space between species
    xlab(bquote(pO[2]~"(kPa)")) +
    ylab("Species")+
    scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    coord_cartesian(xlim=c(0, 10))+
    scale_x_continuous(expand=c(0,0))
}
ggsave(
  paste("output/", output_folder, "/breakpt_est_o2_lab_ests_ribbon_zoomed.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height =11,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

#Save as CSV
#Save as csv
#Round
bp_est_save <- bp_est2
bp_est_save$breakpt<- round(as.numeric(bp_est_save$breakpt), digits=3)
bp_est_save$breakpt_se1<- round(as.numeric(bp_est_save$breakpt_se1), digits=3)
bp_est_save$breakpt_se2<- round(as.numeric(bp_est_save$breakpt_se1), digits=3)
bp_est_save$est_o2<- round(as.numeric(bp_est_save$est_o2), digits=3)
bp_est_save$est_o2_low<- round(as.numeric(bp_est_save$est_o2_low), digits=3)
bp_est_save$est_o2_high<- round(as.numeric(bp_est_save$est_o2_high), digits=3)
bp_est_save$est_o2_lab<- round(as.numeric(bp_est_save$est_o2_lab), digits=3)
write.csv(bp_est_save, file=paste0("output/", output_folder, "/breakpoint_est.csv"))
#-----------------------------------------------------------------------------------------------
#Sequence of temperatures
kelvin = 273.15
boltz = 0.000086173324
tref <- 15

t.range <- seq(min(dat$temperature_C), max(dat$temperature_C), length.out = 100)
t.range2 <-(1 / boltz)  * ( 1 / (t.range+ 273.15) - 1 / (tref + 273.15))

##Calculate pO2 at the range of temperatures and reference body size
bp.2.use <- filter(bp_est, breakpt>0)
for(i in 1:nrow(bp.2.use)){
  test <- bp.2.use[i,]
  #find taxa from species
  taxa.2.use <- taxa$MI_Taxa[taxa$common_name==test$species]
  #body size
  body_size <- 2
  #Model
  model.2.use <- test$model
  
  #calculate pO2 at a reference temperature and body size
  est_o2 <- calc_po2_crit(t.range2,taxa.2.use,test$breakpt,body_size, test$model, fancy=F)
  est_o2_se1 <- calc_po2_crit(t.range2,taxa.2.use,test$breakpt_se1,body_size, test$model, fancy=F)
  est_o2_se2 <- calc_po2_crit(t.range2,taxa.2.use,test$breakpt_se2,body_size, test$model, fancy=F)
  if(i==1){
    bp_est3 <- data.frame(species=test$species, data=test$data, model=test$model, breakpt=test$breakpt, breakpt_se1=test$breakpt_se1,breakpt_se2=test$breakpt_se2, invtemp=t.range2, temp=t.range, po2_crit=est_o2, po2_crit_se1=est_o2_se1, po2_crit_se2=est_o2_se2)
  } else {
    dat3 <- data.frame(species=test$species, data=test$data, model=test$model, breakpt=test$breakpt,  breakpt_se1=test$breakpt_se1,breakpt_se2=test$breakpt_se2, invtemp=t.range2, temp=t.range, po2_crit=est_o2,  po2_crit_se1=est_o2_se1, po2_crit_se2=est_o2_se2)
    bp_est3 <- bind_rows(bp_est3, dat3)
  }
}

##Combine all datasets used
for(i in 1:nrow(bp.2.use)) {
  this_dat <- bp.2.use[i,]
  this_species <- this_dat$species
  dat.2.est <- this_dat$data
  print(this_species)
  best_model <- this_dat$model
  #Pull the data file
  this_datframe <- try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", dat.2.est, "_dat.rds")))
  #Combine data
  this_datframe$model <- paste(best_model)
  this_datframe$data <- paste(dat.2.est)
  this_datframe$species <- paste(this_species)
  if(i==1){
    dat2 <- this_datframe
  } else {
    dat2 <- bind_rows(dat2, this_datframe)
  }
}

if(!iphc){
  bp_est3$data <- factor(bp_est3$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
  dat2$data <- factor(dat2$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
  labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
  names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")
}

##Plot temp vs o2 and po2 crits for all species and models
theme_set(theme_bw(base_size = 30))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             theme(strip.background = element_blank()))

#Add labels
dat2$label <- paste(dat2$species, dat2$data, sep=" ")
bp_est3$label <- paste(bp_est3$species, bp_est3$data, sep=" ")

ggplot(data = dat2,aes(x = temperature_C)) +  
  facet_wrap("label", ncol=4)+
  geom_point(aes(y=po2,colour=depth), alpha=0.5, size=0.25)+
  geom_line(data=bp_est3,mapping=aes(x=temp, y=po2_crit), colour="black")+
  geom_ribbon(data=bp_est3, mapping=aes(x=temp, ymin=(po2_crit_se1), ymax=(po2_crit_se2)), alpha=0.5, fill="lightgrey")+
  xlab("Temperature (C)") +
  ylab(bquote(pO[2]~"(kPa)")) +
  # scale_x_continuous(expand = c(0,0), limits = c(0, 15) ) +
  # scale_y_continuous(expand = c(0,0), limits = c(0,30) )+
  theme(legend.position="right")+
  coord_cartesian(ylim=c(0, 20))+
  scale_x_continuous(expand=c(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=20))+
  # geom_text(aes(label = labels, y=po2, x=temperature_C), data = labels, vjust = 1) +
  scale_colour_viridis(option="mako", guide=guide_colourbar(reverse = TRUE), direction=-1)+
  theme(panel.spacing = unit(0.2, "lines"))

ggsave(
  paste0("output/", output_folder, "/temp_o2_truncated.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 16,
  height =22,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

#Plot separately
for(i in 1:nrow(bp.2.use)){
  test <- bp.2.use[i,]
  this_species <- test$species
  this_dat <- test$data
  dat.2.plot <- filter(dat2, common_name==this_species & data==this_dat)
  dat.2.plot.too <- filter(bp_est3, species==this_species & data==this_dat)
  ggplot(data = dat.2.plot,aes(x = temperature_C)) +
    geom_point(aes(y=po2,colour=depth), alpha=1, size=0.25)+
    geom_line(data=dat.2.plot.too,mapping=aes(x=temp, y=po2_crit), colour="black")+
    geom_ribbon(data=dat.2.plot.too, mapping=aes(x=temp, ymin=(po2_crit_se1), ymax=(po2_crit_se2)), alpha=0.3, fill="lightgrey")+
    xlab("Temperature (C)") +
    ylab(bquote(pO[2]~"(kPa)")) +
    # scale_x_continuous(expand = c(0,0), limits = c(0, 15) ) +
    # scale_y_continuous(expand = c(0,0), limits = c(0,30) )+
    theme(legend.position=c(0.8,0.8))+
    #coord_cartesian(ylim=c(0, 40))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank())+
    scale_colour_viridis(option="mako", guide=guide_colourbar(reverse = TRUE), direction=-1)+
    ggtitle(paste(this_species, " ", this_region))
  
  ggsave(
    paste0("output/", output_folder, "/", this_species, "_", this_dat, ".png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height =8.5,
    units = c("in"),
    bg="white",
    dpi = 600,
    limitsize = TRUE
  )
  
}

#Plot pO2 crit from the different models all on one plot
theme_set(theme_bw(base_size = 15))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             theme(strip.background = element_blank()))
if(!iphc){
  ggplot(bp_est3, aes(x=temp, y=po2_crit))+
    facet_wrap("species", scales="free_y")+
    #geom_ribbon(data=bp_est2, mapping=aes(x=temp, ymin=(po2_crit_se1), ymax=(po2_crit_se2), fill=data), alpha=0.3)+
    geom_line(aes(colour=data), size=1)+
    theme(legend.position="top")+
    theme(legend.title=element_blank(), strip.background = element_blank())+
    theme(text=element_text(size=15))+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    xlab("Temperature (C)") +
    ylab(bquote(pO[2]~"(kPa)"))
  #coord_cartesian(xlim=c(0, 15), ylim=c(0, 30))
}
if(iphc){
  ggplot(bp_est3, aes(x=temp, y=po2_crit))+
    facet_wrap("species", scales="free_y")+
    #geom_ribbon(data=bp_est2, mapping=aes(x=temp, ymin=(po2_crit_se1), ymax=(po2_crit_se2), fill=data), alpha=0.3)+
    geom_line(aes(colour=data), size=1)+
    theme(legend.position="top")+
    theme(legend.title=element_blank(), strip.background = element_blank())+
    theme(text=element_text(size=15))+
    #  scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    xlab("Temperature (C)") +
    ylab(bquote(pO[2]~"(kPa)"))
  #coord_cartesian(xlim=c(0, 15), ylim=c(0, 30))
}

ggsave(
  paste("output", output_folder, "/po2_crit_all_no_se.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height =8.5,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

if(!iphc){
  ggplot(bp_est3, aes(x=temp, y=po2_crit))+
    facet_wrap("species", scales="free_y")+
    geom_ribbon(data=bp_est3, mapping=aes(x=temp, ymin=po2_crit_se1, ymax=po2_crit_se2, fill=data),alpha=0.5)+
    geom_line(aes(colour=data), size=1)+
    theme(legend.position="top")+
    theme(legend.title=element_blank(), strip.background = element_blank())+
    theme(text=element_text(size=15))+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    scale_fill_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    xlab("Temperature (C)") +
    ylab(bquote(pO[2]~"(kPa)"))
  #coord_cartesian(xlim=c(0, 15), ylim=c(0, 40))
}

if(iphc){
  ggplot(bp_est3, aes(x=temp, y=po2_crit))+
    facet_wrap("species", scales="free_y")+
    geom_ribbon(data=bp_est3, mapping=aes(x=temp, ymin=po2_crit_se1, ymax=po2_crit_se2, fill=data),alpha=0.5)+
    geom_line(aes(colour=data), size=1)+
    theme(legend.position="top")+
    theme(legend.title=element_blank(), strip.background = element_blank())+
    theme(text=element_text(size=15))+
    #  scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    #scale_fill_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    xlab("Temperature (C)") +
    ylab(bquote(pO[2]~"(kPa)"))
  #coord_cartesian(xlim=c(0, 15), ylim=c(0, 40))
}
ggsave(
  paste0("output/", output_folder, "/po2_crit_all_with_se.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height =8.5,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

if(iphc){
  #Just IPHC data
  ggplot(filter(bp_est3, grepl("iphc", bp_est3$data)), aes(x=temp, y=po2_crit))+
    facet_wrap("species", scales="free_y")+
    geom_ribbon(data=filter(bp_est3, grepl("iphc", bp_est3$data)), mapping=aes(x=temp, ymin=po2_crit_se1, ymax=po2_crit_se2, fill=data),alpha=0.5)+
    geom_line(aes(colour=data), size=1)+
    theme(legend.position="top")+
    theme(legend.title=element_blank(), strip.background = element_blank())+
    theme(text=element_text(size=15))+
    #  scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    #scale_fill_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    xlab("Temperature (C)") +
    ylab(bquote(pO[2]~"(kPa)"))
  #coord_cartesian(xlim=c(0, 15), ylim=c(0, 40))
  
  ggsave(
    paste("output/plots/po2_crit_all_with_se_iphc_only.png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 10,
    height =8.5,
    units = c("in"),
    bg="white",
    dpi = 600,
    limitsize = TRUE
  )
  #Separate into different plots
  ggplot(filter(bp_est3, grepl("iphc", bp_est3$data)), aes(x=temp, y=po2_crit))+
    facet_wrap("label", scales="free_y")+
    geom_ribbon(data=filter(bp_est3, grepl("iphc", bp_est3$data)), mapping=aes(x=temp, ymin=po2_crit_se1, ymax=po2_crit_se2, fill=data),alpha=0.5)+
    geom_line(aes(colour=data), size=1)+
    theme(legend.position="top")+
    theme(legend.title=element_blank(), strip.background = element_blank())+
    theme(text=element_text(size=13))+
    #  scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    #scale_fill_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    xlab("Temperature (C)") +
    ylab(bquote(pO[2]~"(kPa)"))
  #coord_cartesian(xlim=c(0, 15), ylim=c(0, 40))
  
  ggsave(
    paste("output/plots/po2_crit_all_with_se_iphc_only_sep.png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 10,
    height =8.5,
    units = c("in"),
    bg="white",
    dpi = 600,
    limitsize = TRUE
  )
  
}
####Plot w/ Tim's############################
#Filter for species included
dats50a <- filter(dats50, species %in% unique(bp_est3$species))
dats90a <- filter(dats90, species %in% unique(bp_est3$species))

if(!iphc){
  ggplot(bp_est3, aes(x=temp, y=po2_crit))+
    facet_wrap("species", scales="free_y")+
    geom_ribbon(data=bp_est3, mapping=aes(x=temp, ymin=po2_crit_se1, ymax=po2_crit_se2, fill=data),alpha=0.2)+
    geom_line(aes(colour=data), size=1)+
    geom_polygon(data = dats50a, mapping=aes(x = temp, y = ys), fill = "grey1", color = NA, alpha = 0.5) + 
    geom_polygon(data = dats90a, mapping=aes(x = temp, y = ys), fill = "lightgrey", color = NA, alpha = 0.5) +
    theme(legend.position="top")+
    theme(legend.title=element_blank(), strip.background = element_blank())+
    theme(text=element_text(size=15))+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    scale_fill_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    xlab("Temperature (C)") +
    ylab(bquote(pO[2]~"(kPa)"))
  #coord_cartesian(xlim=c(0, 15), ylim=c(0, 40))+
}
if(iphc){
  ggplot(bp_est3, aes(x=temp, y=po2_crit))+
    facet_wrap("species", scales="free_y")+
    geom_ribbon(data=bp_est3, mapping=aes(x=temp, ymin=po2_crit_se1, ymax=po2_crit_se2, fill=data),alpha=0.2)+
    geom_line(aes(colour=data), size=1)+
    geom_polygon(data = dats50a, mapping=aes(x = temp, y = ys), fill = "grey1", color = NA, alpha = 0.5) + 
    geom_polygon(data = dats90a, mapping=aes(x = temp, y = ys), fill = "lightgrey", color = NA, alpha = 0.5) +
    theme(legend.position="top")+
    theme(legend.title=element_blank(), strip.background = element_blank())+
    theme(text=element_text(size=15))+
    # scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    # scale_fill_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    xlab("Temperature (C)") +
    ylab(bquote(pO[2]~"(kPa)"))
}

ggsave(
  paste0("output/", output_folder, "/po2_crit_with_lab.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height =8.5,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

#Without SE
if(!iphc){
  ggplot(bp_est3, aes(x=temp, y=po2_crit))+
    facet_wrap("species", scales="free_y")+
    #geom_ribbon(data=bp_est3, mapping=aes(x=temp, ymin=po2_crit_se1, ymax=po2_crit_se2, fill=data),alpha=0.5)+
    geom_polygon(data = dats50a, mapping=aes(x = temp, y = ys), fill = "grey1", color = NA, alpha = 0.5) + 
    geom_polygon(data = dats90a, mapping=aes(x = temp, y = ys), fill = "lightgrey", color = NA, alpha = 0.5) +
    geom_line(aes(colour=data), size=1)+
    theme(legend.position="top")+
    theme(legend.title=element_blank(), strip.background = element_blank())+
    theme(text=element_text(size=15))+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    scale_fill_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    xlab("Temperature (C)") +
    ylab(bquote(pO[2]~"(kPa)"))
  #coord_cartesian(xlim=c(0, 15), ylim=c(0, 40))
}
if(iphc){
  ggplot(bp_est3, aes(x=temp, y=po2_crit))+
    facet_wrap("species", scales="free_y")+
    #geom_ribbon(data=bp_est3, mapping=aes(x=temp, ymin=po2_crit_se1, ymax=po2_crit_se2, fill=data),alpha=0.5)+
    geom_polygon(data = dats50a, mapping=aes(x = temp, y = ys), fill = "grey1", color = NA, alpha = 0.5) + 
    geom_polygon(data = dats90a, mapping=aes(x = temp, y = ys), fill = "lightgrey", color = NA, alpha = 0.5) +
    geom_line(aes(colour=data), size=1)+
    theme(legend.position="top")+
    theme(legend.title=element_blank(), strip.background = element_blank())+
    theme(text=element_text(size=15))+
    # scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    # scale_fill_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    xlab("Temperature (C)") +
    ylab(bquote(pO[2]~"(kPa)"))
  #coord_cartesian(xlim=c(0, 15), ylim=c(0, 40))
}
ggsave(
  paste0("output/", output_folder, "/po2_crit_with_lab2.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height =8.5,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

#Plot each separately
if(!iphc){
  for(i in 1:nrow(bp_est3)){
    test <- bp_est3[i,]
    this_species <- test$species
    this_dat <- test$data
    p <- ggplot(bp_est3, aes(x=temp, y=po2_crit))+
      facet_wrap("species", scales="free_y")+
      #geom_ribbon(data=bp_est3, mapping=aes(x=temp, ymin=po2_crit_se1, ymax=po2_crit_se2, fill=data),alpha=0.5)+
      geom_polygon(data = dats50a, mapping=aes(x = temp, y = ys), fill = "grey1", color = NA, alpha = 0.5) + 
      geom_polygon(data = dats90a, mapping=aes(x = temp, y = ys), fill = "lightgrey", color = NA, alpha = 0.5) +
      geom_line(aes(colour=data), size=1)+
      theme(legend.position="top")+
      theme(legend.title=element_blank(), strip.background = element_blank())+
      theme(text=element_text(size=15))+
      scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
      scale_fill_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
      xlab("Temperature (C)") +
      ylab(bquote(pO[2]~"(kPa)"))
    #coord_cartesian(xlim=c(0, 15), ylim=c(0, 40))
    ggtitle(paste(this_species, " ", this_region))
    
    ggsave(
      paste0("output/", output_folder, "/po2_crit_lab", this_species, "_", this_dat, ".png"),
      plot = last_plot(),
      device = NULL,
      path = NULL,
      scale = 1,
      width = 8.5,
      height =8.5,
      units = c("in"),
      bg="white",
      dpi = 600,
      limitsize = TRUE
    )
  }
}

##Plot temp and depth, and where was O2 suitable and not?
##Calculate pO2 breakpoint for each sampling event temp and species (reference body size)
##Do this just for each region, not for the coastwide models
bp_est4 <- filter(bp.2.use, data!="coastwide")
for(i in 1:nrow(bp_est4)){
  dat.2.est <- bp_est4[i,]
  this_species = dat.2.est$species
  print(this_species)
  this_dat <- dat.2.est$data
  this_model <- dat.2.est$model
  best_model <- this_model
  #Pull the data file
  this_datframe <-try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
  #Calculate mean of this_datframe$mi based on best model
  mean_mi <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") mean(this_datframe$mi1) else if(best_model=="model4"|best_model=="model10"|best_model=="model14") mean(this_datframe$mi2) else mean(this_datframe$mi3)
  sd_mi <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") sd(this_datframe$mi1) else if(best_model=="model4"|best_model=="model10"|best_model=="model14") sd(this_datframe$mi2) else sd(this_datframe$mi3)
  #Pull breakpoint and slope
  #find taxa from species
  taxa.2.use <- taxa$MI_Taxa[taxa$common_name==this_species]
  #body size
  body_size <- 1
  #Model
  model.2.use <- best_model
  test <- this_datframe
  test$model <- paste(best_model)
  model<- paste(best_model)
  test$data <- paste(this_dat)
  test$species <- paste(this_species)
  thresh_est <- dat.2.est$breakpt
  thresh_se <- dat.2.est$breakpt_se                                               
  thresh_se1 <- thresh_est-thresh_se
  thresh_se2 <- thresh_est+thresh_se
  #calculate pO2 at a reference temperature and body size
  #calc_po2_crit function across test in an apply function
  invtemp <- test$invtemp
  test$est_o2 <- unlist(lapply(invtemp, calc_po2_crit, taxa.2.use,thresh_est,body_size, model, fancy=T))
  if(!is.null(thresh_se)){
    test$est_o2_se1 <- unlist(lapply(invtemp, calc_po2_crit, taxa.2.use,thresh_se1,body_size, model, fancy=T))
    test$est_o2_se2 <- unlist(lapply(invtemp, calc_po2_crit, taxa.2.use,thresh_se2,body_size, model, fancy=T))
  }
  # test$percentile <- (test$po2-test$est_o2)/test$est_o2_se
  test$unsuitable <- ifelse(test$po2<test$est_o2, "unsuitable", "suitable")
  if(!is.null(thresh_se)){
    test$unsuitable_low <- ifelse(test$po2<test$est_o2_se1, "unsuitable", "suitable")
    test$unsuitable_high <- ifelse(test$po2<test$est_o2_se2,"unsuitable", "suitable")
  }
  #Sum across row to get total suitable
  #test$suitable_total <- rowSums(test[,c("unsuitable", "unsuitable_low", "unsuitable_high")])
  if(i==1){
    dats <- test
  } else {
    dats <- bind_rows(dats, test)
  }
}

#Plot CPUE by depth vs temp, and faded with pO2 below the threshold--
##Plot all together
#Make labels
dats$label <- paste(dats$species, dats$data, sep=" ")

ggplot(data = dats,aes(x = temperature_C, y=depth)) +  
  geom_point(aes(colour=as.factor(unsuitable)), size=0.25, alpha=0.2)+
  facet_wrap("label", ncol=4, scales="free_y")+
  xlab("Temperature (C)") +
  ylab("Depth (m)") +
  scale_y_reverse()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=12))+
  theme(legend.position=c(0.9,0.05))+
  # geom_text(aes(label = labels, y=-Inf, x=temperature_C), data = labels, vjust=1)+
  scale_colour_manual(values=c("blue3", "darkorange"), labels=c("Above pO2 crit", "Below pO2 crit"))+
  guides(color = guide_legend(title="",override.aes = list(size = 7, alpha=1)))

ggsave(
  paste0("output/", output_folder, "/unsuitable_suitable_depth.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height =11,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

if(output_folder=="presence_absence"){
  dats$colour <- case_when(dats$presence==1&dats$unsuitable=="suitable" ~ "True Presence",
                           dats$presence==1&dats$unsuitable=="unsuitable" ~ "False Presence",
                           dats$presence==0&dats$unsuitable=="suitable" ~ "False Absence",
                           dats$presence==0&dats$unsuitable=="unsuitable" ~ "True Absence")
  ggplot(data = dats,aes(x = temperature_C, y=depth)) +  
    geom_point(aes(colour=as.factor(colour)), size=0.25, alpha=0.1)+
    facet_wrap("label", ncol=4, scales="free_y")+
    xlab("Temperature (C)") +
    ylab("Depth (m)") +
    scale_y_reverse()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position=c(0.9,0.05))+
    # geom_text(aes(label = labels, y=-Inf, x=temperature_C), data = labels, vjust=1)+
    scale_colour_manual(values=c("salmon1", "salmon3", "darkseagreen2", "springgreen4"))+
    guides(color = guide_legend(title="",override.aes = list(size = 7, alpha=1)))
  
  ggsave(
    paste0("output/", output_folder, "/unsuitable_suitable_depth.png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height =11,
    units = c("in"),
    bg="white",
    dpi = 600,
    limitsize = TRUE
  )
}
##Shaded by CPUE
ggplot(data = filter(dats, unsuitable==0),aes(x = temperature_C, y=depth)) +  
  geom_point(aes(colour=log(catch_weight+1)), size=0.25, alpha=0.4)+
  scale_colour_distiller(name="Above pO2 crit", type="seq",palette="Blues", direction=1)+
  new_scale_colour() +
  geom_point(data=filter(dats, unsuitable==1), mapping=aes(x = temperature_C, y=depth, colour=log(catch_weight+1)), size=0.25, alpha=0.4)+
  scale_colour_distiller(name="Below pO2 crit", type="seq",palette="Oranges", direction=1)+
  facet_wrap("label", ncol=4, scales="free_y")+
  xlab("Temperature (C)") +
  ylab("Depth (m)") +
  scale_y_reverse()+
  #  coord_cartesian(ylim=c(1500, -100))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=12))+
  theme(legend.position="none")+
  theme(legend.text=element_text(size=10), legend.title=element_text(size=10))
#geom_text(aes(label = labels, y=-Inf, x=temperature_C), data = labels, vjust=1)

ggsave(
  paste("output/plots/unsuitable_suitable_depth_cpue.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height =8.5,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

##save separate for each species, colored by cpue
for(i in 1:nrow(bp_est4)){
  test <- bp_est4[i,]
  this_species <- test$species
  this_dat <- test$data
  dat.2.plot <- filter(dats, common_name==this_species & data==this_dat)
  p <- ggplot(data = filter(dat.2.plot, unsuitable==0),aes(x = temperature_C, y=depth)) +  
    geom_point(aes(colour=log(catch_weight+1)), size=0.25, alpha=0.4)+
    scale_colour_distiller(name="Above pO2 crit", type="seq",palette="Blues", direction=1)+
    new_scale_colour() +
    geom_point(data=filter(dats, unsuitable==1), mapping=aes(x = temperature_C, y=depth, colour=log(catch_weight+1)), size=0.25, alpha=0.4)+
    scale_colour_distiller(name="Below pO2 crit", type="seq",palette="Oranges", direction=1)+
    xlab("Temperature (C)") +
    ylab("Depth (m)") +
    scale_y_reverse()+
    #  coord_cartesian(ylim=c(1500, -100))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          legend.text=element_text(size=10), 
          legend.position="top")+
    ggtitle(paste(this_species, " ", this_region))
  
  ggsave(
    paste("output/plots/po2_obs_cpue/po2_obs_cpue_", this_species, "_", this_dat, ".png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height =8.5,
    units = c("in"),
    bg="white",
    dpi = 600,
    limitsize = TRUE
  )
}

##save separate for each species, with faded by CPUE
for(i in 1:nrow(bp_est4)){
  test <- bp_est4[i,]
  this_species <- test$species
  this_dat <- test$data
  dat.2.plot <- filter(dats, common_name==this_species & data==this_dat)
  p <- ggplot(data = filter(dat.2.plot, unsuitable==0),aes(x = temperature_C, y=depth)) +  
    geom_point(aes(colour=log(catch_weight+1)), size=0.25, alpha=0.4)+
    scale_colour_distiller(name="Above pO2 crit", type="seq",palette="Blues", direction=1)+
    new_scale_colour() +
    geom_point(data=filter(dats, unsuitable==1), mapping=aes(x = temperature_C, y=depth, colour=log(catch_weight+1)), size=0.25, alpha=0.4)+
    scale_colour_distiller(name="Below pO2 crit", type="seq",palette="Oranges", direction=1)+
    xlab("Temperature (C)") +
    ylab("Depth (m)") +
    scale_y_reverse()+
    #  coord_cartesian(ylim=c(1500, -100))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          legend.text=element_text(size=10), 
          legend.position="top")+
    ggtitle(paste(this_species, " ", this_region))
  
  ggsave(
    paste("output/plots/po2_obs_cpue2/po2_obs_cpue_", this_species, "_", this_dat, ".png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height =8.5,
    units = c("in"),
    bg="white",
    dpi = 600,
    limitsize = TRUE
  )
}

###Plot with the SE
ggplot(data = dats,aes(x = temperature_C, y=depth)) +  
  geom_point(aes(colour=as.factor(suitable_total)), size=0.25, alpha=0.2)+
  facet_wrap("label", ncol=4, scales="free_y")+
  xlab("Temperature (C)") +
  ylab("Depth (m)") +
  scale_y_reverse()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=12),
        legend.text=element_text(size=12))+
  theme(legend.position=c(0.9,0.1))+
  # geom_text(aes(label = labels, y=-Inf, x=temperature_C), data = labels, vjust=1)+
  scale_colour_manual(values=c("lightblue", "orange", "darkorange2"), labels=c("Above pO2 crit", "Below pO2 crit+SE", "Below pO2 crit"))+
  guides(color = guide_legend(title="",override.aes = list(size = 7, alpha=1)))

ggsave(
  paste("output/plots/po2_obs_range.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height =8.5,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

###Map#########
map_data <- rnaturalearth::ne_countries(scale = "large",
                                        returnclass = "sf",
                                        continent = "North America")

us_coast_proj <- sf::st_transform(map_data, crs = 32610)

###Map of data available
species <- unique(dats$common_name)
dats$unsuitable <- factor(dats$unsuitable, levels=c("suitable","unsuitable"))

for(i in 1:length(species)){
  species2plot <- species[i]
  dat2plot <- filter(dats, common_name==species2plot)
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(dat2plot, mapping=aes(x=X*1000, y=Y*1000,colour=unsuitable), size=0.1)+
    xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
    ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
    facet_wrap("year", ncol=5)+
    theme_minimal(base_size=12)+
    xlab("Longitude")+
    ylab("Latitude")+
    ggtitle(paste(unique(dat2plot$common_name)))+
    theme(axis.text.x=element_blank())+
    scale_colour_manual(values=c("lightblue", "orange3"), labels=c("Above pO2 crit", "Below pO2 crit"), drop=F)+
    guides(color = guide_legend(title="",override.aes = list(size = 7, alpha=1)))+
    theme(legend.position="none")
  
  
  ggsave(
    paste("output/plots/habitat_mapping/map_",species[i],".png"),
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

#Plot residuals?
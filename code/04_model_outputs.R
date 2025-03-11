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

set.seed(9876)

setwd("~/Dropbox/GitHub/wsg-choke-species")

#Load functions
source("code/helper_funs.R")

##Create sequence of metabolic index values for marginal effects, so same for all
#Load data
files <- list.files(path = "data/processed_data/fish", pattern = ".rds", full.names=T)
dat <- map(files,readRDS)
dat <- bind_rows(dat)

#Remove NAs (and remove IPHC by removing cpue)
dat <- dat  %>%
  drop_na(depth,year, mi1,mi2,mi3, X, Y, cpue_kg_km2)
#Sequence of MI values
#Remove NAs (and remove IPHC by removing cpue)
dat <- dat  %>%
  drop_na(depth,year, mi1,mi2,mi3, X, Y, cpue_kg_km2)
mi1_pred = seq(min(dat$mi1), max(dat$mi1), length.out = 300)
mi2_pred = seq(min(dat$mi2), max(dat$mi2), length.out = 300)
mi3_pred = seq(min(dat$mi3), max(dat$mi3), length.out = 300)


#Load AIC table for model output
#Data type comparison
data_type <- T
region_comp <- F
if(data_type){
aic <- as.data.frame(read.csv("output/data_type/aic_table_data_type_priors_goodonly.csv"))
output_folder <- "data_type"
}
if(region_comp){
aic <- as.data.frame(read.csv("output/region_comp/aic_table_region_comp_priors_goodonly.csv"))
output_folder <- "region_comp"
}

species <- unique(aic$species)

##Loop to pull each set of comparisons that converged, pull the best-fitting model if includes MI, predict to a marginal effect prediction grid, and make a dataframe with them
#Create dataframe
marg_effects_preds = as.data.frame(matrix(NA, 1, 18))
for(i in 1:length(species)) {
  this_species = species[i]
  print(this_species)
  this_aic <- as.data.frame(filter(aic, species==this_species))
  dat_names <- unique(this_aic$data.type)
  for(h in 1:length(dat_names)){
    this_dat <- dat_names[h]
    this_data <- as.data.frame(filter(this_aic, data.type==this_dat))
    #Identify if MI is best-fitting model and pull that model
    mi_models <- this_data[c("model3", "model4", "model5")]
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
      mi_best <- if(best_model=="model3") mi1_pred else if(best_model=="model4") mi2_pred else mi3_pred
      pred_year <- unique(this_datframe$year)[2]
      nd_po2 <- data.frame(mi1_s = mi1_s_pred,
                           mi2_s =mi2_s_pred,
                           mi3_2 =mi3_s_pred,
                           mi1=mi1_pred,
                           mi2=mi2_pred,
                           mi3=mi3_pred,
                           mi_best=mi_best,
                temp_scaled = 0,
                temp_scaled2 = 0,
                log_depth_scaled = 0,
                log_depth_scaled2 = 0,
                year = pred_year,
                survey="wcbts",
                region="cc")
          #nd_po2 <- convert_class(nd_po2)
          p1 <- predict(fit, newdata = nd_po2, se_fit = TRUE, re_form = NA)
          p1$model <- paste(best_model)
          p1$data <- paste(this_dat)
          p1$species <- paste(this_species)
          p1 <-  p1 %>%
            mutate(est_sc= exp(est)/max(exp(est), na.rm=T),
                   est_se_sc1 = est_sc-(exp(est_se)/max(exp(est), na.rm=T)),
                   est_se_sc2 = est_sc+(exp(est_se)/max(exp(est), na.rm=T)))
                    # (exp(est)-exp(est_se))/max(exp(est), na.rm=T),
                   #est_se_sc2 = (exp(est)+exp(est_se))/max(exp(est), na.rm=T))
        #  mutate(est_sc= (exp(est)/max((exp(est)+exp(est_se)), na.rm=T)),
         #        est_se_sc1 = (exp(est)-exp(est_se))/max((exp(est)+exp(est_se)), na.rm=T),
           #      est_se_sc2 = (exp(est)+exp(est_se))/max((exp(est)+exp(est_se)), na.rm=T))
          marg_effects_preds <- bind_rows(marg_effects_preds, p1)
    }
  }
}

#Remove first row
marg_effects_preds <- marg_effects_preds[-1,]
#Remove first 18 columns
marg_effects_preds <- marg_effects_preds[,19:ncol(marg_effects_preds)]
#Remove space in species
marg_effects_preds <- as.data.frame(marg_effects_preds)

##Plot marginal effects
ggplot(marg_effects_preds, aes(mi_best, y=est_sc))+
  facet_wrap("species", scales="free_y")+
  geom_line(aes(colour=data))+
  geom_ribbon(aes(ymin = est_se_sc1, ymax = est_se_sc2, fill=data), alpha=0.4)+
  #xlim(0,200)+
  labs(x = bquote('Metabolic Index'), y = bquote('Population Density'~(kg~km^-2)))+
  theme_minimal()+
  theme(legend.position="top")+
  theme(text=element_text(size=15))



###Alternative
marg_effects_preds = as.data.frame(matrix(NA, 1, 18))
for(i in 1:length(species)) {
  this_species = species[i]
  print(this_species)
  this_aic <- as.data.frame(filter(aic, species==this_species))
  dat_names <- unique(this_aic$data.type)
  for(h in 1:length(dat_names)){
    this_dat <- dat_names[h]
    this_data <- as.data.frame(filter(this_aic, data.type==this_dat))
    #Identify if MI is best-fitting model and pull that model
    mi_models <- this_data[c("model3", "model4", "model5")]
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
      nd_po2 <- data.frame(mi_s=mi_best_s)
      
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
      
      nd_po2$est <- sapply(nd_po2$o2_s,breakpoint_calc, slope$estimate,thresh$estimate)
      marg_effects_preds <- bind_rows(marg_effects_preds, nd_po2)
    }
  }
}




##Plot historical temp & o2 with metabolic index from model predictions and from Tim's

##Make model predictions from data

##Plot historical hotspots in limited density (marginal effects) (or O2 below )

##Plot residuals
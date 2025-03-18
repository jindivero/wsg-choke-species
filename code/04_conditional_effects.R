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
data_type <- F
region_comp <- T
if(data_type){
aic <- as.data.frame(read.csv("output/data_type/aic_table_data_type_priors_goodonly.csv"))
output_folder <- "data_type"
}
if(region_comp){
aic <- as.data.frame(read.csv("output/region_comp/aic_table_region_comp_priors_goodonly.csv"))
output_folder <- "region_comp"
}

#Species to run
species <- unique(aic$species)

#Taxa lookup
taxa <- read_excel("data/species_table.xlsx")
taxa$MI_Taxa <- tolower(taxa$MI_Taxa)
taxa$common_name <- tolower(taxa$common_name)

##Loop to pull each set of comparisons that converged, pull the best-fitting model if includes MI, predict to a marginal effect prediction grid, and make a dataframe with them
#Create dataframe
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
      mi_best_s <- if(best_model=="model3") mi1_s_pred else if(best_model=="model4") mi2_s_pred else mi3_s_pred
      pred_year <- unique(this_datframe$year)[1]
      log_depth_scaled_mean <- mean(this_datframe$log_depth_scaled)
      log_depth_scaled_mean2 <- log_depth_scaled_mean^2
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
                 est_se_sc1 = (exp(est)-exp(est_se))/max(exp(est), na.rm=T),
                 est_se_sc2 = (exp(est)+exp(est_se))/max(exp(est), na.rm=T))
          cond_effects_preds <- bind_rows(cond_effects_preds, p1)
    }
  }
}

##Clean up dataframe
#Remove first row
cond_effects_preds <- cond_effects_preds[-1,]
#Remove first 18 columns
cond_effects_preds <- cond_effects_preds[,19:ncol(cond_effects_preds)]
cond_effects_preds$data <- factor(cond_effects_preds$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")

##Plot conditional effects
ggplot(cond_effects_preds, aes(mi_best, y=est_sc))+
  facet_wrap("species")+
  geom_line(aes(colour=data))+
  #geom_ribbon(aes(ymin = est_se_sc1, ymax = est_se_sc2, fill=data), alpha=0.4)+
  scale_x_continuous(limits=c(0,15, by=5))+
  labs(x = bquote('Metabolic Index'), y = bquote('Conditional Effect Population Density'~(kg~km^-2)))+
  theme_minimal()+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(labels=labs, values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
  theme(legend.key.height = unit(2, "lines"))+
  theme(panel.spacing = unit(1, "lines"))

ggsave(
  paste("output/plots/cond_effects_region_mi_scaled_no_se.png"),
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
  theme_minimal()+
  theme(legend.position="top")+
  theme(text=element_text(size=15))

#Not good

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
  theme_minimal()+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(labels=labs, values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
  theme(legend.key.height = unit(2, "lines"))+
  theme(panel.spacing = unit(1, "lines"))

#Breakpoint estimates of MI into line plots
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
      #Calculate mean of this_datframe$mi based on best model
      mean_mi <- if(best_model=="model3") mean(this_datframe$mi1) else if(best_model=="model4") mean(this_datframe$mi2) else mean(this_datframe$mi3)
      sd_mi <-   if(best_model=="model3") sd(this_datframe$mi1) else if(best_model=="model4") sd(this_datframe$mi2) else sd(this_datframe$mi3)
      #Pull breakpoint and slope
      pars <- as.data.frame(tidy(fit, effects="fixed", conf.int=T))
      slope <- filter(pars, grepl("slope", term))
      thresh <- filter(pars,grepl("breakpt", term))
      thresh$est <- (thresh$est*sd_mi)+mean_mi
      thresh$std.error <- (thresh$std.error*sd_mi)+mean_mi 
      if(i==1 & h==1){
          bp_est <- data.frame(species=this_species, data=this_dat, model=best_model, breakpt=thresh$est, breakpt_se=thresh$std.error)
      } else {
          bp_est <- bind_rows(bp_est, data.frame(species=this_species, model=best_model, data=this_dat, breakpt=thresh$est, breakpt_se=thresh$std.error))
      }
    }
  }
}

##Plot line plot of breakpoint estimates
bp_est$data <- factor(bp_est$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")

ggplot(bp_est, aes(y=species, x=breakpt, colour=data, shape=model))+
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=data, shape=model), size=3, position=ggstance::position_dodgev(height=0.4))+
  geom_linerange(aes(xmin = breakpt-breakpt_se, xmax = breakpt+breakpt_se, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
  theme_minimal()+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  #xlim(-1,10)+
  scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
  guides(colour=guide_legend(nrow=2,byrow=TRUE),shape=guide_legend(nrow=2,byrow=TRUE))+
  theme(panel.spacing = unit(3, "lines"))+
  scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  xlab("Metabolic Index Breakpoint Estimate")+
  ylab("Species")

ggsave(
  paste("output/plots/breakpt_est_all.png"),
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

ggplot(bp_est, aes(y=species, x=breakpt, colour=data, shape=model))+
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=data, shape=model), size=3, position=ggstance::position_dodgev(height=0.4))+
  geom_linerange(aes(xmin = breakpt-breakpt_se, xmax = breakpt+breakpt_se, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
  theme_minimal()+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
  guides(colour=guide_legend(nrow=2,byrow=TRUE),shape=guide_legend(nrow=2,byrow=TRUE))+
  theme(panel.spacing = unit(3, "lines"))+
  scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  xlab("Metabolic Index Breakpoint Estimate")+
  ylab("Species")+
  coord_cartesian(xlim=c(0, 15))

ggsave(
  paste("output/plots/breakpt_est_truncated.png"),
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


##Calculate pO2 at a reference temperature and body size
#Reference temp (12 deg C)
ref_temp <- 12
#Body size at 1kg
body_size <- 1

for(i in 1:nrow(bp_est)){
temp <- bp_est[i,]
#find taxa from species
taxa.2.use <- taxa$MI_Taxa[taxa$common_name==temp$species]
#Calc invtemp
invtemp.2.use <- 1/(273.15+ref_temp)
#Model
model.2.use <- temp$model
#calculate pO2 at a reference temperature and body size
temp$est_o2 <- calc_po2_crit(taxa.2.use,temp$breakpt,invtemp.2.use,body_size, temp$model)
temp$est_o2_se <- calc_po2_crit(taxa.2.use,temp$breakpt_se,invtemp.2.use,body_size, temp$model)
if(i==1){
  bp_est2 <- temp
} else {
  bp_est2 <- bind_rows(bp_est2, temp)
}
}

#Plot
##Plot line plot of breakpoint estimates
bp_est$data <- factor(bp_est$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")

ggplot(bp_est2, aes(y=species, x=est_o2, colour=data))+
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=data), size=2, position=ggstance::position_dodgev(height=0.4))+
  geom_linerange(aes(xmin = est_o2-est_o2_se, xmax = est_o2+est_o2_se, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
  theme_minimal()+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
 # xlim(0,50)+
  scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
  guides(colour=guide_legend(nrow=1,byrow=TRUE),shape=guide_legend(nrow=2,byrow=TRUE))+
  theme(panel.spacing = unit(3, "lines"))+
  xlab(bquote(pO[2]~"(kPa)")) +
  ylab("Species")

ggsave(
  paste("output/plots/breakpt_est_o2_all.png"),
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
ggplot(bp_est2, aes(y=species, x=est_o2, colour=data))+
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=data), size=2, position=ggstance::position_dodgev(height=0.4))+
  geom_linerange(aes(xmin = est_o2-est_o2_se, xmax = est_o2+est_o2_se, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
  theme_minimal()+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  # xlim(0,50)+
  scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
  guides(colour=guide_legend(nrow=1,byrow=TRUE),shape=guide_legend(nrow=2,byrow=TRUE))+
  theme(panel.spacing = unit(3, "lines"))+
  xlab(bquote(pO[2]~"(kPa)")) +
  ylab("Species")+
  coord_cartesian(xlim=c(0, 20))

ggsave(
  paste("output/plots/breakpt_est_o2_truncated.png"),
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
remove.packages("sdmTMB")
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE,  ref="newbreakpt")
library(sdmTMB)
library(sp)
library(dplyr)
library(tidyr)
library(pals)
library(purrr)
library(ggplot2)
library(ggstance)
library(readxl)
library(ggh4x)
library(viridis)

set.seed(9876)

#Set working directory
setwd("~/Dropbox/GitHub/wsg-choke-species")

#Load functions
source("code/helper_funs.R")

#Load AIC table for model output
#Data type comparison
data_type <- F
region_comp <- T
if(data_type){
  aic <- as.data.frame(read.csv("output/data_type/aic_table_data_type_priors_goodonly.csv"))
  output_folder <- "data_type"
}
if(region_comp){
  aic <- as.data.frame(read_excel("output/region_comp/aic_table_region_comp_priors_goodonly.xlsx"))
  output_folder <- "region_comp"
}

#Species to run
species <- unique(aic$species)

#Taxa lookup
taxa <- read_excel("data/species_table.xlsx")
taxa$MI_Taxa <- tolower(taxa$MI_Taxa)
taxa$common_name <- tolower(taxa$common_name)

#Range of temps
##Create sequence of metabolic index values for marginal effects, so same for all
#Load data
files <- list.files(path = "data/processed_data/fish2", pattern = ".rds", full.names=T)
dat <- map(files,readRDS)
dat <- bind_rows(dat)

#Remove NAs (and remove IPHC by removing cpue)
#Remove NAs (keep IPHC)
dat <- dat  %>%
  drop_na(depth,year, mi1,mi2,mi3, X, Y)

#Sequence of temperatures
kelvin = 273.15
boltz = 0.000086173324
tref <- 12

t.range <- seq(min(dat$temperature_C), max(dat$temperature_C), length.out = 100)
t.range2 <-(1 / boltz)  * ( 1 / (t.range+ 273.15) - 1 / (tref + 273.15))

#Breakpoint estimates of MI
for(i in 1:length(species)) {
  this_species = species[i]
  print(this_species)
  this_aic <- as.data.frame(filter(aic, species==this_species))
  this_aic$data.type <- this_aic[,"data type"]
  dat_names <- unique(this_aic$data.type)
  for(h in 1:length(dat_names)){
    this_dat <- dat_names[h]
    this_data <- as.data.frame(filter(this_aic, data.type==this_dat))
    #Identify if MI is best-fitting model and pull that model
    mi_models <- this_data[c("model9", "model10", "model11")]
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
      mean_mi <- if(best_model=="model9") mean(this_datframe$mi1) else if(best_model=="model10") mean(this_datframe$mi2) else mean(this_datframe$mi3)
      sd_mi <-   if(best_model=="model9") sd(this_datframe$mi1) else if(best_model=="model10") sd(this_datframe$mi2) else sd(this_datframe$mi3)
      #Pull breakpoint and slope
      pars <- as.data.frame(tidy(fit, effects="fixed", conf.int=T))
      slope <- filter(pars, grepl("slope", term))
      thresh <- filter(pars,grepl("breakpt", term))
      thresh$est <- (thresh$estimate*sd_mi)+mean_mi
      if(thresh$std.error!="NaN"){
      thresh$std.error <- (thresh$std.error*sd_mi)+mean_mi 
      } else {
        thresh$std.error <- NA
      }
      if(thresh$est>0 & (thresh$std.error<10|is.na(thresh$std.error))){
      if(i==1 & h==1){
        bp_est <- data.frame(species=this_species, data=this_dat, model=best_model, breakpt=thresh$est, breakpt_se=thresh$std.error)
      } else {
        bp_est <- bind_rows(bp_est, data.frame(species=this_species, model=best_model, data=this_dat, breakpt=thresh$est, breakpt_se=thresh$std.error))
      }
    }
  }
}
}


##Calculate pO2 at a reference temperature and body size
for(i in 1:nrow(bp_est)){
  test <- bp_est[i,]
  #find taxa from species
  taxa.2.use <- taxa$MI_Taxa[taxa$common_name==test$species]
  #body size
  body_size <- 2
  #Model
  model.2.use <- test$model
  test$breakpt_se1 <- test$breakpt-test$breakpt_se
  test$breakpt_se2 <- test$breakpt+test$breakpt_se
  
  #calculate pO2 at a reference temperature and body size
  est_o2 <- calc_po2_crit(t.range2,taxa.2.use,test$breakpt,body_size, test$model)
  est_o2_se1 <- calc_po2_crit(t.range2,taxa.2.use,test$breakpt_se1,body_size, test$model)
  est_o2_se2 <- calc_po2_crit(t.range2,taxa.2.use,test$breakpt_se2,body_size, test$model)
  if(i==1){
    bp_est2 <- data.frame(species=test$species, data=test$data, model=test$model, breakpt=test$breakpt, breakpt_se=test$breakpt_se, invtemp=t.range2, temp=t.range, po2_crit=est_o2, po2_crit_se1=est_o2_se1, po2_crit_se2=est_o2_se2)
  } else {
   dat3 <- data.frame(species=test$species, data=test$data, model=test$model, breakpt=test$breakpt, breakpt_se=test$breakpt_se, invtemp=t.range2, temp=t.range, po2_crit=est_o2,  po2_crit_se1=est_o2_se1, po2_crit_se2=est_o2_se2)
   bp_est2 <- bind_rows(bp_est2, dat3)
  }
}

##Combine all datasets used
for(i in 1:nrow(bp_est)) {
  dat.2.est <- bp_est[i,]
  this_species = bp_est$species[i]
  print(this_species)
  this_aic <- as.data.frame(filter(aic, species==this_species))
  dat_names <- unique(this_aic$data.type)
  this_dat <- bp_est$data[i]
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
      #Combine data
      this_datframe$model <- paste(best_model)
      this_datframe$data <- paste(this_dat)
      this_datframe$species <- paste(this_species)
    if(i==1){
      dat2 <- this_datframe
    } else {
      dat2 <- bind_rows(dat2, this_datframe)
    }
    }
  }

  
bp_est2$data <- factor(bp_est2$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
dat2$data <- factor(dat2$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")

##Plot temp vs o2 and po2 crits for all species and models
theme_set(theme_bw(base_size = 20))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             theme(strip.background = element_blank()))

#Add labels
dat2$label <- paste(dat2$species, dat2$data, sep=" ")
bp_est2$label <- paste(bp_est2$species, bp_est2$data, sep=" ")

ggplot(data = dat2,aes(x = temperature_C)) +  
  facet_wrap("label", ncol=4)+
  geom_point(aes(y=po2,colour=depth), alpha=0.5, size=0.25)+
  geom_line(data=bp_est2,mapping=aes(x=temp, y=po2_crit), colour="black")+
  geom_ribbon(data=bp_est2, mapping=aes(x=temp, ymin=(po2_crit_se1), ymax=(po2_crit_se2)), alpha=0.5, fill="lightgrey")+
  xlab("Temperature (C)") +
  ylab(bquote(pO[2]~"(kPa)")) +
 # scale_x_continuous(expand = c(0,0), limits = c(0, 15) ) +
 # scale_y_continuous(expand = c(0,0), limits = c(0,30) )+
  theme(legend.position="none")+
  coord_cartesian(ylim=c(0, 40))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=12))+
 # geom_text(aes(label = labels, y=po2, x=temperature_C), data = labels, vjust = 1) +
  scale_colour_viridis(option="mako", guide=guide_colourbar(reverse = TRUE), direction=-1)

ggsave(
  paste("output/plots/temp_o2_truncated.png"),
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

#Plot separately
for(i in 1:nrow(labels)){
  test <- labels[i,]
  this_species <- test$species
  this_dat <- test$data
  dat.2.plot <- filter(dat2, common_name==this_species & data==this_dat)
  dat.2.plot.too <- filter(bp_est2, species==this_species & data==this_dat)
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
    paste("output/plots/obs_crit/obs_crit_", this_species, "_", this_dat, ".png"),
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
ggplot(bp_est2, aes(x=temp, y=po2_crit))+
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

ggsave(
  paste("output/plots/po2_crit_all_no_se.png"),
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

ggplot(bp_est2, aes(x=temp, y=po2_crit))+
  facet_wrap("species", scales="free_y")+
  geom_ribbon(data=bp_est2, mapping=aes(x=temp, ymin=po2_crit_se1, ymax=po2_crit_se2, fill=data),alpha=0.5)+
  geom_line(aes(colour=data), size=1)+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
  scale_fill_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
  xlab("Temperature (C)") +
  ylab(bquote(pO[2]~"(kPa)"))
  #coord_cartesian(xlim=c(0, 15), ylim=c(0, 40))

ggsave(
  paste("output/plots/po2_crit_all_with_se.png"),
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

####Plot w/ Tim's



##Plot temp and depth, and where was O2 suitable and not?
##Calculate pO2 breakpoint for each sampling event temp and species (reference body size)
##Do this just for each region, not for the coastwide models
bp_est3 <- filter(bp_est, data!="coastwide")
for(i in 1:nrow(bp_est3)){
  dat.2.est <- bp_est3[i,]
  this_species = dat.2.est$species
  print(this_species)
  this_dat <- dat.2.est$data
  this_model <- dat.2.est$model
  best_model <- this_model
  #Pull the data file
  this_datframe <-try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
  #Calculate mean of this_datframe$mi based on best model
  mean_mi <- if(best_model=="model3") mean(this_datframe$mi1) else if(best_model=="model4") mean(this_datframe$mi2) else mean(this_datframe$mi3)
  sd_mi <-   if(best_model=="model3") sd(this_datframe$mi1) else if(best_model=="model4") sd(this_datframe$mi2) else sd(this_datframe$mi3)
  #Pull breakpoint and slope
  #find taxa from species
  taxa.2.use <- taxa$MI_Taxa[taxa$common_name==this_species]
  #body size
  body_size <- 2
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
      test$est_o2 <- unlist(lapply(invtemp, calc_po2_crit, taxa.2.use,thresh_est,body_size, model))
      test$est_o2_se1 <- unlist(lapply(invtemp, calc_po2_crit, taxa.2.use,thresh_se1,body_size, model))
      test$est_o2_se2 <- unlist(lapply(invtemp, calc_po2_crit, taxa.2.use,thresh_se2,body_size, model))
     # test$percentile <- (test$po2-test$est_o2)/test$est_o2_se
      test$unsuitable <- ifelse(test$po2<test$est_o2, 1, 0)
      test$unsuitable_low <- ifelse(test$po2<test$est_o2_se1, 1, 0)
      test$unsuitable_high <- ifelse(test$po2<test$est_o2_se2, 1, 0)
      #Sum across row to get total suitable
      test$suitable_total <- rowSums(test[,c("unsuitable", "unsuitable_low", "unsuitable_high")])
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
  scale_colour_manual(values=c("lightblue", "orange3"), labels=c("Above pO2 crit", "Below pO2 crit"))+
  guides(color = guide_legend(title="",override.aes = list(size = 7, alpha=1)))

ggsave(
  paste("output/plots/unsuitable_suitable_depth.png"),
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

##Shaded by CPUE
ggplot(data = filter(dats, unsuitable==0),aes(x = temperature_C, y=depth)) +  
  geom_point(aes(colour=log(cpue_kg_km2+1)), size=0.25, alpha=0.4)+
  scale_colour_distiller(name="Above pO2 crit", type="seq",palette="Blues", direction=1)+
  new_scale_colour() +
 geom_point(data=filter(dats, unsuitable==1), mapping=aes(x = temperature_C, y=depth, colour=log(cpue_kg_km2+1)), size=0.25, alpha=0.4)+
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

##save separate for each. species, with faded by CPUE
for(i in 1:nrow(bp_est3)){
  test <- bp_est3[i,]
  this_species <- test$species
  this_dat <- test$data
  dat.2.plot <- filter(dats, common_name==this_species & data==this_dat)
  p <- ggplot(data = filter(dat.2.plot, unsuitable==0),aes(x = temperature_C, y=depth)) +  
    geom_point(aes(colour=log(cpue_kg_km2+1)), size=0.25, alpha=0.4)+
    scale_colour_distiller(name="Above pO2 crit", type="seq",palette="Blues", direction=1)+
    new_scale_colour() +
    geom_point(data=filter(dats, unsuitable==1), mapping=aes(x = temperature_C, y=depth, colour=log(cpue_kg_km2+1)), size=0.25, alpha=0.4)+
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

##save separate for each. species, with faded by CPUE
for(i in 1:nrow(bp_est3)){
  test <- bp_est3[i,]
  this_species <- test$species
  this_dat <- test$data
  dat.2.plot <- filter(dats, common_name==this_species & data==this_dat)
  p <- ggplot(data = dat.2.plot, aes(x = temperature_C, y=depth)) +
  geom_point(dat.2.plot,mapping=aes(colour=log(cpue_kg_km2+1), alpha=as.factor(unsuitable)), size=0.5)+
  xlab("Temperature (C)") +
  ylab("Depth (m)") +
  scale_y_reverse()+
  # scale_x_continuous(expand = c(0,0), limits = c(0, 15) ) +
  # scale_y_continuous(expand = c(0,0), limits = c(0,30) )+
  scale_colour_viridis_c(option="viridis")+
  scale_alpha_manual(values=c(1,0.2))+
  guides(alpha="none", colour="none")+
theme(legend.title=element_blank(), legend.position=c(0.8,0.2))+
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

# setup up mapping ####
map_data <- rnaturalearth::ne_countries(scale = "large",
                                        returnclass = "sf",
                                        continent = "North America")

us_coast_proj <- sf::st_transform(map_data, crs = 32610)

###Map of data available
ggplot(us_coast_proj) + geom_sf() +
    xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
    ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+


##Make model predictions from data and marginal effects


##Plot residuals
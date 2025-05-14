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
both <- T

####Plot data fit in models
##Data names
if(!iphc){
dat_names <- c("cc", "bc", "goa", "ebs", "coastwide")
species <- read_excel("data/species_table.xlsx")
species$common_name <- tolower(species$common_name)
species$scientific_name <- tolower(species$scientific_name)
species_table <- species
species <- unique(species_table$common_name)
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
  species_table <- species
  species <- unique(species$common_name)
  species <- unique(species_table$common_name)
  species_iphc <- c("sablefish", "pacific cod", "pacific halibut", "yelloweye rockfish", "longnose skate", "big skate", "spiny dogfish", "rougheye rockfish")
  
}

#Output folder
output_folder <- "region_comp"

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
      geom_point(filter(dat2plot,catch>0),mapping=aes(x=X*1000, y=Y*1000,colour=survey), size=0.1)+
      xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
      ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
      facet_wrap("year", ncol=5)+
      theme_minimal(base_size=12)+
      xlab("Longitude")+
      ylab("Latitude")+
      ggtitle(paste(this_species, this_dat, sep=" "))+
      theme(axis.text.x=element_blank())
    
    
    ggsave(
      paste0(output_folder, "/", "plots/data_fit_mapping/map_",this_species,"_", this_dat,".png"),
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
##Sequence of MI values--to use later, so same for all
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
      #Calculate max MI in data
      max_mi <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") max(this_datframe$mi1, na.rm=T) else if(best_model=="model4"|best_model=="model10"|best_model=="model14") max(this_datframe$mi2, na.rm=T) else max(this_datframe$mi3, na.rm=T)
      #Pull breakpoint and slope
      pars <- as.data.frame(tidy(fit, effects="fixed", conf.int=T))
      slope <- filter(pars, grepl("slope", term))
      thresh <- filter(pars,grepl("breakpt", term))
      thresh$est <- (thresh$estimate*sd_mi)+mean_mi
      thresh$low <- ((thresh$estimate-thresh$std.error)*sd_mi)+mean_mi 
      thresh$high <- ((thresh$estimate+thresh$std.error)*sd_mi)+mean_mi
      thresh$bp_unscaled <- thresh$estimate
      thresh$bp_se_unscaled <- thresh$std.error
      if(!exists("bp_est")){
          bp_est <- data.frame(species=this_species, data=this_dat, model=best_model, breakpt=thresh$est, breakpt_se1=thresh$low, breakpt_se2=thresh$high, max_mi=max_mi, slope=slope$estimate, slope_se=slope$std.error, bp_unscaled=thresh$bp_unscaled, bp_se_unscaled=thresh$bp_se_unscaled)
      } else {
          bp_est <- bind_rows(bp_est, data.frame(species=this_species, model=best_model, data=this_dat, breakpt=thresh$est, breakpt_se1=thresh$low, breakpt_se2=thresh$high, max_mi=max_mi, slope=slope$estimate, slope_se=slope$std.error, bp_unscaled=thresh$bp_unscaled, bp_se_unscaled=thresh$bp_se_unscaled))
      }
    }
  }
}

##Plot line plot of breakpoint estimates--all
if(!iphc && !both){
bp_est$data <- factor(bp_est$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")

ggplot(bp_est, aes(y=species, x=breakpt, colour=data))+
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=data, shape=model), size=3, position=ggstance::position_dodgev(height=0.4))+
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
  scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  xlab(" Breakpoint Estimate")+
  ylab("Species")
 # geom_vline(xintercept=0, linetype="dashed")+
  #geom_vline(xintercept=1, linetype="dashed")
}

if(both){
  #Add column for IPHC versus not
  bp_est$type <- ifelse(grepl("iphc", bp_est$data), "bottom trawl & IPHC", "bottom trawl only")
  #Create region column
  bp_est$region <- ifelse(grepl("iphc", bp_est$data), gsub(" _iphc", "", bp_est$data), bp_est$data)
  bp_est$region <- factor(bp_est$region, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
  labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
  names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")
  
  ggplot(bp_est, aes(y=species, x=breakpt, colour=data))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=region, shape=type), size=3, position=ggstance::position_dodgev(height=0.4))+
    #Can add back shape
    geom_linerange(aes(xmin = breakpt_se1, xmax = breakpt_se2, colour=region),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
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
    ylab("Species")
    #geom_vline(xintercept=0, linetype="dashed")+
   # geom_vline(xintercept=1, linetype="dashed")
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

##Remove estimates below zero, and with SE greater than 10, and above
if(!iphc && !both){
ggplot(filter(bp_est, (breakpt_se2<max_mi & breakpt_se1>0 &(breakpt_se2-breakpt<5))), aes(y=species, x=breakpt, colour=data))+ #can add back shape
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=data, shape=model), size=3, position=ggstance::position_dodgev(height=0.4))+
    geom_linerange(aes(xmin = breakpt_se1, xmax = breakpt_se2, colour=data),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5, drop=FALSE)+
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
  ggplot(filter(bp_est, (breakpt_se2<max_mi & breakpt_se1>0 &(breakpt_se2-breakpt<5))), aes(y=species, x=breakpt, colour=data))+ #can add back shape
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=data, shape=model), size=3, position=ggstance::position_dodgev(height=0.4), drop=FALSE)+
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

if(both){
  ggplot(filter(bp_est, (breakpt_se2<max_mi & breakpt_se1>0 &((breakpt_se2-breakpt)<5))), aes(y=species, x=breakpt, colour=region))+ #can add back shape
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=region, shape=type), size=3, position=ggstance::position_dodgev(height=0.4))+
  geom_linerange(aes(xmin = breakpt_se1, xmax = breakpt_se2, colour=region),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
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
  scale_shape_manual(values=c(17,20))+
  theme(legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(0.25, "lines"), #Minimize legend space
        panel.spacing = unit(5, "lines"))+ #Make more space between species
  #scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  xlab("Metabolic Index Breakpoint Estimate")+
  ylab("Species")
#geom_vline(xintercept=0, linetype="dashed")+
# geom_vline(xintercept=1, linetype="dashed")
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

#Just IPHC species
if(both){
ggplot(filter(bp_est, (breakpt_se2<max_mi & breakpt_se1>0 &((breakpt_se2-breakpt)<5)&species %in% species_iphc)), aes(y=species, x=breakpt, colour=region))+ #can add back shape
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=region, shape=type), size=3, position=ggstance::position_dodgev(height=0.1))+
  geom_linerange(aes(xmin = breakpt_se1, xmax = breakpt_se2, colour=region),  position=ggstance::position_dodgev(height=0.1), size=1, alpha=0.5)+
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
  scale_shape_manual(values=c(17,20))+
  theme(legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(0.25, "lines"), #Minimize legend space
        panel.spacing = unit(5, "lines"))+ #Make more space between species
  #scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  xlab("Metabolic Index Breakpoint Estimate")+
  ylab("Species")

ggsave(
  paste0("output/", output_folder, "/breakpt_est_truncated_iphc_only.png"),
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

#Bottom trawl only
ggplot(filter(bp_est, (breakpt_se2<max_mi & breakpt_se1>0 &((breakpt_se2-breakpt)<5)&type=="bottom trawl only")), aes(y=species, x=breakpt, colour=region))+ #can add back shape
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=region, shape=model), size=3, position=ggstance::position_dodgev(height=0.1))+
  geom_linerange(aes(xmin = breakpt_se1, xmax = breakpt_se2, colour=region),  position=ggstance::position_dodgev(height=0.1), size=1, alpha=0.5)+
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
  scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  theme(legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(0.25, "lines"), #Minimize legend space
        panel.spacing = unit(5, "lines"))+ #Make more space between species
  #scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  xlab("Metabolic Index Breakpoint Estimate")+
  ylab("Species")

ggsave(
  paste0("output/", output_folder, "/breakpt_est_truncated_bottom_trawl_only.png"),
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

#No shape for Eo
ggplot(filter(bp_est, (breakpt_se2<max_mi & breakpt_se1>0 &((breakpt_se2-breakpt)<5)&type=="bottom trawl only")), aes(y=species, x=breakpt, colour=region))+ #can add back shape
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=region), size=3, position=ggstance::position_dodgev(height=0.1))+
  geom_linerange(aes(xmin = breakpt_se1, xmax = breakpt_se2, colour=region),  position=ggstance::position_dodgev(height=0.1), size=1, alpha=0.5)+
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
  scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  theme(legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(0.25, "lines"), #Minimize legend space
        panel.spacing = unit(5, "lines"))+ #Make more space between species
  #scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  xlab("Metabolic Index Breakpoint Estimate")+
  ylab("Species")

ggsave(
  paste0("output/", output_folder, "/breakpt_est_truncated_bottom_trawl_only_noEoshape.png"),
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

##pO2 crit: Calculate pO2 at a reference temperature and body size
#Reference temp (15 deg C)
ref_temp <- 15
#Body size at 1kg
body_size <- 2

##Filter for breakpoints that are less than max_mi, and greater than 0, and se <5
bp_est <- filter(bp_est, (breakpt_se2<max_mi & breakpt_se1>0 &(breakpt_se2-breakpt<5)))

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
if(!iphc && !both){
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
  bp_est2$region <- factor(bp_est2$region, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
  labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
  names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")
  
  ggplot(bp_est2, aes(y=species, x=est_o2, colour=data))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=region, shape=type), size=2, position=ggstance::position_dodgev(height=0.4))+
    #Can add shape back
    geom_linerange(aes(xmin = est_o2_low, xmax = est_o2_high, colour=region),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    # xlim(0,50)+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
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
    xlab(bquote(pO[2]~"(kPa) at 15 C")) +
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

#
if(both){
  #IPHC only
  ggplot(filter(bp_est2, species %in% species_iphc), aes(y=species, x=est_o2, colour=region))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=region, shape=type), size=2, position=ggstance::position_dodgev(height=0.1))+
    #Can add shape back
    geom_linerange(aes(xmin = est_o2_low, xmax = est_o2_high, colour=region),  position=ggstance::position_dodgev(height=0.1), size=1, alpha=0.5)+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    # xlim(0,50)+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
    scale_shape_manual(values=c(17,20))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(1, "lines"))+ #Make more space between species
    # scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    xlab(bquote("Critical"~pO[2]~"(kPa) at 15 C")) +
    ylab("Species")

ggsave(
  paste0("output/", output_folder, "/breakpt_est_o2_iphc_only.png"),
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

#Bottom trawl only
ggplot(filter(bp_est2, type=="bottom trawl only"), aes(y=species, x=est_o2, colour=region))+
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=region), size=2, position=ggstance::position_dodgev(height=0.1))+
  #Can add shape back
  geom_linerange(aes(xmin = est_o2_low, xmax = est_o2_high, colour=region),  position=ggstance::position_dodgev(height=0.1), size=1, alpha=0.5)+
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
        panel.spacing = unit(1, "lines"))+ #Make more space between species
  # scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  xlab(bquote("Critical"~pO[2]~"(kPa) at 15 C")) +
  ylab("Species")

ggsave(
  paste0("output/", output_folder, "/breakpt_est_o2_bottom_trawl_only.png"),
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

#Add laboratory estimates (just single points)
##Add laboratory lines
#Only show reasonable ones
if(!iphc && !both){
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
    xlab(bquote("Critical"~pO[2]~"(kPa) at 15 C")) +
  ylab("Species")+
  scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
 # coord_cartesian(xlim=c(0, 75))+
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
    xlab(bquote("Critical"~pO[2]~"(kPa) at 15 C")) +
    ylab("Species")+
    scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    coord_cartesian(xlim=c(0, 75))+
    scale_x_continuous(expand=c(0,0))
}
if(both){
ggplot(filter(bp_est2, type=="bottom trawl only"),aes(y=species, x=est_o2, colour=region))+
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=region), size=2, position=ggstance::position_dodgev(height=0.4))+
  geom_point(mapping=aes(x=est_o2_lab), colour="black", shape="triangle")+
  geom_linerange(aes(xmin = est_o2_low, xmax = est_o2_high, colour=region),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
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
    xlab(bquote("Critical"~pO[2]~"(kPa) at 15 C")) +
  ylab("Species")+
  scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  # coord_cartesian(xlim=c(0, 75))+
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
if(!iphc && !both){
ggplot(bp_est5, aes(y=species, x=est_o2, colour=data))+
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=data), size=2, position=ggstance::position_dodgev(height=0.4))+
 # geom_point(mapping=aes(x=est_o2_lab), colour="black", shape="triangle")+
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
 scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6", "black"), drop=FALSE)+
  # guides(colour=guide_legend(nrow=1,byrow=TRUE),shape=guide_legend(nrow=1,byrow=TRUE))+
  guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
  theme(legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(0.25, "lines"), #Minimize legend space
        panel.spacing = unit(5, "lines"))+ #Make more space between species
    xlab(bquote("Critical"~pO[2]~"(kPa) at 15 C")) +
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

if(both){
  ggplot(bp_est5, aes(y=species, x=est_o2, colour=region))+
    #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
    geom_point(aes(colour=region), size=2, position=ggstance::position_dodgev(height=0.1))+
    # geom_point(mapping=aes(x=est_o2_lab), colour="black", shape="triangle")+
    geom_point(mapping=aes(x=est_o2_lab), colour="black", shape="triangle")+
    geom_linerange(filter(bp_est5, data!="laboratory phylogenetic imputation"),mapping=aes(xmin = est_o2_low, xmax = est_o2_high, colour=region),  position=ggstance::position_dodgev(height=0.1), size=1, alpha=0.5)+
    geom_linerange(filter(bp_est5, data=="laboratory phylogenetic imputation"),mapping=aes(xmin = est_o2_low, xmax = est_o2_high), colour="black", alpha=0.5, linetype="dashed")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=12))+
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    theme(text=element_text(size=15))+
    # xlim(0,50)+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6", "black"),labels=c("cc", "bc", "goa", "ebs", "coastwide", "lab-derived"), drop=FALSE)+
    # guides(colour=guide_legend(nrow=1,byrow=TRUE),shape=guide_legend(nrow=1,byrow=TRUE))+
    guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
    theme(legend.box = "vertical",
          legend.spacing.y = unit(0, "pt"),
          legend.key.height = unit(0.25, "lines"), #Minimize legend space
          panel.spacing = unit(5, "lines"))+ #Make more space between species
    xlab(bquote("Critical"~pO[2]~"(kPa) at 15 C")) +
    ylab("Species")+
    scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
    #coord_cartesian(xlim=c(0, 75))+
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

#Save as CSV
#Save as csv
#Round
bp_est_save <- bp_est2
bp_est_save$breakpt<- round(as.numeric(bp_est_save$breakpt), digits=3)
bp_est_save$breakpt_se1<- round(as.numeric(bp_est_save$breakpt_se1), digits=3)
bp_est_save$breakpt_se2<- round(as.numeric(bp_est_save$breakpt_se2), digits=3)
bp_est_save$est_o2<- round(as.numeric(bp_est_save$est_o2), digits=3)
bp_est_save$est_o2_low<- round(as.numeric(bp_est_save$est_o2_low), digits=3)
bp_est_save$est_o2_high<- round(as.numeric(bp_est_save$est_o2_high), digits=3)
bp_est_save$est_o2_lab<- round(as.numeric(bp_est_save$est_o2_lab), digits=3)
write.csv(bp_est_save, file=paste0("output/", output_folder, "/breakpoint_est.csv"))

##Loop to pull each set of comparisons that converged, pull the best-fitting model if includes MI, predict to a marginal effect prediction grid, and make a dataframe with them
#Create dataframe
cond_effects_preds = as.data.frame(matrix(NA, 1, 18))

for(i in 1:nrow(bp_est2)) {
  this_species = bp_est2$species[i]
  print(this_species)
  this_aic <- as.data.frame(filter(aic, species==this_species))
  this_aic$data.type <- this_aic$"data type"
  this_dat <- bp_est2$data[i]
  print(this_dat)
  this_data <- as.data.frame(filter(this_aic, data.type==this_dat&species==this_species))
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
        ##Add calculation of pO2 at reference temperature
       # p1$po2_15 <- calc_po2_crit(t.range2,taxa.2.use,test$breakpt,body_size, test$model, fancy=F)
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

##Clean up dataframe
#Remove first row
cond_effects_preds <- cond_effects_preds[-1,]
#Remove first 18 columns
cond_effects_preds <- cond_effects_preds[,19:ncol(cond_effects_preds)]


if(!iphc && !both){
  cond_effects_preds$data <- factor(cond_effects_preds$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
  labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
  names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")
  
  ##Plot conditional effects
  ggplot(cond_effects_preds, aes(mi_best, y=est_sc))+
    facet_wrap("species", ncol=4)+
    geom_line(aes(colour=data))+
    #geom_ribbon(aes(ymin = est_se_sc1, ymax = est_se_sc2, fill=data), alpha=0.4)+
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
  #Add column for IPHC versus not
  cond_effects_preds$type <- ifelse(grepl("iphc", cond_effects_preds$data), "bottom trawl & IPHC", "bottom trawl only")
  #Create region column
  cond_effects_preds$region <- ifelse(grepl("iphc", cond_effects_preds$data), gsub(" _iphc", "", cond_effects_preds$data), cond_effects_preds$data)
  cond_effects_preds$region <- factor(cond_effects_preds$region, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
  labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
  names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")
  
  ##Plot conditional effects
  #Bottom trawl only
  ggplot(filter(cond_effects_preds, type=="bottom trawl only"),aes(mi_best, y=est_sc))+
    facet_wrap("species", ncol=4)+
    geom_line(aes(colour=region))+
    # geom_ribbon(aes(ymin = est_se_sc1, ymax = est_se_sc2, fill=data), alpha=0.4)+
    scale_x_continuous(limits=c(0,15, by=5))+
    labs(x = bquote("Temp-Corrected"~pO[2]~"(kPa)"), y = bquote('Conditional Effect Population Density'~(kg~km^-2)))+
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
  
}

ggsave(
  paste0("output/", output_folder, "/cond_effects_region_mi_scaled_no_se_bottom_trawl_only.png"),
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

if(!iphc && !both){
bp_est3$data <- factor(bp_est3$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
dat2$data <- factor(dat2$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")
}
if(both){
  #Add column for IPHC versus not
  bp_est3$type <- ifelse(grepl("iphc", bp_est3$data), "bottom trawl & IPHC", "bottom trawl only")
  #Create region column
  bp_est3$region <- ifelse(grepl("iphc", bp_est3$data), gsub(" _iphc", "", bp_est3$data), bp_est3$data)
  bp_est3$region <- factor(bp_est3$region, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
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
if(!iphc &&!both){
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

if(both){
  ggplot(filter(bp_est3, type=="bottom trawl only"), aes(x=temp, y=po2_crit))+
    facet_wrap("species", scales="free_y")+
    #geom_ribbon(data=bp_est2, mapping=aes(x=temp, ymin=(po2_crit_se1), ymax=(po2_crit_se2), fill=data), alpha=0.3)+
    geom_line(aes(colour=region), size=1)+
    theme(legend.position="top")+
    theme(legend.title=element_blank(), strip.background = element_blank())+
    theme(text=element_text(size=15))+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"), drop=FALSE)+
    xlab("Temperature (C)") +
    ylab(bquote(pO[2]~"(kPa)"))
}
ggsave(
  paste0("output/", output_folder, "/po2_crit_all_no_se.png"),
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

if(!iphc && !both){
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

if(both){
  ggplot(bp_est3, aes(x=temp, y=po2_crit))+
    facet_wrap("species", scales="free_y")+
    geom_ribbon(data=bp_est3,mapping=aes(x=temp, ymin=po2_crit_se1, ymax=po2_crit_se2, fill=region),alpha=0.5)+
    geom_line(aes(colour=region), size=1)+
    theme(legend.position="top")+
    theme(legend.title=element_blank(), strip.background = element_blank())+
    theme(text=element_text(size=15))+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"), drop=FALSE)+
    scale_fill_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"), drop=FALSE)+
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

if(!iphc && !both){
ggplot(bp_est3, aes(x=temp, y=po2_crit))+
  facet_wrap("species", scales="free_y")+
    geom_polygon(data = dats50a, mapping=aes(x = temp, y = ys), fill = "grey1", color = NA, alpha = 0.2) + 
    geom_polygon(data = dats90a, mapping=aes(x = temp, y = ys), fill = "lightgrey", color = NA, alpha = 0.2) +
  geom_ribbon(data=bp_est3, mapping=aes(x=temp, ymin=po2_crit_se1, ymax=po2_crit_se2, fill=data),alpha=0.2)+
  geom_line(aes(colour=data), size=1)+
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
if(both){
  ggplot(bp_est3, aes(x=temp, y=po2_crit))+
    facet_wrap("species", scales="free_y")+
    geom_ribbon(bp_est3, mapping=aes(x=temp, ymin=po2_crit_se1, ymax=po2_crit_se2, fill=region),alpha=0.2)+
    geom_line(aes(colour=region), size=1)+
    geom_polygon(data = dats50a, mapping=aes(x = temp, y = ys), fill = "grey1", color = NA, alpha = 0.3) + 
    geom_polygon(data = dats90a, mapping=aes(x = temp, y = ys), fill = "grey", color = NA, alpha = 0.3) +
    theme(legend.position="top")+
    theme(legend.title=element_blank(), strip.background = element_blank())+
    theme(text=element_text(size=15))+
    scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"), drop=FALSE)+
    scale_fill_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"), drop=FALSE)+
    xlab("Temperature (C)") +
    ylab(bquote(pO[2]~"(kPa)"))
  #coord_cartesian(xlim=c(0, 15), ylim=c(0, 40))+
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

saveRDS(bp_est4, file=paste0("output/", output_folder, "/bp_est4.rds"))
bp_est4 <- readRDS(file=paste0("output/", output_folder, "/bp_est4.rds"))

##In data, how much was biomass reduced by oxygen at each point?
bp_est4 <- bp.2.use
for(i in 1:nrow(bp_est4)){
  dat.2.est <- bp_est4[i,]
  this_species = dat.2.est$species
  print(this_species)
  this_dat <- dat.2.est$data
  this_model <- dat.2.est$model
  best_model <- dat.2.est$model
  slope.2.use <- dat.2.est$slope
  breakpt.2.use <- dat.2.est$bp_unscaled
  breakpt.2.use_se1 <- breakpt.2.use-dat.2.est$bp_se_unscaled
  breakpt.2.use_se2 <- breakpt.2.use+dat.2.est$bp_se_unscaled
  taxa_table <- filter(species_table, common_name==this_species)
  
  #Pull the data file 
  this_datframe <-try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
  #Which MI to use?
  mi.2.use <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") this_datframe$mi1_s else if(best_model=="model4"|best_model=="model10"|best_model=="model14") this_datframe$mi2_s else this_datframe$mi3_s
  #Calculate raw effect
  this_datframe$est_effect_raw <- sapply(mi.2.use,breakpoint_calc, slope.2.use,breakpt.2.use)
  this_datframe$est_effect_raw_se1<- sapply(mi.2.use,breakpoint_calc, slope.2.use,breakpt.2.use_se1)
  this_datframe$est_effect_raw_se2 <- sapply(mi.2.use,breakpoint_calc, slope.2.use,breakpt.2.use_se2)
  
  #Exponentiate effect
  this_datframe$est_effect_raw <- exp(this_datframe$est_effect_raw)
  this_datframe$est_effect_raw_se1 <- exp(this_datframe$est_effect_raw_se1)
  this_datframe$est_effect_raw_se2 <- exp(this_datframe$est_effect_raw_se2)
  
  #Calculate max effect
  max_effect <- exp(breakpoint_calc(breakpt.2.use, slope.2.use, breakpt.2.use))
  max_effect_se1 <- exp(breakpoint_calc(breakpt.2.use_se1, slope.2.use, breakpt.2.use))
  max_effect_se2 <- exp(breakpoint_calc(breakpt.2.use_se2, slope.2.use, breakpt.2.use))
  this_datframe$max_effect <- max_effect
  this_datframe$max_effect_se1 <- max_effect_se1
  this_datframe$max_effect_se2 <- max_effect_se2
  
  #Raw biomass reduction
  this_datframe$biomass_reduction <- max_effect-this_datframe$est_effect_raw
  this_datframe$biomass_reduction_se1 <- max_effect_se1-this_datframe$est_effect_raw_se1
  this_datframe$biomass_reduction_se2 <- max_effect_se2-this_datframe$est_effect_raw_se2
  
  #Proportional biomass reduction
  this_datframe$est_effect_prop <- (max_effect-this_datframe$est_effect_raw)/max_effect
  this_datframe$est_effect_prop_se1 <- (max_effect_se1-this_datframe$est_effect_raw_se1)/max_effect_se1
  this_datframe$est_effect_prop_se2 <- (max_effect_se2-this_datframe$est_effect_raw_se2)/max_effect_se2
  
  #Other columns
  this_datframe$model <- this_model
  this_datframe$data <- this_dat 
  this_datframe$depth_limit <- taxa_table$depth
  this_datframe$min_lat <- taxa_table$min_lat
  this_datframe$max_lat <- taxa_table$max_lat
  
  if(i==1){
    effects <- this_datframe
  }
  if(i>0){
    effects <- bind_rows(effects, this_datframe)
  }
}

##Overall (not spatial)
ggplot(filter(effects, biomass_reduction>0&!grepl("iphc", data)), aes(x=biomass_reduction))+ 
  geom_density(aes(group=region, colour=region))+
  facet_wrap("common_name")+
  coord_cartesian(xlim=c(0,0.5))
ggplot(filter(effects, biomass_reduction>0&grepl("iphc", data)), aes(x=biomass_reduction))+ 
  geom_density(aes(group=region, colour=region))+
  facet_wrap("common_name")+
  coord_cartesian(xlim=c(0,0.5))

#Plot spatially, hexagon
map_data <- rnaturalearth::ne_countries(scale = "large",
                                        returnclass = "sf",
                                        continent = "North America")

us_coast_proj <- sf::st_transform(map_data, crs = 32610)

for(i in 1:nrow(bp_est4)) {
  this_species <- bp_est4$species[i]
  print(this_species)
  taxa_table <- filter(species_table, common_name==this_species)
  this_dat <- bp_est4$data[i]
  dat2plot <- filter(effects, common_name==this_species & data==this_dat)
  dat2plot <- filter(dat2plot, depth<depth_limit)
  dat2plot <- filter(dat2plot, (latitude>(taxa_table$min_lat-0.5))&(latitude<(taxa_table$max_lat+0.5)))
  
  #Restrict dataset to species ranges
  #if(this_species=="canary rockfish"|this_species=="silvergray rockfish"|this_species=="yellowtail rockfish"){
#    dat2plot <- filter(dat2plot, region!="ebs")
 # }
 # if(this_species=="pacific cod"|this_species=="pacific halibut"|this_species=="walleye pollock"){
  #  min_lat <- 34
   # dat2plot <- filter(dat2plot, latitude>min_lat)
 # }

   ggplot(us_coast_proj) + geom_sf() +
        geom_hex(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000), size=0.1, bins=30)+
        xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
        ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
        #facet_wrap("year", ncol=5)+
        theme_minimal(base_size=12)+
        xlab("Longitude")+
        ylab("Latitude")+
        ggtitle(paste(this_species, this_dat, sep=" "))+
        theme(axis.text.x=element_blank())+
        scale_fill_viridis()
      
      ggsave(
        paste0("output/", output_folder, "/", "plots/counts/map_",this_species,"_", this_dat,".png"),
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

##Plotting each point
for(i in 1:nrow(bp_est4)) {
  this_species <- bp_est4$species[i]
  print(this_species)
  taxa_table <- filter(species_table, common_name==this_species)
  this_dat <- bp_est4$data[i]
  dat2plot <- filter(effects, common_name==this_species & data==this_dat)
  dat2plot <- filter(dat2plot, depth<depth_limit)
  dat2plot <- filter(dat2plot, (latitude>(taxa_table$min_lat-0.5))&(latitude<(taxa_table$max_lat+0.5)))

  #Restrict dataset to species ranges
  #if(this_species=="canary rockfish"|this_species=="silvergray rockfish"|this_species=="yellowtail rockfish"){
  #    dat2plot <- filter(dat2plot, region!="ebs")
  # }
  # if(this_species=="pacific cod"|this_species=="pacific halibut"|this_species=="walleye pollock"){
  #  min_lat <- 34
  # dat2plot <- filter(dat2plot, latitude>min_lat)
  # }
  
  ggplot(us_coast_proj) + geom_sf() +
    #geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="#0D0887FF", size=0.5, alpha=0.1)+
    geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.5, alpha=0.1)+
    geom_point(filter(dat2plot,est_effect_prop>0),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.9)+
    #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
    ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
    #facet_wrap("year", ncol=5)+
    theme_minimal(base_size=18)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(axis.text.x=element_blank(), legend.position=c(0.8,0.15))+
    ggtitle(paste(this_species, this_dat, sep=" "))+
    #scale_colour_viridis(name="proportional biomass reduction \n from oxygen")
    scale_colour_viridis(name="proportional biomass \nreduction from oxygen", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish, option="plasma")
  
  ggsave(
    paste0("output/", output_folder, "/", "plots/points/map_",this_species,"_", this_dat,".png"),
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

###With entire oxygen dataset
bp_est4 <- bp.2.use

#Get all insitu data
insitu <- readRDS("data/processed_data/o2/insitu_combined.rds")
#Just columns of interest
dat <- insitu[,c("survey", "depth", "year", "date", "latitude", "longitude", "X", "Y", "temperature_C", "do_mlpL", "salinity_psu", "sigma0_kgm3", "O2_umolkg", "event_id", "doy")]
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

#Log depth
dat$depth_ln <- log(dat$depth)

#Remove ai
dat <- filter(dat, region!="ai")

###Loop through for each species for all data (entire dataset--not limited to range)
for(i in 1:nrow(bp_est4)){
  dat.2.est <- dat
  this_data <- bp_est4[i,]
  this_species = this_data$species
  print(this_species)
  this_dat <- this_data$data
  this_model <- this_data$model
  best_model <- this_data$model
  slope.2.use <- this_data$slope
  breakpt.2.use <- this_data$bp_unscaled
  breakpt.2.use_se1 <- breakpt.2.use-this_data$bp_se_unscaled
  breakpt.2.use_se2 <- breakpt.2.use+this_data$bp_se_unscaled
  #Pull the data file 
  this_datframe <-try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_dat, "_dat.rds")))

  #Calculate MI
  ###Calculate Metabolic index 
  ##Species parameters from Tim's paper
  #Read table
  taxa_table <- filter(species_table, common_name==this_species)
  taxa <- taxa_table$MI_Taxa
  mi_pars <- read.csv("data/taxa_table.csv")
  mi_pars$Group <- tolower(mi_pars$Group)
  mi_pars <- filter(mi_pars,Group==taxa)
  
  V <- mi_pars$logV
  n <- mi_pars$n
  Eo <- c(mi_pars$Eolow, mi_pars$Eo, mi_pars$Eohigh)
  Ao <- 1/exp(V)
  Eo <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") Eo[1] else if(best_model=="model4"|best_model=="model10"|best_model=="model14") Eo[2] else Eo[3]
  #Abbreviated equation
  #dat$mi1 <- dat$po2*exp(Eo1* dat$invtemp)
  #dat$mi2 <- dat$po2*exp(Eo2* dat$invtemp)
  #dat$mi2 <- dat$po2*exp(Eo2* dat$invtemp)
  fancy <- FALSE
  #Calculate MI for grid
  #Calc invtemp
  kelvin = 273.15
  boltz = 0.000086173324
  tref <- 15
  dat.2.est$invtemp <- (1 / boltz)  * ( 1 / (dat.2.est$temperature_C + 273.15) - 1 / (dat.2.est$temperature_C + 273.15))
  dat.2.est$mi <- calc_mi(Eo=Eo, po2=dat.2.est$po2, inv.temp=dat.2.est$invtemp, fancy=fancy)
  
  #MI from data for scaling
  mi.2.use <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") this_datframe$mi1 else if(best_model=="model4"|best_model=="model10"|best_model=="model14") this_datframe$mi2 else this_datframe$mi3
  mi_mean <- mean(mi.2.use, na.rm=TRUE)
  mi_sd <- sd(mi.2.use, na.rm=TRUE)
  #Scale MI calculated in grid to the fitting dataset
  dat.2.est$mi_s <- (dat.2.est$mi-mi_mean)/mi_sd
  mi_s <- dat.2.est$mi_s
  #Calculate MI effect
  #Calculate raw effect
  dat.2.est$est_effect_raw <- sapply(mi_s,breakpoint_calc, slope.2.use,breakpt.2.use)
  dat.2.est$est_effect_raw_se1<- sapply(mi_s,breakpoint_calc, slope.2.use,breakpt.2.use_se1)
  dat.2.est$est_effect_raw_se2 <- sapply(mi_s,breakpoint_calc, slope.2.use,breakpt.2.use_se2)
  
  #Exponentiate effect
  dat.2.est$est_effect_raw <- exp(dat.2.est$est_effect_raw)
  dat.2.est$est_effect_raw_se1 <- exp(dat.2.est$est_effect_raw_se1)
  dat.2.est$est_effect_raw_se2 <- exp(dat.2.est$est_effect_raw_se2)
  
  #Calculate max effect
  max_effect <- exp(breakpoint_calc(breakpt.2.use, slope.2.use, breakpt.2.use))
  max_effect_se1 <- exp(breakpoint_calc(breakpt.2.use_se1, slope.2.use, breakpt.2.use_se1))
  max_effect_se2 <- exp(breakpoint_calc(breakpt.2.use_se2, slope.2.use, breakpt.2.use_se2))
  dat.2.est$max_effect <- max_effect
  dat.2.est$max_effect_se1 <- max_effect_se1
  dat.2.est$max_effect_se2 <- max_effect_se2
  
  #Raw biomass reduction
  dat.2.est$biomass_reduction <- max_effect-dat.2.est$est_effect_raw
  dat.2.est$biomass_reduction_se1 <- max_effect_se1-dat.2.est$est_effect_raw_se1
  dat.2.est$biomass_reduction_se2 <- max_effect_se2-dat.2.est$est_effect_raw_se2
  
  #Proportional biomass reduction
  dat.2.est$est_effect_prop <- (max_effect-dat.2.est$est_effect_raw)/max_effect
  dat.2.est$est_effect_prop_se1 <- (max_effect_se1-dat.2.est$est_effect_raw_se1)/max_effect_se1
  dat.2.est$est_effect_prop_se2 <- (max_effect_se2-dat.2.est$est_effect_raw_se2)/max_effect_se2
  
  #Other columns
  dat.2.est$model <- this_model
  dat.2.est$data <- this_dat 
  dat.2.est$common_name <- this_species
  dat.2.est$depth_limit <- taxa_table$depth
  dat.2.est$min_lat <- taxa_table$min_lat
  dat.2.est$max_lat <- taxa_table$max_lat

  if(i==1){
    effects_full <- dat.2.est
  }
  if(i>0){
    effects_full <- bind_rows(effects_full, dat.2.est)
  }
}
effects_full <- filter(effects_full, region!="ai")

##Plot
for(i in 1:nrow(bp_est4)) {
  this_species <- bp_est4$species[i]
  print(this_species)
  this_dat <- bp_est4$data[i]
  taxa_table <- filter(species_table, common_name==this_species)
  dat2plot <- filter(effects_full, common_name==this_species & data==this_dat)
  dat2plot <- filter(dat2plot, depth<depth_limit)
  dat2plot <- filter(dat2plot, (latitude>(taxa_table$min_lat-0.5))&(latitude<(taxa_table$max_lat+0.5)))
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.5, alpha=0.1)+
    geom_point(filter(dat2plot,est_effect_prop>0),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.9)+
    #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
   #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
    ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
    #facet_wrap("year", ncol=5)+
    theme_minimal(base_size=12)+
    xlab("Longitude")+
    ylab("Latitude")+
    ggtitle(paste(this_species, this_dat, sep=" "))+
    theme(axis.text.x=element_blank())+
    scale_colour_viridis(name="proportional biomass reduction \n from oxygen", option="magma")
    #scale_colour_viridis(name="proportional biomass reduction \n from oxygen", limits=c(0,1), breaks=c(0,0.25, 0.5, 0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish)
  
  ggsave(
    paste0("output/", output_folder, "/", "plots/points_full_freescales/map_",this_species,"_", this_dat,".png"),
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
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.5, alpha=0.1)+
    geom_point(filter(dat2plot,est_effect_prop>0),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.9)+
    #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
    ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
    #facet_wrap("year", ncol=5)+
    theme_minimal(base_size=12)+
    xlab("Longitude")+
    ylab("Latitude")+
    ggtitle(paste(this_species, this_dat, sep=" "))+
    theme(axis.text.x=element_blank())+
    #scale_colour_viridis(name="proportional biomass reduction \n from oxygen")
    scale_colour_viridis(name="proportional biomass \nreduction from oxygen", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish, option="plasma")
  
  ggsave(
    paste0("output/", output_folder, "/", "plots/points_full_fixedscales/map_",this_species,"_", this_dat,".png"),
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

#Plot just coastwide, on one plot
bp_est8 <- filter(bp_est4, grepl("coastwide", data)) 
#If species in species_iphc, keep if data contains iphc
bp_est8 <- filter(bp_est8, !(species %in% species_iphc & !grepl("iphc", data)))
#Filter effects_full by species and data in bp_est8
effects_full_sub <- filter(effects_full, common_name %in% bp_est8$species & data %in% bp_est8$data)
effects_full_sub <- filter(effects_full_sub, region!="ai")

dat2plot <- filter(effects_full_sub, (depth<depth_limit+200))
dat2plot <- filter(dat2plot, (latitude>(min_lat-0.5))&(latitude<(max_lat+0.5)))

ggplot(us_coast_proj) + geom_sf() +
  #geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="#0D0887FF", size=0.5, alpha=0.1)+
  geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.5, alpha=0.1)+
  geom_point(filter(dat2plot,est_effect_prop>0),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.9)+
   #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  facet_wrap("common_name")+
  xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
  ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
  #facet_wrap("year", ncol=5)+
  theme_minimal(base_size=18)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(axis.text.x=element_blank(), legend.position=c(0.8,0.15))+
  #scale_colour_viridis(name="proportional biomass reduction \n from oxygen")
  scale_colour_viridis(name="proportional biomass \nreduction from oxygen", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish, option="plasma")

ggsave(
  paste0("output/", output_folder, "/", "plots/map_coastwide_species_combined.png"),
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

##Plot annually
for(i in 1:nrow(bp_est4)) {
  this_species <- bp_est4$species[i]
  print(this_species)
  this_dat <- bp_est4$data[i]
  taxa_table <- filter(species_table, common_name==this_species)
  dat2plot <- filter(effects_full, common_name==this_species & data==this_dat)
  dat2plot <- filter(dat2plot, depth<depth_limit)
  dat2plot <- filter(dat2plot, (latitude>(taxa_table$min_lat-0.5))&(latitude<(taxa_table$max_lat+0.5)))
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.5, alpha=0.1)+
    geom_point(filter(dat2plot,est_effect_prop>0),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.9)+
    #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    facet_wrap("year")+
    xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
    ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
    #facet_wrap("year", ncol=5)+
    theme_minimal(base_size=12)+
    xlab("Longitude")+
    ylab("Latitude")+
    ggtitle(paste(this_species, this_dat, sep=" "))+
    theme(axis.text.x=element_blank(), legend.position="top")+
    scale_colour_viridis(name="proportional biomass \nreduction from oxygen", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish, option="plasma")
  
  #  limits=c(0,1), breaks=c(0,0.25, 0.5, 0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish)
  
  ggsave(
    paste0("output/", output_folder, "/", "plots/points_full_annual/map_",this_species,"_", this_dat,".png"),
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

##proportion of observations in each region each year where at least 10% decline in biomass from O2
#Filter to depth limit
test <- filter(effects_full, depth<(depth_limit+200))
test <- filter(test, (latitude>(min_lat-0.5))&(latitude<(max_lat+0.5)))
#Filter to just coastwide species, and IPHC for halibut and cod
test <- filter(test, common_name %in% bp_est8$species & data %in% bp_est8$data)

effects_full_annual_summary <- test %>%
  group_by(year, region, common_name) %>%
  summarise(prop_below = sum(est_effect_prop>0)/n()) %>%
  ungroup()
effects_full_annual_summary2 <- test %>%
  group_by(year, region, common_name) %>%
  summarise(prop_below = sum(est_effect_prop_se1>0)/n()) %>%
  ungroup()
effects_full_annual_summary3 <- test %>%
  group_by(year, region, common_name) %>%
  summarise(prop_below = sum(est_effect_prop_se2>0)/n()) %>%
  ungroup()
effects_full_annual_summary$se1 <- effects_full_annual_summary2$prop_below
effects_full_annual_summary$se2 <- effects_full_annual_summary3$prop_below

#plot
ggplot(effects_full_annual_summary, aes(x=year, y=prop_below))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  facet_wrap("common_name")+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"))+
  scale_fill_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"))+
  xlab("Year") +
  ylab("Proportion of observations below threshold")

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_obs_10_perc_thresh_annually.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 8.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##Evaluate compared to depth
ggplot(filter(effects_full, grepl("coastwide", data)&!grepl("iphc", data)&est_effect_prop>0&depth<depth_limit),aes(x=depth, y=est_effect_prop))+
  geom_point(aes(colour=region), alpha=0.1)+
  facet_wrap("common_name", scales="free")+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"),drop=FALSE)+
  scale_fill_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"), drop=FALSE)+
  xlab("Depth") +
  ylab("Estimated conditional effect of pO2'")

ggsave(
  paste0("output/", output_folder, "/", "plots/depth_vs_effect_coastwide.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 8.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

#Proportion below limit by depth
test$depth_bin <- cut(test$depth, labels=FALSE, breaks=seq(0, 1500, by=50), include.lowest = TRUE)
test$depth_bin <- test$depth_bin*50
effects_full_annual_summary <- test %>%
  group_by(depth_bin, region, common_name) %>%
  summarise(prop_below = sum(est_effect_prop>0)/n()) %>%
  ungroup()
effects_full_annual_summary2 <- test %>%
  group_by(depth_bin, region, common_name) %>%
  summarise(prop_below = sum(est_effect_prop_se1>0)/n()) %>%
  ungroup()
effects_full_annual_summary3 <- test %>%
  group_by(depth_bin, region, common_name) %>%
  summarise(prop_below = sum(est_effect_prop_se2>0)/n()) %>%
  ungroup()
effects_full_annual_summary$se1 <- effects_full_annual_summary2$prop_below
effects_full_annual_summary$se2 <- effects_full_annual_summary3$prop_below

ggplot(effects_full_annual_summary, aes(x=depth_bin, y=prop_below))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  facet_wrap("common_name", scales="free_x")+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"))+
  scale_fill_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"))+
  xlab("Depth (m)") +
  ylab("Proportion of observations below threshold")+
  geom_vline(test, mapping=aes(xintercept=depth_limit), linetype="dashed", colour="black")

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_below_by_depth.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 8.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

ggplot(filter(effects_full, grepl("coastwide", data)&grepl("iphc", data)&est_effect_prop>0&depth<depth_limit),aes(x=depth))+
  geom_density(aes(colour=region), alpha=0.1)+
  facet_wrap("common_name", scales="free")+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"),drop=FALSE)+
  scale_fill_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"), drop=FALSE)+
  xlab("Depth") +
  ylab("Density of Observations Below pO2' threshold")

ggsave(
  paste0("output/", output_folder, "/", "plots/below_limit_by_depth_iphc.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 8.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

ggplot(filter(effects_full, grepl("coastwide", data)&!grepl("iphc", data)&est_effect_prop>0&depth<depth_limit),aes(x=depth))+
  geom_density(aes(colour=region), alpha=0.1)+
  facet_wrap("common_name", scales="free")+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"),drop=FALSE)+
  scale_fill_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"), drop=FALSE)+
  xlab("Depth") +
  ylab("Density of Observations Below pO2' threshold")

ggsave(
  paste0("output/", output_folder, "/", "plots/below_limit_by_depth_bottomtrawl.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 8.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##Limit to just regions where there is data fit
##Plot
for(i in 1:nrow(bp_est4)) {
  this_species <- bp_est4$species[i]
  print(this_species)
  this_dat <- bp_est4$data[i]
  taxa_table <- filter(species_table, common_name==this_species)
  dat2plot <- filter(effects_full, common_name==this_species & data==this_dat)
  dat2plot <- filter(dat2plot, depth<depth_limit)
  dat2plot <- filter(dat2plot, (latitude>(taxa_table$min_lat-0.5))&(latitude<(taxa_table$max_lat+0.5)))
  
  #Restrict to region
  if(grepl("cc", this_dat)){
    dat2plot <- filter(dat2plot, region=="cc")
  }
  if(grepl("bc", this_dat)){
    dat2plot <- filter(dat2plot, region=="bc")
  }
  if(grepl("goa", this_dat)){
    dat2plot <- filter(dat2plot, region=="goa")
  }
  if(grepl("ebs", this_dat)){
    dat2plot <- filter(dat2plot, region=="ebs")
  }
  if(grepl("coastwide", this_dat)){
    dat2plot <- dat2plot
  }
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.5, alpha=0.1)+
    geom_point(filter(dat2plot,est_effect_prop>0),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.9)+
    #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
    ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
    #facet_wrap("year", ncol=5)+
    theme_minimal(base_size=12)+
    xlab("Longitude")+
    ylab("Latitude")+
    ggtitle(paste(this_species, this_dat, sep=" "))+
    theme(axis.text.x=element_blank(), legend.position="top")+
    scale_colour_viridis(name="proportional biomass \nreduction from oxygen", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish, option="plasma")
  
  ggsave(
    paste0("output/", output_folder, "/", "plots/points_full_region/map_",this_species,"_", this_dat,".png"),
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

##Annually
for(i in 1:nrow(bp_est4)) {
  this_species <- bp_est4$species[i]
  print(this_species)
  this_dat <- bp_est4$data[i]
  taxa_table <- filter(species_table, common_name==this_species)
  dat2plot <- filter(effects_full, common_name==this_species & data==this_dat)
  dat2plot <- filter(dat2plot, depth<depth_limit)
  dat2plot <- filter(dat2plot, (latitude>(taxa_table$min_lat-0.5))&(latitude<(taxa_table$max_lat+0.5)))
  
  #Restrict to region
  if(grepl("cc", this_dat)){
    dat2plot <- filter(dat2plot, region=="cc")
  }
  if(grepl("bc", this_dat)){
    dat2plot <- filter(dat2plot, region=="bc")
  }
  if(grepl("goa", this_dat)){
    dat2plot <- filter(dat2plot, region=="goa")
  }
  if(grepl("ebs", this_dat)){
    dat2plot <- filter(dat2plot, region=="ebs")
  }
  if(grepl("coastwide", this_dat)){
    dat2plot <- dat2plot
  }
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.5, alpha=0.1)+
    geom_point(filter(dat2plot,est_effect_prop>0),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.9)+
    #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    facet_wrap("year")+
    xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
    ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
    #facet_wrap("year", ncol=5)+
    theme_minimal(base_size=12)+
    xlab("Longitude")+
    ylab("Latitude")+
    ggtitle(paste(this_species, this_dat, sep=" "))+
    theme(axis.text.x=element_blank(), legend.position="top")+
    scale_colour_viridis(name="proportional biomass \n reduction from oxygen", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish, option="plasma")
  
  ggsave(
    paste0("output/", output_folder, "/", "plots/points_full_region_annual/map_",this_species,"_", this_dat,".png"),
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

###Species max depth by estimated effect--this doesn't really look that nice
#Combine depth
bp_est4_species <- left_join(bp_est4, species_table[,c("common_name", "depth")], by=c("species"="common_name"))
#Plot
ggplot(bp_est4_species, aes(x=depth, y=exp(bp_unscaled*slope)))+geom_point(aes(colour=species))
#Add species without oxygen limitation
species_O2 <- unique(bp_est4$species)
species_no_o2 <- filter(species_table, !common_name %in% species_O2)

species_no_o2$effect <- 0

ggplot(bp_est4_species, aes(x=depth, y=exp(bp_unscaled*slope)))+geom_point(aes(colour=species))+
geom_point(species_no_o2, mapping=aes(x=depth, y=effect), colour="black")

###Predict to grid
bc_grid <- readRDS("data/processed_data/o2/bc_predictions_grid.rds")
cc_grid <- readRDS("data/processed_data/o2/cc_predictions_grid.rds")

##Calculated estimated conditional effect in each grid
#Abundance data 
bp_est5 <- filter(bp_est4, grepl("iphc", data))
bp_est5 <- filter(bp_est5, species!="pacific halibut")
species_test <- unique(bp_est5$species)

#pick which ones to use
for(i in 1:length(species_test)){
  this_data <- filter(bp_est5, species==species_test[i])
  if(nrow(this_data)>1){
    this_data <- filter(this_data, grepl("coastwide", data))
  }
  dat.2.est <- bind_rows(cc_grid, bc_grid)
  this_species = this_data$species
  print(this_species)
  this_dat <- this_data$data
  this_model <- this_data$model
  best_model <- this_data$model
  slope.2.use <- this_data$slope
  breakpt.2.use <- this_data$bp_unscaled
  breakpt.2.use_se1 <- breakpt.2.use-this_data$bp_se_unscaled
  breakpt.2.use_se2 <- breakpt.2.use+this_data$bp_se_unscaled
  #Pull the data file 
  this_datframe <-try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
  
  #Calculate MI
  ###Calculate Metabolic index 
  ##Species parameters from Tim's paper
  #Read table
  taxa_table <- filter(species_table, common_name==this_species)
  taxa <- taxa_table$MI_Taxa
  mi_pars <- read.csv("data/taxa_table.csv")
  mi_pars$Group <- tolower(mi_pars$Group)
  mi_pars <- filter(mi_pars,Group==taxa)
  
  V <- mi_pars$logV
  n <- mi_pars$n
  Eo <- c(mi_pars$Eolow, mi_pars$Eo, mi_pars$Eohigh)
  Ao <- 1/exp(V)
  Eo <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") Eo[1] else if(best_model=="model4"|best_model=="model10"|best_model=="model14") Eo[2] else Eo[3]
  #Abbreviated equation
  #dat$mi1 <- dat$po2*exp(Eo1* dat$invtemp)
  #dat$mi2 <- dat$po2*exp(Eo2* dat$invtemp)
  #dat$mi2 <- dat$po2*exp(Eo2* dat$invtemp)
  fancy <- FALSE
  #Calculate MI for grid
  #Calc invtemp
  kelvin = 273.15
  boltz = 0.000086173324
  tref <- 15
  dat.2.est$invtemp <- (1 / boltz)  * ( 1 / (dat.2.est$temperature_C + 273.15) - 1 / (dat.2.est$temperature_C + 273.15))
  dat.2.est$mi <- calc_mi(Eo=Eo, po2=dat.2.est$po2, inv.temp=dat.2.est$invtemp, fancy=fancy)
  
  #MI from data for scaling
  mi.2.use <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") this_datframe$mi1 else if(best_model=="model4"|best_model=="model10"|best_model=="model14") this_datframe$mi2 else this_datframe$mi3
  mi_mean <- mean(mi.2.use, na.rm=TRUE)
  mi_sd <- sd(mi.2.use, na.rm=TRUE)
  #Scale MI calculated in grid to the fitting dataset
  dat.2.est$mi_s <- (dat.2.est$mi-mi_mean)/mi_sd
  mi_s <- dat.2.est$mi_s
  #Calculate MI effect
  #Calculate raw effect
  dat.2.est$est_effect_raw <- sapply(mi_s,breakpoint_calc, slope.2.use,breakpt.2.use)
  dat.2.est$est_effect_raw_se1<- sapply(mi_s,breakpoint_calc, slope.2.use,breakpt.2.use_se1)
  dat.2.est$est_effect_raw_se2 <- sapply(mi_s,breakpoint_calc, slope.2.use,breakpt.2.use_se2)
  
  #Exponentiate effect
  dat.2.est$est_effect_raw <- exp(dat.2.est$est_effect_raw)
  dat.2.est$est_effect_raw_se1 <- exp(dat.2.est$est_effect_raw_se1)
  dat.2.est$est_effect_raw_se2 <- exp(dat.2.est$est_effect_raw_se2)
  
  #Calculate max effect
  max_effect <- exp(breakpoint_calc(breakpt.2.use, slope.2.use, breakpt.2.use))
  max_effect_se1 <- exp(breakpoint_calc(breakpt.2.use_se1, slope.2.use, breakpt.2.use_se1))
  max_effect_se2 <- exp(breakpoint_calc(breakpt.2.use_se2, slope.2.use, breakpt.2.use_se2))
  dat.2.est$max_effect <- max_effect
  dat.2.est$max_effect_se1 <- max_effect_se1
  dat.2.est$max_effect_se2 <- max_effect_se2
  
  #Raw biomass reduction
  dat.2.est$biomass_reduction <- max_effect-dat.2.est$est_effect_raw
  dat.2.est$biomass_reduction_se1 <- max_effect_se1-dat.2.est$est_effect_raw_se1
  dat.2.est$biomass_reduction_se2 <- max_effect_se2-dat.2.est$est_effect_raw_se2
  
  #Proportional biomass reduction
  dat.2.est$est_effect_prop <- (max_effect-dat.2.est$est_effect_raw)/max_effect
  dat.2.est$est_effect_prop_se1 <- (max_effect_se1-dat.2.est$est_effect_raw_se1)/max_effect_se1
  dat.2.est$est_effect_prop_se2 <- (max_effect_se2-dat.2.est$est_effect_raw_se2)/max_effect_se2
  
  #Other columns
  dat.2.est$model <- this_model
  dat.2.est$data <- this_dat 
  dat.2.est$common_name <- this_species
  dat.2.est$depth_limit <- taxa_table$depth
  dat.2.est$min_lat <- taxa_table$min_lat
  dat.2.est$max_lat <- taxa_table$max_lat
  
  if(i==1){
    grid_abundance <- dat.2.est
  }
  if(i>0){
    grid_abundance <- bind_rows(grid_abundance, dat.2.est)
  }
}

#Biomass data 
bp_est6 <- filter(bp_est4, (!grepl("iphc", data)&species!="pacific halibut"))
phal <- filter(bp_est4,(grepl("iphc", data)&species=="pacific halibut"))
bp_est6 <- filter(bp_est6, !grepl("goa", data))
bp_est6 <- bind_rows(bp_est6, phal)
species_test <- unique(bp_est6$species)

for(i in 1:length(species_test)){
  this_data <- filter(bp_est6, species==species_test[i])
  if(nrow(this_data)>1){
    this_data <- filter(this_data, grepl("coastwide", data))
  }
  dat.2.est <- bind_rows(cc_grid, bc_grid)
  this_species = this_data$species
  print(this_species)
  this_dat <- this_data$data
  this_model <- this_data$model
  best_model <- this_data$model
  slope.2.use <- this_data$slope
  breakpt.2.use <- this_data$bp_unscaled
  breakpt.2.use_se1 <- breakpt.2.use-this_data$bp_se_unscaled
  breakpt.2.use_se2 <- breakpt.2.use+this_data$bp_se_unscaled
  #Pull the data file 
  this_datframe <-try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
  
  #Calculate MI
  ###Calculate Metabolic index 
  ##Species parameters from Tim's paper
  #Read table
  taxa_table <- filter(species_table, common_name==this_species)
  taxa <- taxa_table$MI_Taxa
  mi_pars <- read.csv("data/taxa_table.csv")
  mi_pars$Group <- tolower(mi_pars$Group)
  mi_pars <- filter(mi_pars,Group==taxa)
  
  V <- mi_pars$logV
  n <- mi_pars$n
  Eo <- c(mi_pars$Eolow, mi_pars$Eo, mi_pars$Eohigh)
  Ao <- 1/exp(V)
  Eo <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") Eo[1] else if(best_model=="model4"|best_model=="model10"|best_model=="model14") Eo[2] else Eo[3]
  #Abbreviated equation
  #dat$mi1 <- dat$po2*exp(Eo1* dat$invtemp)
  #dat$mi2 <- dat$po2*exp(Eo2* dat$invtemp)
  #dat$mi2 <- dat$po2*exp(Eo2* dat$invtemp)
  fancy <- FALSE
  #Calculate MI for grid
  #Calc invtemp
  kelvin = 273.15
  boltz = 0.000086173324
  tref <- 15
  dat.2.est$invtemp <- (1 / boltz)  * ( 1 / (dat.2.est$temperature_C + 273.15) - 1 / (dat.2.est$temperature_C + 273.15))
  dat.2.est$mi <- calc_mi(Eo=Eo, po2=dat.2.est$po2, inv.temp=dat.2.est$invtemp, fancy=fancy)
  
  #MI from data for scaling
  mi.2.use <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") this_datframe$mi1 else if(best_model=="model4"|best_model=="model10"|best_model=="model14") this_datframe$mi2 else this_datframe$mi3
  mi_mean <- mean(mi.2.use, na.rm=TRUE)
  mi_sd <- sd(mi.2.use, na.rm=TRUE)
  #Scale MI calculated in grid to the fitting dataset
  dat.2.est$mi_s <- (dat.2.est$mi-mi_mean)/mi_sd
  mi_s <- dat.2.est$mi_s
  #Calculate MI effect
  #Calculate raw effect
  dat.2.est$est_effect_raw <- sapply(mi_s,breakpoint_calc, slope.2.use,breakpt.2.use)
  dat.2.est$est_effect_raw_se1<- sapply(mi_s,breakpoint_calc, slope.2.use,breakpt.2.use_se1)
  dat.2.est$est_effect_raw_se2 <- sapply(mi_s,breakpoint_calc, slope.2.use,breakpt.2.use_se2)
  
  #Exponentiate effect
  dat.2.est$est_effect_raw <- exp(dat.2.est$est_effect_raw)
  dat.2.est$est_effect_raw_se1 <- exp(dat.2.est$est_effect_raw_se1)
  dat.2.est$est_effect_raw_se2 <- exp(dat.2.est$est_effect_raw_se2)
  
  #Calculate max effect
  max_effect <- exp(breakpoint_calc(breakpt.2.use, slope.2.use, breakpt.2.use))
  max_effect_se1 <- exp(breakpoint_calc(breakpt.2.use_se1, slope.2.use, breakpt.2.use_se1))
  max_effect_se2 <- exp(breakpoint_calc(breakpt.2.use_se2, slope.2.use, breakpt.2.use_se2))
  dat.2.est$max_effect <- max_effect
  dat.2.est$max_effect_se1 <- max_effect_se1
  dat.2.est$max_effect_se2 <- max_effect_se2
  
  #Raw biomass reduction
  dat.2.est$biomass_reduction <- max_effect-dat.2.est$est_effect_raw
  dat.2.est$biomass_reduction_se1 <- max_effect_se1-dat.2.est$est_effect_raw_se1
  dat.2.est$biomass_reduction_se2 <- max_effect_se2-dat.2.est$est_effect_raw_se2
  
  #Proportional biomass reduction
  dat.2.est$est_effect_prop <- (max_effect-dat.2.est$est_effect_raw)/max_effect
  dat.2.est$est_effect_prop_se1 <- (max_effect_se1-dat.2.est$est_effect_raw_se1)/max_effect_se1
  dat.2.est$est_effect_prop_se2 <- (max_effect_se2-dat.2.est$est_effect_raw_se2)/max_effect_se2
  
  #Other columns
  dat.2.est$model <- this_model
  dat.2.est$data <- this_dat 
  dat.2.est$common_name <- this_species
  dat.2.est$depth_limit <- taxa_table$depth
  dat.2.est$min_lat <- taxa_table$min_lat
  dat.2.est$max_lat <- taxa_table$max_lat
  
  if(i==1){
    grid_biomass <- dat.2.est
  }
  if(i>0){
    grid_biomass <- bind_rows(grid_biomass, dat.2.est)
  }
}

#Plot just one year of data, all species on same plot
#Combine biomass and abundance
dat2plot <- filter(grid_biomass, depth_m<(depth_limit+200))
dat2plot <- filter(dat2plot, (latitude>(taxa_table$min_lat-0.5))&(latitude<(taxa_table$max_lat+0.5)))
dat2plot <- filter(dat2plot, year==2021)

dat2plot2 <- filter(grid_abundance, depth_m<(depth_limit+200))
dat2plot2 <- filter(dat2plot2, (latitude>(taxa_table$min_lat-0.5))&(latitude<(taxa_table$max_lat+0.5)))
dat2plot2 <- filter(dat2plot2, year==2021)

dat2plot <- bind_rows(dat2plot, dat2plot2)

ggplot(us_coast_proj) + geom_sf() +
  geom_point(dat2plot,mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
  ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
  facet_wrap("common_name", ncol=5)+
  theme_minimal(base_size=12)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position="top")+
  theme(axis.text.x=element_blank())+
  scale_colour_viridis(name="proportional density \n reduction from oxygen", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish, option="plasma")

ggsave(
  paste0("output/", output_folder, "/", "plots/grid_combined_2021.png"),
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

#Plot just one year of data, all species on same plot, zoom to WA
#Combine biomass and abundance
dat2plot <- filter(grid_biomass, depth_m<(depth_limit+200))
dat2plot <- filter(dat2plot, (latitude>(taxa_table$min_lat-0.5))&(latitude<(taxa_table$max_lat+0.5)))
dat2plot <- filter(dat2plot, year==2021)

dat2plot2 <- filter(grid_abundance, depth_m<(depth_limit+200))
dat2plot2 <- filter(dat2plot2, (latitude>(taxa_table$min_lat-0.5))&(latitude<(taxa_table$max_lat+0.5)))
dat2plot2 <- filter(dat2plot2, year==2021)

dat2plot <- bind_rows(dat2plot, dat2plot2)

ggplot(us_coast_proj) + geom_sf() +
  geom_point(dat2plot,mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
  ylim(100*1000, 300*1000)+
  facet_wrap("common_name", ncol=5)+
  theme_minimal(base_size=12)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position="top")+
  theme(axis.text.x=element_blank())+
  scale_colour_viridis(name="proportional density \n reduction from oxygen", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish, option="plasma")

ggsave(
  paste0("output/", output_folder, "/", "plots/grid_combined_2021_WA.png"),
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

#Plot all years, each species separately
to_do <- unique(grid_abundance[c("common_name", "data")])
                
for(i in 1:nrow(to_do)) {
  this_species <- to_do$common_name[i]
  print(this_species)
  this_dat <- to_do$data[i]
  taxa_table <- filter(species_table, common_name==this_species)
  dat2plot <- filter(grid_abundance, common_name==this_species & data==this_dat)
  dat2plot <- filter(dat2plot, depth_m<(depth_limit+200))
  dat2plot$latitude <- dat2plot$latitute
  dat2plot <- filter(dat2plot, (latitude>(taxa_table$min_lat-0.5))&(latitude<(taxa_table$max_lat+0.5)))
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(dat2plot,mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
    ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
    facet_wrap("year", ncol=5)+
    theme_minimal(base_size=16)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(legend.position=c(0.8,0.1))+
    #ggtitle(paste(this_species, this_dat, sep=" "))+
    theme(axis.text.x=element_blank())+
    scale_colour_viridis(name="proportional density \n reduction from oxygen", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish, option="plasma")
  #scale_colour_viridis(name="proportional biomass reduction \n from oxygen", limits=c(0,1), breaks=c(0,0.25, 0.5, 0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish)
  
  ggsave(
    paste0("output/", output_folder, "/", "plots/grid_annual/map_",this_species,"_", this_dat,".png"),
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

to_do <- unique(grid_biomass[c("common_name", "data")])

for(i in 1:nrow(to_do)) {
  this_species <- to_do$common_name[i]
  print(this_species)
  this_dat <- to_do$data[i]
  taxa_table <- filter(species_table, common_name==this_species)
  dat2plot <- filter(grid_biomass, common_name==this_species & data==this_dat)
  dat2plot <- filter(dat2plot, depth_m<(depth_limit+200))
  dat2plot$latitude <- dat2plot$latitute
  dat2plot <- filter(dat2plot, (latitude>(taxa_table$min_lat-0.5))&(latitude<(taxa_table$max_lat+0.5)))
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(dat2plot,mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
    ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
    facet_wrap("year", ncol=5)+
    theme_minimal(base_size=16)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(legend.position=c(0.8,0.1))+
    #ggtitle(paste(this_species, this_dat, sep=" "))+
    theme(axis.text.x=element_blank())+
    scale_colour_viridis(name="proportional density \n reduction from oxygen", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish, option="plasma")
  #scale_colour_viridis(name="proportional biomass reduction \n from oxygen", limits=c(0,1), breaks=c(0,0.25, 0.5, 0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish)
  
  ggsave(
    paste0("output/", output_folder, "/", "plots/grid_annual/map_",this_species,"_", this_dat,".png"),
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

##Plot just one year of data, separate species
#Add

##Calculate total area of reduced habitat (and proportion of habitat)
#Combine
grid_full <- bind_rows(grid_abundance, grid_biomass)
test <- filter(grid_full, depth_m<(depth_limit+200))
test$latitude <- test$latitute
test <- filter(test, (latitude>(min_lat-0.5))&(latitude<max_lat+0.5))

test2 <- filter(test, est_effect_prop>0)
test3 <- filter(test, est_effect_prop_se1>0)
test4 <- filter(test, est_effect_prop_se2>0)

effects_full_annual_summary <- test2 %>%
  group_by(year, region, common_name) %>%
  summarise(area_sum=sum(area)) %>%
  ungroup()
effects_full_annual_summary2 <- test3 %>%
  group_by(year, region, common_name) %>%
  summarise(area_sum=sum(area)) %>%
  ungroup()
effects_full_annual_summary3 <- test4 %>%
  group_by(year, region, common_name) %>%
  summarise(area_sum=sum(area)) %>%
  ungroup()

effects_full_summary <- full_join(effects_full_annual_summary, effects_full_annual_summary2, by=c("year", "region", "common_name"))
effects_full_summary$se1 <- effects_full_summary$area_sum.y
effects_full_summary$area_sum_med <- effects_full_summary$area_sum.x
effects_full_summary$area_sum.y <- NULL
effects_full_summary$area_sum.x <- NULL
effects_full_summary <- full_join(effects_full_summary, effects_full_annual_summary3, by=c("year", "region", "common_name"))
effects_full_summary$se2 <- effects_full_summary$area_sum
effects_full_summary$area_sum <- NULL

#Total area
cc_total <- filter(cc_grid, year==2023)
sum_cc <- sum(cc_total$area)

bc_total <- filter(bc_grid, year==2023)
sum_bc <- sum(bc_total$area)

effects_full_summary$total <- ifelse(effects_full_summary$region=="cc", sum_cc, sum_bc)

#plot
ggplot(effects_full_summary, aes(x=year, y=area_sum_med/total))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1/total, ymax=se2/total, fill=region), alpha=0.2)+
  facet_wrap("common_name", scales="free_y")+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"))+
  scale_fill_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"))+
  xlab("Year") +
  ylab("Proportion of area below oxygen threshold")

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_area_below_grid.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 8.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

#plot
ggplot(effects_full_summary, aes(x=year, y=area_sum_med))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  facet_wrap("common_name", scales="free_y")+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"))+
  scale_fill_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"))+
  xlab("Year") +
  ylab("Total area below oxygen threshold (ha)")

ggsave(
  paste0("output/", output_folder, "/", "plots/raw_area_below_grid.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 8.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)
##Calculate total biomass reduced
test2$reduced <- test2$biomass_reduction*test2$area
test2$reduced2 <- test2$biomass_reduction_se1*test2$area
test2$reduced3 <- test2$biomass_reduction_se2*test2$area

effects_full_annual_summary <- test2 %>%
  group_by(year, region, common_name) %>%
  summarise(biomass=sum(reduced), se1=sum(reduced2), se2=sum(reduced3)) %>%
  ungroup()

#plot
ggplot(effects_full_annual_summary, aes(x=year, y=biomass))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  facet_wrap("common_name", scales="free_y")+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"))+
  scale_fill_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"))+
  xlab("Year") +
  ylab("Total Reduction in Density")

ggsave(
  paste0("output/", output_folder, "/", "plots/density_reduced_anually_grid_raw.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 8.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

#Scaled to maximum
maxs <- effects_full_annual_summary %>%
  group_by(common_name, region) %>%
  summarise(max=max(se2)) %>%
  ungroup()
effects_full_annual_summary <- left_join(effects_full_annual_summary, maxs, by=c("common_name", "region"))

#plot
ggplot(effects_full_annual_summary, aes(x=year, y=biomass/max))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1/max, ymax=se2/max, fill=region), alpha=0.2)+
  facet_wrap("common_name", scales="free_y")+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"))+
  scale_fill_manual(values=c("#F8766D","#7CAE00", "#C77CFF", "#00BFC4"))+
  xlab("Year") +
  ylab("Total Reduction in Density")

ggsave(
  paste0("output/", output_folder, "/", "plots/density_reduced_anually_grid_scaled.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 8.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)
#Washington coast, bin by latitude
test2$reduced <- test2$biomass_reduction*test2$area
test2$reduced2 <- test2$biomass_reduction_se1*test2$area
test2$reduced3 <- test2$biomass_reduction_se2*test2$area

effects_full_annual_summary <- test2 %>%
  group_by(year, region, common_name) %>%
  summarise(biomass=sum(reduced), se1=sum(reduced2), se2=sum(reduced3)) %>%
  ungroup()

##Heat map--test just for WA coast w/ data
#Group grid into quadrants
plot_all<-ggplot(comp4) + geom_tile(aes(x = Quad, y = Date, fill = q50)) + 
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red", 
                       midpoint = 0.5, limits = c(0, 1)) + 
  scale_x_continuous(breaks = seq(0, 500, by = 50), expand = c(0, 0)) + 
  scale_y_date(expand = c(0,0)) + 
  facet_grid(~region, scales = "free_x", space = "free_x", 
             labeller = labeller(region = label_wrap_gen(4))) + 
  theme_minimal() + 
  labs(y = "Year", x = "Whale Museum Quadrant", fill = "Average relative \n SRKW probability") + 
  theme(panel.spacing.x = unit(0, "lines"), 
        strip.background = element_rect(fill = "lavender"), 
        legend.position = "bottom", 
        strip.text = element_text(size = 8))

#Predict full biomass to get proportion reduced from O2?
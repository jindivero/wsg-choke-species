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
library(viridis)
library(ggh4x)
library(ggnewscale)
library(sf)
library(mapview)
library(openxlsx)
library(parallel)
library(ggpubr)

set.seed(9876)

setwd("~/Dropbox/GitHub/wsg-choke-species")

#ggplot themes
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

#Load functions
source("code/helper_funs.R")

#List of regions, species
dat_names <- c("cc", "bc", "goa", "ebs", "coastwide", "cc _iphc", "bc _iphc", "goa _iphc", "ebs _iphc", "coastwide _iphc")
species <- read_excel("data/species_table.xlsx")
species$common_name <- tolower(species$common_name)
species$scientific_name <- tolower(species$scientific_name)
species_table <- species
species <- unique(species$common_name)
species <- unique(species_table$common_name)
species_iphc <- c("sablefish", "pacific cod", "pacific halibut", "yelloweye rockfish", "longnose skate", "big skate", "spiny dogfish", "rougheye rockfish")

#Output folder
output_folder <- "region_comp"

#Metabolic index models to use
mi_models2use <- c("model13", "model14", "model15")

#Taxa lookup
taxa <- read_excel("data/species_table.xlsx")
taxa$MI_Taxa <- tolower(taxa$MI_Taxa)
taxa$common_name <- tolower(taxa$common_name)

#Load AIC table for model output
aic <- as.data.frame(read_excel(paste0("output/",output_folder, "/aic_table.xlsx")))
#Make first 6 columns in AIC numeric
aic[,1:(ncol(aic)-5)] <- sapply(aic[,1:(ncol(aic)-5)], as.numeric)

##Set up mapping
map_data <- rnaturalearth::ne_countries(scale = "large",
                                        returnclass = "sf",
                                        continent = "North America")

us_coast_proj <- sf::st_transform(map_data, crs = 32610)

###Conditional effect plot
#Manually calculating
##Create dataset used by all
#Load all fish data
files <- list.files(path = "data/processed_data/fish2", pattern = ".rds", full.names=T)
dat <- map(files,readRDS)
dat <- bind_rows(dat)

#Remove any rows with necessary data missing
dat <- dat %>%
  drop_na(depth, mi1, temperature_C, salinity_psu, X, Y, year)

#Remove weird depths
dat <- filter(dat, depth>0)

#Remove oxygen outliers
dat <- filter(dat, O2_umolkg<1500)

##Make prediction grid for conditional effects
#Sequence of oxygen values
dat_pred <- as.data.frame(seq(0,30, length.out=100))
colnames(dat_pred) <- "po2"

#Add other columns for conditional effect fitting
dat_pred$temperature_C <- 12
#Calc invtemp
kelvin = 273.15
boltz = 0.000086173324
tref <- 12
dat_pred$invtemp <- (1 / boltz)  * ( 1 / (tref + 273.15) - 1 / (tref + 273.15))
dat_pred$pred_id <- 1:100

#----Conditional Effect: Manually w/ Monte Carlo parameter draws
for(i in 1:length(species)){
  this_species <- species[i]
  print(this_species)
  preds2 <- dat_pred
  #Calculate metabolic index from correct taxa parameters
  #Pull correct 
  mi_pars <- read.csv("data/taxa_table.csv")
  mi_pars$Group <- tolower(mi_pars$Group)
  this_taxa <- taxa$MI_Taxa[taxa$common_name==this_species]
  mi_pars <- filter(mi_pars,Group==this_taxa)
  
  Eo <- c(mi_pars$Eolow, mi_pars$Eo, mi_pars$Eohigh)
  
  #Calculate metabolic index for species
  preds2$mi1 <- calc_mi(po2=preds2$po2, inv.temp=preds2$invtemp, Eo=Eo[1],fancy=F)
  preds2$mi2 <- calc_mi(po2=preds2$po2, inv.temp=preds2$invtemp, Eo=Eo[2],fancy=F)
  preds2$mi3 <- calc_mi(po2=preds2$po2, inv.temp=preds2$invtemp, Eo=Eo[3],fancy=F)
  
  #List of data types to pull
  this_aic <- filter(aic, species==this_species)
  datnames <- unique(this_aic$`data type`)
  #For each data type, pull all of the models that did fit
  for (j in 1:length(datnames)){
    preds <- preds2
    this_dat <- datnames[j]
    print(this_dat)
    #Pull the data file for scaling
    this_datframe <- try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
    #Mean and sd for scaling temp
    mean_temp <- mean(this_datframe$temperature_C, na.rm=T)
    sd_temp <- sd(this_datframe$temperature_C, na.rm=T)
    #Mean and sd for scaling MIs
    mean_mi1 <- mean(this_datframe$mi1, na.rm=T)
    sd_mi1 <- sd(this_datframe$mi1, na.rm=T)
    mean_mi2 <- mean(this_datframe$mi2, na.rm=T)
    sd_mi2 <- sd(this_datframe$mi2, na.rm=T)
    mean_mi3 <- mean(this_datframe$mi3, na.rm=T)
    sd_mi3 <- sd(this_datframe$mi3, na.rm=T)
    #Scale prediction data
    preds$temp_scaled <- (preds$temperature_C-mean_temp)/sd_temp
    preds$temp_scaled2 <- (preds$temp_scaled)^2
    preds$mi1_s <- (preds$mi1-mean_mi1)/sd_mi1
    preds$mi2_s <- (preds$mi2-mean_mi2)/sd_mi2
    preds$mi3_s <- (preds$mi3-mean_mi3)/sd_mi3
    
    #Filter AIC table to the datatype to figure out which models to pull
    this_aic <- filter(aic, species==this_species, `data type`==this_dat)
    #Just the first 5 columns
    this_aic <- this_aic[,1:(ncol(this_aic)-5)]
    #Only the columns that are not NAs (which means they didn't pass sanity checks)
    this_aic <- this_aic[,colSums(is.na(this_aic))<nrow(this_aic)]
    
    #Get list of these columns (these are the model fits to pull for model averaging)
    models <- colnames(this_aic)
    if(length(models)>0){
      #Pull the model output files
      model_fits <- list()
      for(h in 1:length(models)){
        fit <- try(readRDS(file = paste0("output/", output_folder,"/", this_species,"_",this_dat,"_", models[h], ".rds")))
        model_fits[[h]] <- fit
      }
        #Get model weights
        aics <- list()
        for (k in 1:length(models)){
          aic_models <- AIC(model_fits[[k]])
          aics[[k]] <- aic_models
        }
        aics <- unlist(aics)
        delta_aic <- aics - min(aics)
        weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
        set.seed(459384) # for reproducibility and consistency
        for(g in 1:length(model_fits)){
          if(models[g] %in% mi_models2use){
            fit <- model_fits[[g]]
            pars <- as.data.frame(tidy(fit, effects="fixed", conf.int=T))
            slope <- filter(pars, grepl("slope", term))
            thresh <- filter(pars,grepl("breakpt", term))
            #Pull 100 threshold estimates or slope estimates
            bp_pars <- rnorm(mean=thresh$estimate, sd=thresh$std.error, n=100)
            slope_pars <- rnorm(mean=slope$est, sd=slope$std.error, n=100)
            slope <- slope$estimate
            #Calculate the breakpoint effect
            if(models[g]=="model13"){
              mi_s <- preds$mi1_s
            }
            if(models[g]=="model14"){
              mi_s <- preds$mi2_s
            }
            if(models[g]=="model15"){
              mi_s <- preds$mi3_s
            }
            for(l in 1:length(bp_pars)){
            pred <- as.data.frame(sapply(mi_s,breakpoint_calc, slope_pars[l],bp_pars[l]))
            colnames(pred) <- "pred"
            pred$weight <- weights[g]
            pred$sim <- l
            pred$pred_id <-1:100
            #Combine together
            if(l==1){
              pred_model <- pred
            }
            if(l>1){
              pred_model <- bind_rows(pred_model, pred)
            }
            }
      
          } else {
            #For models without pO2 term, conditional effect will be 0--add 100 iterations of this for each prediction point
            pred_model <- as.data.frame(matrix(-1000000, nrow=length(preds$po2)*100))
            colnames(pred_model) <- "pred"
            pred_model$weight <- weights[g]
            pred_model$pred_id <- rep(1:100, 100)
            pred_model$sim <- rep(1:100, each=100)
          }
          if(g==1){
            pred_all <- pred_model
          }
          if(g>1){
            pred_all <- bind_rows(pred_all, pred_model)
          }
        }     
      #Calculate weighted average
      ens_preds <- pred_all %>% 
        group_by(pred_id, sim) %>% 
        summarise(weighted_mean=weighted.mean(exp(pred),weight))  %>% 
        ungroup()
      
      #Ensemble w/ SD
      ens_preds2 <- ens_preds %>% 
        group_by(pred_id) %>% 
        summarise(ensemble_mean=mean(weighted_mean, na.rm=T), ensemble_sd=sd(weighted_mean, na.rm=T)) %>%
        ungroup()
      ens_preds2 <- ens_preds2 %>% 
        mutate(ensemble_lower=(ensemble_mean-ensemble_sd), ensemble_upper=(ensemble_mean+ensemble_sd))%>% 
        ungroup()
      
      #Scaled
      ens_preds3 <- ens_preds2 %>%
        mutate(ensemble_mean_sc= (ensemble_mean/max(ensemble_upper, na.rm=T)),
               ensemble_mean_lower_sc = (ensemble_lower)/max(ensemble_upper, na.rm=T),
               ensemble_mean_upper_sc = (ensemble_upper)/max(ensemble_upper, na.rm=T))
      #Add prediction dataframe back
      preds5 <- left_join(preds, ens_preds3, by="pred_id")
      preds5$species <- this_species
      preds5$data <- this_dat
      if(!grepl("iphc", this_dat)){
      preds5$region <- this_dat
      preds5$data_type <- "bottom trawl only"
      }
      if(grepl("iphc", this_dat)){
      preds5$region <- gsub(" _iphc", "", this_dat)
      preds5$data_type <- "bottom trawl & IPHC"
      }
      if(j==1){
        all_preds <- preds5}
      if(j>1){
        all_preds <- bind_rows(all_preds, preds5)
      }
    }
  }
      
  if(i==1){
    all_preds2 <- all_preds}
  if(i>1){
    all_preds2 <- bind_rows(all_preds2, all_preds)
  }
  }
  
#Save all_preds2
saveRDS(all_preds2, file = paste0("output/", output_folder, "/conditional_effects_data_all.rds"))

##Plot all with SE
#Set region
all_preds2$region <- factor(all_preds2$region, levels=c("ebs", "goa", "bc", "cc", "coastwide"))
labs <- c("Eastern Bering Sea", "Gulf of Alaska", "British Columbia", "California Current", "Coastwide")
names(labs) <- c("ebs", "goa", "bc", "cc", "coastwide")

ggplot(all_preds2, aes(x=po2, y=ensemble_mean_sc))+
  geom_line(aes(colour=region, linetype=data_type))+
 geom_ribbon(aes(ymin=ensemble_mean_lower_sc, ymax=ensemble_mean_upper_sc, fill=region), alpha=0.2)+
  facet_wrap("species", ncol=4, scales="free_y", labeller=labeller(species=label_wrap_gen(15)))+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE, guide="none")+
  scale_fill_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE, labels=labs )+
  scale_linetype_manual(values=c("dashed", "solid"))+
  guides(fill = guide_legend(nrow = 3, override.aes = list(alpha = 1)), linetype=guide_legend(nrow=2))+
  xlab("Dissolved Oxygen (kPa)")+
  ylab("Effect on Fish Density")

ggsave(
  paste0("output/", output_folder, "/plots_final/cond_effect_ensemble_all.png"),
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

#Filtering protocol step 1:  aic to only if a 0 in any of columns named in mi_models2use
aic_full <- aic
aic <- filter(aic, aic[,mi_models2use[1]]==0 | aic[,mi_models2use[2]]==0 | aic[,mi_models2use[3]]==0)
species_o2 <- unique(aic$species)

##Filtering protocol step 2: check reasonableness of breakpoint estimates
#Pull 100 draws of parameter estimates (slope and breakpoint) and calculate ensemble
for(i in 1:length(species_o2)){
  this_species <- species_o2[i]
  print(this_species)
  #Calculate metabolic index from correct taxa parameters
  #Pull correct MI values
  mi_pars <- read.csv("data/taxa_table.csv")
  mi_pars$Group <- tolower(mi_pars$Group)
  this_taxa <- taxa$MI_Taxa[taxa$common_name==this_species]
  mi_pars <- filter(mi_pars,Group==this_taxa)
  
  Eo <- c(mi_pars$Eolow, mi_pars$Eo, mi_pars$Eohigh)
  
  #List of data types to pull
  this_aic <- filter(aic, species==this_species)
  datnames <- unique(this_aic$`data type`)
  #For each data type, pull all of the models that did fit
  for (j in 1:length(datnames)){
    this_dat <- datnames[j]
    print(this_dat)
    #Pull the data file for scaling
    this_datframe <- try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
    #Mean and sd for unscaling MIs
    mean_mi1 <- mean(this_datframe$mi1, na.rm=T)
    sd_mi1 <- sd(this_datframe$mi1, na.rm=T)
    mean_mi2 <- mean(this_datframe$mi2, na.rm=T)
    sd_mi2 <- sd(this_datframe$mi2, na.rm=T)
    mean_mi3 <- mean(this_datframe$mi3, na.rm=T)
    sd_mi3 <- sd(this_datframe$mi3, na.rm=T)
    
    #Filter AIC table to the datatype to figure out which models to pull
    this_aic <- filter(aic, species==this_species, `data type`==this_dat)
    #Just the first 5 columns
    this_aic <- this_aic[,1:(ncol(this_aic)-5)]
    #Only the columns that are not NAs (which means they didn't pass sanity checks)
    this_aic <- this_aic[,colSums(is.na(this_aic))<nrow(this_aic)]
    
    #Get list of these columns (these are the model fits to pull for model averaging)
    models <- colnames(this_aic)
    #Filter to just the MI models
    models <- models[models %in% mi_models2use]
    if(length(models)>0){
      #Pull the model output files
      model_fits <- list()
      for(h in 1:length(models)){
        fit <- try(readRDS(file = paste0("output/", output_folder,"/", this_species,"_",this_dat,"_", models[h], ".rds")))
        model_fits[[h]] <- fit
      }
      #Get model weights
      aics <- list()
      for (k in 1:length(models)){
        aic_models <- AIC(model_fits[[k]])
        aics[[k]] <- aic_models
      }
      aics <- unlist(aics)
      delta_aic <- aics - min(aics)
      weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
      set.seed(459384) # for reproducibility and consistency
      for(g in 1:length(model_fits)){
       # if(models[g] %in% mi_models2use){
          fit <- model_fits[[g]]
          pars <- as.data.frame(tidy(fit, effects="fixed", conf.int=T))
          thresh <- filter(pars,grepl("breakpt", term))
          slope <- filter(pars, grepl("slope", term))
          #Pull 100 threshold estimates or slope estimates
          bp_par <- rnorm(mean=thresh$estimate, sd=thresh$std.error, n=100)
          slope_est <- exp(rnorm(mean=slope$estimate, sd=slope$std.error, n=100))
          pars <- data.frame(bp_par=bp_par, slope_est=slope_est)
          if(models[g]=="model13"){
           mean_mi <- mean_mi1
           sd_mi <- sd_mi1
          }
          if(models[g]=="model14"){
            mean_mi <- mean_mi2
            sd_mi <- sd_mi2
          }
          if(models[g]=="model15"){
            mean_mi <- mean_mi3
            sd_mi <- sd_mi3
          }
          pars$bp_pars <- (pars$bp_par*sd_mi)+mean_mi
          pars$sim <- 1:100
          pars$weight <- weights[g]
        # } else {
        #   #For models without pO2 term, conditional effect will be 0--add 100 iterations of this for each prediction point
        #   pars <- as.data.frame(matrix(0, nrow=100))
        #   colnames(pars) <- "bp_pars"
        #   pars$slope_pars <- 0
        #   pars$weight <- weights[g]
        #   pars$sim <- 1:100
       # }
        if(g==1){
          pars_all <- pars
        }
        if(g>1){
          pars_all <- bind_rows(pars_all, pars)
        }
      }
     # }     
      #Calculate weighted average
      pars_weighted <- pars_all %>% 
        group_by(sim) %>% 
        summarise(weighted_mean_bp=weighted.mean(bp_pars,weight, na.rm=T), weighted_mean_slope=weighted.mean(slope_est, weight, na.rm=T))  %>% 
        ungroup()
      
      #Ensemble w/ SD
      pars_ensemble <- pars_weighted %>% 
        summarise(bp_ensemble_mean=mean(weighted_mean_bp, na.rm=T), bp_ensemble_sd=sd(weighted_mean_bp, na.rm=T), slope_ensemble_mean=mean(weighted_mean_slope, na.rm=T), slope_ensemble_sd=sd(weighted_mean_slope, na.rm=T) ) %>%
        ungroup()
      pars_ensemble <- pars_ensemble %>% 
        mutate(ensemble_lower=(bp_ensemble_mean-bp_ensemble_sd), ensemble_upper=(bp_ensemble_mean+bp_ensemble_sd), slope_lower=(slope_ensemble_mean-slope_ensemble_sd), slope_upper=(slope_ensemble_mean+slope_ensemble_sd))%>% 
        ungroup()
      
      #Add other info
      pars_ensemble$species <- this_species
      pars_ensemble$data <- this_dat
      pars_ensemble$max_mi3 <- max(this_datframe$mi3, na.rm=T)
      pars_ensemble$max_mi2 <- max(this_datframe$mi2, na.rm=T)
      pars_ensemble$max_mi1 <- max(this_datframe$mi1, na.rm=T)
      
      if(!grepl("iphc", this_dat)){
        pars_ensemble$region <- this_dat
        pars_ensemble$data_type <- "bottom trawl only"
      }
      if(grepl("iphc", this_dat)){
        pars_ensemble$region <- gsub(" _iphc", "", this_dat)
        pars_ensemble$data_type <- "bottom trawl & IPHC"
      }
      if(j==1){
        all_preds <- pars_ensemble}
      if(j>1){
        all_preds <- bind_rows(all_preds, pars_ensemble)
      }
    }
  }

  if(i==1){
    bp_est <- all_preds}
  if(i>1){
    bp_est<- bind_rows(bp_est, all_preds)
  }
}

##Filter out species that don't have reasonable estimates
bp_est <- unique(bp_est)
bp_est <- filter(bp_est, (bp_ensemble_mean<max_mi1 & ensemble_lower>0 & bp_ensemble_sd<10 & slope_lower>0))

#Creat ID
bp_est$id <- paste(bp_est$species, bp_est$data, sep="_")
all_preds2$id <- paste(all_preds2$species, all_preds2$data, sep="_")

#conditional effects not negative or NA standard errors
all_preds3 <- all_preds2 %>% filter(id %in% bp_est$id)
ids2use <- all_preds3 %>% 
  group_by(id) %>%
  summarize(min=min(ensemble_lower))  %>%
  filter(!is.na(min)&min>0)
bp_est <- filter(bp_est, id %in% ids2use$id)

#Save
saveRDS(bp_est, file = paste0("output/", output_folder, "/breakpoint_estimates.rds"))

#Save a streamlined csv for a table
bp_save <- select(bp_est, species, region, data_type, bp_ensemble_mean, bp_ensemble_sd, slope_ensemble_mean,slope_ensemble_sd)
write.csv(bp_save, file = paste0("output/", output_folder, "/breakpoint_estimates.csv"), row.names=F)

#Combine with AIC and save
aic_full$data_type <- ifelse(grepl("iphc", aic_full$`data type`), "bottom trawl & IPHC", "bottom trawl only")
aic_full$region <- case_when(
  grepl("ebs", aic_full$`data type`) ~ "ebs",
  grepl("goa", aic_full$`data type`) ~ "goa",
  grepl("bc", aic_full$`data type`) ~ "bc",
  grepl("cc", aic_full$`data type`) ~ "cc",
  grepl("coastwide", aic_full$`data type`) ~ "coastwide"
)
bp_aic <- left_join(aic_full, bp_save, by=c("species", "data_type", "region"))
bp_aic$"data type" <- NULL
#Round model7-model15, bp_ensemble_mean through slope_ensemble_sd
bp_aic[,c("model7", "model8", "model13", "model14", "model15", "bp_ensemble_mean", "bp_ensemble_sd", "slope_ensemble_mean","slope_ensemble_sd")] <- round(bp_aic[,c("model7", "model8", "model13", "model14", "model15", "bp_ensemble_mean", "bp_ensemble_sd", "slope_ensemble_mean","slope_ensemble_sd")], 3)
#Replace NA with --
bp_aic[is.na(bp_aic)] <- "--"
bp_aic$N_obs <- bp_aic$'N obs'
bp_aic$N_region <- bp_aic$'N region'
bp_aic$N_years <- bp_aic$'N years'
#Reorder
bp_aic <- bp_aic %>%
  select(species, region, data_type, model7, model8, model13, model14, model15, bp_ensemble_mean, bp_ensemble_sd, slope_ensemble_mean, slope_ensemble_sd, N_obs, N_region, N_years)

#save as excel file
write.xlsx(bp_aic, file = paste0("output/", output_folder, "/breakpoint_estimates_aic.xlsx"), rowNames=F)

#Plot, filtered out
ggplot(filter(all_preds2, id %in% bp_est$id), aes(x=po2, y=ensemble_mean_sc))+
  geom_line(aes(colour=region, linetype=data_type))+
  geom_ribbon(aes(ymin=ensemble_mean_lower_sc, ymax=ensemble_mean_upper_sc, fill=region), alpha=0.2)+
  facet_wrap("species", ncol=4, scales="free_y", labeller=labeller(species=label_wrap_gen(15)))+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE, labels=labs)+
  scale_fill_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE, labels=labs)+
  scale_linetype_manual(values=c("dashed", "solid"))+
  guides(fill = guide_legend(nrow = 2, labels=labs), color=guide_legend(nrow=2, labels=labs, override.aes=list(size=4)), linetype=guide_legend(nrow=2))+
  xlab("Oxygen (kPa) at 12 C")+
  ylab("Effect on Fish Density")

ggsave(
  paste0("output/", output_folder, "/plots_final/species_filtered_alldat.png"),
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

#Plot line plot
#Set region
bp_est$region <- factor(bp_est$region, levels=c("ebs", "goa", "bc", "cc", "coastwide"))
labs <- c("Eastern Bering Sea", "Gulf of Alaska", "British Columbia", "California Current", "Coastwide")
names(labs) <- c("ebs", "goa", "bc", "cc", "coastwide")

ggplot(bp_est, aes(y=species, x=bp_ensemble_mean, colour=region))+
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=region, shape=data_type), size=3, position=ggstance::position_dodgev(height=0.4))+
  #Can add back shape
  geom_linerange(aes(xmin = ensemble_lower, xmax = ensemble_upper, colour=region),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=12))+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  #xlim(-1,10)+
  scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE)+
  guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
  theme(legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(0.25, "lines"), #Minimize legend space
        panel.spacing = unit(5, "lines"))+ #Make more space between species
  #scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  xlab("Temperature-Corrected Oxygen Breakpoint Estimate (kPa)")+
  ylab("")
#geom_vline(xintercept=0, linetype="dashed")+
# geom_vline(xintercept=1, linetype="dashed")

ggsave(
  paste0("output/", output_folder, "/plots_final/breakpoint_estimates.png"),
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


##Plot only one data type per species
#IPHC if available, coastwide if available
for(i in 1:length(unique(bp_est$species))){
  this_dat <- filter(bp_est, species==unique(bp_est$species)[i])
  this_species <- unique(bp_est$species)[i]
  if(any(grepl("coastwide", this_dat$data))){
    this_dat <- filter(this_dat, grepl("coastwide", data))
    if(this_species %in% species_iphc){
      this_dat <- filter(this_dat, grepl("iphc", data))
    } 
  }
  if(i==1){
    bp_est2 <- this_dat
  } else {
    bp_est2 <- bind_rows(bp_est2, this_dat)
  }
}

saveRDS(bp_est2, file = paste0("output/", output_folder, "/breakpoint_estimates_filtered.rds"))

#Plot
ggplot(filter(all_preds2, id %in% bp_est2$id), aes(x=po2, y=ensemble_mean_sc))+
  geom_line(aes(colour=region, linetype=data_type))+
  geom_ribbon(aes(ymin=ensemble_mean_lower_sc, ymax=ensemble_mean_upper_sc, fill=region), alpha=0.2)+
  facet_wrap("species", ncol=4, scales="free_y")+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE, labels=labs )+
  scale_fill_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE,labels=labs )+
  scale_linetype_manual(values=c("dashed", "solid"))+
  guides(fill = guide_legend(nrow = 2, labels=labs), color=guide_legend(nrow=2, labels=labs, override.aes=list(size=4)), linetype=guide_legend(nrow=2))+
  xlab("Oxygen (kPa) at 12 C")+
  ylab("Effect on Fish Density")

ggsave(
  paste0("output/", output_folder, "/plots_final/species_filtered_onedat.png"),
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

##-----Historical observations and conditional effects
bp_est2 <- readRDS(paste0("output/", output_folder, "/breakpoint_estimates_filtered.rds"))
#Get all insitu data
#insitu <- readRDS("data/processed_data/o2/insitu_combined.rds")
insitu <- as.data.frame(readRDS("data/processed_data/o2/all_o2_dat_filtered.rds"))
#Just columns of interest
dat <- insitu[,c("survey", "depth", "year", "latitude", "longitude", "X", "Y", "temp", "salinity_psu","sigma0", "region", "o2")]
dat$temperature_C <- dat$temp
dat$O2_umolkg <- dat$o2
#Drop if NA in O2
dat <- drop_na(dat, latitude, longitude, O2_umolkg, temperature_C, salinity_psu, year)
#Remove weird depths
dat <- filter(dat, depth>0)

#Remove oxygen outliers
dat <- filter(dat, o2<1500)

#Set minimum sigma
minsigma0 <- 24
dat$sigma0[dat$sigma0 <= minsigma0] <- minsigma0

# remove older (earlier than 2000) data
dat <- dplyr::filter(dat, year >=2000)

#Calculate pO2
#Calculate pO2 from umol kg
dat$po2 <- calc_po2_sat(salinity=dat$salinity_psu, temp=dat$temperature_C, depth=dat$depth, oxygen=dat$O2_umolkg, lat=dat$latitude, long=dat$longitude, umol_m3=T, ml_L=F)

#Log depth
dat$depth_ln <- log(dat$depth)

#Remove ai
dat <- filter(dat, region!="ai")

#Calculate inverse temperature
#Calc invtemp
kelvin = 273.15
boltz = 0.000086173324
tref <- 12
dat$invtemp <- (1 / boltz)  * ( 1 / (dat$temperature_C + 273.15) - 1 / (tref + 273.15))

#Add an ID
dat$pred_id <- 1:nrow(dat)

###Loop through for each species for all data (entire dataset--not limited to range)
#----Conditional Effect: Manually w/ Monte Carlo parameter draws
for(i in 1:nrow(bp_est2)){
  this_species <- bp_est2$species[i]
  this_data <- bp_est2$data[i]
  print(this_species)
  #Calculate metabolic index from correct taxa parameters
  #Pull correct 
  species_tab <- filter(species_table, common_name==this_species)
  #Calculate metabolic index from correct taxa parameters
  #Pull correct Eo parameter
  mi_pars <- read.csv("data/taxa_table.csv")
  mi_pars$Group <- tolower(mi_pars$Group)
  this_taxa <- species_tab$MI_Taxa[species_tab$common_name==this_species]
  mi_pars <- filter(mi_pars,Group==this_taxa)
  
  Eo <- c(mi_pars$Eolow, mi_pars$Eo, mi_pars$Eohigh)
  
  #New dataframe
  preds <- dat
  
  #Calculate metabolic index for species
  preds$mi1 <- calc_mi(po2=preds$po2, inv.temp=preds$invtemp, Eo=Eo[1],fancy=F)
  preds$mi2 <- calc_mi(po2=preds$po2, inv.temp=preds$invtemp, Eo=Eo[2],fancy=F)
  preds$mi3 <- calc_mi(po2=preds$po2, inv.temp=preds$invtemp, Eo=Eo[3],fancy=F)
  
  #Filter depth and latitude for species
  preds <- filter(preds, depth<(species_tab$depth+200))
  northern_limit <- species_tab$northern_limit
  southern_limit <- species_tab$southern_limit
  if(!is.na(northern_limit)){
    preds<- filter(preds, latitude<northern_limit)
  }
  if(!is.na(southern_limit)){
    preds <- filter(preds, latitude>southern_limit)
  }
  
  #Pull the data file for scaling
  this_datframe <- try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_data, "_dat.rds")))
  #Mean and sd for unscaling MIs
  mean_mi1 <- mean(this_datframe$mi1, na.rm=T)
  sd_mi1 <- sd(this_datframe$mi1, na.rm=T)
  mean_mi2 <- mean(this_datframe$mi2, na.rm=T)
  sd_mi2 <- sd(this_datframe$mi2, na.rm=T)
  mean_mi3 <- mean(this_datframe$mi3, na.rm=T)
  sd_mi3 <- sd(this_datframe$mi3, na.rm=T)
  
  #Calculate
  preds$mi1_s <- (preds$mi1-mean_mi1)/sd_mi1
  preds$mi2_s <- (preds$mi2-mean_mi2)/sd_mi2
  preds$mi3_s <- (preds$mi3-mean_mi3)/sd_mi3

  #List of data types to pull
  this_aic <- filter(aic, species==this_species& `data type`==this_data)
  #For each data type, pull all of the models that did fit
  this_dat <- this_data
  print(this_dat)
  #Just the first 5 columns
  this_aic <- this_aic[,1:(ncol(this_aic)-5)]
  #Only the columns that are not NAs (which means they didn't pass sanity checks)
  this_aic <- this_aic[,colSums(is.na(this_aic))<nrow(this_aic)]
    
  #Get list of these columns (these are the model fits to pull for model averaging)
  models <- colnames(this_aic)
  
  #Pull the model output files
  model_fits <- list()
  for(h in 1:length(models)){
    fit <- try(readRDS(file = paste0("output/", output_folder,"/", this_species,"_",this_dat,"_", models[h], ".rds")))
    model_fits[[h]] <- fit
  }
  #Get model weights
  aics <- list()
  for (k in 1:length(models)){
  aic_models <- AIC(model_fits[[k]])
  aics[[k]] <- aic_models
  }
 aics <- unlist(aics)
 delta_aic <- aics - min(aics)
 weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
 set.seed(459384) # for reproducibility and consistency
 for(g in 1:length(model_fits)){
        if(models[g] %in% mi_models2use){
          fit <- model_fits[[g]]
          pars <- as.data.frame(tidy(fit, effects="fixed", conf.int=T))
          slope <- filter(pars, grepl("slope", term))
          thresh <- filter(pars,grepl("breakpt", term))
          #Pull 100 threshold estimates or slope estimates
          bp_pars <- rnorm(mean=thresh$estimate, sd=thresh$std.error, n=100)
          slope_pars <- rnorm(mean=slope$estimate, sd=slope$std.error, n=100)
          #Set slope--slope <- slope$estimate
          #Calculate the breakpoint effect
          if(models[g]=="model13"){
            mi_s <- preds$mi1_s
          }
          if(models[g]=="model14"){
            mi_s <- preds$mi2_s
          }
          if(models[g]=="model15"){
            mi_s <- preds$mi3_s
          }
          for(l in 1:length(bp_pars)){
            pred <- as.data.frame(sapply(mi_s,breakpoint_calc, slope_pars[l],bp_pars[l]))
            colnames(pred) <- "pred"
            #repeat this value the length of nrow(dat)
            max <- rep((bp_pars[l]+1), nrow(preds))
            pred$max_effect <- sapply(max,breakpoint_calc, slope_pars[l],bp_pars[l])
            pred$est_effect_prop <- (exp(pred$max_effect)-exp(pred$pred))/exp(pred$max_effect)
            pred$weight <- weights[g]
            pred$sim <- l
            pred$pred_id <-1:nrow(preds)
            
            #Combine together
            if(l==1){
              pred_model <- pred
            }
            if(l>1){
              pred_model <- bind_rows(pred_model, pred)
            }
          }
          
        } else {
          #For models without pO2' term, conditional effect will be 0--add 100 iterations of this for each prediction point
          pred_model <- as.data.frame(matrix(0, nrow=nrow(preds)*100))
          colnames(pred_model) <- "pred"
          pred_model$est_effect_prop <- 0
          pred_model$weight <- weights[g]
          pred_model$pred_id <- rep(1:nrow(preds), 100)
          pred_model$sim <- rep(1:100, each=nrow(preds))
        }
        if(g==1){
          pred_all <- pred_model
        }
        if(g>1){
          pred_all <- bind_rows(pred_all, pred_model)
        }
 }     
    
  #Calculate weighted average
      ens_preds <- pred_all %>% 
        group_by(pred_id, sim) %>% 
        summarise(weighted_mean=weighted.mean(est_effect_prop,weight, na.rm=T))  %>% 
        ungroup()
      
      #Ensemble w/ SD
      ens_preds2 <- ens_preds %>% 
        group_by(pred_id) %>% 
        summarise(ensemble_mean=mean(weighted_mean, na.rm=T), ensemble_sd=sd(weighted_mean, na.rm=T)) %>%
        ungroup()
      ens_preds2 <- ens_preds2 %>% 
        mutate(ensemble_lower=(ensemble_mean-ensemble_sd), ensemble_upper=(ensemble_mean+ensemble_sd))%>% 
        ungroup()
      
      #Add prediction dataframe back
      preds5 <- left_join(preds, ens_preds2, by="pred_id")
      #Other columns
      preds5$data <- this_data 
      preds5$common_name <- this_species
      preds5$depth_limit <- species_tab$depth
      preds5$min_lat <- species_tab$northern_limit
      preds5$max_lat <- species_tab$southern_limit
      if(!grepl("iphc", this_dat)){
        preds5$data_type <- "bottom trawl only"
      }
      if(grepl("iphc", this_dat)){
        preds5$data_type <- "bottom trawl & IPHC"
      }
      if(i==1){
        all_obs <- preds5}
      if(i>1){
        all_obs <- bind_rows(all_obs, preds5)
      }
}

saveRDS(all_obs, file=paste0("output/", output_folder, "/", "conditional_effects_obs_points_all.rds"))
all_obs <- readRDS(paste0("output/", output_folder, "/", "conditional_effects_obs_points_all"))

##Plot just coastwide species, on one plot
dat2plot <- filter(all_obs, grepl("coastwide", data))
#Depth buffer (200m than typical depth, and 0.5 degree latitude north and south)
#dat2plot <- filter(dat2plot, (depth<depth_limit+200))
#dat2plot <- filter(dat2plot, (latitude>(min_lat-0.5))&(latitude<(max_lat+0.5)))

ggplot(us_coast_proj) + geom_sf() +
  #geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="#0D0887FF", size=0.5, alpha=0.1)+
  geom_point(filter(dat2plot, ensemble_mean==0),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.5, alpha=0.1)+
  geom_point(filter(dat2plot,ensemble_mean>0),mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.9)+
   #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  facet_wrap("common_name", labeller=labeller(common_name=label_wrap_gen(20)))+
  scale_x_continuous(breaks=c(-150,-135,-120), limits=c(min(dat2plot$X)*1000, max(dat2plot$X)*1000))+
  ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
  theme_minimal(base_size=18)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position="top", legend.text = element_text(angle = 45, hjust = 1),legend.justification="center",
        legend.box.spacing = unit(0, "pt"), panel.spacing = unit(1.2, "lines"))+
  #scale_colour_viridis(name="Reduction in biomass")
  scale_colour_viridis(name="Reduction in \nbiomass ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

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

##Plot just species w/ Alaska
temp <- filter(all_obs, region=="goa"&(common_name=="sablefish"|common_name=="rougheye rockfish"|common_name=="walleye pollock"))
dat2plot <- filter(dat2plot, (depth<depth_limit+200))
#temp <- filter(temp, (latitude>(min_lat-0.5))&(latitude<(max_lat+0.5)))
ggplot(us_coast_proj) + geom_sf() +
  #geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="#0D0887FF", size=0.5, alpha=0.1)+
  geom_point(filter(temp, ensemble_mean==0),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.5, alpha=0.1)+
  geom_point(filter(temp,ensemble_mean>0),mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.9)+
  #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  scale_x_continuous(breaks=c(-160,-150,-140), limits=c(min(temp$X)*1000, max(temp$X)*1000))+
  facet_wrap("common_name", ncol=1)+
  scale_y_continuous(breaks=c(46,48,50,52), limits=c(min(temp$Y)*1000, max(temp$Y)*1000))+
  #facet_wrap("year", ncol=5)+
  theme_minimal(base_size=18)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position="top", legend.text = element_text(angle = 45, hjust = 1), legend.justification="center",
        legend.box.spacing = unit(0, "pt"))+
  #scale_colour_viridis(name="Reduction in biomass")
  scale_colour_viridis(name="Reduction in \nbiomass ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")
#  limits=c(0,1), breaks=c(0,0.25, 0.5, 0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish)

ggsave(
  paste0("output/", output_folder, "/plots/AK_species.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 11,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##Plot annually
toplot <- unique(dat2plot$common_name)
for(i in 1:length(toplot)) {
  this_species <- toplot[i]
  temp <- filter(dat2plot, common_name==this_species&year>2009)
  ggplot(us_coast_proj) + geom_sf() +
    #geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="#0D0887FF", size=0.5, alpha=0.1)+
    geom_point(filter(temp, ensemble_mean==0),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.5, alpha=0.1)+
    geom_point(filter(temp,ensemble_mean>0),mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.9)+
    #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    facet_wrap("year", ncol=4)+
    scale_x_continuous(breaks=c(-150,-130), limits=c(min(temp$X)*1000, max(temp$X)*1000))+
    ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
    #facet_wrap("year", ncol=5)+
    theme_minimal(base_size=18)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(legend.position=c(0.92,0.1), panel.spacing = unit(1, "lines"))+
    ggtitle(this_species)+
    #scale_colour_viridis(name="Reduction in biomass")
    scale_colour_viridis(name="Reduction in \nbiomass ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")
  #  limits=c(0,1), breaks=c(0,0.25, 0.5, 0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish)
  
  ggsave(
    paste0("output/", output_folder, "/", "plots/points_full_annual/map_",this_species,".png"),
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

#Plot annually with line where observation below
dat2use <- filter(all_obs, common_name=="walleye pollock"|common_name=="sablefish")
dat2use$region <- factor(dat2use$region, levels=c("ebs", "goa", "bc", "cc"))
labs <- c("Eastern Bering Sea", "Gulf of Alaska", "British Columbia", "California Current")
names(labs) <- c("ebs", "goa", "bc", "cc")
bp_est2$common_name <- bp_est2$species
ggplot(dat2use, aes(x=po2))+
  geom_density(aes(colour=region))+
  facet_wrap("common_name", ncol=1)+
  theme(legend.position="top")+
  xlab("Dissolved Oxygen (pO2)")+
guides(color = guide_legend(nrow = 2))+
  geom_vline(filter(bp_est2, species %in% unique(dat2use$common_name)&region=="goa"), mapping=aes(xintercept=bp_ensemble_mean), linetype="dashed")+
  scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677"), drop=FALSE, labels=labs )

ggsave(
  paste0("output/", output_folder, "/", "plots/density_obs_example.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5,
  height = 6,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)


##proportion of observations in each region with decline
dat2plot <- drop_na(dat2plot, ensemble_mean, ensemble_lower, ensemble_upper)
effects_full_annual_summary <- dat2plot %>%
  group_by(year, region, common_name) %>%
  summarise(prop_below = sum(ensemble_mean>.1)/n(), N=n()) %>%
  ungroup()
effects_full_annual_summary2 <- dat2plot %>%
  group_by(year, region, common_name) %>%
  summarise(prop_below = sum(ensemble_lower>.1)/n()) %>%
  ungroup()
effects_full_annual_summary3 <- dat2plot %>%
  group_by(year, region, common_name) %>%
  summarise(prop_below = sum(ensemble_upper>.1)/n()) %>%
  ungroup()

effects_full_annual_summary$se1 <- effects_full_annual_summary2$prop_below
effects_full_annual_summary$se2 <- effects_full_annual_summary3$prop_below

#Factor region
effects_full_annual_summary$region <- factor(effects_full_annual_summary$region, levels=c("ebs", "goa", "bc", "cc"))
labs <- c("Eastern Bering Sea", "Gulf of Alaska", "British Columbia", "California Current")
names(labs) <- c("ebs", "goa", "bc", "cc")

#At least 50 observations in region
effects_full_annual_summary <- filter(effects_full_annual_summary, N>50)
#plot
ggplot(effects_full_annual_summary, aes(x=year, y=prop_below))+
  geom_line(aes(colour=region))+
  #xlim(2010,2024)+
  geom_ribbon(aes(ymin=se2, ymax=se1, fill=region), alpha=0.2)+
  facet_wrap("common_name")+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=20))+
 scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677"), drop=FALSE, labels=labs )+
 scale_fill_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  xlab("Year") +
  ylab("Proportion of observations w/ >10% biomass reduction")

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_below_obs_annual.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height = 8.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

####Proportion below threshold by depth
dat2plot <- drop_na(dat2plot, ensemble_mean, ensemble_lower, ensemble_upper)
dat2plot$depth_bin <- cut(dat2plot$depth, labels=FALSE, breaks=seq(0, 1500, by=50), include.lowest = TRUE)
dat2plot$depth_bin <- dat2plot$depth_bin*50
effects_full_annual_summary <- dat2plot %>%
  group_by(depth_bin, region, common_name) %>%
  summarise(prop_below = sum(ensemble_mean>0.1)/n(), N=n()) %>%
  ungroup()
effects_full_annual_summary2 <- dat2plot %>%
  group_by(depth_bin, region, common_name) %>%
  summarise(prop_below = sum(ensemble_lower>0.1)/n()) %>%
  ungroup()
effects_full_annual_summary3 <- dat2plot %>%
  group_by(depth_bin, region, common_name) %>%
  summarise(prop_below = sum(ensemble_upper>0.1)/n()) %>%
  ungroup()
effects_full_annual_summary$se1 <- effects_full_annual_summary2$prop_below
effects_full_annual_summary$se2 <- effects_full_annual_summary3$prop_below


#Factor region
effects_full_annual_summary$region <- factor(effects_full_annual_summary$region, levels=c("ebs", "goa", "bc", "cc"))
labs <- c("Eastern Bering Sea", "Gulf of Alaska", "British Columbia", "California Current")
names(labs) <- c("ebs", "goa", "bc", "cc")
effects_full_annual_summary <- filter(effects_full_annual_summary, N>50)
ggplot(effects_full_annual_summary, aes(x=depth_bin, y=prop_below))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  facet_wrap("common_name", scales="free_x")+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=20))+
  scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  scale_fill_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  xlab("Depth (m)") +
  ylab("Proportion of observations w/ >10% biomass reduction")+
  geom_vline(dat2plot, mapping=aes(xintercept=depth_limit), linetype="dashed", colour="black")

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_below_by_depth.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height = 8.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

###Predict to grid
bc_grid <- readRDS("data/processed_data/o2/bc_predictions_grid.rds")
bc_grid$latitude <- bc_grid$latitute
cc_grid <- readRDS("data/processed_data/o2/cc_predictions_grid.rds")
cc_grid$latitude <- cc_grid$latitute

#Calculate inverse temp
#Calc invtemp
kelvin = 273.15
boltz = 0.000086173324
tref <- 12
cc_grid$invtemp <- (1 / boltz)  * ( 1 / (cc_grid$temperature_C + 273.15) - 1 / (tref + 273.15))
bc_grid$invtemp <- (1 / boltz)  * ( 1 / (bc_grid$temperature_C + 273.15) - 1 / (tref + 273.15))
cc_grid2 <- filter(cc_grid, survey=="NWFSC.Combo")

#Plot oxygen data


#Pull correct Eo parameter
mi_pars <- read.csv("data/taxa_table.csv")
mi_pars$Group <- tolower(mi_pars$Group)

#combine grid
grid <- bind_rows(bc_grid, cc_grid)
grid$pred_id <- 1:nrow(grid)

options(dplyr.summarise.inform = FALSE)

bp_est2 <- readRDS(paste0("output/", output_folder, "/breakpoint_estimates_filtered.rds"))

#Function to calculate for each grid cell
calc_grid <- function(mi1_s, mi2_s, mi3_s, pred_id, pars_list, models, model_fits, mi_models2use, weights,n) {
  # Define locally within the worker
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
  
  for(g in 1:length(model_fits)){
    if(models[g] %in% mi_models2use){
      #Calculate the breakpoint effect
      if(models[g]=="model13"){
        mi_s <- mi1_s
      }
      if(models[g]=="model14"){
        mi_s <- mi2_s
      }
      if(models[g]=="model15"){
        mi_s <- mi3_s
      }
      #Pull bp_pars and slope_pars from list
      bp_pars <- pars_list[[g]]$bp_pars
      slope_pars <- pars_list[[g]]$slope_pars
      
      for(l in 1:length(bp_pars)){
        pred <- sapply(mi_s,breakpoint_calc, slope_pars[l],bp_pars[l])
        #repeat this value the length of nrow(dat)
        max <- bp_pars[l]+1
        max_effect <- sapply(max,breakpoint_calc, slope_pars[l],bp_pars[l])
        est_effect_prop <- (exp(max_effect)-exp(pred))/exp(max_effect)
        weight <- weights[g]
        sim <- l
        pred_id <- pred_id
        
        #Combine together
        if(l==1){
          est_effect_prop2 <- est_effect_prop
          weight2 <- weight
          sim2 <- sim
          pred_id2 <- pred_id
        }
        
        if(l>1){
          est_effect_prop2 <- c(est_effect_prop2, est_effect_prop)
          weight2 <- c(weight2, weight)
          sim2 <- c(sim2, sim)
          pred_id2 <- c(pred_id2, pred_id)
        }
      }
      
    } else {
      #For models without pO2' term, conditional effect will be 0--add 100 iterations of this for each prediction point
      est_effect_prop2 <- rep(0, 100)
      weight2 <- rep(weights[g], 100)
      pred_id2 <- rep(pred_id, 100)
      sim2 <- rep(1:100, each=1)
    }
    if(g==1){
    df <- data.frame(est_effect_prop=est_effect_prop2, weight=weight2, sim=sim2, pred_id=pred_id2)
    }
    if(g>1){
    df <- bind_rows(df, data.frame(est_effect_prop=est_effect_prop2, weight=weight2, sim=sim2, pred_id=pred_id2))
    }
  }     
  
  #Calculate weighted average
  ens_preds <- df %>% 
    group_by(pred_id, sim) %>% 
    summarise(weighted_mean=weighted.mean(est_effect_prop,weight, na.rm=T))  %>% 
    ungroup()
  
  #Ensemble w/ SD
ens_preds2 <- ens_preds %>% 
    group_by(pred_id) %>% 
    summarise(ensemble_mean=mean(weighted_mean, na.rm=T), ensemble_sd=sd(weighted_mean, na.rm=T)) %>%
    ungroup()
  ens_preds2 <- ens_preds2 %>% 
    mutate(ensemble_lower=(ensemble_mean-ensemble_sd), ensemble_upper=(ensemble_mean+ensemble_sd))%>% 
    ungroup()
  return(ens_preds2)
}

process_species <- function(species2do){
  this_species <- species2do
  bp_est3 <- filter(bp_est2, species==this_species)
  this_data <- bp_est3$data
  if(this_species=="sablefish"){
    this_data <- "cc"
  }
  print(this_species)
  #Calculate metabolic index from correct taxa parameters
  #Pull correct 
  species_tab <- filter(species_table, common_name==this_species)
  #Calculate metabolic index from correct taxa parameters
  this_taxa <- species_tab$MI_Taxa[species_tab$common_name==this_species]
  mi_pars2 <- filter(mi_pars,Group==this_taxa)
  
  Eo <- c(mi_pars2$Eolow, mi_pars2$Eo, mi_pars2$Eohigh)
  
  #New dataframe
  preds <- grid
  
  #Filter by depth and latitude
  #Depth buffer (200m than typical depth, and 0.5 degree latitude north and south)
  #Filter depth and latitude for species
  preds <- filter(preds, depth_m<(species_tab$depth+200))
  northern_limit <- species_tab$northern_limit
  southern_limit <- species_tab$southern_limit
  if(!is.na(northern_limit)){
    preds<- filter(preds, latitude<northern_limit)
  }
  if(!is.na(southern_limit)){
    preds <- filter(preds, latitude>southern_limit)
  }
  
  #Calculate metabolic index for species
  preds$mi1 <- calc_mi(po2=preds$po2, inv.temp=preds$invtemp, Eo=Eo[1],fancy=F)
  preds$mi2 <- calc_mi(po2=preds$po2, inv.temp=preds$invtemp, Eo=Eo[2],fancy=F)
  preds$mi3 <- calc_mi(po2=preds$po2, inv.temp=preds$invtemp, Eo=Eo[3],fancy=F)
  
  #Pull the data file for scaling
  this_datframe <- try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_data, "_dat.rds")))
  #Mean and sd for unscaling MIs
  mean_mi1 <- mean(this_datframe$mi1, na.rm=T)
  sd_mi1 <- sd(this_datframe$mi1, na.rm=T)
  mean_mi2 <- mean(this_datframe$mi2, na.rm=T)
  sd_mi2 <- sd(this_datframe$mi2, na.rm=T)
  mean_mi3 <- mean(this_datframe$mi3, na.rm=T)
  sd_mi3 <- sd(this_datframe$mi3, na.rm=T)
  
  #Calculate
  preds$mi1_s <- (preds$mi1-mean_mi1)/sd_mi1
  preds$mi2_s <- (preds$mi2-mean_mi2)/sd_mi2
  preds$mi3_s <- (preds$mi3-mean_mi3)/sd_mi3
  mi1_s <- preds$mi1_s
  mi2_s <- preds$mi2_s
  mi3_s <- preds$mi3_s
  pred_id <- preds$pred_id

  #Which AIC/data to pull
  this_aic <- filter(aic, species==this_species& `data type`==this_data)
  this_dat <- this_data
  print(this_dat)
  #Just the first 5 columns
  this_aic <- this_aic[,1:(ncol(this_aic)-5)]
  #Only the columns that are not NAs (which means they didn't pass sanity checks)
  this_aic <- this_aic[,colSums(is.na(this_aic))<nrow(this_aic)]
  
  #Get list of these columns (these are the model fits to pull for model averaging)
  models <- colnames(this_aic)
  
  #Pull the model output files
  model_fits <- list()
  for(h in 1:length(models)){
    fit <- try(readRDS(file = paste0("output/", output_folder,"/", this_species,"_",this_dat,"_", models[h], ".rds")))
    model_fits[[h]] <- fit
  }
  #Get model weights
  aics <- list()
  for (k in 1:length(models)){
    aic_models <- AIC(model_fits[[k]])
    aics[[k]] <- aic_models
  }
  aics <- unlist(aics)
  delta_aic <- aics - min(aics)
  weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
  set.seed(459384) # for reproducibility and consistency
  #Get lists of parameters to use
  pars_list <- list()
  for(g in 1:length(model_fits)){
    if(models[g] %in% mi_models2use){
      fit <- model_fits[[g]]
      pars <- as.data.frame(tidy(fit, effects="fixed", conf.int=T))
      slope <- filter(pars, grepl("slope", term))
      thresh <- filter(pars,grepl("breakpt", term))
      #Pull 100 threshold estimates or slope estimates
      bp_pars <- rnorm(mean=thresh$estimate, sd=thresh$std.error, n=100)
      slope_pars <- rnorm(mean=slope$estimate, sd=slope$std.error, n=100)
      par <- list(bp_pars=bp_pars, slope_pars=slope_pars)
      pars_list[[g]] <- par
    } else {
      pars_list[[g]] <- NA
    }
  }

  #Calculate for each grid cell, parallelized
  #use mcmapply to parallelize
  n <- nrow(preds)
  results <- mcmapply(
    FUN = function(z) {
      calc_grid(
        mi1_s = mi1_s[z],
        mi2_s = mi2_s[z],
        mi3_s = mi3_s[z],
        pred_id=pred_id[z],
        pars_list = pars_list,
        models = models,
        model_fits = model_fits,
        mi_models2use = mi_models2use,
        weights = weights,
        n=n
      )
    },
    z = 1:n,
    mc.cores = 8,
    SIMPLIFY = FALSE
  )
  #convert results to single dataframe
  preds5 <- as.data.frame(bind_rows(results))

  #Add prediction dataframe back
  preds5 <- left_join(preds, preds5, by="pred_id")
  #Other columns
  preds5$data <- this_data 
  preds5$common_name <- this_species
  preds5$depth_limit <- species_tab$depth
  preds5$min_lat <- species_tab$southern_limit
  preds5$max_lat <- species_tab$northern_limit
  if(!grepl("iphc", this_dat)){
    preds5$data_type <- "bottom trawl only"
  }
  if(grepl("iphc", this_dat)){
    preds5$data_type <- "bottom trawl & IPHC"
  }
  saveRDS(preds5, file = paste0("output/", output_folder, "/", this_species, "_", "grid.rds"))
  return(preds5)
}

#Apply for all the species
species_list <- filter(bp_est2, grepl("coastwide|bc|cc", data))
species_list <- unique(species_list$species)
#apply to all species (if need to run)
grids <- list()
for(i in 1:length(species_list)){
grids_x <- process_species(species_list[i])
grids[[i]] <- grids_x
}

#Or load in
grids <- list()
for(i in 1:length(species_list)){
this_species <- species_list[i]
grid_x <- readRDS(paste0("output/", output_folder, "/", this_species, "_", "grid.rds"))
grids[[i]] <- grid_x
}

for(i in 1:length(grids)){
  grid_x <- grids[[i]]
  grid_x <- filter(grid_x, region=="bc"|(region=="cc"&survey=="NWFSC.Combo"))
  grids[[i]] <- grid_x
 saveRDS(grid_x,paste0("output/", output_folder, "/", this_species, "_", "grid.rds"))
  
}

#Filter to one year
for(i in 1:length(grids)){
grids2use <- grids[[i]]
grids2use <- filter(grids2use, year==2021)
#grids2use <- filter(grids2use,depth_m<(depth_limit+200))
#grids2use <- filter(grids2use, (latitude>(min_lat-0.5))&(latitude<(max_lat+0.5)))

if(i==1){
  grids2 <- grids2use
}
if(i>1){
  grids2 <- bind_rows(grids2use, grids2)
}
}

theme_set(theme_bw(base_size = 25))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

#Plot one year, each species, full region
ylims <- c(min(grids2$Y)*1000, max(grids2$Y)*1000)
xlims <- c(min(grids2$X)*1000, max(grids2$X)*1000)

ggplot(us_coast_proj) + geom_sf() +
  geom_point(filter(grids2, ensemble_mean==0),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.5, alpha=0.1)+
  geom_point(filter(grids2,ensemble_mean>0),mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.9)+
  #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  scale_x_continuous(breaks=c(-130,-120), limits=c(xlims))+
  ylim(ylims)+
  facet_wrap("common_name", ncol=5, labeller=labeller(common_name=label_wrap_gen(10)))+
  theme_minimal(base_size=16)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position=c(0.93,0.1), panel.spacing = unit(1.5, "lines"), strip.text=element_text(size=13))+
  scale_colour_viridis(name="Reduction in \nbiomass ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

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

#Plot one year, each species, WA state
###WASHINGTON COAST
#Plot just one year of data, all species on same plot, zoom to WA
dat2plot_wa <- filter(grids2, latitude>46&latitude<48.5)
ggplot(us_coast_proj) + geom_sf() +
  geom_point(dat2plot_wa,mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  scale_x_continuous(breaks=c(-125,-120), limits=c(min(dat2plot_wa$X)*1000, max(dat2plot_wa$X)*1000))+
  ylim(min(dat2plot_wa$Y)*1000, max(dat2plot_wa$Y)*1000)+
  facet_wrap("common_name", ncol=5, labeller=labeller(common_name=label_wrap_gen(10)))+
  theme_minimal(base_size=18)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position=c(0.94,0.1), panel.spacing=unit(1,"lines"))+
  scale_colour_viridis(name="Reduction in \nbiomass ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

ggsave(
  paste0("output/", output_folder, "/", "plots/grid_combined_2021_WA.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 9,
  height = 11,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##Compare two years, in WA state
#Filter to one year
for(i in 1:length(grids)){
  grids2use <- grids[[i]]
  grids2use <- filter(grids2use, year==2010|year==2021)
#  grids2use <- filter(grids2use,depth_m<(depth_limit+200))
#  grids2use <- filter(grids2use, (latitude>(min_lat-0.5))&(latitude<(max_lat+0.5)))
  grids2use <- filter(grids2use, latitude>46&latitude<48.5)
  if(i==1){
    grids3 <- grids2use
  }
  if(i>1){
    grids3 <- bind_rows(grids2use, grids3)
  }
}

temp <- filter(grids3, common_name=="canary rockfish"|common_name=="shortspine thornyhead"|common_name=="redbanded rockfish"|common_name=="pacific hake")

ggplot(us_coast_proj) + geom_sf() +
  geom_point(temp,mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.5)+
  #geom_point(filter(dat2plot_wa, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot_wa, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  scale_x_continuous(breaks=c(-125,-120), limits=c(min(temp$X)*1000, max(temp$X)*1000))+
  ylim(min(temp$Y)*1000, max(temp$Y)*1000)+
  facet_nested(~common_name~year)+
  theme_minimal(base_size=18)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position="top", legend.text = element_text(angle = 45, hjust = 1), legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  scale_colour_viridis(name="Reduction in \nbiomass ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

ggsave(
  paste0("output/", output_folder, "/", "plots/grid_year_comp_WA.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 11,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

#Plot all years, each species separately
for(i in 1:length(grids)) {
  to_do <- grids[[i]]
  this_species <- unique(to_do$common_name)
  print(this_species)
  this_dat <- unique(to_do$data)
 # to_do <- filter(to_do, depth_m<(depth_limit+200))
# to_do<- filter(to_do, (latitude>(min_lat-0.5))&(latitude<(max_lat+0.5)))
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(to_do,mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    scale_x_continuous(breaks=c(-125,-120), limits=c(min(to_do$X)*1000, max(to_do$X)*1000))+
    ylim(min(to_do$Y)*1000, max(to_do$Y)*1000)+
    facet_wrap("year", ncol=5)+
    theme_minimal(base_size=16)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(legend.position=c(0.8,0.1))+
    #ggtitle(paste(this_species, this_dat, sep=" "))
    scale_colour_viridis(name="Reduction in \nbiomass ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

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

#Plot all years, each species separately, WA only
for(i in 1:length(grids)) {
  to_do <- grids[[i]]
  this_species <- unique(to_do$common_name)
  print(this_species)
  this_dat <- unique(to_do$data)
  #to_do <- filter(to_do, depth_m<(depth_limit+200))
  #to_do<- filter(to_do, (latitude>(min_lat-0.5))&(latitude<(max_lat+0.5)))
  to_do <- filter(to_do, latitude>46&latitude<48.5)
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(to_do,mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    scale_x_continuous(breaks=c(-125,-120), limits=c(min(to_do$X)*1000, max(to_do$X)*1000))+
    ylim(min(to_do$Y)*1000, max(to_do$Y)*1000)+
    facet_wrap("year", ncol=5)+
    theme_minimal(base_size=16)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(legend.position=c(0.8,0.1))+
    #ggtitle(paste(this_species, this_dat, sep=" "))+
    scale_colour_viridis(name="Reduction in \nbiomass ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")
  
  ggsave(
    paste0("output/", output_folder, "/", "plots/grid_annual_WA/map_",this_species,"_", this_dat,".png"),
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

#Annual proportion of habitat by year
for(i in 1:length(grids)){
test <- grids[[i]]
#test <- filter(test, depth_m<(depth_limit+200))
#test <- filter(test, (latitude>(min_lat-0.5))&(latitude<max_lat+0.5))

area_sum2 <- filter(test, year==2021)
areas <- area_sum2 %>%
  group_by(region) %>%
  summarise(area_sum=sum(area))  %>%
  ungroup()

test <- left_join(test, areas, by=c("region"))

test2 <- filter(test, ensemble_mean>0.1)
test3 <- filter(test, ensemble_lower>0.1)
test4 <- filter(test, ensemble_upper>0.1)

effects_full_annual_summary <- test2 %>%
  group_by(year, region, common_name) %>%
  summarise(area_sum=sum(area)/area_sum) %>%
  distinct()%>%
  ungroup()

effects_full_annual_summary2 <- test3 %>%
  group_by(year, region, common_name) %>%
  summarise(area_sum=sum(area)/area_sum) %>%
  distinct()%>%
  ungroup()
effects_full_annual_summary3 <- test4 %>%
  group_by(year, region, common_name) %>%
  summarise(area_sum=sum(area)/area_sum) %>%
  distinct()%>%
  ungroup()

effects_full_summary <- full_join(effects_full_annual_summary, effects_full_annual_summary2, by=c("year", "region", "common_name"))
effects_full_summary$se1 <- effects_full_summary$area_sum.y
effects_full_summary$area_sum_med <- effects_full_summary$area_sum.x
effects_full_summary$area_sum.y <- NULL
effects_full_summary$area_sum.x <- NULL
effects_full_summary <- full_join(effects_full_summary, effects_full_annual_summary3, by=c("year", "region", "common_name"))
effects_full_summary$se2 <- effects_full_summary$area_sum
effects_full_summary$area_sum <- NULL
effects_full_summary$species <- unique(test$common_name)

if(i==1){
effects <- effects_full_summary
}
if(i>1){
effects <- bind_rows(effects, effects_full_summary)
}
}

#Regions
effects$region <- factor(effects$region, levels=c("bc", "cc"))
labs <- c("British Columbia", "California Current")
names(labs) <- c("bc", "cc")

#plot
ggplot(effects, aes(x=year, y=area_sum_med))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  facet_wrap("common_name", scales="free_y",ncol=4,labeller=labeller(common_name=label_wrap_gen(15)))+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank(), text=element_text(size=20),
  legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  scale_fill_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  scale_colour_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  xlab("Year") +
  ylab("Proportion of area w/ >10% biomass reduction")

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_area_below_grid.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 13,
  height = 13,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

#####Annual proportion of habitat by depth
for(i in 1:length(grids)){
  test <- grids[[i]]
 #test <- filter(test, depth_m<(depth_limit+200))
  #test <- filter(test, (latitude>(min_lat-0.5))&(latitude<max_lat+0.5))
test$depth_bin <- cut(test$depth_m, labels=FALSE, breaks=seq(0, max(test$depth_m), by=25), include.lowest = TRUE)
test$depth_bin <- test$depth_bin*25

test <- filter(test, year==2021)

areas <- test %>%
  group_by(depth_bin, region) %>%
  summarise(total_area=sum(area))  %>%
ungroup()

test <- left_join(test, areas, by=c("region", "depth_bin"))

test2 <- filter(test, ensemble_mean>0.1)
test3 <- filter(test, ensemble_lower>0.1)
test4 <- filter(test, ensemble_upper>0.1)

effects_full_annual_summary <- test2 %>%
  group_by(depth_bin, region, total_area) %>%
  summarize(area_sum = sum(area)) %>%
  mutate(prop=area_sum/total_area) %>%
  ungroup()
effects_full_annual_summary2 <- test3 %>%
  group_by(depth_bin, region, total_area) %>%
  summarise(area_sum = sum(area)) %>%
  mutate(prop=area_sum/total_area) %>%
  ungroup()
effects_full_annual_summary3 <- test4 %>%
  group_by(depth_bin, region, total_area) %>%
  summarise(area_sum = sum(area)) %>%
  mutate(prop=area_sum/total_area) %>%
  ungroup()

effects_full_summary <- full_join(effects_full_annual_summary, effects_full_annual_summary2, by=c("depth_bin", "region"))
effects_full_summary$se1 <- effects_full_summary$prop.y
effects_full_summary$area_sum_med <- effects_full_summary$prop.x
effects_full_summary$prop.y <- NULL
effects_full_summary$prop.x <- NULL
effects_full_summary <- full_join(effects_full_summary, effects_full_annual_summary3, by=c("depth_bin", "region"))
effects_full_summary$se2 <- effects_full_summary$prop
effects_full_summary$prop <- NULL
effects_full_summary$species <- unique(test$common_name)
effects_full_summary$depth_limit <- unique(test$depth_limit)


if(i==1){
  effects2 <- effects_full_summary
}
if(i>1){
  effects2 <- bind_rows(effects2, effects_full_summary)
}
}

#plot
ggplot(effects2, aes(x=depth_bin, y=area_sum_med))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  facet_wrap("species",scales="free_x", ncol=4, labeller=labeller(species=label_wrap_gen(15)))+
  theme(legend.position="top")+
ggh4x::facetted_pos_scales(x=list(scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)), scale_x_continuous(limits=c(0,800), breaks=c(100,300,500, 700)), NULL, scale_x_continuous(breaks=c(100,300,500, 700), limits=c(0,600)), 
                             scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)), scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)),scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)), scale_x_continuous(breaks=c(100,300,500, 700), limits=c(0,700)), 
                             scale_x_continuous(breaks=c(100,300,500), limits=c(0,600)), scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)), scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)), scale_x_continuous(breaks=c(100,300,500), limits=c(0,600)), 
                             scale_x_continuous(breaks=c(100,300,500, 700), limits=c(0,700)), scale_x_continuous(breaks=c(100,500,900), limits=c(0,1200)), scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)), scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)), 
                             scale_x_continuous(breaks=c(100,300,500), limits=c(0,600)), scale_x_continuous(breaks=c(100,300,500), limits=c(0,600))))+
  theme(legend.title=element_blank(), strip.background = element_blank(), text=element_text(size=16),
  legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  scale_fill_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  scale_colour_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  xlab("Depth (m)") +
  ylab("Proportion of area w/ >10% biomass reduction")+
  geom_vline(aes(xintercept=depth_limit), linetype="dashed")

ggsave(
  paste0("output/", output_folder, "/", "plots/grid_area_depth_bins.png"),
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

##Annual with states
#Add states
for(i in 1:length(grids)){
  test <- grids[[i]]
#  test <- filter(test, depth_m<(depth_limit+200))
 # test <- filter(test, (latitude>(min_lat-0.5))&(latitude<max_lat+0.5))
  #add state
  test$state <- case_when(test$latitude<42~"California",
                          test$latitude >42 & test$latitude<46~"Oregon",
                          test$latitude>46& test$region=="cc"~"Washington",
                          test$latitude>46& test$region=="bc"~"British Columbia")
  
  
  area_sum2 <- filter(test, year==2021)
  areas <- area_sum2 %>%
    group_by(state) %>%
    summarise(total_area=sum(area))  %>%
    ungroup()
  
  test <- left_join(test, areas, by=c("state"))
  
  test2 <- filter(test, ensemble_mean>0.1)
  test3 <- filter(test, ensemble_lower>0.1)
  test4 <- filter(test, ensemble_upper>0.1)
  
  effects_full_annual_summary <- test2 %>%
    group_by(year, state, common_name) %>%
    summarise(area_sum=sum(area)/total_area) %>%
    distinct()%>%
    ungroup()
  
  effects_full_annual_summary2 <- test3 %>%
    group_by(year, state, common_name) %>%
    summarise(area_sum=sum(area)/total_area) %>%
    distinct()%>%
    ungroup()
  effects_full_annual_summary3 <- test4 %>%
    group_by(year, state, common_name) %>%
    summarise(area_sum=sum(area)/total_area) %>%
    distinct()%>%
    ungroup()
  
  effects_full_summary <- full_join(effects_full_annual_summary, effects_full_annual_summary2, by=c("year", "state", "common_name"))
  effects_full_summary$se1 <- effects_full_summary$area_sum.y
  effects_full_summary$area_sum_med <- effects_full_summary$area_sum.x
  effects_full_summary$area_sum.y <- NULL
  effects_full_summary$area_sum.x <- NULL
  effects_full_summary <- full_join(effects_full_summary, effects_full_annual_summary3, by=c("year", "state", "common_name"))
  effects_full_summary$se2 <- effects_full_summary$area_sum
  effects_full_summary$area_sum <- NULL
  effects_full_summary$species <- unique(test$common_name)
  
  if(i==1){
    effects <- effects_full_summary
  }
  if(i>1){
    effects <- bind_rows(effects, effects_full_summary)
  }
}

#plot
effects$state <- factor(effects$state, levels=c("British Columbia", "Washington", "Oregon", "California"))

ggplot(effects, aes(x=year, y=area_sum_med))+
  geom_line(aes(colour=state))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=state), alpha=0.2)+
  facet_wrap("common_name", scales="free_y",ncol=4,labeller=labeller(common_name=label_wrap_gen(18)))+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=16),   legend.justification="center",
        legend.box.spacing = unit(0, "pt"))+
  scale_fill_manual(values=c("red", "darkgreen", "royalblue","gold"), drop=FALSE)+
  scale_colour_manual(values=c("red", "darkgreen", "royalblue","gold"), drop=FALSE)+
  scale_x_continuous(breaks=c(2008,2016,2024))+
  xlab("Year") +
  ylab("Proportion of area w/ >10% biomass reduction")

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_area_below_grid_states.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =9,
  height = 11,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

#Persistence mapping--how many years was a grid point more than a 10%, or 50%, decline--or average?
for(i in 1:length(grids)){
test <- grids[[i]]
test <- unique(test)
#test <- filter(test, depth_m<(depth_limit+200))
#test <- filter(test, (latitude>(min_lat-0.5))&(latitude<max_lat+0.5))
test <- filter(test, ensemble_mean>0.1)

count <- test%>%
  group_by(X, Y, latitude, longitude, region) %>%
  summarise(count=n())

count <- test %>%
  group_by(X, Y, latitude, longitude) %>%
summarise(counts=n()) %>%
distinct()%>%
ungroup()

count$species <- unique(test$common_name)

if(i==1){
  count_total <- count
}
if(i>1){
  count_total <- bind_rows(count_total, count)
}
}

xlims <- c(min(count_total$X)*1000, max(count_total$X)*1000)
ylims <- c(min(count_total$Y)*1000, max(count_total$Y)*1000)
ggplot(us_coast_proj) + geom_sf() +
  geom_point(count_total,mapping=aes(x=X*1000, y=Y*1000, colour=counts), size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  scale_x_continuous(breaks=c(-130,-120), limits=c(xlims))+
  ylim(ylims)+
  facet_wrap("species", ncol=5, labeller=labeller(species=label_wrap_gen(10)))+
  theme_minimal(base_size=18)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position=c(0.92,0.1), panel.spacing = unit(1, "lines"))+
scale_colour_viridis(name="Reduction in \nbiomass ",option="plasma")

ggsave(
  paste0("output/", output_folder, "/", "plots/grid_persistence.png"),
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

#Persistence mapping--how many years was a grid point more than a 50%, or 50%, decline--or average?
for(i in 1:length(grids)){
  test <- grids[[i]]
  test <- unique(test)
  test <- filter(test, depth_m<(depth_limit+200))
  test <- filter(test, (latitude>(min_lat-0.5))&(latitude<max_lat+0.5))
  test <- filter(test, ensemble_mean>0.5)
  
  count <- test%>%
    group_by(X, Y, latitude, longitude, region) %>%
    summarise(count=n())
  
  count <- test %>%
    group_by(X, Y, latitude, longitude) %>%
    summarise(counts=n()) %>%
    distinct()%>%
    ungroup()
  
  count$species <- unique(test$common_name)
  
  if(i==1){
    count_total <- count
  }
  if(i>1){
    count_total <- bind_rows(count_total, count)
  }
}

#Add all grid points
grid1 <- filter(grids[[1]], year==2021)

count_total2 <- filter(count_total, latitude>46&latitude<48.5)
xlims <- c(min(count_total2$X)*1000, max(count_total2$X)*1000)
ylims <- c(min(count_total2$Y)*1000, max(count_total2$Y)*1000)
ggplot(us_coast_proj) + geom_sf() +
  geom_point(grid1,mapping=aes(x=X*1000, y=Y*1000), color="grey", size=0.9)+
  geom_point(count_total2,mapping=aes(x=X*1000, y=Y*1000, colour=counts), size=0.9)+
  #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  scale_x_continuous(breaks=c(-130,-120), limits=c(xlims))+
  ylim(ylims)+
  facet_wrap("species", ncol=5, labeller=labeller(species=label_wrap_gen(10)))+
  theme_minimal(base_size=18)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position="top")+
  scale_colour_viridis(name="Number of years \n w/ >50% reduction",option="plasma")

ggsave(
  paste0("output/", output_folder, "/", "plots/grid_persistence_0.5_wa.png"),
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

xlims <- c(min(count_total$X)*1000, max(count_total$X)*1000)
ylims <- c(min(count_total$Y)*1000, max(count_total$Y)*1000)
ggplot(us_coast_proj) + geom_sf() +
  geom_point(grid1,mapping=aes(x=X*1000, y=Y*1000), color="grey", size=0.9)+
  geom_point(count_total,mapping=aes(x=X*1000, y=Y*1000, colour=counts), size=0.9)+
  #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  scale_x_continuous(breaks=c(-130,-120), limits=c(xlims))+
  ylim(ylims)+
  facet_wrap("species", ncol=5, labeller=labeller(species=label_wrap_gen(10)))+
  theme_minimal(base_size=18)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position="top")+
  scale_colour_viridis(name="Number of years \n w/ >50% reduction",option="plasma")

ggsave(
  paste0("output/", output_folder, "/", "plots/grid_persistence_0.5.png"),
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
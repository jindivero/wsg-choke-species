library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(openxlsx2)
library(ggpubr)

setwd("~/Dropbox/GitHub/wsg-choke-species")

#Load functions
source("code/helper_funs.R")

#Output folder
output_folder <- "region_comp"

#Plot themes
theme_set(theme_bw(base_size = 15))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

#Load data
dat <- list.files(path = "data/processed_data/fish2", pattern = ".rds", full.names=T) %>%
  map(readRDS) %>% 
  bind_rows()

#Remove AI
dat <- filter(dat, region!="ai")
dat <- filter(dat, region!="nbs")

#Combine IPHC and bottom trawl catch data
dat$catch_weight_combined <- ifelse(is.na(dat$catch_weight), dat$cpue_weight, dat$catch_weight)
dat$catch_count_combined <- ifelse(is.na(dat$catch_numbers), dat$cpue_count, dat$catch_numbers)

#Clean up
dat$catch_weight_combined <- replace(dat$catch_weight_combined, dat$catch_weight_combined == "Inf", NA)
dat$catch_count_combined <- replace(dat$catch_count_combined, dat$catch_count_combined == "Inf", NA)
dat$catch_weight_combined <- replace(dat$catch_weight_combined, dat$catch_weight_combined == "NaN", NA)
dat$catch_count_combined <- replace(dat$catch_count_combined, dat$catch_count_combined == "NaN", NA)

##Get depth habitat for filtering
species_table <- read_excel("data/species_table.xlsx")
species_table$common_name <- tolower(species_table$common_name)
species <- unique(species_table$common_name)
#IPHC species that use counts and not weights
species_iphc <- c("sablefish", "pacific cod", "yelloweye rockfish", "longnose skate", "big skate", "spiny dogfish", "rougheye rockfish")

bottom_trawl_only <- F

for(i in 1:length(species)){
  this_species <- species[i]
  dat_depth2 <- filter(dat, common_name==this_species)
  print(this_species)
  
  # Sort by depth
  dat_depth2 <- dat_depth2[order(dat_depth2$depth), ]
  
  #Calculate the cumulative sum of catch by depth
  if(this_species %in% species_iphc){
    dat_depth2$catch2use <- dat_depth2$catch_count_combined
  } else {
    dat_depth2$catch2use <- dat_depth2$catch_weight_combined
  } 
  if(bottom_trawl_only){
    dat_depth2$catch2use <- dat_depth2$catch_weight
  }
  
  dat_depth2 <- drop_na(dat_depth2, catch2use)
  dat_depth2$cumsum_catch <- cumsum(dat_depth2$catch2use)

  #Calculate the proportional cumulative sum
  dat_depth2$prop_cumsum_catch <- dat_depth2$cumsum_catch / sum(dat_depth2$catch2use, na.rm=T)
  
  # Find the index of the closest value to 99% (0.99) in prop_cumsum_var1
  closest_index <- which.min(abs(dat_depth2$prop_cumsum_catch - 0.99))
  
  # Get the depth at 99% catch
  closest_value <- dat_depth2$depth[closest_index]
  
  #Add to dataframe
  dat_depth2$filtered_depth <- closest_value
  
  #Dataframe
  depth2use <- data.frame(common_name=this_species, depth=closest_value)
  
  if(i==1){
    depths <- depth2use
    dat_depths <- dat_depth2
  } else {
    depths <- bind_rows(depths, depth2use)
    dat_depths <- bind_rows(dat_depths, dat_depth2)
  }
}

#Add depths to species table
species_table <- left_join(species_table, depths, by="common_name")
write_xlsx(species_table, "data/species_table.xlsx")

#Plot cumulative sum by depth
ggplot(dat_depths, aes(x=depth, y=prop_cumsum_catch))+
  geom_line()+
  geom_vline(species_table, mapping=aes(xintercept=depth), linetype="dashed")+
  facet_wrap("common_name", ncol=4, labeller=labeller(common_name=label_wrap_gen(15)))+
  ylab("Cumulative Sum of Catch")+
  xlab("Depth (m)")

ggsave(
  paste0("output/", output_folder, "/plots/cumulative_depth.png"),
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

###Filter for latitude
#Regions
regions <- c("cc", "bc", "goa", "ebs")

#Latitude
for(i in 1:length(species)){
  this_species <- species[i]
  dat2 <- filter(dat, common_name==this_species)
  print(this_species)
  
  #Select catch data to use
  if(this_species %in% species_iphc){
    dat2$catch2use <- dat2$catch_count_combined
  } else {
    dat2$catch2use <- dat2$catch_weight_combined
  } 
  if(bottom_trawl_only){
    dat2$catch2use <- dat2$catch_weight
  }
  
  #Calculate the cumulative sum of catch by latitude in each region
  for(j in 1:length(regions)){
  print(regions[j])
  dat_depth2 <- filter(dat2, region==regions[j])
  
  # Sort by latitude
  if(regions[j]=="cc"){
  dat_depth2 <- dat_depth2[order(dat_depth2$latitude, decreasing=T), ]
  } else {
    dat_depth2 <- dat_depth2[order(dat_depth2$latitude, decreasing=F),]
  }
  #drop NA
  dat_depth2 <- drop_na(dat_depth2, catch2use)
  #Calculate cumulative sum
  dat_depth2$cumsum_catch <- cumsum(dat_depth2$catch2use)
  
  #If species is just completely not present in a region, skip
  pos_catch <- filter(dat_depth2, catch2use>0)
  if(sum(dat_depth2$cumsum_catch)>0 &nrow(pos_catch)>50){
  #Calculate the proportional cumulative sum
  dat_depth2$prop_cumsum_catch <- dat_depth2$cumsum_catch / sum(dat_depth2$catch2use, na.rm=T)
  
  # Find the index of the closest value to 9% (0.99) in prop_cumsum_var1
  closest_index <- which.min(abs(dat_depth2$prop_cumsum_catch - 0.99))
  
  # Get the depth at 99% catch
  closest_value <- dat_depth2$latitude[closest_index]
  
  #Add to dataframe
  dat_depth2$filtered_lat <- closest_value
  } else {
    closest_value <- NA
    dat_depth2$filtered_lat <- NA
  }
  
  #Dataframe
  depth2use <- data.frame(common_name=this_species, latitude=closest_value, region=regions[j])
  
  if(j==1){
    lats <- depth2use
    dat_lats <- dat_depth2
  } else {
    lats <- bind_rows(lats, depth2use)
    dat_lats <- bind_rows(dat_lats, dat_depth2)
  }
  }
  if(i==1){
    lats_filtered <- lats
    dat_lats2 <- dat_lats
  } else {
    lats_filtered <- bind_rows(lats_filtered, lats)
    dat_lats2 <- bind_rows(dat_lats2, dat_lats)
  }
}

#Plot
dat_lats2$region <- factor(dat_lats2$region, levels=c("ebs", "goa", "bc", "cc"))
lats_filtered$region <- factor(lats_filtered$region, levels=c("ebs", "goa", "bc", "cc"))
labs <- c("Eastern Bering Sea", "Gulf of Alaska", "British Columbia", "California Current")
names(labs) <- c("ebs", "goa", "bc", "cc")

ggplot(dat_lats2, aes(x=latitude, y=prop_cumsum_catch))+
  geom_line()+
  geom_vline(lats_filtered, mapping=aes(xintercept=latitude), linetype="dashed")+
  facet_grid(common_name~region, labeller=labeller(common_name=label_wrap_gen(15)))+
  ylab("Cumulative Sum of Catch")+
  xlab("Latitude")

ggsave(
  paste0("output/", output_folder, "/plots/cumulative_latitude.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 25,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##Just CC
ggplot(filter(dat_lats2,region=="cc"), aes(x=latitude, y=prop_cumsum_catch))+
  scale_x_reverse()+
  geom_line()+
  geom_vline(filter(lats_filtered, region=="cc"), mapping=aes(xintercept=latitude), linetype="dashed")+
  facet_wrap("common_name", ncol=4, labeller=labeller(common_name=label_wrap_gen(15)))+
  ylab("Cumulative Sum of Catch")+
  xlab("Latitude")+
  ggtitle("California Current")

ggsave(
  paste0("output/", output_folder, "/plots/cumulative_latitude_CC.png"),
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

##Just EBS
lats_filtered2 <- filter(lats_filtered, region=="ebs"&!is.na(latitude))
dat_lats3 <- filter(dat_lats2, region=="ebs" & (common_name %in% unique(lats_filtered2$common_name)))

ggplot(dat_lats3, aes(x=latitude, y=prop_cumsum_catch))+
  geom_line()+
  geom_vline(lats_filtered2, mapping=aes(xintercept=latitude), linetype="dashed")+
  facet_wrap("common_name", ncol=4, labeller=labeller(common_name=label_wrap_gen(15)))+
  ylab("Cumulative Sum of Catch")+
  xlab("Latitude")+
  ggtitle("Eastern Bering Sea")

ggsave(
  paste0("output/", output_folder, "/plots/cumulative_latitude_EBS.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 7,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

#Add to species table
cc <- filter(lats_filtered, region=="cc")
species_table$southern_limit <- cc$latitude
ebs <- filter(lats_filtered, region=="ebs")
species_table$northern_limit <- ebs$latitude

#Manually remove ones without
no_s_limits <- c('longspine thornyhead', 'shortspine thornyhead', 'slender sole', 'spotted ratfish', 'southern rock sole', 'blackbelly eelpout','dover sole', 'english sole', 'lingcod', 'longnose skate', 'pacific hake', 'pacific sanddab', 'petrale sole', 'sablefish')
species_table$southern_limit <- ifelse(species_table$common_name %in% no_s_limits, NA, species_table$southern_limit)

#Remove ebs limit from 
no_n_limits <- c("longspine thornyhead", "walleye pollock", "arrowtooth flounder", "flathead sole", "pacific cod", "pacific halibut", "sandpaper skate")
species_table$northern_limit <- ifelse(species_table$common_name %in% no_n_limits, NA, species_table$northern_limit)

#Save
write_xlsx(species_table, "data/species_table.xlsx")

##Plot only ones with clear limits
lats_filtered2 <- filter(species_table, !is.na(northern_limit))
dat_lats3 <- filter(dat_lats2, region=="ebs" & (common_name %in% unique(lats_filtered2$common_name)))

a <- ggplot(dat_lats3, aes(x=latitude, y=prop_cumsum_catch))+
  geom_line()+
  geom_vline(lats_filtered2, mapping=aes(xintercept=northern_limit), linetype="dashed")+
  facet_wrap("common_name", ncol=4, labeller=labeller(common_name=label_wrap_gen(15)))+
  ylab("Cumulative Sum of Catch")+
  xlab("Latitude")+
  ggtitle(" Species with Northern Limits (Eastern Bering Sea)")

ggsave(
  paste0("output/", output_folder, "/plots/cumulative_latitude_EBS_only.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 4,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##Southern limits Plot only ones with clear limits
lats_filtered2 <- filter(species_table, !is.na(southern_limit))
dat_lats3 <- filter(dat_lats2, region=="cc" & (common_name %in% unique(lats_filtered2$common_name)))

b <- ggplot(dat_lats3, aes(x=latitude, y=prop_cumsum_catch))+
  geom_line()+
scale_x_reverse()+
  geom_vline(lats_filtered2, mapping=aes(xintercept=southern_limit), linetype="dashed")+
  facet_wrap("common_name", ncol=4, labeller=labeller(common_name=label_wrap_gen(15)))+
  ylab("Cumulative Sum of Catch")+
  xlab("Latitude")+
  ggtitle("Species with Southern Limits (California Current)")

ggsave(
  paste0("output/", output_folder, "/plots/cumulative_latitude_CC_only.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 7,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

#Combine
ggarrange(a,b,labels=c("A", "B"), nrow=2, heights=c(0.6,1), widths=c(0.8,1))
ggsave(
  paste0("output/", output_folder, "/plots/cumulative_latitude_combined.png"),
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

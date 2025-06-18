library(sdmTMB)
library(dplyr)
library(Metrics)
library(ggplot2)
library(tidyr)
library(visreg)
library(ggpubr)
library(purrr)
library(readxl)
library(openxlsx2)
library(ggpattern)
setwd("~/Dropbox/GitHub/wsg-choke-species")

#Load functions
source("code/helper_funs.R")

#Load data
dat <- list.files(path = "data/processed_data/fish2", pattern = ".rds", full.names=T) %>%
  map(readRDS) %>% 
  bind_rows()

##Clean up data for plotting

#Remove any rows with necessary data missing
dat <- dat %>%
  drop_na(depth, mi1, temperature_C, salinity_psu, X, Y, year)

#Remove weird depths
dat <- filter(dat, depth>0)

#Remove oxygen outliers
dat <- filter(dat, O2_umolkg<1500)

#Add survey type
dat$survey_type <- ifelse(dat$survey=="iphc", "IPHC longline", "bottom trawl survey")

#Remove aleutian islands
dat <- filter(dat, region!="ai")

#Remove NBS
dat <- filter(dat, region!="nbs")

#Remove duplicates
dat <- unique(dat)

###Sanity checks
##Just bottom trawl surveys
test <- subset(dat, dat$survey_type!="IPHC longline")

# Check that there are some zeros
out2 = tapply( test[,'catch_weight'], INDEX=list(test[,'region'],test[,'common_name']), FUN=length )
#PASSES

## For IPHC
test <- subset(dat, dat$survey_type=="IPHC longline")
out2 = tapply( test[,'cpue_weight'], INDEX=list(test[,'region'],test[,'common_name']), FUN=length )
#PASSES

#Region labels
#Set region
dat$region <- factor(dat$region, levels=c("ebs", "goa", "bc", "cc"))
labs <- c("Eastern Bering Sea", "Gulf of Alaska", "British Columbia", "California Current")
names(labs) <- c("ebs", "goa", "bc", "cc")

#Plot themes
theme_set(theme_bw(base_size = 15))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())
##Plot barplot of all data available--for diagnostics
ggplot(filter(dat, survey!="iphc"&catch_weight>0), aes(x=year))+
  stat_count(aes(fill=region))+  
  facet_wrap("common_name", ncol=4, scales="free_y")+
  scale_x_continuous(breaks=c(2000,2010,2020), limits=c(2000,2027))+
  xlab("Year")+
  ylab("Number of Observations")+
  theme(legend.position="top")+
  guides(fill = guide_legend(nrow = 1))

ggplot(filter(dat, survey=="iphc"), aes(x=year))+
  stat_count(aes(fill=region))+  
  facet_wrap("common_name", ncol=2, scales="free_y")+
  scale_x_continuous(breaks=c(2000,2010,2020), limits=c(2000,2027))+
  xlab("Year")+
  ylab("Number of Observations")+
  theme(legend.position="top")+
  guides(fill = guide_legend(nrow = 1))

##Combine IPHC and bottom trawl surveys, positive catch only
##Plot barplot with IPHC and bottom trawl
ggplot(filter(dat, (survey!="iphc" & catch_weight>0)|(survey=="iphc"&(cpue_weight>0|cpue_count>0))), aes(x=year, fill=region, pattern=survey_type))+
  #stat_count(aes(fill=region, pattern=survey_type))+
  facet_wrap("common_name", ncol=4, scales="free_y")+# labeller=labeller(common_name=label_wrap_gen(10)))+
  scale_x_continuous(breaks=c(2008,2016,2024), limits=c(2008,2026))+
  geom_bar_pattern(
    stat = "count",
    colour = "black",            # Outline color
    pattern_fill = "black",      # Pattern stripe color
    pattern_angle = 45,
    pattern_density = 0.2,
    pattern_spacing = 0.05,
    pattern_size = 0.1)+
  scale_pattern_manual(
    name="data type",
    values = c("none", "stripe"))+
  xlab("Year")+
  ylab("Number of Observations")+
  theme(legend.position="top")+
 guides(fill = guide_legend(nrow = 2), pattern=guide_legend(nrow=2,override.aes = list(pattern = c("none", "stripe"))))+
  scale_fill_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677"), drop=FALSE, labels=labs)

ggsave(
  paste("output/plots/dat_availability/barplot_positiveonly_all.png"),
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

##With zeroes
ggplot(filter(dat, (survey!="iphc" & catch_weight>0)|(survey=="iphc"& (cpue_count>0|cpue_weight>0))), aes(x=year))+
  stat_count(dat, mapping=aes(x=year), fill="grey", alpha=0.3)+
  stat_count(aes(alpha=survey_type, fill=region))+
  facet_wrap("common_name", ncol=4, scales="free_y")+
  scale_x_continuous(breaks=c(2008,2016,2024), limits=c(2008,2026))+
  xlab("Year")+
  ylab("Number of Observations")+
  theme(legend.position="top")+
  guides(fill = guide_legend(nrow = 1), alpha="none")+
  scale_alpha_manual(values=c(0.5,1))+
  theme(strip.text=element_text(size=10))


ggsave(
  paste("output/plots/dat_availability/barplot_all_zeroes.png"),
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

##With zeroes--bottom trawl only
ggplot(filter(dat, (survey!="iphc" & catch_weight>0)), aes(x=year))+
  stat_count(filter(dat, (survey!="iphc")), mapping=aes(x=year), fill="grey", alpha=0.3)+
  stat_count(aes(fill=region))+
  facet_wrap("common_name", ncol=4, scales="free_y")+
  scale_x_continuous(breaks=c(2008,2016,2024), limits=c(2008,2026))+
  xlab("Year")+
  ylab("Number of Observations")+
  theme(legend.position="top")+
  guides(fill = guide_legend(nrow = 1), alpha="none")+
  scale_alpha_manual(values=c(0.5,1))


ggsave(
  paste("output/plots/dat_availability/barplot_all_zeroes_bottom_trawl_only.png"),
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

###Map of all data available
# setup up mapping ####
map_data <- rnaturalearth::ne_countries(scale = "large",
                                        returnclass = "sf",
                                        continent = "North America")

us_coast_proj <- sf::st_transform(map_data, crs = 32610)

###Map of data available
species <- unique(dat$common_name)
for(i in 1:length(species)){
dat2plot <- filter(dat, common_name==species[i]&((survey!="iphc" & catch_weight>0)|(survey=="iphc"& (cpue_count>0|cpue_weight>0))))
ggplot(us_coast_proj) + geom_sf() +
  xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
  ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
  theme(axis.text.x=element_blank())+
  geom_point(dat2plot, mapping=aes(x=X*1000, y=Y*1000,colour=survey))+
  facet_wrap("year", ncol=7)+
  theme_minimal(base_size=12)+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle(paste(unique(dat2plot$common_name)))
ggsave(
  paste("output/plots/dat_availability/map_",species[i],".png"),
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

##All combined, example
example <- filter(dat, common_name=="pacific halibut")
ggplot(example, aes(x=year))+
  stat_count(aes(alpha=survey_type, fill=region))+
  scale_x_continuous(breaks=c(2008,2016,2024), limits=c(2008,2024))+
  xlab("Year")+
  ylab("Number of Observations")+
  theme(legend.position="top")+
  guides(fill = guide_legend(nrow = 1))+
  scale_alpha_manual(values=c(0.4,1))+
guides(fill = guide_legend(nrow = 1), alpha="none")

ggsave(
  paste("output/plots/data_available_overall.png"),
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

##Make a table
#Number of positive catches per species per region
table <- dat %>% filter(survey_type!="iphc" & catch_weight>0)%>% 
  group_by(common_name, region) %>% 
  summarize(N_pos = length(catch_weight))
table <- pivot_wider(table, names_from=region, id_cols=common_name, values_from=c(N_pos))
#Species table
species <- read_excel("data/species_table.xlsx")
species$common_name <- tolower(species$common_name)
#Combine
species <- left_join(species, table)

##IPHC data
table2 <- dat %>% filter(survey=="iphc"& (cpue_count>0|cpue_weight>0))%>% 
  mutate(cpue=ifelse(is.na(cpue_count), cpue_weight, cpue_count))%>% 
  group_by(common_name, region) %>% 
  summarize(N_pos_IPHC = length(cpue))
  
table2 <- pivot_wider(table2, names_from=region, id_cols=common_name, values_from=c(N_pos_IPHC))
colnames(table2) <- c("common_name", "iphc_cc", "iphc_bc", "iphc_goa", "iphc_ebs")

species <- left_join(species, table2)

#save
write.csv(species, file="output/species_summary.csv")

##Get depth habitat for filtering
species <- read_excel("data/species_table.xlsx")
species$common_name <- tolower(species$common_name)
species <- unique(species$common_name)
for(i in 1:length(species)){
  this_species <- species[i]
  dat_depth2 <- filter(dat, common_name==this_species)
  dat_depth2 <- dat_depth2  %>% drop_na(depth, catch_weight)
  
  # Sort by depth
  dat_depth2 <- dat_depth2[order(dat_depth2$depth), ]
  
  #Calculate the cumulative sum of catch by depth
  dat_depth2$cumsum_catch <- cumsum(dat_depth2$catch_weight)
  
  #Calculate the proportional cumulative sum
  dat_depth2$prop_cumsum_catch <- dat_depth2$cumsum_catch / sum(dat_depth2$catch_weight, na.rm=T)
  
  # Find the index of the closest value to 98% (0.97) in prop_cumsum_var1
  closest_index <- which.min(abs(dat_depth2$prop_cumsum_catch - 0.99))
  
  # Get the depth at 97% catch
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

#Species table
species <- read_excel("data/species_table.xlsx")
species$common_name <- tolower(species$common_name)
species <- left_join(species, depths)
write_xlsx(species, "data/species_table.xlsx")

#Plot cumulative sum by depth
ggplot(dat_depths, aes(x=depth, y=prop_cumsum_catch))+
  geom_line()+
  geom_vline(species, mapping=aes(xintercept=depth), linetype="dashed")+
  facet_wrap("common_name", ncol=4)+
  ylab("Cumulative Sum of Catch")+
  xlab("Depth (m)")

ggsave(
  paste("output/plots/dat_availability/cumulative_depth.png"),
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

##Models fit to data
#Plot
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

##Combine IPHC and bottom trawl surveys, positive catch only
##Plot barplot with IPHC and bottom trawl
ggplot(filter(dat, (survey!="iphc" & catch_weight>0)|(survey=="iphc"&(cpue_weight>0|cpue_count>0))), aes(x=year, fill=region, pattern=survey_type))+
  #stat_count(aes(fill=region, pattern=survey_type))+
  facet_wrap("common_name", ncol=4, scales="free_y")+# labeller=labeller(common_name=label_wrap_gen(10)))+
  scale_x_continuous(breaks=c(2009,2016,2023), limits=c(2008,2023))+
  geom_bar_pattern(
    stat = "count",
    colour = "black",            # Outline color
    pattern_fill = "black",      # Pattern stripe color
    pattern_angle = 45,
    pattern_density = 0.2,
    pattern_spacing = 0.05,
    pattern_size = 0.1)+
  scale_pattern_manual(
    name="data type",
    values = c("none", "stripe"))+
  xlab("Year")+
  ylab("Number of Observations")+
  theme(legend.position="top")+
  guides(fill = guide_legend(nrow = 2), pattern=guide_legend(nrow=2,override.aes = list(pattern = c("none", "stripe"))))+
  scale_fill_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677"), drop=FALSE, labels=labs)

ggsave(
  paste("output/plots/dat_availability/barplot_positiveonly_all.png"),
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

##All combined, example
example <- filter(dat, common_name=="sablefish")

ggplot(example, aes(x=year, fill=region, pattern=survey_type))+
  geom_bar_pattern(
    stat = "count",
    colour = "black",            # Outline color
    pattern_fill = "black",      # Pattern stripe color
    pattern_angle = 45,
    pattern_density = 0.2,
    pattern_spacing = 0.02,
    pattern_size = 0.1)+
  scale_pattern_manual(
    name="data type",
    values = c("none", "stripe"))+
  scale_x_continuous(breaks=c(2008,2016,2024))+
  xlab("Year")+
  ylab("Number of Observations")+
  theme(legend.position="top")+
  guides(fill = guide_legend(nrow = 2), pattern=guide_legend(nrow=2,override.aes = list(pattern = c("none", "stripe"))))+
  scale_fill_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677"), drop=FALSE, labels=labs)

ggsave(
  paste("output/plots/data_available_overall.png"),
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


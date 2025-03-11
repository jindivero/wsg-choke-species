library(sdmTMB)
library(dplyr)
library(Metrics)
library(ggplot2)
library(tidyr)
library(visreg)
library(ggpubr)
library(purrr)
library(readxl)
setwd("~/Dropbox/GitHub/wsg-choke-species")

#Load functions
source("code/helper_funs.R")

#Load data
dat <- list.files(path = "data/processed_data/fish", pattern = ".rds", full.names=T) %>%
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
dat$survey_type <- ifelse(dat$survey=="iphc", "IPHC longline", "NOAA/DFO bottom trawl survey")

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
out2 = tapply( test[,'cpue_kg_km2'], INDEX=list(test[,'region'],test[,'common_name']), FUN=length )
#PASSES

# Check that there's some zeros
out = tapply( test[,'cpue_kg_km2'], INDEX=list(test[,'region'],test[,'common_name']), FUN=function(x){sum(x==0)} )

## For IPHC
test <- subset(dat, dat$survey_type=="IPHC longline")
out2 = tapply( test[,'cpue_kg_km2'], INDEX=list(test[,'region'],test[,'common_name']), FUN=length )
#PASSES

#Region labels
labs <- c("British Columbia", "California Current", "Eastern Bering Sea", "Gulf of Alaska")
names(labs) <- c("bc", "cc", "ebs", "goa")
dat$region <- factor(dat$region, levels=c("cc", "bc", "goa", "ebs")) 

#Plot themes
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())
##Plot barplot of all data available--for diagnostics
ggplot(filter(dat, survey!="iphc"&cpue_kg_km2>0), aes(x=year))+
  stat_count(aes(fill=region))+  
  facet_wrap("common_name", ncol=2, scales="free_y")+
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
##Plot barplot with IPHC and non
ggplot(filter(dat, (survey!="iphc" & cpue_kg_km2>0)|(survey=="iphc"&(cpue_weight>0|cpue_count>0))), aes(x=year))+
  stat_count(aes(alpha=survey_type, fill=region))+
  facet_wrap("common_name", ncol=4, scales="free_y")+
  scale_x_continuous(breaks=c(2008,2016,2024), limits=c(2008,2026))+
  xlab("Year")+
  ylab("Number of Observations")+
  theme(legend.position="top")+
 guides(fill = guide_legend(nrow = 1), alpha="none")+
  scale_alpha_manual(values=c(0.4,1))

ggsave(
  paste("output/plots/dat_availability/barplot_all.png"),
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
ggplot(filter(dat, (survey!="iphc" & cpue_kg_km2>0)|(survey=="iphc"& (cpue_count>0|cpue_weight>0))), aes(x=year))+
  stat_count(dat, mapping=aes(x=year), fill="grey", alpha=0.3)+
  stat_count(aes(alpha=survey_type, fill=region))+
  facet_wrap("common_name", ncol=4, scales="free_y")+
  scale_x_continuous(breaks=c(2008,2016,2024), limits=c(2008,2026))+
  xlab("Year")+
  ylab("Number of Observations")+
  theme(legend.position="top")+
  guides(fill = guide_legend(nrow = 1), alpha="none")+
  scale_alpha_manual(values=c(0.5,1))


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
ggplot(filter(dat, (survey!="iphc" & cpue_kg_km2>0)), aes(x=year))+
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
dat2plot <- filter(dat, common_name==species[i]&((survey!="iphc" & cpue_kg_km2>0)|(survey=="iphc"& (cpue_count>0|cpue_weight>0))))
ggplot(us_coast_proj) + geom_sf() +
  xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
  ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
  geom_point(dat2plot, mapping=aes(x=X*1000, y=Y*1000,colour=survey))+
  facet_wrap("year", ncol=7)+
  theme_minimal(base_size=20)+
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
table <- dat %>% filter(survey_type!="iphc" & cpue_kg_km2>0)%>% 
  group_by(common_name, region) %>% 
  summarize(N_pos = length(cpue_kg_km2))
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

#Range
temp_ranges <- dplyr::group_by(dat, species, region) %>%
  dplyr::filter(cpue_kg_km2 > 0) %>%
  dplyr::summarize(min_temp = min(temp,na.rm=T), max_temp = max(temp,na.rm=T))
saveRDS(temp_ranges, "outputs/MI_ranges.rds")
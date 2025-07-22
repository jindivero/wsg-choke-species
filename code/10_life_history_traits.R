# look up size and other life history traits 
install.packages("rfishbase")
devtools::install_github("james-thorson/FishLife", dep=TRUE)
library(FishLife)
library(rfishbase)
library(readxl)
library(dplyr)
library(ggplot2)

#Species list
setwd("~/Dropbox/GitHub/wsg-choke-species")

#ggplot themes
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

#Read in table with species names
species_table <- read_excel("data/species_table.xlsx")
sp_list <- species_table$scientific_name

# fix some names
sp_list <- gsub("sebastes ruberrimus", "Sebastes ruberrimus" , sp_list)
sp_list <- gsub("Raja binoculata", "Beringraja binoculata" , sp_list)
sp_list <- gsub("Raja rhina", "Beringraja rhina" , sp_list)
sp_list <- gsub("sebastes aleutianus", "Sebastes aleutianus" , sp_list)


## fishlife

# get species names and tip order from FishLife
tip_order_fl = FishBase_and_Morphometrics$tree$tip.label 

tips <- match(sp_list, tip_order_fl)

# get trait values
traits <- FishBase_and_Morphometrics$beta_gv[tips,]
species <- rownames(traits)

traits <- bind_cols(species, traits) |> 
  rename(species = '...1') 

# looks like Bathyraja kincaidii is missing from FishLife
# add back in for later merging with Fishbase
traits$species[is.na(traits$species)] <- "Bathyraja kincaidii"


# can select traits/columns of interest (e.g, 'log(growth coefficient)')

# aspect ratio is in here - could be interesting to see if there is anything there


## fishbase ####

# Get Winf, Linf, K

# get species codes first
#Get Fishbase species codes

sp_codes <- fb_tbl("species") |>
  select(SpecCode, Genus, Species)

sp_codes <- sp_codes |>
  mutate(genus_species = paste0(Genus, " ", Species))

sp_codes <- sp_codes |>
  filter(genus_species %in% sp_list)

setdiff(sp_list, sp_codes$genus_species)

# Beringraja rhina not in here...

growth <- fb_tbl("popgrowth") |>
  filter(SpecCode %in% sp_codes$SpecCode)

length(unique(growth$SpecCode))

# 28 species here

# which are missing?
growth <- left_join(growth, sp_codes)

setdiff(sp_list, growth$genus_species)

# filter by FAO area
fao_areas <- faoareas(sp_codes$genus_species) |>
  select(AreaCode, StockCode, SpecCode)

growth <- left_join(growth, fao_areas) |>
  filter(AreaCode %in% c(67, 77)) |> # NEP FAO areas
  select(StockCode, SpecCode, genus_species,
         Loo, Winfinity, K, Locality, AreaCode) 

# mean by species
LH_sum <- growth |>
  group_by(genus_species) |>
  summarise(mean_K =   mean(K, na.rm = TRUE),
            mean_Loo = mean(Loo,  na.rm = TRUE),
            mean_Woo = mean(Winfinity,  na.rm = TRUE)
  )

# which species have NA for Winfinity?
missing_sp <- LH_sum %>%
  filter(is.nan(mean_Woo)) %>%
  pull(genus_species)

missing_sp


# put together

FL_traits <- traits |>
  select(
    species,
    `log(age_max)`,`log(aspect_ratio)`,
    `log(fecundity)`, `log(growth_coefficient)`,
    `log(length_infinity)`, `log(length_max)`,
    `log(length_maturity)`, `log(age_maturity)`,
    `log(natural_mortality)`, `log(weight_infinity)`) |>
  mutate(across(where(is.numeric), exp)) |>
  rename(
    max_age =           `log(age_max)`,
    aspect_ratio =      `log(aspect_ratio)`,
    fecundity =         `log(fecundity)`,
    growth_coefficient =`log(growth_coefficient)`,
    length_infinity =   `log(length_infinity)`,
    max_length =        `log(length_max)`,
    length_maturity =   `log(length_maturity)`,
    age_maturity =      `log(age_maturity)`,
    natural_mortality = `log(natural_mortality)`,
    weight_infinity =   `log(weight_infinity)`
  )

FB_traits <- LH_sum |>
  rename(species = genus_species)

all_traits <- left_join(FL_traits, FB_traits)
	
	#Add oxygen threshold estimates
	bp_est2 <- readRDS(paste0("output/region_comp/breakpoint_estimates_filtered.rds"))
	#Add scientific name
	bp_est2 <- left_join(species_table, bp_est2, by=c("common_name"="species"))
	
	#Combine
	summary <- left_join(all_traits, bp_est2, by=c("species"="scientific_name"))
	summary$threshold <- ifelse(is.na(summary$bp_ensemble_mean), 0, 1)
	
a <-	ggplot(summary, aes(x = as.factor(threshold), y=aspect_ratio)) +
	  geom_boxplot(aes(colour=as.factor(threshold)))+
  geom_point(aes(colour=as.factor(threshold)))+
	  labs(colour="Oxygen Limitation Threshold")+
    scale_colour_discrete(labels=c("No Threshold", "Threshold"))+
	  scale_x_discrete(labels=c("No Threshold", "Threshold"))+
	  xlab("")
	
b <- 	ggplot(summary, aes(x = as.factor(threshold), y=weight_infinity)) +
	  geom_boxplot(aes(colour=as.factor(threshold)))+
  geom_point(aes(colour=as.factor(threshold)))+
	  labs(colour="Oxygen Limitation Threshold")+
	  scale_x_discrete(labels=c("No Threshold", "Threshold"))+
  scale_colour_discrete(labels=c("No Threshold", "Threshold"))+
	  xlab("")+
  coord_cartesian(ylim=c(0, 15000))
	
c <- 	ggplot(summary, aes(x = as.factor(threshold), y=max_age)) +
	  geom_boxplot(aes(colour=as.factor(threshold)))+
  geom_point(aes(colour=as.factor(threshold)))+
	  labs(colour="Oxygen Limitation Threshold")+
	  scale_x_discrete(labels=c("No Threshold", "Threshold"))+
  scale_colour_discrete(labels=c("No Threshold", "Threshold"))+
	  xlab("")
	
d <- 	ggplot(summary, aes(x = as.factor(threshold), y=length_infinity)) +
    geom_boxplot(aes(colour=as.factor(threshold)))+
  geom_point(aes(colour=as.factor(threshold)))+
	  labs(colour="Oxygen Limitation Threshold")+
	  scale_x_discrete(labels=c("No Threshold", "Threshold"))+
  scale_colour_discrete(labels=c("No Threshold", "Threshold"))+
	  xlab("")

e <- 	ggplot(summary, aes(x = as.factor(threshold), y=max_length)) +
  geom_boxplot(aes(colour=as.factor(threshold)))+
  geom_point(aes(colour=as.factor(threshold)))+
  labs(colour="Oxygen Limitation Threshold")+
  scale_x_discrete(labels=c("No Threshold", "Threshold"))+
  scale_colour_discrete(labels=c("No Threshold", "Threshold"))+
  xlab("")

f <- 	ggplot(summary, aes(x = as.factor(threshold), y=growth_coefficient)) +
  geom_boxplot(aes(colour=as.factor(threshold)))+
  geom_point(aes(colour=as.factor(threshold)))+
  labs(colour="Oxygen Limitation Threshold")+
  scale_x_discrete(labels=c("No Threshold", "Threshold"))+
  scale_colour_discrete(labels=c("No Threshold", "Threshold"))+
  xlab("")

g <- 	ggplot(summary, aes(x = as.factor(threshold), y=fecundity)) +
  geom_boxplot(aes(colour=as.factor(threshold)))+
  geom_point(aes(colour=as.factor(threshold)))+
  labs(colour="Oxygen Limitation Threshold")+
  scale_x_discrete(labels=c("No Threshold", "Threshold"))+
  scale_colour_discrete(labels=c("No Threshold", "Threshold"))+
  xlab("")

h <- 	ggplot(summary, aes(x = as.factor(threshold), y=age_maturity)) +
  geom_boxplot(aes(colour=as.factor(threshold)))+
  geom_point(aes(colour=as.factor(threshold)))+
  labs(colour="Oxygen Limitation Threshold")+
  scale_x_discrete(labels=c("No Threshold", "Threshold"))+
  scale_colour_discrete(labels=c("No Threshold", "Threshold"))+
  xlab("")

i <- 	ggplot(summary, aes(x = as.factor(threshold), y=natural_mortality)) +
  geom_boxplot(aes(colour=as.factor(threshold)))+
  geom_point(aes(colour=as.factor(threshold)))+
  labs(colour="Oxygen Limitation Threshold")+
  scale_x_discrete(labels=c("No Threshold", "Threshold"))+
  scale_colour_discrete(labels=c("No Threshold", "Threshold"))+
  xlab("")

k <- 	ggplot(summary, aes(x = as.factor(threshold), y=mean_Woo)) +
  geom_boxplot(aes(colour=as.factor(threshold)))+
  geom_point(aes(colour=as.factor(threshold)))+
  labs(colour="Oxygen Limitation Threshold")+
  scale_x_discrete(labels=c("No Threshold", "Threshold"))+
  scale_colour_discrete(labels=c("No Threshold", "Threshold"))+
  xlab("")+
  coord_cartesian(ylim=c(0, 15000))


ggarrange(a, b, c, d, common.legend=T)
	
#save
ggsave(
  paste0("output/region_comp/plots/life_history_traits.png"),
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
	
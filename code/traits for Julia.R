# look up size and other life history traits 

  # devtools::install_github("james-thorson/FishLife", dep=TRUE)

  library(fishlife)
  library(fishbase)
  
  # species list
  path <- "/Users/jenniferbigman/Library/CloudStorage/Dropbox/Students/Indivero/species_list.rds"
  sp_list <- readRDS(path)
  
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
	
	sp_codes <- fb_tbl("species") |>
	  select(SpecCode, Genus, Species)
	
	sp_codes <- sp_codes |>
	  mutate(genus_species = paste0(Genus, " ", Species))
	
	sp_codes <- sp_codes |>
	  filter(genus_species %in% sp_list)
	
	setdiff(sp_list, sp_codes2$genus_species)
	
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
  
  
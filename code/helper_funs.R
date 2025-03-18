##Unit conversions and equations
calc_po2_sat <- function(salinity, temp, depth, oxygen, lat, long, umol_m3, ml_L) {
  # Input:       S = Salinity (pss-78)
  #              T = Temp (deg C) ! use potential temp
  #depth is in meters
  
  #Pena et al. ROMS and GLORYS are in mmol per m^3 (needs to be converted to umol per kg) (o2 from trawl data was in mL/L, so had to do extra conversions)
  gas_const = 8.31
  partial_molar_vol = 0.000032
  kelvin = 273.15
  boltz = 0.000086173324
  
  #this was dumb because it is actually equivalent to umol/kg
  #convert mmol to umol
  umol_m3 <- oxygen*1000
  #convert m3 to l
  umol_l <- umol_m3/1000
  #convert from molality (moles per volume) to molarity (moles per mass)
  #1 L of water = 1 kg of water, so no equation needed?? Right??
  o2_umolkg <- umol_l * 1/1
  
  SA = gsw_SA_from_SP(salinity,depth,long,lat) #absolute salinity for pot T calc
  pt = gsw_pt_from_t(SA,temp,depth) #potential temp at a particular depth
  #this is for if using data that has oxygen in ml/L
  #CT = gsw_CT_from_t(SA,temp,depth) #conservative temp
  #sigma0 = gsw_sigma0(SA,CT)
  # o2_umolkg = oxygen*44660/(sigma0+1000) 
  
  O2_Sat0 = gsw_O2sol_SP_pt(salinity,pt)
  
  #= o2satv2a(sal,pt) #uses practical salinity and potential temp - solubity at p =1 atm
  press = exp(depth*10000*partial_molar_vol/gas_const/(temp+kelvin))
  O2_satdepth = O2_Sat0*press
  
  #solubility at p=0
  sol0 = O2_Sat0/0.209
  sol_Dep = sol0*press
  po2 = o2_umolkg/sol_Dep
  po2 <- po2 * 101.325 # convert to kPa
  return(po2)
  
}

calc_mi <- function(Eo, Ao, W, n,po2, inv.temp) {
  mi = W^n*Ao*po2 *exp(Eo * inv.temp)
  return(mi)
}

calc_po2_crit <- function(inv.temp, taxa, mi, body_size, model) {
  W <- body_size
  inv.temp <- inv.temp
  ###Calculate Metabolic index 
  ##Species parameters from Tim's paper
  #Gadidae (lower, median, upper)
  if(taxa=="gadidae") {
    Eo <- c(-0.03064428,0.1883451,0.40414)
    V <- c(2.852883, 1.716267, 8.796332)
    n <- c(-0.1416557, 0.1883451, 0.05982743)
  }
  if(taxa=="perciformes"){
    Eo <- c(0.0267307, 0.3102310, 0.5859761)
    V <- c(1.734103, 1.59323, 8.836855)
    n <- c(-0.1720744, 0.3102311, 0.03412396)
  }
  if(taxa=="elasmobranchii"){
    #ie. Squalidae
    Eo <- c(-0.07659791, 0.2661157, 0.6013813)
    V <- c(0.7086717, 1.502317, 9.654109)
    n <-c(-0.2180308, 0.2661157, 0.05064831)
  }
  if(taxa=="teleostei"){
    #Orders Scorpaenidae, Soleida, Pleuronectida)
    Eo <- c(-0.009823544, 0.2554435,0.5315111)
    V <- c(1.639916,1.565438,8.815259)
    n <-c(-0.1562307,0.2554435,0.03221273)
  }
  Ao <- 1/exp(V[2])
  n <- n[2]
  Eo <-  if(model=="model3") Eo[1] else if(model=="model4") Eo[2] else Eo[3]
  po2 = mi/(W^n*Ao*exp(Eo * inv.temp))
  return(po2)
}

# calc o2 solubility, relies on o2 in umol/kg
gsw_O2sol_SP_pt <- function(sal,pt) {
  x = sal
  pt68 = pt*1.00024
  y = log((298.15 - pt68)/(273.15 + pt68))
  
  a0 =  5.80871
  a1 =  3.20291
  a2 =  4.17887
  a3 =  5.10006
  a4 = -9.86643e-2
  a5 =  3.80369
  b0 = -7.01577e-3
  b1 = -7.70028e-3
  b2 = -1.13864e-2
  b3 = -9.51519e-3
  c0 = -2.75915e-7
  
  O2sol = exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y)))) + x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x))
  return(O2sol)
}

back.convert <- function(x, mean_orig, sd_orig) {
  x* sd_orig+mean_orig
}

convert_class <- function(x) {
  for (i in 1:ncol(x)) x[,i] <- as(x[,i], Class = "matrix")
  return(x)
}

###DATA PROCESSING ###

print_species <- function(type){
  biomass <- combine_all(type)
  species <- biomass[1:2]
  return(species)
}

length_expand_nwfsc <- function(spc, sci_name) {
  # load, clean, and join data
  bio <- readRDS("data/fish_raw/NOAA/nwfsc_bio.rds")
  catch <- readRDS("data/fish_raw/NOAA/nwfsc_catch.rds")
  names(catch) = tolower(names(catch))
  names(bio) = tolower(names(bio))
  bio$scientific_name <- tolower(bio$scientific_name)
  bio$common_name <- tolower(bio$common_name)
  catch$common_name <- tolower(catch$common_name)
  
  bio$trawl_id = as.character(bio$trawl_id)
  catch$trawl_id=as.character(catch$trawl_id)
  
  #fix dogfish
  catch$common_name <- ifelse(catch$common_name=="pacific spiny dogfish", "spiny dogfish", catch$common_name)
  bio$common_name <- ifelse(bio$common_name=="pacific spiny dogfish", "spiny dogfish", bio$common_name)
  
  #haul$sampling_end_hhmmss = as.numeric(haul$sampling_end_hhmmss)
  #haul$sampling_start_hhmmss = as.numeric(haul$sampling_start_hhmmss)
  
  #Combine data
  dat = dplyr::left_join(filter(bio[,c("trawl_id", "year", "scientific_name", "common_name", "weight", "ageing_lab", "length_cm", "width_cm", "sex", "age")], !is.na(length_cm)), 
                         catch[,c("trawl_id","common_name", "subsample_count","area_swept_ha","longitude_dd", "latitude_dd","subsample_wt_kg","total_catch_numbers","total_catch_wt_kg","cpue_kg_km2")],
                         relationship = "many-to-many")
  
  
  # filter out species of interest from joined (catch/haul/bio) dataset
  dat_sub = dplyr::filter(dat, common_name == spc)
  
  # fit length-weight regression by year to predict fish weights that have lengths only.
  # note a rank-deficiency warning may indicate there is insufficient data for some year/sex combinations (likely for unsexed group)
  #Create one set of data with only year/survey combinations with at least some weight data:
  
  if(nrow(dat_sub)>0) {
    #Add column counting number of length observations per catch
    dat_sub$is_length <- ifelse(!is.na(dat_sub$length_cm), 1,0) 
    dat_sub <- group_by(dat_sub, trawl_id) %>% mutate(nlength=sum(is_length)) %>% ungroup()
    dat_sub$trawl_id <- as.numeric(dat_sub$trawl_id)
    
    fitted <-  filter(dat_sub, !is.na(length_cm)) %>%
      group_by(year) %>%
      mutate(sum = sum(weight, na.rm=T)) %>%
      filter(sum>0) %>%
      ungroup()
    fitted$weight <- ifelse(fitted$weight==0, NA, fitted$weight)
    #And one for those without any weight data:
    not_fitted <-  filter(dat_sub, !is.na(length_cm)) %>%
      group_by(year) %>%
      mutate(sum = sum(weight, na.rm=T)) %>%
      filter(sum==0) %>%
      ungroup()
    
    if(nrow(fitted)>0){
      fitted <-  
        group_nest(fitted, year) %>%
        mutate(
          model = purrr::map(data, ~ lm(log(weight) ~ log(length_cm), data = .x)),
          tidied = purrr::map(model, broom::tidy),
          augmented = purrr::map(model, broom::augment),
          predictions = purrr::map2(data, model, modelr::add_predictions)
        )
      
      # replace missing weights with predicted weights
      dat_pos2 = fitted %>%
        tidyr::unnest(predictions) %>%
        dplyr::select(-data, -model, -tidied, -augmented) %>%
        dplyr::mutate(weight = ifelse(is.na(weight), exp(pred), weight))
      
      dat_pos <- bind_rows(dat_pos2, not_fitted)
    }
    
    if(nrow(fitted)==0 & nrow(not_fitted)>0){
      dat_pos <- not_fitted
    }
    
    if(nrow(fitted)>0 | nrow(not_fitted)>0){
      #If there is no weight data available to get weight-length empirical interpolation, use fishbase to fill in
      # find length-weight relationship parameters for one species
      pars <- rfishbase::length_weight(
        #convert scientific name to first letter capitalized for fishbase
        species_list = stringr::str_to_sentence(sci_name, locale = "en"))
      #Get mean
      a <- mean(pars$a)
      b <- mean(pars$b)
      #Make NAs where text
      dat_pos$weight <- ifelse(dat_pos$weight=="NaN", NA, dat_pos$weight)
      #Remove any zero weights
      dat_pos$weight <- ifelse(dat_pos$weight==0, NA, dat_pos$weight)
      #Calculate weight (and convert from g to kg)
      dat_pos <- dplyr::mutate(dat_pos, weight = ifelse(is.na(weight), ((a*length_cm^b)*0.001), weight))
      #Add column getting the mean individual weight 
      dat_test <- group_by(dat_pos, trawl_id) %>% mutate(haul_weight=mean(weight)) %>% ungroup()
      
      trawlids <- unique(dat_pos$trawl_id)
      if(length(trawlids!=0)){
        p <- data.frame(trawl_id = trawlids,
                        p1 = 0,
                        p2 = 0,
                        p3 = 0,
                        p4 = 0)
        
        sizethresholds <- quantile(dat_pos$weight, c(0.15, 0.5, 0.85, 1), na.rm = T)
        for (i in 1:length(trawlids)) {
          haul_sample<- dplyr::filter(dat_pos, trawl_id == trawlids[i])
          if(nrow(haul_sample) > 0 | var(haul_sample$weight >0)) {
            # fit kernel density to weight frequency
            smoothed_w <- KernSmooth::bkde(haul_sample$weight, range.x = c(min(dat_pos$weight), max(dat_pos$weight)), bandwidth = 2)
            # make sure smoother predicts positive or zero density
            smoothed_w$y[smoothed_w$y<0] <- 0
            # calculate proportion by biomass and by number
            p_w_byweight <- smoothed_w$y * smoothed_w$x / sum(smoothed_w$x*smoothed_w$y)
            
            p_w_byweight[p_w_byweight<0] <- 0
            #p_w_bynum[p_w_bynum<0] <- 0
            
            p1 <- sum(p_w_byweight[smoothed_w$x<=sizethresholds[1]])
            p2 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[1] & smoothed_w$x <=sizethresholds[2]])
            p3 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[2] & smoothed_w$x <=sizethresholds[3]])
            p4 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[3]])
            
            p[i,2:5] <- c(p1, p2, p3, p4)
            
          }
          else {
            indx <- which(sizethresholds>haul_sample$weight)
            p[i, min(indx)+1] <- 1
          }
        }
        # add hauls with zero catch back in
        absent = filter(dat_sub, cpue_kg_km2 == 0)
        trawlids <- unique(absent$trawl_id)
        #  absent.df <- data.frame(trawl_id = trawlids,
        #      p1 = 0,
        #     p2 = 0,
        #     p3 = 0,
        #     p4 = 0)
        all_hauls <- p
        # all_hauls <- rbind(p, absent.df)
        all_hauls$trawl_id <- as.numeric(all_hauls$trawl_id)
        dat_sub$median_weight <- median(dat_sub$weight, na.rm=T)
        nlengths <- unique(dat_sub[,c("trawl_id","nlength", "median_weight")])
        meanweight <- unique(dat_test[,c("trawl_id","haul_weight")])
        all_hauls2 <- left_join(all_hauls, nlengths)
        all_hauls2 <- left_join(all_hauls2, meanweight)
        return(all_hauls2)
      }
    }
  }
  if(nrow(dat_sub)>0){
    if(nrow(fitted)==0 & nrow(not_fitted)==0){
      trawlids <- unique(dat_sub$trawl_id)
      absent.df <- data.frame(trawl_id = trawlids,
                              p1 = NA,
                              p2 = NA,
                              p3 = NA,
                              p4 = NA,
                              nlength=0,
                              haul_weight=NA, 
                              median_weight=NA)
      return(absent.df)
    }
  }
  if(nrow(dat_sub)>0){
    if(length(trawlids)==0){
      trawlids <- unique(dat_sub$trawl_id)
      absent.df <- data.frame(trawl_id = trawlids,
                              p1 = NA,
                              p2 = NA,
                              p3 = NA,
                              p4 = NA,
                              nlength=0,
                              haul_weight=NA)
      return(absent.df)
    }
  }
if(nrow(dat_sub)==0){
  return(warning("species not present in data"))
}
}


load_data_nwfsc <- function(spc,sci_name, dat.by.size, length) {
  dat <- readRDS("data/fish_raw/NOAA/nwfsc_catch.rds")
  names(dat) = tolower(names(dat))
  dat$common_name <- tolower(dat$common_name)
  dat$scientific_name <- tolower(dat$scientific_name)
  #Fix dogfish
  dat$common_name <- ifelse(dat$common_name=="pacific spiny dogfish", "spiny dogfish", dat$common_name)
  dat = dplyr::filter(dat, scientific_name == sci_name)
  if(is.data.frame(dat.by.size)){
    dat.by.size$trawl_id <- as.character(dat.by.size$trawl_id)
    dat <- left_join(dat, dat.by.size, by = "trawl_id")
    # remove tows where there was positive catch but no length measurements
  }
  if(!is.data.frame(dat.by.size)){
    dat$nlength <- NA
    dat$median_weight <- NA
    dat$haul_weight <- NA
    dat$p1 <- NA
    dat$p2 <- NA
    dat$p3 <- NA
    dat$p4 <- NA
  }
  # remove tows where there was positive catch but no length measurements
  if(length){
    dat <- dplyr::filter(dat, !is.na(p1))
  }
  # analyze or years and hauls with adequate oxygen and temperature data, within range of occurrence
  
  #O2 from trawl data is in ml/l 
  # just in case, remove any missing or nonsense values from sensors
  # dat <- dplyr::filter(dat, !is.na(o2), !is.na(sal), !is.na(temp), is.finite(sal))
  # dat <- calc_po2_mi(dat)
  # dat <- dplyr::filter(dat, !is.na(temp), !is.na(mi))
  
  # prepare data and models -------------------------------------------------
  dat$longitude <- dat$longitude_dd
  dat$latitude <- dat$latitude_dd
  dat$event_id <- dat$trawl_id
  dat$date <- as.POSIXct(as.Date(as.POSIXct("1970-01-01")+as.difftime(dat$date,units="days")))
  # get julian day
  dat$julian_day <- rep(NA, nrow(dat))
  for (i in 1:nrow(dat)){ 
    dat$julian_day[i] <- as.POSIXlt(dat$date[i], format = "%Y-%b-%d")$yday
  }
  dat <- dplyr::select(dat, event_id, common_name, project, vessel, tow, year, date, longitude_dd, latitude_dd, longitude, latitude, cpue_kg_km2,
                       depth_m, julian_day, nlength, median_weight, haul_weight, pass, p1, p2, p3, p4)
  
  # UTM transformation
  dat_ll = dat
  sp::coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
  sp::proj4string(dat_ll) <- sp::CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = sp::spTransform(dat_ll, 
                            sp::CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # convert back from sp object to data frame
  dat = as.data.frame(dat_utm)
  dat = dplyr::rename(dat, X = coords.x1,
                      Y = coords.x2)
  dat$scientific_name <- sci_name
  dat$depth <- dat$depth_m
  dat$depth_m <- NULL
  dat$survey <- "nwfsc"
  return(dat)
}

length_expand_afsc <- function(sci_name) {
  # load, clean, and join data
  bio2 <-readRDS("data/fish_raw/NOAA/ak_bts_goa_ebs_nbs_indivero_all_levels.RDS")
  catch2 <- readRDS("data/fish_raw/NOAA/ak_bts_goa_ebs_nbs_indivero_cpue_zerofilled.RDS")
  
  #Isolate necessary parts of full data to get specimen weights/lengths per haul
  haul <- bio2$haul
  specimen <- bio2$specimen
  species <- bio2$species
  size <- bio2$size
  
  #make lowercase
  names(haul) <- tolower(names(haul))
  names(specimen) <- tolower(names(specimen))
  names(species) <- tolower(names(species))
  names(size) <- tolower(names(bio2$size))
  
  species$species_name <- tolower(species$species_name)
  
  #Combine size and species
  lengths <- dplyr::left_join(size, species)
  lengths <- dplyr::left_join(lengths, haul)
  
  #Expand to make separate row for each measurement
  lengths <- dplyr::filter(lengths, !is.na(frequency))
  lengths2 <- tidyr::uncount(lengths, weights=frequency)
  
  #Combine
  bio <- dplyr::left_join(specimen, species)
  bio <- dplyr::left_join(bio, haul)
  
  #Combine specimen and length data
  bio3 <- dplyr::bind_rows(bio, lengths)
  
  #Convert length from mm to cm
  bio3$length_cm <- bio3$length*0.1
  
  #Convert weight from g to kg
  bio3$weight <- bio3$weight*0.001
  
  #Combine catch data with species data
  names(catch2) <- tolower(names(catch2))
  catch <- dplyr::left_join(catch2, species)
  
  # filter out species of interest from joined (catch/haul/bio) dataset
  catch_sub = dplyr::filter(catch, species_name == sci_name)
  bio_sub = dplyr::filter(bio3, species_name == sci_name)
  
  #Select only necessary columns for joining
  catch4 <- catch_sub[,c("hauljoin", "survey", "year", "depth_m", "latitude_dd_start", "longitude_dd_start", "cpue_kgkm2", "species_name", "common_name")]
  bio4 <- dplyr::filter(bio_sub[,c("hauljoin", "performance", "species_name", "common_name", "length_cm", "weight", "sex", "age")], !is.na(length_cm))
  #Combine data
  dat <-dplyr::left_join(bio4, catch4, relationship = "many-to-many")
  dat <- dplyr::mutate(dat, trawl_id=hauljoin)
  #According to the codebook https://repository.library.noaa.gov/view/noaa/50147, 0 means Good performance, and the other numbers are for "Satisfactory, and then a "but"..."; negative numbers are Unsatisfactory
  #Dataset already includes only Good and Satisfactory hauls, Unsatisfactory are removed
  
  #If years=T in the function, this will subset data to 1999 onward to remove possibly funky data prior to 1999
  years <- F
  if(years){
    dat_sub = dplyr::filter(dat, year>1999)
  }
  
  if(!years){
    dat_sub = dat
  }
  
  #Add column counting number of length observations per catch
  dat_sub$is_length <- ifelse(!is.na(dat_sub$length_cm), 1,0) 
  dat_sub <- group_by(dat_sub, trawl_id) %>% mutate(nlength=sum(is_length)) %>% ungroup()
  dat_sub$trawl_id <- as.numeric(dat_sub$trawl_id)
  
  # fit length-weight regression by year to predict fish weights that have lengths only.
  # note a rank-deficiency warning may indicate there is insufficient data for some year/sex combinations (likely for unsexed group)
  
  if(nrow(dat_sub)>0){
    fitted = dat_sub
    # #If region=T in function, this will do the length-weight regression for each broad geographic region (EBS, NBS, and GOA)
    # if(region){
    #   # dplyr::select(trawl_id,year,
    #   #               subsample_wt_kg, total_catch_wt_kg, area_swept_ha_der, cpue_kg_km2,
    #   #               individual_tracking_id, sex, length_cm, weight) %>%
    #   
    # #Create one set of data with only year/survey combinations with at least some weight data:
    #   fitted <-  filter(fitted, !is.na(length_cm)) %>%
    #     group_by(year,survey) %>%
    #     mutate(sum = sum(weight, na.rm=T)) %>%
    #     filter(sum>0) %>%
    #     ungroup()
    #   
    #   #Remove erroneous zero weight observations because cause error in lm()
    #   fitted$weight <- ifelse(fitted$weight==0, NA, fitted$weight)
    #   
    #  #And one for those without any weight data:
    #   not_fitted <-  filter(dat_sub, !is.na(length_cm)) %>%
    #     group_by(year,survey) %>%
    #     mutate(sum = sum(weight, na.rm=T)) %>%
    #     filter(sum==0) %>%
    #     ungroup()
    #   
    # #Fit regression model for years with data
    #   
    #   fitted <-
    #   group_nest(fitted, year,survey) %>%
    #   mutate(
    #     model = purrr::map(data, ~ lm(log(weight) ~ log(length_cm), data = .x)),
    #     tidied = purrr::map(model, broom::tidy),
    #     augmented = purrr::map(model, broom::augment),
    #     predictions = purrr::map2(data, model, modelr::add_predictions)
    #   )
    # }
    
    #If region=F in function, this will do length-weight regression for all of Alaska combined
    # dplyr::select(trawl_id,year,
    #               subsample_wt_kg, total_catch_wt_kg, area_swept_ha_der, cpue_kg_km2,
    #               individual_tracking_id, sex, length_cm, weight) %>%
    
    #Create one set of data with only year/survey combinations with at least some weight data:
    fitted <-  filter(fitted, !is.na(length_cm)) %>%
      group_by(year) %>%
      mutate(sum = sum(weight, na.rm=T)) %>%
      filter(sum>0) %>%
      ungroup()
    #And one for those without any weight data:
    not_fitted <-  filter(dat_sub, !is.na(length_cm)) %>%
      group_by(year,survey) %>%
      mutate(sum = sum(weight, na.rm=T)) %>%
      filter(sum==0) %>%
      ungroup()
    
    if(nrow(fitted)>0){
      #Fit regression model for years with data
      #Remove weird zero weights
      fitted$weight <- ifelse(fitted$weight==0, NA, fitted$weight)
      fitted <-
        group_nest(fitted, year) %>%
        mutate(
          model = purrr::map(data, ~ lm(log(weight) ~ log(length_cm), data = .x)),
          tidied = purrr::map(model, broom::tidy),
          augmented = purrr::map(model, broom::augment),
          predictions = purrr::map2(data, model, modelr::add_predictions)
        )
      # replace missing weights with predicted weights for years with data
      dat_pos2 = fitted %>%
        tidyr::unnest(predictions) %>%
        dplyr::select(-data, -model, -tidied, -augmented) %>%
        try(dplyr::mutate(weight = ifelse(is.na(weight), exp(pred), weight)))
      
      #combine back with data for years without data
      dat_pos <- bind_rows(dat_pos2, not_fitted)
    }
    if(nrow(fitted)==0){
      dat_pos <- not_fitted
    }
    ##Checked that this works by checking number of rows--looks like it passes!
    
    #If there is no weight data available to get weight-length empirical interpolation, use fishbase to fill in
    # find length-weight relationship parameters for one species
    pars <- rfishbase::length_weight(
      #convert scientific name to first letter capitalized for fishbase
      species_list = stringr::str_to_sentence(sci_name, locale = "en"))
    #Get mean
    a <- mean(pars$a)
    b <- mean(pars$b)
    #For longspine thornyhead, which is not in database for some reason, got these values from fishbase direct page
    if(sci_name=="sebastolobus altivelis"){
      a <- 0.00912
      b=3.09
    }
    #Calculate weight (and convert from g to kg)
    dat_pos <- dplyr::mutate(dat_pos, weight = ifelse(is.na(weight), ((a*length_cm^b)*0.001), weight))
    #convert cm-g units
    #haul level weight
    dat_test <- group_by(dat_pos, trawl_id) %>% mutate(haul_weight=mean(weight)) %>% ungroup()
    
    #make column of trawl_id, which is called hauljoin originally in the AFSC data
    trawlids <- unique(dat_pos$trawl_id)
    if(length(trawlids!=0)){
      p <- data.frame(trawl_id = trawlids,
                      p1 = 0,
                      p2 = 0,
                      p3 = 0,
                      p4 = 0)
      
      sizethresholds <- quantile(dat_pos$weight, c(0.15, 0.5, 0.85, 1), na.rm = T)
      for (i in 1:length(trawlids)) {
        haul_sample<- dplyr::filter(dat_pos, trawl_id == trawlids[i])
        if(nrow(haul_sample) > 0 | var(haul_sample$weight >0)) {
          # fit kernel density to weight frequency
          smoothed_w <- KernSmooth::bkde(haul_sample$weight, range.x = c(min(dat_pos$weight), max(dat_pos$weight)), bandwidth = 2)
          # make sure smoother predicts positive or zero density
          smoothed_w$y[smoothed_w$y<0] <- 0
          # calculate proportion by biomass and by number
          p_w_byweight <- smoothed_w$y * smoothed_w$x / sum(smoothed_w$x*smoothed_w$y)
          
          
          p_w_byweight[p_w_byweight<0] <- 0
          #p_w_bynum[p_w_bynum<0] <- 0
          
          p1 <- sum(p_w_byweight[smoothed_w$x<=sizethresholds[1]])
          p2 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[1] & smoothed_w$x <=sizethresholds[2]])
          p3 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[2] & smoothed_w$x <=sizethresholds[3]])
          p4 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[3]])
          
          p[i,2:5] <- c(p1, p2, p3, p4)
          
        }
        else {
          indx <- which(sizethresholds>haul_sample$weight)
          p[i, min(indx)+1] <- 1
        }
      }
      
      # add hauls with zero catch back in
      absent = filter(dat_sub, cpue_kgkm2 == 0)
      if(nrow(absent)>0){
        trawlids <- unique(absent$trawl_id)
        absent.df <- data.frame(trawl_id = trawlids,
                                p1 = 0,
                                p2 = 0,
                                p3 = 0,
                                p4 = 0)
        
        all_hauls <- rbind(p, absent.df)
      }
      if(nrow(absent)==0){
        all_hauls <- p
      }
      all_hauls$trawl_id <- as.numeric(all_hauls$trawl_id)
      dat_sub$median_weight <- median(dat_sub$weight, na.rm=T)
      nlengths <- unique(dat_sub[,c("trawl_id","nlength", "median_weight")])
      meanweight <- unique(dat_test[,c("trawl_id","haul_weight")])
      all_hauls2 <- left_join(all_hauls, nlengths)
      all_hauls2 <- left_join(all_hauls2, meanweight)
      return(all_hauls2)
    }
  }
  if(nrow(dat_sub)>0){
    if(nrow(fitted)==0 & nrow(not_fitted)==0){
      trawlids <- unique(dat_sub$trawl_id)
      absent.df <- data.frame(trawl_id = trawlids,
                              p1 = NA,
                              p2 = NA,
                              p3 = NA,
                              p4 = NA,
                              nlength=0,
                              haul_weight=NA, 
                              median_weight=NA)
      return(absent.df)
    }
  }
  #If there are no hauls at all with any length measurements, do this instead (because caused an error in the kernel density function otherwise)
  if(nrow(dat_sub)>0){
    trawlids <- unique(dat_sub$trawl_id)
    if(length(trawlids)==0){
      absent.df <- data.frame(trawl_id = trawlids,
                              p1 = NA,
                              p2 = NA,
                              p3 = NA,
                              p4 = NA,
                              nlength=0,
                              haul_weight=NA)
      return(absent.df)
    }
  }
if(nrow(dat_sub)==0){
  
  return(warning("species not present in data"))
}
  }


load_data_afsc <- function(sci_name, spc, dat.by.size, length=T) {
  bio2 <-readRDS("data/fish_raw/NOAA/ak_bts_goa_ebs_nbs_indivero_all_levels.RDS")
  dat <-readRDS("data/fish_raw/NOAA/ak_bts_goa_ebs_nbs_indivero_cpue_zerofilled.RDS")
  species <- bio2$species
  names(dat) = tolower(names(dat))
  names(species) =tolower(names(species))
  species$species_name <- tolower(species$species_name)
  species$common_name <- tolower(species$common_name)
  dat <- dplyr::left_join(dat, species)
  dat.by.size$trawl_id <- as.character(dat.by.size$trawl_id)
  dat$trawl_id <-as.character(dat$hauljoin)
  dat = dplyr::filter(dat, species_name ==sci_name)
  dat <- dplyr::left_join(dat, dat.by.size, by = "trawl_id")
  # remove tows where there was positive catch but no length measurements
  if(length){
    dat <- dplyr::filter(dat, !is.na(p1))
  }
  # analyze or years and hauls with adequate oxygen and temperature data, within range of occurrence
  
  # get julian day
  dat$julian_day <- rep(NA, nrow(dat))
  haul <- bio2$haul
  names(haul) <- tolower(names(haul))
  haul <- haul[,c("hauljoin", "start_time")]
  dat <- left_join(dat, haul)
  for (i in 1:nrow(dat)) dat$julian_day[i] <- as.POSIXlt(dat$start_time[i], format = "%Y-%b-%d")$yday
  
  #O2 from trawl data is in ml/l 
  # just in case, remove any missing or nonsense values from sensors
  # dat <- dplyr::filter(dat, !is.na(o2), !is.na(sal), !is.na(temp), is.finite(sal))
  # dat <- calc_po2_mi(dat)
  # dat <- dplyr::filter(dat, !is.na(temp), !is.na(mi))
  
  # prepare data and models -------------------------------------------------
  dat$longitude_dd <- dat$longitude_dd_start
  dat$latitude_dd <- dat$latitude_dd_start
  dat$longitude <- dat$longitude_dd
  dat$latitude <- dat$latitude_dd
  dat$scientific_name <- dat$species_name
  dat$cpue_kg_km2 <- dat$cpue_kgkm2
  dat$project <- dat$survey
  dat$event_id <- dat$trawl_id
  dat$date <- as.POSIXct(as.Date(dat$start_time, format = "%Y-%b-%d"))
  dat <- dplyr::select(dat, event_id, common_name, scientific_name, project, survey, year, date, bottom_temperature_c, longitude_dd, latitude_dd, longitude, latitude, cpue_kg_km2,
                       depth_m, julian_day, nlength, median_weight, haul_weight, p1, p2, p3, p4)
  
  
  # UTM transformation
  dat_ll = dat
  sp::coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
  sp::proj4string(dat_ll) <- sp::CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = sp::spTransform(dat_ll, 
                            sp::CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # convert back from sp object to data frame
  dat = as.data.frame(dat_utm)
  dat = dplyr::rename(dat, X = coords.x1,
                      Y = coords.x2)
  dat$depth <- dat$depth_m
  dat$depth_m <- NULL
  return(dat)
}

length_expand_bc <- function(sci_name, spc) {
  # load, clean, and join data
  itis <- readRDS("data/fish_raw/BC/species-table.rds")
  haul <- readRDS("data/fish_raw/BC/pbs-haul.rds")
  catch <- readRDS("data/fish_raw/BC/pbs-catch.rds")
  bio2 <- readRDS("data/fish_raw/BC/pbs-bio-samples.rds")
  
  catch$common_name <- catch$species_common_name
  catch$species_science_name <- NULL
  catch$species_common_name <- NULL
  
  #Rename dogfish
  catch$common_name <- ifelse(catch$common_name=="north pacific spiny dogfish", "spiny dogfish", catch$common_name)
  
  #Merge the official BC bio data and the official BC haul data to get metadata (from Sean, or here: https://open.canada.ca/data/en/dataset/86af7918-c2ab-4f1a-ba83-94c9cebb0e6c)
  bio <- dplyr::full_join(bio2, haul, by="event_id", relationship="many-to-many")
  
  #Combine with species data
  bio <- dplyr::full_join(bio, itis, by="itis", relationship="many-to-many")
  catch <- dplyr::full_join(catch, itis, relationship="many-to-many")

  #Put bio data in the same format as the NOAA bio data
  #Convert g to kg
  bio$weight <- bio$weight*0.001
  
  #rename columnns
  bio$length_cm <- bio$length
  
  #Clean catch data
  names(catch) = tolower(names(catch))
  
  #Merge with ITIS information to get scientific name
  catch <- left_join(catch,itis)
  
  #haul$sampling_end_hhmmss = as.numeric(haul$sampling_end_hhmmss)
  #haul$sampling_start_hhmmss = as.numeric(haul$sampling_start_hhmmss)
  
  #Combine catch data with haul data
  dat <- dplyr::left_join(catch, haul, relationship = "many-to-many")
  
  #Rename missing species
  if(sci_name=="sebastes aleutianus"){
    bio$scientific_name <- ifelse(str_detect(bio$scientific_name, "sebastes aleutianus"), "sebastes aleutianus", bio$scientific_name)
    bio$common_name <- ifelse(str_detect(bio$common_name, "rougheye"), "rougheye rockfish", bio$common_name)
  }
  #Combine bio/haul data with catch data
  dat <- dplyr::left_join(dat, filter(bio[,c("event_id", "age","length_cm", "weight", "scientific_name", "common_name")], !is.na(length_cm)), relationship = "many-to-many")
  
  # filter out species of interest from joined (catch/haul/bio) dataset
  dat_sub = dplyr::filter(dat, scientific_name==sci_name)

  dat_sub$event_id <- as.numeric(dat_sub$event_id)
  trawlids <- unique(dat_sub$event_id)
  
  if(nrow(dat_sub)>0) {
    #Add column counting number of length observations per catch
    dat_sub$is_length <- ifelse(!is.na(dat_sub$length_cm), 1,0) 
    dat_sub <- group_by(dat_sub, event_id) %>% mutate(nlength=sum(is_length)) %>% ungroup()
    
    # fit length-weight regression by year to predict fish weights that have lengths only.
    # note a rank-deficiency warning may indicate there is insufficient data for some year/sex combinations (likely for unsexed group)
    
    #Create one set of data with only year/survey combinations with at least some weight data:
    fitted <-  filter(dat_sub, !is.na(length_cm)) %>%
      group_by(year) %>%
      mutate(sum = sum(weight, na.rm=T)) %>%
      filter(sum>0) %>%
      ungroup()
    fitted$weight <- ifelse(fitted$weight==0, NA, fitted$weight)
    #And one for those without any weight data:
    not_fitted <-  filter(dat_sub, !is.na(length_cm)) %>%
      group_by(year) %>%
      mutate(sum = sum(weight, na.rm=T)) %>%
      filter(sum==0) %>%
      ungroup()
    
    #Fit regression model for years with data
    fitted <-
      group_nest(fitted, year) %>%
      mutate(
        model = purrr::map(data, ~ lm(log(weight) ~ log(length_cm), data = .x)),
        tidied = purrr::map(model, broom::tidy),
        augmented = purrr::map(model, broom::augment),
        predictions = purrr::map2(data, model, modelr::add_predictions)
      )
    
    # replace missing weights with predicted weights
    dat_pos2 = fitted %>%
      tidyr::unnest(predictions) %>%
      dplyr::select(-data, -model, -tidied, -augmented) %>%
      try(dplyr::mutate(weight = ifelse(is.na(weight), exp(pred), weight)))
    
    #combine back with data for years without data
    dat_pos <- bind_rows(dat_pos2, not_fitted)
    
    ##Checked that this works by checking number of rows--looks like it passes!
    
    #If there is no weight data available to get weight-length empirical interpolation, use fishbase to fill in
    # find length-weight relationship parameters for one species
    pars <- rfishbase::length_weight(
      #convert scientific name to first letter capitalized for fishbase
      species_list = stringr::str_to_sentence(sci_name, locale = "en"))
    #Get mean
    a <- mean(pars$a)
    b <- mean(pars$b)
    #For longspine thornyhead, which is not in database for some reason, got these values from fishbase direct page
    if(sci_name=="sebastolobus altivelis"){
      a <- 0.00912
      b=3.09
    }
    #Calculate weight (and convert from g to kg)
    dat_pos <- dplyr::mutate(dat_pos, weight = ifelse(is.na(weight), ((a*length_cm^b)*0.001), weight))
    dat_test <- group_by(dat_pos, event_id) %>% mutate(haul_weight=mean(weight)) %>% ungroup()
    
    trawlids <- unique(dat_pos$event_id)
    if(length(trawlids!=0)){
      p <- data.frame(event_id = trawlids,
                      p1 = 0,
                      p2 = 0,
                      p3 = 0,
                      p4 = 0)
      
      sizethresholds <- quantile(dat_pos$weight, c(0.15, 0.5, 0.85, 1), na.rm = T)
      for (i in 1:length(trawlids)) {
        haul_sample<- dplyr::filter(dat_pos, event_id == trawlids[i])
        if(nrow(haul_sample) > 0 | var(haul_sample$weight >0)) {
          # fit kernel density to weight frequency
          smoothed_w <- KernSmooth::bkde(haul_sample$weight, range.x = c(min(dat_pos$weight), max(dat_pos$weight)), bandwidth = 2)
          # make sure smoother predicts positive or zero density
          smoothed_w$y[smoothed_w$y<0] <- 0
          # calculate proportion by biomass and by number
          p_w_byweight <- smoothed_w$y * smoothed_w$x / sum(smoothed_w$x*smoothed_w$y)
          
          p_w_byweight[p_w_byweight<0] <- 0
          #p_w_bynum[p_w_bynum<0] <- 0
          
          p1 <- sum(p_w_byweight[smoothed_w$x<=sizethresholds[1]])
          p2 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[1] & smoothed_w$x <=sizethresholds[2]])
          p3 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[2] & smoothed_w$x <=sizethresholds[3]])
          p4 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[3]])
          
          p[i,2:5] <- c(p1, p2, p3, p4)
          
        }
        else {
          indx <- which(sizethresholds>haul_sample$weight)
          p[i, min(indx)+1] <- 1
        }
      }
      
      # add hauls with zero catch back in
      absent = filter(dat_sub, catch_weight == 0)
      trawlids <- unique(absent$event_id)
      absent.df <- data.frame(event_id = trawlids,
                              p1 = 0,
                              p2 = 0,
                              p3 = 0,
                              p4 = 0)
      
      all_hauls <- rbind(p, absent.df)
      all_hauls$event_id <- as.numeric(all_hauls$event_id)
      dat_sub$median_weight <- median(dat_sub$weight, na.rm=T)
      nlengths <- unique(dat_sub[,c("event_id","nlength", "median_weight")])
      meanweight <- unique(dat_test[,c("event_id","haul_weight")])
      all_hauls2 <- left_join(all_hauls, nlengths)
      all_hauls2 <- left_join(all_hauls2, meanweight)
      return(all_hauls2)
    }
  }
  if(nrow(dat_sub)>0){
    if(nrow(fitted)==0 & nrow(not_fitted)==0){
      trawlids <- unique(dat_sub$event_id)
      absent.df <- data.frame(event_id = trawlids,
                              p1 = NA,
                              p2 = NA,
                              p3 = NA,
                              p4 = NA,
                              nlength=0,
                              haul_weight=NA, 
                              median_weight=NA)
      return(absent.df)
    }
  }
  if(nrow(dat_sub)>0){
    if(length(trawlids)==0){
      absent.df <- data.frame(event_id = trawlids,
                              p1 = NA,
                              p2 =NA,
                              p3 = NA,
                              p4 = NA,
                              nlength=0,
                              haul_weight=NA)
      return(absent.df)
    }
  }
  if(nrow(dat_sub)==0){
    
    return(warning("species not present in data"))
  }
  }

load_data_bc <- function(sci_name,dat.by.size, length=T, spc) {
  catch <- readRDS("data/fish_raw/BC/pbs-catch.rds")
  haul <- readRDS("data/fish_raw/BC/pbs-haul.rds")
  itis <- readRDS("data/fish_raw/BC/species-table.rds")
  #Rename dogfish
  catch$common_name <- ifelse(catch$species_common_name=="north pacific spiny dogfish", "spiny dogfish", catch$species_common_name)
  catch$species_common_name <- NULL
  catch <- left_join(catch,itis)
  
  #haul$sampling_end_hhmmss = as.numeric(haul$sampling_end_hhmmss)
  #haul$sampling_start_hhmmss = as.numeric(haul$sampling_start_hhmmss)
  
  #Combine catch data with haul data
  dat <- dplyr::left_join(catch, haul, relationship = "many-to-many")
  
  # dat.by.size$event_id <- as.character(dat.by.size$event_id)
  dat = dplyr::filter(dat, scientific_name == sci_name)
  dat <- left_join(dat, dat.by.size, by = "event_id")
  # remove tows where there was positive catch but no length measurements
  if(length){
    dat <- dplyr::filter(dat, !is.na(p1))
  }
  # analyze or years and hauls with adequate oxygen and temperature data, within range of occurrence
  
  # get julian day
  dat$julian_day <- rep(NA, nrow(dat))
  for (i in 1:nrow(dat)) dat$julian_day[i] <- as.POSIXlt(dat$date[i], format = "%Y-%b-%d")$yday
  
  
  #O2 from trawl data is in ml/l 
  # just in case, remove any missing or nonsense values from sensors
  # dat <- dplyr::filter(dat, !is.na(o2), !is.na(sal), !is.na(temp), is.finite(sal))
  # dat <- calc_po2_mi(dat)
  # dat <- dplyr::filter(dat, !is.na(temp), !is.na(mi))
  
  # prepare data and models -------------------------------------------------
  dat$cpue_kg_km2 <- dat$catch_weight
  dat$longitude_dd <- dat$lon_start
  dat$latitude_dd <- dat$lat_start
  dat$longitude <- dat$lon_start
  dat$latitude <- dat$lat_start
  dat$event_id <- as.character(dat$event_id)
  dat$year <- substr(dat$date, start=1, stop=4)
  dat$year <- as.integer(dat$year)
  dat$date <- as.POSIXct(as.Date(dat$date, format = "%Y-%b-%d"))
  dat$project <- dat$survey_name
  dat$salinity_psu <- dat$salinity_PSU
  dat$salinity_PSU <- NULL
  dat <- dplyr::select(dat, event_id, scientific_name, project, year, date, longitude_dd, latitude_dd, longitude, latitude, cpue_kg_km2,
                       depth_m, julian_day, nlength,median_weight, haul_weight, pass, p1, p2, p3, p4, temperature_C, do_mlpL, salinity_psu)
  dat <- filter(dat, !is.na(latitude_dd))
  
  
  # UTM transformation
  dat_ll = dat
  sp::coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
  sp::proj4string(dat_ll) <- sp::CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = sp::spTransform(dat_ll, 
                            sp::CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # convert back from sp object to data frame
  dat = as.data.frame(dat_utm)
  dat = dplyr::rename(dat, X = coords.x1,
                      Y = coords.x2)
  dat$common_name <- spc
  dat$depth <- dat$depth_m
  dat$depth_m <- NULL
  dat$survey <- "dfo"
  
  #convert oxygen mg/L to umol_kg
  SA = gsw_SA_from_SP(dat$salinity_psu,dat$depth,dat$longitude,dat$latitude) #absolute salinity for pot T calc
  pt = gsw_pt_from_t(SA,dat$temperature_C,dat$depth) #potential temp at a particular depth
  CT = gsw_CT_from_t(SA,dat$temperature_C,dat$depth) #conservative temp
  dat$sigma0_kgm3 = gsw_sigma0(SA,CT)
  dat$O2_umolkg = dat$do_mlpL*44660/(dat$sigma0_kgm3+1000) 
  
  return(dat)
}

IPHC <- function (catch, adjustment, spc, sci_name) {
  #make lowercase
  colnames(catch) <- tolower(colnames(catch))
  colnames(adjustment) <- tolower(colnames(adjustment))
  
  if(spc!="pacific halibut"){
    catch$scientific_name <- tolower(catch$'scientific name')
    catch$common_name<- paste(spc)
    catch <- subset(catch, catch$scientific_name==sci_name)
    catch$cpue_count <- (catch$'number observed'/catch$hooksobserved)*catch$hooksretrieved
    catch$cpue_weight <- NA
    #Get CTD data and sample metadata to join
    dat <-  read_excel("data/fish_raw/IPHC/Set and Pacific halibut data.xlsx")
    colnames(dat) <- tolower(colnames(dat))
    dat <- dat[,c("year","stlkey","setno","date", 'midlat fished','midlon fished', "avgdepth (fm)", "temp c", "salinity psu", "oxygen_ml")]
    #Join
    catch <- left_join(catch, dat)
    ##Add zeros
    #All sets
    sets <- unique(dat)
    #Which are not in the catch data (positive catches)
    missing <- anti_join(sets, catch)
    #Add species columns back
    missing$common_name <- paste(spc)
    missing$scientific_name <- paste(sci_name)
    #Add weight and count
    missing$cpue_weight <- 0
    missing$cpue_count <- 0
    #Add to dataset
    catch <- bind_rows(catch, missing)
     }
  if(spc=="pacific halibut"){
  #Calculate CPUE as U32 and O32 count and weight divided by effective skates hauled
  catch$cpue_O32_count <- catch$`o32 pacific halibut count`/catch$`effective skates hauled`
  catch$cpue_U32_count <- catch$`u32 pacific halibut count`/catch$`effective skates hauled`
  catch$cpue_O32_weight <- catch$`o32 pacific halibut weight`/catch$`effective skates hauled`
  catch$cpue_U32_weight <- catch$`u32 pacific halibut weight`/catch$`effective skates hauled`
  #Sum weight classes together
  catch$cpue_count <- catch$cpue_O32_count+catch$cpue_U32_count
  catch$cpue_weight <- catch$cpue_O32_weight+catch$cpue_U32_weight
  catch$common_name <- "pacific halibut"
  catch$scientific_name <-"hippoglossus stenolepis"
  }
  #join
  adjustment$stlkey <- as.character(adjustment$stlkey)
  data <- left_join(catch, adjustment, by="stlkey")
  data$h.adj <- as.numeric(data$h.adj)
  
  #Calculate ajustment factor
  data$cpue_count <- data$cpue_count * data$h.adj
  data$cpue_weight <- data$cpue_weight * data$h.adj
  #Convert biomass to kg
  data$cpue_weight <- data$cpue_weight*0.453592
  #Extract columns of interest
  dat <- data[,c("year.x", "date.x", "common_name", "scientific_name", "midlat fished", "midlon fished", "avgdepth (fm)","temp c", "salinity psu", "oxygen_ml", "cpue_weight", "cpue_count")]
  #Re-name columns to match the NOAA and BC data
  colnames(dat) <- c("year", "date", "common_name", "scientific_name", "latitude", "longitude", "depth", "temperature_C", "salinity_psu", "do_mlpL", "cpue_weight", "cpue_count")

  #Add columns to match the NOAA and BC data
  dat$survey <- "iphc"
  
  ##Add region
  # load regional polygons
  regions.hull <- readRDS("data/processed_data/regions_hull.rds")
  #make dataframe an sf object
  dat <- drop_na(dat, latitude, longitude)
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
  #Date and month in right format
  dat$month <- case_when(grepl("May",dat$date) ~5,
                         grepl("Jun",dat$date)  ~6,
                         grepl("Jul",dat$date)  ~7,
                         grepl("Aug",dat$date)  ~8,
                         grepl("Sep",dat$date)  ~9,
                         grepl("Oct",dat$date)  ~10)
  dat$day <- as.numeric(substr(dat$date, 1,2))
  dat$date <-  as.POSIXct(as.Date(with(dat,paste(year,month,day,sep="-")),"%Y-%m-%d"))
  dat$doy <- as.POSIXlt(dat$date, format = "%Y-%b-%d")$yday
  dat$year <- as.numeric(dat$year)
  dat$day <- NULL
  #convert fathoms to m
  dat$depth <- dat$depth*1.8288
  
  #convert oxygen mg/L to umol_kg
  SA = gsw_SA_from_SP(dat$salinity_psu,dat$depth,dat$longitude,dat$latitude) #absolute salinity for pot T calc
  pt = gsw_pt_from_t(SA,dat$temperature_C,dat$depth) #potential temp at a particular depth
  CT = gsw_CT_from_t(SA,dat$temperature_C,dat$depth) #conservative temp
  dat$sigma0_kgm3 = gsw_sigma0(SA,CT)
  dat$O2_umolkg = dat$do_mlpL*44660/(dat$sigma0_kgm3+1000) 
  
  #Convert coordinates
  dat <- subset(dat, !is.na(latitude))
  dat <- dat %>%
    st_as_sf(coords=c('longitude','latitude'),crs=4326,remove = F) %>%  
    st_transform(crs = "+proj=utm +zone=10 +datum=WGS84 +units=km") %>% 
    mutate(X=st_coordinates(.)[,1],Y=st_coordinates(.)[,2]) 
 
  return(dat)
}

combine_all <- function(type){
  #BC
  #BC raw data
  bio_qc <- read.csv("data/fish_raw/BC/QCS_biology.csv")
  bio_vi <- read.csv("data/fish_raw/BC/WCVI_biology.csv")
  bio_hs <- read.csv("data/fish_raw/BC/HS_biology.csv")
  bio_hg <- read.csv("data/fish_raw/BC/WCHG_biology.csv")
  itis  <- readRDS("data/fish_raw/BC/itis_bc.rds")
  
  haul <- readRDS("data/fish_raw/BC/pbs-haul.rds")
  catch <- readRDS("data/fish_raw/BC/pbs-catch.rds")
  
  haul_qc <- read.csv("data/fish_raw/BC/QCS_effort.csv")
  haul_vi <- read.csv("data/fish_raw/BC/WCVI_effort.csv")
  haul_hs <- read.csv("data/fish_raw/BC/HS_effort.csv")
  haul_hg <- read.csv("data/fish_raw/BC/WCHG_effort.csv")
  
  #Combine BC bio and haul data to combine together
  bio2 <- rbind(bio_hg, bio_hs, bio_qc, bio_vi)
  haul2 <- rbind(haul_hg, haul_hs, haul_qc, haul_vi)
  
  #Merge the official BC bio data and the official BC haul data to get metadata (from here: https://open.canada.ca/data/en/dataset/86af7918-c2ab-4f1a-ba83-94c9cebb0e6c)
  bio2$Set.number <- bio2$Tow.number
  bio <- dplyr::left_join(bio2, haul2, relationship="many-to-many")
  
  names(bio) = tolower(names(bio))
  names(haul) = tolower(names(haul))
  
  #Put bio data in the same format as the NOAA bio data
  bio$scientific_name <- tolower(bio$scientific.name)
  bio$common_name <- tolower(bio$english.common.name)
  bio$weight <- bio$weight..g.*0.001
  bio$length_cm <- bio$fork.length..mm.*0.1
  bio$length_cm <- ifelse(is.na(bio$length_cm), (bio$total.length..mm.*0.1), bio$length_cm)
  
  #Make column names consistent so can join to surveyjoin haul data to get event_id
  haul$year <- as.character(substr(haul$date, start=1, stop=4))
  bio$year <- as.character(bio$survey.year)
  haul$set.date <- substr(haul$date, start=1, stop=10)
  bio$lat_start <- bio$start.latitude
  bio$lat_end <- bio$end.latitude
  bio$lon_start <- bio$start.longitude
  bio$lon_end <- bio$end.longitude
  
  #Combine BC bio data and haul metadata with surveyjoin metadata
  bio3 <- dplyr::left_join(bio[,c("age", "sex", "length_cm", "weight", "scientific_name", "common_name", "lat_start", "lon_start", "lat_end", "lon_end", "set.date")], haul, relationship="many-to-many")
  bio <- bio3
  
  #Clean catch data
  names(catch) = tolower(names(catch))
  
  catch <- dplyr::left_join(catch,itis)
  
  #haul$sampling_end_hhmmss = as.numeric(haul$sampling_end_hhmmss)
  #haul$sampling_start_hhmmss = as.numeric(haul$sampling_start_hhmmss)
  
  #Combine catch data with haul data
  dat_bc <- dplyr::left_join(catch, haul, relationship = "many-to-many")
  #Combine bio/haul data with catch data
  #Only positive catches
  dat_bc <- subset(dat_bc, catch_weight>0)
  dat_bc$cpue_kg_km2 <- dat_bc$catch_weight
  #Biomass
  if(type=="hauls"){
    counts_bc <- aggregate(cpue_kg_km2~scientific_name, dat_bc, FUN=length)
  }
  if(type=="biomass"){
    counts_bc <- aggregate(cpue_kg_km2~scientific_name, dat_bc, FUN=sum)
  }
  
  #Add common name
  counts_bc <- unique(dplyr::left_join(counts_bc, bio[,c("scientific_name", "common_name")]))
  counts_bc <- dplyr::mutate(counts_bc, bc=cpue_kg_km2)
  counts_bc$cpue_kg_km2 <- NULL
  rm(list=setdiff(ls(), "counts_bc"))
  
  #NWFSC
  # load, clean, and join data
  bio <- readRDS("data/fish_raw/NOAA/nwfsc_bio.rds")
  load("data/fish_raw/NOAA/nwfsc_haul.rda")
  haul <- nwfsc_haul
  catch <- readRDS("data/fish_raw/NOAA/nwfsc_catch.rds")
  names(catch) = tolower(names(catch))
  names(bio) = tolower(names(bio))
  names(haul) = tolower(names(haul))
  bio$scientific_name <- tolower(bio$scientific_name)
  bio$common_name <- tolower(bio$common_name)
  catch$common_name <- tolower(catch$common_name)
  
  bio$trawl_id = as.character(bio$trawl_id)
  haul$trawl_id = as.character(haul$event_id)
  catch$trawl_id=as.character(catch$trawl_id)
  haul$year <- as.character(substr(haul$date, start=1, stop=4))
  bio$year <- as.character(bio$year)
  catch$date <- NULL
  bio$date <- NULL
  bio$year <- NULL
  
  #haul$sampling_end_hhmmss = as.numeric(haul$sampling_end_hhmmss)
  #haul$sampling_start_hhmmss = as.numeric(haul$sampling_start_hhmmss)
  
  #Combine data
  dat = dplyr::left_join(catch[,c("trawl_id","common_name", "subsample_count","area_swept_ha","longitude_dd", "latitude_dd",
                                  "subsample_wt_kg","total_catch_numbers","total_catch_wt_kg","cpue_kg_km2")], haul, relationship = "many-to-many") %>%
    dplyr::left_join(filter(bio[,c("trawl_id", "scientific_name", "common_name", "weight", "ageing_lab", "oto_id", "length_cm", "width_cm", "sex", "age")], !is.na(length_cm)), relationship = "many-to-many") %>%
    filter(performance == "Satisfactory")
  
  #Only positive catches
  dat_nw <- subset(dat, cpue_kg_km2>0)
  #Biomass
  if(type=="hauls"){
    counts_nw <- aggregate(cpue_kg_km2~scientific_name, dat_nw, FUN=length)
  }
  if(type=="biomass"){
    counts_nw <- aggregate(cpue_kg_km2~scientific_name, dat_nw, FUN=sum)
  }
  counts_nw <- unique(left_join(counts_nw, bio[,c("scientific_name", "common_name")]))
  counts_nw <- dplyr::mutate(counts_nw, nw=cpue_kg_km2)
  counts_nw$cpue_kg_km2 <- NULL
  
  rm(list=setdiff(ls(), c("counts_bc", "counts_nw")))
  
  ##Alaska
  bio2 <-readRDS("data/fish_raw/NOAA/ak_bts_goa_ebs_nbs_all_levels.RDS")
  catch2 <- readRDS("data/fish_raw/NOAA/ak_bts_goa_ebs_nbs_cpue_zerofilled.RDS")
  
  #Isolate necessary parts of full data to get specimen weights/lengths per haul
  haul <- bio2$haul
  specimen <- bio2$specimen
  species <- bio2$species
  size <- bio2$size
  
  #make lowercase
  names(haul) <- tolower(names(haul))
  names(specimen) <- tolower(names(specimen))
  names(species) <- tolower(names(species))
  names(size) <- tolower(names(bio2$size))
  
  species$species_name <- tolower(species$species_name)
  
  #Combine catch data with species data
  names(catch2) <- tolower(names(catch2))
  catch <- dplyr::left_join(catch2, species)
  
  #Select only necessary columns for joining
  catch4 <- catch[,c("hauljoin", "survey", "year", "depth_m", "latitude_dd_start", "longitude_dd_start", "cpue_kgkm2", "species_name", "common_name")]
  dat <- dplyr::filter(catch4, year>1998)
  #bio4 <- dplyr::filter(bio4, year>1998)
  
  #Combine data
  dat <- dplyr::mutate(dat, trawl_id=hauljoin)
  #According to the codebook https://repository.library.noaa.gov/view/noaa/50147, 0 means Good performance, and the other numbers are for "Satisfactory, and then a "but"..."; negative numbers are Unsatisfactory
  #Dataset already includes only Good and Satisfactory hauls, Unsatisfactory are removed
  
  #Only positive catches
  dat$cpue_kg_km2 <- dat$cpue_kgkm2
  dat_ak <- subset(dat, cpue_kg_km2>0)
  dat_ak$scientific_name <- dat_ak$species_name
  #Biomass
  if(type=="hauls"){
    counts_ak <- aggregate(cpue_kg_km2~scientific_name, dat_ak, FUN=length)
  }
  if(type=="biomass"){
    counts_ak <- aggregate(cpue_kg_km2~scientific_name, dat_ak, FUN=sum)
  }
  counts_ak <- dplyr::mutate(counts_ak, ak=cpue_kg_km2)
  counts_ak$cpue_kg_km2 <- NULL
  counts_ak <- dplyr::left_join(counts_ak, species[,c("common_name", "species_name")], by=c("scientific_name"="species_name"))
  counts_ak$common_name <- tolower(counts_ak$common_name)
  
  #counts_ak <- unique(left_join(counts_nw, bio4[,c("scientific_name", "common_name")]))
  
  #Combine
  counts <- full_join(counts_nw, counts_bc)
  counts <- full_join(counts, counts_ak)
  return(counts)
}

###LOAD DATA ###
prepare_data <- function(spc,sci_name,mi,iphc, file_name, taxa){
  dat.by.size <- try(length_expand_bc(sci_name))
  gc()
  if(is.data.frame(dat.by.size)){
    dat3 <- try(load_data_bc(sci_name = sci_name, spc=spc, dat.by.size = dat.by.size, length=F))
    dat3 <- try(unique(dat3))
  }
  gc()
  rm(dat.by.size)
  
  dat.by.size <- try(length_expand_nwfsc(spc=spc, sci_name=sci_name))
  gc()
  if(is.data.frame(dat.by.size)){
    dat2 <- try(load_data_nwfsc(spc= spc, sci_name=sci_name, dat.by.size = dat.by.size, length=F))
    dat2 <- try(unique(dat2))
  }
  
  if(!is.data.frame(dat.by.size)){
    catch <- readRDS("data/fish_raw/NOAA/nwfsc_catch.rds")
    names(catch) <- tolower(names(catch))
    catch$common_name <- tolower(catch$common_name)
    catch$scientific_name <- tolower(catch$scientific_name)
    catch <- dplyr::filter(catch, scientific_name == sci_name)
    if(length(catch)>0){
      dat.by.size <- matrix(data=NA, nrow=1,ncol=1)
      dat2 <- try(load_data_nwfsc(spc= spc, sci_name=sci_name, dat.by.size = dat.by.size, length=F))
      dat2 <- try(unique(dat2))
    }
  }
  
  gc()
  rm(dat.by.size)
  
  dat.by.size <- try(length_expand_afsc(sci_name))
  gc()
  if(is.data.frame(dat.by.size)){
    dat5 <- try(load_data_afsc(sci_name = sci_name, spc=spc, dat.by.size = dat.by.size, length=F))
    dat5 <- try(unique(dat5))
  }
  
  gc()
  
  #All regions present
  if(exists("dat3") & exists("dat2") & exists("dat5")){
    dat4 <- bind_rows(dat3, dat2, dat5)
  }
  #Only BC
  if(exists("dat3") & !exists("dat2") & !exists("dat5")){
    dat4 <- dat3
  }
  #Only AK
  if(!exists("dat3") & !exists("dat2") & exists("dat5")){
    dat4 <- dat5
  }
  
  #Only NWFSC
  if(!exists("dat3") & exists("dat2") & !exists("dat5")){
    dat4 <- dat2
    
  }
  #BC & NW
  if(exists("dat3") & exists("dat2") & !exists("dat5")){
    dat4 <- bind_rows(dat3, dat2)
  }
  #BC & AK
  if(exists("dat3") & !exists("dat2")& exists("dat5")){
    dat4 <- bind_rows(dat3, dat5)
  }
  #AK & NW
  if(!exists("dat3") & exists("dat2")& exists("dat5")){
    dat4 <- bind_rows(dat2, dat5)
  }
  gc()

  #  adjustment <- read_excel("~/Dropbox/choke species/code/choke-species-data/data/fish_raw/IPHC/iphc-2023-fiss-hadj-20231031.xlsx")
   # dat_IPHC <- IPHC(catch, adjustment)
  #  dat4 <- bind_rows(dat4, dat_IPHC)
#  }
  #Make a duplicate copy object
  dat <- dat4
  
  ###Add in situ data####
    #Isolate just NWFSC and EBS/NBS/GOA data for joining, because  dfo are already in the dataset
    insitu <- readRDS("data/processed_data/o2/insitu_combined.rds")
    insitu <- filter(insitu, survey!="iphc")
    insitu <- filter(insitu, survey!="dfo")
    insitu$event_id <- as.character(insitu$event_id)
    #Northwest
    dat6 <- filter(dat4, survey=="nwfsc")
    #Just columns of interest
    dat6 <- left_join(dat6[,c("event_id", "date", "year", "survey", "latitude", "longitude", "depth", "X", "Y", "cpue_kg_km2", "julian_day", "nlength", "median_weight", "haul_weight", "pass", "p1", "p2", "p3", "p4", "common_name", "scientific_name", "depth", "vessel", "tow", "bottom_temperature_c")], insitu, by=c("date", "latitude", "longitude"))
    #Edit columns
    dat6$event_id <- dat6$event_id.x
    dat6$event_id.x <- NULL
    dat6$event_id.y <- NULL
    dat6$year <- dat6$year.x
    dat6$year.y <- NULL
    dat6$year.x <- NULL
    dat6$depth <- dat6$depth.x
    dat6$depth.x <- NULL
    dat6$depth.y <- NULL
    dat6$X <- dat6$X.x
    dat6$Y <- dat6$Y.x
    dat6$X.x <- NULL
    dat6$X.y <- NULL
    dat6$Y.x <- NULL
    dat6$Y.y <- NULL
    dat6$depth.1 <- NULL
    dat6$month <- NULL
    dat6$doy <- NULL
    dat6$survey <- dat6$survey.x
    dat6$survey.x <- NULL
    dat6$survey.y <- NULL
    dat6$time <- NULL
    
    dat7 <- dat6
    
    #GOA, EBS & NBS 
    dat6 <- filter(dat4, survey=="EBS"|survey=="NBS"|survey=="GOA")
    dat6 <- left_join(dat6[,c("event_id", "date", "year", "survey", "latitude", "longitude", "depth", "X", "Y", "cpue_kg_km2", "julian_day", "nlength", "median_weight", "haul_weight", "pass", "p1", "p2", "p3", "p4", "common_name", "scientific_name", "depth", "vessel", "tow", "bottom_temperature_c")], insitu[,c("event_id", "temperature_C", "do_mlpL", "salinity_psu", "sigma0_kgm3", "O2_umolkg")], by=c("event_id"))
    
    #Edit columns
    dat6$depth.1 <- NULL

    #Recombine back with DFO data
    if(exists("dat3")){
      dat <- bind_rows(dat6, dat7,dat3)
    }
    
    if(!exists("dat3")){
      dat <- bind_rows(dat6, dat7)
    }
  
    #Some other columns
    dat$doy <- dat$julian_day
    dat$julian_day <- NULL
    dat$cpue_kg_km2_sub <- dat$cpue_kg_km2 * (dat$p2+dat$p3)
    dat$cpue_kg_km2_sub <- ifelse(dat$cpue_kg_km2==0, 0, dat$cpue_kg_km2_sub)
    dat$p1 <- ifelse(dat$cpue_kg_km2==0, 0, dat$p1)
    dat$p2 <- ifelse(dat$cpue_kg_km2==0, 0, dat$p2)
    dat$p3 <- ifelse(dat$cpue_kg_km2==0, 0, dat$p3)
    dat$p4 <- ifelse(dat$cpue_kg_km2==0, 0, dat$p4)
    dat$survey <- tolower(dat$survey)
    dat$survey <- ifelse(dat$survey=="nwfsc", "wcbts", dat$survey)
    dat$region<- case_when(dat$survey=="nbs"~"nbs",
                           dat$survey=="wcbts"~"cc",
                           dat$survey=="ebs"~"ebs",
                           dat$survey=="goa"~"goa",
                           dat$survey=="dfo"~"bc")
     
    #Combine with IPHC if needed
    if(iphc){
      if(spc=="pacific halibut"){
        catch <-  read_excel("data/fish_raw/IPHC/Set and Pacific halibut data.xlsx")
      }
      if(spc!="pacific halibut"){
        catch <-  read_excel("data/fish_raw/IPHC/Non_Pacific_halibut_data.xlsx")
      }
      adjustment <- read_excel("data/fish_raw/IPHC/iphc-2024-fiss-hadj.xlsx")
      dat_IPHC <- IPHC(catch, adjustment, spc, sci_name)
      dat <- bind_rows(dat, dat_IPHC)
    }
    
    #Some other columns
    dat$log_depth_scaled <- scale(log(dat$depth))
    dat$log_depth_scaled2 <- with(dat, log_depth_scaled ^ 2)
    dat$month <- month(dat$date)    #Some other columns
    
    if(mi==T){
    ##Calculate metabolic index
    #Calculate pO2 from umol kg
    dat$po2 <- calc_po2_sat(salinity=dat$salinity_psu, temp=dat$temperature_C, depth=dat$depth, oxygen=dat$O2_umolkg, lat=dat$latitude, long=dat$longitude, umol_m3=T, ml_L=F)
    
    ##Calculate inverse temp
    kelvin = 273.15
    boltz = 0.000086173324
    tref <- 12
    dat$invtemp <- (1 / boltz)  * ( 1 / (dat$temperature_C + 273.15) - 1 / (tref + 273.15))
    
    ###Calculate Metabolic index 
    ##Species parameters from Tim's paper
    #Gadidae (lower, median, upper)
    if(taxa=="gadidae") {
    Eo <- c(-0.03064428,0.1883451,0.40414)
    V <- c(2.852883, 1.716267, 8.796332)
    n <- c(-0.1416557, 0.1883451, 0.05982743)
    }
    if(taxa=="perciformes"){
    Eo <- c(0.0267307, 0.3102310, 0.5859761)
    V <- c(1.734103, 1.59323, 8.836855)
    n <- c(-0.1720744, 0.3102311, 0.03412396)
    }
    if(taxa=="elasmobranchii"){
    #ie. Squalidae
    Eo <- c(-0.07659791, 0.2661157, 0.6013813)
    V <- c(0.7086717, 1.502317, 9.654109)
    n <-c(-0.2180308, 0.2661157, 0.05064831)
    }
    if(taxa=="teleostei"){
    #Orders Scorpaenidae, Soleida, Pleuronectida)
    Eo <- c(-0.009823544, 0.2554435,0.5315111)
    V <- c(1.639916,1.565438,8.815259)
    n <-c(-0.1562307,0.2554435,0.03221273)
    }
    
    #Abbreviated equation
    #dat$mi1 <- dat$po2*exp(Eo1* dat$invtemp)
    #dat$mi2 <- dat$po2*exp(Eo2* dat$invtemp)
    #dat$mi2 <- dat$po2*exp(Eo2* dat$invtemp)
    
    ##If using full equation and weight
    #Could link to Tim's estimates
    #MI_pars <- readRDS("metabolic_index-main/MI_pars.rds")
    #Use average weight across hauls of the species if haul is missing the weight
    dat$haul_weight <- ifelse(is.na(dat$haul_weight), mean(dat$median_weight, na.rm=T), dat$haul_weight)
    #If no weights at all for the species, use an average mid-size weight from the sablefish paper (# this works for the adult class (0.5 - 6 kg).  for the large adult, adjust)
    dat$haul_weight <- ifelse(is.na(dat$haul_weight), 3, dat$haul_weight)
    
    #Calculate metabolic index with each Eo
    #Calculate Ao from V above (using median)
    Ao <- 1/exp(V[2])
    n <- n[2]
    
    dat$mi1 <- calc_mi(Eo[1], Ao, dat$haul_weight,  n, dat$po2, dat$invtemp)
    dat$mi2 <- calc_mi(Eo[2], Ao, dat$haul_weight, n, dat$po2, dat$invtemp)
    dat$mi3 <- calc_mi(Eo[3], Ao, dat$haul_weight, n, dat$po2, dat$invtemp)
    
    #For average body size from each haul
    #dat$mi5 <- calc_mi(Eo1, Ao, dat$haul_weight,  n, dat$po2, dat$invtemp)
    
    dat$po2_s <- (scale(dat$po2))
    dat$mi1_s <-(scale(dat$mi1))
    dat$mi2_s <-(scale(dat$mi2))
    dat$mi3_s <-(scale(dat$mi3))
    #dat$mi5_s <- scale(dat$mi5)
    }
    
    #Reorder columns
    if(iphc==F){
      dat <- relocate(dat, scientific_name, common_name, survey, year, date, doy, month, depth, longitude, latitude, cpue_kg_km2, cpue_kg_km2_sub, salinity_psu, temperature_C, sigma0_kgm3,do_mlpL, O2_umolkg, invtemp, po2, mi1,mi2,mi3,log_depth_scaled, log_depth_scaled2, X, Y, p1,p2,p3,p4,median_weight, haul_weight, nlength, pass, vessel, tow, bottom_temperature_c)
    }
    if(iphc==T){
      dat <- relocate(dat, scientific_name, common_name, survey, year, date, doy, month, depth, longitude, latitude, cpue_kg_km2, cpue_kg_km2_sub,cpue_weight, cpue_count, salinity_psu, temperature_C, sigma0_kgm3, do_mlpL, O2_umolkg, invtemp, po2,mi1,mi2,mi3,log_depth_scaled, log_depth_scaled2, X, Y, p1,p2,p3,p4,median_weight, haul_weight, nlength, pass, vessel, tow, bottom_temperature_c)
      
    }
#Remove duplicates
    dat <- unique(dat)
    #Remove geometry
    dat$geometry <- NULL
    saveRDS(dat, file = paste("data/processed_data/fish/dat_", file_name, ".rds", sep=""))
  try(return(dat))
}

load_fish <- function(species, test_region){
  file <- list.files("data/processed_data", pattern=paste(species))
  file <- paste("data/processed_data/", file, sep="")
  dat <- readRDS(file)
  
  #Filter region
  if(test_region=="cc"){
    test_survey <- "nwfsc"
  }
  if(test_region=="bc"){
    test_survey <- "dfo"
  }
  dat <- as.data.frame(dat)
  #Filter to region
  dat <- filter(dat, survey==test_survey)
  return(dat)
}

###MODEL FITTING ###
model_comparison <- function(species, test_region, gloryswd, basewd, GOBH, remove_outlier, filter_size, new_breakpt, filter_years, combine_pred){
  #Pull species data
  #Set seed
  # set.seed(5)
  
  #Fish data
  file <- list.files("data/processed_data", pattern=paste(species))
  file <- paste("data/processed_data/", file, sep="")[1]
  dat <- readRDS(file)
  dat <- as.data.frame(dat)
  dat$geometry <- NULL
  
  #Filter region
  if(test_region=="cc"){
    test_survey <- "nwfsc"
    #Filter to region
    dat <- filter(dat, survey==test_survey)
  }
  if(test_region=="bc"){
    test_survey <- "dfo"
    #Filter to region
    dat <- filter(dat, survey==test_survey)
  }
  dat <- as.data.frame(dat)
  
  #Merge temp data
  dat$temperature_C <- ifelse((is.na(dat$temperature_C)& !is.na(dat$bottom_temperature_c)), dat$bottom_temperature_c, dat$temperature_C)
  dat$temp <- dat$temperature_C
  #Remove NAs
  dat <- dat  %>%
    drop_na(depth,year, temp, X, Y)
  if(remove_outlier==T){
    #Remove outliers = catch > 10 sd above the mean
    dat$cpue_s <- scale(dat$cpue_kg_km2)
    dat <- dplyr::filter(dat, cpue_s <=20)
  }
  
  if(filter_size==T){
    dat$cpue_kg_km2 <- dat$cpue_kg_km2_sub
  }
  
  #Format columns
  dat$depth_ln <- log(dat$depth)
  dat$log_depth_scaled <- scale(log(dat$depth))
  dat$log_depth_scaled2 <- dat$log_depth_scaled^2
  dat$temp_s <- scale(dat$temp)
  
  ###Syoptic in situ only###
  dat_syn <- dat %>%
    drop_na(O2_umolkg)
  dat_syn$o2 <- dat_syn$O2_umolkg
  dat_syn$o2_s <- scale(dat_syn$o2)
  dat_syn$log_depth_scaled <- scale(log(dat_syn$depth))
  dat_syn$log_depth_scaled2 <- dat_syn$log_depth_scaled^2
  dat_syn$o2_s <- scale(dat_syn$o2)
  
  ###Fit model and predict to data###
  #Load oxygen data
  dat_o2 <- as.data.frame(readRDS("data/processed_data/all_o2_dat_filtered.rds"))
  
  #Filter test region
  dat_o2 <- filter(dat_o2, region==test_region)
  
  #Remove any rows with missing dat_o2a
  dat_o2 <- dat_o2 %>%
    drop_na(depth, o2, temp,doy, X, Y, year)
  
  #Remove oxygen outliers
  dat_o2 <- filter(dat_o2, o2<1500)
  
  #Set minimum sigma
  minsigma0 <- 24
  dat_o2$sigma0[dat_o2$sigma0 <= minsigma0] <- minsigma0
  
  # remove older (earlier than 2000) dat_o2a
  dat_o2 <- dplyr::filter(dat_o2, year >=2000)
  
  #Log depth
  dat_o2$depth_ln <- log(dat_o2$depth)
  
  ##Fit in situ model ##
  #Extra years
  missing_years <- setdiff(unique(dat$year), unique(dat_o2$year))
  train_years <- seq(from=min(dat_o2$year), to=max(dat_o2$year), by=1)
  train_years2 <- setdiff(train_years, unique(dat_o2$year))
  extra_years <- append(missing_years, train_years2)
  extra_years <- extra_years[!is.na(extra_years)]
  if(length(extra_years)==0) {extra_years = NULL}  
  
  spde <- make_mesh(data = dat_o2,
                    xy_cols = c("X", "Y"),
                    cutoff = 45)
  
  #Scale
  dat_o2$o2 <- dat_o2$o2/100
  dat_o2 <- as.data.frame(dat_o2)
  
  print("fitting integrated model")
  m <- sdmTMB(
    formula = o2  ~ 1+s(depth_ln) + s(doy)+s(temp),
    mesh = spde,
    dat = dat_o2,
    family = gaussian(),
    time = "year",
    spatial = "on",
    spatiotemporal  = "ar1",
    extra_time=c(extra_years))
  
  #Add predictions of O2
  dat_pred <- predict(m, newdata = dat)
  dat_pred$o2 <- dat_pred$est*100
  #Replace with zero
  dat_pred$o2 <- ifelse(dat_pred$o2<0, 0, dat_pred$o2)
  dat_pred$o2_s <- scale(dat_pred$o2)
  
  if(GOBH!=T){
    dat_gobh <- dat
    dat_gobh$o2 <- dat_gobh$o2_gobh
    dat_gobh$o2 <- ifelse(dat_gobh$o2<0, 0, dat_gobh$o2)
    dat_gobh$o2_s <- scale(dat_gobh$o2)
    dat_gobh <- dat_gobh %>%
      drop_na(o2_s, temp, X, Y, depth, year)
  }
  
  if(GOBH==T){
    ##GOBH predictions
    #Set do threshold level for GLORYS data
    do_threshold <- 0
    
    #Filter days of GLORYS data to every 10th day?
    filter_time <- T
    
    #Filter GLORYS data to bottom depth? ()
    filter_depth <- F
    
    print("fitting and predicting GOBH model")
    #Apply
    dat_gobh <- gobh_interpolate(dat, "cc", filter_depth, filter_time, gloryswd, basewd, do_threshold)
    #Combine
    dat_gobh <- bind_rows(dat_gobh)
    
    #Get columns set up
    dat_gobh$o2 <- dat_gobh$est*100
    dat_gobh$o2 <- ifelse(dat_gobh$o2<0, 0, dat_pred$o2)
    dat_gobh$o2_s <- scale(dat_gobh$o2)
    
  }
  
  #Filter years to the same years as synoptic?
  if(filter_years==T){
    dat_pred <- dat_pred %>% drop_na(O2_umolkg)
    dat_gobh <- dat_gobh %>% drop_na(O2_umolkg)
  }
  #Merge predicted and in situ data--use in situ when available, and predicted only when not
  if(combine_pred==T){
    dat_pred$o2 <- ifelse(is.na(dat_pred$O2_umolkg), dat_pred$o2, dat_pred$O2_umolkg)
  }
  
  #List the dataframes to test
  dats <- list(dat_syn, dat_pred, dat_gobh)
  #List data names
  dat_names <- c("Synoptic", "Predicted", "GOBH")
  model_names <- c("null", "temp","o2", "temp+o2", "temp+o2+temp*02", "breakpt(o2)", "breakpt(o2)+temp")
  
  #Create summary matrix
  aic <- as.data.frame(matrix(nrow=length(model_names), ncol=length(dat_names)))
  par_summary <- list()
  preds <- list()
  models <- list()
  
 # if(new_breakpt==T){
  #  remove.packages("sdmTMB")
  #  remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE, ref="newlogistic")
  #  library(sdmTMB)
#  }
  
  print("fitting SDMs")
  for (i in 1:length(dats)){
    ###Test alternative fish SDMs with the three different oxygen data
    dat <- as.data.frame(dats[i])
    paste(dat_names[i])
    extra_time <- setdiff(min(dat$year):max(dat$year),unique(dat$year))
    
    spde <- make_mesh(data = dat,
                      xy_cols = c("X", "Y"),
                      cutoff = 45)
    print("fitting SDM m1")
    m1 <- try(sdmTMB(
      formula = cpue_kg_km2  ~ 1+as.factor(year)+log_depth_scaled+log_depth_scaled2,
      mesh = spde,
      data = dat,
      family = tweedie(link="log"),
      time=NULL,
      spatial = "on",
      spatiotemporal  = "off",
      anisotropy=TRUE
    ))
    print("fitting SDM m2")
    m2 <- try(sdmTMB(
      formula = cpue_kg_km2  ~ 1+as.factor(year)+log_depth_scaled+log_depth_scaled2+temp_s,
      mesh = spde,
      data = dat,
      family = tweedie(link="log"),
      time=NULL,
      spatial = "on",
      spatiotemporal  = "off",
      anisotropy=TRUE
    ))
    print("fitting SDM m3")
    m3 <- try(sdmTMB(
      formula = cpue_kg_km2  ~ 1+as.factor(year)+log_depth_scaled+log_depth_scaled2+o2_s,
      mesh = spde,
      data = dat,
      family = tweedie(link="log"),
      time=NULL,
      spatial = "on",
      spatiotemporal  = "off",
      anisotropy=TRUE
    ))
    print("fitting SDM m4")
    m4 <- try(sdmTMB(
      formula = cpue_kg_km2  ~ 1+as.factor(year)+log_depth_scaled+log_depth_scaled2+temp_s+o2_s,
      mesh = spde,
      data = dat,
      family = tweedie(link="log"),
      time=NULL,
      spatial = "on",
      spatiotemporal  = "off",
      anisotropy=TRUE
    ))
    print("fitting SDM m5")
    m5 <- try(sdmTMB(
      formula = cpue_kg_km2  ~ 1+as.factor(year)+log_depth_scaled+log_depth_scaled2+temp_s+o2_s+temp_s*o2_s,
      mesh = spde,
      data = dat,
      family = tweedie(link="log"),
      time=NULL,
      spatial = "on",
      spatiotemporal  = "off",
      anisotropy=TRUE
    ))
    print("fitting SDM m6")
    if(test_region=="bc"){
      start <- matrix(0, nrow = 2, ncol = 1)
      start[1,1] <- 1
      start[2,1] <- -1
      lower <- matrix(0, nrow = 2, ncol = 1)
      lower[1,1] <- 0
      lower[2,1] <- -Inf
      upper <- matrix(0, nrow = 2, ncol = 1)
      upper[1,1] <- Inf
      upper[2,1] <- -Inf
      m6 <- try(sdmTMB(
        formula = cpue_kg_km2  ~ 1+as.factor(year)+log_depth_scaled+log_depth_scaled2+breakpt(o2_s),
        mesh = spde,
        data = dat,
        family = tweedie(link="log"),
        time=NULL,
        spatial = "on",
        spatiotemporal  = "off",
        anisotropy=TRUE,
        control = sdmTMBcontrol(start = list(b_threshold = start),
                                #lower = list(b_threshold = lower),
                                #  upper=list(b_threshold = upper)
        )))
    }
    if(test_region=="cc"){
      m6 <- try(sdmTMB(
        formula = cpue_kg_km2  ~ 1+as.factor(year)+log_depth_scaled+log_depth_scaled2+breakpt(o2_s),
        mesh = spde,
        data = dat,
        family = tweedie(link="log"),
        time=NULL,
        spatial = "on",
        spatiotemporal  = "off",
        anisotropy=TRUE))
    }
    print("fitting SDM m7")
    if(test_region=="bc"){
      m7 <- try(sdmTMB(
        formula = cpue_kg_km2  ~ 1+as.factor(year)+log_depth_scaled+log_depth_scaled2+breakpt(o2_s)+temp_s,
        mesh = spde,
        data = dat,
        family = tweedie(link="log"),
        time=NULL,
        spatial = "on",
        spatiotemporal  = "off",
        anisotropy=TRUE,
        control = sdmTMBcontrol(start = list(b_threshold = start),
                                #lower = list(b_threshold = lower),
                                #  upper=list(b_threshold = upper)
        )))
    }
    if(test_region=="cc"){
      m7 <- try(sdmTMB(
        formula = cpue_kg_km2  ~ 1+as.factor(year)+log_depth_scaled+log_depth_scaled2+breakpt(o2_s)+temp_s,
        mesh = spde,
        data = dat,
        family = tweedie(link="log"),
        time=NULL,
        spatial = "on",
        spatiotemporal  = "off",
        anisotropy=TRUE))
    }
    model_list <- list(m1, m2,m3, m4,m5,m6,m7)
    names(model_list) <- model_names
    
    #Other model ideas: smoother on depth, bivariate spline on depth and o2
    
    ##Add AIC to summary table
    temp <- matrix(nrow=7, ncol=2)
    temp[1,1] <- try(AIC(m1))
    temp[2,1] <- try(AIC(m2))
    temp[3,1] <- try(AIC(m3))
    temp[4,1] <- try(AIC(m4))
    temp[5,1] <- try(AIC(m5))
    temp[6,1] <- try(AIC(m6))
    temp[7,1] <- try(AIC(m7))
    
    temp[,2] <- try(abs(min(temp[,1], na.rm=T)-(temp[,1])))
    
    aic[,i] <- temp[,2]
    
    ##Extract each parameter table
    pars <- mapply(extract_pars, model_list, model_names, SIMPLIFY=F)
    pars <- try(bind_rows(pars))
    pars$data <- paste(dat_names[i])
    par_summary[[i]] <- pars
    
    ##Predict fish density
    temp <- list()
    for(j in 1:length(model_list)){
      dat <- dats[[i]]
      if(test_survey=="dfo"){pred_year <- 2019L}
      if(test_survey=="nwfsc"){pred_year <- 2012L}
      nd_po2 <- data.frame(o2_s = seq(min(dat$o2_s), max(dat$o2_s), length.out = 300), 
                           temp_s = 0,
                           log_depth_scaled = 0,
                           log_depth_scaled2 = 0,
                           year = pred_year)
      nd_po2 <- convert_class(nd_po2)
      p1 <- predict(model_list[[j]], newdata = nd_po2, se_fit = TRUE, re_form = NA)
      p1$model <- paste(model_names[j])
      p1$data <- paste(dat_names[i])
      temp[[j]] <- p1
    }
    names(temp) <- model_names
    preds[[i]] <- temp
    models[[i]] <- model_list
  }
  par_summary <- as.data.frame(bind_rows(par_summary))
  colnames(aic) <- dat_names
  aic <- as.data.frame(aic)
  aic$data <- model_names
  names(models) <- dat_names
  names(preds) <- dat_names
  names(dats) <- dat_names
  output <- list(dats,models, aic, par_summary, preds)
  names(output) <- c("data", "models", "aic_table", "parameter_estimates", "predictions")
  
  return(output)
}

###MODEL EVALUATION ###
extract_pars <- function(model, model_name){
  model_name <- paste(model_name)               
  if(is.list(model)){
    p1 <- try(tidy(model, effects="fixed", conf.int=T))
    p2 <- try(tidy(model, effects="ran_pars", conf.int=T))
    p <- try(as.data.frame(bind_rows(p1,p2)))
    p$model <- try(paste(model_name))
  }
  return(p)
  if(!is.list(model)){
    p <- data.frame(ncol=7, nrow=1)
    p$model <- paste(model_name)
  }
}

plot_marginal2 <- function(output, model){
  pred <- output[["predictions"]][["Synoptic"]][[print(model)]]
  pred2 <- output[["predictions"]][["Predicted"]][[print(model)]]
  pred3 <- output[["predictions"]][["GOBH"]][[print(model)]]
  dat <- output[["data"]][["Synoptic"]]
  pred$center <- attr(dat$o2_s, "scaled:center")
  pred$scale <- attr(dat$o2_s, "scaled:scale")
  dat <- output[["data"]][["Predicted"]]
  pred2$center <- attr(dat$o2_s, "scaled:center")
  pred2$scale <- attr(dat$o2_s, "scaled:scale")
  dat <- output[["data"]][["GOBH"]]
  pred3$center <- attr(dat$o2_s, "scaled:center")
  pred3$scale <- attr(dat$o2_s, "scaled:scale")
  preds <- bind_rows(pred, pred2, pred3)
  preds$unscaled <-  preds$o2_s* preds$scale+preds$center
  plot <- ggplot(preds, aes(unscaled, y=exp(est)))+
    geom_line(aes(colour=data))+
    geom_ribbon(aes(ymin = (exp(est) - exp(est_se)), ymax = (exp(est) + exp(est_se)), fill=data), alpha=0.4)+
    #scale_y_continuous(limits = c(100, 400), expand = expansion(mult = c(0, 0.0))) +
    xlim(0,200)+
    labs(x = bquote('Oxygen'~mu~"mol"~(kg^-1)), y = bquote('Population Density'~(kg~km^-2)))+
    theme_minimal()+
    theme(text=element_text(size=15))
  return(plot)
}

plot_marginal3 <- function(output, data_type, model){
  toplot <- list()
  
  for (i in 1:length(outputs)){
    output <- outputs[[i]]
    pred <- output[["predictions"]][["Synoptic"]][[print(model)]]
    pred2 <- output[["predictions"]][["Predicted"]][[print(model)]]
    pred3 <- output[["predictions"]][["GOBH"]][[print(model)]]
    dat <- output[["data"]][["Synoptic"]]
    pred$center <- attr(dat$o2_s, "scaled:center")
    pred$scale <- attr(dat$o2_s, "scaled:scale")
    dat <- output[["data"]][["Predicted"]]
    pred2$center <- attr(dat$o2_s, "scaled:center")
    pred2$scale <- attr(dat$o2_s, "scaled:scale")
    dat <- output[["data"]][["GOBH"]]
    pred3$center <- attr(dat$o2_s, "scaled:center")
    pred3$scale <- attr(dat$o2_s, "scaled:scale")
    preds <- bind_rows(pred, pred2, pred3)
    preds$unscaled <-  preds$o2_s* preds$scale+preds$center
    preds$region <- paste(region[i])
    preds$species <- paste(species[i])
    toplot[[i]] <- preds
  }
  preds <- bind_rows(toplot)
  preds$species <- factor(preds$species, levels=c("Sablefish", "Dover sole")) 
  preds$region <- factor(preds$region, levels=c("California Current", "British Columbia")) 
  
  plot <- ggplot(preds, aes(unscaled, y=exp(est)))+
    facet_wrap(region~species, scales="free_y")+
    geom_line(aes(colour=data))+
    geom_ribbon(aes(ymin = (exp(est) - exp(est_se)), ymax = (exp(est) + exp(est_se)), fill=data), alpha=0.4)+
    # scale_y_continuous(limits=c(0,50), expand = expansion(mult = c(0, 0.0))) +
    xlim(0,200)+
    labs(x = bquote('Oxygen'~mu~"mol"~(kg^-1)), y = bquote('Population Density'~(kg~km^-2)))+
    theme_minimal()+
    theme(legend.position="top")+
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank())+
    theme(text=element_text(size=15))
  return(plot)
}


plot_marginal <- function(output, data_type, model){
  pred <- output[["predictions"]][[print(data_type)]][[print(model)]]
  dat <- output[["data"]][[print(data_type)]]
  center <- attr(dat$o2_s, "scaled:center")
  scale <- attr(dat$o2_s, "scaled:scale")
  m <- output[["models"]][[print(data_type)]][[print(model)]]
  plot <- ggplot(pred, aes(back.convert(pred$o2_s, center, scale), y=exp(est)))+
    geom_line()+
    geom_ribbon(pred, mapping=aes(ymin = (exp(est) - exp(est_se)), ymax = (exp(est) + exp(est_se))), alpha=0.4)+
    #scale_y_continuous(limits = c(100, 400), expand = expansion(mult = c(0, 0.0))) +
    xlim(0,200)+
    labs(x = bquote('Oxygen'~mu~"mol"~(kg^-1)), y = bquote('Population Density'~(kg~km^-2)))+
    theme_minimal()+
    theme(text=element_text(size=15))
  return(plot)
}

par_table <- function(output, species, region){
  pars <- output[["parameter_estimates"]]
  pars <- filter(pars, (model=="breakpt(o2)")& (term=="o2_s-slope"|term=="o2_s-breakpt"))
  dat1 <- output[["data"]][["Synoptic"]]
  dat2 <- output[["data"]][["Predicted"]]
  dat3 <- output[["data"]][["GOBH"]]
  pars$center <- case_when(pars$data=="Synoptic"~attr(dat1$o2_s, "scaled:center"),
                           pars$data=="Predicted"~attr(dat2$o2_s, "scaled:center"),
                           pars$data=="GOBH"~attr(dat3$o2_s, "scaled:center"))
  pars$scale <- case_when(pars$data=="Synoptic"~attr(dat1$o2_s, "scaled:scale"),
                          pars$data=="Predicted"~attr(dat2$o2_s, "scaled:scale"),
                          pars$data=="GOBH"~attr(dat3$o2_s, "scaled:scale"))
  #Un-scale
  pars$estimate <-  round(pars$estimate* pars$scale+pars$center, 2)
  pars$std.error <-  round(pars$std.error* pars$scale+pars$center, 2)
  #Columns needed
  pars <- select(pars, data, term,estimate, std.error)
  #Combine into one
  pars$est <- paste(pars$estimate,"+-",pars$std.error)
  #Columns needed
  pars <- select(pars, data, term, est)
  #Wide
  pars <- pivot_wider(pars,names_from=data, values_from=est)
  write.csv(pars, file=paste("outputs/table_pars", species, region, sep="_"), row.names=F)
  pars <- as.data.frame(pars)
  return(pars)
}

plot_dat <- function(output){
  dat1 <- output[["data"]][[1]]
  dat1$type <- "Synoptic"
  dat2 <- output[["data"]][[2]]
  dat2$type <- "Integrated \n Prediction"
  dat3 <- output[["data"]][[3]]
  dat3$type <- "GOBH"
  all_dat <- bind_rows(dat1, dat2,dat3)
  
  ggplot(all_dat, aes(x=o2, y=-depth))+geom_point(aes(colour=log(cpue_kg_km2+1)))+
    facet_wrap("type")+                                          
    theme_minimal()+
    theme(text=element_text(size=15))+
    labs(x=bquote('Oxygen'~mu~"mol"~(kg^-1)))+
    ylab("Depth (m)")+
    scale_colour_viridis_c(
      limits = c(0, 10),
      oob = scales::squish,
      breaks = c(0, 5,10),
      name=bquote('log(CPUE'~"kg"~km^-2~")"))
}

par_table_raw <- function(output, species, region){
  pars <- output[["parameter_estimates"]]
  pars <- filter(pars, (model=="breakpt(o2)")& (term=="o2_s-slope"|term=="o2_s-breakpt"))
  dat1 <- output[["data"]][["Synoptic"]]
  dat2 <- output[["data"]][["Predicted"]]
  dat3 <- output[["data"]][["GOBH"]]
  pars$center <- case_when(pars$data=="Synoptic"~attr(dat1$o2_s, "scaled:center"),
                           pars$data=="Predicted"~attr(dat2$o2_s, "scaled:center"),
                           pars$data=="GOBH"~attr(dat3$o2_s, "scaled:center"))
  pars$scale <- case_when(pars$data=="Synoptic"~attr(dat1$o2_s, "scaled:scale"),
                          pars$data=="Predicted"~attr(dat2$o2_s, "scaled:scale"),
                          pars$data=="GOBH"~attr(dat3$o2_s, "scaled:scale"))
  #Un-scale
  pars$estimate <-  round(pars$estimate* pars$scale+pars$center, 2)
  pars$std.error <-  round(pars$std.error* pars$scale+pars$center, 2)
  #Columns needed
  # pars <- select(pars, data, term,estimate, std.error)
  #Wide
  #  pars <- pivot_wider(pars,names_from=data, values_from=est)
  # write.csv(pars, file=paste("outputs/table_pars", species, region, sep="_"), row.names=F)
  pars <- as.data.frame(pars)
  return(pars)
}
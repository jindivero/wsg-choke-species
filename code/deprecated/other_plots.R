bp_est4 <- bp.2.use
for(i in 1:nrow(bp_est4)){
  dat.2.est <- bp_est4[i,]
  this_species = dat.2.est$species
  print(this_species)
  this_dat <- dat.2.est$data
  this_model <- dat.2.est$model
  best_model <- this_model
  #Pull the data file
  this_datframe <-try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
  #Calculate mean of this_datframe$mi based on best model
  mean_mi <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") mean(this_datframe$mi1) else if(best_model=="model4"|best_model=="model10"|best_model=="model14") mean(this_datframe$mi2) else mean(this_datframe$mi3)
  sd_mi <- if(best_model=="model3"|best_model=="model9"|best_model=="model13") sd(this_datframe$mi1) else if(best_model=="model4"|best_model=="model10"|best_model=="model14") sd(this_datframe$mi2) else sd(this_datframe$mi3)
  #Pull breakpoint and slope
  #find taxa from species
  taxa.2.use <- taxa$MI_Taxa[taxa$common_name==this_species]
  #body size
  body_size <- 1
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
  test$est_o2 <- unlist(lapply(invtemp, calc_po2_crit, taxa.2.use,thresh_est,body_size, model, fancy=T))
  if(!is.null(thresh_se)){
    test$est_o2_se1 <- unlist(lapply(invtemp, calc_po2_crit, taxa.2.use,thresh_se1,body_size, model, fancy=T))
    test$est_o2_se2 <- unlist(lapply(invtemp, calc_po2_crit, taxa.2.use,thresh_se2,body_size, model, fancy=T))
  }
  # test$percentile <- (test$po2-test$est_o2)/test$est_o2_se
  test$unsuitable <- ifelse(test$po2<test$est_o2, "unsuitable", "suitable")
  if(!is.null(thresh_se)){
    test$unsuitable_low <- ifelse(test$po2<test$est_o2_se1, "unsuitable", "suitable")
    test$unsuitable_high <- ifelse(test$po2<test$est_o2_se2,"unsuitable", "suitable")
  }
  #Sum across row to get total suitable
  #test$suitable_total <- rowSums(test[,c("unsuitable", "unsuitable_low", "unsuitable_high")])
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
  scale_colour_manual(values=c("blue3", "darkorange"), labels=c("Above pO2 crit", "Below pO2 crit"))+
  guides(color = guide_legend(title="",override.aes = list(size = 7, alpha=1)))

ggsave(
  paste0("output/", output_folder, "/unsuitable_suitable_depth.png"),
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

if(output_folder=="presence_absence"){
  dats$colour <- case_when(dats$presence==1&dats$unsuitable=="suitable" ~ "True Presence",
                           dats$presence==1&dats$unsuitable=="unsuitable" ~ "False Presence",
                           dats$presence==0&dats$unsuitable=="suitable" ~ "False Absence",
                           dats$presence==0&dats$unsuitable=="unsuitable" ~ "True Absence")
  ggplot(data = dats,aes(x = temperature_C, y=depth)) +  
    geom_point(aes(colour=as.factor(colour)), size=0.25, alpha=0.1)+
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
    # scale_colour_manual(values=c("salmon1", "salmon3", "darkseagreen2", "springgreen4"))+
    guides(color = guide_legend(title="",override.aes = list(size = 7, alpha=1)))
  
  ggsave(
    paste0("output/", output_folder, "/unsuitable_suitable_depth.png"),
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
}
##Shaded by CPUE
ggplot(data = filter(dats, unsuitable==0),aes(x = temperature_C, y=depth)) +  
  geom_point(aes(colour=log(catch_weight+1)), size=0.25, alpha=0.4)+
  scale_colour_distiller(name="Above pO2 crit", type="seq",palette="Blues", direction=1)+
  new_scale_colour() +
  geom_point(data=filter(dats, unsuitable==1), mapping=aes(x = temperature_C, y=depth, colour=log(catch_weight+1)), size=0.25, alpha=0.4)+
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

##save separate for each species, colored by cpue
for(i in 1:nrow(bp_est4)){
  test <- bp_est4[i,]
  this_species <- test$species
  this_dat <- test$data
  dat.2.plot <- filter(dats, common_name==this_species & data==this_dat)
  p <- ggplot(data = filter(dat.2.plot, unsuitable==0),aes(x = temperature_C, y=depth)) +  
    geom_point(aes(colour=log(catch_weight+1)), size=0.25, alpha=0.4)+
    scale_colour_distiller(name="Above pO2 crit", type="seq",palette="Blues", direction=1)+
    new_scale_colour() +
    geom_point(data=filter(dats, unsuitable==1), mapping=aes(x = temperature_C, y=depth, colour=log(catch_weight+1)), size=0.25, alpha=0.4)+
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

##save separate for each species, with faded by CPUE
for(i in 1:nrow(bp_est4)){
  test <- bp_est4[i,]
  this_species <- test$species
  this_dat <- test$data
  dat.2.plot <- filter(dats, common_name==this_species & data==this_dat)
  p <- ggplot(data = filter(dat.2.plot, unsuitable==0),aes(x = temperature_C, y=depth)) +  
    geom_point(aes(colour=log(catch_weight+1)), size=0.25, alpha=0.4)+
    scale_colour_distiller(name="Above pO2 crit", type="seq",palette="Blues", direction=1)+
    new_scale_colour() +
    geom_point(data=filter(dats, unsuitable==1), mapping=aes(x = temperature_C, y=depth, colour=log(catch_weight+1)), size=0.25, alpha=0.4)+
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

###Map#########
map_data <- rnaturalearth::ne_countries(scale = "large",
                                        returnclass = "sf",
                                        continent = "North America")

us_coast_proj <- sf::st_transform(map_data, crs = 32610)

###Map of data available
species <- unique(dats$common_name)
dats$unsuitable <- factor(dats$unsuitable, levels=c("suitable","unsuitable"))

for(i in 1:length(species)){
  species2plot <- species[i]
  dat2plot <- filter(dats, common_name==species2plot)
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(dat2plot, mapping=aes(x=X*1000, y=Y*1000,colour=unsuitable), size=0.1)+
    xlim(min(dat2plot$X)*1000, max(dat2plot$X)*1000)+
    ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
    facet_wrap("year", ncol=5)+
    theme_minimal(base_size=12)+
    xlab("Longitude")+
    ylab("Latitude")+
    ggtitle(paste(unique(dat2plot$common_name)))+
    theme(axis.text.x=element_blank())+
    scale_colour_manual(values=c("lightblue", "orange3"), labels=c("Above pO2 crit", "Below pO2 crit"), drop=F)+
    guides(color = guide_legend(title="",override.aes = list(size = 7, alpha=1)))+
    theme(legend.position="none")
  
  
  ggsave(
    paste("output/plots/habitat_mapping/map_",species[i],".png"),
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

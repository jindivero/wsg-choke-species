plot_predict <- function(test_predict_O2, test_year, us_coast_proj) {
  # Plot ####
  data_plot <- ggplot(us_coast_proj) + geom_sf() +
    geom_point(
      data = test_predict_O2,
      aes(
        x = X * 1000,
        y = Y * 1000,
        col = do
      ),
      size = 1.0,
      alpha = 1.0
    ) +
    scale_x_continuous(breaks = c(-125, -120), limits = xlimits) +
    ylim(ylimits[1], ylimits[2]) +
    scale_colour_viridis_c(
      limits = c(0, 200),
      oob = scales::squish,
      name = bquote(O[2]),
      breaks = c(0, 100, 200)
    ) +
    labs(x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank()
      ,
      panel.grid.minor = element_blank()
      ,
      panel.border = element_blank()
      ,
      strip.background = element_blank()
      ,
      strip.text = element_blank()
    ) +
    theme(axis.line = element_line(color = "black")) +
    theme(axis.text = element_text(size = 12)) +
    theme(axis.title = element_text(size = 14)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.position = "bottom") +
    guides(colour = guide_colourbar(title.position = "top", title.hjust =
                                      0.5))
  
  predict_plot <- ggplot(us_coast_proj) + geom_sf() +
    geom_point(
      data = test_predict_O2,
      aes(
        x = X * 1000,
        y = Y * 1000,
        col = (est)
      ),
      size = 1.0,
      alpha = 1.0
    ) +
    scale_x_continuous(breaks = c(-125, -120), limits = xlimits) +
    ylim(ylimits[1], ylimits[2]) +
    scale_colour_viridis_c(
      limits = c(0, 200),
      oob = scales::squish,
      name = bquote(O[2]),
      breaks = c(0, 100, 200)
    ) +
    labs(x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank()
      ,
      panel.grid.minor = element_blank()
      ,
      panel.border = element_blank()
      ,
      strip.background = element_blank()
      ,
      strip.text = element_blank()
    ) +
    theme(axis.line = element_line(color = "black")) +
    theme(axis.text = element_text(size = 12)) +
    theme(axis.title = element_text(size = 14)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.position = "bottom") +
    guides(colour = guide_colourbar(title.position = "top", title.hjust =
                                      0.5))
  
  residual_plot <- ggplot(us_coast_proj) + geom_sf() +
    geom_point(
      data = test_predict_O2,
      aes(
        x = X * 1000,
        y = Y * 1000,
        col = residual
      ),
      size = 1.0,
      alpha = 1.0
    ) +
    scale_x_continuous(breaks = c(-125, -120), limits = xlimits) +
    ylim(ylimits[1], ylimits[2]) +
    scale_colour_distiller(palette = "RdBu", limits = c(-50, 50)) +
    #, limits = c(-40, 40), oob = scales::squish, name = bquote(O[2]), breaks = c(-40, 0, 40)) +
    labs(x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank()
      ,
      panel.grid.minor = element_blank()
      ,
      panel.border = element_blank()
      ,
      strip.background = element_blank()
      ,
      strip.text = element_blank()
    ) +
    theme(axis.line = element_line(color = "black")) +
    theme(axis.text = element_text(size = 12)) +
    theme(axis.title = element_text(size = 14)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.position = "bottom") +
    guides(colour = guide_colourbar(title.position = "top", title.hjust =
                                      0.5))
  
  ## put all plots in one ####
  grid.arrange(data_plot, predict_plot, residual_plot, ncol = 3)
  # plot residuals vs. prediction
  resid_vs_pred <- ggplot(data = test_predict_O2, aes(x = (est), y = residual, col = Y)) +
    geom_point() +
    scale_colour_viridis_c(
      limits = c(31, 50),
      oob = scales::squish,
      name = "latitude",
      breaks = c(35, 40, 45)
    ) +
    ggtitle(test_year) +
    labs(x = "Predicted", y = "Residual") +
    theme(legend.position = "none")
  pred_vs_actual <- ggplot(data = test_predict_O2, aes(x = o2, y = est, col = Y)) +
    geom_point() +
    scale_colour_viridis_c(
      limits = c(31, 50),
      oob = scales::squish,
      name = "latitude",
      breaks = c(35, 40, 45)
    ) +
    ggtitle(test_year) +
    labs(x = "Observed", y = "Predicted") +
    geom_abline(intercept = 0, slope = 1) +
    theme(legend.position = "none")
  return(grid.arrange(pred_vs_actual, resid_vs_pred, ncol = 2))
}

calc_rmse <- function(rmse_list, n){
  rmse2 <- as.data.frame(rmse_list)
  rmse2$n <- n
  rmse2 <- filter(rmse2, n>50)
  rmse2$rmse2 <- rmse2$rmse_list ^ 2
  rmse2$xminusxbarsq <- rmse2$n * rmse2$rmse2
  rmse2 <- drop_na(rmse2, xminusxbarsq)
  rmse_total<- sqrt(sum(rmse2$xminusxbarsq, na.rm=T) / sum(rmse2$n, na.rm=T))
  return(rmse_total)
}

plot_simple <- function(output, dat.2.use){
  #Separate test and training data, predictions, and models from output list
  train_data <- output[[1]]
  test_data <- output[[2]]
  preds <- output[[3]]
  models <- output[[4]]
  #Set latitude and longitude
  xlims <- c(min(dat.2.use$X)*1000, max(dat.2.use$X)*1000)
  ylims <- c(min(dat.2.use$Y)*1000, max(dat.2.use$Y)*1000)
  lats <- c(round(min(dat.2.use$latitude)),  round(max(dat.2.use$latitude)))
  lons <- c(round(min(dat.2.use$longitude)+2), round(max(dat.2.use$longitude)))
  data_map <- 
    ggplot(us_coast_proj) + geom_sf() +
    geom_point(train_data, mapping=aes(x=X*1000, y=Y*1000, col=o2),
               # data = train_data,
               #aes(
               # x = X * 1000,
               # y = Y * 1000,
               #  col = o2
               #  ),
               size = 1.0,
               alpha = 1.0
    ) +
    ylim(ylims)+
    scale_x_continuous(breaks=lons, limits=xlims)+
    scale_colour_viridis_c(
      limits = c(0, 200),
      oob = scales::squish,
      name = bquote(O[2]),
      breaks = c(0, 100, 200)
    ) +
    labs(x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank()
      ,
      panel.grid.minor = element_blank()
      ,
      panel.border = element_blank()
      ,
      strip.background = element_blank()
      ,
      strip.text = element_blank()
    ) +
    theme(axis.line = element_line(color = "black")) +
    theme(axis.text = element_text(size = 11)) +
    theme(axis.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 11)) +
    theme(legend.position = "none")
  #guides(colour = guide_colourbar(title.position = "top", title.hjust =
  #   0.5))
  
  dat_available <- ggplot(train_data, aes(x=year))+
    stat_count(aes(fill=survey))+
    theme_bw() +
    theme(
      panel.grid.major = element_blank()
      ,
      panel.grid.minor = element_blank()
      ,
      panel.border = element_blank()
      ,
      strip.background = element_blank()
      ,
      strip.text = element_blank()
    ) +
    theme(axis.line = element_line(color = "black")) +
    theme(axis.text = element_text(size = 11)) +
    theme(axis.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 11)) +
    theme(legend.position="bottom")+
    xlab("Year")+
    ylab("Number of \n observations")
  
  model_names <- c("Persistent Spatial", "Persistent Spatial + Year", "Year+Temp+Salinity", "Temp+Salinity+Spatio-temporal")
  test_year <- unique(test_data$year)
  test_region <- unique(test_data$region)
  pred_plots <- list()
  resid_plots <- list()
  pred_obs <- list()
  for (i in 1:length(preds)){
    try(pred_plots[[i]] <-  ggplot(us_coast_proj) + geom_sf() +
          geom_point(preds[[i]], mapping=aes(x=X*1000, y=Y*1000, col=o2),
                     size = 1.0,
                     alpha = 1.0
          ) +
          ylim(ylims)+
          scale_x_continuous(breaks=lons, limits=xlims)+
          scale_colour_viridis_c(
            limits = c(0, 200),
            oob = scales::squish,
            name = bquote(O[2]~Predictions),
            breaks = c(0, 100, 200)
          ) +
          labs(x = "Longitude", y = "Latitude") +
          theme_bw() +
          theme(
            panel.grid.major = element_blank()
            ,
            panel.grid.minor = element_blank()
            ,
            panel.border = element_blank()
            ,
            strip.background = element_blank()
            ,
            strip.text = element_blank()
          ) +
          theme(axis.line = element_line(color = "black")) +
          theme(axis.text = element_text(size = 11)) +
          theme(axis.title = element_text(size = 12)) +
          theme(legend.text = element_text(size = 11)) +
          theme(legend.position = "bottom") +
          guides(colour = guide_colourbar(title.position = "top", title.hjust =
                                            0.5)))
    
    try(resid_plots[[i]] <-  ggplot(us_coast_proj) + geom_sf() +
          geom_point(preds[[i]], mapping=aes(x=X*1000, y=Y*1000, col=residual),
                     size = 1.0,
                     alpha = 1.0
          ) +
          scale_colour_distiller(palette = "RdBu", limits = c(-50, 50)) +
          ylim(ylims)+
          scale_x_continuous(breaks=lons, limits=xlims)+
          #, limits = c(-40, 40), oob = scales::squish, name = bquote(O[2]), breaks = c(-40, 0, 40)) +
          labs(x = "Longitude", y = "Latitude") +
          theme_bw() +
          theme(
            panel.grid.major = element_blank()
            ,
            panel.grid.minor = element_blank()
            ,
            panel.border = element_blank()
            ,
            strip.background = element_blank()
            ,
            strip.text = element_blank()
          ) +
          theme(axis.line = element_line(color = "black")) +
          theme(axis.text = element_text(size = 11)) +
          theme(axis.title = element_text(size = 12)) +
          theme(legend.text = element_text(size = 11)) +
          theme(legend.position = "bottom") +
          guides(colour = guide_colourbar(title.position = "top", title.hjust =
                                            0.5)))
    
    try(pred_obs[[i]] <- ggplot(data = preds[[i]], aes(x = o2, y = est, col = latitude)) +
          geom_point() +
          scale_colour_distiller(
            # limits = c(31, 50),
            #oob = scales::squish,
            name = "latitude",
            palette="Greys"
            # breaks = c(35, 40, 45)
          ) +
          theme(legend.position = "bottom") +
          theme_bw() +
          theme(
            panel.grid.major = element_blank()
            ,
            panel.grid.minor = element_blank()
            ,
            panel.border = element_blank()
            ,
            strip.background = element_blank()
            ,
            strip.text = element_blank()
          ) +
          theme(axis.line = element_line(color = "black")) +
          theme(axis.text = element_text(size = 11)) +
          theme(axis.title = element_text(size = 12)) +
          theme(legend.text = element_text(size = 11)) +
          labs(x = "Observed", y = "Predicted") +
          geom_abline(intercept = 0, slope = 1)+
          theme(legend.position="bottom"))
  }
  # plot residuals vs. prediction
  #  resid_vs_pred <- ggplot(data = test_predict_O2, aes(x = (est), y = residual, col = Y)) +
  # geom_point() +
  #scale_colour_viridis_c(
  #  limits = c(31, 50),
  #  oob = scales::squish,
  #  name = "latitude",
  #  breaks = c(35, 40, 45)
  # ) +
  # ggtitle(test_year) +
  # labs(x = "Predicted", y = "Residual") +
  # theme(legend.position = "none")
  
  figure1 <- ggarrange(data_map, dat_available, labels=c("A", "B"),
                       ncol = 2, nrow = 1)
  figure1 <- annotate_figure(figure1, left="Training Data", fig.lab.size=14, fig.lab.face="bold")
  figure2 <- ggarrange(pred_plots[[1]], resid_plots[[1]], pred_obs[[1]], ncol=3, nrow=1, legend="none", labels=c("C", "D", "E"))
  figure2 <- annotate_figure(figure2, left=paste(model_names[1]), fig.lab.size=14, fig.lab.face="bold")
  figure3 <- ggarrange(pred_plots[[2]], resid_plots[[2]], pred_obs[[2]], ncol=3, nrow=1, legend="none")
  figure3 <- annotate_figure(figure3, left=paste(model_names[2]), fig.lab.size=14, fig.lab.face="bold")
  figure4 <- ggarrange(pred_plots[[3]], resid_plots[[3]], pred_obs[[3]], ncol=3, nrow=1, legend="none")
  figure4 <- annotate_figure(figure3, left=paste(model_names[3]), fig.lab.size=14, fig.lab.face="bold")
  figure5 <- try(ggarrange(pred_plots[[4]], resid_plots[[4]], pred_obs[[4]], ncol=3, nrow=1))
  figure5 <- try(annotate_figure(figure5, left=paste(model_names[4]), fig.lab.size=14, fig.lab.face="bold"))
  if(!is.list(figure5)){
    figure4 <- ggarrange(ggplot(), ggplot(), ggplot())
  }
  
  figure <- ggarrange(figure1, figure2, figure3, figure4, figure5, ncol=1, nrow=5, heights=c(1,1,1,1, 1.25), align="h")
  annotate_figure(figure, top=paste(test_year), fig.lab.size=18, fig.lab.face="bold")
  
  ggsave(
    paste("outputs/plots/", test_region, "/preds/preds_", test_year, ".pdf", sep=""),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height = 11,
    units = c("in"),
    dpi = 600,
    limitsize = TRUE
  )
  return(figure)
}

plot_glorys <- function(preds, dat.2.use){
  #Separate test and training data, predictions, and models from output list
  #Set latitude and longitude
  xlims <- c(min(dat.2.use$X)*1000, max(dat.2.use$X)*1000)
  ylims <- c(min(dat.2.use$Y)*1000, max(dat.2.use$Y)*1000)
  lats <- c(round(min(dat.2.use$latitude)),  round(max(dat.2.use$latitude)))
  lons <- c(round(min(dat.2.use$longitude)+2), round(max(dat.2.use$longitude)))
  test_region <- unique(preds$region)
  test_year <- unique(preds$year)
  
  try(pred_plot <-  ggplot(us_coast_proj) + geom_sf() +
        geom_point(preds, mapping=aes(x=X*1000, y=Y*1000, col=o2),
                   size = 1.0,
                   alpha = 1.0
        ) +
        ylim(ylims)+
        scale_x_continuous(breaks=lons, limits=xlims)+
        scale_colour_viridis_c(
          limits = c(0, 200),
          oob = scales::squish,
          name = bquote(O[2]~Predictions),
          breaks = c(0, 100, 200)
        ) +
        labs(x = "Longitude", y = "Latitude") +
        theme_bw() +
        theme(
          panel.grid.major = element_blank()
          ,
          panel.grid.minor = element_blank()
          ,
          panel.border = element_blank()
          ,
          strip.background = element_blank()
          ,
          strip.text = element_blank()
        ) +
        theme(axis.line = element_line(color = "black")) +
        theme(axis.text = element_text(size = 11)) +
        theme(axis.title = element_text(size = 12)) +
        theme(legend.text = element_text(size = 11)) +
        theme(legend.position = "bottom") +
        guides(colour = guide_colourbar(title.position = "top", title.hjust =
                                          0.5)))+
    ggtitle("Predictions")
  
  try(resid_plot <-  ggplot(us_coast_proj) + geom_sf() +
        geom_point(preds, mapping=aes(x=X*1000, y=Y*1000, col=residual),
                   size = 1.0,
                   alpha = 1.0
        ) +
        scale_colour_distiller(palette = "RdBu", limits = c(-50, 50)) +
        ylim(ylims)+
        scale_x_continuous(breaks=lons, limits=xlims)+
        #, limits = c(-40, 40), oob = scales::squish, name = bquote(O[2]), breaks = c(-40, 0, 40)) +
        labs(x = "Longitude", y = "Latitude") +
        theme_bw() +
        theme(
          panel.grid.major = element_blank()
          ,
          panel.grid.minor = element_blank()
          ,
          panel.border = element_blank()
          ,
          strip.background = element_blank()
          ,
          strip.text = element_blank()
        ) +
        theme(axis.line = element_line(color = "black")) +
        theme(axis.text = element_text(size = 11)) +
        theme(axis.title = element_text(size = 12)) +
        theme(legend.text = element_text(size = 11)) +
        theme(legend.position = "bottom") +
        guides(colour = guide_colourbar(title.position = "top", title.hjust =
                                          0.5)))+
    ggtitle("Prediction Residuals")
  
  try(pred_obs <- ggplot(data = preds, aes(x = o2, y = est, col = latitude)) +
        geom_point() +
        scale_colour_distiller(
          # limits = c(31, 50),
          #oob = scales::squish,
          name = "latitude",
          palette="Greys"
          # breaks = c(35, 40, 45)
        ) +
        theme(legend.position = "bottom") +
        theme_bw() +
        theme(
          panel.grid.major = element_blank()
          ,
          panel.grid.minor = element_blank()
          ,
          panel.border = element_blank()
          ,
          strip.background = element_blank()
          ,
          strip.text = element_blank()
        ) +
        theme(axis.line = element_line(color = "black")) +
        theme(axis.text = element_text(size = 11)) +
        theme(axis.title = element_text(size = 12)) +
        theme(legend.text = element_text(size = 11)) +
        labs(x = "Observed", y = "Predicted") +
        geom_abline(intercept = 0, slope = 1)+
        theme(legend.position="bottom"))
  
  # plot residuals vs. prediction
  #  resid_vs_pred <- ggplot(data = test_predict_O2, aes(x = (est), y = residual, col = Y)) +
  # geom_point() +
  #scale_colour_viridis_c(
  #  limits = c(31, 50),
  #  oob = scales::squish,
  #  name = "latitude",
  #  breaks = c(35, 40, 45)
  # ) +
  # ggtitle(test_year) +
  # labs(x = "Predicted", y = "Residual") +
  # theme(legend.position = "none")
  
  figure <- ggarrange(pred_plot, resid_plot, pred_obs, ncol=3, nrow=1, labels=c("A", "B", "C"), widths=c(1,1,2), heights=c(1,1,0.3))
  annotate_figure(figure, top=paste(test_year), fig.lab.size=18, fig.lab.face="bold")
  
  ggsave(
    paste("outputs/plots/", test_region, "/glorys/_glorys", test_year, ".pdf", sep=""),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height = 5.5,
    units = c("in"),
    dpi = 600,
    limitsize = TRUE
  )
  return(figure)
}

##Function to plot marginal effects of spatio-temporal model
plot_marginal_effects <- function(models,preds, dat.2.use, i){
  m <- try(models[[i]])
  if(is.list(m)){ #don't do all this if the model failed anyway
    xlims <- c(min(dat.2.use$X)*1000, max(dat.2.use$X)*1000)
    ylims <- c(min(dat.2.use$Y)*1000, max(dat.2.use$Y)*1000)
    lats <- c(round(min(dat.2.use$latitude)),  round(max(dat.2.use$latitude)))
    lons <- c(round(min(dat.2.use$longitude)+2), round(max(dat.2.use$longitude)))
    #Prediction dataframe for epsilon and omega
    preds <- try(preds[[i]])
    test_region <- unique(preds$region)
    test_year <- unique(preds$year)
    #marginal effects
    pdf(paste("outputs/plots/", test_region, "/margeffects/margeffects_", test_year,".pdf", sep=""))
    par(mfrow=c(2,2))
    visreg(m, "sigma0")
    visreg(m, "temp")
    visreg(m, "doy")
    visreg(m, "depth_ln")
    dev.off()
    
    omega <- ggplot(us_coast_proj) + geom_sf() +
      geom_point(preds, mapping=aes(x=X*1000, y=Y*1000, col=omega_s),
                 size = 1.0,
                 alpha = 1.0
      ) +
      scale_colour_distiller(palette = "RdBu") +
      ylim(ylims)+
      scale_x_continuous(breaks=lons, limits=xlims)+
      #, limits = c(-40, 40), oob = scales::squish, name = bquote(O[2]), breaks = c(-40, 0, 40)) +
      labs(x = "Longitude", y = "Latitude") +
      theme_bw() +
      theme(
        panel.grid.major = element_blank()
        ,
        panel.grid.minor = element_blank()
        ,
        panel.border = element_blank()
        ,
        strip.background = element_blank()
        ,
        strip.text = element_blank()
      ) +
      theme(axis.line = element_line(color = "black")) +
      theme(axis.text = element_text(size = 11)) +
      theme(axis.title = element_text(size = 12)) +
      theme(legend.text = element_text(size = 11)) +
      theme(legend.position = "bottom") +
      guides(colour = guide_colourbar(title.position = "top", title.hjust =
                                        0.5))
    try(epsilon <- ggplot(us_coast_proj) + geom_sf() +
          geom_point(preds, mapping=aes(x=X*1000, y=Y*1000, col=epsilon_st),
                     size = 1.0,
                     alpha = 1.0
          ) +
          scale_colour_distiller(palette = "RdBu") +
          ylim(ylims)+
          scale_x_continuous(breaks=lons, limits=xlims)+
          #, limits = c(-40, 40), oob = scales::squish, name = bquote(O[2]), breaks = c(-40, 0, 40)) +
          labs(x = "Longitude", y = "Latitude") +
          theme_bw() +
          theme(
            panel.grid.major = element_blank()
            ,
            panel.grid.minor = element_blank()
            ,
            panel.border = element_blank()
            ,
            strip.background = element_blank()
            ,
            strip.text = element_blank()
          ) +
          theme(axis.line = element_line(color = "black")) +
          theme(axis.text = element_text(size = 11)) +
          theme(axis.title = element_text(size = 12)) +
          theme(legend.text = element_text(size = 11)) +
          theme(legend.position = "bottom") +
          guides(colour = guide_colourbar(title.position = "top", title.hjust =
                                            0.5)))
    
    ggsave(
      paste("outputs/plots/", test_region, "/spatiotemp_var/spatiotemp_var_", test_year, ".pdf", sep=""),
      plot = last_plot(),
      device = NULL,
      path = NULL,
      scale = 1,
      width = 8.5,
      height = 11,
      units = c("in"),
      dpi = 600,
      limitsize = TRUE
    )
    
  }
  if(!is.list(m)){
    return(paste("model not fit"))
  }
}

#Convert GLORYS from .nc format 
convert_glorys <- function(file_name, do_threshold, filter_depth, filter_time) {
  nc <- tidync(file_name)
  if(filter_depth){
    nc_df <- nc %>%
      hyper_tibble %>%
      group_by(longitude, latitude) %>%
      filter(depth == max(depth)) %>%
      ungroup() 
  }
  if(!filter_depth){
    nc_df <- nc %>%
      hyper_tibble %>%
      group_by(longitude, latitude) %>%
      ungroup() 
    
  }
  # replace DO below threshold with the threshold level
  nc_df <- nc_df %>%
    mutate(o2 = case_when(
      o2 < do_threshold ~ do_threshold,
      TRUE ~ o2  # Keep other values unchanged
    ))
  
  #Make columns numeric
  nc_df$longitude <- as.numeric(nc_df$longitude)
  nc_df$latitude <- as.numeric(nc_df$latitude)
  nc_df$depth <- as.numeric(nc_df$depth)
  
  # remove large list from memory
  rm(nc)
  nc_df$time <- as.Date(substr(nc_df$time, 1, 10))
  if(filter_time){
    # nc_df$time <- round(nc_df$time)
    days <- unique(nc_df$time)
    n_days <- length(days)
    first_day <- days[1]
    last_day <- days[n_days]
    days.2.use <- seq(first_day, last_day, by = 10)
    nc_df <- nc_df %>%
      filter(time %in% days.2.use)
  }
  
  #nc_df$time <- round(nc_df$time)
  #nc_df$time <- (as_datetime("1950-01-01")+hours(nc_df$time))
  nc_df <- nc_df %>%
    st_as_sf(coords=c('longitude','latitude'),crs=4326,remove = F) %>%  
    st_transform(crs = "+proj=utm +zone=10 +datum=WGS84 +units=km") %>% 
    mutate(X=st_coordinates(.)[,1],Y=st_coordinates(.)[,2]) %>% 
    st_set_geometry(NULL)
  nc_df$doy <- as.POSIXlt(nc_df$time, format = "%Y-%b-%d")$yday
  nc_df$year <- year(nc_df$time)
  nc_df$month <- month(nc_df$time)
  return(nc_df)
}


isolate_preds <- function(x, region){
  region <- region
  if(region=="bc"|region=="cc"){
    preds <- try(x[["predictions"]][[6]])
  }
  if(region=="goa"|region=="ebs"|region=="ai"){
    preds <- try(x[["predictions"]][[3]])
  }
  return(preds)
}

isolate_preds2 <- function(x){
  preds <- try(x[["predictions"]])
  return(preds)
}


combine_preds <- function(x, glorys, region){
  if(!glorys){
    preds <- lapply(x, isolate_preds, region)
  }
  if(glorys){
    preds <- lapply(x, isolate_preds2)
  }
  row_lt2 <- which(sapply(preds, is.data.frame))
  preds <- preds[row_lt2]
  preds <- bind_rows(preds)
  return(preds)
}


gobh_interpolate <- function(dat, test_region, filter_depth, filter_time, gloryswd, basewd, do_threshold){
  predictions <- list()
  yearlist <- sort(unique(dat$year))
  for (i in 1:length(yearlist)) {
    test_year <- yearlist[i]
    print(test_year)
    ##Trawl testing data
    data <- dat %>%
      filter(year==test_year)
    data <- as.data.frame(data)
    ##Pull GLORYS data
    setwd(gloryswd)
    if(test_region=="cc"){
      files <- list.files("wc_o2", pattern=paste(test_year))
      files <- paste("wc_o2/", files, sep="")
    }
    
    if(test_region=="bc"){
      files <- list.files("bc_o2", pattern=paste(test_year))
      files <- paste("bc_o2/", files, sep="")
    }
    
    if(test_region=="ebs"|test_region=="goa"|test_region=="ai"){
      files <- list.files("alaska_o2_combined", pattern=paste(test_year))
      files <- paste("alaska_o2_combined/", files, sep="")
    }
    
    #Convert GLOYRS file from .nc format to dataframe in right format
    print("converting GLORYS")
    if(length(files)==1){
      glorys <- convert_glorys(files, do_threshold, filter_depth, filter_time)
    }
    
    #Use function to get data for all files if multiple files in year (this is mostly a thing for Alaska)
    if(length(files)>1){
      # get first year of glorys
      glorys <- convert_glorys(files[1], do_threshold, filter_depth, filter_time)
      # get remaining years of glorys and combine into single data frame
      for (j in 2:length(files)) {
        tmp_glorys <- convert_glorys(files[j], do_threshold, filter_depth, filter_time)
        glorys <- rbind(glorys, tmp_glorys)
      }
    }
    
    ##Pull in survey extent polygon from GLORYS data
    # Regional polygon
    poly <- filter(regions.hull, region==test_region)
    #Convert GLORYS to an sf
    glorys_sf <-  st_as_sf(glorys, coords = c("longitude", "latitude"), crs = st_crs(4326))
    # pull out observations within each region
    region_dat  <- st_filter(glorys_sf, poly)
    region_dat <- as.data.frame(region_dat)
    
    #Log depth
    region_dat$depth_ln <- log(region_dat$depth)
    
    #Set working directory back to base for saving
    setwd(basewd)
    
    #Mesh
    spde <- make_mesh(data = region_dat,
                      xy_cols = c("X", "Y"),
                      cutoff = 45)
    
    region_dat$o2 <- region_dat$o2/100
    
    #Fit model
    print("fitting GLOBH model intercept")
    m <- try(sdmTMB(formula = o2 ~ 1 +s(depth_ln) + s(doy),
                    mesh = spde,
                    data = region_dat, 
                    family = gaussian(),
                    spatial = "on",
                    spatiotemporal  = "off"))
    preds <- try(predict(m, data))
    if(is.list(m)){
      predictions[[i]] <- preds
    }
  }
  return(predictions)
}

load_all_hauls <- function() {
  install.packages("remotes")
  remotes::install_github("nwfsc-assess/nwfscSurvey")
  haul = nwfscSurvey::PullHaul.fn(SurveyName = "NWFSC.Combo")
  haul <- plyr::rename(haul, replace=c("salinity_at_gear_psu_der" = "sal", 
                                       "temperature_at_gear_c_der" = "temp", 
                                       "o2_at_gear_ml_per_l_der" = "o2",
                                       "depth_hi_prec_m" = "depth"))
  
  # read in the grid cell data from the survey design
  grid_cells = readxl::read_excel("data/Selection Set 2018 with Cell Corners.xlsx")
  grid_cells = dplyr::mutate(grid_cells,
                             depth_min = as.numeric(unlist(strsplit(grid_cells$Depth.Range,"-"))[1]),
                             depth_max = as.numeric(unlist(strsplit(grid_cells$Depth.Range,"-"))[2]))
  
  # convert grid_cells to sp object
  grid = SpatialPoints(cbind(grid_cells$Cent.Long,grid_cells$Cent.Lat),
                       proj4string = CRS("+proj=longlat +datum=WGS84"))
  r = raster::rasterize(x=grid, y = raster(nrow=length(unique(grid_cells$Cent.Lat)),
                                           ncol=length(unique(grid_cells$Cent.Long))))
  rasterToPoints(r)
  
  raster = aggregate(r, fact = 2)
  raster = projectRaster(raster, crs = "+proj=tmerc +lat_0=31.96 +lon_0=-121.6 +k=1 +x_0=390000 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  # create matrix of point data with coordinates and depth from raster
  grid = as.data.frame(rasterToPoints(raster))
  
  # Figure out the grid cell corresponding to each tow location
  haul$Cent.Lat = NA
  haul$Cent.Lon = NA
  haul$Cent.ID = NA
  for(i in 1:nrow(haul)) {
    indx = which(grid_cells$NW.LAT > haul$latitude_dd[i] &
                   grid_cells$SW.LAT < haul$latitude_dd[i] &
                   grid_cells$NW.LON < haul$longitude_dd[i] &
                   grid_cells$NE.LON > haul$longitude_dd[i])
    if(length(indx) > 0) {
      haul$Cent.ID[i] = grid_cells$Cent.ID[indx]
      haul$Cent.Lat[i] = grid_cells$Cent.Lat[indx]
      haul$Cent.Lon[i] = grid_cells$Cent.Long[indx]
    }
  }
  
  # project lat/lon to UTM, after removing missing values and unsatisfactory hauls
  haul = haul %>% filter(!is.na(Cent.Lon), performance == "Satisfactory")
  
  haul_trans = haul
  coordinates(haul_trans) <- c("Cent.Lon", "Cent.Lat")
  proj4string(haul_trans) <- CRS("+proj=longlat +datum=WGS84")
  newproj = paste("+proj=utm +zone=10 ellps=WGS84")
  haul_trans <- spTransform(haul_trans, CRS(newproj))
  haul_trans = as.data.frame(haul_trans)
  haul_trans$Cent.Lon = haul_trans$Cent.Lon/10000
  haul_trans$Cent.Lat = haul_trans$Cent.Lat/10000
  haul_trans$year = as.numeric(substr(haul_trans$date_yyyymmdd,1,4))
  
  haul$X = haul_trans$Cent.Lon
  haul$Y = haul_trans$Cent.Lat
  haul$year = haul_trans$year
  #haul$year_centered = haul$year - mean(unique(haul$year))
  
  return(haul)
}
library(ggplot2)

setwd("~/Dropbox/GitHub/wsg-choke-species")

#ggplot themes
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

#Load functions
source("code/helper_funs.R")

est_pa <- read.csv("output/presence_absence/breakpoint_est.csv")
est_pa$type <-"presence_absence"
est_density <- read.csv("output/region_comp/breakpoint_est.csv")
est_density$type <- "density"

bp_est2 <- bind_rows(est_pa, est_density)

bp_est2$data <- factor(bp_est2$data, levels=c("cc", "bc", "goa", "ebs", "coastwide"))
labs <- c("British Columbia", "California Current", "Gulf of Alaska", "Eastern Bering Sea", "Coastwide")
names(labs) <- c("bc", "cc", "goa", "ebs", "coastwide")

ggplot(bp_est2, aes(y=species, x=breakpt, colour=data))+
  geom_point(aes(colour=type, shape=data), size=2, position=ggstance::position_dodgev(height=0.5))+
  #Can add shape back
  geom_linerange(aes(xmin = breakpt_se1, xmax = breakpt_se2, colour=type),  position=ggstance::position_dodgev(height=0.5), size=1, alpha=0.5)+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  # xlim(0,50)+
#  scale_colour_manual(values=c("#F8766D","#7CAE00", "#00BFC4",  "#C77CFF", "#00B0F6"))+
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
  xlab("Estimated Breakpoint Metabolic Index") +
  ylab("Species")

ggsave(
  paste("output/density_presence_comparison.png"),
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


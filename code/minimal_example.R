remove.packages("sdmTMB")
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE,  ref="newbreakpt")
library(sdmTMB)

library(ggeffects)
library(visreg)
library(ggplot2)

set.seed(9876)

dat <- readRDS("output/region_comp/rex sole_coastwide_dat.rds")

#Make mesh
bnd <- INLA::inla.nonconvex.hull(cbind(dat$X, dat$Y), 
                                   convex = -0.05)
inla_mesh <- INLA::inla.mesh.2d(
      boundary = bnd,
      max.edge = c(150, 1000),
      offset = -0.1, # default -0.1
      cutoff = 50,
      min.angle = 5 # default 21
    )

spde <- make_mesh(dat, c("X", "Y"), mesh = inla_mesh)

#Priors
priors <- sdmTMBpriors(
      matern_s = pc_matern(
        range_gt = 50, range_prob = 0.05, #A value one expects the range is greater than with 1 - range_prob probability.
        sigma_lt = 25, sigma_prob = 0.05 #A value one expects the marginal SD (sigma_O or sigma_E internally) is less than with 1 - sigma_prob probability.
      ),
        matern_st = pc_matern(
        range_gt = 50, range_prob = 0.05,
        sigma_lt = 25, sigma_prob = 0.05
      ),
      #  ar1_rho = normal(0.7,0.1),
      #tweedie_p = normal(1.5,0.2)
    )

# refactor to avoid identifiability errors
dat$region <- as.factor(as.character(dat$region))
dat$year <- as.factor(as.character(dat$year))
  
#Model formula
formula = "catch ~ -1 + year+ region +breakpt(mi1_s)+ +temp_scaled + temp_scaled2+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
##Add survey for species/ w/ IPHC data
formula = "catch ~ -1 + year+ region + survey+breakpt(mi1_s)+ +temp_scaled + temp_scaled2+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"

#Fit model
m <-try(sdmTMB(
    formula = as.formula(formula),
    mesh = spde,
    time = "year",
    family = tweedie(link = "log"),
    data = dat,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = "iid",
    control = sdmTMBcontrol(normalize = TRUE,
                            multiphase = TRUE,
                            newton_loops = 3,
                            nlminb_loops = 2)
  ))

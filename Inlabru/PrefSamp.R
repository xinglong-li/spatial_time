library(dplyr)
library(sf)
library(sp)
library(reshape2)
library(INLA)

library(ggplot2)
library(inlabru)


# Prepare variables for the model ==================================================================

Data$time <- (Data$year - min(Data$year)) / (max(Data$year) - min(Data$year))
Data$locs <- coordinates(Data[, c("east", "north")])
Data$site_number <- as.numeric(as.factor(Data$site))

# Site-selection indicator
Data$R <- as.numeric(!is.na(Data$bsmoke))
Data$R_lag <- c(rep(NA, no_sites), Data$R[1:(dim(Data)[1]-no_sites)])

# Response variable for the auxiliary model
Data$zero <- 0

# Compute Euclidean distances between all the sites
dists <- spDists(cbind(BS$east, BS$north))

r <- 0.1 # The maximum radius of interest to be 10 km
Data$repulsion_ind <- 0

counter <- no_sites + 1
for (i in sort(unique(Data$year))[-1]) {
  # First extract the data at time i
  data_i <- Data[Data$year == i, ]
  # Compute the repulsion indicator. Was there a site at year i-1 within radius r of it?
  Data$repulsion_ind[counter:(counter+(no_sites - 1))] = rowSums(dists[, which(data_i$R_lag==1)] < r) > 0
  counter <- counter + no_sites
}

# Create new mesh
edg_in = 0.3 # The range set to be [0.03, 0.15]
edg_out = 2 * edg_in
mesh <- fm_mesh_2d_inla(loc = cbind(Data$east, Data$north),
                        boundary = bord,
                        offset = c(2*edg_in, edg_out),
                        max.edge = 2*c(edg_in, edg_out),
                        cutoff = edg_in,
                        min.angle = 21)
mesh$n 
ggplot(Data) + gg(mesh) + geom_point(aes(x = east, y = north), color = "blue") + coord_fixed()

# Maybe we should consider set the PC prior using data infor
spde_obj <- inla.spde2.pcmatern(mesh = mesh,
                                alpha = 2,
                                prior.range = c(edg_in, 0.01),
                                prior.sigma = c(sqrt(mean(var_annual$var_bs)), 0.01),
                                constr = T)

# Build inlabru model ==============================================================================

# Joint independent model --------------------------------------------------------------------------

comp_joint_indep <- ~ Intercept_obs(1) + # Components for observation model
  Time_obs_1(time) +
  Time_obs_2(time^2) +
  Random_obs_0(site_number, model = "iid2d", n = no_sites*2, constr=TRUE) +
  Random_obs_1(site_number, weights = time, copy = "Random_obs_0") +
  Spatial_obs_0(locs, model = spde_obj) +
  Spatial_obs_1(locs, weights = time, model = spde_obj) +
  Spatial_obs_2(locs, weights = time^2, model = spde_obj) +
  
  # Components for site selection model
  Intercept_slc(1) + Time_slc_1(time) +
  Time_slc_2(time^2) +
  R_lag_slc(R_lag) +
  Repuls_slc(repulsion_ind) +
  AR_slc(year, model='ar1', hyper=list(theta1=list(prior="pcprec",param=c(2, 0.01)))) +
  Spatial_slc(locs, model = spde_obj) 

like_obs <- like(
  formula = bsmoke ~ Intercept_obs + Time_obs_1 + Time_obs_2 + 
    Random_obs_0 + Random_obs_1 + 
    Spatial_obs_0 + Spatial_obs_1 + Spatial_obs_2,
  family = "gaussian",
  data = Data
)

like_slc <- like(
  formula = R ~ Intercept_slc + Time_slc_1 + Time_slc_2 + 
    R_lag_slc + Repuls_slc + AR_slc + Spatial_slc,
  family = "binomial",
  Ntrials = rep(1, times = length(Data$R)),
  data = Data
)

bru_options_reset()
bru_options_set(bru_max_iter = 5,
                control.inla = list(strategy = "gaussian", int.strategy = 'eb'),
                bru_verbose = T)

fit_bru_joint_indep <- bru(comp_joint_indep, like_obs, like_slc)


# Joint preferential sampling model ----------------------------------------------------------------

comp_aux <- ~ Intercept_obs(1) + # Components for observation model
  Time_obs_1(time) +
  Time_obs_2(time^2) +
  Random_obs_0(site_number, model = "iid2d", n = no_sites*2, constr=TRUE) +
  Random_obs_1(site_number, weights = time, copy = "Random_obs_0") +
  Spatial_obs_0(locs, model = spde_obj) +
  Spatial_obs_1(locs, weights = time, model = spde_obj) +
  Spatial_obs_2(locs, weights = time^2, model = spde_obj) +
  
  # Components for site selection model
  Intercept_slc(1) + Time_slc_1(time) +
  Time_slc_2(time^2) +
  R_lag_slc(R_lag) +
  Repuls_slc(repulsion_ind) +
  AR_slc(year, model='ar1', hyper=list(theta1=list(prior="pcprec",param=c(2, 0.01)))) +
  Spatial_slc(locs, model = spde_obj) +
  Comp_share1(site_number, copy = 'Comp_aux1', fixed = FALSE) + 
  Comp_share2(site_number, copy = 'Comp_aux2', fixed = FALSE) +
  
  # Components for 1st auxiliary model
  Random_aux1_0(site_number, copy = "Random_obs_0", fixed = TRUE) +
  Random_aux1_1(site_number, weights = time, copy = "Random_obs_1", fixed = TRUE) +
  Comp_aux1(site_number, model = 'iid', hyper = list(prec = list(initial = -20, fixed=TRUE))) +
  
  # Components for 2nd auxiliary model
  Spatial_aux2_0(locs, copy = "Spatial_obs_0", fixed = TRUE) +
  Spatial_aux2_1(locs, weights = time, copy = "Spatial_obs_1", fixed = TRUE) +
  Spatial_aux2_2(locs, weights = time^2, copy = "Spatial_obs_2", fixed = TRUE) +
  Comp_aux2(site_number, model = 'iid', hyper = list(prec = list(initial = -20, fixed=TRUE)))


# All likelihoods
like_obs <- like(
  formula = bsmoke ~ Intercept_obs + Time_obs_1 + Time_obs_2 + 
    Random_obs_0 + Random_obs_1 + 
    Spatial_obs_0 + Spatial_obs_1 + Spatial_obs_2,
  family = "gaussian",
  data = Data
)

like_slc_share <- like(
  formula = R ~ Intercept_slc + Time_slc_1 + Time_slc_2 + 
    R_lag_slc + Repuls_slc + AR_slc + Spatial_slc + 
    Comp_share1 + Comp_share2,
  family = "binomial",
  Ntrials = rep(1, times = length(Data$R)),
  data = Data
)

like_aux_1 <- like(
  formula = zero ~ Random_aux1_0 + Random_aux1_1 + Comp_aux1,
  family = "gaussian",
  data = Data
)

like_aux_2 <- like(
  formula = zero ~  Spatial_aux2_0 +  Spatial_aux2_1 + Spatial_aux2_2 + Comp_aux2,
  family = "gaussian",
  data = Data
)

bru_options_reset()
bru_options_set(bru_max_iter = 20,
                control.inla = list(strategy = "gaussian", int.strategy = 'eb'),#)# h = 1e-5),
                # control.mode = list(theta = c(fit_bru_aux_1$mode$theta, 0),
                #                     restart = TRUE),
                control.family = list(
                  list(),
                  list(),
                  list(hyper = list(prec = list(initial = 20, fixed=TRUE))),
                  list(hyper = list(prec = list(initial = 20, fixed=TRUE)))),
                bru_verbose = T)

start_time_aux <- Sys.time()
fit_bru_aux <- bru(comp_aux, like_obs, like_slc_share, like_aux_1, like_aux_2)
end_time_aux <- Sys.time()
runtime_init <- end_time_aux - start_time_aux


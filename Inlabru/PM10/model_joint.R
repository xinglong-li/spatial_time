library(dplyr)
library(sf)
library(sp)
library(reshape2)
library(INLA)
library(ggplot2)
library(inlabru)



# Prepare the data =================================================================================

# Load the (scaled) locations of sites and the border of California(the unit is 100km)
# The data also excludes sites that contain outlier records (site_number = 0030 and 0022).
path <- "./PM10/"
PM10s0 <- readRDS(sprintf("%sPM10s_CA_summary_scaled.rds", path))
print(sprintf("Total number of annual measurements: %s", dim(PM10s0)[1]))
print(sprintf("Total number of sites: %s", length(unique(PM10s0$site_number))))
print(sprintf("The sd of East coordinates is: %s", sd(PM10s0$east)))
print(sprintf("The sd of North coordinates is: %s", sd(PM10s0$north)))

CA_border <- readRDS(sprintf("%sCA_border_scaled.rds", path))

ggplot(PM10s0) + gg(CA_border) + geom_point(aes(x = east, y = north)) + coord_equal()

mean_annually <- group_by(PM10s0, year) %>%
  summarise(mean_pm = mean(annual_mean)) 
mean_pm_annually <- mean(mean_annually$mean_pm)
PM10s0$annual_mean <- log(PM10s0$annual_mean) - log(mean_pm_annually)

no_sites = length(unique(PM10s0$site_number))
no_T = length(unique(PM10s0$year))

# Reform the data so that each site has num_of_year rows, the empty records are filled with NAs
PM10s_flat <- dcast(PM10s0, site_number + north + east ~ year, value.var = "annual_mean")
PM10s = melt(PM10s_flat, id.vars = c(1,2,3), variable.name = 'year', value.name = 'annual_mean')
PM10s$year <- as.numeric(as.character(factor(PM10s$year, labels = 1985:2022)))

stopifnot(dim(PM10s)[1] == no_sites * no_T)

subsmaple_plotting_data = PM10s[which(PM10s$site_number %in% unique(PM10s$site_number)[
  sample.int(111, size = 30)]), ]
means_plot = ggplot(data = subsmaple_plotting_data, aes(x = year, y = annual_mean)) +
  geom_line(aes(group = site_number, colour = site_number)) + xlab("Year") + ylab("log(PM10)")
means_plot

var_annual <- group_by(PM10s, year) %>% 
  summarise(var_pm = var(annual_mean, na.rm = T))
variance_plot =  ggplot(data = var_annual, aes(x = year, y = var_pm)) +
  geom_line() + geom_smooth() + xlab("Year") + ylab("Variance log(PM10)") 
ggtitle('A plot of the variance of the log annual means with fitted smoother')
variance_plot

# Prepare variables for the model ==================================================================

PM10s$time <- (PM10s$year - min(PM10s$year)) / (max(PM10s$year) - min(PM10s$year))
PM10s$locs <- coordinates(PM10s[, c("east", "north")])
PM10s$site_number <- as.numeric(as.factor(PM10s$site_number))

# Site-selection indicator
PM10s$R <- as.numeric(!is.na(PM10s$annual_mean))
PM10s$R_lag <- c(rep(NA, no_sites), PM10s$R[1:(dim(PM10s)[1]-no_sites)])

# Response variable for the auxiliary model
PM10s$zero <- 0

# Compute Euclidean distances between all the sites
dists <- spDists(cbind(PM10s$east[1:no_sites], PM10s$north[1:no_sites]))

r <- 0.1 # The maximum radius of interest to be 10 km
PM10s$repulsion_ind <- 0

counter <- no_sites + 1
for (i in sort(unique(PM10s$year))[-1]) {
  # First extract the data at time i
  data_i <- PM10s[PM10s$year == i, ]
  # Compute the repulsion indicator. Was there a site at year i-1 within radius r of it?
  PM10s$repulsion_ind[counter:(counter+(no_sites - 1))] = rowSums(dists[, which(data_i$R_lag==1)] < r) > 0
  counter <- counter + no_sites
}

# Create mesh
edge_in = 0.15 # 20km
edge_out = 2 * edge_in
mesh = fm_mesh_2d_inla(loc = cbind(PM10s$east, PM10s$north),
                       boundary = CA_border,
                       offset = c(2*edge_in, edge_out), 
                       max.edge = 2*c(edge_in, edge_out),
                       cutoff = edge_in,
                       min.angle = 21)
mesh$n
ggplot(PM10s) + gg(mesh) + geom_point(aes(x = east, y = north), color = "blue") + coord_fixed()

# Maybe we should consider set the PC prior using data info
spde_obj <- inla.spde2.pcmatern(mesh = mesh, 
                                alpha = 2, 
                                prior.range = c(edge_in, 0.01),
                                prior.sigma = c(sqrt(mean(var_annual$var_pm)), 0.1),
                                constr = T)

# Build inlabru model ==============================================================================

# Model for the observation process ----------------------------------------------------------------

comp_obs <- ~ Intercept_obs(1) +
  Time_obs_1(time) +
  Time_obs_2(time^2) +
  Random_obs_0(site_number, model = "iid2d", n = no_sites*2, constr=TRUE) +
  Random_obs_1(site_number, weights = time, copy = "Random_obs_0") +
  Spatial_obs_0(locs, model = spde_obj) +
  Spatial_obs_1(locs, weights = time, model = spde_obj) +
  Spatial_obs_2(locs, weights = time^2, model = spde_obj)

like_obs <- like(
  formula = annual_mean ~ Intercept_obs + Time_obs_1 + Time_obs_2 + 
    Random_obs_0 + Random_obs_1 + 
    Spatial_obs_0 + Spatial_obs_1 + Spatial_obs_2,
  family = "gaussian",
  data = PM10s
)

bru_options_reset()
bru_options_set(bru_max_iter = 5,
                control.inla = list(strategy = "gaussian", int.strategy = 'eb'),
                bru_verbose = T)

fit_bru_obs <- bru(comp_obs, like_obs)

# Model for the selection process ------------------------------------------------------------------

comp_slc <- ~ Intercept_slc(1) + 
  Time_slc_1(time) +
  Time_slc_2(time^2) +
  R_lag_slc(R_lag) +
  Repuls_slc(repulsion_ind) +
  AR_slc(year, model='ar1', hyper=list(theta1=list(prior="pcprec",param=c(2, 0.01)))) +
  Spatial_slc(locs, model = spde_obj) 

like_slc <- like(
  formula = R ~ Intercept_slc + Time_slc_1 + Time_slc_2 + 
    R_lag_slc + Repuls_slc + AR_slc + Spatial_slc,
  family = "binomial",
  Ntrials = rep(1, times = length(PM10s$R)),
  data = PM10s
)

bru_options_reset()
bru_options_set(bru_max_iter = 5,
                control.inla = list(strategy = "gaussian", int.strategy = 'eb'),
                bru_verbose = T)

fit_bru_slc <- bru(comp_slc, like_slc)

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
  formula = annual_mean ~ Intercept_obs + Time_obs_1 + Time_obs_2 + 
    Random_obs_0 + Random_obs_1 + 
    Spatial_obs_0 + Spatial_obs_1 + Spatial_obs_2,
  family = "gaussian",
  data = PM10s
)

like_slc <- like(
  formula = R ~ Intercept_slc + Time_slc_1 + Time_slc_2 + 
    R_lag_slc + Repuls_slc + AR_slc + Spatial_slc,
  family = "binomial",
  Ntrials = rep(1, times = length(PM10s$R)),
  data = PM10s
)

bru_options_reset()
bru_options_set(bru_max_iter = 10,
                control.inla = list(strategy = "gaussian", int.strategy = 'eb'),
                # control.mode = list(theta = c(fit_bru_obs$mode$theta, fit_bru_slc$mode$theta),
                # restart = TRUE),
                bru_verbose = T)

fit_bru_joint_indep <- bru(comp_joint_indep, like_obs, like_slc)

# Joint Auxiliary model 1 --------------------------------------------------------------------------

comp_aux_1 <- ~ Intercept_obs(1) + # Components for observation model
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
  
  # Components for 1st auxiliary model
  Random_aux1_0(site_number, copy = "Random_obs_0", fixed = TRUE) +
  Random_aux1_1(site_number, weights = time, copy = "Random_obs_1", fixed = TRUE) +
  Comp_aux1(site_number, model = 'iid', hyper = list(prec = list(initial = -20, fixed=TRUE)))


# All likelihoods
like_obs <- like(
  formula = annual_mean ~ Intercept_obs + Time_obs_1 + Time_obs_2 + 
    Random_obs_0 + Random_obs_1 + 
    Spatial_obs_0 + Spatial_obs_1 + Spatial_obs_2,
  family = "gaussian",
  data = PM10s
)

like_slc_share_1 <- like(
  formula = R ~ Intercept_slc + Time_slc_1 + Time_slc_2 + 
    R_lag_slc + Repuls_slc + AR_slc + Spatial_slc + Comp_share1,
  family = "binomial",
  Ntrials = rep(1, times = length(PM10s$R)),
  data = PM10s
)

like_aux_1 <- like(
  formula = zero ~ Random_aux1_0 + Random_aux1_1 + Comp_aux1,
  family = "gaussian",
  data = PM10s
)

bru_options_reset()
bru_options_set(bru_max_iter = 10,
                control.inla = list(strategy = "gaussian", int.strategy = 'eb'),#)# h = 1e-5),
                # control.mode = list(theta = c(fit_bru_obs$mode$theta, fit_bru_slc$mode$theta),
                #                     restart = TRUE),
                control.family = list(
                  list(),
                  list(),
                  list(hyper = list(prec = list(initial = 20, fixed=TRUE)))),
                bru_verbose = T)

start_time_aux1 <- Sys.time()
fit_bru_aux_1 <- bru(comp_aux_1, like_obs, like_slc_share_1, like_aux_1)
end_time_aux1 <- Sys.time()
runtime_aux_1 <- end_time_aux1 - start_time_aux1

# Joint Auxiliary model 2 --------------------------------------------------------------------------

comp_aux_2 <- ~ Intercept_obs(1) + # Components for observation model
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
  Comp_share2(site_number, copy = 'Comp_aux2', fixed = FALSE) +
  
  # Components for 2nd auxiliary model
  Spatial_aux2_0(locs, copy = "Spatial_obs_0", fixed = TRUE) +
  Spatial_aux2_1(locs, weights = time, copy = "Spatial_obs_1", fixed = TRUE) +
  Spatial_aux2_2(locs, weights = time^2, copy = "Spatial_obs_2", fixed = TRUE) +
  Comp_aux2(site_number, model = 'iid', hyper = list(prec = list(initial = -20, fixed=TRUE)))

# All likelihoods
like_obs <- like(
  formula = annual_mean ~ Intercept_obs + Time_obs_1 + Time_obs_2 + 
    Random_obs_0 + Random_obs_1 + 
    Spatial_obs_0 + Spatial_obs_1 + Spatial_obs_2,
  family = "gaussian",
  data = PM10s
)

like_slc_share_2 <- like(
  formula = R ~ Intercept_slc + Time_slc_1 + Time_slc_2 + 
    R_lag_slc + Repuls_slc + AR_slc + Spatial_slc + Comp_share2,
  family = "binomial",
  Ntrials = rep(1, times = length(PM10s$R)),
  data = PM10s
)

like_aux_2 <- like(
  formula = zero ~ Spatial_aux2_0 + Spatial_aux2_1 + Spatial_aux2_2 + Comp_aux2,
  family = "gaussian",
  data = PM10s
)

bru_options_reset()
bru_options_set(bru_max_iter = 10,
                control.inla = list(strategy = "gaussian", int.strategy = 'eb'),#)# h = 1e-5),
                # control.mode = list(theta = c(fit_bru_obs$mode$theta, fit_bru_slc$mode$theta),
                #                     restart = TRUE),
                control.family = list(
                  list(),
                  list(),
                  list(hyper = list(prec = list(initial = 20, fixed=TRUE)))),
                bru_verbose = T)

fit_bru_aux_2 <- bru(comp_aux_2, like_obs, like_slc_share_2, like_aux_2)

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
  formula = annual_mean ~ Intercept_obs + Time_obs_1 + Time_obs_2 + 
    Random_obs_0 + Random_obs_1 + 
    Spatial_obs_0 + Spatial_obs_1 + Spatial_obs_2,
  family = "gaussian",
  data = PM10s
)

like_slc_share <- like(
  formula = R ~ Intercept_slc + Time_slc_1 + Time_slc_2 + 
    R_lag_slc + Repuls_slc + AR_slc + Spatial_slc + 
    Comp_share1 + Comp_share2,
  family = "binomial",
  Ntrials = rep(1, times = length(PM10s$R)),
  data = PM10s
)

like_aux_1 <- like(
  formula = zero ~ Random_aux1_0 + Random_aux1_1 + Comp_aux1,
  family = "gaussian",
  data = PM10s
)

like_aux_2 <- like(
  formula = zero ~ Spatial_aux2_0 + Spatial_aux2_1 + Spatial_aux2_2 + Comp_aux2,
  family = "gaussian",
  data = PM10s
)

bru_options_reset()
bru_options_set(bru_max_iter = 20,
                control.inla = list(strategy = "gaussian", int.strategy = 'eb'),#)# h = 1e-5),
                # control.mode = list(theta = c(fit_bru_obs$mode$theta, fit_bru_slc$mode$theta),
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

# Predict at grid ==================================================================================


# Posterior sample at the original sites ----------

pred_bru <- generate(fit_bru_obs, 
                     PM10s, 
                     ~ exp(Intercept_obs + Time_obs_1 + Time_obs_2 + Random_obs_0 + Random_obs_1 + 
                             Spatial_obs_0 + Spatial_obs_1 + Spatial_obs_2),
                     n.samples = 1000) %>%
  as.data.frame() %>%
  `*`(offsets)

pred_bru$year <- BS2$year
pred_bru$R <- BS2$R

# Posterior mean of sites -------------------------

annual_quantiles <- function(ann, ...){
  ann_mean <- colMeans(ann)
  qs <- quantile(ann_mean, probs=c(0.025, 0.975))
  data.frame('ann_mean' = mean(ann_mean), 
             'q_low' = qs[1],
             'q_up' = qs[2])
}

pred <- group_by(pred_bru, year) %>%
  group_modify(annual_quantiles) 

pred_1 <- filter(pred_bru, R == 1) %>%
  group_by(year) %>% 
  group_modify(annual_quantiles)

pred_0 <- filter(pred_bru, R == 0) %>%
  group_by(year) %>%
  group_modify(annual_quantiles)

# Annual mean of observations ---------------------

bs_summary <- group_by(BS2, year) %>%
  summarise(ann_mean_exp = mean(exp(bsmoke), na.rm=TRUE)*offsets)

ggplot(pred) +
  geom_line(aes(x = year, y = bs_summary$ann_mean_exp), col='black') +
  geom_line(aes(x = year, y = ann_mean), col='blue') +
  geom_ribbon(aes(x = year, ymin = q_low, ymax = q_up), fill='blue', alpha = 0.5) +
  geom_line(aes(x = year, y = pred_1$ann_mean), col='green') +
  geom_ribbon(aes(x = year, ymin = pred_1$q_low, ymax = pred_1$q_up), fill='green', alpha = 0.5) +
  geom_line(aes(x = year, y = pred_0$ann_mean), col='red') +
  geom_ribbon(aes(x = year, ymin = pred_0$q_low, ymax = pred_0$q_up), fill='red', alpha = 0.5) +
  xlab("Year") +
  ylab("BlackSmoke")


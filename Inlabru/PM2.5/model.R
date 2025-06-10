library(dplyr)
library(ggplot2)
library(sf)
library(sp)
#library(INLA)
library(inlabru)
library(reshape2)


PM25 <- readRDS("PM2.5/PM25_with_scaled_boundary.rds")
CA_boundary <- readRDS("PM2.5/CA_scaled_boundary.rds")

# Fit the Preferential sampling model ==============================================================

# Data preprocessing -------------------------------------------------------------------------------
print(sprintf("Total number of annual measurements: %s", dim(PM25)[1]))
print(sprintf("Total number of sites: %s", length(unique(PM25$site_number))))
print(sprintf("The sd of East coordinates is: %s", sd(PM25$east)))
print(sprintf("The sd of North coordinates is: %s", sd(PM25$north)))

mean_annually <- group_by(PM25, year) %>%
  summarise(mean_pm = mean(annual_mean)) 
mean_pm_annually <- mean(mean_annually$mean_pm)
PM25$annual_mean <- log(PM25$annual_mean) - log(mean_pm_annually)

no_sites = length(unique(PM25$site_number))
no_T = length(unique(PM25$year))

# Reform the data so that each site has num_of_year rows, the empty records are filled with NAs
PM25_flat <- dcast(PM25, site_number + north + east ~ year, value.var = "annual_mean")
PM25 = melt(PM25_flat, id.vars = c(1,2,3), variable.name = 'year', value.name = 'annual_mean')
PM25$year <- as.numeric(as.character(factor(PM25$year, labels = 1999:2025)))

stopifnot(dim(PM25)[1] == no_sites * no_T)

means_plot = ggplot(data = PM25, aes(x = year, y = annual_mean)) +
  geom_line(aes(group = site_number, colour = site_number)) + xlab("Year") + ylab("log(PM2.5)")
means_plot

var_annual <- group_by(PM25, year) %>% 
  summarise(var_pm = var(annual_mean, na.rm = T))
variance_plot =  ggplot(data = var_annual, aes(x = year, y = var_pm)) +
  geom_line() + geom_smooth() + xlab("Year") + ylab("Variance of log(PM2.5)") 
ggtitle('A plot of the variance of the log annual means with fitted smoother')
variance_plot

# Prepare variables for the model ==================================================================

# Rescale time to the range of [0, 1]
PM25$time <- (PM25$year - min(PM25$year)) / (max(PM25$year) - min(PM25$year))
PM25$locs <- coordinates(PM25[, c("east", "north")])

PM25 <- select(PM25, !c("east", "north"))
PM25$site_number <- as.numeric(as.factor(PM25$site_number))

# Site-selection indicator
PM25$R <- as.numeric(!is.na(PM25$annual_mean))
PM25$R_lag <- c(rep(NA, no_sites), PM25$R[1:(dim(PM25)[1]-no_sites)])

# Response variable for the auxiliary model
PM25$zero <- 0

# Compute Euclidean distances between all the sites
dists <- spDists(cbind(PM25_flat$east, PM25_flat$north))

r <- 0.1 # The maximum radius of interest to be 1 km
PM25$repulsion_ind <- 0

counter <- no_sites + 1
for (i in sort(unique(PM25$year))[-1]) {
  # First extract the data at time i
  data_i <- PM25[PM25$year == i, ]
  # Compute the repulsion indicator. Was there a site at year i-1 within radius r of it?
  PM25$repulsion_ind[counter:(counter+(no_sites - 1))] = rowSums(dists[, which(data_i$R_lag==1)] < r) > 0
  counter <- counter + no_sites
}

###################################################################################################
#                           Tune Meshgrid Size                                                    
# The size of the meshgrid length is crutial to make sure that the model fitting function "bru"   
# converges. Too small size will lead to too much unknown variables and make the fitting function 
# hard to converge with little observations. This is especially for preferential sampling model on
# the expanded dataset
###################################################################################################

# Create mesh
edge_in = 0.5 # 5 km
edge_out = 2 * edge_in
mesh = fm_mesh_2d_inla(loc = cbind(PM25$east, PM25$north),
                       boundary = SOCAB_bord,
                       offset = c(2*edge_in, edge_out), 
                       max.edge = 2*c(edge_in, edge_out),
                       cutoff = edge_in,
                       min.angle = 21)
mesh$n
ggplot(PM25_flat) + gg(mesh) + geom_point(aes(x = east, y = north), color = "blue") + 
  xlab("East / 10km") +
  ylab("North / 10km") +
  coord_fixed()

# Maybe we should consider set the PC prior using data info
spde_obj <- inla.spde2.pcmatern(mesh = mesh, 
                                alpha = 2, 
                                prior.range = c(2.5*edge_in, 0.01),
                                prior.sigma = c(1.5*sqrt(mean(var_annual$var_pm)), 0.1),
                                constr = T)

# Build inlabru model ==============================================================================

# Joint independent model --------------------------------------------------------------------------

comp_joint_indep <- ~ 
  # Components for observation model
  Intercept_obs(1) +
  Time_obs_1(time) +
  Time_obs_2(time^2) +
  Random_obs_0(site_number, model = "iid2d", n = no_sites*2, constr=TRUE) +
  Random_obs_1(site_number, weights = time, copy = "Random_obs_0") +
  Spatial_obs_0(locs, model = spde_obj) +
  Spatial_obs_1(locs, weights = time, model = spde_obj) +
  Spatial_obs_2(locs, weights = time^2, model = spde_obj) +

  # Components for site selection model
  Intercept_slc(1) + 
  Time_slc_1(time) +
  # Time_slc_2(time^2) +   # The fitted model shows that t^2 and repuls_slc are not significant, so we revomve them and refit the model one by one
  R_lag_slc(R_lag) +
  # Repuls_slc(repulsion_ind) +
  AR_slc(year, model='ar1', hyper=list(theta1=list(prior="pcprec",param=c(2, 0.01)))) +
  Spatial_slc(locs, model = spde_obj)

like_obs <- like(
  formula = annual_mean ~ Intercept_obs + Time_obs_1 + Time_obs_2 +
    Random_obs_0 + Random_obs_1 +
    Spatial_obs_0 + Spatial_obs_1 + Spatial_obs_2,
  family = "gaussian",
  data = PM25
)

like_slc <- like(
  formula = R ~ Intercept_slc + Time_slc_1 + #Time_slc_2 +
    # R_lag_slc + Repuls_slc + AR_slc + Spatial_slc,
    R_lag_slc + AR_slc + Spatial_slc,
  family = "binomial",
  Ntrials = rep(1, times = length(PM25$R)),
  data = PM25
)

bru_options_reset()
bru_options_set(bru_max_iter = 10,
                control.inla = list(strategy = "gaussian", int.strategy = 'eb'),
                control.predictor=list(link=1),
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
  Intercept_slc(1) + 
  Time_slc_1(time) +
  # Time_slc_2(time^2) +
  R_lag_slc(R_lag) +
  # Repuls_slc(repulsion_ind) +
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
  data = PM25
)

like_slc_share <- like(
  formula = R ~ Intercept_slc + Time_slc_1 + # Time_slc_2 + 
    # R_lag_slc + Repuls_slc + AR_slc + Spatial_slc + 
    R_lag_slc + AR_slc + Spatial_slc +
    Comp_share1 + Comp_share2,
  family = "binomial",
  Ntrials = rep(1, times = length(PM25$R)),
  data = PM25
)

like_aux_1 <- like(
  formula = zero ~ Random_aux1_0 + Random_aux1_1 + Comp_aux1,
  family = "gaussian",
  data = PM25
)

like_aux_2 <- like(
  formula = zero ~ Spatial_aux2_0 + Spatial_aux2_1 + Spatial_aux2_2 + Comp_aux2,
  family = "gaussian",
  data = PM25
)

bru_options_reset()
bru_options_set(bru_max_iter = 20,
                control.inla = list(strategy = "gaussian", int.strategy = 'eb'),
                control.predictor=list(link=1),
                control.family = list(
                  list(),
                  list(),
                  list(hyper = list(prec = list(initial = 20, fixed=TRUE))),
                  list(hyper = list(prec = list(initial = 20, fixed=TRUE)))),
                bru_verbose = T)

start_time_aux <- Sys.time()
fit_bru_aux <- bru(comp_aux, like_obs, like_slc_share, like_aux_1, like_aux_2)
end_time_aux <- Sys.time()
runtime_aux <- end_time_aux - start_time_aux


# The model with extended data set =================================================================

# Once we have created the mesh, we expand the data and treat all mesh nodes inside the border as pseudo sites
PM25_flat <- dcast(PM250, site_number + north + east ~ year, value.var = "annual_mean")
nodes_locs <- data.frame("N" = mesh$loc[, 2]*10, "E" = mesh$loc[, 1]*10) %>%
  st_as_sf(coords = c("E", "N"), crs = crs_utm_km)
idx_in <- st_contains(SOCAB_border, nodes_locs) %>% as.matrix() %>% c()

sites_locs <- data.frame("east" = PM25_flat$east, "north" = PM25_flat$north)
nodes_in_locs <- data.frame("east" = mesh$loc[idx_in, 1], "north" = mesh$loc[idx_in, 2])
pseudo_sites <- setdiff(nodes_in_locs, sites_locs) 

# Number of pseudo sites is a little bit larger than no_nodes_in - no_sites because there are sites too close
# and these sites are not on top of mesh nodes.
PM25_flat_pseudo <- data.frame(matrix(data=NA, nrow=dim(pseudo_sites)[1], ncol=dim(PM25_flat)[2])) %>%
  `colnames<-`(colnames(PM25_flat))
PM25_flat_pseudo$north <- pseudo_sites$north
PM25_flat_pseudo$east <- pseudo_sites$east
PM25_flat_expand <- rbind(PM25_flat, PM25_flat_pseudo)
# Create site number for all pseudo and real sites
PM25_flat_expand$site_number <- 1:dim(PM25_flat_expand)[1]

ggplot(PM25_flat_expand) + gg(mesh) + 
  geom_point(aes(x = east, y = north), 
             color = c(rep(2, dim(PM25_flat)[1]), rep(3, dim(PM25_flat_pseudo)[1]))) + 
  xlab("East / 10km") +
  ylab("North / 10km") + 
  coord_fixed()

# Transform it into long format 
PM25_expand = melt(PM25_flat_expand, id.vars = c(1,2,3), variable.name = 'year', value.name = 'annual_mean')
PM25_expand$year <- as.numeric(as.character(factor(PM25_expand$year, labels = 1985:2022)))

no_sites_expand = dim(PM25_flat_expand)[1]
stopifnot(dim(PM25_expand)[1] == no_sites_expand * no_T)

PM25_expand$time <- (PM25_expand$year - min(PM25_expand$year)) / (max(PM25_expand$year) - min(PM25_expand$year))
PM25_expand$locs <- coordinates(PM25_expand[, c("east", "north")])
PM25_expand <- select(PM25_expand, !c("east", "north"))
PM25_expand$site_number <- as.numeric(as.factor(PM25_expand$site_number))

# Site-selection indicator
PM25_expand$R <- as.numeric(!is.na(PM25_expand$annual_mean))
PM25_expand$R_lag <- c(rep(NA, no_sites_expand), PM25_expand$R[1:(dim(PM25_expand)[1]-no_sites_expand)])

# Response variable for the auxiliary model
PM25_expand$zero <- 0

# Compute Euclidean distances between all the sites
dists_expand <- spDists(cbind(PM25_flat_expand$east, PM25_flat_expand$north))

r <- 0.1 # The maximum radius of interest to be 1 km
PM25_expand$repulsion_ind <- 0

counter <- no_sites_expand + 1
for (i in sort(unique(PM25_expand$year))[-1]) {
  # First extract the data at time i
  data_i <- PM25_expand[PM25_expand$year == i, ]
  # Compute the repulsion indicator. Was there a site at year i-1 within radius r of it?
  PM25_expand$repulsion_ind[counter:(counter+(no_sites_expand - 1))] = rowSums(dists_expand[, which(data_i$R_lag==1)] < r) > 0
  counter <- counter + no_sites_expand
}


# Joint preferential sampling model ----------------------------------------------------------------

comp_aux_expand <- ~ 
  # Components for observation model
  Intercept_obs(1) + 
  Time_obs_1(time) +
  Time_obs_2(time^2) +
  Random_obs_0(site_number, model = "iid2d", n = no_sites_expand*2, constr=TRUE) +
  Random_obs_1(site_number, weights = time, copy = "Random_obs_0") +
  Spatial_obs_0(locs, model = spde_obj) +
  Spatial_obs_1(locs, weights = time, model = spde_obj) +
  Spatial_obs_2(locs, weights = time^2, model = spde_obj) +
  
  # Components for site selection model
  Intercept_slc(1) + Time_slc_1(time) +
  # Time_slc_2(time^2) +
  R_lag_slc(R_lag) +
  # Repuls_slc(repulsion_ind) +
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
  data = PM25_expand
)

like_slc_share <- like(
  formula = R ~ Intercept_slc + Time_slc_1 + #Time_slc_2 + 
    R_lag_slc + AR_slc + Spatial_slc + 
    # R_lag_slc + Repuls_slc + AR_slc + Spatial_slc + 
    Comp_share1 + Comp_share2,
  family = "binomial",
  Ntrials = rep(1, times = length(PM25_expand$R)),
  data = PM25_expand
)

like_aux_1 <- like(
  formula = zero ~ Random_aux1_0 + Random_aux1_1 + Comp_aux1,
  family = "gaussian",
  data = PM25_expand
)

like_aux_2 <- like(
  formula = zero ~ Spatial_aux2_0 + Spatial_aux2_1 + Spatial_aux2_2 + Comp_aux2,
  family = "gaussian",
  data = PM25_expand
)

bru_options_reset()
bru_options_set(bru_max_iter = 20,
                control.inla = list(strategy = "gaussian", int.strategy = 'eb'),
                control.predictor=list(link=1),
                control.family = list(
                  list(),
                  list(),
                  list(hyper = list(prec = list(initial = 20, fixed=TRUE))),
                  list(hyper = list(prec = list(initial = 20, fixed=TRUE)))),
                bru_verbose = T)

start_time_aux_expand <- Sys.time()
fit_bru_aux_expand <- bru(comp_aux_expand, like_obs, like_slc_share, like_aux_1, like_aux_2)
end_time_aux_expand <- Sys.time()
runtime_aux_expand <- end_time_aux_expand - start_time_aux_expand



# Predict at grid ==================================================================================
# uses the original dataset ========================================================================

# 3 different models that are fitted 
model_joint <- fit_bru_aux            # the joint model on the original dataset
model_enlarged <- fit_bru_aux_expand  # the joint model on the expanded dataset

# make prediction on the original 
pred_bru <- generate(model_joint, 
                     PM25, 
                     ~ exp(Intercept_obs + Time_obs_1 + Time_obs_2 + Random_obs_0 + Random_obs_1 + 
                           Spatial_obs_0 + Spatial_obs_1 + Spatial_obs_2),
                     n.samples = 1000) %>%
  as.data.frame() %>%
  `*`(mean_pm_annually)

pred_bru$year <- PM25$year
pred_bru$R <- PM25$R

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

bs_summary <- group_by(PM25, year) %>%
  summarise(ann_mean_exp = mean(exp(annual_mean), na.rm=TRUE)*mean_pm_annually)

ggplot(pred) +
  geom_line(aes(x = year, y = bs_summary$ann_mean_exp), col='black') +
  geom_line(aes(x = year, y = ann_mean), col='blue') +
  geom_ribbon(aes(x = year, ymin = q_low, ymax = q_up), fill='blue', alpha = 0.5) +
  geom_line(aes(x = year, y = pred_1$ann_mean), col='green') +
  geom_ribbon(aes(x = year, ymin = pred_1$q_low, ymax = pred_1$q_up), fill='green', alpha = 0.5) +
  geom_line(aes(x = year, y = pred_0$ann_mean), col='red') +
  geom_ribbon(aes(x = year, ymin = pred_0$q_low, ymax = pred_0$q_up), fill='red', alpha = 0.5) +
  xlab("Year") +
  ylab("PM25 (Micrograms/cubic meter (25 C))")


# # Plot the marginal pdf of individual effect, fixed or random.
# # To check the names of all effects:
# # names(fit_bru$marginals.fixed) or names(fit_bru$marginals.random)
# 
# plot(model_joint, "Time_obs_1")
# plot(model_joint, "Time_obs_2")
# # plot(model_joint, "Random_obs_0")
# 
# # What we are interested in is the range and variance of the Matern covariance funcion, 
# # which are functions of the parameters internally used in inlabru.
# # We can look at the posterior distributions of the range parameter and the log of the variance parameters.
# 
# spde.range <- spde.posterior(model_joint, "Spatial_obs_0", what = "range")
# spde.logvar <- spde.posterior(model_joint, "Spatial_obs_0", what = "log.variance")
# 
# range.plot <- plot(spde.range)
# var.plot <- plot(spde.logvar)
# multiplot(range.plot, var.plot)
# 
# 
# # We can look at the posterior distributions of the Matern correlation and covariance functions as follows:
# plot(spde.posterior(model_joint, "Spatial_obs_0", what = "matern.correlation"))
# plot(spde.posterior(model_joint, "Spatial_obs_0", what = "matern.covariance"))

# Predict at grid ==================================================================================
# uses the enlarged dataset ========================================================================

# make prediction on the original dataset but use the model fitted on the enlarged dataset
pred_bru_expand_initial <- generate(model_enlarged, 
                            PM25, 
                            ~ exp(Intercept_obs + Time_obs_1 + Time_obs_2 + Random_obs_0 + Random_obs_1 + 
                                    Spatial_obs_0 + Spatial_obs_1 + Spatial_obs_2),
                            n.samples = 1000) %>%
  as.data.frame() %>%
  `*`(mean_pm_annually)

pred_bru_expand_initial$year <- PM25$year
pred_bru_expand_initial$R <- PM25$R

# Posterior mean of sites -------------------------

annual_quantiles <- function(ann, ...){
  ann_mean <- colMeans(ann)
  qs <- quantile(ann_mean, probs=c(0.025, 0.975))
  data.frame('ann_mean' = mean(ann_mean), 
             'q_low' = qs[1],
             'q_up' = qs[2])
}

pred_expand_initial <- group_by(pred_bru_expand_initial, year) %>%
  group_modify(annual_quantiles) 

pred_expand_initial_1 <- filter(pred_bru_expand_initial, R == 1) %>%
  group_by(year) %>% 
  group_modify(annual_quantiles)

pred_expand_initial_0 <- filter(pred_bru_expand_initial, R == 0) %>%
  group_by(year) %>%
  group_modify(annual_quantiles)


# make prediction on the enlarged dataset -----------------------------------
pred_bru_expand <- generate(model_enlarged, 
                            PM25_expand, 
                            ~ exp(Intercept_obs + Time_obs_1 + Time_obs_2 + Random_obs_0 + Random_obs_1 + 
                                    Spatial_obs_0 + Spatial_obs_1 + Spatial_obs_2),
                            n.samples = 1000) %>%
  as.data.frame() %>%
  `*`(mean_pm_annually)

pred_bru_expand$year <- PM25_expand$year
pred_bru_expand$R <- PM25_expand$R

# Posterior mean of sites -------------------------

annual_quantiles <- function(ann, ...){
  ann_mean <- colMeans(ann)
  qs <- quantile(ann_mean, probs=c(0.025, 0.975))
  data.frame('ann_mean' = mean(ann_mean), 
             'q_low' = qs[1],
             'q_up' = qs[2])
}

pred_expand <- group_by(pred_bru_expand, year) %>%
  group_modify(annual_quantiles) 

pred_expand_1 <- filter(pred_bru_expand, R == 1) %>%
  group_by(year) %>% 
  group_modify(annual_quantiles)

pred_expand_0 <- filter(pred_bru_expand, R == 0) %>%
  group_by(year) %>%
  group_modify(annual_quantiles)

# Annual mean of observations ---------------------

bs_summary <- group_by(PM25_expand, year) %>%
  summarise(ann_mean_exp = mean(exp(annual_mean), na.rm=TRUE)*mean_pm_annually)

ggplot(pred_expand) +
  geom_line(aes(x = year, y = bs_summary$ann_mean_exp), col='black') +
  geom_line(aes(x = year, y = ann_mean), col='blue') +
  geom_ribbon(aes(x = year, ymin = q_low, ymax = q_up), fill='blue', alpha = 0.5) +
  geom_line(aes(x = year, y = pred_expand_1$ann_mean), col='green') +
  geom_ribbon(aes(x = year, ymin = pred_expand_1$q_low, ymax = pred_expand_1$q_up), fill='green', alpha = 0.5) +
  geom_line(aes(x = year, y = pred_expand_initial_0$ann_mean), col='red') +
  geom_ribbon(aes(x = year, ymin = pred_expand_initial_0$q_low, ymax = pred_expand_initial_0$q_up), fill='red', alpha = 0.5) +
  xlab("Year") +
  ylab("PM25 (Micrograms/cubic meter (25 C))")


# Plot the marginal pdf of individual effect, fixed or random.
# To check the names of all effects:
# names(fit_bru$marginals.fixed) or names(fit_bru$marginals.random)

plot(model_enlarged, "Time_obs_1")
plot(model_enlarged, "Time_obs_2")
# plot(model_enlarged, "Random_obs_0")
# 
# # What we are interested in is the range and variance of the Matern covariance funcion, 
# # which are functions of the parameters internally used in inlabru.
# # We can look at the posterior distributions of the range parameter and the log of the variance parameters.
# 
# spde.range <- spde.posterior(model_enlarged, "Spatial_obs_0", what = "range")
# spde.logvar <- spde.posterior(model_enlarged, "Spatial_obs_0", what = "log.variance")
# 
# range.plot <- plot(spde.range)
# var.plot <- plot(spde.logvar)
# multiplot(range.plot, var.plot)
# 
# 
# # We can look at the posterior distributions of the Matern correlation and covariance functions as follows:
# plot(spde.posterior(model_enlarged, "Spatial_obs_0", what = "matern.correlation"))
# plot(spde.posterior(model_enlarged, "Spatial_obs_0", what = "matern.covariance"))
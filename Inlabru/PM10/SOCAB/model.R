library(dplyr)
library(ggplot2)
library(sf)
library(sp)
library(INLA)
library(inlabru)
library(reshape2)


# The raw data downloaded from EPA
PM10s_CA_raw <- readRDS("./PM10/PM10s_CA_raw.rds")

# The data after filtering out the extreme events
PM10s_CA_selected <- dplyr::select(PM10s_CA_raw, c("site_number",
                                                   "poc",
                                                   "latitude",
                                                   "longitude",
                                                   "datum",
                                                   "year",
                                                   "event_type",
                                                   "arithmetic_mean")) %>% 
  filter(event_type %in% c("No Events",
                           "Events Excluded",
                           "Concurred Events Excluded"))

# Preprocess the data and border of SOCAB ==========================================================

# Project the sites to UTM CRS ---------------------------------------------------------------------
# Locations with CRS NAD83
PM10s_nad83 <- filter(PM10s_CA_selected, datum == "NAD83")
loc_nad83 <- dplyr::select(PM10s_nad83, c("longitude", "latitude"))

# Locations with CRS WGS84
PM10s_wgs84 <- filter(PM10s_CA_selected, datum == "WGS84")
loc_wgs84 <- dplyr::select(PM10s_wgs84, c("longitude", "latitude"))

stopifnot(
  dim(loc_nad83)[1] + dim(loc_wgs84)[1] == dim(PM10s_CA_selected)[1]
)

# Project the map into UTM zone 10 with unit kilometer
crs_utm_km <- CRS("+proj=utm +zone=10 + ellps=WGS84 +units=km")

loc_nad83_spt <- SpatialPoints(loc_nad83, proj4string=CRS("+init=epsg:4296"))
loc_nad83_to_utm <- spTransform(loc_nad83_spt, crs_utm_km) %>% coordinates()
loc_wgs84_spt <- SpatialPoints(loc_wgs84, proj4string=CRS("+init=epsg:4326"))
loc_wgs84_to_utm <- spTransform(loc_wgs84_spt, crs_utm_km) %>% coordinates()

PM10s_nad83_utm <- dplyr::mutate(PM10s_nad83, "N" = loc_nad83_to_utm[, 2],
                                 "E" = loc_nad83_to_utm[, 1])
PM10s_wgs84_utm <- dplyr::mutate(PM10s_wgs84, "N" = loc_wgs84_to_utm[, 2],
                                 "E" = loc_wgs84_to_utm[, 1])

PM10s_CA_utm <- bind_rows(PM10s_nad83_utm, PM10s_wgs84_utm) %>% 
  dplyr::select(-c("latitude", "longitude", "datum"))

# Load the border of SOCAB (South Coast Air Basin) -------------------------------------------------
CA_AirBasins <- st_read(dsn = "./Data/Air_Basins_SCAG_Region", layer = "Air_Basins_SCAG_Region")
SOCAB_border <- CA_AirBasins[CA_AirBasins$AirBasins == "South Coast Air Basin",] %>%
  st_transform(crs = crs_utm_km) %>%
  # st_polygonize() %>%
  st_union() %>%
  st_geometry()

# Filter the sites in SOCAB ------------------------------------------------------------------------
site_locs <- data.frame("N" = PM10s_CA_utm$N, "E" = PM10s_CA_utm$E) %>%
  st_as_sf(coords = c("E", "N"), crs = crs_utm_km)

xy_in <- st_contains(SOCAB_border, site_locs) %>% as.matrix() %>% c()

PM10s_SOCAB_utm <- PM10s_CA_utm[xy_in, ]

ggplot(PM10s_SOCAB_utm) + gg(as(SOCAB_border, "Spatial")) + 
  geom_point(aes(x = E, y = N), color = "blue") + coord_equal()

# Aggregate by site --------------------------------------------------------------------------------
# (Should we aggregate monitors belong to the same sites together?)
PM10s_SOCAB_summary <- group_by(PM10s_SOCAB_utm, site_number, year) %>%
  summarise(annual_mean = mean(arithmetic_mean),
            north = mean(N),
            east = mean(E)) %>%
  mutate(north = mean(north[site_number == site_number]),
         east = mean(east[site_number == site_number]))

ggplot(PM10s_SOCAB_summary) + gg(as(SOCAB_border, "Spatial")) + 
  geom_point(aes(x = east, y = north), color = "blue") + coord_equal()

# Rescale the coordinates of sites and border, the new unit is 10km --------------------------------

sd_north <- sd(PM10s_SOCAB_summary$north) # = 20.45
sd_east <- sd(PM10s_SOCAB_summary$east) # = 41.68

PM10s_SOCAB_scaled <- PM10s_SOCAB_summary
PM10s_SOCAB_scaled$north <- PM10s_SOCAB_summary$north / 10
PM10s_SOCAB_scaled$east <- PM10s_SOCAB_summary$east / 10

SOCAB_border_scaled <- (SOCAB_border / matrix(data = c(10, 10), ncol = 2)) %>%
  as("Spatial")

ggplot(PM10s_SOCAB_scaled) + gg(SOCAB_border_scaled) + 
  geom_point(aes(x = east, y = north), color = "blue") + coord_equal()


# Fit the Preferential sampling model ==============================================================

# Data preprocessing -------------------------------------------------------------------------------

PM10s0 <- PM10s_SOCAB_scaled
SOCAB_bord <- SOCAB_border_scaled

print(sprintf("Total number of annual measurements: %s", dim(PM10s0)[1]))
print(sprintf("Total number of sites: %s", length(unique(PM10s0$site_number))))
print(sprintf("The sd of East coordinates is: %s", sd(PM10s0$east)))
print(sprintf("The sd of North coordinates is: %s", sd(PM10s0$north)))

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

means_plot = ggplot(data = PM10s, aes(x = year, y = annual_mean)) +
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
PM10s <- select(PM10s, !c("east", "north"))
PM10s$site_number <- as.numeric(as.factor(PM10s$site_number))

# Site-selection indicator
PM10s$R <- as.numeric(!is.na(PM10s$annual_mean))
PM10s$R_lag <- c(rep(NA, no_sites), PM10s$R[1:(dim(PM10s)[1]-no_sites)])

# Response variable for the auxiliary model
PM10s$zero <- 0

# Compute Euclidean distances between all the sites
dists <- spDists(cbind(PM10s_flat$east, PM10s_flat$north))

r <- 0.1 # The maximum radius of interest to be 1 km
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
edge_in = 0.3 # 1.5km
edge_out = 2 * edge_in
mesh = fm_mesh_2d_inla(loc = cbind(PM10s$east, PM10s$north),
                       boundary = SOCAB_bord,
                       offset = c(2*edge_in, edge_out), 
                       max.edge = 2*c(edge_in, edge_out),
                       cutoff = edge_in,
                       min.angle = 21)
mesh$n
ggplot(PM10s_flat) + gg(mesh) + geom_point(aes(x = east, y = north), color = "blue") + coord_fixed()

# Maybe we should consider set the PC prior using data info
spde_obj <- inla.spde2.pcmatern(mesh = mesh, 
                                alpha = 2, 
                                prior.range = c(1.5*edge_in, 0.01),
                                prior.sigma = c(1.5*sqrt(mean(var_annual$var_pm)), 0.1),
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
                control.inla = list(strategy = "gaussian", int.strategy = 'eb'),
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

pred_bru <- generate(fit_bru_aux, 
                     PM10s, 
                     ~ exp(Intercept_obs + Time_obs_1 + Time_obs_2 + Random_obs_0 + Random_obs_1 + 
                             Spatial_obs_0 + Spatial_obs_1 + Spatial_obs_2),
                     n.samples = 2000) %>%
  as.data.frame() %>%
  `*`(mean_pm_annually)

pred_bru$year <- PM10s$year
pred_bru$R <- PM10s$R

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

bs_summary <- group_by(PM10s, year) %>%
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
  ylab("PM10s")







# The model with extended data set =================================================================

# Once we have created the mesh, we expand the data and treat all mesh nodes inside the border as pseudo sites
PM10s_flat <- dcast(PM10s0, site_number + north + east ~ year, value.var = "annual_mean")
nodes_locs <- data.frame("N" = mesh$loc[, 2]*10, "E" = mesh$loc[, 1]*10) %>%
  st_as_sf(coords = c("E", "N"), crs = crs_utm_km)
idx_in <- st_contains(SOCAB_border, nodes_locs) %>% as.matrix() %>% c()

sites_locs <- data.frame("east" = PM10s_flat$east, "north" = PM10s_flat$north)
nodes_in_locs <- data.frame("east" = mesh$loc[idx_in, 1], "north" = mesh$loc[idx_in, 2])
pseudo_sites <- setdiff(nodes_in_locs, sites_locs) 

# Number of pseudo sites is a little bit larger than no_nodes_in - no_sites because there are sites too close
# and these sites are not on top of mesh nodes.
PM10s_flat_pseudo <- data.frame(matrix(data=NA, nrow=dim(pseudo_sites)[1], ncol=dim(PM10s_flat)[2])) %>%
  `colnames<-`(colnames(PM10s_flat))
PM10s_flat_pseudo$north <- pseudo_sites$north
PM10s_flat_pseudo$east <- pseudo_sites$east
PM10s_flat_expand <- rbind(PM10s_flat, PM10s_flat_pseudo)
# Create site number for all pseudo and real sites
PM10s_flat_expand$site_number <- 1:dim(PM10s_flat_expand)[1]

ggplot(PM10s_flat_expand) + gg(mesh) + 
  geom_point(aes(x = east, y = north), 
             color = c(rep(2, dim(PM10s_flat)[1]), rep(3, dim(PM10s_flat_pseudo)[1]))) + coord_fixed()

# Transform it into long format 
PM10s_expand = melt(PM10s_flat_expand, id.vars = c(1,2,3), variable.name = 'year', value.name = 'annual_mean')
PM10s_expand$year <- as.numeric(as.character(factor(PM10s_expand$year, labels = 1985:2022)))

no_sites_expand = dim(PM10s_flat_expand)[1]
stopifnot(dim(PM10s_expand)[1] == no_sites_expand * no_T)

PM10s_expand$time <- (PM10s_expand$year - min(PM10s_expand$year)) / (max(PM10s_expand$year) - min(PM10s_expand$year))
PM10s_expand$locs <- coordinates(PM10s_expand[, c("east", "north")])
PM10s_expand <- select(PM10s_expand, !c("east", "north"))
PM10s_expand$site_number <- as.numeric(as.factor(PM10s_expand$site_number))

# Site-selection indicator
PM10s_expand$R <- as.numeric(!is.na(PM10s_expand$annual_mean))
PM10s_expand$R_lag <- c(rep(NA, no_sites_expand), PM10s_expand$R[1:(dim(PM10s_expand)[1]-no_sites_expand)])

# Response variable for the auxiliary model
PM10s_expand$zero <- 0

# Compute Euclidean distances between all the sites
dists_expand <- spDists(cbind(PM10s_flat_expand$east, PM10s_flat_expand$north))

r <- 0.1 # The maximum radius of interest to be 1 km
PM10s_expand$repulsion_ind <- 0

counter <- no_sites_expand + 1
for (i in sort(unique(PM10s_expand$year))[-1]) {
  # First extract the data at time i
  data_i <- PM10s_expand[PM10s_expand$year == i, ]
  # Compute the repulsion indicator. Was there a site at year i-1 within radius r of it?
  PM10s_expand$repulsion_ind[counter:(counter+(no_sites_expand - 1))] = rowSums(dists_expand[, which(data_i$R_lag==1)] < r) > 0
  counter <- counter + no_sites_expand
}

# Build inlabru model ==============================================================================

# Joint independent model --------------------------------------------------------------------------

comp_joint_indep <- ~ Intercept_obs(1) + # Components for observation model
  Time_obs_1(time) +
  Time_obs_2(time^2) +
  Random_obs_0(site_number, model = "iid2d", n = no_sites_expand*2, constr=TRUE) +
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
  data = PM10s_expand
)

like_slc <- like(
  formula = R ~ Intercept_slc + Time_slc_1 + Time_slc_2 + 
    R_lag_slc + Repuls_slc + AR_slc + Spatial_slc,
  family = "binomial",
  Ntrials = rep(1, times = length(PM10s_expand$R)),
  data = PM10s_expand
)

bru_options_reset()
bru_options_set(bru_max_iter = 10,
                control.inla = list(strategy = "gaussian", int.strategy = 'eb'),
                bru_verbose = T)

fit_bru_joint_indep <- bru(comp_joint_indep, like_obs, like_slc)


# Joint preferential sampling model ----------------------------------------------------------------

comp_aux <- ~ Intercept_obs(1) + # Components for observation model
  Time_obs_1(time) +
  Time_obs_2(time^2) +
  Random_obs_0(site_number, model = "iid2d", n = no_sites_expand*2, constr=TRUE) +
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
                control.inla = list(strategy = "gaussian", int.strategy = 'eb'),
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

pred_bru <- generate(fit_bru_aux, 
                     PM10s, 
                     ~ exp(Intercept_obs + Time_obs_1 + Time_obs_2 + Random_obs_0 + Random_obs_1 + 
                             Spatial_obs_0 + Spatial_obs_1 + Spatial_obs_2),
                     n.samples = 2000) %>%
  as.data.frame() %>%
  `*`(mean_pm_annually)

pred_bru$year <- PM10s$year
pred_bru$R <- PM10s$R

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

bs_summary <- group_by(PM10s, year) %>%
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
  ylab("PM10s")





library(dplyr)
library(sf)
library(sp)
library(reshape2)
library(INLA)
library(ggplot2)
library(inlabru)



# Prepare the data =================================================================================

# Load data 
load("./Reproducibility/Data2Joe.RData")
# Load the indicator telling us if the mesh vertices lie in GB
xy_in = readRDS("./Data/Reproducibility/xy_in.rds") 

# Standardize the location
BS <- BlackSmokePrefData
sd_x <- sd(BS[, 2])
BS[, c(2,3)] <- BS[, c(2,3)] / sd_x

ggplot(BS) + geom_point(aes(x = east, y = north)) + coord_equal()

# Reshape the data with one observation per row (required by INLA)
BS2 <- melt(BS, id.vars = c(1, 2, 3), variable.name = 'year', value.name = 'bsmoke')
BS2$year = as.numeric(as.character(factor(BS2$year, labels = 66:96)))

ggplot(BS2) + geom_histogram(aes(x = bsmoke), bins = 40) 
# Right skew - take natural log and make it unitless by minus log of mean value
BS2$bsmoke = log(BS2$bsmoke / mean(colMeans(BS[, 4:34], na.rm = T)))
ggplot(BS2) + geom_histogram(aes(x = bsmoke), bins = 40) 

no_sites = length(unique(BS2$site))
no_T = length(unique(BS2$year))
stopifnot(dim(BS2)[1] == no_sites * no_T)

subsmaple_plotting_data = BS2[which(BS2$site %in% levels(BS2$site)[
  sample.int(1466, size = 30)]), ]

ggplot(data = subsmaple_plotting_data, aes(x = year, y = bsmoke)) +
  geom_line(aes(group = site, colour = site))

var_annual <- group_by(BS2, year) %>% 
  summarise(var_bs = var(bsmoke, na.rm = T))
ggplot(data = var_annual, aes(x = year, y = var_bs)) +
  geom_line() + geom_smooth() + 
  ggtitle('A plot of the variance of the log annual means with fitted smoother')


# Prepare the GB border shape file =================================================================

# Read in the map & make use it is in the same CRS: OSGB36, with reference number 27700 in sf package
uk_border <- st_read(dsn = "/Users/Xinglong/Downloads/infuse_uk_2011", layer = "infuse_uk_2011")
uk_border <- st_transform(uk_border$geometry, 27700) 

# Rescale the coordinates
bord <- (uk_border / matrix(data = c(sd_x, sd_x), ncol = 2)) %>%
  as("Spatial")

ggplot(BS2) + gg(bord) + geom_point(aes(x = east, y = north), color = "blue") + coord_equal()


# Prepare variables for the model ==================================================================

# # The mesh grid used by Joe 
# mesh_jw <- readRDS('./Data/Reproducibility/mesh_5_7.rds')
# mesh_jw$n  # equals 8695
# ggplot(BS2) + gg(mesh_jw) + geom_point(aes(x = east, y = north), color = "blue") + coord_fixed()


# Create new mesh
edg_in = 0.15 # The range set to be [0.03, 0.15]
edg_out = 4 * edg_in
mesh <- fm_mesh_2d_inla(loc = cbind(BS2$east, BS2$north),
                        boundary = bord,
                        offset = c(2*edg_in, edg_out),
                        max.edge = 4*c(edg_in, edg_out),
                        cutoff = edg_in)
mesh$n
ggplot(BS2) + gg(mesh) + geom_point(aes(x = east, y = north), color = "blue") + coord_fixed()

# Maybe we should consider set the PC prior using data infor
spde_obj <- inla.spde2.pcmatern(mesh = mesh,
                                alpha = 2,
                                prior.range = c(0.1, 0.01),
                                prior.sigma = c(1, 0.01),
                                constr = T)

BS2$time <- (BS2$year - min(BS2$year)) / (max(BS2$year) - min(BS2$year))
BS2$locs <- coordinates(BS2[, c("east", "north")])
BS2$site_number <- as.numeric(as.factor(BS2$site))

comp <- bsmoke ~ Intercept(1) + Time_1(time) + Time_2(time^2) +
  Random_0(site_number, model = "iid2d", n = no_sites*2, constr=TRUE) + 
  Random_1(site_number, weights = time, copy = "Random_0") +
  Spatial_0(locs, model = spde_obj) + 
  Spatial_1(locs, weights = time, model = spde_obj) + 
  Spatial_2(locs, weights = time^2, model = spde_obj)

# theta.ini <- fit_bru$mode$theta
# bru_options_set(control.mode = list(theta = theta.ini, restart = TRUE))
bru_options_set(bru_max_iter = 3, control.inla = list(strategy = "gaussian", int.strategy = 'eb'))

fit_bru <- bru(comp, family = "gaussian", data = BS2)
fit_bru$mode$theta

# Predict at grid ==================================================================================

pred_bru <- predict(fit_bru, 
                    BS2, 
                    ~ exp(Intercept + Time_1 + Time_2 + Random_0 + Random_1 + Spatial_0 + Spatial_1 + Spatial_2), 
                    n.samples = 1000)

pred_bru$year <- BS2$year
pred_summary <- group_by(pred_bru, year) %>%
  summarise(annual_mean = mean(mean),
            annual_sd = mean(sd))

bs_summary <- group_by(BS2, year) %>%
  summarise(annual_mean_exp = mean(exp(bsmoke), na.rm=TRUE))

ggplot(pred_summary) +
  geom_line(aes(x = year, y = annual_mean)) +
  geom_line(aes(x = year, y = bs_summary$annual_mean_exp), col='blue') +
  geom_ribbon(aes(x = year, ymin = annual_mean - annual_sd, ymax = annual_mean + annual_sd), fill = "grey70", alpha = 0.5) +
  xlab("Year") +
  ylab("BlackSmoke")


# Plot the marginal pdf of individual effect, fixed or random.
# To check the names of all effects:
# names(fit_bru$marginals.fixed) or names(fit_bru$marginals.random)

plot(fit_bru, "Intercept")
plot(fit_bru, "Time_1")
plot(fit_bru, "Time_2")

plot(fit_bru, "Random_0")

# What we are interested in is the range and variance of the Matern covariance funcion, 
# which are functions of the parameters internally used in inlabru.
# We can look at the posterior distributions of the range parameter and the log of the variance parameters.

spde.range <- spde.posterior(fit_bru, "Spatial_0", what = "range")
spde.logvar <- spde.posterior(fit_bru, "Spatial_0", what = "log.variance")

range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)
multiplot(range.plot, var.plot)


# We can look at the posterior distributions of the Matern correlation and covariance functions as follows:
plot(spde.posterior(fit_bru, "Spatial_0", what = "matern.correlation"))
plot(spde.posterior(fit_bru, "Spatial_0", what = "matern.covariance"))




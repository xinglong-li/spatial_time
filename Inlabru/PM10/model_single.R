library(dplyr)
library(sp)
library(rgdal)
library(reshape2)
library(INLA)
library(ggplot2)
library(inlabru)

path <- "./PM10/"

# Load the (scaled) locations of sites and the border of California(the unit is 100km) =============
# The data also excludes sites that contain outlier records (site_number = 0030 and 0022).

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

subsmaple_plotting_data = PM10s[which(PM10s$site_number %in% unique(PM10s$site_number)[
  sample.int(111, size = 30)]), ]
means_plot = ggplot(data = subsmaple_plotting_data, aes(x = year, y = annual_mean)) +
  geom_line(aes(group = site_number, colour = site_number)) + xlab("Year") + ylab("log(PM10)")
means_plot

var_annually <- group_by(PM10s, year) %>% 
  summarise(var_pm = var(annual_mean, na.rm = T))
variance_plot =  ggplot(data = var_annually, aes(x = year, y = var_pm)) +
  geom_line() + geom_smooth() + xlab("Year") + ylab("Variance log(PM10)") 
  ggtitle('A plot of the variance of the log annual means with fitted smoother')
variance_plot

# Build the INLA mesh ==============================================================================

cutoff_dist = 0.2 # 20km
cutoff_outer = 2 * cutoff_dist

mesh = inla.mesh.2d(loc = cbind(PM10s$east, PM10s$north),
                    boundary = CA_border,
                    offset = c(0.1, 0.2), max.edge = c(cutoff_dist, cutoff_outer),
                    cutoff = c(cutoff_dist, cutoff_outer),
                    min.angle = 26)

# mesh <- fm_mesh_2d_inla(loc = cbind(PM10s$east, PM10s$north), 
#                         boundary = CA_border,
#                         offset = c(0.1, 0.2), max.edge = c(cutoff_dist, cutoff_outer),
#                         cutoff = cutoff_dist,
#                         min.angle = 26)

ggplot(PM10s) + gg(mesh) + geom_point(aes(x = east, y = north)) + coord_fixed()

mesh$n

spde_obj <- inla.spde2.pcmatern(mesh = mesh, 
                                alpha = 2, 
                                prior.range = c(0.04, 0.05),
                                prior.sigma = c(1, 0.01),
                                constr = T)
# alpha = 2 implies 1st order smoothness, i.e., 1 times differentiable
# PC prior says we believe the lower 1st percentile of range if 3.4 km?
# 99th percentile for the GRF's standard deviation is 1. We don't believe sd higher.


# Build inlabru model ==============================================================================

# comp <- ~ Intercept(1) + Time(time) + Time_square(I(time^2)) + 
#   Random_0(site_number, model = "iid2d") + 
#   Random_1(site_number, weights = time, copy = Random_Intercept) + 
#   Spatial_0(mesh, model = spde_obj) + 
#   Spatial_1(mesh, weights = time, model = spde_obj) + 
#   Spatial_2(mesh, weights = I(time^2), model = spde_obj)

time <- (PM10s$year - min(PM10s$year)) / (max(PM10s$year) - min(PM10s$year))

km_proj <- CRS("+proj=utm +zone=10 + ellps=WGS84 +units=km")
locs <- SpatialPoints(PM10s[, c("east", "north")], km_proj)

annual_mean <- PM10s$annual_mean

site_number <- as.numeric(as.factor(PM10s$site_number))

N <- dim(PM10s)[1]

# comp <- annual_mean ~ Intercept(1) + Time(time) + 
#   Random_0(site_number, model = "iid") +
#   Spatial_0(locs, model = spde_obj)

comp <- annual_mean ~ Intercept(1) + Time(time) +
  Random_0(site_number, time, model = "iid", constr=TRUE)

comp <- annual_mean ~ Intercept(1) + Time(time) +
  Random_0(site_number, model = "iid2d", n = no_sites*2, constr=TRUE) + 
  Random_1(no_sites + site_number, time, copy = "Random_0") +
  Spatial_0(locs, model = spde_obj) #+ Spatial_1(locs, weights = time, model = spde_obj)

fit_bru <- bru(comp, family = "gaussian")

# Predict at grid ==================================================================================

data_4pred <- list(time = time,
                   locs = locs,
                   site_number = site_number)
data_4pred <- as.data.frame(data_4pred)

pred_bru <- predict(fit_bru, 
                    data_4pred, 
                    ~ exp(Intercept + Time + Random_0 + Spatial_0), 
                    n.samples = 1000)

pred_bru$year <- PM10s$year
pred_summary <- group_by(pred_bru, year) %>%
  summarise(annual_mean = mean(mean),
            annual_sd = mean(sd))

ggplot(pred_summary) +
  geom_line(aes(x = year, y = annual_mean)) +
  geom_ribbon(aes(x = year, ymin = annual_mean - annual_sd, ymax = annual_mean + annual_sd), fill = "grey70", alpha = 0.5) +
  xlab("Year") +
  ylab("PM10")


# Plot the marginal pdf of individual effect, fixed or random.
# To check the names of all effects:
# names(fit_bru$marginals.fixed) or names(fit_bru$marginals.random)

plot(fit_bru, "Intercept")

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



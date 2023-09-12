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
print(PM10s0, n=30)
print(sprintf("Total number of annual measurements: %s", dim(PM10s0)[1]))
print(sprintf("Total number of sites: %s", length(unique(PM10s0$site_number))))
print(sprintf("The sd of East coordinates is: %s", sd(PM10s0$east)))
print(sprintf("The sd of North coordinates is: %s", sd(PM10s0$north)))

CA_border <- readRDS(sprintf("%sCA_border_scaled.rds", path))

# plot(CA_border)
# points(PM10s0$east, PM10s0$north, pch=24, col='blue', cex=0.6)
ggplot(PM10s) + gg(CA_border) + geom_point(aes(x = east, y = north)) + coord_equal()

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

# Build the INLA mesh ==============================================================================

cutoff_dist = 0.1 # 10km
mesh = inla.mesh.2d(loc = cbind(PM10s$east, PM10s$north),
                    boundary = CA_border,
                    offset = c(0.1, 0.2), max.edge = c(cutoff_dist, 0.3),
                    cutoff = c(cutoff_dist, 0.3),
                    min.angle = 26)

# plot(mesh, asp = 1)
# points(x = PM10s$east, y = PM10s$north, col = 'red')

ggplot(PM10s) + gg(mesh) + geom_point(aes(x = east, y = north)) + coord_equal()

mesh$n


matern <- inla.spde2.pcmatern(mesh = mesh, 
                              alpha = 2, 
                              prior.range = c(0.04, 0.05),
                              prior.sigma = c(1, 0.01),
                              constr = T)
# alpha = 2 implies 1st order smoothness, i.e., 1 times differentialble
# PC prior says we believe the lower 1st percentile of range if 3.4 km ?
# 99th percentile for the GRF's standard deviation is 1. We don't believe sd higher.


time <- (1:no_T) / no_T

comp <- Intercept(1) + Time(time) + Time_square(I(time^2)) + 
  Random_0(site_number, model = "iid2d") + 
  Random_1(site_number, weights = time, copy = Random_Intercept) + 
  Spatial_0(mesh, model = spde_obj) + 
  Spatial_1(mesh, weights = time, model = spde_obj) + 
  Spatial_2(mesh, weights = I(time^2), model = spde_obj)

formula_obs <- annual_mean ~ 

like_obs <- like(formula = formula_obs,
                 data = PM10s,
                 family = "gaussian"
)

fit_bru <- bru(comp,
               like_obs,
               options = list(control.inla = list(int.strategy = "eb"),
                              bru_max_iter = 1)
)

# Predict at grid ==================================================================================

grid_4pred <- data.frame(x = xs)

pred_bru <- predict(fit, 
                    grid_4pred, 
                    x ~ exp(field + Intercept), 
                    n.samples = 1000)

ggplot(PM10s) +
  gg(pred_bru) +
  geom_point(data = cd, aes(x = x, y = count / exposure), cex = 2) +
  geom_point(data = true.lambda, aes(x, y), pch = "_", cex = 9, col = "blue") +
  coord_cartesian(xlim = c(0, 55), ylim = c(0, 6)) +
  xlab("x") +
  ylab("Intensity")


# Plot the marginal pdf of individual effect, fixed or random.
# To check the names of all effects:
# names(fit_bru$marginals.fixed) or names(fit_bru$marginals.random)

plot(fit_bru, "Intercept")

# What we are interested in is the range and variance of the Matern covariance funcion, 
# which are functions of the parameters internally used in inlabru.
# We can look at the posterior distributions of the range parameter and the log of the variance parameters.

spde.range <- spde.posterior(fit2.bru, "Spatial_0", what = "range")
spde.logvar <- spde.posterior(fit2.bru, "Spatial_0", what = "log.variance")

range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)
multiplot(range.plot, var.plot)


# We can look at the posterior distributions of the Matern correlation and covariance functions as follows:
plot(spde.posterior(fit2.bru, "field", what = "matern.correlation"))
plot(spde.posterior(fit2.bru, "field", what = "matern.covariance"))



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

cutoff_dist = 0.3 # 20km
cutoff_outer = 2 * cutoff_dist

mesh = inla.mesh.2d(loc = cbind(PM10s$east, PM10s$north),
                    boundary = CA_border,
                    offset = c(0.1, 0.2), max.edge = c(cutoff_dist, cutoff_outer),
                    cutoff = c(cutoff_dist, cutoff_outer),
                    min.angle = 26)

ggplot(PM10s) + gg(mesh) + geom_point(aes(x = east, y = north)) + coord_fixed()

mesh$n

# Maybe we should consider set the PC prior using data infor
spde_obj <- inla.spde2.pcmatern(mesh = mesh, 
                                alpha = 2, 
                                prior.range = c(0.1, 0.01),
                                prior.sigma = c(1, 0.01),
                                constr = T)

# Build inlabru model ==============================================================================

PM10s$time <- (PM10s$year - min(PM10s$year)) / (max(PM10s$year) - min(PM10s$year))
PM10s$locs <- coordinates(PM10s[, c("east", "north")])
PM10s$site_number <- as.numeric(as.factor(PM10s$site_number))

# When adding 'weights' in building component to introduce random slopes,
# the model can be fitted but the summary function does not work and pops up error message?

comp <- annual_mean ~ Intercept(1) + Time_1(time) + Time_2(time^2) +
  Random_0(site_number, model = "iid2d", n = no_sites*2, constr=TRUE) + 
  Random_1(site_number, weights = time, copy = "Random_0") +
  Spatial_0(locs, model = spde_obj) + 
  Spatial_1(locs, weights = time, model = spde_obj) + 
  Spatial_2(locs, weights = time^2, model = spde_obj)

# theta.ini <- fit_bru$mode$theta
# bru_options_set(control.mode = list(theta = theta.ini, restart = TRUE))
bru_options_set(bru_max_iter = 5)

fit_bru <- bru(comp, family = "gaussian", data = PM10s)
fit_bru$mode$theta

# Predict at grid ==================================================================================

pred_bru <- predict(fit_bru, 
                    PM10s, 
                    ~ exp(Intercept + Time_1 + Time_2 + Random_0 + Random_1 + Spatial_0 + Spatial_1 + Spatial_2), 
                    n.samples = 1000)

pred_bru$year <- PM10s$year
pred_summary <- group_by(pred_bru, year) %>%
  summarise(annual_mean = mean(mean),
            annual_sd = mean(sd))

pm_summary <- group_by(PM10s, year) %>%
  summarise(annual_mean_exp = mean(exp(annual_mean), na.rm=TRUE))

ggplot(pred_summary) +
  geom_line(aes(x = year, y = annual_mean)) +
  geom_line(aes(x = year, y = pm_summary$annual_mean_exp), col='blue') +
  geom_ribbon(aes(x = year, ymin = annual_mean - annual_sd, ymax = annual_mean + annual_sd), fill = "grey70", alpha = 0.5) +
  xlab("Year") +
  ylab("PM10")


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



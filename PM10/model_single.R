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

# Build the INLA mesh =============================================================================

cutoff_dist = 0.2 # 20km
cutoff_outer = 2 * cutoff_dist

mesh = inla.mesh.2d(loc = cbind(PM10s$east, PM10s$north),
                    boundary = CA_border,
                    offset = c(0.1, 0.2), max.edge = c(cutoff_dist, cutoff_outer),
                    cutoff = c(cutoff_dist, cutoff_outer),
                    min.angle = 26)

ggplot(PM10s) + gg(mesh) + geom_point(aes(x = east, y = north)) + coord_fixed()

mesh$n

# Create the Matern SPDE object ====================================================================

spde_obj <- inla.spde2.pcmatern(mesh = mesh, 
                                alpha = 2, 
                                prior.range = c(0.04, 0.05),
                                prior.sigma = c(1, 0.01),
                                constr = T)
# alpha = 2 implies 1st order smoothness, i.e., 1 times differentialble
# PC prior says we believe the lower 1st percentile of range if 3.4 km ?
# 99th percentile for the GRF's standard deviation is 1. We don't believe sd higher.

# Create the projector matrix ----------------------------------------------------------------------

A_proj <- inla.spde.make.A(mesh = mesh,
                           loc = as.matrix(cbind(PM10s$east, PM10s$north)),
                           group = PM10s$year - 1984,
                           n.group = no_T)

stopifnot(dim(A_proj)[1] == no_sites * no_T & dim(A_proj)[2] == mesh$n * no_T)

# Create the stack =================================================================================

s_index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = spde_obj$n.spde,
                                n.group = no_T)

s_index_copy <- s_index
names(s_index_copy) <- c("spatial.field.copy", "spatial.field.group.copy", "spatial.field.repl.copy")

s_index_copy2 <- s_index
names(s_index_copy2) <- c("spatial.field.copy2", "spatial.field.group.copy2", "spatial.field.repl.copy2")

time = (1:no_T) / no_T
time2 = time^2

cov_y <- data.frame(
  year = rep(time, each = no_sites),
  year_2 = rep(time2, each = no_sites),
  spatial_ind = rep(1:no_sites, times = no_T),
  spatial_ind2 = no_sites + rep(1:no_sites, times = no_T)
  )

stack_y_est <- inla.stack(
  data = list(y = PM10s$annual_mean),
  A = list(A_proj, A_proj, A_proj, 1),
  effects= list(c(s_index, list(Intercept = 1)),
                c(s_index_copy, list(Intercept_copy = 1)),
                c(s_index_copy2, list(Intercept_copy2 = 1)),
                cov_y),
  tag = 'y_est')

# The formula ======================================================================================

formula_naive <- y ~ -1 + Intercept +
  f(spatial.field, model = spde_obj) +
  f(spatial.field.copy, I(spatial.field.group.copy/no_T), model = spde_obj) +
  # f(spatial.field.copy2, I((spatial.field.group.copy2/no_T)^2), model = spde_obj) +
  I(spatial.field.group/no_T) #+ 
  # I((spatial.field.group/no_T)^2) +
  # f(spatial_ind, model = "iid2d", n = no_sites*2, constr = TRUE) +
  # f(spatial_ind2, year, copy = "spatial_ind")

# Fit the model ====================================================================================

# theta.ini = c(1.597900, -1.277423, -0.443820, -1.441220, 0.036510, -1.441336, 0.016919,
#               4.462918, 1.437147, 4, 4, 4)

out.naive = inla(formula_naive, family = 'gaussian',
                 data = inla.stack.data(stack_y_est),
                 control.predictor = list(A = inla.stack.A(stack_y_est), compute = F),
                 control.compute = list(dic=F, config = T, cpo = F),
                 control.fixed = list(mean = list(Intercept = 1, default = 0),
                                      prec = list(Intercept = 0.25, default = 0.001)),
                 # control.results = list(return.marginals.random = F,
                 #                        return.marginals.predictor = F),
                 control.inla = list(h = 0.00001, strategy = "gaussian", int.strategy = 'eb'),
                 # control.mode = list(theta = theta.ini, restart=T),
                 verbose = T, num.threads=20)

summary(out.naive)








# Project the map into UTM zone 10 with unit kilometer
km_proj <- CRS("+proj=utm +zone=10 + ellps=WGS84 +units=km")

locs <- SpatialPoints(PM10s[, c(3, 2)], km_proj) %>% coordinates()
loc_wgs84_spt <- SpatialPoints(loc_wgs84, proj4string=CRS("+init=epsg:4326"))
loc_wgs84_to_utm <- spTransform(loc_wgs84_spt, km_proj) %>% coordinates()

PM10s_nad83_utm <- dplyr::mutate(PM10s_nad83, "N" = loc_nad83_to_utm[, 2],
                                 "E" = loc_nad83_to_utm[, 1])
PM10s_wgs84_utm <- dplyr::mutate(PM10s_wgs84, "N" = loc_wgs84_to_utm[, 2],
                                 "E" = loc_wgs84_to_utm[, 1])

PM10s_utm <- bind_rows(PM10s_nad83_utm, PM10s_wgs84_utm) %>% 
  dplyr::select(-c("latitude", "longitude", "datum"))


annual_mean <- PM10s$annual_mean
comp <- annual_mean ~ Intercept(1) + Spatial_0(locs, model = spde_obj)
fit_bru <- bru(comp, 
              family = "gaussian")


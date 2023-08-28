library(dplyr)
library(sp)
library(rgdal)
library(reshape2)
library(INLA)
library(ggplot2)

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

plot(CA_border)
points(PM10s0$east, PM10s0$north, pch=24, col='blue', cex=0.6)

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

# Build the INLA mesh =============================================================================

cutoff_dist = 0.1 # 10km
mesh = inla.mesh.2d(loc = cbind(PM10s$east, PM10s$north),
                    boundary = CA_border,
                    offset = c(0.1, 0.2), max.edge = c(cutoff_dist, 0.3),
                    cutoff = c(cutoff_dist, 0.3),
                    min.angle = 26)
plot(mesh, asp = 1)
points(x = PM10s$east, y = PM10s$north, col = 'red')
mesh$n

# Create the projector matrix ----------------------------------------------------------------------

A_proj <- inla.spde.make.A(mesh = mesh,
                           loc = as.matrix(cbind(PM10s$east, PM10s$north)),
                           group = PM10s$year - 1984,
                           n.group = no_T)

stopifnot(dim(A_proj)[1] == no_sites * no_T & dim(A_proj)[2] == mesh$n * no_T)


# Create the Matern SPDE object ====================================================================

spde_obj <- inla.spde2.pcmatern(mesh = mesh, 
                                alpha = 2, 
                                prior.range = c(0.04, 0.05),
                                prior.sigma = c(1, 0.01),
                                constr = T)
# alpha = 2 implies 1st order smoothness, i.e., 1 times differentialble
# PC prior says we believe the lower 1st percentile of range if 3.4 km ?
# 99th percentile for the GRF's standard deviation is 1. We don't believe sd higher.

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
  spatial_ind2 = no_sites + rep(1:no_sites, times = no_T))

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
  f(spatial.field.copy2, I((spatial.field.group.copy2/no_T)^2), model = spde_obj) +
  I(spatial.field.group/no_T) + I((spatial.field.group/no_T)^2) +
  f(spatial_ind, model = "iid2d", n = no_sites*2, constr = TRUE) +
  f(spatial_ind2, year, copy = "spatial_ind")

# Fit the model ====================================================================================

out.naive = inla(formula_naive, family = 'gaussian',
                 data = inla.stack.data(stack_y_est),
                 control.predictor = list(A = inla.stack.A(stack_y_est), compute = F),
                 control.compute = list(dic=F, config = T, cpo = F),
                 control.fixed = list(mean = list(Intercept = 1, default = 0),
                                      prec = list(Intercept = 0.25, default = 0.001)),
                 control.results = list(return.marginals.random = F,
                                        return.marginals.predictor = F),
                 control.inla = list(h = 0.00001, strategy = "gaussian", int.strategy = 'eb'),
                 control.inla = list(strategy = "gaussian", int.strategy = 'eb'),
                 control.mode = list(theta = theta.ini, restart=T),
                 verbose = T, num.threads=20)

summary(out.naive)


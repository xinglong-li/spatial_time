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

hist(PM10s0$annual_mean, breaks=20)

mean_annually <- group_by(PM10s0, year) %>%
  summarise(mean_pm = mean(annual_mean)) 
mean_pm_annually <- mean(mean_annually$mean_pm)
PM10s0$annual_mean <- log(PM10s0$annual_mean) - log(mean_pm_annually)
hist(PM10s0$annual_mean, breaks=20)

subsmaple_plotting_data = PM10s0[which(PM10s0$site_number %in% unique(PM10s0$site_number)[
  sample.int(111, size = 30)]),]
means_plot = ggplot(data = subsmaple_plotting_data, aes(x = year, y = annual_mean)) +
  geom_line(aes(group = site_number, colour = site_number))
means_plot

var_annually <- group_by(PM10s0, year) %>% # Here the annual_mean of PM10s are already in standardized log scale
  summarise(var_pm = var(annual_mean))
variance_plot =  ggplot(data = var_annually, aes(x = year, y = var_pm)) +
  geom_line() + geom_smooth() + ggtitle('A plot of the variance of the log annual means with fitted smoother')
variance_plot

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

mesh_loc <- data.frame("north" = mesh$loc[, 1], "east" = mesh$loc[, 2])
coordinates(mesh_loc) <- ~ north + east
proj4string(mesh_loc) <- proj4string(CA_border)

xy_in <- over(mesh_loc, CA_border)$REGION
xy_in <- ! is.na(xy_in)
plot(mesh, asp = 1)
points(x = mesh_loc$north[xy_in], y = mesh_loc$east[xy_in], col = 'red')

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


# Prepare the variables (Follow Joe Watson's work) =================================================

r <- 0.1 # The maximum radius of interest to be 10 km

R <- as.numeric(!is.na(PM10s_flat[, -c(1,2,3)]))
PM10s$R <- R
PM10s$R_lag <- c(rep(NA, no_sites), PM10s$R[1:(dim(PM10s)[1]-no_sites)])

# Compute Euclidean distances between all the sites
Dists <- spDists(cbind(PM10s$east, PM10s$north))

# Observe the locations of the sites through time
PM10s$repulsion_ind <- 0
counter <- no_sites + 1
for (i in sort(unique(PM10s$year))[-1])
{
  # First extract the data at time i
  Data_i =  PM10s[PM10s$year == i,]
  
  # print(plot(x = PM10s$east[PM10s$year==i & !is.na(PM10s$annual_mean)],
  #            y = PM10s$north[PM10s$year==i & !is.na(PM10s$annual_mean)],
  #            xlab = i, xlim = range(PM10s$east), ylim = range(PM10s$north)))
  
  # Compute the repulsion indicator. Was there a site at year i-1 within radius r of it
  PM10s$repulsion_ind[counter:(counter+(no_sites - 1))] = rowSums(Dists[, which(Data_i$R_lag==1)] < r) > 0
  counter = counter + no_sites
}

# Find the mesh indices that correspond to R = 1, then we can try to predict at mesh locations elsewhere
# by pretending R = 0 to see that effect is.
# First obtain the mesh locations associated with each observation
mesh_obs_ind <- which(A_proj > 0, arr.ind = TRUE)[, 2]
number_unvisited_nodes <- mesh$n * no_T - length(unique(mesh_obs_ind))
unvisited_nodes <- which(! (1:(mesh$n*no_T) %in% unique(mesh_obs_ind)))

stopifnot(! (unvisited_nodes[1] %in% unique(mesh_obs_ind)))
stopifnot(mesh_obs_ind[1] %in% unique(mesh_obs_ind))

# Now expand the A_proj matrix by number_unvisited_nodes columns, choose weight 1 at column unvisited_nodes[1] for row i.
Ind_mat <- Matrix(0, nrow = number_unvisited_nodes, ncol = mesh$n * no_T)
for (i in 1:number_unvisited_nodes) {
  Ind_mat[i, unvisited_nodes[[i]]] = 1
}
A_proj_expand <- rbind(A_proj, Ind_mat)

# Expand the Y with NAs and R with number_unvisited_nodes 0s.
number_unvisited_nodes_per_year <- number_unvisited_nodes / no_T
data_expand <- matrix(NA, nrow = number_unvisited_nodes, ncol = dim(PM10s)[2])
data_expand = as.array(data_expand)
colnames(data_expand) <- names(PM10s)
data_expand[, 'east'] <- rep(mesh$loc[unvisited_nodes[1:number_unvisited_nodes_per_year], 1], times = no_T)
data_expand[, 'north'] <- rep(mesh$loc[unvisited_nodes[1:number_unvisited_nodes_per_year], 2], times = no_T)
data_expand[, 'repulsion_ind'] <- 0
# Put a 0 in the unvisited nodes if and only if they lie inside the GB boundary. Else put NA
data_expand[, 'R'] <- rep(ifelse(xy_in == 1, 0, NA)[unvisited_nodes[1:number_unvisited_nodes_per_year]], 
                          times = no_T)
data_expand[, 'R_lag'] <- c(rep(NA, number_unvisited_nodes_per_year), 
                            rep(0, number_unvisited_nodes_per_year*(no_T - 1)) )
# Map the repulsion onto the nodes 
data_expand[, "repulsion_ind"] <- 0
for(i in sort(unique(PM10s$year))[-1])
{
  Dists <- spDists(x = data_expand[, c('east','north')], 
                   y = PM10s[PM10s$R_lag == 1 & PM10s$year == i, c('east', 'north')])
  data_expand[, 'repulsion_ind'] <- as.numeric(rowSums(Dists < r) > 0)
}
data_expand[, "year"] <- rep(sort(unique(PM10s$year)), each = number_unvisited_nodes_per_year)
data_expand <- rbind(PM10s, data_expand)
data_expand <- as.data.frame(data_expand)


# Create the stack =================================================================================

# The stack for estimating observation process y ---------------------------------------------------

s_index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = spde_obj$n.spde,
                                n.group = no_T)

s_index_copy <- s_index
names(s_index_copy) <- c("spatial.field.copy", "spatial.field.group.copy", "spatial.field.repl.copy")

s_index_copy2 <- s_index
names(s_index_copy2) <- c("spatial.field.copy2", "spatial.field.group.copy2", "spatial.field.repl.copy2")

cov_y_expand <- data.frame(
  year = (data_expand$year - 1984) / no_T,
  year_2 = ((data_expand$year - 1984) / no_T)^2,
  spatial_ind = c(rep(1:no_sites, times = no_T), rep(NA, times = dim(data_expand)[1] - length(PM10s$R_lag))),
  spatial_ind2 = no_sites + c(rep(1:no_sites, times = no_T), 
                              rep(NA, times = dim(data_expand)[1] - length(PM10s$R_lag)))
  )

stack_y_est_expand <- inla.stack(
  data = list(y = data_expand$annual_mean,
              alldata = cbind(data_expand$annual_mean, NA, NA, NA),
              Ntrials = rep(0, times = length(data_expand$annual_mean))),
  A = list(A_proj_expand, A_proj_expand, A_proj_expand, 1),
  effects = list(c(s_index, list(Intercept = 1)),
                 c(s_index_copy, list(Incercept_copy = 1)),
                 c(s_index_copy2, list(Incercept_copy2 =1)),
                 cov_y_expand),
  tag = "y_est_expand"
  )

# The stack of R -----------------------------------------------------------------------------------

R_s_index = s_index
names(R_s_index) = c("R.spatial.filed", "R.spatial.field.group", "R.spatial.field.repl")
R_s_index$R.spatial.field = 1:(mesh$n * no_T)

R_s_index_copy = s_index
names(R_s_index_copy) = c("R.spatial.field.copy", "R.spatial.field.group.copy", "R.spatial.field.repl.copy")
R_s_index_copy2 = s_index
names(R_s_index_copy2) = c("R.spatial.field.copy2", "R.spatial.field.group.copy2", "R.spatial.field.repl.copy2")

cov_R_expand = data.frame(
  R_lag = data_expand$R_lag,
  repulsion_ind = data_expand$repulsion_ind,
  R_year = (data_expand$year - 1984) / no_T,
  yearR = (data_expand$year - 1984) / no_T,
  yearR_2 = ((data_expand$year - 1984) / no_T)^2,
  R_copyind = c(1:length(PM10s$annual_mean), rep(NA, times = dim(data_expand)[1] - length(PM10s$annual_mean)))
  )

stack_R_est_expand <- inla.stack(
  data = list(R = data_expand$R,
              alldata = cbind(NA, data_expand$R, NA, NA),
              Ntrials = rep(1, times = length(data_expand$R))),
  A = list(A_proj_expand, A_proj_expand, A_proj_expand, 1),
  effects = list(c(R_s_index, list(R.Intercept = 1)),
                 c(R_s_index_copy, list(R.Intercept = 1)),
                 c(R_s_index_copy2, list(R.Intercept = 1)),
                 cov_R_expand),
  tag = "R_est_expand"
  )

# The stack for dummy variable ---------------------------------------------------------------------

# Create dummy variable for copying linear predictor across models
PM10s$dummy = 0

s_index_dummy = s_index
names(s_index_dummy) = c("spatial.field.dummy", "spatial.field.group.dummy", "spatial.field.repl.dummy")
s_index_dummy$spatial.field.dummy = 1:(mesh$n * no_T)

s_index_copy_dummy = s_index
names(s_index_copy_dummy) = c(
  "spatial.field.copy.dummy", "spatial.field.group.copy.dummy", "spatial.field.repl.copy.dummy")
# add mesh$n 1's to the beginning as site selection is initially at time 1 #
s_index_copy_dummy$spatial.field.group.copy.dummy = c(
  rep(1, mesh$n), s_index_copy_dummy$spatial.field.group.copy.dummy[1:(mesh$n) * (no_T-1)])

s_index_copy_dummy2 = s_index
names(s_index_copy_dummy2) = c(
  "spatial.field.copy.dummy2", "spatial.field.group.copy.dummy2", "spatial.field.repl.copy.dummy2")
s_index_copy_dummy2$spatial.field.group.copy.dummy2 = c(
  rep(1, mesh$n), s_index_copy_dummy2$spatial.field.group.copy.dummy2[1:((mesh$n) * (no_T-1))])

s_index_copy_dummy3 = s_index
names(s_index_copy_dummy3) = c(
  "spatial.field.copy.dummy3", "spatial.field.group.copy.dummy3", "spatial.field.repl.copy.dummy3")
s_index_copy_dummy3$spatial.field.group.copy.dummy3 = c(
  rep(1, mesh$n), s_index_copy_dummy3$spatial.field.group.copy.dummy3[1:((mesh$n) * (no_T-1))])

cov_dummy_est = data.frame(
  spatial_ind_dummy = rep(1:no_sites, times = no_T),
  spatial_ind_dummy2 = no_sites + rep(1:no_sites, times = no_T),
  time_dummy = c(rep(sort(unique((data_expand$year - 1984) / no_T))[1], times = no_sites), 
                 rep(sort(unique((data_expand$year - 1984) / no_T))[1:(no_T - 1)], each = no_sites)),
  spatial.field.dummy2 = 1:length(PM10s$annual_mean),
  spatial.field.repl.dummy2 = rep(1, times = length(PM10s$annual_mean))
  )

stack_dummy_est = inla.stack(
  data = list(dummy2 = rep(0, times = length(PM10s$annual_mean)),
              alldata = cbind(rep(NA, times = length(PM10s$annual_mean)),
                              rep(NA, times = length(PM10s$annual_mean)),
                              rep(NA, times = length(PM10s$annual_mean)),
                              rep(0, times = length(PM10s$annual_mean))),
              Ntrials = rep(0, times = length(PM10s$annual_mean))),
  A = list(1),
  effects= list(cov_dummy_est),
  tag = 'dummy_est')

# Stack for dummy expand ---------------------------------------------------------------------------

A_dummy = Diagonal(mesh$n * no_T) # Matrix(0, nrow = length(sim_data_expand$Y_obs)-(mesh$n*nyears)), ncol = mesh$n*nyears)
stack_dummy_est_expand = inla.stack(
  data = list(dummy = rep(0, times = mesh$n * no_T),#,rep(NA, times = (length(sim_data$Y_obs)-(mesh$n*nyears)) )), #single model
              alldata = cbind(rep(NA, times = mesh$n * no_T),
                              rep(NA, times = mesh$n * no_T),
                              rep(0, times = mesh$n * no_T),
                              rep(NA, times = mesh$n * no_T)),
              Ntrials = rep(0, times = mesh$n*no_T)), #joint model
  A = list(A_dummy, A_dummy, A_dummy, A_dummy),
  effects= list(c(s_index_dummy),
                c(s_index_copy_dummy),
                c(s_index_copy_dummy2),
                c(s_index_copy_dummy3)),
  tag = 'dummy_est_expand')

# The joint stack ----------------------------------------------------------------------------------

stack_joint2 <- inla.stack(stack_y_est_expand, 
                           stack_R_est_expand,
                           stack_dummy_est_expand,
                           stack_dummy_est)


# The formula ======================================================================================

formula_joint2 = alldata ~ -1 + Intercept +
  f(spatial.field, model = spde_obj) + #spatial.field.group) +
  f(R.spatial.field, copy = "spatial.field.dummy", fixed = F) +
  f(R.spatial.field.copy, model = spde_obj) +
  f(R_copyind, copy = "spatial.field.dummy2", fixed = F) +
  f(spatial.field.copy, I(spatial.field.group.copy / no_T), model = spde_obj) +
  f(spatial.field.copy2, I((spatial.field.group.copy2 / no_T)^2), model = spde_obj) +
  R_lag + year + year_2 + yearR + yearR_2 + repulsion_ind +
  f(R_year, model = "ar1", hyper = list(theta1 = list(prior = "pcprec", param=c(2, 0.01)))) +
  I(R_year == 1 / no_T) +
  f(spatial.field.dummy, I(-spatial.field.repl.dummy), model = "iid", 
    hyper = list(prec = list(initial = -20, fixed = TRUE))) +
  f(spatial.field.dummy2, I(-spatial.field.repl.dummy2), model = "iid", 
    hyper = list(prec = list(initial = -20, fixed = TRUE))) +
  f(spatial.field.copy.dummy, copy = "spatial.field", fixed = T) + #, replicate = I(pmax(1,spatial.field.group.copy.dummy-1))
  f(spatial.field.copy.dummy2, I(spatial.field.group.copy.dummy2 / no_T), copy = "spatial.field.copy", fixed = T) +
  f(spatial.field.copy.dummy3, I((spatial.field.group.copy.dummy3 / no_T)^2), copy = "spatial.field.copy2", fixed = T) +
  f(spatial_ind, model = "iid2d", n = no_sites * 2, constr = TRUE) + #random site-specific intercepts
  f(spatial_ind2, year, copy = "spatial_ind") + # random site-specific slopes
  f(spatial_ind_dummy, copy = "spatial_ind", fixed = T) +
  f(spatial_ind_dummy2, time_dummy, copy = 'spatial_ind', fixed = T)


# Fit the model ====================================================================================

theta.ini2 <- c(2.95712439, -0.39504940, -0.02055101, -2.30589883, -2.52425155, -2.22719145, -0.52338512,
                -2.22312867, -0.38897875,  0.93644135,  0.63899640,  0.51169121, -0.10588640, -2.10631419,
                2.77008655,  0.11748417)

out.joint2 <- inla(formula_joint2,
                   family = c("gaussian", "binomial", "gaussian", "gaussian"),
                   data = inla.stack.data(stack_joint2),
                   Ntrials = inla.stack.data(stack_joint2)$Ntrials,
                   control.predictor = list(A = inla.stack.A(stack_joint2), compute = F),
                   control.compute = list(dic = F, config = T, cpo = F),
                   control.fixed = list(mean = list(Intercept = 1, default = 0),
                                        prec = list(Intercept = 0.25, default = 0.001)),
                   # control.results = list(return.marginals.random = F,
                   #                        return.marginals.predictor = F),
                   control.inla = list(h = 0.00001, int.strategy = "eb"),
                   # control.mode = list(theta = theta.ini2, restart = T),
                   control.family = list(list(),
                                         list(),
                                         list(hyper = list(prec = list(initial = 20, fixed = TRUE))),
                                         list(hyper = list(prec = list(initial = 20, fixed = TRUE)))),
                   verbose = TRUE,
                   num.threads = 16)

summary(out.joint2)

saveRDS(out.joint2, sprintf("%sjoint_model2.rds", path))

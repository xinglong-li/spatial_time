library(INLA)
library(reshape2)
library(ggplot2)
library(sp)


################################
#       Data pre-process       #
################################

# Load the Black Smoke Pref Data 
BS <- get(load("./Data2Joe.RData"))

# The columns are: site, east, north (both in British National Grid),
# black smoke by year (mu gm/3) where NA means it wasn’t measured.
sd_x <- sd(BS$east)
sd_y <- sd(BS$north)

BS[, c(2,3)] <- BS[, c(2,3)] / sd_x # standardize x,y coords (dividing by 86021)
plot(x = BS$east, y = BS$north)

# Reshape the data with one observation per row (required by INLA)
BS2 <- melt(BS, id.vars = c(1,2,3), variable.name = 'year', value.name = 'bsmoke')
BS2$year <- as.numeric(as.character(factor(BS2$year, labels = 66:96)))
hist(BS2$bsmoke) #right skew 

# Take natural log of measurement
BS2$bsmoke <- log(BS2$bsmoke / mean(colMeans(BS[, 4:34], na.rm = T)))
hist(BS2$bsmoke) # more bell-shaped

# Plot samples of data
subsmaple_plot_data <- BS2[which(BS2$site %in% levels(BS2$site)[sample.int(1466, size = 30)]), ]
means_plot <- ggplot(data = subsmaple_plot_data, aes(x = year, y = bsmoke)
                     ) + geom_line(aes(group = site, colour = site))
means_plot

subsample_plot_data2 <- data.frame(variance = 0, year = 66:96)
subsample_plot_data2$variance <- apply(log(BS[,4:34] / 30.7), 2 , var, na.rm=T )
variance_plot <-  ggplot(data = subsample_plot_data2, aes(x = year, y = variance)
                         ) + ggtitle('A plot of the variance of the log annual means with fitted smoother'
                         ) + geom_line() + geom_smooth()
variance_plot

# Extract selection variable (1 if site selected at time t)
BS2$R <- as.numeric(!is.na(BS2$bsmoke))


############################
#     Prepare for INLA     #
############################

# Grid for projection
ncol <- 100 
nrow <- 100
num_of_grid <- nrow * ncol
east_grid <- seq(from = min(BS$east, na.rm=T), to = max(BS$east, na.rm=T), length.out = ncol)
north_grid <- seq(from = min(BS$north, na.rm=T), to = max(BS$north, na.rm=T), length.out = nrow)

# Form the regular mesh
# Create a rough convex boundary for the UK
# Form the grid independently from the sites to avoid preferential grid selection
UK_domain <- cbind(c(2, 7.7, 7.7, 6, 4, 1), c(0.5, 0.5, 6, 13.5, 13.5, 12))
hull <- inla.nonconvex.hull(cbind(BS2$east, BS2$north))

# Set cutoff distance (minimum edge) and max.edge to less than the process range
# 17km was found to be the smallest range for the spatial fields fit in Shaddick and Zidek (2014)
# Set max.edge and cutoff to be 10km
cutoff_dist = 16000 / sd_x # 10km < 17km as min edge
max.edge = 100000 / sd_x # 100km max edge as in Shaddick and Zidek

mesh <- inla.mesh.2d(loc = cbind(BS2$east, BS2$north),
                     boundary = hull,
                     offset = c(0.1, 0.2), max.edge = c(cutoff_dist, 0.5),
                     cutoff = c(cutoff_dist, 0.5),
                     min.angle = 26)

plot(mesh, asp = 1)
points(x = BS2$east, y = BS2$north, col = 'red')
mesh$n

# Compute Euclidean distances between all the sites #
Dists <- spDists(cbind(BS$east, BS$north))

# Set the maximum radius of interest
r <- 0.116 # Note that this is ~ 10km

# Compute R_lag - the indicator for if the site was selected at time t-1
num_of_sites <- as.numeric(length(unique(BS2$site)))
num_of_T <- as.numeric(length(unique(BS2$year)))
BS2$R_lag <- c(rep(NA, num_of_sites), BS2$R[1:(dim(BS2)[1]-num_of_sites)])

# Observe the locations of the sites through time
BS2$repulsion_ind <- 0
counter <- num_of_sites + 1
for (i in sort(unique(BS2$year))[-1]){
  # First extract the data at time i
  Data_i <-  BS2[BS2$year == i,]
  
  print(plot(x = BS2$east[BS2$year==i & !is.na(BS2$bsmoke)],
             y = BS2$north[BS2$year==i & !is.na(BS2$bsmoke)],
             xlab = i, xlim = range(BS2$east), ylim = range(BS2$north)))
  
  # Compute the repulsion indicator. Was there a site at year i-1 within radius r of it
  BS2$repulsion_ind[counter:(counter+(num_of_sites-1))] <- rowSums(Dists[,which(Data_i$R_lag==1)] < r) > 0
  counter <- counter + num_of_sites
}


##############################
#        Build Model         #
##############################

######################################################
#        Create the Matern spde object for Y_grf     #
######################################################

spde_obj <- inla.spde2.pcmatern(mesh = mesh, alpha = 2,
                                prior.range = c(0.04, 0.05),
                                prior.sigma = c(1, 0.01),
                                constr = T) 
# alpha = 2 implies 1st order smoothness i.e. 1 times differentiable (same as exponential)
# PC prior says we believe the lower 1st percentile of range is 3.4km 
# (a fifth of the min range found by Shaddick and Zidek)
# 99th percentile for the GRF's standard deviation is 1. We don't believe sd higher.

# Create the projector matrix
# nrows = number of sites * number of time points
# ncols = no.nodes * number of time points
A_proj <- inla.spde.make.A(mesh = mesh,
                           loc = as.matrix(cbind(BS2$east, BS2$north)),
                           group = BS2$year-65, # group membership needs to be 1:num_of_T
                           n.group = num_of_T)

s_index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = spde_obj$n.spde,
                                n.group = num_of_T)


#######################
#     Naive model     #
#######################

# Create the stack object for estimating observation process y

# Variables:
time <- (1:num_of_T) / num_of_T
time2 <- time^2
cov_y <- data.frame(year = rep(time, each = num_of_sites),
                    year_2 = rep(time2, each = num_of_sites),
                    # site-specific random intercepts
                    spatial_ind = rep(1:num_of_sites, times = num_of_T),
                    spatial_ind2 = num_of_sites + rep(1:num_of_sites, times = num_of_T))

s_index_copy <- s_index
names(s_index_copy) <- c('spatial.field.copy',"spatial.field.group.copy", "spatial.field.repl.copy")

s_index_copy2 <- s_index
names(s_index_copy2) <- c('spatial.field.copy2',"spatial.field.group.copy2", "spatial.field.repl.copy2")


stack_y_est <- inla.stack(data = list(y = BS2$bsmoke, #single model
                                      alldata = cbind(BS2$bsmoke, NA, NA),
                                      Ntrials = rep(0,times = length(BS2$bsmoke))), #joint model
                         A = list(A_proj, A_proj, A_proj, 1),
                         effects = list(c(s_index, list(Intercept = 1)),
                                        c(s_index_copy, list(Intercept_copy = 1)),
                                        c(s_index_copy2, list(Intercept_copy2 = 1)),
                                        cov_y),
                         tag = 'y_est')

# Find the MESH indices that correspond to R = 1
# Then we can try to predict at mesh locations elsewhere by pretending R = 0 to see what effect is
# First obtain the mesh location/s associated with each observation
mesh_obs_idx <- which(A_proj > 0, arr.ind = T)[, 2]
num_unvisited_nodes <- mesh$n * num_of_T - length(unique(mesh_obs_idx))
unvisited_nodes <- which( ! ((1:(mesh$n * num_of_T)) %in% unique(mesh_obs_idx)) )

# Check this is correct #
unvisited_nodes[1] %in% unique(mesh_obs_idx) # Should be false
mesh_obs_idx[1] %in% unique(mesh_obs_idx) # Should be true

# Formula for naive model
formula_naive <- y ~ -1 + Intercept +
  f(spatial.field, model = spde_obj) +
  f(spatial.field.copy, I(spatial.field.group.copy/num_of_T), model = spde_obj) +
  f(spatial.field.copy2, I((spatial.field.group.copy2/num_of_T)^2), model = spde_obj) +
  I(spatial.field.group/no_T) + I((spatial.field.group/no_T)^2) +
  f(year, model = 'ar1') +
  f(spatial_ind, model = "iid2d", n = num_of_sites*2, constr = TRUE) + #random site-specific intercepts
  f(spatial_ind2, year, copy = "spatial_ind") # random site-specific slopes
#f(spatial_ind, model = 'iid', hyper=list(theta=list(prior="loggamma",param=c(2,0.04)))) +

theta.ini = c(1.597900, -1.277423, -0.443820, -1.441220, 0.036510, -1.441336, 0.016919,
              4.462918, 1.437147, 4, 4, 4)
out.naive = inla(formula_naive, family = 'gaussian',
                 data = inla.stack.data(stack_y_est),
                 control.predictor = list(A = inla.stack.A(stack_y_est), compute = F),
                 control.compute = list(dic=F, config = T, cpo = F),
                 control.inla = list(strategy = "gaussian", int.strategy = 'eb'),
                 control.mode = list(theta = theta.ini, restart=T),
                 verbose = T, num.threads=2)
summary(out.naive)

saveRDS(out.naive,file="finalmod_naive.rds")
# rm(out.naive)


library(INLA)
library(reshape2)
library(ggplot2)
library(sp)


################################
#       Data pre-process       #
################################

# Load the Black Smoke Pref Data 
BS <- get(load("~/git_local/spacial_time/Data2Joe.RData"))

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

# Create the Matern spde object for Y_grf
spde_obj <- inla.spde2.pcmatern(mesh = mesh, alpha = 2,
                                prior.range = c(0.04, 0.05),
                                prior.sigma = c(1, 0.01),
                                constr = T) 
# alpha = 2 implies 1st order smoothness i.e. 1 times differentiable (same as exponential)
# PC prior says we believe the lower 1st percentile of range is 3.4km 
# (a fifth of the min range found by Shaddick and Zidek)
# 99th percentile for the GRF's standard deviation is 1. We don't believe sd higher.

# Create the projector matrix
A_proj <- inla.spde.make.A(mesh = mesh,
                           loc = as.matrix(cbind(BS2$east, BS2$north)),
                           group = BS2$year-65, # group membership needs to be 1:num_of_T
                           n.group = num_of_T)
dim(A_proj) 
# nrows = number of sites * number of time points, ncols = no.nodes * no time points #
no_sites * no_T
mesh$n * no_T

s_index = inla.spde.make.index(name = "spatial.field",
                               n.spde = spde_obj$n.spde,
                               n.group = no_T)

time = (1:no_T)/no_T
time2 = time^2
# create the stack object for estimating observation process y #
cov_y = data.frame(year = rep(time, each = no_sites),
                   year_2 = rep(time2, each = no_sites),
                   spatial_ind = rep(1:no_sites, times = no_T),
                   spatial_ind2 = no_sites + rep(1:no_sites, times = no_T)) # site-specific random intercepts


s_index_copy = s_index
names(s_index_copy) = c('spatial.field.copy',"spatial.field.group.copy", "spatial.field.repl.copy")

s_index_copy2 = s_index
names(s_index_copy2) = c('spatial.field.copy2',"spatial.field.group.copy2", "spatial.field.repl.copy2")


stack_y_est = inla.stack(data = list(y = BlackSmokePrefData2$bsmoke, #single model
                                     alldata = cbind(BlackSmokePrefData2$bsmoke, NA, NA),
                                     Ntrials = rep(0,times = length(BlackSmokePrefData2$bsmoke))), #joint model
                         A = list(A_proj,A_proj,A_proj,1),
                         effects= list(c(s_index, list(Intercept = 1)),
                                       c(s_index_copy, list(Intercept_copy = 1)),
                                       c(s_index_copy2, list(Intercept_copy2 = 1)),
                                       cov_y),
                         tag = 'y_est')

# Find the MESH indices that correspond to R = 1 #
# Then we can try to predict at mesh locations elsewhere by pretending R = 0 to see what effect is #
# First obtain the mesh location/s associated with each observation #
mesh_obs_ind = which(A_proj > 0, arr.ind = T)[,2]
number_unvisited_nodes = mesh$n * no_T - length(unique(mesh_obs_ind))
unvisited_nodes = which( ! ((1:(mesh$n * no_T)) %in% unique(mesh_obs_ind)) )

# Check this is correct #
unvisited_nodes[1] %in% unique(mesh_obs_ind) # Should be false
mesh_obs_ind[1] %in% unique(mesh_obs_ind) # Should be true

# Now expand the A_proj matrix by number_unvisited_nodes columns - choosing weight 1 at column unvisited_nodes[i] for row i
Ind_mat = Matrix(0, nrow = number_unvisited_nodes, ncol = mesh$n * no_T)
for(i in 1:number_unvisited_nodes)
{
  Ind_mat[i,unvisited_nodes[i]]=1
}
A_proj_expand = rbind(A_proj, Ind_mat)

# Expand the Y with NAs and R with number_unvisited_nodes 0s
# number_unvisited_nodes per year
number_unvisited_nodes_per_year = number_unvisited_nodes / no_T
data_expand = matrix(NA, nrow = number_unvisited_nodes, ncol = dim(BlackSmokePrefData2)[2])
data_expand = as.array(data_expand)
colnames(data_expand) = names(BlackSmokePrefData2)
data_expand[,'east'] = rep(mesh$loc[unvisited_nodes[1:number_unvisited_nodes_per_year],1], times = no_T)
data_expand[,'north'] = rep(mesh$loc[unvisited_nodes[1:number_unvisited_nodes_per_year],2], times = no_T)
data_expand[,'repulsion_ind'] = 0
data_expand[,'R'] = 0
data_expand[,'R_lag'] = c(rep(NA, number_unvisited_nodes_per_year), rep(0, number_unvisited_nodes_per_year*(no_T-1)) )
# Map the repulsion onto the nodes #
for(i in sort(unique(BlackSmokePrefData2$year))[-1])
{
  Dists = spDists(x = data_expand[,c('east','north')], y = BlackSmokePrefData2[BlackSmokePrefData2$R_lag==1 & BlackSmokePrefData2$year == i,c('east','north')])
  data_expand[,'repulsion_ind'] = as.numeric(rowSums(Dists < r) > 0)
}

data_expand[,'year'] = rep(sort(unique(BlackSmokePrefData2$year)), each = number_unvisited_nodes_per_year)
data_expand = rbind(BlackSmokePrefData2, data_expand)
data_expand = as.data.frame(data_expand)
# Next choose only those corresponding to an observed value #
#mesh_obs_ind = mesh_obs_ind[sim_data$R_lag == 1]

# create a vector of size mesh$n * nyears with all zeros except observed indices #
#R_zero = rep(0, times = mesh$n * nyears)
#R_zero[mesh_obs_ind] = 1

# create the stack object for estimating observation process y #
cov_y_expand = data.frame(year = data_expand$year/no_T,
                          year_2 = (data_expand$year/no_T)^2,
                          spatial_ind = c(rep(1:no_sites, times = no_T),no_sites + rep(1:number_unvisited_nodes_per_year, times = no_T)), # Only estimate random effects for site locations
                          spatial_ind2 = dim(data_expand)[1] + c(rep(1:no_sites, times = no_T), no_sites+rep(1:number_unvisited_nodes_per_year, times = no_T)))

#diag_model = diag(mesh$n * nyears)
stack_y_est_expand = inla.stack(data = list(y = data_expand$bsmoke, #single model
                                            alldata = cbind(data_expand$bsmoke, NA,NA,NA),
                                            Ntrials = rep(0,times = length(data_expand$bsmoke))), #joint model
                                A = list(A_proj_expand,A_proj_expand,A_proj_expand,1),
                                effects= list(c(s_index, list(Intercept = 1)),
                                              c(s_index_copy, list(Intercept_copy = 1)),
                                              c(s_index_copy2, list(Intercept_copy2 = 1)),
                                              cov_y_expand),
                                tag = 'y_est_expand')

# create the stack object for estimating selection process R #
cov_R = data.frame(R_lag = BlackSmokePrefData2$R_lag,
                   yearR = rep(time, each = no_sites),
                   yearR_2 = rep(time2, each = no_sites),
                   repulsion_ind = BlackSmokePrefData2$repulsion_ind,
                   R_year = BlackSmokePrefData2$year)
#spatial_ind2 = sim_data$spatial_ind

R_s_index = s_index
names(R_s_index) = c('R.spatial.field',"R.spatial.field.group", "R.spatial.field.repl")
# change the R.spatial.field to 1:lengthlength(s_index_dummy$spatial.field.dummy) instead of rep(1:mesh$n, nyears)
R_s_index$R.spatial.field = 1:(mesh$n * no_T)

R_s_index_copy = s_index
names(R_s_index_copy) = c('R.spatial.field.copy',"R.spatial.field.group.copy", "R.spatial.field.repl.copy")
R_s_index_copy2 = s_index
names(R_s_index_copy2) = c('R.spatial.field.copy2',"R.spatial.field.group.copy2", "R.spatial.field.repl.copy2")
#R_s_index$R.spatial.field.group = pmax(1, R_s_index$R.spatial.field.group-1) #lag the time by 1

stack_R_est = inla.stack(data = list(R = BlackSmokePrefData2$R, #for single model
                                     alldata = cbind(NA,BlackSmokePrefData2$R, NA),
                                     Ntrials = rep(1,times = length(BlackSmokePrefData2$R))), #for joint model
                         A = list(A_proj,A_proj,A_proj,1),
                         effects= list(c(R_s_index, list(R.Intercept = 1)),
                                       c(R_s_index_copy, list(R.Intercept_copy = 1)),
                                       c(R_s_index_copy2, list(R.Intercept_copy2 = 1)),
                                       cov_R),
                         tag = 'R_est')


cov_R_expand = data.frame(R_lag = data_expand$R_lag,
                          repulsion_ind = data_expand$repulsion_ind,
                          R_year = data_expand$year,
                          yearR = data_expand$year/no_T,
                          yearR_2 = (data_expand$year/no_T)^2,
                          R_copyind = c(rep(1:no_sites, times = no_T),no_sites + rep(1:number_unvisited_nodes_per_year, times = no_T)))

stack_R_est_expand = inla.stack(data = list(R = data_expand$R, #for single model
                                            alldata = cbind(NA,data_expand$R,NA,NA),
                                            Ntrials = rep(1,times = length(data_expand$R))), #for joint model
                                A = list(A_proj_expand,A_proj_expand,A_proj_expand,1),
                                effects= list(c(R_s_index, list(R.Intercept = 1)),
                                              c(R_s_index_copy, list(R.Intercept_copy = 1)),
                                              c(R_s_index_copy2, list(R.Intercept_copy2 = 1)),
                                              cov_R_expand),
                                tag = 'R_est_expand')# Form the dummy variable with the zero observations #

# Create dummy variable for copying linear predictor across models
BlackSmokePrefData2$dummy = 0

s_index_dummy = s_index
names(s_index_dummy) = c('spatial.field.dummy',"spatial.field.group.dummy", "spatial.field.repl.dummy")
# change the spatial.field.dummy to 1:lengthlength(s_index_dummy$spatial.field.dummy) instead of rep(1:mesh$n, nyears)
s_index_dummy$spatial.field.dummy = 1:(mesh$n * no_T)

s_index_copy_dummy = s_index
names(s_index_copy_dummy) = c('spatial.field.copy.dummy',"spatial.field.group.copy.dummy", "spatial.field.repl.copy.dummy")
# add mesh$n 1's to the beginning as site selection is initially at time 1 #
s_index_copy_dummy$spatial.field.group.copy.dummy = c(rep(1,mesh$n), s_index_copy_dummy$spatial.field.group.copy.dummy[1:((mesh$n)*(no_T-1))])

s_index_copy_dummy2 = s_index
names(s_index_copy_dummy2) = c('spatial.field.copy.dummy2',"spatial.field.group.copy.dummy2", "spatial.field.repl.copy.dummy2")
s_index_copy_dummy2$spatial.field.group.copy.dummy2 = c(rep(1,mesh$n), s_index_copy_dummy2$spatial.field.group.copy.dummy2[1:((mesh$n)*(no_T-1))])

s_index_copy_dummy3 = s_index
names(s_index_copy_dummy3) = c('spatial.field.copy.dummy3',"spatial.field.group.copy.dummy3", "spatial.field.repl.copy.dummy3")
s_index_copy_dummy3$spatial.field.group.copy.dummy3 = c(rep(1,mesh$n), s_index_copy_dummy3$spatial.field.group.copy.dummy3[1:((mesh$n)*(no_T-1))])


#w = c(w, rep(-1, nyears*mesh$n))
#dummy_var = 1:(nyears*mesh$n)


cov_dummy_est = data.frame(spatial_ind_dummy = rep(1:no_sites, times = no_T),
                           spatial_ind_dummy2 = no_sites + rep(1:no_sites, times = no_T),
                           time_dummy = c(rep(1, times = no_sites), rep(1:(no_T-1), each = no_sites)))

cov_dummy2_expand = data.frame(spatial_ind_dummy = rep(1:no_sites, times = no_T),
                               spatial_ind_dummy2 = no_sites + rep(1:no_sites, times = no_T),
                               time_dummy = c(rep(1, times = no_sites), rep(1:(no_T-1), each = no_sites)),
                               spatial.field.dummy2 = 1:(no_sites*no_T),
                               spatial.field.repl.dummy2 = rep(1, times = no_sites*no_T))

# cov_dummy2_expand = data.frame(spatial_ind_dummy = c(rep(1:no_sites, times = no_T),rep(NA, times = dim(data_expand)[1] - (no_sites*no_T))), # Only estimate random effects for site locations
#                               spatial_ind_dummy2 = no_sites + c(rep(1:no_sites, times = no_T),rep(NA, times = dim(data_expand)[1] - (no_sites*no_T))),
#                               time_dummy = c(rep(1, times = no_sites), rep(1:(no_T-1), each = no_sites),rep(NA, times = dim(data_expand)[1] - (no_sites*no_T))))

stack_dummy_est = inla.stack(data = list(dummy = rep(0, times = length(BlackSmokePrefData2$bsmoke)), #single model
                                         alldata = cbind(rep(NA, times = length(BlackSmokePrefData2$bsmoke)),
                                                         rep(NA, times = length(BlackSmokePrefData2$bsmoke)),
                                                         rep(0, times = length(BlackSmokePrefData2$bsmoke))),
                                         Ntrials = rep(0,times = length(BlackSmokePrefData2$bsmoke))), #joint model
                             A = list(A_proj,A_proj,A_proj, A_proj,1),
                             effects= list(c(s_index_dummy),
                                           c(s_index_copy_dummy),
                                           c(s_index_copy_dummy2),
                                           c(s_index_copy_dummy3),
                                           cov_dummy_est),
                             tag = 'dummy_est')

# This is for copying the random intercepts and slopes across
stack_dummy2_expand = inla.stack(data = list(dummy2 = rep(0, times = length(BlackSmokePrefData2$bsmoke)), #single model
                                             alldata = cbind(rep(NA, times = length(BlackSmokePrefData2$bsmoke)),
                                                             rep(NA, times = length(BlackSmokePrefData2$bsmoke)),
                                                             rep(NA, times = length(BlackSmokePrefData2$bsmoke)),
                                                             rep(0, times = length(BlackSmokePrefData2$bsmoke))),
                                             Ntrials = rep(0,times = length(BlackSmokePrefData2$bsmoke))), #joint model
                                 A = list(1),
                                 effects= list(cov_dummy2_expand),
                                 tag = 'dummy2_expand')

A_dummy = Diagonal(mesh$n*no_T)#, Matrix(0, nrow = length(sim_data_expand$Y_obs)-(mesh$n*nyears)), ncol = mesh$n*nyears)
stack_dummy_est_expand = inla.stack(data = list(dummy = rep(0, times = mesh$n*no_T),#,rep(NA, times = (length(sim_data$Y_obs)-(mesh$n*nyears)) )), #single model
                                                alldata = cbind(rep(NA, times = mesh$n*no_T),
                                                                rep(NA, times = mesh$n*no_T),
                                                                rep(0, times = mesh$n*no_T),
                                                                rep(NA, times = mesh$n*no_T)),
                                                Ntrials = rep(0,times = mesh$n*no_T)), #joint model
                                    A = list(A_dummy,A_dummy,A_dummy,A_dummy),
                                    effects= list(c(s_index_dummy),
                                                  c(s_index_copy_dummy),
                                                  c(s_index_copy_dummy2),
                                                  c(s_index_copy_dummy3)),
                                    tag = 'dummy_est_expand')


#stack_y_combine = inla.stack(stack_y_est, stack_y_pred)
#stack_combine_est = inla.stack(stack_y_est, stack_y_pred, stack_R_est, stack_R_pred, stack_dummy_est, stack_dummy_pred)
stack_joint1 = inla.stack(stack_y_est, stack_R_est, stack_dummy_est)
stack_joint2 = inla.stack(stack_y_est_expand, stack_R_est_expand, stack_dummy_est_expand, stack_dummy2_expand)

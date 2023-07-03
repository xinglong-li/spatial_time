# R Script to re-fit the models used in Watson et al 2021
# Black smoke analysis #
library(INLA)
library(mvtnorm)
library(boot)
library(reshape2)
library(sp)
library(ggplot2)

# load data #
load("./Reproducibility/Data2Joe.RData")
# xy_in = readRDS("./Reproducibility/xy_in.rds") # load the indicator telling us if the mesh vertices lie in GB

# The columns are: site name, x, y (both in British National Grid),
# black smoke by year (mu gm/3) where NA means it wasn’t measured.
sd_x = sd(BlackSmokePrefData[,c(2)])
sd_y = sd(BlackSmokePrefData[,c(3)])

BlackSmokePrefData[,c(2,3)] = BlackSmokePrefData[,c(2,3)]/sd_x # standardize x,y coords (dividing by 86021)
plot(x = BlackSmokePrefData$east, y = BlackSmokePrefData$north)

# subset the data into 1 x 1 square for quicker calculation - optional #
#BlackSmokePrefData = BlackSmokePrefData[BlackSmokePrefData$east>0 & BlackSmokePrefData$east<0.5,]
#BlackSmokePrefData = BlackSmokePrefData[BlackSmokePrefData$north>0 & BlackSmokePrefData$north<0.5,]
#plot(x = BlackSmokePrefData$east, y = BlackSmokePrefData$north)

# reshape the data with one observation per row (required by INLA)
BlackSmokePrefData2 = melt(BlackSmokePrefData,id.vars = c(1,2,3), variable.name = 'year', value.name = 'bsmoke')
BlackSmokePrefData2$year = as.numeric(as.character(factor(BlackSmokePrefData2$year, labels =66:96 )))

hist(BlackSmokePrefData2$bsmoke) #right skew - take natural log
BlackSmokePrefData2$bsmoke = log(BlackSmokePrefData2$bsmoke / mean(colMeans(BlackSmokePrefData[,4:34], na.rm = T)))
# Divided by 30.7 first - the mean of the annual means across all sites - to make it unitless
hist(BlackSmokePrefData2$bsmoke) # more bell-shaped

subsmaple_plotting_data = BlackSmokePrefData2[which(BlackSmokePrefData2$site %in% levels(BlackSmokePrefData2$site)[sample.int(1466, size = 30)]),]

means_plot = ggplot(data = subsmaple_plotting_data, aes(x = year, y = bsmoke)) +
  geom_line(aes(group = site, colour = site))
means_plot

subsample_plotting_data2 = data.frame(variance = 0, year = 66:96)
subsample_plotting_data2$variance = apply(log(BlackSmokePrefData[,4:34] / 30.7), 2 , var, na.rm=T )
variance_plot =  ggplot(data = subsample_plotting_data2, aes(x = year, y = variance)) +
  geom_line() + geom_smooth() + ggtitle('A plot of the variance of the log annual means with fitted smoother')
variance_plot

# Set the maximum radius of interest #
r = 0.116
#Note that this is ~ 10km

# extract key information from data #
no_sites = as.numeric(length(unique(BlackSmokePrefData2$site))) # number of sites
no_T = as.numeric(length(unique(BlackSmokePrefData2$year))) # number of years
ind_select = matrix(0, nrow = no_sites, ncol = no_T)
ind_select = as.matrix(!is.na(BlackSmokePrefData[,-c(1,2,3)]))
R = as.numeric(!is.na(BlackSmokePrefData[,-c(1,2,3)]))
BlackSmokePrefData2$R = R # selection variable (1 if site selected at time t)

ncol = 100 # grid for projection
nrow = 100
L = nrow * ncol # number of grid sites
east_grid = seq(from = min(BlackSmokePrefData$east, na.rm=T),
                to = max(BlackSmokePrefData$east, na.rm=T), length.out = ncol)
north_grid = seq(from = min(BlackSmokePrefData$north, na.rm=T),
                 to = max(BlackSmokePrefData$north, na.rm=T), length.out = nrow)

# UK_domain <- cbind(c(2, 7.7, 7.7, 6, 4, 1), c(0.5, 0.5, 6, 13.5, 13.5, 12))
hull = inla.nonconvex.hull(cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north),
                           convex = -0.02)
# cutoff_dist = 16000/sd_x # 10km < 17km as min edge
cutoff_dist = 6000/sd_x # 10km < 17km as min edge
# max.edge = 100000/sd_x # 100km max edge as in Shaddick and Zidek
mesh = inla.mesh.2d(loc = cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north),
                    boundary = hull,
                    offset = c(0.1, 0.2), max.edge = c(cutoff_dist, 0.5),
                    cutoff = c(cutoff_dist, 0.5),
                    min.angle = 26)
plot(mesh, asp = 1)
points(x = BlackSmokePrefData$east, y = BlackSmokePrefData$north, col = 'red')
mesh$n #2432 vertices

# # Load the regular mesh #
# mesh <- readRDS('./Reproducibility/mesh_5_7.rds')
# plot(mesh, asp = 1)
# points(x = BlackSmokePrefData$east, y = BlackSmokePrefData$north, col = 'red')
# mesh$n #2432 vertices
# # scale the mesh onto the transformed scale
# mesh$loc = mesh$loc / sd_x

# Compute Euclidean distances between all the sites #
Dists = spDists(cbind(BlackSmokePrefData$east, BlackSmokePrefData$north))

# Compute R_lag - the indicator for if the site was selected at time t-1
BlackSmokePrefData2$R_lag = c(rep(NA, no_sites), BlackSmokePrefData2$R[1:(dim(BlackSmokePrefData2)[1]-no_sites)])

# Observe the locations of the sites through time #
BlackSmokePrefData2$repulsion_ind = 0
counter = no_sites+1
for (i in sort(unique(BlackSmokePrefData2$year))[-1])
{
  # First extract the data at time i
  Data_i =  BlackSmokePrefData2[BlackSmokePrefData2$year == i,]
  
  print(plot(x = BlackSmokePrefData2$east[BlackSmokePrefData2$year==i & !is.na(BlackSmokePrefData2$bsmoke)],
             y = BlackSmokePrefData2$north[BlackSmokePrefData2$year==i & !is.na(BlackSmokePrefData2$bsmoke)],
             xlab = i, xlim = range(BlackSmokePrefData2$east), ylim = range(BlackSmokePrefData2$north)))
  
  # Compute the repulsion indicator. Was there a site at year i-1 within radius r of it
  BlackSmokePrefData2$repulsion_ind[counter:(counter+(no_sites-1))] = rowSums(Dists[,which(Data_i$R_lag==1)] < r) > 0
  counter = counter + no_sites
}

# create the Matern spde object for Y_grf #
spde_obj = inla.spde2.pcmatern(mesh=mesh, alpha = 2,
                               prior.range = c(0.04,0.05),
                               prior.sigma = c(1,0.01),
                               constr = T) # alpha = 2 implies 1st order smoothness ie 1 times differentiable (same as exponential)
# PC prior says we believe the lower 1st percentile of range is 3.4km (a fifth of the min range found by Shaddick and Zidek)
# 99th percentile for the GRF's standard deviation is 1. We don't believe sd higher.

# create the projector matrix #
A_proj = inla.spde.make.A(mesh = mesh,
                          loc = as.matrix(cbind( BlackSmokePrefData2$east, BlackSmokePrefData2$north )),
                          group = BlackSmokePrefData2$year-65, # group membership needs to be 1:no_T
                          n.group = no_T)
dim(A_proj) # nrows = number of sites * number of time points, ncols = no.nodes * no time points #
no_sites * no_T
mesh$n * no_T

s_index = inla.spde.make.index(name = "spatial.field",
                               n.spde = spde_obj$n.spde,
                               n.group = no_T)

time = (1:no_T)/no_T
time2 = time^2
# # create the stack object for estimating observation process y #
# cov_y = data.frame(year = rep(time, each = no_sites),
#                    year_2 = rep(time2, each = no_sites),
#                    spatial_ind = rep(1:no_sites, times = no_T),
#                    spatial_ind2 = no_sites + rep(1:no_sites, times = no_T)) # site-specific random intercepts


s_index_copy = s_index
names(s_index_copy) = c('spatial.field.copy',"spatial.field.group.copy", "spatial.field.repl.copy")

s_index_copy2 = s_index
names(s_index_copy2) = c('spatial.field.copy2',"spatial.field.group.copy2", "spatial.field.repl.copy2")


# stack_y_est = inla.stack(data = list(y = BlackSmokePrefData2$bsmoke, #single model
#                                      alldata = cbind(BlackSmokePrefData2$bsmoke, NA, NA),
#                                      Ntrials = rep(0,times = length(BlackSmokePrefData2$bsmoke))), #joint model
#                          A = list(A_proj,A_proj,A_proj,1),
#                          effects= list(c(s_index, list(Intercept = 1)),
#                                        c(s_index_copy, list(Intercept_copy = 1)),
#                                        c(s_index_copy2, list(Intercept_copy2 = 1)),
#                                        cov_y),
#                          tag = 'y_est')

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
# # Put a 0 in the unvisited nodes if and only if they lie inside the GB boundary. Else put NA.
# data_expand[,'R'] = rep(ifelse(xy_in == 1, 0, NA)[unvisited_nodes[1:number_unvisited_nodes_per_year]], times = no_T)
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
cov_y_expand = data.frame(year = (data_expand$year-65)/no_T,
                          year_2 = ((data_expand$year-65)/no_T)^2,
                          spatial_ind = c(rep(1:no_sites, times = no_T), rep(NA, times = dim(data_expand)[1] - length(BlackSmokePrefData2$R_lag))), # Only estimate random effects for site locations
                          spatial_ind2 = no_sites + c(rep(1:no_sites, times = no_T), rep(NA, times = dim(data_expand)[1] - length(BlackSmokePrefData2$R_lag))))

#diag_model = diag(mesh$n * nyears)
stack_y_est_expand = inla.stack(data = list(y = data_expand$bsmoke, #single model
                                            alldata = cbind(data_expand$bsmoke, NA,NA,NA),
                                            Ntrials = rep(0,times = length(data_expand$bsmoke))), #joint model
                                A = list(A_proj_expand, A_proj_expand, A_proj_expand, 1),
                                effects= list(c(s_index, list(Intercept = 1)),
                                              c(s_index_copy, list(Intercept_copy = 1)),
                                              c(s_index_copy2, list(Intercept_copy2 = 1)),
                                              cov_y_expand),
                                tag = 'y_est_expand')
print('stack y est expand complete')

# # create the stack object for estimating selection process R #
# cov_R = data.frame(R_lag = BlackSmokePrefData2$R_lag,
#                    yearR = rep(time, each = no_sites),
#                    yearR_2 = rep(time2, each = no_sites),
#                    repulsion_ind = BlackSmokePrefData2$repulsion_ind,
#                    R_year = BlackSmokePrefData2$year)
# #spatial_ind2 = sim_data$spatial_ind
# 
# R_s_index = s_index
# names(R_s_index) = c('R.spatial.field',"R.spatial.field.group", "R.spatial.field.repl")
# # change the R.spatial.field to 1:lengthlength(s_index_dummy$spatial.field.dummy) instead of rep(1:mesh$n, nyears)
# R_s_index$R.spatial.field = 1:(mesh$n * no_T)
# 
# R_s_index_copy = s_index
# names(R_s_index_copy) = c('R.spatial.field.copy',"R.spatial.field.group.copy", "R.spatial.field.repl.copy")
# R_s_index_copy2 = s_index
# names(R_s_index_copy2) = c('R.spatial.field.copy2',"R.spatial.field.group.copy2", "R.spatial.field.repl.copy2")
# #R_s_index$R.spatial.field.group = pmax(1, R_s_index$R.spatial.field.group-1) #lag the time by 1

# stack_R_est = inla.stack(data = list(R = BlackSmokePrefData2$R, #for single model
#                                      alldata = cbind(NA,BlackSmokePrefData2$R, NA),
#                                      Ntrials = rep(1,times = length(BlackSmokePrefData2$R))), #for joint model
#                          A = list(A_proj,A_proj,A_proj,1),
#                          effects= list(c(R_s_index, list(R.Intercept = 1)),
#                                        c(R_s_index_copy, list(R.Intercept_copy = 1)),
#                                        c(R_s_index_copy2, list(R.Intercept_copy2 = 1)),
#                                        cov_R),
#                          tag = 'R_est')


cov_R_expand = data.frame(R_lag = data_expand$R_lag,
                          repulsion_ind = data_expand$repulsion_ind,
                          R_year = (data_expand$year-65)/no_T,
                          yearR = (data_expand$year-65)/no_T,
                          yearR_2 = ((data_expand$year-65)/no_T)^2,
                          R_copyind = c(1:length(BlackSmokePrefData2$bsmoke), rep(NA, times = dim(data_expand)[1] - length(BlackSmokePrefData2$bsmoke))))

stack_R_est_expand = inla.stack(data = list(R = data_expand$R, #for single model
                                            alldata = cbind(NA, data_expand$R, NA, NA),
                                            Ntrials = rep(1, times = length(data_expand$R))), #for joint model
                                A = list(A_proj_expand, A_proj_expand, A_proj_expand, 1),
                                effects= list(c(R_s_index, list(R.Intercept = 1)),
                                              c(R_s_index_copy, list(R.Intercept_copy = 1)),
                                              c(R_s_index_copy2, list(R.Intercept_copy2 = 1)),
                                              cov_R_expand),
                                tag = 'R_est_expand')# Form the dummy variable with the zero observations #
print('stack R est expand complete')

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

# cov_dummy2_expand = data.frame(spatial_ind_dummy = c(rep(1:no_sites, times = no_T), no_sites + rep(1:number_unvisited_nodes_per_year, times = no_T)), # Only estimate random effects for site locations
#                                spatial_ind_dummy2 = dim(data_expand)[1]/no_T + c(rep(1:no_sites, times = no_T), no_sites+rep(1:number_unvisited_nodes_per_year, times = no_T)),
#                                time_dummy = c(((data_expand$year-65)/no_T)[1:(dim(data_expand)[1]/no_T)], ((data_expand$year-65)/no_T)[1:((dim(data_expand)[1]/no_T)*(no_T-1))]),
#                                spatial.field.dummy2 = 1:dim(data_expand)[1],
#                                spatial.field.repl.dummy2 = rep(1, times = dim(data_expand)[1]))


cov_dummy_est = data.frame(spatial_ind_dummy = rep(1:no_sites, times = no_T),
                           spatial_ind_dummy2 = no_sites + rep(1:no_sites, times = no_T),
                           time_dummy = c(rep(sort(unique((data_expand$year-65)/no_T))[1], times = no_sites), rep(sort(unique((data_expand$year-65)/no_T))[1:(no_T-1)], each = no_sites)),
                           spatial.field.dummy2 = 1:length(BlackSmokePrefData2$bsmoke),
                           spatial.field.repl.dummy2 = rep(1, times = length(BlackSmokePrefData2$bsmoke)))

# cov_dummy_expand = data.frame(spatial_ind_dummy = c(rep(1:no_sites, times = no_T),rep(NA, times = dim(data_expand)[1] - (no_sites*no_T))), # Only estimate random effects for site locations
#                               spatial_ind_dummy2 = no_sites + c(rep(1:no_sites, times = no_T),rep(NA, times = dim(data_expand)[1] - (no_sites*no_T))),
#                               time_dummy = c(rep(sort(unique((data_expand$year-65)/no_T))[1], times = no_sites), rep(sort(unique((data_expand$year-65)/no_T))[1:(no_T-1)], each = no_sites),rep(sort(unique((data_expand$year-65)/no_T))[c(1,1:(no_T-1))], each = (dim(data_expand)[1] - (no_sites*no_T))/no_T )))

stack_dummy_est = inla.stack(data = list(dummy2 = rep(0, times = length(BlackSmokePrefData2$bsmoke)), #single model
                                         alldata = cbind(rep(NA, times = length(BlackSmokePrefData2$bsmoke)),
                                                         rep(NA, times = length(BlackSmokePrefData2$bsmoke)),
                                                         rep(NA, times = length(BlackSmokePrefData2$bsmoke)),
                                                         rep(0, times = length(BlackSmokePrefData2$bsmoke))),
                                         Ntrials = rep(0,times = length(BlackSmokePrefData2$bsmoke))), #joint model
                             A = list(1),
                             effects= list(cov_dummy_est),
                             tag = 'dummy_est')

# # This is for copying the random intercepts and slopes across
# stack_dummy2_expand = inla.stack(data = list(dummy2 = rep(0, times = dim(data_expand)[1]), #single model
#                                              alldata = cbind(rep(NA, times = dim(data_expand)[1]),
#                                                              rep(NA, times = dim(data_expand)[1]),
#                                                              rep(NA, times = dim(data_expand)[1]),
#                                                              rep(0, times = dim(data_expand)[1])),
#                                              Ntrials = rep(0,times = dim(data_expand)[1])), #joint model
#                                  A = list(1),
#                                  effects= list(cov_dummy2_expand),
#                                  tag = 'dummy2_expand')
# print('stack dummy2 expand complete')

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
print('stack dummy est expand complete')

#stack_y_combine = inla.stack(stack_y_est, stack_y_pred)
#stack_combine_est = inla.stack(stack_y_est, stack_y_pred, stack_R_est, stack_R_pred, stack_dummy_est, stack_dummy_pred)
#stack_joint1 = inla.stack(stack_y_est, stack_R_est, stack_dummy_est)
#print('stack joint1 complete')

stack_joint2 = inla.stack(stack_y_est_expand, stack_R_est_expand, stack_dummy_est_expand, stack_dummy_est)
print('stack joint2 complete')

# Step 1 - fit the naive model and save object to file #
# formula_naive = y ~ -1 + Intercept +
#   f(spatial.field, model = spde_obj) +
#   f(spatial.field.copy, I(spatial.field.group.copy/no_T), model = spde_obj) +
#   f(spatial.field.copy2, I((spatial.field.group.copy2/no_T)^2), model = spde_obj) +
#   I(spatial.field.group/no_T) + I((spatial.field.group/no_T)^2) +
#   f(year, model = 'ar1') +
#   f(spatial_ind, model = "iid2d", n = no_sites*2, constr = TRUE) + #random site-specific intercepts
#   f(spatial_ind2, year, copy = "spatial_ind") # random site-specific slopes
# #f(spatial_ind, model = 'iid', hyper=list(theta=list(prior="loggamma",param=c(2,0.04)))) +
#
# theta.ini = c(1.597900, -1.277423, -0.443820, -1.441220, 0.036510, -1.441336, 0.016919,
#               4.462918, 1.437147, 4, 4, 4)
# out.naive = inla(formula_naive, family = 'gaussian',
#                  data = inla.stack.data(stack_y_est),
#                  control.predictor = list(A = inla.stack.A(stack_y_est), compute = F),
#                  control.compute = list(dic=F, config = T, cpo = F),
#                  control.fixed = list(mean = list(Intercept = 1, default = 0),
#                                       prec = list(Intercept = 0.25, default = 0.001)),
#                  control.results = list(return.marginals.random = F,
#                                         return.marginals.predictor = F),
#                  control.inla = list(h = 0.00001, strategy = "gaussian", int.strategy = 'eb'),
#                  control.inla = list(strategy = "gaussian", int.strategy = 'eb'),
#                  control.mode = list(theta = theta.ini, restart=T),
#                  verbose = T, num.threads=20)
# summary(out.naive)
#
# saveRDS(out.naive,file="finalmod_naive2_extended_approx.rds")
# rm(out.naive)
#
# # Fit joint model version 1 - without zeroes over whole grid #
#  formula_joint1 = alldata ~ -1 + Intercept +
#    f(spatial.field, model = spde_obj)+ #spatial.field.group) +
#    f(R.spatial.field, copy = "spatial.field.dummy", fixed = F) +
#    f(R.spatial.field.copy, model = spde_obj) +
#    f(spatial.field.copy, I(spatial.field.group.copy/no_T), model = spde_obj) +
#    f(spatial.field.copy2, I((spatial.field.group.copy2/no_T)^2), model = spde_obj) +
#    R_lag + year + year_2 + yearR + yearR_2 + repulsion_ind +
#    f(R_year, model="ar1",hyper=list(theta1=list(prior="pcprec",param=c(2,0.01))))+
#    I(R_year==1) +
#    f(spatial.field.dummy, I(-spatial.field.repl.dummy), model="iid", hyper = list(prec = list(initial = -20, fixed=TRUE))) +
#    f(spatial.field.copy.dummy, copy = "spatial.field", fixed = T) + #, replicate = I(pmax(1,spatial.field.group.copy.dummy-1))
#    f(spatial.field.copy.dummy2, I(spatial.field.group.copy.dummy2/no_T), copy = "spatial.field.copy", fixed = T) +
#    f(spatial.field.copy.dummy3, I((spatial.field.group.copy.dummy3/no_T)^2), copy = "spatial.field.copy2", fixed = T) +
#    f(spatial_ind, model = "iid2d", n = no_sites*2, constr = TRUE) + #random site-specific intercepts
#    f(spatial_ind2, year, copy = "spatial_ind") + # random site-specific slopes
#    f(spatial_ind_dummy, copy = "spatial_ind", fixed = T) +
#    f(spatial_ind_dummy2, time_dummy, copy = 'spatial_ind', fixed = T)
#
#I(R.spatial.field.group/30) + I((R.spatial.field.group/30)^2) + f(spatial_ind2, copy = 'spatial_ind', fixed=F)
#f(spatial_ind, model = 'iid',hyper=list(theta=list(prior="loggamma",param=c(2,0.04)))) +
#theta.ini2  = c(1.5600280,-1.3272852,-0.5072330,-2.4681359,-2.7030130,-1.4166110,
#                 -0.0262280, -1.3215250, -0.0538500,  0.8493710,  0.3394635, 4, 4, 4, 0.3733780, 0.37)
theta.ini2 = c(2.95712439, -0.39504940, -0.02055101, -2.30589883, -2.52425155, -2.22719145, -0.52338512,
               -2.22312867, -0.38897875,  0.93644135,  0.63899640,  0.51169121, -0.10588640, -2.10631419,
               2.77008655,  0.11748417)

# out.joint1 = inla(formula_joint1, family = c('gaussian','binomial','gaussian'),
#                  data = inla.stack.data(stack_joint1),
#                  Ntrials=inla.stack.data(stack_joint1)$Ntrials,
#                  control.predictor = list(A = inla.stack.A(stack_joint1), compute = F),
#                  control.compute = list(dic=F, config = T, cpo = F),
#                  control.results = list(return.marginals.random = F,
#                                         return.marginals.predictor = F),
#                  control.fixed = list(mean = list(Intercept = 1, default = 0),
#                                       prec = list(Intercept = 0.25, default = 0.001)),
#                  control.inla = list(h = 0.00001, strategy = "gaussian", int.strategy = 'eb'),
#                  control.mode = list(theta = theta.ini2, restart=T),
#                  control.family = list(
#                    list(),
#                    list(),
#                    list(hyper = list(prec = list(initial = 20, fixed=TRUE)))),
#                  verbose = T, num.threads=20)
# summary(out.joint1)
#
# saveRDS(out.joint1,file="finalmod_joint1_extended_approx.rds")
# theta.ini3 = out.joint1$mode$theta
# rm(out.joint1)

# Fit joint model version 2 #
formula_joint2 = alldata ~ -1 + Intercept +
  f(spatial.field, model = spde_obj)+ #spatial.field.group) +
  f(R.spatial.field, copy = "spatial.field.dummy", fixed = F) +
  f(R.spatial.field.copy, model = spde_obj) +
  f(R_copyind, copy = "spatial.field.dummy2", fixed = F) +
  f(spatial.field.copy, I(spatial.field.group.copy/no_T), model = spde_obj) +
  f(spatial.field.copy2, I((spatial.field.group.copy2/no_T)^2), model = spde_obj) +
  R_lag + year + year_2 + yearR + yearR_2 + repulsion_ind +
  f(R_year, model="ar1",hyper=list(theta1=list(prior="pcprec",param=c(2,0.01))))+
  I(R_year==1/no_T) +
  f(spatial.field.dummy, I(-spatial.field.repl.dummy), model="iid", hyper = list(prec = list(initial = -20, fixed=TRUE))) +
  f(spatial.field.dummy2, I(-spatial.field.repl.dummy2), model="iid", hyper = list(prec = list(initial = -20, fixed=TRUE))) +
  f(spatial.field.copy.dummy, copy = "spatial.field", fixed = T) + #, replicate = I(pmax(1,spatial.field.group.copy.dummy-1))
  f(spatial.field.copy.dummy2, I(spatial.field.group.copy.dummy2/no_T), copy = "spatial.field.copy", fixed = T) +
  f(spatial.field.copy.dummy3, I((spatial.field.group.copy.dummy3/no_T)^2), copy = "spatial.field.copy2", fixed = T) +
  f(spatial_ind, model = "iid2d", n = no_sites*2, constr = TRUE) + #random site-specific intercepts
  f(spatial_ind2, year, copy = "spatial_ind") + # random site-specific slopes
  f(spatial_ind_dummy, copy = "spatial_ind", fixed = T) +
  f(spatial_ind_dummy2, time_dummy, copy = 'spatial_ind', fixed = T)

out.joint2 = inla(formula_joint2, family = c('gaussian','binomial','gaussian','gaussian'),
                  data = inla.stack.data(stack_joint2),
                  Ntrials=inla.stack.data(stack_joint2)$Ntrials,
                  control.predictor = list(A = inla.stack.A(stack_joint2), compute = F),
                  control.compute = list(dic=F, config = T, cpo = F),
                  control.fixed = list(mean = list(Intercept = 1, default = 0),
                                       prec = list(Intercept = 0.25, default = 0.001)),
                  # control.results = list(return.marginals.random = F,
                  #                        return.marginals.predictor = F),
                  control.inla = list(h = 0.00001, int.strategy = 'eb'),
                  control.mode = list(theta = theta.ini2, restart=T),
                  control.family = list(
                    list(),
                    list(),
                    list(hyper = list(prec = list(initial = 20, fixed=TRUE))),
                    list(hyper = list(prec = list(initial = 20, fixed=TRUE)))),
                  verbose = T, num.threads=8)
summary(out.joint2)

# saveRDS(out.joint2,
#         file="finalmod_joint2_reproduce.rds")


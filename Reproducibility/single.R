library(dplyr)
library(sp)
library(rgdal)
library(reshape2)
library(INLA)
library(ggplot2)
library(inlabru)


# load data #
load("./Reproducibility/Data2Joe.RData")
xy_in = readRDS("./Data/Reproducibility/xy_in.rds") # load the indicator telling us if the mesh vertices lie in GB

# The columns are: site name, x, y (both in British National Grid),
# black smoke by year (mu gm/3) where NA means it wasn’t measured.
sd_x = sd(BlackSmokePrefData[, c(2)])
sd_y = sd(BlackSmokePrefData[, c(3)])

BlackSmokePrefData[,c(2,3)] = BlackSmokePrefData[,c(2,3)]/sd_x # standardize x,y coords (dividing by 86021)
plot(x = BlackSmokePrefData$east, y = BlackSmokePrefData$north)

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
ind_select = as.matrix(!is.na(BlackSmokePrefData[,-c(1,2,3)]) )
R = as.numeric(!is.na(BlackSmokePrefData[,-c(1,2,3)]) )
BlackSmokePrefData2$R = R # selection variable (1 if site selected at time t)

ncol = 100 # grid for projection
nrow = 100
L = nrow * ncol # number of grid sites
east_grid = seq(from = min(BlackSmokePrefData$east, na.rm=T),
                to = max(BlackSmokePrefData$east, na.rm=T), length.out = ncol)
north_grid = seq(from = min(BlackSmokePrefData$north, na.rm=T),
                 to = max(BlackSmokePrefData$north, na.rm=T), length.out = nrow)


# Load the regular mesh #
mesh <- readRDS('./Data/Reproducibility/mesh_5_7.rds')

plot(mesh, asp = 1)
points(x = BlackSmokePrefData$east*sd_x, y = BlackSmokePrefData$north*sd_x, col = 'red')

mesh$n #2432 vertices

# scale the mesh onto the transformed scale
mesh$loc = mesh$loc / sd_x


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
  
  # print(plot(x = BlackSmokePrefData2$east[BlackSmokePrefData2$year==i & !is.na(BlackSmokePrefData2$bsmoke)],
  #            y = BlackSmokePrefData2$north[BlackSmokePrefData2$year==i & !is.na(BlackSmokePrefData2$bsmoke)],
  #            xlab = i, xlim = range(BlackSmokePrefData2$east), ylim = range(BlackSmokePrefData2$north)))
  
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


# Step 1 - fit the naive model and save object to file #
formula_naive = y ~ -1 + Intercept +
  f(spatial.field, model = spde_obj) +
  f(spatial.field.copy, I(spatial.field.group.copy/no_T), model = spde_obj) +
  f(spatial.field.copy2, I((spatial.field.group.copy2/no_T)^2), model = spde_obj) +
  I(spatial.field.group/no_T) + I((spatial.field.group/no_T)^2) +
  f(year, model = 'ar1') +
  f(spatial_ind, model = "iid2d", n = no_sites*2, constr = TRUE) + #random site-specific intercepts
  f(spatial_ind2, year, copy = "spatial_ind") # random site-specific slopes
#f(spatial_ind, model = 'iid', hyper=list(theta=list(prior="loggamma",param=c(2,0.04)))) +

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
                 # control.inla = list(h = 0.00001, strategy = "gaussian", int.strategy = 'eb'),
                 control.inla = list(strategy = "gaussian", int.strategy = 'eb'),
                 # control.mode = list(theta = theta.ini, restart=T),
                 verbose = T, num.threads=20)
summary(out.naive)
#
# saveRDS(out.naive,file="finalmod_naive2_extended_approx.rds")
# rm(out.naive)


# The INLA naive model does not converge without setting initial values.



# Analysis of the final non-PS joint model #

# Read in UK shapefile #
# library(raster)
library(terra)
require("rgdal")
require("rgeos")
require("dplyr")
library(INLA)
#INLA:::inla.dynload.workaround() #make it compatible with older linux
#library(RandomFields)
library(mvtnorm)
library(boot)
library(geoR)
library(reshape2)
library(sp)
library(ggplot2)


load("./Data2Joe.RData")

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
ind_select = as.matrix(!is.na(BlackSmokePrefData[,-c(1,2,3)]) )
R = as.numeric(!is.na(BlackSmokePrefData[,-c(1,2,3)]) )
BlackSmokePrefData2$R = R # selection variable (1 if site selected at time t)

# Generate the repulsion_ind variable for analysis #
Dists = spDists(cbind(BlackSmokePrefData$east, BlackSmokePrefData$north))
repulsion_ind = matrix(0, nrow = no_sites, ncol = no_T)
counter3 = 2
for(i in sort(unique(BlackSmokePrefData2$year))[-1])
{
  repulsion_ind[,counter3] = rowSums(Dists[,which(BlackSmokePrefData2[BlackSmokePrefData2$year == (i-1),]$R == 1)] < r) > 0
  counter3 = counter3 + 1
}

ncol = 100 # grid for projection
nrow = 100
L = nrow * ncol # number of grid sites
east_grid = seq(from = min(BlackSmokePrefData$east, na.rm=T),
                to = max(BlackSmokePrefData$east, na.rm=T), length.out = ncol)
north_grid = seq(from = min(BlackSmokePrefData$north, na.rm=T),
                 to = max(BlackSmokePrefData$north, na.rm=T), length.out = nrow)


# Form the regular mesh #
# Create a rough convex boundary for the UK #
# Form the grid independently from the sites to avoid preferential grid selection #
UK_domain = cbind(c(2, 7.7, 7.7, 6, 4, 1), c(0.5, 0.5, 6, 13.5, 13.5, 12))
hull = inla.nonconvex.hull(cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north))

# Set cutoff distance (minimum edge) and max.edge to less than the process range #
# 17km was found to be the smallest range for the spatial fields fit in Shaddick and Zidek (2014)
# set max.edge and cutoff to be 10km #

######################
# Regular grid satisfying the min and max edge lengths < practical range is too large
# Need to form irregular mesh and use inla.project to project onto a regular lattice afterwards
# Will affect the de-biasing effect over the whole grid
# mesh = inla.mesh.2d(boundary = hull,
#                     offset = c(0.1,0.2), max.edge = cutoff_dist,
#                     cutoff = c(cutoff_dist, cutoff_dist*3),
#                     min.angle = 30)
# plot(mesh)
#######################

cutoff_dist = 10000/sd_x # 10km < 17km as min edge
max.edge = 100000/sd_x # 100km max edge as in Shaddick and Zidek

mesh = inla.mesh.2d(loc = cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north),
                    boundary = hull,
                    offset = c(0.1,0.2), max.edge = max.edge,
                    cutoff = c(cutoff_dist, 0.5),
                    min.angle = 26)
#max.edge = c(0.05,0.09),
#cutoff = 0.5)
plot(mesh, asp = 1)
points(x = BlackSmokePrefData$east, y = BlackSmokePrefData$north, col = 'red')

mesh$n #1243 vertices


# download.file("https://borders.ukdataservice.ac.uk/ukborders/easy_download/prebuilt/shape/infuse_uk_2011.zip",
#               destfile = "lad-region-lookup.zip")
# unzip("lad-region-lookup.zip", exdir = ".")
region <- readOGR(".", "infuse_uk_2011")
UK_polygon = region@polygons

plot(region)

#####################################
### A matrix for projecting SPDEs ###
#####################################
# Define lattice square side length in kilometers - e.g 5km
stepsize = 5000/sd_x
nxy = round(c(diff(range(BlackSmokePrefData2$east)),
              diff(range(BlackSmokePrefData2$north)))/stepsize)
proj_grid = inla.mesh.projector(mesh, xlim = range(BlackSmokePrefData2$east),
                                ylim = range(BlackSmokePrefData2$north),
                                dims = nxy)

# Find all locations in the inla.projector lattice that fall outside of mainland Great Britain
# download.file("https://borders.ukdataservice.ac.uk/ukborders/easy_download/prebuilt/shape/infuse_gb_2011.zip",
#               destfile = "lad-region-lookup2.zip")
# unzip("lad-region-lookup2.zip", exdir = ".")
region_GB <- readOGR(".", "infuse_gb_2011")
GB_polygons = polygons(region_GB)
GB_polygons@polygons[["Polygon"]]

points <- SpatialPoints(coords = proj_grid$lattice$loc*sd_x,
                       proj4string=CRS(proj4string(region_GB)))
table(xy.in <- gContains(region_GB,points, byid=T))

############################################################################
### Sampling from Posterior distribution or load in pre-compiled Samples ###
############################################################################
# Creating Samples
m_samples = 1000
#mod2 = readRDS('LOCATION_OF_MODEL.rds')
#samp2 <- inla.posterior.sample(n=m_samples, mod2)

# Load pre-compiled samples
samp2 = readRDS('LOCATION_OF_SAMPLES')

# reduce the size of samp by removing the reduntant information #
for (i in 1:1000)
{
  samp2[[i]]$latent = samp2[[i]]$latent[-c(grep('Predictor', rownames(samp2[[i]]$latent), fixed = T)),]
}

# Creating the A projector Matrix
A.grid2 = inla.spde.make.A(mesh,
                         loc = cbind(BlackSmokePrefData$east, BlackSmokePrefData$north))

# 45466 Apredictors and 50252 Predictors
tmp2 = matrix(0, nrow = no_T, ncol = mesh$n)
year = (1:no_T)/no_T
year2 = year^2
Post_Array2 = array(0, dim = c(m_samples, no_T, mesh$n))
Post_Sites_Array2 = array(0, dim = c(m_samples, no_T, dim(BlackSmokePrefData)[1]))
Post_R_Array2 = array(0, dim = c(m_samples, no_T, mesh$n))
Post_R_Sites_Array2 = array(0, dim = c(m_samples, no_T, dim(BlackSmokePrefData)[1]))
Post_R_Sites_Array3 = array(0, dim = c(m_samples, no_T, dim(BlackSmokePrefData)[1])) # for the HT estimator
residuals_array = array(0, dim = c(m_samples, no_T, dim(BlackSmokePrefData)[1]))
RI_array = array(0, dim = c(m_samples, no_T, dim(BlackSmokePrefData)[1])) # random intercepts
RS_array = array(0, dim = c(m_samples, no_T, dim(BlackSmokePrefData)[1])) # random slopes

beta_gamma = rnorm(m_samples, mean = 0.66, sd = 0.1748)#rep( samp2[[1]]$hyperpar['Beta for R.spatial.field'], m_samples) #
beta_b = rnorm(m_samples, mean = 0.061, sd = 0.0429)#rep( samp2[[1]]$hyperpar['Beta for R_copyind'], m_samples)#
# Looping for each sample
for (i in 1:m_samples){
  # Creating predictions for Y process
  print(i)
  print('from')
  print(m_samples)
  counter = 0
  counter2 = 0
  for(j in 1:no_T) # loop over the years
  {
    Post_Array2[i,j,] <- samp2[[i]]$latent['Intercept'] + # Intercept
    samp2[[i]]$latent['year'] * year[j] + # linear fixed effect
    samp2[[i]]$latent['year_2'] * year2[j] + # quadratic fixed effect
    as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial.field:',x = names(samp2[[i]]$latent)))]) + # Spatial Intercept
    as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial.field.copy:',x = names(samp2[[i]]$latent)))]) * year[j] + # Spatial linear slope
    as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial.field.copy2:',x = names(samp2[[i]]$latent)))]) * year2[j] #+ # Spatial quadratic slope
    #samp[[i]]$latent[grep('year:',rownames(samp[[i]]$latent), fixed = TRUE)[j],1] # AR1 zero mean time effect

    # Repeat to obtain the selection probabilities
    if(j == 1)
    {
      Post_R_Array2[i,j,] <- samp2[[i]]$latent['I(R_year == 1/31)TRUE'] + # Intercept
      samp2[[i]]$latent['yearR'] * year[j] + # linear fixed effect
      samp2[[i]]$latent['yearR_2'] * year2[j] + # quadratic fixed effect
      beta_gamma[j] * as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'R.spatial.field:',x = names(samp2[[i]]$latent)))[counter + c(1:mesh$n)]]) + # copied linear combination of latent processes
      as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'R.spatial.field.copy:',x = names(samp2[[i]]$latent)))]) + # R spatial residual field
      as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'R_year:',x = names(samp2[[i]]$latent)))[j]])  # AR(1) process
    }
    if(j!=1)
    {
      Post_R_Array2[i,j,] <- samp2[[i]]$latent['I(R_year == 1/31)FALSE'] + # Intercept
        samp2[[i]]$latent['yearR'] * year[j] + # linear fixed effect
        samp2[[i]]$latent['yearR_2'] * year2[j] + # quadratic fixed effect
        beta_gamma[j] * as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'R.spatial.field:',x = names(samp2[[i]]$latent)))[counter + c(1:mesh$n)]]) + # copied linear combination of latent processes
        as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'R.spatial.field.copy:',x = names(samp2[[i]]$latent)))]) + # R spatial residual field
        as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'R_year:',x = names(samp2[[i]]$latent)))[j]])  # AR(1) process
    }


    # Project onto the site locations and add random site-specific effects #
    Post_Sites_Array2[i,j,] = as.numeric(A.grid2 %*% Post_Array2[i,j,]) +
                              as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial_ind:',x = names(samp2[[i]]$latent)))[1:no_sites]]) + # random Intercept
                              as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial_ind2:',x = names(samp2[[i]]$latent)))[no_sites+(1:no_sites)]]) * year[j] # random slope

    # Project selection field onto the site locations and add random site-specific effects and R covariates #
    if(j == 1)
    {
      Post_R_Sites_Array2[i,j,] = as.numeric(A.grid2 %*% Post_R_Array2[i,j,]) +
      as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'R_copyind:',x = names(samp2[[i]]$latent)))[counter2 + c(1:no_sites)]]) # copied linear combinations of site-specific random effects

      Post_R_Sites_Array3[i,j,] = as.numeric(A.grid2 %*% Post_R_Array2[i,j,]) +
      as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'R_copyind:',x = names(samp2[[i]]$latent)))[counter2 + c(1:no_sites)]]) # copied linear combinations of site-specific random effects
    }

    if(j != 1)
    {
      Post_R_Sites_Array2[i,j,] = as.numeric(A.grid2 %*% Post_R_Array2[i,j,]) +
      beta_b[j] * as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'R_copyind:',x = names(samp2[[i]]$latent)))[counter2 + c(1:no_sites)]]) + # copied linear combinations of site-specific random effects
      samp2[[i]]$latent['R_lag'] * R[(counter2-no_sites) + c(1:no_sites)] +
      samp2[[i]]$latent['repulsion_ind']*repulsion_ind[,j]

      Post_R_Sites_Array3[i,j,] = as.numeric(A.grid2 %*% Post_R_Array2[i,j,]) +
      as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'R_copyind:',x = names(samp2[[i]]$latent)))[counter2 + c(1:no_sites)]]) # copied linear combinations of site-specific random effects
    }

    # Extract the random intercepts #
    RI_array[i,j,] = as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial_ind:',x = names(samp2[[i]]$latent)))[1:no_sites]])

    # Extract the random slopes #
    RS_array[i,j,] = as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial_ind2:',x = names(samp2[[i]]$latent)))[1:no_sites]])

    # Extract the residuals #
    residuals_array[i,j,] = log(BlackSmokePrefData[,(3+j)] / mean(colMeans(BlackSmokePrefData[,4:34], na.rm = T))) - Post_Sites_Array2[i,j,]

    counter = counter + mesh$n
    counter2 = counter2 + no_sites
  }
}

Post_R_Sites_Array2 = inv.logit(Post_R_Sites_Array2)
Post_R_Sites_Array3 = inv.logit(Post_R_Sites_Array3)

# Summarising predictions
Post_mean2 <- t(apply(Post_Array2, c(2,3), function(x) {mean(x, na.rm=TRUE)}))
Post_sd2 <- t(apply(Post_Array2, c(2,3), function(x) {sd(x, na.rm=TRUE)}))
Post_R_mean2 <- t(apply(Post_R_Array2, c(2,3), function(x) {mean(x, na.rm=TRUE)}))
Post_R_sd2 <- t(apply(Post_R_Array2, c(2,3), function(x) {sd(x, na.rm=TRUE)}))
Post_LCL2 <- apply((apply(Post_Array2, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.025))
Post_UCL2 <- apply((apply(Post_Array2, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.975))
RI_array <- t(apply(RI_array, c(2,3), function(x) {mean(x, na.rm=TRUE)}))
RS_array <- t(apply(RS_array, c(2,3), function(x) {mean(x, na.rm=TRUE)}))
residuals_array <- t(apply(residuals_array, c(2,3), function(x) {mean(x, na.rm=TRUE)}))

# Project onto the grid using inla.projector
post_Matrix_grid_mean2 = array(NA, dim = c(nxy[1], nxy[2], no_T))
post_Matrix_grid_sd2 = array(NA, dim = c(nxy[1], nxy[2], no_T))
post_Matrix_grid_R_mean2 = array(NA, dim = c(nxy[1], nxy[2], no_T))
post_Matrix_grid_R_sd2 = array(NA, dim = c(nxy[1], nxy[2], no_T))

for(i in 1:no_T)
{
  post_Matrix_grid_mean2[,,i] = inla.mesh.project(proj_grid, Post_mean2[,i])
  post_Matrix_grid_mean2[,,i][!xy.in] <- NA
  post_Matrix_grid_sd2[,,i] = inla.mesh.project(proj_grid, Post_sd2[,i])
  post_Matrix_grid_sd2[,,i][!xy.in] <- NA

  post_Matrix_grid_R_mean2[,,i] = inla.mesh.project(proj_grid, Post_R_mean2[,i])
  post_Matrix_grid_R_mean2[,,i][!xy.in] <- NA
  post_Matrix_grid_R_sd2[,,i] = inla.mesh.project(proj_grid, Post_R_sd2[,i])
  post_Matrix_grid_R_sd2[,,i][!xy.in] <- NA
}

par(mfrow = c(2,2))
image(post_Matrix_grid_mean2[,,1], asp=1, main = 'Plot of the expected pollution in 1966 (red is low)')
image(post_Matrix_grid_mean2[,,30], asp=1, main = 'Plot of the expected pollution in 1996')
image(post_Matrix_grid_sd2[,,1], asp=1, main = 'Plot of the se of the expected pollution in 1966 (red is low)')
image(post_Matrix_grid_sd2[,,30], asp=1, main = 'Plot of the se of the expected pollution in 1996')

par(mfrow = c(2,2))
image(post_Matrix_grid_R_mean2[,,1], asp=1, main = 'Plot of the selection probability field in 1966 (red is low)')
image(post_Matrix_grid_R_mean2[,,30], asp=1, main = 'Plot of the selection probability field in 1996')
image(post_Matrix_grid_R_sd2[,,1], asp=1, main = 'Plot of the se of the selection probability field in 1966 (red is low)')
image(post_Matrix_grid_R_sd2[,,30], asp=1, main = 'Plot of the se of the selection probability field in 1996')

# # Step one - Make results tables for fixed parameter estimates
fixed_results = as.data.frame(mod1$summary.fixed)

# # Step two - Make results tables for hyperparameter estimates
hyper_results = as.data.frame(mod1$summary.hyperpar)

# # Step three - Make results tables for selected (R=1) sites estimates at each t
site_results2 = matrix(0, nrow = no_T, ncol = 4) # prediction means, sd, LCL, UCL
colnames(site_results2) = c('predciction mean', 'prediction sd', 'LCL', 'UCL')

R1_results2 = matrix(NA, nrow = no_T, ncol = 4) # prediction means, sd, LCL, UCL
colnames(R1_results2) = c('predciction mean', 'prediction sd', 'LCL', 'UCL')

# # Step four - Make results tables for unselected (R=0) sites estimates at each t
R0_results2 = matrix(NA, nrow = no_T, ncol = 4) # prediction means, sd, LCL, UCL
colnames(R0_results2) = c('predciction mean', 'prediction sd','LCL', 'UCL')

# Step 5 - Make results table for consistent sites vs those that are dropped
library(plyr)
temp_df = ddply(BlackSmokePrefData2, .(site), summarize, nR = sum(R))
consistent_ind = which(temp_df$nR == 31)
inconsistent_ind = which(temp_df$nR != 31)

# Step 6 - Make results tables for the Horvitz-Thompson Estimated mean #
HT_results = matrix(NA, nrow = no_T, ncol = 4)
colnames(HT_results) = c('predciction mean', 'prediction sd', 'LCL', 'UCL')

consistent_results2 = matrix(NA, nrow = no_T, ncol = 4) # prediction means, sd, LCL, UCL
colnames(consistent_results2) = c('predciction mean', 'prediction sd', 'LCL', 'UCL')

inconsistent_results2 = matrix(NA, nrow = no_T, ncol = 4) # prediction means, sd, LCL, UCL
colnames(inconsistent_results2) = c('predciction mean', 'prediction sd', 'LCL', 'UCL')

simple_ind_site = 1:dim(BlackSmokePrefData2)[1]
R_1_ind = which(BlackSmokePrefData2$R==1)
R_0_ind = which(BlackSmokePrefData2$R==0)

site_results2[,3] = apply((apply(Post_Sites_Array2, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.025))
site_results2[,4] = apply((apply(Post_Sites_Array2, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.975))

  for(j in 1:no_T)
  {
    site_results2[j,1] = mean(Post_Sites_Array2[,j,], na.rm=T)
    site_results2[j,2] = sd(apply(Post_Sites_Array2[,j,], 1, mean, na.rm=T), na.rm=T)

    R1_results2[j,1] = mean(Post_Sites_Array2[,j,which(!is.na(BlackSmokePrefData[,3+j]))], na.rm=T)
    R1_results2[j,2] = sd(apply(Post_Sites_Array2[,j,which(!is.na(BlackSmokePrefData[,3+j]))], c(1), mean, na.rm=T))
    R1_results2[j,3] = quantile(apply(Post_Sites_Array2[,j,which(!is.na(BlackSmokePrefData[,3+j]))], 1, mean, na.rm=T), probs = 0.025)
    R1_results2[j,4] = quantile(apply(Post_Sites_Array2[,j,which(!is.na(BlackSmokePrefData[,3+j]))], 1, mean, na.rm=T), probs = 0.975)

    R0_results2[j,1] = mean(Post_Sites_Array2[,j,which(is.na(BlackSmokePrefData[,3+j]))], na.rm=T)
    R0_results2[j,2] = sd(apply(Post_Sites_Array2[,j,which(is.na(BlackSmokePrefData[,3+j]))], c(1), mean, na.rm=T))
    R0_results2[j,3] = quantile(apply(Post_Sites_Array2[,j,which(is.na(BlackSmokePrefData[,3+j]))], 1, mean, na.rm=T), probs = 0.025)
    R0_results2[j,4] = quantile(apply(Post_Sites_Array2[,j,which(is.na(BlackSmokePrefData[,3+j]))], 1, mean, na.rm=T), probs = 0.975)

    # consistent_results2[j,1] = mean(Post_Sites_Array2[,j,consistent_ind], na.rm=T)
    # consistent_results2[j,2] = sd(apply(Post_Sites_Array2[,j,consistent_ind], c(1), mean, na.rm=T))
    # consistent_results2[j,3] = quantile(apply(Post_Sites_Array2[,j,consistent_ind], 1, mean, na.rm=T), probs = 0.025)
    # consistent_results2[j,4] = quantile(apply(Post_Sites_Array2[,j,consistent_ind], 1, mean, na.rm=T), probs = 0.975)
    #
    # inconsistent_results2[j,1] = mean(Post_Sites_Array2[,j,inconsistent_ind], na.rm=T)
    # inconsistent_results2[j,2] = sd(apply(Post_Sites_Array2[,j,inconsistent_ind], c(1), mean, na.rm=T))
    # inconsistent_results2[j,3] = quantile(apply(Post_Sites_Array2[,j,inconsistent_ind], 1, mean, na.rm=T), probs = 0.025)
    # inconsistent_results2[j,4] = quantile(apply(Post_Sites_Array2[,j,inconsistent_ind], 1, mean, na.rm=T), probs = 0.975)
    #
    # HT_results[j,1] = mean(apply(Post_R_Sites_Array3[,j,], 1, FUN = function(x) sum(BlackSmokePrefData[,j+3]/x, na.rm=T)/no_sites))
    # HT_results[j,2] = sd(apply(Post_R_Sites_Array3[,j,], 1, FUN = function(x) sum(BlackSmokePrefData[,j+3]/x, na.rm=T)/no_sites))
    # HT_results[j,3] = quantile(apply(Post_R_Sites_Array3[,j,], 1, FUN = function(x) sum(BlackSmokePrefData[,j+3]/x, na.rm=T)/no_sites), probs = 0.025)
    # HT_results[j,4] = quantile(apply(Post_R_Sites_Array3[,j,], 1, FUN = function(x) sum(BlackSmokePrefData[,j+3]/x, na.rm=T)/no_sites), probs = 0.975)
  }

# # Step five - Make results table for whole grid prediction at each t
Grid_results2 = matrix(NA, nrow = no_T, ncol = 4) # prediction means, sd, LCL, UCL
colnames(Grid_results2) = c('predciction mean', 'prediction sd', 'LCL', 'UCL')

for(i in 1:no_T)
{
  Grid_results2[i,1] = mean(post_Matrix_grid_mean2[,,i], na.rm=T)
  Grid_results2[i,2] = mean(post_Matrix_grid_sd2[,,i], na.rm=T)
}
Grid_results2[,3] = Post_LCL2
Grid_results2[,4] = Post_UCL2

Grid_df_2 = data.frame(y = c(Grid_results2[,1],R1_results2[,1],R0_results2[,1]),
                     x = rep(66:96, times = 3),
                     ymin = c(Post_LCL2, R1_results2[,3], R0_results2[,3]),
                     ymax = c(Post_UCL2, R1_results2[,4], R0_results2[,4]),
                     group = c(rep('Whole UK', times = 31), rep('R1', times = 31), rep('R0', times = 31)))
Grid_plot_2 = ggplot(aes(x = x, y = y, ymin = ymin, ymax = ymax, alpha = 0, linetype=group), data = Grid_df_2) +
  geom_ribbon(aes(colour = group, alpha = 0, fill = group)) +
  xlab('year') + ylab('posterior mean bs on transformed scale') +
  ggtitle('Mod2 posterior means of BS at the selected, unselected sites and across the UK (R1, R0 and Whole UK resp.). ')
Grid_plot_2

# On the original scale
Grid_df2_2 = data.frame(y = c(exp(c(Grid_results2[,1],R1_results2[,1],R0_results2[,1]))*30.7),
                     x = rep(66:96, times = 3),
                     ymin = c(exp(c(Post_LCL2, R1_results2[,3], R0_results2[,3]))*30.7),
                     ymax = c(exp(c(Post_UCL2, R1_results2[,4], R0_results2[,4]))*30.7),
                     group = c(rep('Whole GB', times = 31), rep('R1', times = 31), rep('R0', times = 31)))
Grid_plot2_2 = ggplot(aes(x = x, y = y, ymin = ymin, ymax = ymax, linetype=group), data = Grid_df2_2) +
  geom_ribbon(aes(colour = group, fill = group), alpha = 0, size = 1) +
  scale_linetype_manual(values=c('solid',"dashed", "dotted")) +
  xlab('year') + ylab('posterior mean bs on original scale') +
  ggtitle('')
Grid_plot2_2

# Plot sites that were consistent vs sites that were dropped #
consistent_df2 = data.frame(y = exp(c(consistent_results2[,1],inconsistent_results2[,1]))*30.7,
                      x = rep(66:96, times = 2),
                      ymin = exp(c(consistent_results2[,3], inconsistent_results2[,3]))*30.7,
                      ymax = exp(c(consistent_results2[,4], inconsistent_results2[,4]))*30.7,
                      group = c(rep('consistent', times = 31), rep('inconsistent', times = 31)))
consistent_plot2 = ggplot(aes(x = x, y = y, ymin = ymin, ymax = ymax), data = consistent_df2) +
  geom_ribbon(aes(colour = group, alpha = 0.05, fill = group)) +
  xlab('year') + ylab('posterior mean bs on original scale') +
  ggtitle('Mod2 posterior mean BS with 95% credible intervals of the consistent (not dropped) vs. inconsistent sites.')
consistent_plot2

# Plot the AR1 time effect for R process #
AR_df_2 = matrix(0, nrow = m_samples, ncol = no_T)
for (i in 1:m_samples)
{
  AR_df_2[i,] = samp2[[i]]$latent[grep('R_year',names(samp2[[i]]$latent), fixed = TRUE)[1:no_T]] # AR1 zero mean time effect
}
AR_df2_2 = data.frame(mean = apply(AR_df_2, 2, mean, na.rm=T),
                    LCL = apply(AR_df_2, 2, quantile, probs = 0.025, na.rm=T),
                    UCL = apply(AR_df_2, 2, quantile, probs = 0.975, na.rm=T),
                    x = 66:96)

AR_plot2 = ggplot(data = AR_df2_2, aes(x = x, y = mean, ymax = UCL, ymin = LCL)) +
  geom_ribbon(alpha = 0.3, colour = 'blue') +
  geom_abline(intercept = 0, slope = 0) +
  xlab('year') + ylab('posterior mean increase in the average logit') +
  ggtitle('A plot of the AR(1) process in the R process in mod2 reflecting public mood')
AR_plot2
#

###################### Diagnostics #############################

# First check the normality assumptions of the residuals, the random intercepts and the random slopes #
# Also plot the residuals by year to assess temporal correlation #

RI_plot_df = data.frame(y = RI_array[,1])
RS_plot_df = data.frame(y = RS_array[,1])
residual_plot_df = data.frame(x = rep(66:96,each = dim(BlackSmokePrefData)[1]), y = as.vector(residuals_array))

par(mfrow = c(2,2))
hist(RI_plot_df$y, main = 'Histogram of random site-specific intercepts')
qqnorm(RI_plot_df$y)
qqline(RI_plot_df$y)
hist(RS_plot_df$y, main = 'Histogram of random site-specific slopes')
qqnorm(RS_plot_df$y)
qqline(RS_plot_df$y)
par(mfrow = c(1,1))

residuals_plot = ggplot(data = residual_plot_df, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth() +
  xlab('year') + ylab('residual') + ggtitle('A plot of the residuals vs year')
residuals_plot

par(mfrow=c(1,1))
qqPlot(residual_plot_df$y, main = 'Normal QQ-plot of residuals', ylab = 'observed quantiles')
qqline(residual_plot_df$y)

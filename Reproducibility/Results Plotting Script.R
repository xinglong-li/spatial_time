# This is an R script for generating all the plots seen in Watson et al 2021

# Load the packages
library(raster)
require("rgdal")
require("rgeos")
require("dplyr")
library(INLA)
library(mvtnorm)
library(boot)
library(reshape2)
library(sp)
library(ggplot2)
library(MASS)

# Load in all the data used to fit the models
load('./Reproducibility/final_mod_joint1_helper2.RData')
# load a higher resolution mesh
mesh <- readRDS('./Reproducibility/mesh_5_7.rds')
# the stack_X_expand objects contain the dummy points situated throughout GB
# the stack_X_est objects are restricted to the sampled locations
# The stack_joint1 object contains only the X_'est' data objects and does not
# adjust for PS fully (corresponds to population 1 in the paper)

# We can quickly see this by running the following
table(stack_joint1$data$data$R)
table(is.na(stack_joint1$data$data$y)) # Note the FALSE match the number of R=1
table(stack_R_est$data$data$R) # matches

# We can instead look to the P2 population
table(stack_R_est_expand$data$data$R) # Notice the much higher number of 0's
table(is.na(stack_y_est_expand$data$data$y)) # Notice the match

# Load the P2 model corresponding to the previously loaded mesh
mod1 <- readRDS('./Reproducibility/Mod3_V5_superhighresmesh_regular_eblaplace.rds')
# Or load the model fit using the script Fit Models.R
#mod1 <- readRDS('finalmod_joint2_reproduce.rds')


############################################
### Sampling from Posterior distribution ###
############################################
# Creating Samples
m_samples = 200
# The paper uses 1000, but we suggest running only 100 to save memory

# Do we want to sample from the model now, or load precompiled samples?
# Note the samples come from the model saved in file: Mod3_V5_superhighresmesh_regular_eblaplace.rds
sample <- T

if(sample)
{
  samp2 = inla.posterior.sample(mod1, n=m_samples, selection=list(Apredictor=1, Predictor=1))
}
if(!sample)
{
  # To save time, we load our samples
  samp2 = readRDS('./Reproducibility/Mod3_V5_samples_superhighresmesh_regular_eblaplace1.rds')
}

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
points <- SpatialPoints(coords = proj_grid$lattice$loc*sd_x,
                        proj4string = CRS(proj4string(region_GB)))
table(xy.in <- gContains(region_GB, points, byid=T))

# Creating the A projector Matrix 
A.grid2 = inla.spde.make.A(mesh,
                           loc = cbind(BlackSmokePrefData$east, BlackSmokePrefData$north))

# 45466 Apredictors and 50252 Predictors
tmp2 = matrix(0, nrow = no_T, ncol = mesh$n)
year = (1:no_T)/no_T
year2 = year^2
Post_Array2 = array(0, dim = c(m_samples, no_T, mesh$n))
Post_Sites_Array2 = array(0, dim = c(m_samples, no_T, dim(BlackSmokePrefData)[1]))
residuals_array = array(0, dim = c(m_samples, no_T, dim(BlackSmokePrefData)[1]))
RI_array = array(0, dim = c(m_samples, no_T, dim(BlackSmokePrefData)[1])) # random intercepts
RS_array = array(0, dim = c(m_samples, no_T, dim(BlackSmokePrefData)[1])) # random slopes

# Looping for each sample 
for (i in 1:m_samples){
  # Creating predictions for Y process
  print(i)
  counter = 1
  for(j in 1:no_T) # loop over the years
  {
    # If old naming scheme
    if(sum(names(samp2[[2]]$latent)=='Intercept')>0)
    {
      Post_Array2[i,j,] <- samp2[[i]]$latent['Intercept'] + # Intercept 
        samp2[[i]]$latent['year'] * year[j] + # linear fixed effect
        samp2[[i]]$latent['year_2'] * year2[j] + # quadratic fixed effect
        as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial.field:',x = names(samp2[[i]]$latent)))]) + # Spatial Intercept
        as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial.field.copy:',x = names(samp2[[i]]$latent)))]) * year[j] + # Spatial linear slope
        as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial.field.copy2:',x = names(samp2[[i]]$latent)))]) * year2[j] #+ # Spatial quadratic slope
      #samp[[i]]$latent[grep('year:',rownames(samp[[i]]$latent), fixed = TRUE)[j],1] # AR1 zero mean time effect
      
      # Project onto the site locations and add random site-specific effects #
      Post_Sites_Array2[i,j,] = as.numeric(A.grid2 %*% Post_Array2[i,j,]) +
        as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial_ind:',x = names(samp2[[i]]$latent)))[1:no_sites]]) + # random Intercept
        as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial_ind2:',x = names(samp2[[i]]$latent)))[no_sites+(1:no_sites)]]) * year[j] # random slope
      
      # Extract the random intercepts #
      RI_array[i,j,] = as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial_ind:',x = names(samp2[[i]]$latent)))[1:no_sites]])
      
      # Extract the random slopes #
      RS_array[i,j,] = as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial_ind2:',x = names(samp2[[i]]$latent)))[no_sites+(1:no_sites)]])
      
      # Extract the residuals #
      residuals_array[i,j,] = log(BlackSmokePrefData[,(3+j)] / mean(colMeans(BlackSmokePrefData[,4:34], na.rm = T))) - Post_Sites_Array2[i,j,]
    }
    # If new naming scheme
    if(sum(rownames(samp2[[2]]$latent)=='Intercept:1')>0)
    {
      Post_Array2[i,j,] <- samp2[[i]]$latent[which(rownames(samp2[[i]]$latent)=='Intercept:1')] + # Intercept 
        samp2[[i]]$latent[which(rownames(samp2[[i]]$latent)=='year:1')] * year[j] + # linear fixed effect
        samp2[[i]]$latent[which(rownames(samp2[[i]]$latent)=='year_2:1')] * year2[j] + # quadratic fixed effect
        as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial.field:',x = rownames(samp2[[i]]$latent)))]) + # Spatial Intercept
        as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial.field.copy:',x = rownames(samp2[[i]]$latent)))]) * year[j] + # Spatial linear slope
        as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial.field.copy2:',x = rownames(samp2[[i]]$latent)))]) * year2[j] #+ # Spatial quadratic slope
      #samp[[i]]$latent[grep('year:',rownames(samp[[i]]$latent), fixed = TRUE)[j],1] # AR1 zero mean time effect
      
      # Project onto the site locations and add random site-specific effects #
      Post_Sites_Array2[i,j,] = as.numeric(A.grid2 %*% Post_Array2[i,j,]) +
        as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial_ind:',x = rownames(samp2[[i]]$latent)))[1:no_sites]]) + # random Intercept
        as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial_ind2:',x = rownames(samp2[[i]]$latent)))[no_sites+(1:no_sites)]]) * year[j] # random slope
      
      # Extract the random intercepts #
      RI_array[i,j,] = as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial_ind:',x = rownames(samp2[[i]]$latent)))[1:no_sites]])
      
      # Extract the random slopes #
      RS_array[i,j,] = as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial_ind2:',x = rownames(samp2[[i]]$latent)))[no_sites+(1:no_sites)]])
      
      # Extract the residuals #
      residuals_array[i,j,] = log(BlackSmokePrefData[,(3+j)] / mean(colMeans(BlackSmokePrefData[,4:34], na.rm = T))) - Post_Sites_Array2[i,j,]
    }

    counter = counter + mesh$n
  }
}

# Summarising predictions 
Post_mean2 <- t(apply(Post_Array2, c(2,3), function(x) {mean(x, na.rm=TRUE)}))
Post_sd2 <- t(apply(Post_Array2, c(2,3), function(x) {sd(x, na.rm=TRUE)}))
Post_LCL2 <- apply((apply(Post_Array2, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.025))
Post_UCL2 <- apply((apply(Post_Array2, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.975))
RI_array <- t(apply(RI_array, c(2,3), function(x) {mean(x, na.rm=TRUE)}))
RS_array <- t(apply(RS_array, c(2,3), function(x) {mean(x, na.rm=TRUE)}))
residuals_array <- t(apply(residuals_array, c(2,3), function(x) {mean(x, na.rm=TRUE)}))

# Project onto the grid using inla.projector 
post_Matrix_grid_mean2 = array(NA, dim = c(nxy[1], nxy[2], no_T))
post_Matrix_grid_sd2 = array(NA, dim = c(nxy[1], nxy[2], no_T))

for(i in 1:no_T)
{
  post_Matrix_grid_mean2[,,i] = inla.mesh.project(proj_grid, Post_mean2[,i])
  post_Matrix_grid_mean2[,,i][!xy.in] <- NA
  post_Matrix_grid_sd2[,,i] = inla.mesh.project(proj_grid, Post_sd2[,i])
  post_Matrix_grid_sd2[,,i][!xy.in] <- NA
}

par(mfrow = c(2,2))
image(post_Matrix_grid_mean2[,,1], asp=1, main = 'Plot of the expected pollution in 1966 (red is low)')
image(post_Matrix_grid_mean2[,,30], asp=1, main = 'Plot of the expected pollution in 1996')
image(post_Matrix_grid_sd2[,,1], asp=1, main = 'Plot of the se of the expected pollution in 1966 (red is low)')
image(post_Matrix_grid_sd2[,,30], asp=1, main = 'Plot of the se of the expected pollution in 1996') 

# # Step one - Make results tables for fixed parameter estimates
fixed_results = as.data.frame(mod1$summary.fixed)

# # Step two - Make results tables for hyperparameter estimates
hyper_results = as.data.frame(mod1$summary.hyperpar)
# In the paper we used the Monte Carlo samples to get the quantities
# It appears that some of the sampling was done in error in an old version of INLA 
# (possibly due to empirical Bayes being used) leading to incorrect posterior SD of some hyperparameters
# The conclusions don't change qualitatively
# The sampling seems to have been fixed now in the later versions of INLA
if(sample)
{
  mod1_hyperparameter_samples <-
    inla.hyperpar.sample(mod1, n=m_samples) 
  # Per the author, Havard Rue, "If the task is to sample the hyperparameters, 
  # the you can use inla.hyperpar.sample(), which has a different behaviour than inla.posterior.sample
  
}
if(!sample)
{
  # To save time, we load our samples
  mod1_hyperparameter_samples <- readRDS("EB_Mod3_V5_hyperpar_samples.rds")
  # Notice the posterior SD estimates - they appear to low
}
apply(mod1_hyperparameter_samples, 2, mean)
apply(mod1_hyperparameter_samples, 2, sd)

# # Step three - Make results tables for selected (R=1) sites estimates at each t
site_results2 = matrix(0, nrow = no_T, ncol = 4) # prediction means, sd, LCL, UCL
colnames(site_results2) = c('predciction mean', 'prediction sd', 'LCL', 'UCL')

R1_results2 = matrix(NA, nrow = no_T, ncol = 4) # prediction means, sd, LCL, UCL
colnames(R1_results2) = c('predciction mean', 'prediction sd', 'LCL', 'UCL')

# # Step four - Make results tables for unselected (R=0) sites estimates at each t
R0_results2 = matrix(NA, nrow = no_T, ncol = 4) # prediction means, sd, LCL, UCL
colnames(R0_results2) = c('predciction mean', 'prediction sd','LCL', 'UCL')

# Step 5 - Make results table for consistent sites vs those that are dropped
# library(plyr)
# temp_df = ddply(BlackSmokePrefData2, .(site), summarize, nR = sum(R))
# consistent_ind = which(temp_df$nR == 31)
# inconsistent_ind = which(temp_df$nR != 31)
# 
# consistent_results2 = matrix(NA, nrow = no_T, ncol = 4) # prediction means, sd, LCL, UCL
# colnames(consistent_results2) = c('predciction mean', 'prediction sd', 'LCL', 'UCL')
# 
# inconsistent_results2 = matrix(NA, nrow = no_T, ncol = 4) # prediction means, sd, LCL, UCL
# colnames(inconsistent_results2) = c('predciction mean', 'prediction sd', 'LCL', 'UCL')

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
Grid_plot_2 = ggplot(aes(x = x, y = y, ymin = ymin, ymax = ymax, alpha = 0.05), data = Grid_df_2) +
  geom_ribbon(aes(colour = group, alpha = 0.05, fill = group)) +
  xlab('year') + ylab('posterior mean bs on transformed scale') +
  ggtitle('Mod2 posterior means of BS at the selected, unselected sites and across the UK (R1, R0 and Whole UK resp.). ')
Grid_plot_2

# Improved plot style
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

# On the original scale 
Grid_df2_2 = data.frame(y = exp(c(Grid_results2[,1],R1_results2[,1],R0_results2[,1]))*30.7,
                        x = rep(66:96, times = 3),
                        ymin = exp(c(Post_LCL2, R1_results2[,3], R0_results2[,3]))*30.7,
                        ymax = exp(c(Post_UCL2, R1_results2[,4], R0_results2[,4]))*30.7,
                        group = c(rep('Whole UK', times = 31), rep('R1', times = 31), rep('R0', times = 31)))
Grid_plot2_2 = ggplot(aes(x = x, y = y, ymin = ymin, ymax = ymax, alpha = 0.05), data = Grid_df2_2) +
  geom_ribbon(aes(colour = group, fill = group)) + 
  xlab('year') + ylab('posterior mean bs on original scale') +
  ggtitle('Posterior means of BS at the selected, unselected sites and across the UK (R1, R0 and Whole UK resp.). ')
Grid_plot2_2

# Plot sites that were consistent vs sites that were dropped #
# consistent_df2 = data.frame(y = exp(c(consistent_results2[,1],inconsistent_results2[,1]))*30.7,
#                             x = rep(66:96, times = 2),
#                             ymin = exp(c(consistent_results2[,3], inconsistent_results2[,3]))*30.7,
#                             ymax = exp(c(consistent_results2[,4], inconsistent_results2[,4]))*30.7,
#                             group = c(rep('consistent', times = 31), rep('inconsistent', times = 31)))
# consistent_plot2 = ggplot(aes(x = x, y = y, ymin = ymin, ymax = ymax), data = consistent_df2) +
#   geom_ribbon(aes(colour = group, alpha = 0.05, fill = group)) +
#   xlab('year') + ylab('posterior mean bs on original scale') +
#   ggtitle('Mod2 posterior mean BS with 95% credible intervals of the consistent (not dropped) vs. inconsistent sites.')
# consistent_plot2

# Plot the AR1 time effect for R process #
AR_df_2 = matrix(0, nrow = m_samples, ncol = no_T)
for (i in 1:m_samples)
{
  # If old naming scheme
  if(sum(names(samp2[[2]]$latent)=='Intercept')>0)
  {
    AR_df_2[i,] = samp2[[i]]$latent[grep('R_year',names(samp2[[i]]$latent), fixed = TRUE)[1:no_T]] # AR1 zero mean time effect
  }
  # If new naming scheme
  if(sum(rownames(samp2[[2]]$latent)=='Intercept:1')>0)
  {
    AR_df_2[i,] = samp2[[i]]$latent[grep('R_year',rownames(samp2[[i]]$latent), fixed = TRUE)[1:no_T]] # AR1 zero mean time effect
  }
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

qqnorm(residual_plot_df$y, main = 'Normal QQ-plot of residuals')
qqline(residual_plot_df$y)

## Now we plot the population density curves

# create the population density raster
pop_dens = raster('./Reproducibility/UK_residential_population_2011_1_km.asc')
pop_dens = as(pop_dens, 'SpatialPixelsDataFrame')
pop_dens@proj4string = GB_polygons_simp@proj4string
pop_dens@data$UK_residential_population_2011_1_km[is.na(pop_dens@data$UK_residential_population_2011_1_km)] = 0
pop_dens = pop_dens[!is.na(sp::over(pop_dens, GB_polygons_simp)),]

# create the projection matrix onto the population density raster
A_popden = inla.spde.make.A(mesh, loc = pop_dens@coords/sd_x)

# sample random multivariate normal values for the 1km x 1km pixels of pop_dens
var_b1 = 1 / ( mod1$summary.hyperpar['Precision for spatial_ind (component 1)','mean'] )
var_b2 = 1 / ( mod1$summary.hyperpar['Precision for spatial_ind (component 2)','mean'] )
rho_b = mod1$summary.hyperpar['Rho1:2 for spatial_ind','mean']
Sigma_b = matrix(c(var_b1, rho_b, rho_b, var_b2), byrow = T, nrow=2, ncol=2)

# Create objects to project onto a regular lattice grid covering GB
# Project onto the grid using inla.projector 
post_Matrix_grid = array(NA, dim = c(m_samples, nxy[1], nxy[2], no_T))

# Create objects to project onto the population density pixels grid covering GB
# Project onto the grid using inla.projector 
# first create a matrix to store the average annual exposures per captia
post_popdens_avg = array(NA, dim = c(m_samples, no_T))

# next look at proportion of individuals who are exposed to bs levels above the EU standard
post_popdens_prop = array(NA, dim = c(m_samples, no_T))

pop_dens_data = as.numeric(pop_dens@data$UK_residential_population_2011_1_km)
# normalize to sum to 1 so that we form an average
pop_dens_data = pop_dens_data / sum(pop_dens_data)
n_pop_dens = length(pop_dens_data)

# This is the EU exceedance value of interest during the time period on the log scale
exceedance_val = log(34 / mean(colMeans(BlackSmokePrefData[,4:34], na.rm = T)))

set.seed(12345)
#beta_gamma = rep(0.6240, times = 1000)   #rnorm(m_samples, mean = 1.974, sd = 0.1570)
#beta_b = rep(0.0640, times = 1000)  #rnorm(m_samples, mean = 0.0894, sd = 0.0408)
#nugget = 1/rep(19.2331, times = 1000)  #
# Looping for each sample 
for (i in 1:m_samples){
  # Creating predictions for Y process
  print(i)
  print('from')
  print(m_samples)
  counter = 0
  counter2 = 0
  
  # sample random multivariate normal values for the 1km x 1km pixels of pop_dens
  RI_s_pop_dens = mvrnorm(n = n_pop_dens, mu = c(0,0), Sigma = Sigma_b)
  
  for(j in 1:no_T) # loop over the years
  {
    #print(j)
    
    # Post_Array2[i,j,] <- samp2[[i]]$latent['Intercept'] + # Intercept 
    #    samp2[[i]]$latent['year'] * year[j] + # linear fixed effect
    #    samp2[[i]]$latent['year_2'] * year2[j] + # quadratic fixed effect
    #    as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial.field:',x = names(samp2[[i]]$latent)))]) + # Spatial Intercept
    #    as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial.field.copy:',x = names(samp2[[i]]$latent)))]) * year[j] + # Spatial linear slope
    #    as.numeric(samp2[[i]]$latent[which(startsWith(prefix = 'spatial.field.copy2:',x = names(samp2[[i]]$latent)))]) * year2[j] #+ # Spatial quadratic slope

    # map to the grid
    post_Matrix_grid[i,,,j] = inla.mesh.project(proj_grid, Post_Array2[i,j,])
    post_Matrix_grid[i,,,j][!xy.in] <- NA
    
    # create the linear predictor addon for the pop_dens values
    linear_pred_addon = RI_s_pop_dens[,1] + RI_s_pop_dens[,2] * year[j]
    
    # compute population average exposures
    post_popdens_avg[i,j] = sum( ( (A_popden %*% Post_Array2[i,j,]) + linear_pred_addon ) * pop_dens_data )
    
    # compute the proportion of population exposed to levels of BS above EU threshold
    post_popdens_prop[i,j] = sum( ( ( (A_popden %*% Post_Array2[i,j,]) + linear_pred_addon ) > exceedance_val) * pop_dens_data )
    
    counter = counter + mesh$n
    counter2 = counter2 + no_sites
  }
}

# Project onto the grid using inla.projector 
post_Matrix_grid_mean2 = apply(post_Matrix_grid, c(2,3,4), function(x) {mean(x, na.rm=TRUE)})
post_Matrix_grid_sd2 = apply(post_Matrix_grid, c(2,3,4), function(x) {sd(x, na.rm=TRUE)})

# compute the exceedance value on our scale
exceedance_val = log(34 / mean(colMeans(BlackSmokePrefData[,4:34], na.rm = T)))
post_Matrix_grid_exceedance = apply(post_Matrix_grid, c(2,3,4), function(x) {mean(x>exceedance_val, na.rm = T)})

post_grid_exceedance_annual_proportion_means = apply(post_Matrix_grid, c(1,4), function(x) {mean(x>exceedance_val, na.rm = T)})
post_grid_exceedance_annual_proportion_mean = apply(post_grid_exceedance_annual_proportion_means, c(2), mean, na.rm=T)
post_grid_exceedance_annual_proportion_UCL95 = apply(post_grid_exceedance_annual_proportion_means, c(2), function(x){quantile(x, probs = c(0.975), na.rm=T)})
post_grid_exceedance_annual_proportion_LCL95 = apply(post_grid_exceedance_annual_proportion_means, c(2), function(x){quantile(x, probs = c(0.025), na.rm=T)})
post_grid_exceedance_annual_proportion_SD = apply(post_grid_exceedance_annual_proportion_means, c(2), sd, na.rm=T)

# now for population density
post_popdens_mean_annual_exposure = apply(post_popdens_avg, 2, mean, na.rm=T)
post_popdens_UCL95_annual_exposure = apply(post_popdens_avg, 2, function(x){quantile(x, probs = c(0.975), na.rm=T)})
post_popdens_LCL95_annual_exposure = apply(post_popdens_avg, 2, function(x){quantile(x, probs = c(0.025), na.rm=T)})
post_popdens_SD_annual_exposure = apply(post_popdens_avg, 2, sd, na.rm=T)

# what proportion of GB's population has BS exposure levels exceeding EU limits
post_popdens_mean_annual_exceedance = apply(post_popdens_prop, 2, mean, na.rm=T)
post_popdens_UCL95_annual_exceedance = apply(post_popdens_prop, 2, function(x){quantile(x, probs = c(0.975), na.rm=T)})
post_popdens_LCL95_annual_exceedance = apply(post_popdens_prop, 2, function(x){quantile(x, probs = c(0.025), na.rm=T)})
post_popdens_SD_annual_exceedance = apply(post_popdens_prop, 2, sd, na.rm=T)

rm(samp2)

results = list(post_Matrix_grid_mean2 = post_Matrix_grid_mean2,
               post_Matrix_grid_sd2 = post_Matrix_grid_sd2,
               exceedance_val = exceedance_val,
               post_Matrix_grid_exceedance = post_Matrix_grid_exceedance,
               post_grid_exceedance_annual_proportion = data.frame(Post_mean = post_grid_exceedance_annual_proportion_mean,
                                                                   Post_SD = post_grid_exceedance_annual_proportion_SD,
                                                                   Post_LCL95 = post_grid_exceedance_annual_proportion_LCL95,
                                                                   Post_UCL95 = post_grid_exceedance_annual_proportion_UCL95),
               pop_dens_annual_avg = data.frame(Post_mean = post_popdens_mean_annual_exposure,
                                                Post_SD = post_popdens_SD_annual_exposure,
                                                Post_LCL95 = post_popdens_LCL95_annual_exposure,
                                                Post_UCL95 = post_popdens_UCL95_annual_exposure),
               pop_dens_annual_exceedanceprob = data.frame(Post_mean = post_popdens_mean_annual_exceedance,
                                                           Post_SD = post_popdens_SD_annual_exceedance,
                                                           Post_LCL95 = post_popdens_LCL95_annual_exceedance,
                                                           Post_UCL95 = post_popdens_UCL95_annual_exceedance))


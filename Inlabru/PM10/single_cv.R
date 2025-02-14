library(dplyr)
library(sp)
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
PM10s0$annual_mean <-
  log(PM10s0$annual_mean) - log(mean_pm_annually)

no_sites = length(unique(PM10s0$site_number))
no_T = length(unique(PM10s0$year))

# Reform the data so that each site has num_of_year rows, the empty records are filled with NAs
PM10s_flat <-
  dcast(PM10s0, site_number + north + east ~ year, value.var = "annual_mean")
PM10s = melt(PM10s_flat, id.vars = c(1, 2, 3), variable.name = 'year', value.name = 'annual_mean')
PM10s$year <- as.numeric(as.character(factor(PM10s$year, labels = 1985:2022)))

stopifnot(dim(PM10s)[1] == no_sites * no_T)

subsmaple_plotting_data = PM10s[which(PM10s$site_number %in% unique(PM10s$site_number)[sample.int(111, size = 30)]), ]
means_plot = ggplot(data = subsmaple_plotting_data, aes(x = year, y = annual_mean)) +
  geom_line(aes(group = site_number, colour = site_number)) + xlab("Year") + ylab("log(PM10)")
means_plot

var_annual <- group_by(PM10s, year) %>%
  summarise(var_pm = var(annual_mean, na.rm = T))
variance_plot =  ggplot(data = var_annual, aes(x = year, y = var_pm)) +
  geom_line() + geom_smooth() + xlab("Year") + ylab("Variance log(PM10)")
ggtitle('A plot of the variance of the log annual means with fitted smoother')
variance_plot

PM10s$time <-
  (PM10s$year - min(PM10s$year)) / (max(PM10s$year) - min(PM10s$year))
PM10s$locs <- coordinates(PM10s[, c("east", "north")])
PM10s$site_number <- as.numeric(as.factor(PM10s$site_number))

# Use K-fold cross validation ======================================================================

cross_validation <- function(K, cutoff_dist, idx_sites, hyper_ini=NULL){
  # K          : K-fold cross validation.
  # cutoff_dist: maximum innner cut-off distance of the mesh grid.
  # idx_sites  : random permumation of indeces of sites, used to split K subsets.
  # hyper_ini  : list of hyper parameters, one for each of the K subsets. The hyper parameters are 
  #              posterior modes of the hyperparameters fitted using coerse grid on the same subset. 
  
  bru_options_set(bru_max_iter = 1,
                  control.inla = list(strategy = "gaussian", int.strategy = 'eb'))
  
  cutoff_outer <- 2 * cutoff_dist
  
  mesh <- fm_mesh_2d_inla(loc = cbind(PM10s$east, PM10s$north),
                          boundary = CA_border,
                          offset = c(0.1, 0.2), 
                          max.edge = c(cutoff_dist, cutoff_outer),
                          cutoff = cutoff_dist,
                          min.angle = 26)
  
  spde_obj <- inla.spde2.pcmatern(mesh = mesh, 
                                  alpha = 2, 
                                  prior.range = c(0.1, 0.01),
                                  prior.sigma = c(1, 0.01),
                                  constr = T)
  
  comp <- annual_mean ~ Intercept(1) + Time_1(time) + Time_2(time^2) +
    Random_0(site_number, model = "iid2d", n = no_sites*2, constr=TRUE) + 
    Random_1(site_number, weights = time, copy = "Random_0") +
    Spatial_0(locs, model = spde_obj) + 
    Spatial_1(locs, weights = time, model = spde_obj) + 
    Spatial_2(locs, weights = time^2, model = spde_obj)

  no_site_per_group <- no_sites %/% K
  cv_errs <- NULL
  hyper_modes <- NULL
  
  for(k in 1:K){
    # Split the dataset into training and prediction set.
    print(sprintf("Subset %s out of %s", k, K))
    if(k == K){
      i_pred <- (1+(K-1)*no_site_per_group) : no_sites
      site_pred <- idx_sites[i_pred]
      site_fit <- idx_sites[-i_pred]
    }else{
      i_pred <- (1+(k-1)*no_site_per_group) : (k*no_site_per_group)
      site_pred <- idx_sites[i_pred]
      site_fit <- idx_sites[-i_pred]
    }
    
    data_fit <- filter(PM10s, site_number %in% site_fit)
    data_pred <- filter(PM10s, site_number %in% site_pred)
    
    # Initialize the model using the fit of a model with coarse grid
    if(!is.null(hyper_ini)){
      bru_options_set(control.mode = list(theta = hyper_ini[[k]], restart = TRUE))
    }
    
    # Fit model
    fit_bru <- bru(comp, family = "gaussian", data = data_fit)

    # Save the posterior modes of hyper parameters for later initialization
    hyper_modes <- c(hyper_modes, list(fit_bru$mode$theta))
    
    # Predict
    pred_bru <- predict(fit_bru, 
                        data_pred, 
                        ~ exp(Intercept + Time_1 + Time_2 + Random_0 + Random_1 + Spatial_0 + Spatial_1 + Spatial_2), 
                        n.samples = 1000)
    cv_err_k <- mean(pred_bru$mean - data_pred$annual_mean, na.rm = TRUE)
    cv_errs <- c(cv_errs, cv_err_k)
  }
  list(pred_err = cv_errs, 
       post_modes = hyper_modes)
}


mesh_select <- function(K, cutoff_dist_0){
  # cutoff_dist_0: starting value of cut-off distance for inner grid
  
  cutoff_dist <- cutoff_dist_0
  idx_site_permute <- sample(no_sites, no_sites)
  
  hyper_ini <- NULL
  cv_errs <- NULL
  t <- 0
  while(TRUE){
    print(sprintf("Cutoff distance is %.3f", cutoff_dist))
    t <- t + 1
    cv_fit <- cross_validation(K, cutoff_dist, idx_site_permute, hyper_ini)
    cv_errs <- c(cv_errs, mean(cv_fit$pred_err))
    print(sprintf("The cross validation prediction error is %.3f", mean(cv_fit$pred_err)))

    hyper_ini <- cv_fit$post_modes
    cutoff_dist <- cutoff_dist * 0.95
    
    if(t > 2 && cv_errs[t] > cv_errs[t-1] && cv_errs[t-1] > cv_errs[t-2]) {
      print("CV error starts increasing. STOP")
    }
  }
  cv_errs
}

# Train ============================================================================================

cutoff_dist_0 <- 0.3
K <- 20

set.seed(2)
mesh_select(K, cutoff_dist_0)

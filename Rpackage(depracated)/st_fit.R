library(INLA)


spatio_temporal <- function(formula_observation,
                            data,
                            sites,
                            locations,
                            times = NULL,
                            boundary = NULL,
                            formula_selection = NULL,
                            cov_observation = NULL,
                            cov_selection = NULL,
                            family = "gaussian") {
  y ~ x1 + x2 + (1 | x3) + (1 + x1 | x3)
  cov_observation = ["iid", "dependent", "Mattern"]
  
  # Args:
  #   formula_observation: Formula for observation model. Use the syntax of the package lmer.
  #   formula_selection:   Formula for site selection model. Use the syntax of the package lmer.
  #   cov_observation:     List of covariance matrix for random effects in the observation model.
  #                        One for each set of random effects specified in formula_observation.
  #                        If 'Mattern' used, then spatial SPDE mattern covariance for the 
  #                        complete mesh built for INLA (not only on observed sites) is used.
  #                        If not speficied, the diagonal covariance matrix is used for each
  #                        set of random effects.
  #   cov_selection:       List of covariance matrix for random effects in the site selection model.
  #                        One for each set of random effects specified in formula_selection.
  #                        If 'Mattern' used, then spatial SPDE mattern covariance for the 
  #                        complete mesh built for INLA (not only on observed sites) is used.
  #                        If not speficied, the diagonal covariance matrix is used for each
  #   data:                The data frame contains observation variable and covariates.
  #                        Each row is the observation of one site at one time point.
  #   sites:               The site names used as indicator of each site. 
  #                        Specified by the name of the column in data.
  #   times:               Time points used as indicator of the time points of each observation.
  #                        Specified by the name of the column in data.
  #                        If not specified, the whole observation is treated as from one time point,
  #                        and only a spatial model (no temporal trend) can be fitted.
  #   locations:           Locations of sites. It can be either a SpatialPoints object from 
  #                        package 'sp', or a data frame with each row contains site_name, loc_N, loc_E.
  #   boundary:            Boundary of sites, must be a SpatialPolygons object from the package 'sp'.
  #   family:              String that specifies the distribution of the observation.
  #                        The default family for site selection is always 'Bernoulli'.
  #                        Check ```names(inla.models()$likelihood)``` for the full list of 
  #                        acceptable families.
  
  # Prepare standardized variables/factors and for INLA formula ------------------------------------

  data_formated <- format_data(data, sites, times)
  locations_formated <- format_locations(locations)
  boudary_formated <- format_boundary(boundary)
  formula_obs <- parse_formula(formula_observation)
  formula_slc <- parse_formula(formula_selection)
  
  model <- build_model(formula_obs,
                       formula_slc,
                       data_formated,
                       locations_formated,
                       boundary_formated)
  
  # Model fitting ----------------------------------------------------------------------------------
  
  inla_obj <- inla(formula = model$formula,
                   family = family,
                   data = inla.stack.data(model$stack_joint,
                                          spde = st_spde),
                   control.predictor = list(A = inla.stack.A(stack_joint),
                                            compute = TRUE),
                   control.family = list(list(hyper=precprior),
                                         list(hyper=precprior)),
                   control.compute = list(dic = TRUE, 
                                          config = TRUE,
                                          cpo = TRUE),
                   control.fixed = list(mean = list(Intercept = 1, default = 0),
                                        prec = list(Intercept = 0.25, default = 0.001)),
                   control.results = list(return.marginals.random = FALSE,
                                          return.marginals.predictor = FALSE),
                   control.inla = list(h = 1e-5, 
                                       strategy = "laplace",
                                       int.strategy = 'eb'),
                   verbose = verbose, 
                   num.threads = 8)
  
  # Return results ---------------------------------------------------------------------------------
  
  st_obj <- list(inla_mesh = st_mesh,
                 inla_obj = inla_obj)
  class(st_obj) <- "spatial_temporal"

  return(st_obj)
}



summary.spatial_temporal <- function(obj) {
  # S3 method for printing a summary of analysis results.
  #
  # Args:
  #  obj: 
}


print.spatial_temporal <- function(obj) {
  # S3 method for printing a summary of analysis results.
  #
  # Args:
}

library(inla)

build_model <- function(formula_obs,
                        formula_slc,
                        data_formated,
                        locations_formated,
                        boundary_formated) {
  # Construct the INLA model.
  #
  # Args:
  #   formula_obs: 
  #
  # Returns:
  #  
  
  # Create components.
  
  model_obs <- build_observation_model(formula_obs, 
                                       data_formated)
  model_slc <- build_selection_model(formula_slc)
  
  formula_joint <- model_obs$formula + model_slc$formula
  stack_joint <- inla.stack(model_obs$stack,
                            model_slc$stack)
  model_joint <- list(formula = formula_joint,
                      stack = stack_joint)
}


build_observation_model <- function() {
  
}


build_selection_model <- function() {
  
}








build_model <- function() {
  # Prepare the mesh grid --------------------------------------------------------------------------
  
  st_mesh <- inla.mesh.2d(loc = coodinates,
                          boundary = boundary,
                          max.edge = c(),
                          min.angle = NULL,
                          offset = c(),
                          cutoff = c()
  )
  
  # Construct the SPDE -----------------------------------------------------------------------------
  
  st_spde <- inla.spde2.pcmatern(mesh = st_mesh,
                                 alpha = 2,
                                 prior.range = c(0.3, 0.5), # P(practic.range < 0.3) = 0.5
                                 prior.sigma = c(10, 0.01) # P(sigma > 10) = 0.01
  )
  
  # Extract the projection matrix for estimation ---------------------------------------------------
  
  A_est <- inla.spde.make.A(mesh = st_mesh,
                            loc = loc
  )
  
  grid_index <- inla.spde.make.index(name = "spatial_field",
                                     n.spde = st_spde$n.spde
  )
  
  if (temporal) {
    A_est <- inla.make.A(mesh = st_mesh,
                         loc = coordinates,
                         group = time_points,
                         n.group = n_days)
    grid_index <- inla.spde.make.index(name = "spatial_field",
                                       n.spde = st_mesh$n.spde,
                                       n.groups = n_days)
  }
  
  # Prepare the stack for estimation ---------------------------------------------------------------
  
  A_stack_est <- list(A_est, 1) # we can put all covariates into one data frame and share the '1'
  effects_stack_est <- list(c(grid_index, list(Intercept = 1)),
                            list(covariates)
  )
  stack_est <- inla.stack(data = list(y = NULL),
                          A = A_stack_est,
                          effects = effects_stack_est,
                          tag = "estimation")
  
  if (multiplelikelihoods) {
    stack_response1 <- inla.stack(data = list(y = cbind(x, NA),
                                              alldata = cbind(y ,NA)),
                                  effects = list(xi.filed = 1: spde$n.spde),
                                  A = list(A_x),
                                  tag = "est.x")
    stack_response2 <- inla.stack(data = list(y = cbind(NA, y),
                                              alldata = cbind(NA, z)),
                                  effects = list(
                                    list()
                                  ),
                                  A = list(),
                                  tag = "est.y",
    )
    stack <- inla.stack(stk_x, stk_y)
    # The copy feature has to be defined when specifying the formula
  }
  
  join_stack <- inla.stack(stack_est, stack_val)
  
  st_formula <- y ~ -1 + Intercept + covaraites + f(spatial_field, 
                                                    model = st_spde,
                                                    group = spatial_field.group,
                                                    control.group = list(model = "ar1"))
  precprior <- list(theta = list(param=c(1, 0.1)))
}
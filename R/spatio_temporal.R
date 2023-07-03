library(INLA)


# spatial_temporal <- function(formula = NULL,
#                              data = NULL,
#                              coordicates = NULL,
#                              boundary = NULL,
#                              family = NULL) {

spatial_temporal <- function(st_formula_joint, 
                             family = c('gaussian','binomial','gaussian','gaussian'),
                             data = inla.stack.data(stack_joint2),
                             Ntrials=inla.stack.data(stack_joint2)$Ntrials,
                             control.predictor = list(A = inla.stack.A(stack_joint2), compute = F),
                             verbose = TRUE) {
  # Args:
  #   family: a string that specifies the distribution of the data.
  #           Check \code{names(inla.models()$likelihood)} for the full list acceptable families.
  
  # Prepare standardized variables/factors and for INLA formula ------------------------------------

  st_input <- standardize_input(data, coordinates, boundary)
  loc <- st_input$loc
  boudary <- st_input$boundary
  covariates <- st_input$covariates
  obvervations <- st_input$observations
  
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
  
  # Model fitting ----------------------------------------------------------------------------------
  
  inla_obj <- inla(formula = st_formula_joint,
                   family = family,
                   data = inla.stack.data(stack_joint,
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


# Make spatial prediction ==========================================================================

spatial_predict <- function(st_mesh, 
                            ) {
  # predict the marginal distribution of linear predictor, which is 
  # computationally expensive.
  A_pred <- inla.spde.make.A(mesh = st_mesh,
                             loc = locations_to_be_predicted,
                             group = i_day, # select the time point for prediction
                             n.group = n_days)
  
  stack_pred_A <- list(A_pred, 1)
  stack_pred_effects <- list(c(grid_index, list(Intercept = 1)),
                             list(covariate_matrix))
  stack_pred <- inla.stack(data = list(y = NA),
                           A = stack_pred_A,
                           effects = stack_pred_effects,
                           tag = "prediction")
  stack_joint <- inla.stack(stack_est, stack_pred)
  pred_output <- inla(formula = st_formula,
                      family = "",
                      data = inla.stack.data(stack_join, spde=st_spde),
                      control.predictor = list(A = inla.stack.A(stack_join),
                                               compute = TRUE),
                      control.mode = list(theta = fitted$mode$theta,
                                          restart = FALSE),
                      quantiles = NULL, # save space
                      control.results = list(return.marginals.random = FALSE, # Save space
                                             return.marginals.predictor = FALSE)
                      )
  # To obtain the prediction of the data for the selected day, 
  # We just need to extract the posterior marginals of the linear predictor
  index_pred <- inla.stack.index(stack, "prediction")$data
  lp_marginals <- output$marginals$linear.predictor[index_pred] # or $fitted.values[index_pred]
  lp_mean <- unlist(lapply(lp_marginals,
                           function(x) inla.emarginal(exp, x)))
  lp_grid_mean <- matrix(lp_mean, n_x, n_y, byrow = T)
}


spatial_predict <- function() {
  # If the linear predictor marginal distributions are not necessary,
  # a less computationally heavy solution consists in projecting the latent field 
  # – estimated at the mesh vertices – onto the grid locations.
  A_pred <- inla.spde.make.A(mesh=Swiss.mesh) # Or this is equivalent to A = 1, so that its just on mesh point? page101
  stack_pred_A <- list() 
  stack_pred_effects <- list()
  stack_pred <- inla.stack(data = list(y = NA),
                           effects = stack_pred_effects,
                           tag = "prediction")
  stack_joint <- inla.stack(stack_est, stack_pred)
  pred_output <- inla(formula = st_formula,
                      family = "",
                      data = inla.stack.data(),
                      control.predictor = list(A = inla.stack.A(stack_joint),
                                               compute = TRUE),
                      control.compute = list(config = TRUE),
                      quantiles = NULL,
                      control.results = list(return.marginals.random = FALSE,
                                             return.marginals.predictor = FALSE)
                      )
  index_pred <- inla.stack.index(stack_joint, "prediction")$data
  post_mean_pred <- pred_out$summary.linear.predictor[index.pred, "mean"]
  post_sd_pred <- pred_out$summary.linear.predictor[index.pred, "sd"]
  
  proj_grid <- inla.mesh.projector(st_mesh,
                                   points)
  
  mean_pred <- inla.mesh.project(proj_grid, post_mean_pred)
  sd_pred <- inla.mesh.project(proj_grid, post_sd_pred)
  
  # Taking into account the uncertainty by sampling from its posterior distributions.
  # This can be done useing inla.posterior.sample function which needs the output from the inla
  # function by setting config = TRUE
  
  
  s1 <- inla.posterior.sample(n = 1, result = pred_output)
  
  # Since all the latent field components are stacked, we need to query the index for each component.
  
  # Getting samples from the posterior distribution of the linear predictor at mesh nodes.
  # Then these values are interpolated to the grid and the expected value computed in the response 
  # scale by applying the inverse of the link function.
  pred.nodes <- exp(sapply(s1, function(x) x$latent[index_pred]))
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

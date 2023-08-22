
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

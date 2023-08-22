
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


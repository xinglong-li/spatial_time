library(INLA)


spatial_temporal <- function(formula = NULL,
                             data = NULL,
                             boundary = NULL,
                             family = NULL) {
  # Args:
  #   family: a string that specifies the distribution of the data.
  #           Check \code{names(inla.models()$likelihood)} for the full list acceptable families.
  
  st_mesh <- inla.mesh.2d(loc = ,
                          loc.domain = ,
                          max.edge = ,
                          offset = ,
                          cutoff = )
  st_spde <- inla.spde2.matern(mesh = st_mesh)
  
  A_est <- inla.spde.make.A(mesh = st_mesh,
                             loc = )
  grid_index <- inla.spde.make.index(name = "spatial_field",
                                   n.spde = st_spde$n.spde)
  
  stack_est_A <- list()
  stack_est_effects <- list()
  stack_est <- inla.stack(data = NULL,
                          A = stack_est_A,
                          effects = stack_est_effects,
                          tag = "estimation")
  
  join_stack <- inla.stack(stack_est, stack_val)
  
  st_formula <- NULL
  
  inla_obj <- inla(formula = ,
                   family = "",
                   data = ,
                   control.predictor = list(),
                   control.compute = list())
  st_obj <- list(inla_obj = inla_obj)
  class(st_obj) <- "spatial_temporal"

  return(st_obj)
}


spatial_predict <- function(st_mesh, 
                            ) {
  # predict the marginal distribution of linear predictor, which is 
  # computationally expensive.
  A_pred <- inla.spde.make.A(mesh = st_mesh,
                             )
  
  stack_pred_A <- list()
  stack_pred_effects <- list()
  stack_pred <- inla.stack(data = list(y = NA),
                           A = stack_pred_A,
                           effects = stack_pred_effects,
                           tag = "prediction")
  stack_joint <- inla.stack(stack_est, stack_pred)
  pred_output <- inla(formula = st_formula,
                      family = "",
                      data = inla.stack.data(stack_join, spde=st_spde),
                      control.predictor = list(A = inla.stack.A(stack_join),
                                               compute = TRUE)
                      )
}


spatial_predict <- function() {
  # If the linear predictor marginal distributions are not necessary, 
  # a less computationally heavy solution consists in projecting the latent field 
  # – estimated at the mesh vertices – onto the grid locations.
  A_pred <- inla.spde.make.A(mesh=Swiss.mesh)
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
                                               compute = TRUE)
                      )
  index_pred <- inla.stack.index(stack_joint, "prediction")$data
  post_mean_pred <- pred_out$summary.linear.predictor[index.pred, "mean"]
  post_sd_pred <- pred_out$summary.linear.predictor[index.pred, "sd"]
  
  proj_grid <- inla.mesh.projector(st_mesh,
                                   points)
  
  mean_pred <- inla.mesh.project(proj_grid, post_mean_pred)
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



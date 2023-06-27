
standardize_input <- function (data, boundary) {
  loc <- NULL
  boundary <- NULL
  
  if (is.null(boundary)) {
    boundary <- inla.nonconvex.hull(points = coords)
  }

  return(list(loc = loc,
              boundary = boundary,
              observations = observations,
              covarites = covariates)
         )
}
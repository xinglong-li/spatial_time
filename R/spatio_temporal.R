library(INLA)


spatial_temporal <- function(data = NULL, 
                             ) {
  model_obj <- INLA::inla()
  model_obj <- list()
  class(model_obj) <- "spatial_temporal"

  return(model_obj)
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



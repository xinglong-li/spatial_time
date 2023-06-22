library(INLA)


spatial_temporal <- function(formula = NULL,
                             data = NULL,
                             boundary = NULL,
                             family = NULL) {
  # Args:
  #   family: a string that specifies the distribution of the data.
  #           Check \code{names(inla.models()$likelihood)} for the full list acceptable families.
  
  inla_obj <- INLA::inla()
  st_obj <- list(inla_obj = inla_obj)
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



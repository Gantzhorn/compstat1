#' Estimate density in RCPP with plugin oracle bandwidth selection
#'
#' @param .data The data for estimation
#' @param .x Points to evaluate the density estimate
#' @param .kernel The kernel for density estimation
#' @param .npoints If .x is not entered then the range of .data will be used.
#' @param .bw The bandwidth estimator
#'
#' @return Density estimate object
#' @export
#'
estimate_density <- function(.data, .kernel = "gaussian", .bw = plugin_oracle_bandwidth, .x = NULL, .npoints = 512L) {
  
  if(is.null(.x)) {
    start_point <- min(.data)
    end_point <- max(.data)
    
    .x <- seq(start_point, end_point, length.out = .npoints)
  }
  
  .kernel_number <- switch(
    .kernel,
    "gaussian" = 1L,
    "epanechnikov" = 2L,
    "uniform" = 3L
  )
  
  .bandwidth <- if(is.numeric(.bw)) {
    .bw
  } else {
    .bw(.data, .kernel)
    
  }
  
  .estimate <- calculate_density_cpp(data = .data,
                                    calculation_points = .x,
                                    lambda = .bandwidth,
                                    kernel_number = .kernel_number)
  
  structure(
    list(
      x = .x,
      y = .estimate,
      bandwidth = .bandwidth,
      kernel_name = .kernel
    ),
    class = "densityEstimate"
  )
  
}

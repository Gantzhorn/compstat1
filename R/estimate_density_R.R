source("./R/kernels_R.R")

calculate_density_R <- function(data,
                                    calculation_points,
                                    lambda,
                                    kernel_number) {
  
  kernel_call <- switch(kernel_number,
                        "1" = kernel_gaussian_R,
                        "2" = kernel_epanechnikov_R,
                        "3" = kernel_uniform_R)
  
  out_vector <- numeric(length(calculation_points))
  
  for(i in seq_along(calculation_points)) {
    out_vector[i] <- sum(kernel_call((data - calculation_points[i]) / lambda))
  }
  
  out_vector / (lambda * length(data))
}

estimate_density_R <- function(.data, .kernel = "gaussian", .bw = plugin_oracle_bandwidth, .x = NULL, .npoints = 512L) {
  
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
  
  .estimate <- calculate_density_R(data = .data,
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

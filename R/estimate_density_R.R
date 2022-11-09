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



calculate_density_R_profiler <- function(data,
                                calculation_points,
                                lambda,
                                kernel_number) {
  
  kernel_call <- switch(kernel_number,
                        "1" = kernel_gaussian_R,
                        "2" = kernel_epanechnikov_R,
                        "3" = kernel_uniform_R)
  
  out_vector <- numeric(length(calculation_points))
  
  for(i in seq_along(calculation_points)) {
    step1 <- (data - calculation_points[i])/lambda
    step2 <- kernel_call(step1)
    out_vector[i] <- sum(step2)
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

calculate_density_R_profiler(test, x1, 0.6, 2)

density_estimator_apply <- function(data, m = 512) {
  h <- plugin_oracle_bandwidth(data, .bandwidth_selector = "epanechnikov")
  rg <- range(data)
  data_seq <- seq(rg[1], rg[2], length.out = m)
  dens <- purrr::map_dbl(data_seq, function(x) {
    sum(kernel_epanechnikov_R((x - data) / h))
  }) 
  dens/(h*length(data))
}

density_estimator_apply2 <- function(data, m = 512) {
  h <- plugin_oracle_bandwidth(data, .bandwidth_selector = "epanechnikov")
  rg <- range(data)
  data_seq <- seq(rg[1], rg[2], length.out = m)
  dens <- sapply(data_seq, function(x) {
    sum(kernel_epanechnikov_R((x - data) / h))
  }) 
  dens/(h*length(data))
}

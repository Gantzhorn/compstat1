kernel_constant <- function(x) {
  UseMethod("kernel_constant")
}

kernel_constant.oracle_bandwidth <- function(x) {
  (8 * sqrt(pi) * x$square_norm / (3 * x$sigma))^0.2
}

gaussian_lambda <- function() {
  structure(
    list(
      name = "Gaussian",
      sigma = 1,
      square_norm = 1 / (2 * sqrt(pi))
    ),
    class = "oracle_bandwidth"
  )
}

epanechnikov_lambda <- function() {
  structure(
    list(
      name = "Epanechnikov",
      sigma = 1 / 25,
      square_norm = 3 / 5
    ),
    class = "oracle_bandwidth"
  )
}

uniform_lambda <- function() {
  structure(
    list(
      name = "Uniform",
      sigma = 1 / 9,
      square_norm = 1 / 2
    ),
    class = "oracle_bandwidth"
  )
}

plugin_oracle_bandwidth <- function(.data, .bandwidth_selector = "gaussian") {
  
  kernel_attr <- switch(
    .bandwidth_selector,
    "gaussian" = gaussian_lambda(),
    "epanechnikov" = epanechnikov_lambda(),
    "uniform" = uniform_lambda()
  )
  
  sigma_silverman <- min(sd(.data), IQR(.data) / 1.34)
  
  constant <- kernel_constant(kernel_attr)
  
  constant * sigma_silverman * length(.data)^(-0.2)
  
}

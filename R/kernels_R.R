kernel_epanechnikov_R <- function(x) {
  0.75 * (1 - x^2) * (abs(x) <= 1)
}

kernel_gaussian_R <- function(x) {
  0.3989422804014327 * exp(-0.5 * x^2)
}

kernel_uniform_R <- function(x) {
  0.5 * (abs(x) <= 1)
}

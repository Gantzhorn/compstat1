library(magrittr)
library(ggplot2)
library(microbenchmark)
library(Rcpp)

theme_set(theme_bw())

sourceCpp("./src/estimate_density.cpp")

source("./R/estimate_density_R.R")

test <- rnorm(1000L)
x1 <- seq(min(test), max(test), length.out = 512L)

microbenchmark(calculate_density_cpp(test, x1, 0.6, 2),
               calculate_density_R(test, x1, 0.6, 2))


source("./R/estimate_density_rcpp.R")
source("./R/plugin_oracle_bandwidth.R") # Bandwidth selection using oracle bandwidth

microbenchmark(plugin_oracle_bandwidth(test, "epanechnikov"))

microbenchmark(
  estimate_density(test, .kernel = "epanechnikov"),
  estimate_density_R(test, .kernel = "epanechnikov"),
  density(test, kernel = "epanechnikov", bw = "nrd")
  )

test2 <- rnorm(8192L)

rcpp_density_times <- 
  microbenchmark(
    estimate_density(test2[1:32], .kernel = "epanechnikov"),
    estimate_density(test2[1:64], .kernel = "epanechnikov"),
    estimate_density(test2[1:128], .kernel = "epanechnikov"),
    estimate_density(test2[1:256], .kernel = "epanechnikov"),
    estimate_density(test2[1:512], .kernel = "epanechnikov"),
    estimate_density(test2[1:1024], .kernel = "epanechnikov"),
    estimate_density(test2[1:2048], .kernel = "epanechnikov"),
    estimate_density(test2[1:4096], .kernel = "epanechnikov"),
    estimate_density(test2[1:8192], .kernel = "epanechnikov"),
    unit = "milliseconds"
  ) %>% 
  summary() %>% 
  extract2("median")

r_density_times <- 
  microbenchmark(
    estimate_density_R(test2[1:32], .kernel = "epanechnikov"),
    estimate_density_R(test2[1:64], .kernel = "epanechnikov"),
    estimate_density_R(test2[1:128], .kernel = "epanechnikov"),
    estimate_density_R(test2[1:256], .kernel = "epanechnikov"),
    estimate_density_R(test2[1:512], .kernel = "epanechnikov"),
    estimate_density_R(test2[1:1024], .kernel = "epanechnikov"),
    estimate_density_R(test2[1:2048], .kernel = "epanechnikov"),
    estimate_density_R(test2[1:4096], .kernel = "epanechnikov"),
    estimate_density_R(test2[1:8192], .kernel = "epanechnikov"),
    unit = "milliseconds"
  ) %>% 
  summary() %>% 
  extract2("median")

default_density_times <- 
  microbenchmark(
    density(test2[1:32], kernel = "epanechnikov", bw = "nrd"),
    density(test2[1:64], kernel = "epanechnikov", bw = "nrd"),
    density(test2[1:128], kernel = "epanechnikov", bw = "nrd"),
    density(test2[1:256], kernel = "epanechnikov", bw = "nrd"),
    density(test2[1:512], kernel = "epanechnikov", bw = "nrd"),
    density(test2[1:1024], kernel = "epanechnikov", bw = "nrd"),
    density(test2[1:2048], kernel = "epanechnikov", bw = "nrd"),
    density(test2[1:4096], kernel = "epanechnikov", bw = "nrd"),
    density(test2[1:8192], kernel = "epanechnikov", bw = "nrd"),
    unit = "milliseconds"
  ) %>% 
  summary() %>% 
  extract2("median")

plot_dataframe <- data.frame(times = c(rcpp_density_times, r_density_times, default_density_times),
                             n = rep(c(32, 64, 128, 256, 512, 1024, 2048, 4096, 8192), 3),
                             implementation = rep(c("RCPP", "R", "Default"), each = 9))

plot_dataframe %>% 
  ggplot(aes(x = n, y = times, col = implementation)) +
  geom_line(size = 1) + 
  scale_x_log10() +
  scale_y_log10() +
  xlab("Number of observations") +
  ylab("Runtime [milliseconds]")


# Plots of density estimates for simulated data
density_plot <- estimate_density(test, "epanechnikov")
default_density <- density(test, kernel = "epanechnikov", bw = "nrd")

plot(density_plot)
lines(default_density, col = "red", lwd = 2)

# Plots of density estimates for real data
density_plot_real <- estimate_density(infert$age, "epanechnikov")
default_density_real <- density(infert$age, kernel = "epanechnikov", bw = "nrd")

plot(density_plot_real)
lines(default_density_real, col = "red", lwd = 2)


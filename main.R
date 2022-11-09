library(tidyverse)
library(magrittr)
library(profvis)
library(microbenchmark)
library(Rcpp)
library(xtable)

ggplot2::theme_set(theme_bw())

sourceCpp("./src/estimate_density.cpp")
test <- rnorm(100000L)
x1 <- seq(min(test), max(test), length.out = 512L)

source("./R/estimate_density_R.R")
standard_est <- tibble(x = density(test, kernel = "epanechnikov", bw = 0.6, from = min(x1), to = max(x1))$x,
                       y = density(test, bw = 0.6, kernel = "epanechnikov", from = min(x1), to = max(x1))$y,
                       type = "density")
R_est <- calculate_density_R(test, x1, 0.6, 2) %>%
  tibble(y =., x = x1, type = "R_implementation")

bind_rows(standard_est, R_est) %>% ggplot(aes(x = x, y = y, col = type)) +
  geom_line(size = 1) + ylab("Density estimate") +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size =12),
        strip.text = element_text(face = "bold", size = 12))

bind_rows(standard_est, R_est) %>% pivot_wider(names_from = type, 
                                               values_from = y) %>% 
  mutate(Forskel = (density-R_implementation)^2) %>% 
  ggplot(aes(x= x, y = Forskel)) + 
  geom_line(size = 1.2) + ylab("Squared error") + 
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size =12),
        strip.text = element_text(face = "bold", size = 12))



profvis::profvis(source("./R/estimate_density_R.R"))


initialRvRcpp <- microbenchmark(calculate_density_cpp(test, x1, 0.6, 2),
               calculate_density_R(test, x1, 0.6, 2), times = 50)

initialRvRcpp %>% as_tibble() %>%
  mutate(Implementation = ifelse(str_detect(expr, "cpp"), "Rcpp", "R"),
                                         time = time/1000000) %>% 
  ggplot(aes(x = log(time), y = Implementation)) + geom_violin() +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size =12),
        strip.text = element_text(face = "bold", size = 12))

source("./R/estimate_density_rcpp.R")
source("./R/plugin_oracle_bandwidth.R") # Bandwidth selection using oracle bandwidth

microbenchmark(plugin_oracle_bandwidth(test, "epanechnikov"))

microbenchmark(
  estimate_density(test, .kernel = "epanechnikov"),
  estimate_density_R(test, .kernel = "epanechnikov"),
  density(test, kernel = "epanechnikov", bw = "nrd")
  )

test2 <- rnorm(32768)

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
    estimate_density(test2[1:16384], .kernel = "epanechnikov"),
    estimate_density(test2[1:32768], .kernel = "epanechnikov"),
    unit = "milliseconds"
  ) %>% 
  summary() %>% 
  magrittr::extract2("median")

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
    estimate_density_R(test2[1:16384], .kernel = "epanechnikov"),
    estimate_density_R(test2[1:32768], .kernel = "epanechnikov"),
    unit = "milliseconds"
  ) %>% 
  summary() %>% 
  magrittr::extract2("median")

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
    density(test2[1:16384], kernel = "epanechnikov", bw = "nrd"),
    density(test2[1:32768], kernel = "epanechnikov", bw = "nrd"),
    unit = "milliseconds"
  ) %>% 
  summary() %>% 
  magrittr::extract2("median")

r_apply_times <- 
  microbenchmark(
    density_estimator_apply(test2[1:32]),
    density_estimator_apply(test2[1:64]),
    density_estimator_apply(test2[1:128]),
    density_estimator_apply(test2[1:256]),
    density_estimator_apply(test2[1:512]),
    density_estimator_apply(test2[1:1024]),
    density_estimator_apply(test2[1:2048]),
    density_estimator_apply(test2[1:4096]),
    density_estimator_apply(test2[1:8192]),
    density_estimator_apply(test2[1:16384]),
    density_estimator_apply(test2[1:32768]),
    unit = "milliseconds"
  ) %>% 
  summary() %>% 
  magrittr::extract2("median")

r_apply_times2 <- 
  microbenchmark(
    density_estimator_apply2(test2[1:32]),
    density_estimator_apply2(test2[1:64]),
    density_estimator_apply2(test2[1:128]),
    density_estimator_apply2(test2[1:256]),
    density_estimator_apply2(test2[1:512]),
    density_estimator_apply2(test2[1:1024]),
    density_estimator_apply2(test2[1:2048]),
    density_estimator_apply2(test2[1:4096]),
    density_estimator_apply2(test2[1:8192]),
    density_estimator_apply2(test2[1:16384]),
    density_estimator_apply2(test2[1:32768]),
    unit = "milliseconds"
  ) %>% 
  summary() %>% 
  magrittr::extract2("median")

plot_dataframe <- data.frame(times = c(rcpp_density_times, r_density_times, default_density_times, r_apply_times, r_apply_times2),
                             n = rep(c(32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768), 5),
                             implementation = rep(c("RCPP", "R", "Default", "Purrr::map", "apply"), each = 11))

plot_dataframe %>% 
  ggplot(aes(x = n, y = times, col = implementation)) +
  geom_line(size = 1.2) + geom_point() + 
  scale_x_log10() +
  scale_y_log10() +
  xlab("Number of observations") +
  ylab("Runtime [milliseconds]") +
  scale_color_brewer(palette="Set1") +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size =12),
        strip.text = element_text(face = "bold", size = 12)) +
  scale_fill_discrete(labels=c('sapply_R','Density','purr::map','forloop_R','Rcpp'))


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

# Testing on worse real data
bestand_data <- read_csv("bestanddata.csv")

UDGIFTER <- bestand_data %>% filter(UDGIFT>0)

# bestand_data %>% select(VAERDI, AAR, GRUNDAREAL, EKSTRAAREAL)
realDataEstimate <- estimate_density(bestand_data$VAERDI, "epanechnikov")
realDataEstimate1 <- density(bestand_data$VAERDI, kernel = "epanechnikov", from = min(bestand_data$VAERDI), to = max(bestand_data$VAERDI))

bind_rows(tibble(x = realDataEstimate$x, y = realDataEstimate$y, Implementation = "Rcpp"),
          tibble(x = realDataEstimate1$x, y = realDataEstimate1$y, Implementation = "Density" , bw = "nrd")) %>%
  filter(x<5.0e+06) %>% 
  ggplot(aes(x = x, y = y, col = Implementation)) + geom_line(size = 1) +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size =12),
        strip.text = element_text(face = "bold", size = 12))


bind_rows(tibble(x = realDataEstimate$x, y = realDataEstimate$y, Implementation = "Rcpp"),
          tibble(x = realDataEstimate1$x, y = realDataEstimate1$y, Implementation = "Density")) %>% 
  pivot_wider(values_from = y, names_from = Implementation) %>% 
  filter(x<5.0e+06) %>% 
  mutate(Difference = abs(Rcpp-Density)) %>% ggplot(aes(x = x, y = Difference)) +
  geom_line(size = 1) + theme(axis.title = element_text(face = "bold", size = 14),
                              axis.text = element_text(face = "bold", size = 12),
                              legend.title = element_text(face = "bold", size = 16),
                              legend.text = element_text(face = "bold", size =12),
                              strip.text = element_text(face = "bold", size = 12))

#Fit to bimodial
normaldata1 <- rnorm(10000, mean = 5, sd = 1)
normaldata2 <- rnorm(10000, mean = 10, sd = 1)

bernoulli1 <- rbinom(10000, 1, prob = 0.5)

bimodialdata <- normaldata1*(bernoulli1)+normaldata2*(1-bernoulli1)
  
bimodialEstRcpp <- estimate_density(bimodialdata, "epanechnikov")
bimodialEstDens <- density(bimodialdata, kernel = "epanechnikov", from = min(bimodialdata), to = max(bimodialdata), bw = "nrd")

bind_rows(tibble(x = bimodialEstRcpp$x, y = bimodialEstRcpp$y, Implementation = "Rcpp"),
          tibble(x = bimodialEstDens$x, y = bimodialEstDens$y, Implementation = "Density")) %>%
  ggplot(aes(x = x, y = y, col = Implementation)) + geom_line(size = 1) +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size =12),
        strip.text = element_text(face = "bold", size = 12))


bind_rows(tibble(x = bimodialEstRcpp$x, y = bimodialEstRcpp$y, Implementation = "Rcpp"),
          tibble(x = bimodialEstDens$x, y = bimodialEstDens$y, Implementation = "Density")) %>% 
  pivot_wider(values_from = y, names_from = Implementation) %>% 
  filter(x<5.0e+06) %>% 
  mutate(Difference = Rcpp-Density) %>% ggplot(aes(x = x, y = Difference)) +
  geom_line(size = 1) + theme(axis.title = element_text(face = "bold", size = 14),
                              axis.text = element_text(face = "bold", size = 12),
                              legend.title = element_text(face = "bold", size = 16),
                              legend.text = element_text(face = "bold", size =12),
                              strip.text = element_text(face = "bold", size = 12))

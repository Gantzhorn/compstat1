#include <Rcpp.h>
#include <cmath>

double kernel_epanechnikov(double x) {
  if(std::abs(x) > 1) {
    return 0;
  } else {
    return 0.75 * (1 - x * x);
  }
}

double kernel_gaussian(double x) {
  return 0.3989422804014327 * std::exp(-0.5 * x * x);
}

double kernel_uniform(double x) {
  if(std::abs(x) > 1) {
    return 0;
  } else {
    return 0.5;
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector calculate_density_cpp(Rcpp::NumericVector data,
                                         Rcpp::NumericVector calculation_points,
                                         double lambda,
                                         int kernel_number) {
  
  
  Rcpp::NumericVector result_vector(calculation_points.length());
  
  double (*kernel)(double);
  
  switch(kernel_number) {
  case 1:
    kernel = &kernel_gaussian;
    break;
    
  case 2:
    kernel = &kernel_epanechnikov;
    break;
    
  case 3:
    kernel = &kernel_uniform;
    break;
    
  default:
    kernel = &kernel_gaussian;
  break;
  }
  
  int points = calculation_points.length();
  int data_length = data.length();
  
  for(int i = 0; i < points; ++i) {
    
    result_vector[i] = 0;
    
    for(int j = 0; j < data_length; ++j) {
      result_vector[i] += kernel((data[j] - calculation_points[i]) / lambda);
    }
    
  }
  
  return result_vector / (lambda * data_length);
}

// This header configures R interface wrappers ("wr_*.cpp" files)
// It is only used at the R package level and should only be included in such.

#ifndef ACTIONET_R_CONFIG_H
#define ACTIONET_R_CONFIG_H

// Automatically return all arma::vec and related as Rcpp:NumericVector instead of 1-D matrix
#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#ifndef LIBACTIONET_BUILD_R
    #define LIBACTIONET_BUILD_R
#endif //LIBACTIONET_BUILD_R

// [[Rcpp::interfaces(r)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // For some dumbass reason R won't use Rcpp from the header below.
#include "libactionet.hpp"

#endif //ACTIONET_R_CONFIG_H

// This header configures R interface wrappers ("wr_*.cpp" files)
// It is only used at the R package level and should only be included in such.

#ifndef ACTIONET_R_CONFIG_H
#define ACTIONET_R_CONFIG_H

// Automatically return all arma::vec and related as Rcpp:NumericVector instead of 1-D matrix
#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#ifndef LIBACTIONET_BUILD_R
    #define LIBACTIONET_BUILD_R
#endif //LIBACTIONET_BUILD_R

#include <RcppArmadillo.h> // For some dumbass reason R won't use Rcpp from the header below.
#include "libactionet_config.hpp"
#include "libactionet.hpp"
// [[Rcpp::depends(RcppArmadillo)]]

//' Set the RNG Seed from within Rcpp
//'
//' Within Rcpp, one can set the R session seed without triggering
//' the CRAN rng modifier check.
//' @param seed A \code{double} that is the seed one wishes to use.
//' @return A set RNG scope.
//' @examples
//' set.seed(10)
//' x = rnorm(5,0,1)
//' set_seed(10)
//' y = rnorm(5,0,1)
//' all.equal(x,y, check.attributes = F)
// [[Rcpp::export]]
inline void set_seed(double seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(std::floor(std::fabs(seed)));
}

// TODO: Necessary?
// template <typename T>
// Rcpp::NumericVector arma2vec(const T& x) {
//     return Rcpp::NumericVector(x.begin(), x.end());
// }

#endif //ACTIONET_R_CONFIG_H

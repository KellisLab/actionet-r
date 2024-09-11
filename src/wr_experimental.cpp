// Rcpp interface for `action` module
// Organized by module header in th order imported.
// [[Rcpp::interfaces(r, cpp)]]
#include <visualization/color_map.hpp>
#include "actionet_r_config.h"

// aa ==================================================================================================================

// [[Rcpp::export]]
arma::mat computeNodeColors(arma::mat& coordinates, int thread_no) {
    arma::mat rgb = actionet::computeNodeColors(coordinates, thread_no);
    return rgb;
}




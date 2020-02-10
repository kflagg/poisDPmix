#include <RcppArmadillo.h>
#include "distributions.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec density_mcmcMap(const arma::field<arma::mat>& mu,
                          const arma::field<arma::cube>& Sigma,
                          const arma::field<arma::ivec>& h,
                          const arma::vec& location){

  // Get the number of MC draws.
  int n_iter = mu.n_elem;

  // Initialize a vector to accumulate the density for each draw
  // across components.
  arma::vec dense(n_iter, arma::fill::zeros);

  // Temporary storage number of components.
  int L;
  int sum_h;

  // Loop through draws and components.
  for(int j = 0; j < n_iter; j++){
    L = h(j).n_elem;
    sum_h = 0;

    // Loop through components.
    for(int l = 0; l < L; l++){

      // Add the density for this component.
      if(h(j)(l) > 0){
        dense(j) += h(j)(l) *
                      densmvn(location, mu(j).col(l), Sigma(j).slice(l));
        sum_h += h(j)(l);
      }

    }

    // Scale to integrate to 1.
    dense(j) /= sum_h;
  }

  // Return the density vector.
  return dense;
}

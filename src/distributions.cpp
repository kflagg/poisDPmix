#include <RcppArmadillo.h>

// Generate a single draw from a Wishart(nu, S) distribution.
arma::mat randw(const int& nu, const arma::mat& S){
  arma::mat alpha = arma::sqrtmat_sympd(S) * arma::randn(S.n_rows, nu);
  return alpha * arma::trans(alpha);
}

// Generate a single draw from a multivariate N(mu, S) distribution.
arma::mat randmvn(const arma::vec& mu, const arma::mat& S){
  return arma::sqrtmat_sympd(S) * arma::randn(mu.n_elem) + mu;
}

// Multivariate N(mu, S) density (from Hoff).
double densmvn(const arma::vec& x, const arma::vec& mu, const arma::mat& S){
  int d = mu.n_elem;
  arma::vec dev = x - mu;
  arma::mat expon = arma::trans(dev) * arma::inv_sympd(S) * dev;
  double ldet;
  double sgn;
  arma::log_det(ldet, sgn, S);
  return exp((d / -2.0) * log(2.0 * arma::datum::pi) -
           0.5 * ldet - expon(0, 0) / 2.0);
}

// Location-shifted multivariate t(nu, mu, S) density (from Genz and Bretz).
double densmvt(const arma::vec& x, const double& nu, const arma::vec& mu, const arma::mat& S){
  int d = x.n_elem;
  double nud2 = (nu + d) * 0.5;
  arma::vec dev = x - mu;
  arma::mat ss = arma::trans(dev) * arma::inv_sympd(S) * dev;
  double ldet;
  double sgn;
  arma::log_det(ldet, sgn, S);
  return exp(lgamma(nud2) - lgamma(nu * 0.5) - 0.5 * ldet -
         d * 0.5 * log(nu * arma::datum::pi) - nud2 * log(1 + ss(0, 0) / nu));
}

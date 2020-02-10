#include <RcppArmadillo.h>

arma::mat randw(const int& nu, const arma::mat& S);
arma::mat randmvn(const arma::vec& mu, const arma::mat& S);
double densmvn(const arma::vec& x, const arma::vec& mu, const arma::mat& S);
double densmvt(const arma::vec& x, const double& nu, const arma::vec& mu, const arma::mat& S);

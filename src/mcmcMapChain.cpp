#include <RcppArmadillo.h>
#include <gdal/ogrsf_frmts.h>
#include <gdal/ogr_geometry.h>
#include <time.h>
#include "distributions.h"
#include "inRegion.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Fully observed, unbounded. {{{
// [[Rcpp::export]]
Rcpp::List mcmc_chain_unb(const std::string& chain_id, const arma::mat& x,
                          const double& a, const double& b,
                          const arma::vec& eta, const arma::mat& Psi,
                          const int& nu, const arma::mat& Omega,
                          const double& alpha, const int& n_iter,
                          const int& n_burnin, const int& n_thin,
                          arma::vec& Lambda, arma::imat& k, arma::ivec& m,
                          const arma::mat& mu_initial,
                          const arma::cube& Sigma_initial){

  // Intermediate variables used in computations.
  int n = x.n_cols;
  int d = x.n_rows;
  arma::mat Psi_inv = arma::inv_sympd(Psi);
  arma::vec Psi_inv_eta = Psi_inv * eta;
  arma::mat Omega_inv = arma::inv_sympd(Omega);
  double c_0 = 2.0 * nu - d + 1.0;
  arma::mat B_0 = 2.0 / c_0 * (Psi + Omega); // Check this.

  int h_used = arma::max(k.row(0)) + 1; // Number of components used.
  int i_new; // Label of observation to copy component label from.
  int k0; // Index for accessing k_cdf.
  arma::mat k_quant = arma::randu(n_iter, n);  // Quantiles for choosing k.
  arma::vec k_cdf(n); // Cumulative distribution of k.
  arma::vec q(n);  // Complete conditional mixing weights.
  arma::vec sum_x(d); // Sum of component members.
  arma::mat sum_sq(d, d); //Component sum of squares.
  arma::uvec unused; // Indices of unused component labels.
  arma::uvec members; // Indices of component members.
  arma::vec dev(d); // Devation of an observation from the component center.
  arma::mat mu_cov_inv(d, d); // Current precision of mu.
  arma::mat mu_cov(d, d); // Current covariance of mu.
  arma::cube Sigma_inv; // Current precisions.

  // Initialize storage. (Lambda, k, and m passed by reference.)
  arma::field<arma::mat> mu(n_iter);      // Component means.
  arma::field<arma::cube> Sigma(n_iter);  // Component covariances.
  arma::field<arma::ivec> h(n_iter);      // Component sizes.
  m.zeros();                              // Number of components.

  // In the initial iteration, initialize mu, Sigma, and h based on the
  // initial values provided in k.
  int j = 0;

  // Look like we're doing something.
  Rcpp::Function msg("message");
  msg(chain_id + std::string(": iteration ") + std::to_string(j + 1) +
    std::string(" of ") + std::to_string(n_iter));

  // Create a vector of component sizes.
  h(j).set_size(h_used);

  // Set mu and Sigma to their initial values.
  mu(j) = mu_initial;
  Sigma(j) = Sigma_initial;

  // Create a cube of component precisions.
  Sigma_inv.set_size(d, d, Sigma(j).n_slices);

  // Loop through components.
  for(int l = 0; l < h_used; l++){

    // Update component size.
    members = find(k.row(j) == l);
    h(j)(l) = members.n_elem;

    // Update component precision.
    Sigma_inv.slice(l) = arma::inv_sympd(Sigma(j).slice(l));
  }

  // Find unused labels and update the number of components used.
  unused = find(h(j) < 1);
  m(j) = h_used - unused.n_elem;

  // Loop through iterations.
  for(j = 1; j < n_iter; j++){

    // Look like we're doing something.
    Rcpp::Function msg("message");
    msg(chain_id + std::string(": iteration ") + std::to_string(j + 1) +
      std::string(" of ") + std::to_string(n_iter));

    // Loop through observations.
    for(int i = 0; i < n; i++){

      // Compute mixing weights for conditional posteriors.
      // TODO: Find an efficient way to calculate h_l^(-i) and avoid making
      // n calls to densmvn().
      for(int idx = 0; idx < n; idx++){
        if(idx == i){
          q(idx) = alpha * densmvt(x.col(i), c_0, eta, B_0);
          // To understand the above, think "prior predictive density".
        }else{
          q(idx) = densmvn(x.col(i),
                           mu(j - 1).col(k(j - 1, idx)),
                           Sigma(j - 1).slice(k(j - 1, idx)));
        }
      }
      k_cdf = arma::cumsum(arma::normalise(q, 1));

      // Draw an observation's index and use its component label.
      i_new = 0;
      while(k_quant(j, i) > k_cdf(i_new)){
        i_new++;
      }

      // Update observation's component label.
      if(i_new == i){
       if(unused.n_elem > 0){
         k(j, i) = unused(0);
         unused.shed_row(0);
         Sigma_inv.slice(k(j, i)) = randw(nu, Omega_inv);
       } else {
         k(j, i) = h_used;
         h_used++;
         Sigma_inv = arma::join_slices(Sigma_inv, randw(nu, Omega_inv));
       }
      } else {
        k(j, i) = k(j - (i_new > i), i_new);
      }
    }

    // Create a vector of component sizes.
    h(j).set_size(h_used);

    // Create a matrix of component centers.
    mu(j).set_size(d, h_used);

    // Create a cube of component covariances.
    Sigma(j).set_size(d, d, h_used);

    // Loop through components.
    for(int l = 0; l < h_used; l++){

      // Update component size.
      members = find(k.row(j) == l);
      h(j)(l) = members.n_elem;

      if(h(j)(l) > 0){

        // Update component sum and covariance of center.
        sum_x = sum(x.cols(members), 1);
        mu_cov_inv = Psi_inv + h(j)(l) * Sigma_inv.slice(l);
        mu_cov = arma::inv_sympd(mu_cov_inv);

        // Draw a new center.
        mu(j).col(l) = randmvn(mu_cov * (Psi_inv_eta +
                                 Sigma_inv.slice(l) * sum_x),
                               mu_cov);

        // Update component sum of squares.
        sum_sq.zeros();
        for(int i = 0; i < h(j)(l); i++){
          dev = x.col(members(i)) - mu(j).col(l);
          sum_sq += dev * arma::trans(dev);
        }

        // Draw a new covariance.
        Sigma_inv.slice(l) = randw(nu + h(j)(l),
                                   arma::inv_sympd(Omega + sum_sq));
        Sigma(j).slice(l) = arma::inv_sympd(Sigma_inv.slice(l));

      } else {
        mu(j).col(l).fill(arma::datum::nan);
        Sigma(j).slice(l).fill(arma::datum::nan);
        Sigma_inv.slice(l).fill(arma::datum::nan);
      }
    }

    // Find unused labels and update the number of components used.
    unused = find(h(j) < 1);
    m(j) = h_used - unused.n_elem;

    // Draw from the posterior of Lambda
    Lambda(j) = (arma::randg(1,
                 arma::distr_param(a + n, 1.0 / (b + 1.0))))(0);
  }

  return Rcpp::List::create(Rcpp::Named("mu") = Rcpp::wrap(mu),
                            Rcpp::Named("Sigma") = Rcpp::wrap(Sigma),
                            Rcpp::Named("h") = Rcpp::wrap(h));
}
// }}}

// Fully observed, bounded. {{{
// [[Rcpp::export]]
Rcpp::List mcmc_chain(const std::string& chain_id, const arma::mat& x,
                      const double& a, const double& b,
                      const arma::vec& eta, const arma::mat& Psi,
                      const int& nu, const arma::mat& Omega,
                      const double& alpha, const int& n_iter,
                      const int& n_burnin, const int& n_thin,
                      arma::vec& Lambda, arma::vec& c,
                      arma::imat& k, arma::ivec& m,
                      const arma::mat& mu_initial,
                      const arma::cube& Sigma_initial,
                      const char* site_source, const char* site_layer,
                      const int& n_mcint){

  // Intermediate variables used in computations.
  int n = x.n_cols;
  int d = x.n_rows;
  arma::mat Psi_inv = arma::inv_sympd(Psi);
  arma::vec Psi_inv_eta = Psi_inv * eta;
  arma::mat Omega_inv = arma::inv_sympd(Omega);
  double c_0 = 2.0 * nu - d + 1.0;
  arma::mat B_0 = 2.0 / c_0 * (Psi + Omega); // Check this.

  int h_used = arma::max(k.row(0)) + 1; // Number of components used.
  int i_new; // Label of observation to copy component label from.
  int k0; // Index for accessing k_cdf.
  arma::mat k_quant = arma::randu(n_iter, n);  // Quantiles for choosing k.
  arma::vec k_cdf(n); // Cumulative distribution of k.
  arma::vec q(n);  // Complete conditional mixing weights.
  arma::vec sum_x(d); // Sum of component members.
  arma::mat sum_sq(d, d); //Component sum of squares.
  arma::uvec unused; // Indices of unused component labels.
  arma::uvec members; // Indices of component members.
  arma::vec dev(d); // Devation of an observation from the component center.
  arma::mat mu_cov_inv(d, d); // Current precision of mu.
  arma::mat mu_cov(d, d); // Current covariance of mu.
  arma::cube Sigma_inv; // Current precisions.

  int l_new; // Component label of new point.
  arma::mat x_mc(d, n_mcint); // Points for MC integration.
  arma::vec h_predict; // Component sizes for the predictive distribution.
  arma::vec k_mc_quant; // Label quantile for the MC integration.
  arma::vec in_site(n_mcint); // In-site indicator vector.
  OGRPoint test_x; // Point to test for intersections.

  // Initialize storage. (Lambda, k, and m passed by reference.)
  arma::field<arma::mat> mu(n_iter);      // Component means.
  arma::field<arma::cube> Sigma(n_iter);  // Component covariances.
  arma::field<arma::ivec> h(n_iter);      // Component sizes.
  m.zeros();                              // Number of components.

  // Load site shapefile.
  GDALAllRegister();
  GDALDataset *site_src;
  site_src = (GDALDataset *) GDALOpenEx(site_source, GDAL_OF_VECTOR,
                                        NULL, NULL, NULL);
  OGRLayer *site_shp;
  site_shp = site_src->GetLayerByName(site_layer);
  site_shp->ResetReading();
  OGRFeature *site_feature;
  site_feature = site_shp->GetNextFeature(); // Only using first feature...
  OGRGeometry *site_geom;
  site_geom = site_feature->GetGeometryRef();

  // In the initial iteration, initialize mu, Sigma, and h based on the
  // initial values provided in k.
  int j = 0;

  // Look like we're doing something.
  Rcpp::Function msg("message");
  msg(chain_id + std::string(": iteration ") + std::to_string(j + 1) +
    std::string(" of ") + std::to_string(n_iter));

  // Create a vector of component sizes.
  h(j).set_size(h_used);

  // Set mu and Sigma to their initial values.
  mu(j) = mu_initial;
  Sigma(j) = Sigma_initial;

  // Create a cube of component precisions.
  Sigma_inv.set_size(d, d, Sigma(j).n_slices);

  // Loop through components.
  for(int l = 0; l < h_used; l++){

    // Update component size.
    members = find(k.row(j) == l);
    h(j)(l) = members.n_elem;

    // Update component precision.
    Sigma_inv.slice(l) = arma::inv_sympd(Sigma(j).slice(l));
  }

  // Find unused labels and update the number of components used.
  unused = find(h(j) < 1);
  m(j) = h_used - unused.n_elem;

  // MC integration.
  h_predict = arma::vec(h_used + 1);
  h_predict.head(h_used) = arma::conv_to<arma::vec>::from(h(j));
  h_predict(h_used) = alpha;
  k_cdf = arma::cumsum(arma::normalise(h_predict, 1));
  k_mc_quant = arma::randu(n_mcint);  // Quantiles for choosing k.
  for(int idx = 0; idx < n_mcint; idx++){
    l_new = 0;
    while(k_mc_quant(idx) > k_cdf(l_new)){
      l_new++;
    }
    if(l_new < h_used){
      x_mc.col(idx) = randmvn(mu(j).col(l_new), Sigma(j).slice(l_new));
    } else {
      x_mc.col(idx) = randmvn(
        randmvn(eta, Psi),
        arma::inv_sympd(randw(nu, Omega_inv))
      );
    }
    test_x = OGRPoint(x_mc(0, idx), x_mc(1, idx));
    in_site(idx) = site_geom->Intersects((OGRGeometry *) &test_x);
  }

  c(j) = arma::mean(in_site);

  // Loop through iterations.
  for(j = 1; j < n_iter; j++){

    // Look like we're doing something.
    Rcpp::Function msg("message");
    msg(chain_id + std::string(": iteration ") + std::to_string(j + 1) +
      std::string(" of ") + std::to_string(n_iter));

    // Loop through observations.
    for(int i = 0; i < n; i++){

      // Compute mixing weights for conditional posteriors.
      // TODO: Find an efficient way to calculate h_l^(-i) and avoid making
      // n calls to densmvn().
      for(int idx = 0; idx < n; idx++){
        if(idx == i){
          q(idx) = alpha * densmvt(x.col(i), c_0, eta, B_0);
          // To understand the above, think "prior predictive density".
        }else{
          q(idx) = densmvn(x.col(i),
                           mu(j - 1).col(k(j - 1, idx)),
                           Sigma(j - 1).slice(k(j - 1, idx)));
        }
      }
      k_cdf = arma::cumsum(arma::normalise(q, 1));

      // Draw an observation's index and use its component label.
      i_new = 0;
      while(k_quant(j, i) > k_cdf(i_new)){
        i_new++;
      }

      // Update observation's component label.
      if(i_new == i){
       if(unused.n_elem > 0){
         k(j, i) = unused(0);
         unused.shed_row(0);
         Sigma_inv.slice(k(j, i)) = randw(nu, Omega_inv);
       } else {
         k(j, i) = h_used;
         h_used++;
         Sigma_inv = arma::join_slices(Sigma_inv, randw(nu, Omega_inv));
       }
      } else {
        k(j, i) = k(j - (i_new > i), i_new);
      }
    }

    // Create a vector of component sizes.
    h(j).set_size(h_used);

    // Create a matrix of component centers.
    mu(j).set_size(d, h_used);

    // Create a cube of component covariances.
    Sigma(j).set_size(d, d, h_used);

    // Loop through components.
    for(int l = 0; l < h_used; l++){

      // Update component size.
      members = find(k.row(j) == l);
      h(j)(l) = members.n_elem;

      if(h(j)(l) > 0){

        // Update component sum and covariance of center.
        sum_x = sum(x.cols(members), 1);
        mu_cov_inv = Psi_inv + h(j)(l) * Sigma_inv.slice(l);
        mu_cov = arma::inv_sympd(mu_cov_inv);

        // Draw a new center.
        mu(j).col(l) = randmvn(mu_cov * (Psi_inv_eta +
                                 Sigma_inv.slice(l) * sum_x),
                               mu_cov);

        // Update component sum of squares.
        sum_sq.zeros();
        for(int i = 0; i < h(j)(l); i++){
          dev = x.col(members(i)) - mu(j).col(l);
          sum_sq += dev * arma::trans(dev);
        }

        // Draw a new covariance.
        Sigma_inv.slice(l) = randw(nu + h(j)(l),
                                   arma::inv_sympd(Omega + sum_sq));
        Sigma(j).slice(l) = arma::inv_sympd(Sigma_inv.slice(l));

      } else {
        mu(j).col(l).fill(arma::datum::nan);
        Sigma(j).slice(l).fill(arma::datum::nan);
        Sigma_inv.slice(l).fill(arma::datum::nan);
      }
    }

    // Find unused labels and update the number of components used.
    unused = find(h(j) < 1);
    m(j) = h_used - unused.n_elem;

    // MC integration.
    h_predict = arma::vec(h_used + 1);
    h_predict.head(h_used) = arma::conv_to<arma::vec>::from(h(j));
    h_predict(h_used) = alpha;
    k_cdf = arma::cumsum(arma::normalise(h_predict, 1));
    k_mc_quant = arma::randu(n_mcint);  // Quantiles for choosing k.
    for(int idx = 0; idx < n_mcint; idx++){
      l_new = 0;
      while(k_mc_quant(idx) > k_cdf(l_new)){
        l_new++;
      }
      if(l_new < h_used){
        x_mc.col(idx) = randmvn(mu(j).col(l_new), Sigma(j).slice(l_new));
      } else {
        x_mc.col(idx) = randmvn(
          randmvn(eta, Psi),
          arma::inv_sympd(randw(nu, Omega_inv))
        );
      }
      test_x = OGRPoint(x_mc(0, idx), x_mc(1, idx));
      in_site(idx) = site_geom->Intersects((OGRGeometry *) &test_x);
    }

    // Draw from the posterior of Lambda
    c(j) = arma::mean(in_site);
    Lambda(j) = (arma::randg(1,
                 arma::distr_param(a + n, 1.0 / (b + c(j)))))(0);
  }

  // Close shapefiles.
  // TODO: Figure out how to do this earlier without destroying the geometry.
  // Then move it to another function.
  OGRFeature::DestroyFeature(site_feature);
  GDALClose(site_src);

  return Rcpp::List::create(Rcpp::Named("mu") = Rcpp::wrap(mu),
                            Rcpp::Named("Sigma") = Rcpp::wrap(Sigma),
                            Rcpp::Named("h") = Rcpp::wrap(h));
}
// }}}

// Partially observed, bounded. {{{
// [[Rcpp::export]]
Rcpp::List mcmc_chain_s(const std::string& chain_id, const arma::mat& x,
                        const double& a, const double& b,
                        const arma::vec& eta, const arma::mat& Psi,
                        const int& nu, const arma::mat& Omega,
                        const double& alpha, const int& n_iter,
                        const int& n_burnin, const int& n_thin,
                        arma::vec& Lambda, arma::vec& Lambda_S,
                        arma::ivec& N_mis, arma::vec& c_S, arma::vec& c,
                        arma::imat& k, arma::ivec& m,
                        const arma::mat& mu_initial,
                        const arma::cube& Sigma_initial,
                        const char* site_source, const char* site_layer,
                        const char* sampled_source, const char* sampled_layer,
                        const int& n_mcint, const bool& keep_mis){

  // Timers.
  clock_t tick1, tick2;
  clock_t integrating = 0;
  clock_t imputing = 0;
  clock_t mixing = 0;
  clock_t updating = 0;

  // Wrapper for R's message().
  Rcpp::Function msg("message");

  // Intermediate variables used in computations.
  int n = x.n_cols;
  int d = x.n_rows;
  arma::mat Psi_inv = arma::inv_sympd(Psi);
  arma::vec Psi_inv_eta = Psi_inv * eta;
  arma::mat Omega_inv = arma::inv_sympd(Omega);
  double c_0 = 2.0 * nu - d + 1.0;
  arma::mat B_0 = 2.0 / c_0 * (Psi + Omega); // Check this.

  int h_used = arma::max(k.row(0)) + 1; // Number of components used.
  int h_next = h_used; // Label for next new component.
  int k0; // Index for accessing k_cdf.
  arma::vec k_cdf(n); // Cumulative distribution of k.
  arma::vec q; // Complete conditional mixing weights.
  arma::vec sum_x(d); // Sum of component members.
  arma::mat sum_sq(d, d); //Component sum of squares.
  arma::uvec unused; // Indices of unused component labels.
  arma::uvec members; // Indices of component members.
  arma::uvec members_mis; // Indices of component imputed members.
  arma::vec dev(d); // Devation of an observation from the component center.
  arma::mat mu_cov_inv(d, d); // Current precision of mu.
  arma::mat mu_cov(d, d); // Current covariance of mu.
  arma::cube Sigma_inv; // Current precisions.

  int l_new; // Component label of new point.
  arma::mat x_mc(d, n_mcint); // Points for MC integration.
  arma::vec h_predict; // Component sizes for the predictive distribution.
  double k_mc_quant; // Label quantile for the MC integration.
  arma::vec in_site(n_mcint); // In-site indicator vector.
  arma::vec in_sampled(n_mcint); // In-sampled-region indicator vector.
  OGRPoint test_x; // Point to test for intersections.
  OGRPoint temp_x; // Temporary storage for points.
  arma::vec mu_new(d); // Temporary center.
  arma::mat Sigma_new(d, d); // Temporary precision.
  arma::mat Sigma_inv_new(d, d); // Temporary covariance.
  bool discard; // Indicates if draw should be discarded or not.

  // Tracking for thinning.
  bool thinning = n_thin > 1;
  int buf_size = 2;
  if(thinning){
    buf_size = n_thin;
  }
  int n_keep = (n_iter - n_burnin) / n_thin +
                  ((n_iter - n_burnin) % n_thin > 0);
  int write_pos = 0;

  // Initialize storage. (Lambda, k, and m passed by reference.)
  arma::field<arma::mat> mu(n_keep);      // Component means.
  arma::field<arma::cube> Sigma(n_keep);  // Component covariances.
  arma::field<arma::ivec> h(n_keep);      // Component sizes.
  arma::field<arma::mat> X_mis(n_keep);   // Imputed points.
  arma::field<arma::mat> k_mis(n_keep);   // Imputed points.
  m.zeros();                              // Number of components.

  // Initialize buffers.
  arma::vec Lambda_buf(buf_size);
  arma::vec Lambda_S_buf(buf_size);
  arma::ivec N_mis_buf(buf_size);
  arma::vec c_S_buf(buf_size);
  arma::vec c_buf(buf_size);
  arma::imat k_buf(buf_size, n);
  arma::ivec m_buf(buf_size);
  arma::field<arma::mat> mu_buf(buf_size);
  arma::field<arma::cube> Sigma_buf(buf_size);
  arma::field<arma::ivec> h_buf(buf_size);
  arma::field<arma::mat> X_mis_buf(buf_size);
  arma::field<arma::mat> k_mis_buf(buf_size);
  arma::mat k_quant = arma::randu(buf_size, n); // Quantiles for choosing k.

  // Load site shapefile.
  GDALAllRegister();
  GDALDataset *site_src;
  site_src = (GDALDataset *) GDALOpenEx(site_source, GDAL_OF_VECTOR,
                                        NULL, NULL, NULL);
  OGRLayer *site_shp;
  site_shp = site_src->GetLayerByName(site_layer);
  site_shp->ResetReading();
  OGRFeature *site_feature;
  site_feature = site_shp->GetNextFeature(); // Only using first feature...
  OGRPolygon *site_geom;
  site_geom = (OGRPolygon *) site_feature->GetGeometryRef();
  int site_nholes = site_geom->getNumInteriorRings();
  OGRLinearRing *ring;
  ring = site_geom->getExteriorRing();

  // Create arrays of site exterior ring vertices.
  int site_ext_npolys = 1;
  int site_ext_npt = ring->getNumPoints();
  double *site_ext_x[site_ext_npolys];
  double *site_ext_y[site_ext_npolys];
  double *thisx = new double [site_ext_npt];
  double *thisy = new double [site_ext_npt];
  for(int idx = 0; idx < site_ext_npt; idx++){
    ring->getPoint(idx, &temp_x);
    thisx[idx] = temp_x.getX();
    thisy[idx] = temp_x.getY();
  }
  site_ext_x[0] = thisx;
  site_ext_y[0] = thisy;

  // Create arrays of site interior ring vertices.
  int site_nholes_offset = 0;
  int site_int_npt[site_nholes];
  double *site_int_x[site_nholes];
  double *site_int_y[site_nholes];
  for(int hole = 0; hole < site_nholes; hole++){
    ring = site_geom->getInteriorRing(hole);
    site_int_npt[hole] = ring->getNumPoints();
    thisx = new double [site_int_npt[hole]];
    thisy = new double [site_int_npt[hole]];
    for(int idx = 0; idx < site_int_npt[hole]; idx++){
      ring->getPoint(idx, &temp_x);
      thisx[idx] = temp_x.getX();
      thisy[idx] = temp_x.getY();
    }
    site_int_x[hole] = thisx;
    site_int_y[hole] = thisy;
  }

  GDALDataset *sampled_src;
  sampled_src = (GDALDataset *) GDALOpenEx(sampled_source, GDAL_OF_VECTOR,
                                           NULL, NULL, NULL);
  OGRLayer *sampled_shp;
  sampled_shp = sampled_src->GetLayerByName(sampled_layer);
  sampled_shp->ResetReading();
  OGRFeature *sampled_feature;
  sampled_feature = sampled_shp->GetNextFeature();
  OGRMultiPolygon *sampled_geom;
  sampled_geom = (OGRMultiPolygon *) sampled_feature->GetGeometryRef();
  int sampled_npolys = sampled_geom->getNumGeometries();
  OGRPolygon *sampled_poly;

  // Create arrays of sampled exterior ring vertices.
  int sampled_ext_npt[sampled_npolys];
  int sampled_nholes[sampled_npolys];
  int sampled_nholes_offset[sampled_npolys];
  sampled_nholes_offset[0] = 0;
  int sampled_nholes_cum[sampled_npolys];
  double *sampled_ext_x[sampled_npolys];
  double *sampled_ext_y[sampled_npolys];
  for(int poly = 0; poly < sampled_npolys; poly++){
    sampled_poly = (OGRPolygon *) sampled_geom->getGeometryRef(poly);
    ring = sampled_poly->getExteriorRing();
    sampled_nholes[poly] = sampled_poly->getNumInteriorRings();
    sampled_nholes_cum[poly] = sampled_nholes_offset[poly] +
                                 sampled_nholes[poly];
    if(poly < sampled_npolys){
      sampled_nholes_offset[poly + 1] = sampled_nholes_cum[poly];
    }
    sampled_ext_npt[poly] = ring->getNumPoints();
    thisx = new double [sampled_ext_npt[poly]];
    thisy = new double [sampled_ext_npt[poly]];
    for(int idx = 0; idx < sampled_ext_npt[poly]; idx++){
      ring->getPoint(idx, &temp_x);
      thisx[idx] = temp_x.getX();
      thisy[idx] = temp_x.getY();
    }
    sampled_ext_x[poly] = thisx;
    sampled_ext_y[poly] = thisy;
  }
  int sampled_nholes_tot = sampled_nholes_cum[sampled_npolys - 1];

  // Create arrays of sampled interior ring vertices.
  int sampled_int_npt[sampled_nholes_tot];
  double *sampled_int_x[sampled_nholes_tot];
  double *sampled_int_y[sampled_nholes_tot];
  for(int poly = 0; poly < sampled_npolys; poly++){
    sampled_poly = (OGRPolygon *) sampled_geom->getGeometryRef(poly);
    for(int hole = sampled_nholes_offset[poly];
        hole < sampled_nholes_cum[poly]; hole++){
      ring = sampled_poly->getInteriorRing(hole - sampled_nholes_offset[poly]);
      sampled_int_npt[hole] = ring->getNumPoints();
      thisx = new double [sampled_int_npt[hole]];
      thisy = new double [sampled_int_npt[hole]];
      for(int idx = 0; idx < sampled_int_npt[hole]; idx++){
        ring->getPoint(idx, &temp_x);
        thisx[idx] = temp_x.getX();
        thisy[idx] = temp_x.getY();
      }
      sampled_int_x[hole] = thisx;
      sampled_int_y[hole] = thisy;
    }
  }

  // Close shapefiles.
  OGRFeature::DestroyFeature(site_feature);
  GDALClose(site_src);
  OGRFeature::DestroyFeature(sampled_feature);
  GDALClose(sampled_src);

  // Initialize temporary varibles used by inRegions().
  double x1;
  double y1;
  double x2;
  double y2;
  int site_counter;
  int samp_counter;
  int in_i;
  int in_j;
  int in_pidx;
  int in_hidx;

  // In the initial iteration, initialize mu, Sigma, and h based on the
  // initial values provided in k.
  int j = buf_size - 1;
  int j_prev = buf_size - 1;

  // Look like we're doing something.
  msg(chain_id + std::string(": iteration 1 of ") + std::to_string(n_iter));

  // Create a vector of component sizes.
  h_buf(j).set_size(h_used);

  // Set mu and Sigma to their initial values.
  mu_buf(j) = mu_initial;
  Sigma_buf(j) = Sigma_initial;

  // Create a cube of component precisions.
  Sigma_inv.set_size(d, d, Sigma_buf(j).n_slices);

  // Loop through components.
  for(int l = 0; l < h_used; l++){

    // Update component size.
    members = find(k_buf.row(j) == l);
    h_buf(j)(l) = members.n_elem;

    // Update component precision.
    Sigma_inv.slice(l) = arma::inv_sympd(Sigma_buf(j).slice(l));
  }

  // Find unused labels and update the number of components used.
  unused = find(h_buf(j) < 1);
  m_buf(j) = h_used - unused.n_elem;

  // MC integration.
  h_predict = arma::vec(h_used + 1);
  h_predict.head(h_used) = arma::conv_to<arma::vec>::from(h_buf(j));
  h_predict(h_used) = alpha;
  k_cdf = arma::cumsum(arma::normalise(h_predict, 1));
  for(int idx = 0; idx < n_mcint; idx++){
    k_mc_quant = R::unif_rand();  // Quantile for choosing k.
    l_new = 0;
    while(k_mc_quant > k_cdf(l_new)){
      l_new++;
    }
    if(l_new < h_used){
      x_mc.col(idx) = randmvn(mu_buf(j).col(l_new), Sigma_buf(j).slice(l_new));
    } else {
      x_mc.col(idx) = randmvn(
        randmvn(eta, Psi),
        arma::inv_sympd(randw(nu, Omega_inv))
      );
    }
    site_counter = 0;
    in_site(idx) = inRegions(x_mc(0, idx), x_mc(1, idx),
                             site_ext_x, site_ext_y,
                             site_ext_npolys, &site_ext_npt,
                             site_int_x, site_int_y, site_int_npt,
                             &site_nholes_offset, &site_nholes,
                             x1, y1, x2, y2,
                             site_counter, in_i, in_j, in_pidx, in_hidx);
    samp_counter = 0;
    in_sampled(idx) = inRegions(x_mc(0, idx), x_mc(1, idx),
                                sampled_ext_x, sampled_ext_y,
                                sampled_npolys, sampled_ext_npt,
                                sampled_int_x, sampled_int_y, sampled_int_npt,
                                sampled_nholes_offset, sampled_nholes_cum,
                                x1, y1, x2, y2,
                                samp_counter, in_i, in_j, in_pidx, in_hidx);
  }

  c_buf(j) = arma::mean(in_site);
  c_S_buf(j) = arma::mean(in_sampled);

  // Update Lambda.
  Lambda_S_buf(j) = Lambda_buf(j) * c_S_buf(j);

  // Impute.
  N_mis_buf(j) = R::rpois(Lambda_buf(j) * (1 - c_S_buf(j)));
  if(N_mis_buf(j) < 0){
    N_mis_buf(j) = 0;
  }
  X_mis_buf(j).set_size(d, N_mis_buf(j));
  k_mis_buf(j).set_size(N_mis_buf(j));
  if(N_mis_buf(j) > 0){
    for(int idx = 0; idx < N_mis_buf(j); idx++){
      discard = TRUE;
      while(discard){
        k_mc_quant = R::unif_rand();  // Quantile for choosing k.
        l_new = 0;
        while(k_mc_quant > k_cdf(l_new)){
          l_new++;
        }
        if(l_new < h_used){
          X_mis_buf(j).col(idx) = randmvn(
            mu_buf(j).col(l_new),
            Sigma_buf(j).slice(l_new)
          );
        } else {
          mu_new = randmvn(eta, Psi);
          Sigma_inv_new = randw(nu, Omega_inv);
          Sigma_new = arma::inv_sympd(Sigma_inv_new);
          X_mis_buf(j).col(idx) = randmvn(mu_new, Sigma_new);
        }
        site_counter = 0;
        samp_counter = 0;
        discard = inRegions(X_mis_buf(j)(0, idx), X_mis_buf(j)(1, idx),
                            sampled_ext_x, sampled_ext_y,
                            sampled_npolys, sampled_ext_npt,
                            sampled_int_x, sampled_int_y, sampled_int_npt,
                            sampled_nholes_offset, sampled_nholes_cum,
                            x1, y1, x2, y2,
                            samp_counter, in_i, in_j, in_pidx, in_hidx) ||
                    !(inRegions(X_mis_buf(j)(0, idx), X_mis_buf(j)(1, idx),
                                site_ext_x, site_ext_y,
                                site_ext_npolys, &site_ext_npt,
                                site_int_x, site_int_y, site_int_npt,
                                &site_nholes_offset, &site_nholes,
                                x1, y1, x2, y2,
                                site_counter, in_i, in_j, in_pidx, in_hidx));
      }
      if(l_new == h_used){
        if(unused.n_elem > 0){
          k_mis_buf(j)(idx) = unused(0);
          h_buf(j)(unused(0)) = 1;
          unused.shed_row(0);
          mu_buf(j).col(k_mis_buf(j)(idx)) = mu_new;
          Sigma_buf(j).slice(k_mis_buf(j)(idx)) = Sigma_new;
          Sigma_inv.slice(k_mis_buf(j)(idx)) = Sigma_inv_new;
        } else {
          k_mis_buf(j)(idx) = h_next;
          h_buf(j) = arma::join_cols(h_buf(j), arma::ivec(1, arma::fill::ones));
          mu_buf(j) = arma::join_rows(mu_buf(j), mu_new);
          Sigma_buf(j) = arma::join_slices(Sigma_buf(j), Sigma_new);
          Sigma_inv = arma::join_slices(Sigma_inv, Sigma_inv_new);
          h_next++;
        }
      } else {
        k_mis_buf(j)(idx) = l_new;
        h_buf(j)(l_new)++;
      }
    }
  }

  // Update number of component labels in use.
  h_used = h_buf(j).n_elem;

  // Flush buffer if necessary.
  if(n_burnin < 1){
    Lambda(write_pos) = Lambda_buf(j);
    Lambda_S(write_pos) = Lambda_S_buf(j);
    N_mis(write_pos) = N_mis_buf(j);
    c_S(write_pos) = c_S_buf(j);
    c(write_pos) = c_buf(j);
    k.row(write_pos) = k_buf.row(j);
    m(write_pos) = m_buf(j);
    mu(write_pos) = mu_buf(j);
    Sigma(write_pos) = Sigma_buf(j);
    h(write_pos) = h_buf(j);
    if(keep_mis){
      X_mis(write_pos) = X_mis_buf(j);
      k_mis(write_pos) = k_mis_buf(j);
    }
    write_pos++;
  }

  // Go to beginning of buffer.
  j = 0;

  tick1 = clock();

  // Loop through iterations.
  for(int iter_ctr = 1; iter_ctr < n_iter; iter_ctr++){

    // Look like we're doing something.
    msg(chain_id + std::string(": iteration ") + std::to_string(iter_ctr + 1) +
      std::string(" of ") + std::to_string(n_iter));

    h_next = h_used;

    // Loop through observations.
    for(int i = 0; i < n; i++){

      // Compute mixing weights for conditional posteriors.
      q.set_size(h_used + 1);
      for(int idx = 0; idx < h_used; idx++){
        if(h_buf(j_prev)(idx) > 0){
          q(idx) = (h_buf(j_prev)(idx) - (idx == k_buf(j_prev, i))) *
                     densmvn(x.col(i),
                             mu_buf(j_prev).col(idx),
                             Sigma_buf(j_prev).slice(idx));
        } else {
          q(idx) = 0;
        }
      }
      q(h_used) = alpha * densmvt(x.col(i), c_0, eta, B_0);
      // To understand the above, think "prior predictive density".
      k_cdf = arma::cumsum(arma::normalise(q, 1));

      // Draw a component label.
      l_new = 0;
      while(k_quant(j, i) > k_cdf(l_new)){
        l_new++;
      }

      // Update observation's component label.
      if(l_new == h_used){
        if(unused.n_elem > 0){
          k_buf(j, i) = unused(0);
          unused.shed_row(0);
          Sigma_inv.slice(k_buf(j, i)) = randw(nu, Omega_inv);
        } else {
          k_buf(j, i) = h_next;
          Sigma_inv = arma::join_slices(Sigma_inv, randw(nu, Omega_inv));
          h_next++;
        }
      } else {
        k_buf(j, i) = l_new;
      }
    }
    for(int i = 0; i < N_mis_buf(j_prev); i++){

      // Compute mixing weights for conditional posteriors.
      q.set_size(h_used + 1);
      for(int idx = 0; idx < h_used; idx++){
        if(h_buf(j_prev)(idx) > 0){
          q(idx) = (h_buf(j_prev)(idx) - (idx == k_mis_buf(j_prev)(i))) *
                     densmvn(X_mis_buf(j_prev).col(i),
                             mu_buf(j_prev).col(idx),
                             Sigma_buf(j_prev).slice(idx));
        } else {
          q(idx) = 0;
        }
      }
      q(h_used) = alpha * densmvt(X_mis_buf(j_prev).col(i), c_0, eta, B_0);
      // To understand the above, think "prior predictive density".
      k_cdf = arma::cumsum(arma::normalise(q, 1));

      // Draw a component label.
      l_new = 0;
      k_mc_quant = R::unif_rand();
      while(k_mc_quant > k_cdf(l_new)){
        l_new++;
      }

      // Update observation's component label.
      if(l_new == h_used){
        if(unused.n_elem > 0){
          k_mis_buf(j_prev)(i) = unused(0);
          unused.shed_row(0);
          Sigma_inv.slice(k_mis_buf(j_prev)(i)) = randw(nu, Omega_inv);
        } else {
          k_mis_buf(j_prev)(i) = h_next;
          Sigma_inv = arma::join_slices(Sigma_inv, randw(nu, Omega_inv));
          h_next++;
        }
      } else {
        k_mis_buf(j_prev)(i) = l_new;
      }
    }

    tick2 = clock();
    mixing += tick2 - tick1;

    // Update number of component labels in use.
    h_used = h_next;

    // Create a vector of component sizes.
    h_buf(j).set_size(h_used);

    // Create a matrix of component centers.
    mu_buf(j).set_size(d, h_used);

    // Create a cube of component covariances.
    Sigma_buf(j).set_size(d, d, h_used);

    // Loop through components.
    for(int l = 0; l < h_used; l++){

      // Update component size.
      members = find(k_buf.row(j) == l);
      members_mis = find(k_mis_buf(j_prev) == l);
      h_buf(j)(l) = members.n_elem + members_mis.n_elem;

      if(h_buf(j)(l) > 0){

        // Update component sum and covariance of center.
        sum_x = sum(x.cols(members), 1) +
                  sum(X_mis_buf(j_prev).cols(members_mis), 1);
        mu_cov_inv = Psi_inv + h_buf(j)(l) * Sigma_inv.slice(l);
        mu_cov = arma::inv_sympd(mu_cov_inv);

        // Draw a new center.
        mu_buf(j).col(l) = randmvn(mu_cov * (Psi_inv_eta +
                                     Sigma_inv.slice(l) * sum_x),
                                   mu_cov);

        // Update component sum of squares.
        sum_sq.zeros();
        for(int i = 0; i < members.n_elem; i++){
          dev = x.col(members(i)) - mu_buf(j).col(l);
          sum_sq += dev * arma::trans(dev);
        }
        for(int i = 0; i < members_mis.n_elem; i++){
          dev = X_mis_buf(j_prev).col(members_mis(i)) - mu_buf(j).col(l);
          sum_sq += dev * arma::trans(dev);
        }

        // Draw a new covariance.
        Sigma_inv.slice(l) = randw(nu + h_buf(j)(l),
                               arma::inv_sympd(Omega + sum_sq));
        Sigma_buf(j).slice(l) = arma::inv_sympd(Sigma_inv.slice(l));

      } else {
        mu_buf(j).col(l).fill(arma::datum::nan);
        Sigma_buf(j).slice(l).fill(arma::datum::nan);
        Sigma_inv.slice(l).fill(arma::datum::nan);
      }
    }

    // Find unused labels and update the number of components used.
    unused = find(h_buf(j) < 1);
    m_buf(j) = h_used - unused.n_elem;

    tick1 = clock();
    updating += tick1 - tick2;

    // MC integration.
    h_predict = arma::vec(h_used + 1);
    h_predict.head(h_used) = arma::conv_to<arma::vec>::from(h_buf(j));
    h_predict(h_used) = alpha;
    k_cdf = arma::cumsum(arma::normalise(h_predict, 1));
    for(int idx = 0; idx < n_mcint; idx++){
      k_mc_quant = R::unif_rand();  // Quantile for choosing k.
      l_new = 0;
      while(k_mc_quant > k_cdf(l_new)){
        l_new++;
      }
      if(l_new < h_used){
        x_mc.col(idx) = randmvn(
          mu_buf(j).col(l_new),
          Sigma_buf(j).slice(l_new)
        );
      } else {
        x_mc.col(idx) = randmvn(
          randmvn(eta, Psi),
          arma::inv_sympd(randw(nu, Omega_inv))
        );
      }
      site_counter = 0;
      in_site(idx) = inRegions(x_mc(0, idx), x_mc(1, idx),
                               site_ext_x, site_ext_y,
                               site_ext_npolys, &site_ext_npt,
                               site_int_x, site_int_y, site_int_npt,
                               &site_nholes_offset, &site_nholes,
                               x1, y1, x2, y2,
                               site_counter, in_i, in_j, in_pidx, in_hidx);
      samp_counter = 0;
      in_sampled(idx) = inRegions(x_mc(0, idx), x_mc(1, idx),
                                  sampled_ext_x, sampled_ext_y,
                                  sampled_npolys, sampled_ext_npt,
                                  sampled_int_x, sampled_int_y,
                                  sampled_int_npt,
                                  sampled_nholes_offset, sampled_nholes_cum,
                                  x1, y1, x2, y2,
                                  samp_counter, in_i, in_j, in_pidx, in_hidx);
    }

    c_buf(j) = arma::mean(in_site);
    c_S_buf(j) = arma::mean(in_sampled);

    tick2 = clock();
    integrating += tick2 - tick1;

    // Update Lambda.
    Lambda_buf(j) = (arma::randg(1,
                     arma::distr_param(a + n, 1.0 / (b + c_S_buf(j)))))(0);
    Lambda_S_buf(j) = Lambda_buf(j) * c_S_buf(j);

    // Impute.
    // Note: This is predicting the next point many times,
    // NOT simulating the Dirichlet process. Should it simulate the DP?
    N_mis_buf(j) = R::rpois(Lambda_buf(j) * (1 - c_S_buf(j)));
    if(N_mis_buf(j) < 0){
      N_mis_buf(j) = 0;
    }
    X_mis_buf(j).set_size(d, N_mis_buf(j));
    k_mis_buf(j).set_size(N_mis_buf(j));
    if(N_mis_buf(j) > 0){
      for(int idx = 0; idx < N_mis_buf(j); idx++){
        discard = TRUE;
        while(discard){
          k_mc_quant = R::unif_rand();  // Quantile for choosing k.
          l_new = 0;
          while(k_mc_quant > k_cdf(l_new)){
            l_new++;
          }
          if(l_new < h_used){
            X_mis_buf(j).col(idx) = randmvn(
              mu_buf(j).col(l_new),
              Sigma_buf(j).slice(l_new)
            );
          } else {
            mu_new = randmvn(eta, Psi);
            Sigma_inv_new = randw(nu, Omega_inv);
            Sigma_new = arma::inv_sympd(Sigma_inv_new);
            X_mis_buf(j).col(idx) = randmvn(mu_new, Sigma_new);
          }
          site_counter = 0;
          samp_counter = 0;
          discard = inRegions(X_mis_buf(j)(0, idx), X_mis_buf(j)(1, idx),
                              sampled_ext_x, sampled_ext_y,
                              sampled_npolys, sampled_ext_npt,
                              sampled_int_x, sampled_int_y, sampled_int_npt,
                              sampled_nholes_offset, sampled_nholes_cum,
                              x1, y1, x2, y2,
                              samp_counter, in_i, in_j, in_pidx, in_hidx) ||
                      !(inRegions(X_mis_buf(j)(0, idx), X_mis_buf(j)(1, idx),
                                  site_ext_x, site_ext_y,
                                  site_ext_npolys, &site_ext_npt,
                                  site_int_x, site_int_y, site_int_npt,
                                  &site_nholes_offset, &site_nholes,
                                  x1, y1, x2, y2,
                                  site_counter, in_i, in_j,
                                  in_pidx, in_hidx));
        }
        if(l_new == h_used){
          if(unused.n_elem > 0){
            k_mis_buf(j)(idx) = unused(0);
            unused.shed_row(0);
            mu_buf(j).col(k_mis_buf(j)(idx)) = mu_new;
            Sigma_buf(j).slice(k_mis_buf(j)(idx)) = Sigma_new;
            Sigma_inv.slice(k_mis_buf(j)(idx)) = Sigma_inv_new;
          } else {
            k_mis_buf(j)(idx) = h_next;
            mu_buf(j) = arma::join_rows(mu_buf(j), mu_new);
            Sigma_buf(j) = arma::join_slices(Sigma_buf(j), Sigma_new);
            Sigma_inv = arma::join_slices(Sigma_inv, Sigma_inv_new);
            h_next++;
          }
        } else {
          k_mis_buf(j)(idx) = l_new;
        }
      }
    }

    tick1 = clock();
    imputing += tick1 - tick2;

    // Flush buffer.
    if(!thinning ||
         (iter_ctr >= n_burnin && (iter_ctr - n_burnin) % n_thin == 0)){
      Lambda(write_pos) = Lambda_buf(j);
      Lambda_S(write_pos) = Lambda_S_buf(j);
      N_mis(write_pos) = N_mis_buf(j);
      c_S(write_pos) = c_S_buf(j);
      c(write_pos) = c_buf(j);
      k.row(write_pos) = k_buf.row(j);
      m(write_pos) = m_buf(j);
      mu(write_pos) = mu_buf(j);
      Sigma(write_pos) = Sigma_buf(j);
      h(write_pos) = h_buf(j);
      if(keep_mis){
        X_mis(write_pos) = X_mis_buf(j);
        k_mis(write_pos) = k_mis_buf(j);
      }
      write_pos++;
    }

    // Advance index in buffer.
    j_prev = j;
    j++;

    // Check if we should go back to the beginning of the buffer.
    if(j >= buf_size){
      j = 0;

      // Draw new k quantiles.
      k_quant = arma::randu(buf_size, n);
    }
  }

  msg(chain_id + std::string(": Done\n  ") +
    std::to_string(integrating * 1000 / CLOCKS_PER_SEC) +
    std::string(" ms integrating\n  ") +
    std::to_string(imputing * 1000 / CLOCKS_PER_SEC) +
    std::string(" ms imputing\n  ") +
    std::to_string(mixing * 1000 / CLOCKS_PER_SEC) +
    std::string(" ms updating mixture\n  ") +
    std::to_string(updating * 1000 / CLOCKS_PER_SEC) +
    std::string(" ms updating components"));

  // Free up memory used by vertices.
  delete [] site_ext_x[0];
  delete [] site_ext_y[0];
  for(int hole = 0; hole < site_nholes; hole++){
    delete [] site_int_x[hole];
    delete [] site_int_y[hole];
  }
  for(int poly = 0; poly < sampled_npolys; poly++){
    delete [] sampled_ext_x[poly];
    delete [] sampled_ext_y[poly];
  }
  for(int hole = 0; hole < sampled_nholes_tot; hole++){
    delete [] sampled_int_x[hole];
    delete [] sampled_int_y[hole];
  }

  return Rcpp::List::create(Rcpp::Named("mu") = Rcpp::wrap(mu),
                            Rcpp::Named("Sigma") = Rcpp::wrap(Sigma),
                            Rcpp::Named("h") = Rcpp::wrap(h),
                            Rcpp::Named("X_mis") = Rcpp::wrap(X_mis),
                            Rcpp::Named("k_mis") = Rcpp::wrap(k_mis));
}
// }}}

// Partially observed, bounded, no augmentation. {{{
// [[Rcpp::export]]
Rcpp::List mcmc_chain_sn(const std::string& chain_id, const arma::mat& x,
                         const double& a, const double& b,
                         const arma::vec& eta, const arma::mat& Psi,
                         const int& nu, const arma::mat& Omega,
                         const double& alpha, const int& n_iter,
                         const int& n_burnin, const int& n_thin,
                         arma::vec& Lambda, arma::vec& Lambda_S,
                         arma::vec& c_S, arma::vec& c,
                         arma::imat& k, arma::ivec& m,
                         const arma::mat& mu_initial,
                         const arma::cube& Sigma_initial,
                         const char* site_source, const char* site_layer,
                         const char* sampled_source, const char* sampled_layer,
                         const int& n_mcint){

  // Timers.
  clock_t tick1, tick2;
  clock_t integrating = 0;
  clock_t imputing = 0;
  clock_t mixing = 0;
  clock_t updating = 0;

  // Wrapper for R's message().
  Rcpp::Function msg("message");

  // Intermediate variables used in computations.
  int n = x.n_cols;
  int d = x.n_rows;
  arma::mat Psi_inv = arma::inv_sympd(Psi);
  arma::vec Psi_inv_eta = Psi_inv * eta;
  arma::mat Omega_inv = arma::inv_sympd(Omega);
  double c_0 = 2.0 * nu - d + 1.0;
  arma::mat B_0 = 2.0 / c_0 * (Psi + Omega); // Check this.

  int h_used = arma::max(k.row(0)) + 1; // Number of components used.
  int h_next = h_used; // Label for next new component.
  int k0; // Index for accessing k_cdf.
  arma::vec k_cdf(n); // Cumulative distribution of k.
  arma::vec q; // Complete conditional mixing weights.
  arma::vec sum_x(d); // Sum of component members.
  arma::mat sum_sq(d, d); //Component sum of squares.
  arma::uvec unused; // Indices of unused component labels.
  arma::uvec members; // Indices of component members.
  arma::vec dev(d); // Devation of an observation from the component center.
  arma::mat mu_cov_inv(d, d); // Current precision of mu.
  arma::mat mu_cov(d, d); // Current covariance of mu.
  arma::cube Sigma_inv; // Current precisions.

  int l_new; // Component label of new point.
  arma::mat x_mc(d, n_mcint); // Points for MC integration.
  arma::vec h_predict; // Component sizes for the predictive distribution.
  double k_mc_quant; // Label quantile for the MC integration.
  arma::vec in_site(n_mcint); // In-site indicator vector.
  arma::vec in_sampled(n_mcint); // In-sampled-region indicator vector.
  OGRPoint test_x; // Point to test for intersections.
  OGRPoint temp_x; // Temporary storage for points.
  arma::vec mu_new(d); // Temporary center.
  arma::mat Sigma_new(d, d); // Temporary precision.
  arma::mat Sigma_inv_new(d, d); // Temporary covariance.
  bool discard; // Indicates if draw should be discarded or not.

  // Tracking for thinning.
  bool thinning = n_thin > 1;
  int buf_size = 2;
  if(thinning){
    buf_size = n_thin;
  }
  int n_keep = (n_iter - n_burnin) / n_thin +
                  ((n_iter - n_burnin) % n_thin > 0);
  int write_pos = 0;

  // Initialize storage. (Lambda, k, and m passed by reference.)
  arma::field<arma::mat> mu(n_keep);      // Component means.
  arma::field<arma::cube> Sigma(n_keep);  // Component covariances.
  arma::field<arma::ivec> h(n_keep);      // Component sizes.
  m.zeros();                              // Number of components.

  // Initialize buffers.
  arma::vec Lambda_buf(buf_size);
  arma::vec Lambda_S_buf(buf_size);
  arma::vec c_S_buf(buf_size);
  arma::vec c_buf(buf_size);
  arma::imat k_buf(buf_size, n);
  arma::ivec m_buf(buf_size);
  arma::field<arma::mat> mu_buf(buf_size);
  arma::field<arma::cube> Sigma_buf(buf_size);
  arma::field<arma::ivec> h_buf(buf_size);
  arma::mat k_quant = arma::randu(buf_size, n); // Quantiles for choosing k.

  // Load site shapefile.
  GDALAllRegister();
  GDALDataset *site_src;
  site_src = (GDALDataset *) GDALOpenEx(site_source, GDAL_OF_VECTOR,
                                        NULL, NULL, NULL);
  OGRLayer *site_shp;
  site_shp = site_src->GetLayerByName(site_layer);
  site_shp->ResetReading();
  OGRFeature *site_feature;
  site_feature = site_shp->GetNextFeature(); // Only using first feature...
  OGRPolygon *site_geom;
  site_geom = (OGRPolygon *) site_feature->GetGeometryRef();
  int site_nholes = site_geom->getNumInteriorRings();
  OGRLinearRing *ring;
  ring = site_geom->getExteriorRing();

  // Create arrays of site exterior ring vertices.
  int site_ext_npolys = 1;
  int site_ext_npt = ring->getNumPoints();
  double *site_ext_x[site_ext_npolys];
  double *site_ext_y[site_ext_npolys];
  double *thisx = new double [site_ext_npt];
  double *thisy = new double [site_ext_npt];
  for(int idx = 0; idx < site_ext_npt; idx++){
    ring->getPoint(idx, &temp_x);
    thisx[idx] = temp_x.getX();
    thisy[idx] = temp_x.getY();
  }
  site_ext_x[0] = thisx;
  site_ext_y[0] = thisy;

  // Create arrays of site interior ring vertices.
  int site_nholes_offset = 0;
  int site_int_npt[site_nholes];
  double *site_int_x[site_nholes];
  double *site_int_y[site_nholes];
  for(int hole = 0; hole < site_nholes; hole++){
    ring = site_geom->getInteriorRing(hole);
    site_int_npt[hole] = ring->getNumPoints();
    thisx = new double [site_int_npt[hole]];
    thisy = new double [site_int_npt[hole]];
    for(int idx = 0; idx < site_int_npt[hole]; idx++){
      ring->getPoint(idx, &temp_x);
      thisx[idx] = temp_x.getX();
      thisy[idx] = temp_x.getY();
    }
    site_int_x[hole] = thisx;
    site_int_y[hole] = thisy;
  }

  GDALDataset *sampled_src;
  sampled_src = (GDALDataset *) GDALOpenEx(sampled_source, GDAL_OF_VECTOR,
                                           NULL, NULL, NULL);
  OGRLayer *sampled_shp;
  sampled_shp = sampled_src->GetLayerByName(sampled_layer);
  sampled_shp->ResetReading();
  OGRFeature *sampled_feature;
  sampled_feature = sampled_shp->GetNextFeature();
  OGRMultiPolygon *sampled_geom;
  sampled_geom = (OGRMultiPolygon *) sampled_feature->GetGeometryRef();
  int sampled_npolys = sampled_geom->getNumGeometries();
  OGRPolygon *sampled_poly;

  // Create arrays of sampled exterior ring vertices.
  int sampled_ext_npt[sampled_npolys];
  int sampled_nholes[sampled_npolys];
  int sampled_nholes_offset[sampled_npolys];
  sampled_nholes_offset[0] = 0;
  int sampled_nholes_cum[sampled_npolys];
  double *sampled_ext_x[sampled_npolys];
  double *sampled_ext_y[sampled_npolys];
  for(int poly = 0; poly < sampled_npolys; poly++){
    sampled_poly = (OGRPolygon *) sampled_geom->getGeometryRef(poly);
    ring = sampled_poly->getExteriorRing();
    sampled_nholes[poly] = sampled_poly->getNumInteriorRings();
    sampled_nholes_cum[poly] = sampled_nholes_offset[poly] +
                                 sampled_nholes[poly];
    if(poly < sampled_npolys){
      sampled_nholes_offset[poly + 1] = sampled_nholes_cum[poly];
    }
    sampled_ext_npt[poly] = ring->getNumPoints();
    thisx = new double [sampled_ext_npt[poly]];
    thisy = new double [sampled_ext_npt[poly]];
    for(int idx = 0; idx < sampled_ext_npt[poly]; idx++){
      ring->getPoint(idx, &temp_x);
      thisx[idx] = temp_x.getX();
      thisy[idx] = temp_x.getY();
    }
    sampled_ext_x[poly] = thisx;
    sampled_ext_y[poly] = thisy;
  }
  int sampled_nholes_tot = sampled_nholes_cum[sampled_npolys - 1];

  // Create arrays of sampled interior ring vertices.
  int sampled_int_npt[sampled_nholes_tot];
  double *sampled_int_x[sampled_nholes_tot];
  double *sampled_int_y[sampled_nholes_tot];
  for(int poly = 0; poly < sampled_npolys; poly++){
    sampled_poly = (OGRPolygon *) sampled_geom->getGeometryRef(poly);
    for(int hole = sampled_nholes_offset[poly];
        hole < sampled_nholes_cum[poly]; hole++){
      ring = sampled_poly->getInteriorRing(hole - sampled_nholes_offset[poly]);
      sampled_int_npt[hole] = ring->getNumPoints();
      thisx = new double [sampled_int_npt[hole]];
      thisy = new double [sampled_int_npt[hole]];
      for(int idx = 0; idx < sampled_int_npt[hole]; idx++){
        ring->getPoint(idx, &temp_x);
        thisx[idx] = temp_x.getX();
        thisy[idx] = temp_x.getY();
      }
      sampled_int_x[hole] = thisx;
      sampled_int_y[hole] = thisy;
    }
  }

  // Close shapefiles.
  OGRFeature::DestroyFeature(site_feature);
  GDALClose(site_src);
  OGRFeature::DestroyFeature(sampled_feature);
  GDALClose(sampled_src);

  // Initialize temporary varibles used by inRegions().
  double x1;
  double y1;
  double x2;
  double y2;
  int site_counter;
  int samp_counter;
  int in_i;
  int in_j;
  int in_pidx;
  int in_hidx;

  // In the initial iteration, initialize mu, Sigma, and h based on the
  // initial values provided in k.
  int j = buf_size - 1;
  int j_prev = buf_size - 1;

  // Look like we're doing something.
  msg(chain_id + std::string(": iteration 1 of ") + std::to_string(n_iter));

  // Create a vector of component sizes.
  h_buf(j).set_size(h_used);

  // Set mu and Sigma to their initial values.
  mu_buf(j) = mu_initial;
  Sigma_buf(j) = Sigma_initial;

  // Create a cube of component precisions.
  Sigma_inv.set_size(d, d, Sigma_buf(j).n_slices);

  // Loop through components.
  for(int l = 0; l < h_used; l++){

    // Update component size.
    members = find(k_buf.row(j) == l);
    h_buf(j)(l) = members.n_elem;

    // Update component precision.
    Sigma_inv.slice(l) = arma::inv_sympd(Sigma_buf(j).slice(l));
  }

  // Find unused labels and update the number of components used.
  unused = find(h_buf(j) < 1);
  m_buf(j) = h_used - unused.n_elem;

  // MC integration.
  h_predict = arma::vec(h_used + 1);
  h_predict.head(h_used) = arma::conv_to<arma::vec>::from(h_buf(j));
  h_predict(h_used) = alpha;
  k_cdf = arma::cumsum(arma::normalise(h_predict, 1));
  for(int idx = 0; idx < n_mcint; idx++){
    k_mc_quant = R::unif_rand();  // Quantile for choosing k.
    l_new = 0;
    while(k_mc_quant > k_cdf(l_new)){
      l_new++;
    }
    if(l_new < h_used){
      x_mc.col(idx) = randmvn(mu_buf(j).col(l_new), Sigma_buf(j).slice(l_new));
    } else {
      x_mc.col(idx) = randmvn(
        randmvn(eta, Psi),
        arma::inv_sympd(randw(nu, Omega_inv))
      );
    }
    site_counter = 0;
    in_site(idx) = inRegions(x_mc(0, idx), x_mc(1, idx),
                             site_ext_x, site_ext_y,
                             site_ext_npolys, &site_ext_npt,
                             site_int_x, site_int_y, site_int_npt,
                             &site_nholes_offset, &site_nholes,
                             x1, y1, x2, y2,
                             site_counter, in_i, in_j, in_pidx, in_hidx);
    samp_counter = 0;
    in_sampled(idx) = inRegions(x_mc(0, idx), x_mc(1, idx),
                                sampled_ext_x, sampled_ext_y,
                                sampled_npolys, sampled_ext_npt,
                                sampled_int_x, sampled_int_y, sampled_int_npt,
                                sampled_nholes_offset, sampled_nholes_cum,
                                x1, y1, x2, y2,
                                samp_counter, in_i, in_j, in_pidx, in_hidx);
  }

  c_buf(j) = arma::mean(in_site);
  c_S_buf(j) = arma::mean(in_sampled);

  // Update Lambda.
  Lambda_S_buf(j) = Lambda_buf(j) * c_S_buf(j);

  // Update number of component labels in use.
  h_used = h_buf(j).n_elem;

  // Flush buffer if necessary.
  if(n_burnin < 1){
    Lambda(write_pos) = Lambda_buf(j);
    Lambda_S(write_pos) = Lambda_S_buf(j);
    c_S(write_pos) = c_S_buf(j);
    c(write_pos) = c_buf(j);
    k.row(write_pos) = k_buf.row(j);
    m(write_pos) = m_buf(j);
    mu(write_pos) = mu_buf(j);
    Sigma(write_pos) = Sigma_buf(j);
    h(write_pos) = h_buf(j);
    write_pos++;
  }

  // Go to beginning of buffer.
  j = 0;

  tick1 = clock();

  // Loop through iterations.
  for(int iter_ctr = 1; iter_ctr < n_iter; iter_ctr++){

    // Look like we're doing something.
    msg(chain_id + std::string(": iteration ") + std::to_string(iter_ctr + 1) +
      std::string(" of ") + std::to_string(n_iter));

    h_next = h_used;

    // Loop through observations.
    for(int i = 0; i < n; i++){

      // Compute mixing weights for conditional posteriors.
      q.set_size(h_used + 1);
      for(int idx = 0; idx < h_used; idx++){
        if(h_buf(j_prev)(idx) > 0){
          q(idx) = (h_buf(j_prev)(idx) - (idx == k_buf(j_prev, i))) *
                     densmvn(x.col(i),
                             mu_buf(j_prev).col(idx),
                             Sigma_buf(j_prev).slice(idx));
        } else {
          q(idx) = 0;
        }
      }
      q(h_used) = alpha * densmvt(x.col(i), c_0, eta, B_0);
      // To understand the above, think "prior predictive density".
      k_cdf = arma::cumsum(arma::normalise(q, 1));

      // Draw a component label.
      l_new = 0;
      while(k_quant(j, i) > k_cdf(l_new)){
        l_new++;
      }

      // Update observation's component label.
      if(l_new == h_used){
        if(unused.n_elem > 0){
          k_buf(j, i) = unused(0);
          unused.shed_row(0);
          Sigma_inv.slice(k_buf(j, i)) = randw(nu, Omega_inv);
        } else {
          k_buf(j, i) = h_next;
          Sigma_inv = arma::join_slices(Sigma_inv, randw(nu, Omega_inv));
          h_next++;
        }
      } else {
        k_buf(j, i) = l_new;
      }
    }

    tick2 = clock();
    mixing += tick2 - tick1;

    // Update number of component labels in use.
    h_used = h_next;

    // Create a vector of component sizes.
    h_buf(j).set_size(h_used);

    // Create a matrix of component centers.
    mu_buf(j).set_size(d, h_used);

    // Create a cube of component covariances.
    Sigma_buf(j).set_size(d, d, h_used);

    // Loop through components.
    for(int l = 0; l < h_used; l++){

      // Update component size.
      members = find(k_buf.row(j) == l);
      h_buf(j)(l) = members.n_elem;

      if(h_buf(j)(l) > 0){

        // Update component sum and covariance of center.
        sum_x = sum(x.cols(members), 1);
        mu_cov_inv = Psi_inv + h_buf(j)(l) * Sigma_inv.slice(l);
        mu_cov = arma::inv_sympd(mu_cov_inv);

        // Draw a new center.
        mu_buf(j).col(l) = randmvn(mu_cov * (Psi_inv_eta +
                                     Sigma_inv.slice(l) * sum_x),
                                   mu_cov);

        // Update component sum of squares.
        sum_sq.zeros();
        for(int i = 0; i < members.n_elem; i++){
          dev = x.col(members(i)) - mu_buf(j).col(l);
          sum_sq += dev * arma::trans(dev);
        }

        // Draw a new covariance.
        Sigma_inv.slice(l) = randw(nu + h_buf(j)(l),
                               arma::inv_sympd(Omega + sum_sq));
        Sigma_buf(j).slice(l) = arma::inv_sympd(Sigma_inv.slice(l));

      } else {
        mu_buf(j).col(l).fill(arma::datum::nan);
        Sigma_buf(j).slice(l).fill(arma::datum::nan);
        Sigma_inv.slice(l).fill(arma::datum::nan);
      }
    }

    // Find unused labels and update the number of components used.
    unused = find(h_buf(j) < 1);
    m_buf(j) = h_used - unused.n_elem;

    tick1 = clock();
    updating += tick1 - tick2;

    // MC integration.
    h_predict = arma::vec(h_used + 1);
    h_predict.head(h_used) = arma::conv_to<arma::vec>::from(h_buf(j));
    h_predict(h_used) = alpha;
    k_cdf = arma::cumsum(arma::normalise(h_predict, 1));
    for(int idx = 0; idx < n_mcint; idx++){
      k_mc_quant = R::unif_rand();  // Quantile for choosing k.
      l_new = 0;
      while(k_mc_quant > k_cdf(l_new)){
        l_new++;
      }
      if(l_new < h_used){
        x_mc.col(idx) = randmvn(
          mu_buf(j).col(l_new),
          Sigma_buf(j).slice(l_new)
        );
      } else {
        x_mc.col(idx) = randmvn(
          randmvn(eta, Psi),
          arma::inv_sympd(randw(nu, Omega_inv))
        );
      }
      site_counter = 0;
      in_site(idx) = inRegions(x_mc(0, idx), x_mc(1, idx),
                               site_ext_x, site_ext_y,
                               site_ext_npolys, &site_ext_npt,
                               site_int_x, site_int_y, site_int_npt,
                               &site_nholes_offset, &site_nholes,
                               x1, y1, x2, y2,
                               site_counter, in_i, in_j, in_pidx, in_hidx);
      samp_counter = 0;
      in_sampled(idx) = inRegions(x_mc(0, idx), x_mc(1, idx),
                                  sampled_ext_x, sampled_ext_y,
                                  sampled_npolys, sampled_ext_npt,
                                  sampled_int_x, sampled_int_y,
                                  sampled_int_npt,
                                  sampled_nholes_offset, sampled_nholes_cum,
                                  x1, y1, x2, y2,
                                  samp_counter, in_i, in_j, in_pidx, in_hidx);
    }

    c_buf(j) = arma::mean(in_site);
    c_S_buf(j) = arma::mean(in_sampled);

    tick2 = clock();
    integrating += tick2 - tick1;

    // Update Lambda.
    Lambda_buf(j) = (arma::randg(1,
                     arma::distr_param(a + n, 1.0 / (b + c_S_buf(j)))))(0);
    Lambda_S_buf(j) = Lambda_buf(j) * c_S_buf(j);

    tick1 = clock();
    imputing += tick1 - tick2;

    // Flush buffer.
    if(!thinning ||
         (iter_ctr >= n_burnin && (iter_ctr - n_burnin) % n_thin == 0)){
      Lambda(write_pos) = Lambda_buf(j);
      Lambda_S(write_pos) = Lambda_S_buf(j);
      c_S(write_pos) = c_S_buf(j);
      c(write_pos) = c_buf(j);
      k.row(write_pos) = k_buf.row(j);
      m(write_pos) = m_buf(j);
      mu(write_pos) = mu_buf(j);
      Sigma(write_pos) = Sigma_buf(j);
      h(write_pos) = h_buf(j);
      write_pos++;
    }

    // Advance index in buffer.
    j_prev = j;
    j++;

    // Check if we should go back to the beginning of the buffer.
    if(j >= buf_size){
      j = 0;

      // Draw new k quantiles.
      k_quant = arma::randu(buf_size, n);
    }
  }

  msg(chain_id + std::string(": Done\n  ") +
    std::to_string(integrating * 1000 / CLOCKS_PER_SEC) +
    std::string(" ms integrating\n  ") +
    std::to_string(imputing * 1000 / CLOCKS_PER_SEC) +
    std::string(" ms imputing\n  ") +
    std::to_string(mixing * 1000 / CLOCKS_PER_SEC) +
    std::string(" ms updating mixture\n  ") +
    std::to_string(updating * 1000 / CLOCKS_PER_SEC) +
    std::string(" ms updating components"));

  // Free up memory used by vertices.
  delete [] site_ext_x[0];
  delete [] site_ext_y[0];
  for(int hole = 0; hole < site_nholes; hole++){
    delete [] site_int_x[hole];
    delete [] site_int_y[hole];
  }
  for(int poly = 0; poly < sampled_npolys; poly++){
    delete [] sampled_ext_x[poly];
    delete [] sampled_ext_y[poly];
  }
  for(int hole = 0; hole < sampled_nholes_tot; hole++){
    delete [] sampled_int_x[hole];
    delete [] sampled_int_y[hole];
  }

  return Rcpp::List::create(Rcpp::Named("mu") = Rcpp::wrap(mu),
                            Rcpp::Named("Sigma") = Rcpp::wrap(Sigma),
                            Rcpp::Named("h") = Rcpp::wrap(h));
}
// }}}

// vim: foldmethod=marker:

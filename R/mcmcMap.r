# Runs several Gibbs samplers to fit the inhomogeneous Poisson process model
# assuming a Dirichlet process mixture of Gaussian kernels for the intensity.
mcmcMap <- function(x, a, b, eta, Psi, nu, Omega, alpha = 1, n_chains = 1,
                    chain_ids = paste('Chain', seq_len(n_chains)),
                    k_initial = NULL, mu_initial = NULL, Sigma_initial = NULL,
                    Lambda_initial = NULL, Lambda_S_initial = NULL,
                    inits = NULL, n_iter = 1000, n_burnin = 0, n_thin = 1,
                    site_shp = NULL, sampled_shp = NULL,
                    n_mcint = 10000, sim_mis = TRUE, keep_mis = sim_mis,
                    cl = NULL){

  if(is.null(inits)){
    inits <- lapply(seq_len(n_chains), function(chain){
      return(list(chain_id = chain_ids[chain], k = k_initial[[chain]],
                  mu = mu_initial[[chain]], Sigma = Sigma_initial[[chain]],
                  Lambda = Lambda_initial[[chain]],
                  Lambda_S = Lambda_S_initial[[chain]]))
    })
  }

  if(inherits(cl, 'cluster') & requireNamespace('parallel', quietly = TRUE)){
    out <- parallel::parLapply(cl, inits,
                               mcmcMapChain, x,
                               a = a, b = b, eta = eta,
                               Psi = Psi, nu = nu, Omega = Omega,
                               alpha = alpha,
                               n_iter = n_iter, n_burnin = n_burnin,
                               n_thin = n_thin,
                               site_shp = site_shp, sampled_shp = sampled_shp,
                               n_mcint = n_mcint, sim_mis = sim_mis,
                               keep_mis = keep_mis)
  } else {
    out <- lapply(inits,
                  mcmcMapChain, x,
                  a = a, b = b, eta = eta,
                  Psi = Psi, nu = nu, Omega = Omega,
                  alpha = alpha,
                  n_iter = n_iter, n_burnin = n_burnin,
                  n_thin = n_thin,
                  site_shp = site_shp, sampled_shp = sampled_shp,
                  n_mcint = n_mcint, sim_mis = sim_mis, keep_mis = keep_mis)
  }
  return(structure(out, class = c('mcmcMap', 'mcmc.list')))
}

`[.mcmcMap` <- function(x, ...){
  cls <- class(x)
  x <- NextMethod('[', drop = FALSE)
  return(structure(x, class = cls))
}

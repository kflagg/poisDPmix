# Runs a single Gibbs sampler to fit the inhomogeneous Poisson process model
# assuming a Dirichlet process mixture of Gaussian kernels for the intensity.
mcmcMapChain <- function(init = list(chain_id = 'Chain', k = NULL, mu = NULL,
                                     Sigma = NULL, Lambda = NULL),
                         x, a, b, eta, Psi, nu, Omega, alpha = 1,
                         n_iter = 1000, n_burnin = 0, n_thin = 1,
                         site_shp = NULL, sampled_shp = NULL,
                         n_mcint = 10000, sim_mis = TRUE, keep_mis = sim_mis){

  # Flags to control which model we use.
  UNBOUNDED <- is.null(site_shp)
  SAMPLED <- !is.null(sampled_shp)
  AUGMENT <- SAMPLED & sim_mis

  # Intermediate variables used in computations.
  n <- nrow(x)
  d <- ncol(x)

  if(is.null(init$k)){
    init$k <- rep(NA_integer_, n)
    h <- integer(0)
    for(i in seq_len(n)){
      init$k[i] <- sample.int(length(h) + 1, 1, prob = c(h, alpha))
      if(init$k[i] > length(h)){
        h <- c(h, 1L)
      }else{
        h[init$k[i]] <- h[init$k[i]] + 1L
      }
    }
    init$k <- init$k - 1L
  }
  if(is.null(init$mu)){
    init$mu <- t(chol(Psi)) %*% matrix(rnorm((max(init$k) + 1) * d), nrow = d) +
                 matrix(eta, nrow = d, ncol = max(init$k) + 1)
  }
  if(is.null(init$Sigma)){
    init$Sigma <- structure(apply(rWishart(max(init$k) + 1, nu, solve(Omega)),
                                  3, solve),
                            dim = c(d, d, max(init$k) + 1))
  }
  if(is.null(init$Lambda)){
    init$Lambda <- rgamma(1, shape = a, rate = b)
  }

  # Initialize storage. (mu, Sigma, and h initialized in arma.)
  n_keep <- ceiling((n_iter - n_burnin) / n_thin)
  Lambda <- c(init$Lambda, rep(NA_real_, n_keep - 1)) # Mean count.
  if(SAMPLED){
    Lambda_S <- rep(NA_real_, n_keep)
    N_mis <- rep(NA_integer_, n_keep)  # Number of imputed points.
    c_S <- rep(NA_real_, n_keep)
  }
  c_norm <- rep(NA_real_, n_keep)  # Intensity normalising constant.
  k <- rbind(as.integer(init$k),   # Component labels.
             matrix(NA_integer_, nrow = n_keep - 1, ncol = n))
  m <- rep(NA_integer_, n_keep)    # Component sizes.

  # TODO: Find another way to do this without copying things into a list.
  call_args <- c(
    list(
      paste0('_poisDPmix_mcmc_chain',
             ifelse(UNBOUNDED, '_unb', ''),
             ifelse(SAMPLED, '_s', ''),
             ifelse(SAMPLED & !AUGMENT, 'n', '')),
      PACKAGE = 'poisDPmix',
      init$chain_id, t(x), a, b, eta, Psi, nu, Omega, alpha,
      n_iter, n_burnin, n_thin, Lambda
    ),
    if(SAMPLED) list(Lambda_S),
    if(AUGMENT) list(N_mis),
    if(SAMPLED) list(c_S),
    if(!UNBOUNDED) list(c_norm),
    list(
      k, m, init$mu, init$Sigma
    ),
    if(!UNBOUNDED) list(site_shp[1], site_shp[2]),
    if(SAMPLED) list(sampled_shp[1], sampled_shp[2]),
    if(!UNBOUNDED | SAMPLED) list(n_mcint),
    if(AUGMENT) list(keep_mis)
  )
  chain <- try(do.call(.Call, call_args))
  colnames(k) <- paste0('k_', seq_len(ncol(k)))

  # Assemble a coda::mcmc object with some extra attributes.
  return(structure(cbind(Lambda, Lambda_S = if(SAMPLED) Lambda_S else NULL,
                         N_mis = if(SAMPLED) N_mis else NULL,
                         c_S = if(SAMPLED) c_S else NULL,
                         c = if(UNBOUNDED) NULL else c_norm,
                         m, k),
                   mu = chain$mu,
                   Sigma = chain$Sigma,
                   h = lapply(chain$h,
                              function(x) ifelse(x > 0, x, NA_integer_)),
                   X_mis = if(SAMPLED) chain$X_mis
                     else NULL,
                   k_mis = if(SAMPLED) chain$k_mis
                     else NULL,
                   mcpar = c(n_burnin + 1,
                             n_burnin + 1 + (n_keep - 1) * n_thin,
                             n_thin),
                   class = c('mcmcMapChain', 'mcmc')))
}

`[.mcmcMapChain` <- function(x, i, ...){
  if(missing(i)) i <- seq_len(nrow(x))
  cls <- class(x)
  mu <- attr(x, 'mu')[i, drop = FALSE]
  Sigma <- attr(x, 'Sigma')[i, drop = FALSE]
  h <- attr(x, 'h')[i, drop = FALSE]
  X_mis <- attr(x, 'X_mis')[i, drop = FALSE]
  k_mis <- attr(x, 'k_mis')[i, drop = FALSE]
  mcpar <- attr(x, 'mcpar')
  x <- NextMethod('[')
  mcpar[2] <- mcpar[1] + mcpar[3] * (length(i) - 1)
  return(structure(x, mu = mu, Sigma = Sigma, h = h,
                   X_mis = X_mis, k_mis = k_mis,
                   mcpar = mcpar, class = cls))
}

# Simulates a point process from a fixed mixture, i.e. the number of
# components and mixing weights are specified.
# Arguments:
#  Lambda    Mean number of foreground anomalies.
#  p         Vector of mixing weights for the foreground process.
#  mu        length(p) by d matrix of foreground component means.
#  Sigma     length(p) by d by d array of foreground component covariances.
#  bound     A list of length d containing vectors of length 2 defining
#            intervals that bound the region.
#  win       An indicator function defining the region, if not a
#            (hyper-)rectangle. It must take an n by d matrix and return a
#            logical vector of length n where true means the row correponds to
#            a point inside the region.
#  Lambda_0  Mean number of background events.
simFixedMixture <- function(Lambda, p, mu, Sigma, bound, win, Lambda_0 = 0){

  d <- ncol(mu)
  in_bound <- function(y){return(rep(TRUE, nrow(y)))}
  if(nrow(mu) != length(p)){
    stop('mu must have dimension length(p) by d')
  }
  if(any(dim(Sigma) != c(length(p), d, d))){
    stop('sigma must have dimension length(p) by d by d')
  }
  if(missing(bound)){
    if(Lambda_0 > 0){
      stop('Lambda_0 > 0 but bound is missing')
    }
  } else {
    if(!(is.list(bound) & length(bound) == d)){
      stop('bound must be a list of length d')
    }
    in_bound <- function(y){apply(y, 1, function(y){
        return(all(y > sapply(bound, min) & y < sapply(bound, max)))
      })}
  }
  if(missing(win)){
    win <- function(y){return(rep(TRUE, nrow(y)))}
  }

  # Simulate uniform bg.
  bg_count <- rpois(1, Lambda_0)
  bg_bad <- rep(TRUE, bg_count)
  bg <- matrix(NA_real_, nrow = bg_count, ncol = d)
  while(any(bg_bad)){
    bg[bg_bad,] <- do.call(cbind, lapply(seq_len(d), function(i){
        return(runif(sum(bg_bad), min(bound[[i]]), max(bound[[i]])))
      }))
    bg_bad <- !(win(bg) & in_bound(bg))
  }

  # Simulate mixture fg.
  k <- length(p)
  fg_count <- rpois(1, Lambda)
  fg_labels <- sample.int(k, fg_count, replace = TRUE, prob = p)
  fg_bad <- rep(TRUE, fg_count)
  fg <- matrix(NA_real_, nrow = fg_count, ncol = d)
  while(any(fg_bad)){
    fg[fg_bad,] <- t(sapply(fg_labels[fg_bad], function(i){
      return(mvtnorm::rmvnorm(1, mu[i,], Sigma[i,,]))
    }))
    fg_bad <- !(win(fg) & in_bound(fg))
  }

  return(rbind(bg, fg))
}

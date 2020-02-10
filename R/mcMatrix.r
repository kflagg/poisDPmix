# Methods to  reformat ragged mcmcMap and mcmcMapChain objects
# into rectangular structures.
as.mcMatrix <- function(x, ...) UseMethod('as.mcMatrix')

as.mcMatrix.mcMatrix <- function(x){
  return(x)
}

as.mcMatrix.mcmcMap <- function(x){
  cls <- class(x)
  return(structure(lapply(x, as.mcMatrix), class = cls))
}

as.mcMatrix.mcmcMapChain <- function(x){

  # Get Attributes.
  mcpar <- attr(x, 'mcpar')
  cls <- class(x)
  cls <- cls[cls != 'mcmcMapChain']
  mu <- attr(x, 'mu')
  Sigma <- attr(x, 'Sigma')
  h <- attr(x, 'h')
  n_iter <- nrow(x)
  d <- length(dim(mu[[1]]))

  x <- as.matrix(x)

  # Get the label columns all anomalies.
  label_cols <- colnames(x)
  label_cols <- label_cols[grep('k_', label_cols)]

  # Get the h.
  temp_matrix <- t(sapply(seq_len(n_iter), function(iter){
    return(h[[iter]][x[iter, label_cols] + 1])
  }))
  colnames(temp_matrix) <- sub('k', 'h', label_cols)
  x <- cbind(x, temp_matrix)

  # Get the mu_l.
  for(l in seq_len(d)){
    temp_matrix <- t(sapply(seq_len(n_iter), function(iter){
      return(mu[[iter]][l, x[iter, label_cols] + 1])
    }))
    colnames(temp_matrix) <- sub('k_', paste0('mu_', l, ','), label_cols)
    x <- cbind(x, temp_matrix)
  }

  # Get the Sigma_m,l.
  for(l in seq_len(d)){
    for(m in seq_len(l)){
      temp_matrix <- t(sapply(seq_len(n_iter), function(iter){
        return(Sigma[[iter]][m, l, x[iter, label_cols] + 1])
      }))
      colnames(temp_matrix) <- sub('k_', paste0('Sigma_', m, ',' , l, ','),
                                   label_cols)
      x <- cbind(x, temp_matrix)
    }
  }

  return(structure(x, mcpar = mcpar, class = c('mcMatrix', cls)))
}

`[.mcMatrix` <- function(x, i, ...){
  if(missing(i)) i <- seq_len(nrow(x))
  cls <- class(x)
  mcpar <- attr(x, 'mcpar')
  x <- NextMethod('[')
  mcpar[2] <- mcpar[1] + mcpar[3] * (length(i) - 1)
  return(structure(x, mcpar = mcpar, class = cls))
}

# Posterior central ellipse methods for poisDPmix.
ellipses <- function(x, ...) UseMethod('ellipses')

ellipses.mcmcMap <- function(x, ..., add = TRUE, newplot = !add, type = 'l'){
  # Set up a window if !add.
  if(!add){
    # Choose a sensible default plotting window from the posterior.
    mu_1 <- unlist(lapply(do.call(c, lapply(x, attr, 'mu')),
                          function(m) m[1,]))
    sigma_1 <- unlist(lapply(do.call(c, lapply(x, attr, 'Sigma')),
                             function(s) s[1, 1,]))
    range_x <- c(
      min(mu_1 - 2 * sqrt(sigma_1), na.rm = TRUE),
      max(mu_1 + 2 * sqrt(sigma_1), na.rm = TRUE)
    )
    mu_2 <- unlist(lapply(do.call(c, lapply(x, attr, 'mu')),
                          function(m) m[2,]))
    sigma_2 <- unlist(lapply(do.call(c, lapply(x, attr, 'Sigma')),
                             function(s) s[2, 2,]))
    range_y <- c(
      min(mu_2 - 2 * sqrt(sigma_2), na.rm = TRUE),
      max(mu_2 + 2 * sqrt(sigma_2), na.rm = TRUE)
    )
    plot.default(range_x, range_y, ..., type = 'n')
  }

  for(chain in seq_len(length(x))){
    ellipses(x[[chain]], ..., add = TRUE, newplot = FALSE, type = type)
  }
}

ellipses.mcmcMapChain <- function(x, ..., add = TRUE, newplot = !add,
                                  alpha = 1, col = 'black', type = 'l'){
  # Set up a window if !add.
  if(!add){
    # Choose a sensible default plotting window from the posterior.
    mu_1 <- unlist(lapply(attr(x, 'mu'), function(m) m[1,]))
    sigma_1 <- unlist(lapply(attr(x, 'Sigma'), function(s) s[1, 1,]))
    range_x <- c(
      min(mu_1 - 2 * sqrt(sigma_1), na.rm = TRUE),
      max(mu_1 + 2 * sqrt(sigma_1), na.rm = TRUE)
    )
    mu_2 <- unlist(lapply(attr(x, 'mu'), function(m) m[2,]))
    sigma_2 <- unlist(lapply(attr(x, 'Sigma'), function(s) s[2, 2,]))
    range_y <- c(
      min(mu_2 - 2 * sqrt(sigma_2), na.rm = TRUE),
      max(mu_2 + 2 * sqrt(sigma_2), na.rm = TRUE)
    )
    plot.default(range_x, range_y, ..., type = 'n')
  }

  h <- attr(x, 'h')
  for(j in seq_len(length(h))){
    for(l in seq_len(length(h[[j]]))){
      if(!is.na(h[[j]][l]) & h[[j]][l] > 0){
        maxh <- max(h[[j]], na.rm = TRUE)
        mixtools::ellipse(mu = attr(x, 'mu')[[j]][,l],
                          sigma = attr(x, 'Sigma')[[j]][,,l],
                          newplot = FALSE, type = type,
                          col = paste0(
                            spatstat::col2hex(col),
                            format(
                              as.hexmode(round(min(255, max(1, alpha * 255 *
                                           h[[j]][l] / maxh)))),
                              width = 2, upper.case = TRUE
                            )
                          ), ...)
      }
    }
  }
}

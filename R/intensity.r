# Posterior intensity methods for poisDPmix.
#intensity <- function(x, ...) UseMethod('intensity')
intensityMap <- function(x, ...) UseMethod('intensityMap')

intensity.mcmcMap <- function(x, location){
  return(do.call(rbind, lapply(x, intensity, location)))
}

intensity.mcmcMapChain <- function(x, location){
  return(as.matrix(x[,'Lambda']) *
           .Call('_poisDPmix_density_mcmcMap', PACKAGE = 'poisDPmix',
                 attr(x, 'mu'), attr(x, 'Sigma'), attr(x, 'h'),
                 location))
}

intensityMap.mcmcMap <- function(x, FUN, q, ...,
                                 grid_x, grid_y, n_x = 41, n_y = 41,
                                 win = NULL, add = FALSE, cl = NULL){
  if(missing(FUN)){
    if(!missing(q)){
      FUN <- 'quantile'
    }else{
      FUN <- 'mean'
    }
  }
  if(is.character(FUN)){
    fun_name <- FUN
  } else {
    fun_name <- deparse(substitute(FUN))
  }
  fun_args <- list(...)
  if(fun_name == 'quantile' & !missing(q)){
    fun_args$probs <- q
  }

  # Choose a sensible default plotting window from the posterior or the window.
  if(add){
    current_win <- par('usr')
  }
  if(missing(grid_x)){
    if(add){
      grid_x <- seq(current_win[1], current_win[2], length.out = n_x)
    }else if(is.owin(win)){
      grid_x <- seq(min(win$xrange), max(win$xrange), length.out = n_x)
    }else{
      mu_1 <- unlist(lapply(do.call(c, lapply(x, attr, 'mu')),
                            function(m) m[1,]))
      sigma_1 <- unlist(lapply(do.call(c, lapply(x, attr, 'Sigma')),
                               function(s) s[1, 1,]))
      grid_x <- seq(
        min(mu_1 - 2 * sqrt(sigma_1), na.rm = TRUE),
        max(mu_1 + 2 * sqrt(sigma_1), na.rm = TRUE),
        length.out = n_x
      )
    }
  }
  if(missing(grid_y)){
    if(add){
      grid_y <- seq(current_win[3], current_win[4], length.out = n_y)
    }else if(is.owin(win)){
      grid_y <- seq(min(win$yrange), max(win$yrange), length.out = n_y)
    }else{
      mu_2 <- unlist(lapply(do.call(c, lapply(x, attr, 'mu')),
                            function(m) m[2,]))
      sigma_2 <- unlist(lapply(do.call(c, lapply(x, attr, 'Sigma')),
                               function(s) s[2, 2,]))
      grid_y <- seq(
        min(mu_2 - 2 * sqrt(sigma_2), na.rm = TRUE),
        max(mu_2 + 2 * sqrt(sigma_2), na.rm = TRUE),
        length.out = n_y
      )
    }
  }

  intensity_grid <- expand.grid(x = grid_x, y = grid_y)
  intensity_grid$inside <- TRUE
  if(is.owin(win)){
    intensity_grid$inside <- inside.owin(intensity_grid$x, intensity_grid$y, win)
  }

  # Function to apply at each grid locations.
  f_intensity <- function(location, mcmcobj, f, f_args){
    if(location[3]){
      return(do.call(
        match.fun(f), c(

          # Loop through posterior draws.
          list(intensity(mcmcobj, location[1:2])),
          f_args

        )
      ))
    }
    return(NA_real_)
  }

  # Loop through grid locations.
  if(inherits(cl, 'cluster') & requireNamespace('parallel', quietly = TRUE)){
    parallel::clusterEvalQ(cl, library(poisDPmix, quietly = TRUE))
    z <- parallel::parApply(cl, intensity_grid, 1, f_intensity,
                      mcmcobj = x, f = FUN, f_args = fun_args)
  }else{
    z <- apply(intensity_grid, 1, f_intensity,
               mcmcobj = x, f = FUN, f_args = fun_args)
  }

  return(structure(matrix(z, nrow = length(grid_x), ncol = length(grid_y)),
                   x = grid_x, y = grid_y, add = add, class = c('postGrid')))
}

# Posterior density methods for poisDPmix.
#density <- function(x, ...) UseMethod('intensity')
densityMap <- function(x, ...) UseMethod('densityMap')

density.mcmcMap <- function(x, location){
  return(do.call(rbind, lapply(x, density, location)))
}

density.mcmcMapChain <- function(x, location){
  return(.Call('_poisDPmix_density_mcmcMap', PACKAGE = 'poisDPmix',
               attr(x, 'mu'), attr(x, 'Sigma'), attr(x, 'h'),
               location))
}

densityMap.mcmcMap <- function(x, FUN, q, ...,
                               grid_x, grid_y, n_x = 41, n_y = 41,
                               win = NULL, add = FALSE, cl = NULL){
  if(missing(FUN)){
    if(!missing(q)){
      FUN <- 'quantile'
    }else{
      FUN <- 'mean'
    }
  }
  if(is.character(FUN)){
    fun_name <- FUN
  } else {
    fun_name <- deparse(substitute(FUN))
  }
  fun_args <- list(...)
  if(fun_name == 'quantile' & !missing(q)){
    fun_args$probs <- q
  }

  # Choose a sensible default plotting window from the posterior or the window.
  if(add){
    current_win <- par('usr')
  }
  if(missing(grid_x)){
    if(add){
      grid_x <- seq(current_win[1], current_win[2], length.out = n_x)
    }else if(is.owin(win)){
      grid_x <- seq(min(win$xrange), max(win$xrange), length.out = n_x)
    }else{
      mu_1 <- unlist(lapply(do.call(c, lapply(x, attr, 'mu')),
                            function(m) m[1,]))
      sigma_1 <- unlist(lapply(do.call(c, lapply(x, attr, 'Sigma')),
                               function(s) s[1, 1,]))
      grid_x <- seq(
        min(mu_1 - 2 * sqrt(sigma_1), na.rm = TRUE),
        max(mu_1 + 2 * sqrt(sigma_1), na.rm = TRUE),
        length.out = n_x
      )
    }
  }
  if(missing(grid_y)){
    if(add){
      grid_y <- seq(current_win[3], current_win[4], length.out = n_y)
    }else if(is.owin(win)){
      grid_y <- seq(min(win$yrange), max(win$yrange), length.out = n_y)
    }else{
      mu_2 <- unlist(lapply(do.call(c, lapply(x, attr, 'mu')),
                            function(m) m[2,]))
      sigma_2 <- unlist(lapply(do.call(c, lapply(x, attr, 'Sigma')),
                               function(s) s[2, 2,]))
      grid_y <- seq(
        min(mu_2 - 2 * sqrt(sigma_2), na.rm = TRUE),
        max(mu_2 + 2 * sqrt(sigma_2), na.rm = TRUE),
        length.out = n_y
      )
    }
  }

  density_grid <- expand.grid(x = grid_x, y = grid_y)
  density_grid$inside <- TRUE
  if(is.owin(win)){
    density_grid$inside <- inside.owin(density_grid$x, density_grid$y, win)
  }

  # Function to apply at each grid locations.
  f_density <- function(location, mcmcobj, f, f_args){
    if(location[3]){
      return(do.call(
        match.fun(f), c(

          # Loop through posterior draws.
          list(density(mcmcobj, location[1:2])),
          f_args

        )
      ))
    }
    return(NA_real_)
  }

  # Loop through grid locations.
  if(inherits(cl, 'cluster') & requireNamespace('parallel', quietly = TRUE)){
    parallel::clusterEvalQ(cl, library(poisDPmix, quietly = TRUE))
    z <- parallel::parApply(cl, density_grid, 1, f_density,
                      mcmcobj = x, f = FUN, f_args = fun_args)
  }else{
    z <- apply(density_grid, 1, f_density,
               mcmcobj = x, f = FUN, f_args = fun_args)
  }

  return(structure(matrix(z, nrow = length(grid_x), ncol = length(grid_y)),
                   x = grid_x, y = grid_y, add = add, class = c('postGrid')))
}

contour.postGrid <- function(z, ..., add = attr(z, 'add')){
  contour.default(x = attr(z, 'x'), y = attr(z, 'y'), z = z, ..., add = add)
}

image.postGrid <- function(z, ...){
  image.default(x = attr(z, 'x'), y = attr(z, 'y'), z = z, ...)
}

as.im.postGrid <- function(z, ...){
  imargs <- list(...)
  imargs$mat <- t(z)
  imargs$xcol <- attr(z, 'x')
  imargs$yrow <- attr(z, 'y')
  return(do.call(im, imargs))
}

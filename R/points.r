# Point plotting methods for poisDPmix.
points.mcmcMap <- function(x, ..., add = TRUE){
  attr_mis <- lapply(lapply(x, attr, 'X_mis'), do.call, what = cbind)
  if(!is.null(attr_mis)){
    X_mis <- do.call(cbind, attr_mis)
    if(add){
      points(X_mis[1,], X_mis[2,], ...)
    } else {
      plot.default(X_mis[1,], X_mis[2,], ...)
    }
  }
}
points.mcmcMapChain <- function(x, ..., add = TRUE){
  attr_mis <- attr(x, 'X_mis')
  if(!is.null(attr_mis)){
    X_mis <- do.call(cbind, attr_mis)
    if(add){
      points(X_mis[1,], X_mis[2,], ...)
    } else {
      plot.default(X_mis[1,], X_mis[2,], ...)
    }
  }
}

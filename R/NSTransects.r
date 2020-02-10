# Generate evenly-spaced transects going north-south.
NSTransects <- function(win, width, spacing,
                        offset = runif(1, 0, spacing)){
  first <- min(vertices(Frame(win))$x) + width / 2 + offset
  last <- max(vertices(Frame(win))$x) - width / 2
  bottom <- min(vertices(Frame(win))$y)
  top <- max(vertices(Frame(win))$y)
  x <- seq(first, last, spacing + width)
  return(structure(data.frame(
    x = sort(rep(x, 2)),
    y = rep(c(bottom - width, top + width, top + width, bottom - width),
            ceiling(length(x)/2))[1:(2*length(x))]
  ), width = width, class = c('cog', 'data.frame')))
}

cog <- function(track, width){
  cls <- class(track)
  if(!('cog' %in% cls)){
    cls <- c('cog', cls)
  }
  return(structure(track, width = width, class = cls))
}

as.owin.cog <- function(track, x_col = 1, y_col = 2){
  n_segments <- nrow(track)
  track_vertices <- cbind(
    x0 = track[-n_segments, x_col],
    y0 = track[-n_segments, y_col],
    x1 = track[-1, x_col],
    y1 = track[-1, y_col]
  )
  track_vertices <- track_vertices[
    !(apply(is.na(track_vertices), 1, any)),
  ]
  return(dilation(psp(track_vertices[,'x0'], track_vertices[,'y0'],
                      track_vertices[,'x1'], track_vertices[,'y1'],
                      dilation(owin(range(track[,x_col], na.rm = TRUE),
                                    range(track[,y_col], na.rm = TRUE)),
                               attr(track, 'width') / 2)),
                  attr(track, 'width') / 2))
}

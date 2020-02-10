# Function to calculate coefficient of variation.
cv <- function(x, ...){return(sd(x, ...) / mean(x, ...))}

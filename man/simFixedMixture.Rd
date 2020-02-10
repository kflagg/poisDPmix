\name{simFixedMixture}
\alias{simFixedMixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a Poisson Process with Fixed Gaussian Intensity
}
\description{
This function simulates a Poisson point processes in \deqn{d}-dimensional
space using a mixture of Gaussian kernels as the intensity of the foreground
process. The mixing weights, centers, and scales of the mixture components are
fixed. A homogeneous background component can also be included.
}
\usage{
simFixedMixture(Lambda, p, mu, Sigma, bound, win, Lambda_0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lambda}{
    the mean number of foreground points.
}
  \item{p}{
    a vector of mixing probabilities for the foreground mixture components.
}
  \item{mu}{
    a matrix of mixing component centers, with a row for each component.
}
  \item{Sigma}{
    a 3-dimensional array of mixing component scale matrices, with
    \code{Sigma[i,,]} being the scale matrix for the \deqn{i}th component.
}
  \item{bound}{
    a list of vectors bounding the region in each dimension. Optional if not
    simulating a background component.
}
  \item{win}{
    a function defining the region in which to simulate. Needs to take an
    \deqn{n} by \deqn{d} matrix of locations as an argument and return a
    logical vector of length \deqn{n}, where the \deqn{i}th entry is
    \code{TRUE} if the \deqn{i}th row of the argument is a location inside the
    region. Optional if the region is a Cartesian product.
}
  \item{Lambda_0}{
    the mean number of background points. Defaults to 0.
}
}
\details{
  TODO: Detail about \code{bound} and \code{win}.
}
\value{
  Returns a matrix with a random number of rows where each row is a random
  location in the region of interest.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  K.A. Flagg
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Lambda, p, mu, Sigma, bound, win, Lambda_0 = 0) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{simulation}% use one of  RShowDoc("KEYWORDS")

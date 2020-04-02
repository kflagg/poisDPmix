# `poisDPmix`: R package to fit Bayesian spatial Cox process models with Dirchlet process mixture intensities

This package implements the methods of Flagg, Hoegh, and Borkowski (2020)
"Modeling Partially Surveyed Point Process Data: Inferring Spatial Point
Intensity of Geomagnetic Anomalies", _JABES_ and supports the example code
accompanying that paper.

The package is not under active development so new features are unlikely to
appear. However, if you find bugs or have trouble installing it, please open
an issue on Github.

# Dependencies

R packages:

- `coda`
- `mixtools`
- `mvtnorm`
- `Rcpp`
- `RcppArmadillo`
- `rgdal`
- `spatstat`

Other libraries:

- `libgdal`

# Installation option 1

If you have GNU-style build tools, simply clone this repo and run:
```
make
```

# Installation option 2

Clone this repo. Open an R terminal in the `poisDPmix/` directory. In R, run:
```
Rcpp::compileAttributes(verbose = TRUE)
```
Then, in a command line terminal, go into the `poisDPmix/` directory and run:
```
R CMD INSTALL .
```

all: R/RcppExports.R src/RcppExports.cpp
	R CMD INSTALL .

clean:
	rm -rf R/RcppExports.R src/RcppExports.cpp src/*.o src/*.so

R/RcppExports.R:
src/RcppExports.cpp: src/density_mcmcMap.cpp src/distributions.cpp src/inRegion.cpp
	Rscript -e "Rcpp::compileAttributes(verbose = TRUE)"


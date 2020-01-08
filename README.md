## SpatPCA Package

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/SpatPCA)](https://CRAN.R-project.org/package=SpatPCA)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/SpatPCA)](https://CRAN.R-project.org/package=SpatPCA)
[![Travis-CI Build Status](https://travis-ci.org/egpivo/SpatPCA.svg?branch=master)](https://travis-ci.org/egpivo/SpatPCA)
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%203%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![Research software impact](http://depsy.org/api/package/cran/SpatPCA/badge.svg)](http://depsy.org/package/r/SpatPCA)



### Description
***SpatPCA*** is an R package that facilitates regularized principal component analysis, 

* seeking the dominant patterns (eigenfunctions), which can be smooth and localized
* computing spatial prediction (Kriging) at new locations
* suitable for either regularly or irregularly spaced data, including 1D, 2D, and 3D
* by the alternating direction method of multipliers (ADMM) algorithm


### Installation
To get the current released version from CRAN:

```r
install.packages("SpatPCA")
```

To get the current development version from GitHub:

```r
devtools::install_github("egpivo/SpatPCA")
```
To compile C++ code with the package [`RcppArmadillo`](https://CRAN.R-project.org/package=RcppArmadillo),

 * Windows users require [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/)
 * Mac users require Xcode Command Line Tools, and install the library gfortran by typing the following lines into terminal

```
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```
More details can be found [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/).

### Usage
```{r example}
library(SpatPCA)

### location
x_1D <- as.matrix(seq(-5, 5, length = 50))

###underlying eigenfunction
Phi_1D <- exp(-x_1D^2) / norm(exp(-x_1D^2), "F")

### Realizations on x_1D
set.seed(1234)
Y_1D <- rnorm(n = 100, sd = 3) %*% t(Phi_1D) + matrix(rnorm(n = 100 * 50), 100, 50)

### main function: spatpca()
cv_1D <- spatpca(x = x_1D, Y = Y_1D)

### Plot the estimate
plot(x_1D, cv_1D$eigenfn[,1], type = 'l', main = '1st eigenfunction')
lines(x_1D, svd(Y_1D)$v[,1], col='red')
legend('topleft', c('SpatPCA', 'PCA'), lty = 1:1, col = 1:2)
```
### Author
[Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b "Wen-Ting Wang") and [Hsin-Cheng Huang](http://www.stat.sinica.edu.tw/hchuang/ "Hsin-Cheng Huang")
 
### Maintainer
[Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b "Wen-Ting Wang")

### Reference
Wang, W.-T. and Huang, H.-C. (2017). [Regularized principal component analysis for spatial data](https://arxiv.org/pdf/1501.03221v3.pdf, "Regularized principal component analysis for spatial data"). *Journal of Computational and Graphical Statistics*, **26**, 14-25.
 
### License
  GPL-3

## SpatPCA Package

[![R build status](https://github.com/egpivo/SpatPCA/workflows/R-CMD-check/badge.svg)](https://github.com/egpivo/SpatPCA/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/egpivo/SpatPCA/master.svg)](https://codecov.io/github/egpivo/SpatpCA?branch=master)

### Description
***SpatPCA*** is an R package that facilitates regularized principal component analysis, 

* seeking the dominant patterns (eigenfunctions), which can be smooth and localized
* computing spatial prediction (Kriging) at new locations
* suitable for either regularly or irregularly spaced data, including 1D, 2D, and 3D
* by the alternating direction method of multipliers (ADMM) algorithm


### Installation
To get the current development version from GitHub:
   ```r
   devtools::install_github("egpivo/SpatPCA")
   ```
To compile C++ code with the package [`RcppArmadillo`](https://CRAN.R-project.org/package=RcppArmadillo),

 * Windows users require [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/)
 * Mac users require Xcode Command Line Tools, and install the library gfortran by typing the following lines into terminal
    ```
    brew update
    brew install gcc
    ```
    The detailed solution is describd [here](https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/), or download and install the library [`gfortran`](https://github.com/fxcoudert/gfortran-for-macOS/releases) to solve the error `ld: library not found for -lgfortran`.

### Usage
```r
library(SpatPCA)
spatpca(position, realizations)
```
- Input: realizations with the corresponding position
- Output: return the most dominant eigenfunctions automatically.
- More details can be referred to [Demo](https://egpivo.github.io/SpatPCA/articles/)

### Author
[Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b) and [Hsin-Cheng Huang](https://sites.stat.sinica.edu.tw/hchuang/)
 
### Maintainer
[Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b)

### Reference
Wang, W.-T. and Huang, H.-C. (2017). [Regularized principal component analysis for spatial data](https://arxiv.org/pdf/1501.03221v3.pdf, "Regularized principal component analysis for spatial data"). *Journal of Computational and Graphical Statistics*, **26**, 14-25.
 
### License
  GPL-3

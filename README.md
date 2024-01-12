# SpatPCA: Regularized Principal Component Analysis for Spatial Data

[![License](https://eddelbuettel.github.io/badges/GPL2+.svg)](https://www.gnu.org/licenses/gpl-2.0.html)
[![R build status](https://github.com/egpivo/SpatPCA/workflows/R-CMD-check/badge.svg)](https://github.com/egpivo/SpatPCA/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/egpivo/SpatPCA/master.svg)](https://app.codecov.io/github/egpivo/SpatpCA?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/SpatPCA)](https://CRAN.R-project.org/package=SpatPCA)
[![Downloads (monthly)](https://cranlogs.r-pkg.org/badges/SpatPCA?color=brightgreen)](https://www.r-pkg.org/pkg/SpatPCA)
[![Downloads (total)](https://cranlogs.r-pkg.org/badges/grand-total/SpatPCA?color=brightgreen)](https://www.r-pkg.org/pkg/SpatPCA)
[![JCGS](https://img.shields.io/badge/JCGS-10.18637%2F10618600.2016.1157483-brightgreen)](https://doi.org/10.1080/10618600.2016.1157483)


## Description
**SpatPCA** is an R package designed for efficient regularized principal component analysis, providing the following features:

- Identify dominant spatial patterns (eigenfunctions) with both smooth and localized characteristics.
- Conduct spatial prediction (Kriging) at new locations.
- Adapt to regularly or irregularly spaced data, spanning 1D, 2D, and 3D datasets.
- Implement using the alternating direction method of multipliers (ADMM) algorithm.


## Installation
You can install **SpatPCA** using either of the following methods:

### Install from CRAN

```r
install.packages("SpatPCA")
```
### Install the Development Version from GitHub
```r
remotes::install_github("egpivo/SpatPCA")
```
### Compilation Requirements
To compile C++ code with the required [`RcppArmadillo`](https://CRAN.R-project.org/package=RcppArmadillo) and [`RcppParallel`](https://CRAN.R-project.org/package=RcppParallel)  packages, follow these instructions based on your operating system:


#### For Windows users
Install [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/)

#### For Mac users
1. Install Xcode Command Line Tools
2. install the `gfortran` library. You can achieve this by running the following commands in the terminal:
  ```bash
  brew update
  brew install gcc
  ```

  For a detailed solution, refer to [this link](https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/), or download and install the library [`gfortran`](https://github.com/fxcoudert/gfortran-for-macOS/releases) to resolve the error `ld: library not found for -lgfortran`.

## Usage
To use **SpatPCA**, first load the package:

```r
library(SpatPCA)
```

Then, apply the `spatpca` function with the following syntax:
```r
spatpca(position, realizations)
```
   - Input: Realizations with the corresponding positions.
   - Output: Return the most dominant eigenfunctions automatically.

For more details, refer to the [Demo](https://egpivo.github.io/SpatPCA/articles/).

### Authors
- [Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b) ([GitHub](https://www.github.com/egpivo))
- [Hsin-Cheng Huang](https://sites.stat.sinica.edu.tw/hchuang/)
 
### Maintainer
[Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b) ([GitHub](https://www.github.com/egpivo))

### Reference
Wang, W.-T. and Huang, H.-C. (2017). [Regularized principal component analysis for spatial data](https://arxiv.org/pdf/1501.03221v3.pdf, "Regularized principal component analysis for spatial data"). *Journal of Computational and Graphical Statistics*, **26**, 14-25.
 
## License
GPL (>= 2)

## Citation
- To cite package ‘SpatPCA’ in publications use:
```
  Wang W, Huang H (2023). SpatPCA: Regularized Principal Component Analysis for
  Spatial Data_. R package version 1.3.5,
  <https://CRAN.R-project.org/package=SpatPCA>.
```

- A BibTeX entry for LaTeX users is
```
  @Manual{,
    title = {SpatPCA: Regularized Principal Component Analysis for Spatial Data},
    author = {Wen-Ting Wang and Hsin-Cheng Huang},
    year = {2023},
    note = {R package version 1.3.5},
    url = {https://CRAN.R-project.org/package=SpatPCA},
  }
```

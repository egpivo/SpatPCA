# SpatPCA Package
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/SpatPCA)](https://CRAN.R-project.org/package=SpatPCA)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/SpatPCA)](https://CRAN.R-project.org/package=SpatPCA)
[![R build status](https://github.com/egpivo/SpatPCA/workflows/R-CMD-check/badge.svg)](https://github.com/egpivo/SpatPCA/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/egpivo/SpatPCA/master.svg)](https://app.codecov.io/github/egpivo/SpatpCA?branch=master)

## Description
**SpatPCA** is an R package designed for efficient regularized principal component analysis, providing the following features:

- Identification of dominant spatial patterns (eigenfunctions) with both smooth and localized characteristics.
- Spatial prediction (Kriging) at new locations.
- Adaptability for regularly or irregularly spaced data, spanning 1D, 2D, and 3D datasets.
- Implementation using the alternating direction method of multipliers (ADMM) algorithm.


## Installation
There are two ways to install the package:

1. **Install from CRAN**:
   ```r
   install.packages("SpatPCA")
   ```
2. **Install the current development version from GitHub**:
   ```r
   remotes::install_github("egpivo/SpatPCA")
   ```

For compiling C++ code with the required [`RcppArmadillo`](https://CRAN.R-project.org/package=RcppArmadillo) and [`RcppParallel`](https://CRAN.R-project.org/package=RcppParallel) packages, follow these instructions:

* Windows users: Install [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/)
* Mac users: Install Xcode Command Line Tools, and install the `gfortran` library. You can achieve this by running the following commands in the terminal:
  ```bash
  brew update
  brew install gcc
  ```

For a detailed solution, refer to [this link](https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/), or download and install the library [`gfortran`](https://github.com/fxcoudert/gfortran-for-macOS/releases) to resolve the error `ld: library not found for -lgfortran`.

## Usage
```r
library(SpatPCA)
spatpca(position, realizations)
```

- Input: Realizations with the corresponding positions.
- Output: Return the most dominant eigenfunctions automatically.
- For more details, refer to the [Demo](https://egpivo.github.io/SpatPCA/articles/).

### Author
- [Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b)
- [Hsin-Cheng Huang](https://sites.stat.sinica.edu.tw/hchuang/)
 
### Maintainer
[Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b)

### Reference
Wang, W.-T. and Huang, H.-C. (2017). [Regularized principal component analysis for spatial data](https://arxiv.org/pdf/1501.03221v3.pdf, "Regularized principal component analysis for spatial data"). *Journal of Computational and Graphical Statistics*, **26**, 14-25.
 
## License
GPL-3

## Citation
- To cite package ‘SpatPCA’ in publications use:
```
  Wang W, Huang H (2023). _SpatPCA: Regularized Principal Component Analysis for
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

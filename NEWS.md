## SpatPCA v1.3.6 (Release date: 2024-01-12)
#### Overview 
This release focuses on 
  1. fixing the warnings from the earlier Rcpp pacakge.
  2. adding LICENSE
  3. refining README 


## SpatPCA v1.3.4 (Release date: 2023-11-11)
#### Overview 
This release focuses on refining the package to comply with CRAN's regulatory standards.

#### MAINTENANCE
- Achieve code coverage exceeding 100%
- Upgrade RcppParallel to version 5.1.7
- Correct typos for improved clarity


## SpatPCA v1.3.3.7 (Release date: 2023-01-28)
#### Overview 
In this release, we focus on package maintenance.

#### MAINTENANCE
- Enhance code coverage over `99%`
- Update deprecated R settings and GitHub workflow files
- Fix memory leak errors caused by `tbb` backend in RcppParallel when testing the package on the platform `Debian Linux, R-devel, GCC ASAN/UBSAN`
- Enhance the readability of the resultant plot by adding math expression
- Fix grammar errors


## SpatPCA v1.3.3.0 (Release date: 2021-01-31)
#### Overview 
In this release, we take care of the perspective of software quality by refactoring code for better readability and adding unit-tests. Accordingly, we add multiple features by separating implicit functions from `spatpca`.

#### NEW FEATURES
1. `thinPlateSplineMatrix()` for producing a thin-plane spline matrix
2. `predictEigenfunction()`for estimating K dominant patterns on new sites
3. `predict()` for predicting new target variable `Y` on new sites   
4. `plot()` for `spatpca` objects for plotting M-fold CV results

#### MAINTENANCE
- Add unit tests with code coverage `87%`




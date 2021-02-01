// includes from the plugin
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;
using namespace arma;

struct thinPlateSpline: public RcppParallel::Worker {
  const mat& P;
  mat& L;  
  int p;
  int d;
  thinPlateSpline(const mat &P, mat& L, int p, int d): P(P), L(L), p(p), d(d) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++) {
      for(unsigned j = 0; j < p; ++j) {
        if(j > i) {
          if(d == 1) {
            double r = sqrt(pow(P(i, 0) - P(j, 0), 2));
            L(i, j) = pow(r, 3) / 12;
          }
          else if(d == 2) {
            double r = sqrt(pow(P(i, 0) - P(j, 0), 2) + (pow(P(i, 1) - P(j, 1), 2)));
            L(i, j) = r * r * log(r) / (8.0 * datum::pi);
          }
          else if(d == 3) {
            double r = sqrt(pow(P(i, 0) - P(j, 0), 2) +
                            pow(P(i, 1) - P(j, 1), 2) +
                            pow(P(i, 2) - P(j, 2), 2));
            L(i, j) = -r / (8.0 * datum::pi);
          }
        }
      }  
      
      L(i, p) = 1;
      for(unsigned k = 0; k < d; ++k)
        L(i, p + k + 1) = P(i, k);
    }
  }
};

//' @title Thin-plane spline matrix
//' 
//' @description Produce a thin-plane spline matrix based on a given location matrix
//' 
//' @param location A location matrix
//' @return A thin-plane spline matrix
//' @examples
//' pesudo_sequence <- seq(-5, 5, length = 5)
//' two_dim_location <- as.matrix(expand.grid(x = pesudo_sequence, y = pesudo_sequence))
//' thin_plate_matrix <- thinPlateSplineMatrix(two_dim_location)
// [[Rcpp::export]]
arma::mat thinPlateSplineMatrix(const arma::mat location) {
  int p = location.n_rows, d = location.n_cols;
  int total_size = p + d;
  mat L, Lp, Ip;
  
  L.zeros(total_size + 1, total_size + 1);
  Ip.eye(total_size + 1, total_size + 1);
  thinPlateSpline thin_plate_spline(location, L, p, d);
  parallelFor(0, p, thin_plate_spline);
  L = symmatu(L);
  Lp = inv(L + 1e-8 * Ip);   
  Lp.shed_cols(p, total_size);
  Lp.shed_rows(p, total_size);
  L.shed_cols(p, total_size);
  L.shed_rows(p, total_size);
  mat result = Lp.t() * (L * Lp);
  return(result);
}

//' @title Interpolated Eigen-function
//' 
//' @description Produce Eigen-function values based on new locations
//' 
//' @keywords internal
//' @param new_location A location matrix
//' @param original_location A location matrix
//' @param Phi An eigenvector matrix
//' @return A predictive estimate matrix
//' @examples
//' pesudo_sequence <- seq(-5, 5, length = 2)
//' original_location <- as.matrix(expand.grid(x = pesudo_sequence, y = pesudo_sequence))
//' new_location <- matrix(c(0.1, 0.2), nrow = 1, ncol = 2)
//' Phi <- matrix(c(1, 0, 0, 0), nrow = 4, ncol = 1)
//' thin_plate_matrix <- eigenFunction(new_location, original_location, Phi)
// [[Rcpp::export]]
arma::mat eigenFunction(const arma::mat new_location, const arma::mat original_location, const arma::mat Phi) {
  mat L;
  int p = original_location.n_rows, d = original_location.n_cols, K = Phi.n_cols;
  int total_size = p + d;
  L.zeros(total_size + 1, total_size + 1);
  thinPlateSpline thin_plate_spline(original_location, L, p, d);
  parallelFor(0, p, thin_plate_spline);
  L = L + L.t();
  
  mat Phi_star, para(total_size + 1, K);
  Phi_star.zeros(total_size + 1, K);
  Phi_star.rows(0, p - 1) = Phi;
  para = solve(L, Phi_star);
  int pnew = new_location.n_rows;
  mat eigen_fn(pnew, K);
  double psum, r;

  for(unsigned newi = 0; newi < pnew ; newi++) {
    for(unsigned i = 0; i < K; i++) {
      psum = 0;
      for(unsigned j = 0; j < p; j++) {
        if(d == 1) {  
          r = norm(new_location.row(newi) - original_location.row(j), "f");
          if(r != 0)
            psum += para(j, i) * pow(r, 3) / 12;
        }
        else if(d == 2) {
          r = sqrt(pow(new_location(newi, 0) - original_location(j, 0), 2) + 
            (pow(new_location(newi, 1) - original_location(j, 1), 2)));
          if(r != 0)
            psum += para(j, i) * r * r * log(r) / (8.0 * datum::pi);
        }
        else if(d == 3) {
          double r = sqrt(pow(new_location(newi, 0) - original_location(j, 0), 2) +
                          pow(new_location(newi, 1) - original_location(j, 1), 2) +
                          pow(new_location(newi, 2) - original_location(j, 2), 2));
          if(r != 0)
            psum -= para(j, i) * r / (8.0 * datum::pi);
        }
      }
      if(d == 1)
        eigen_fn(newi, i) = psum + para(p + 1, i) * new_location(newi, 0) + para(p, i);
      else if(d == 2)
        eigen_fn(newi, i) = psum + para(p + 1, i) * new_location(newi, 0) +
          para(p + 2, i) * new_location(newi, 1) + para(p, i);
      else if(d == 3)
        eigen_fn(newi, i) = psum + para(p + 1, i) * new_location(newi, 0) +
          para(p + 2, i) * new_location(newi, 1) + para(p + 3, i) * new_location(newi, 2) + para(p, i); 
    }
  }
  
  return(eigen_fn);
}

// user includes
void spatpcaCore2(
    const mat gram_matrix_Y,
    mat& Phi,
    mat& C,
    mat& Lambda2,
    const mat Omega,
    const double tau1,
    const double rho,
    const int maxit,
    const double tol) {
  int p = Phi.n_rows, K = Phi.n_cols, iter = 0;
  mat Ip, Sigtau1, temp, tempinv, U, V, diff, Cold = C, Lambda2old = Lambda2;
  vec error(2), S;
  Ip.eye(p,p);
  Sigtau1 = tau1 * Omega - gram_matrix_Y;
  tempinv = inv_sympd(2 * Sigtau1 + rho * Ip);

  for(iter = 0; iter < maxit; iter++) {
    Phi = tempinv * ((rho * Cold) - Lambda2old);
    temp = Phi + (Lambda2old / rho);
    svd_econ(U, S, V, temp);
    C = U.cols(0, V.n_cols - 1) * V.t();
    diff = Phi - C;
    Lambda2 = Lambda2old + rho * diff;
    
    error[0] = norm(diff, "fro") / sqrt(p * K / 1.0);
    error[1] = norm(C - Cold, "fro") / sqrt(p * K / 1.0);
    
    if(max(error) <= tol)
      break;
    Cold = C;
    Lambda2old = Lambda2;
  }
  
  iter++;
}

mat spatpcaCore2p(
    const mat gram_matrix_Y,
    mat& C,
    mat& Lambda2,
    const mat Omega,
    const double tau1,
    const double rho,
    const int maxit,
    const double tol) {
  int p = C.n_rows, K = C.n_cols, iter = 0;
  mat Ip, Sigtau1, temp, tempinv, U, V, diff, Phi, Cold = C, Lambda2old = Lambda2;
  vec error(2), S;

  Ip.eye(p, p);
  Sigtau1 = tau1 * Omega - gram_matrix_Y;
  
  tempinv = inv_sympd(2 * Sigtau1 + rho * Ip);
  for(iter = 0; iter < maxit; iter++) {
    Phi = tempinv * ((rho * Cold) - Lambda2old);
    temp = Phi + (Lambda2old / rho);
    svd_econ(U, S, V, temp);
    C = U.cols(0, V.n_cols - 1) * V.t();
    diff = Phi - C;
    Lambda2 = Lambda2old + rho * diff;
    
    error[0] = norm(diff, "fro") / sqrt(p * K / 1.0);
    error[1] = norm(C - Cold, "fro") / sqrt(p * K / 1.0);
    
    if(max(error) <= tol)
      break;
    Cold = C;
    Lambda2old = Lambda2;
  }
  iter++;
  return(Phi);
}

void spatpcaCore3(
    const mat tempinv,
    mat& Phi,
    mat& R,
    mat& C,
    mat& Lambda1,
    mat& Lambda2,
    const double tau2,
    double rho,
    const int maxit,
    const double tol) {
  int p = Phi.n_rows, K = Phi.n_cols, iter = 0;
  mat temp, scaled_tau2, zero, one, U, V, difference_Phi_R, difference_Phi_C;
  mat Phi_old = Phi, Rold = R, Cold = C, Lambda1old = Lambda1, Lambda2old = Lambda2;
  vec error(4), S;
  
  zero.zeros(p, K);
  one.ones(p, K);
  scaled_tau2 = tau2 * one / rho;
  
  
  for(iter = 0; iter < maxit; iter++) {
    Phi = 0.5 * tempinv * (rho * (Rold + Cold) - (Lambda1old + Lambda2old));
    R = sign(((Lambda1old / rho) + Phi)) % max(zero, abs(((Lambda1old / rho) + Phi)) - scaled_tau2);
    temp = Phi + Lambda2old / rho;
    svd_econ(U, S, V, temp);
    C = U.cols(0, V.n_cols - 1) * V.t();
    difference_Phi_R = Phi - R;
    difference_Phi_C = Phi - C;    
    Lambda1 = Lambda1old + rho * difference_Phi_R;
    Lambda2 = Lambda2old + rho * difference_Phi_C;
    
    error[1] = norm(difference_Phi_R, "fro") / sqrt(p / 1.0);
    error[2] = norm(R - Rold, "fro") / sqrt(p / 1.0);
    error[3] = norm(difference_Phi_C, "fro") / sqrt(p / 1.0);
    error[4] = norm(C - Cold, "fro") / sqrt(p / 1.0);
    
    if(max(error) <= tol)
      break;
    Phi_old = Phi;
    Rold = R;
    Cold = C;
    Lambda1old = Lambda1;
    Lambda2old = Lambda2;
  }
  
  iter++;
}

struct spatpcaCVPhi: public RcppParallel::Worker {
  const mat& Y;
  int K;
  const mat& Omega;
  const vec& tau1;
  const vec& nk;
  int maxit;
  double tol;
  mat& output;
  cube& gram_matrix_Y_train;
  cube& Phi_cv;
  cube& Lambd2_cv;
  mat& rho;
  spatpcaCVPhi(
    const mat& Y,
    int K,
    const mat& Omega,
    const vec& tau1,
    const vec& nk,
    int maxit,
    double tol,
    mat& output,
    cube& gram_matrix_Y_train,
    cube& Phi_cv,
    cube& Lambd2_cv,
    mat& rho):
      Y(Y),
      K(K),
      Omega(Omega),
      tau1(tau1),
      nk(nk),
      maxit(maxit),
      tol(tol),
      output(output),
      gram_matrix_Y_train(gram_matrix_Y_train),
      Phi_cv(Phi_cv),
      Lambd2_cv(Lambd2_cv),
      rho(rho)
    {}

  void operator()(std::size_t begin, std::size_t end) {
    mat Ip;
    Ip.eye(Y.n_cols, Y.n_cols);
    for(std::size_t k = begin; k < end; k++) {
      mat svd_U, Phi_old, Phi,C, Lambda2;
      vec singular_value;
      mat Y_train = Y.rows(find(nk != (k + 1)));
      mat Y_validation = Y.rows(find(nk == (k + 1)));
      
      svd_econ(svd_U, singular_value, Phi_old, Y_train, "right");  
      rho(k, 0) = 10 * pow(singular_value[0], 2.0);
      
      Phi_cv.slice(k) = Phi_old.cols(0, K - 1);
      Phi = C = Phi_cv.slice(k);
      Lambda2 = Lambd2_cv.slice(k) = Phi * (diagmat(rho(k, 0) - 1 / (rho(k, 0) - 2 * pow(singular_value.subvec(0, K - 1), 2))));
      output(k, 0) = pow(norm(Y_validation * (Ip - (Phi_cv.slice(k)) * (Phi_cv.slice(k)).t()), "fro"), 2.0); 
      gram_matrix_Y_train.slice(k) = Y_train.t() * Y_train;
      for(uword i = 1; i < tau1.n_elem; i++) {
        spatpcaCore2(gram_matrix_Y_train.slice(k), Phi, C, Lambda2, Omega, tau1[i], rho(k, 0), maxit, tol);
        output(k, i) = pow(norm(Y_validation * (Ip - Phi * Phi.t()), "fro"), 2.0); 
      }
    }
  }
};

struct spatpcaCVPhi2: public RcppParallel::Worker {
  const mat& Y;
  const cube& gram_matrix_Y_train;
  cube& Phi_cv;
  cube& Lambd2_cv;
  const mat& rho;
  int K;
  double tau1;
  const mat& Omega;
  const vec& tau2;
  const vec& nk;
  int maxit;
  double tol;
  mat& output;
  cube& tempinv;
  
  spatpcaCVPhi2(
    const mat& Y,
    const cube& gram_matrix_Y_train,
    cube& Phi_cv,
    cube& Lambd2_cv,
    const mat& rho,
    int K,
    double tau1,
    const mat& Omega,
    const vec& tau2,
    const vec& nk,
    int maxit, 
    double tol,
    mat& output,
    cube& tempinv):
      Y(Y),
      gram_matrix_Y_train(gram_matrix_Y_train),
      Phi_cv(Phi_cv),
      Lambd2_cv(Lambd2_cv),
      rho(rho),
      K(K),
      tau1(tau1),
      Omega(Omega),
      tau2(tau2),
      nk(nk),
      maxit(maxit),
      tol(tol),
      output(output),
      tempinv(tempinv)
    {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t k = begin; k < end; k++) {
      mat Phi, R, C, Lambda1,Ip, Lambda2;
      mat Y_validation = Y.rows(find(nk == (k + 1)));
      Ip.eye(Y.n_cols, Y.n_cols);
      Phi = C = Phi_cv.slice(k);
      Lambda1 = 0 * Phi_cv.slice(k);
      Lambda2 = Lambd2_cv.slice(k);
      if(tau1 != 0)
        spatpcaCore2(gram_matrix_Y_train.slice(k), Phi,C, Lambda2, Omega, tau1, rho(k, 0), maxit, tol);
      R = Phi;
      Phi_cv.slice(k) = Phi;
      Lambd2_cv.slice(k) = Lambda2;
      tempinv.slice(k) = inv_sympd((tau1 * Omega) - gram_matrix_Y_train.slice(k) + (rho(k, 0) * Ip));
      for(uword i = 0; i < tau2.n_elem; i++) {
        spatpcaCore3(tempinv.slice(k), Phi, R, C, Lambda1, Lambda2, tau2[i], rho(k, 0), maxit, tol);
        output(k, i) = pow(norm(Y_validation * (Ip - Phi * Phi.t()), "fro"), 2.0); 
      }
    }
  }
};

struct spatpcaCVPhi3: public RcppParallel::Worker {
  const mat& Y;
  const cube& gram_matrix_Y_train;
  const cube& Phi_cv;
  const cube& Lambd2_cv;
  const mat& rho;
  const cube& tempinv;
  const uword index;
  int K;
  const mat& Omega;
  const double tau1; 
  const vec& tau2;
  const vec& gamma;
  const vec& nk;
  int maxit;
  double tol;
  mat& output;
  
  spatpcaCVPhi3(
    const mat& Y,
    const cube& gram_matrix_Y_train,
    const cube& Phi_cv,
    const cube& Lambd2_cv,
    const mat& rho,
    const cube& tempinv,
    const uword index,
    int K,
    const mat& Omega,
    const double tau1,
    const vec& tau2,
    const vec& gamma,
    const vec& nk,
    int maxit,
    double tol,
    mat& output):
      Y(Y),
      gram_matrix_Y_train(gram_matrix_Y_train),
      Phi_cv(Phi_cv),
      Lambd2_cv(Lambd2_cv),
      rho(rho),
      tempinv(tempinv),
      index(index),
      K(K),
      Omega(Omega),
      tau1(tau1),
      tau2(tau2),
      gamma(gamma),
      nk(nk),
      maxit(maxit),
      tol(tol),
      output(output)
    {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t k = begin; k < end; k++) {
      int p = Y.n_cols, tempL, K;
      mat  Phi, C, Ip, Lambda2;
      mat Vc, Vc2, covariance_mtraix_train, covariance_mtraix_validation, estimated_covariance, eigenvalue;
      mat Y_validation = Y.rows(find(nk == (k + 1)));
      vec Sc, Sc2, Sct, Sctz, cv;
      double total_variance, error, temp, tempSc, tempSc2;

      Ip.eye(Y.n_cols, Y.n_cols);
      Phi = C = Phi_cv.slice(k);
      K = Phi.n_cols; 
      Lambda2 = Lambd2_cv.slice(k);
      
      if(max(tau2) != 0 || max(tau1) != 0) {
        if(tau2.n_elem == 1) {
          spatpcaCore2(gram_matrix_Y_train.slice(k), Phi, C, Lambda2, Omega, tau1, rho(k, 0), maxit, tol);
        }
        else {
          mat R = Phi;
          mat Lambda1 = 0 * Phi_cv.slice(k);
          for(uword i = 0; i <= index; i++) {
            spatpcaCore3(tempinv.slice(k), Phi, R, C, Lambda1, Lambda2, tau2[i], rho(k, 0), maxit, tol);
          }
        }
      }
      else {
        Phi = Phi_cv.slice(k);
      }
 
      covariance_mtraix_train = gram_matrix_Y_train.slice(k) / (Y.n_rows - Y_validation.n_rows);
      covariance_mtraix_validation = trans(Y_validation) * Y_validation / Y_validation.n_rows;
      total_variance = trace(covariance_mtraix_train);
      eig_sym(Sc, Vc, trans(Phi) * covariance_mtraix_train * Phi);
      tempSc2 = accu(Sc);
      Sc2 = sort(Sc,"descend");
      Vc2 = Vc.cols(sort_index(Sc,"descend"));
      Sct.ones(Sc2.n_elem);
      Sctz.zeros(Sc2.n_elem);

      for(unsigned int gj = 0; gj < gamma.n_elem; gj++) {
        tempSc = tempSc2;
        tempL = K;
        if(Sc2[0] > gamma[gj]) {
          error = (total_variance - tempSc + K * gamma[gj]) / (p - tempL);
          temp = Sc2[tempL - 1];
          while(temp - gamma[gj] < error) {
            if(tempL == 1) {
              error = (total_variance - Sc2[0] + gamma[gj]) / (p - 1);
              break;
            }
            tempSc -= Sc2[tempL - 1];
            tempL --;
            error = (total_variance - tempSc + tempL * gamma[gj]) / (p - tempL);
            temp = Sc2[tempL - 1];
          }
          if(Sc2[0] - gamma[gj] < error)
            error = total_variance / p;
        }
        else
          error = total_variance / p;
        eigenvalue = max(Sc2 - (error + gamma[gj]) * Sct, Sctz);
        estimated_covariance =  Phi * Vc2 * diagmat(eigenvalue) * trans(Vc2) * trans(Phi);
        output(k, gj) = pow(norm(covariance_mtraix_validation - estimated_covariance - error * Ip, "fro"), 2.0);
      }
    }
  }
};

//' Internal function: M-fold Cross-validation 
//' @keywords internal
//' @param sxyr A location matrix
//' @param Yr A data matrix
//' @param M The number of folds for CV
//' @param K The number of estimated eigen-functions
//' @param tau1r A range of tau1
//' @param tau2r A range of tau2
//' @param gammar A range of gamma
//' @param nkr A vector of fold numbers
//' @param maxit A maximum number of iteration
//' @param tol A tolerance rate
//' @param l2r A given tau2
//' @return A list of selected parameters
// [[Rcpp::export]]
List spatpcaCV(
    NumericMatrix sxyr,
    NumericMatrix Yr,
    int M,
    int K,
    NumericVector tau1r,
    NumericVector tau2r,
    NumericVector gammar,
    NumericVector nkr,
    int maxit,
    double tol,
    NumericVector l2r) {
  int n = Yr.nrow(), p = Yr.ncol(), d = sxyr.ncol();
  mat Y(Yr.begin(), n, p, false), sxy(sxyr.begin(), p, d, false);
  colvec tau1(tau1r.begin(), tau1r.size(),false);
  colvec tau2(tau2r.begin(), tau2r.size(), false);
  colvec gamma(gammar.begin(), gammar.size(), false);
  colvec nk(nkr.begin(), nkr.size(), false);
  colvec l2(l2r.begin(), l2r.size(), false);
  mat cv(M, tau1.n_elem), cv3(M, gamma.n_elem), cv_score_tau1, cv_score_tau2, cv_score_gamma, Omega, svd_U, svd_V;
  double selected_tau1, selected_tau2 = 0, selected_gamma;
  mat gram_matrix_Y = Y.t() * Y;
  vec singular_value;
  svd_econ(svd_U, singular_value, svd_V, Y, "right");
  double estimated_rho = 10 * pow(singular_value[0], 2.0);
  mat estimated_Phi = svd_V.cols(0, K - 1); 
  mat estimated_C = estimated_Phi;
  mat estimated_Lambda2 = estimated_Phi * (diagmat(estimated_rho - 1 / (estimated_rho - 2 * pow(singular_value.subvec(0, K - 1), 2))));
  mat rho_cv(M, 1);
  cube gram_matrix_Y_train(p, p, M), Phi_cv(p, K, M), Lambd2_cv(p, K, M), tempinv_cv(p, p, M);
  uword index1, index2 = 0, index3;
  mat Ip;
  Ip.eye(Y.n_cols, Y.n_cols);
  if(max(tau1) != 0 || max(tau2) != 0) {
    Omega = thinPlateSplineMatrix(sxy) + 1e-8 * Ip;
  }
  else {
    if(gamma.n_elem > 1) {
      Omega = Ip;
      mat Y_train, Y_validation; 
      mat svd_U, Phi_oldg;
      vec singular_value;
      for(unsigned k = 0; k < M; ++k) {
        Y_train = Y.rows(find(nk != (k + 1)));
        Y_validation = Y.rows(find(nk == (k + 1)));
        svd_econ(svd_U, singular_value, Phi_oldg, Y_train, "right");
        Phi_cv.slice(k) = Phi_oldg.cols(0, K - 1);
        gram_matrix_Y_train.slice(k) = Y_train.t() * Y_train;
      }
    }
  }

  if(tau1.n_elem > 1) {  
    spatpcaCVPhi spatpcaCVPhi(Y, K, Omega, tau1, nk, maxit, tol, cv, gram_matrix_Y_train, Phi_cv, Lambd2_cv, rho_cv);
    RcppParallel::parallelFor(0, M, spatpcaCVPhi);
    (sum(cv, 0)).min(index1);
    selected_tau1 = tau1[index1];  
    if(index1 > 0)
      estimated_Phi = spatpcaCore2p(gram_matrix_Y, estimated_C, estimated_Lambda2, Omega, selected_tau1, estimated_rho, maxit, tol);
    cv_score_tau1 = sum(cv, 0) / M;
  }
  else {
    selected_tau1 = max(tau1);
    if(selected_tau1 != 0 && max(tau2) == 0) {
      estimated_Phi = spatpcaCore2p(gram_matrix_Y, estimated_C, estimated_Lambda2, Omega, selected_tau1, estimated_rho, maxit, tol);
      mat Phigg, Cgg, Lambda2gg;
      mat Y_validation, Y_train;
      mat svd_U, Phi_oldg;
      vec singular_value;

      for(unsigned k = 0; k < M; ++k) {
        Y_train = Y.rows(find(nk != (k + 1)));
        gram_matrix_Y_train.slice(k) = Y_train.t() * Y_train;
        svd_econ(svd_U, singular_value, Phi_oldg, Y_train, "right");
        Phigg = Cgg = Phi_oldg.cols(0, K - 1);
        rho_cv(k, 0) = 10 * pow(singular_value[0], 2.0);
        Lambda2gg = Phigg * (diagmat(rho_cv(k, 0) - 1 / (rho_cv(k, 0) - 2 * pow(singular_value.subvec(0, K - 1), 2))));
        spatpcaCore2(gram_matrix_Y_train.slice(k), Phigg,Cgg, Lambda2gg, Omega, selected_tau1, rho_cv(k, 0), maxit, tol);
        Phi_cv.slice(k) = Phigg;
        Lambd2_cv.slice(k) = Lambda2gg;
        tempinv_cv.slice(k) = inv_sympd((selected_tau1 * Omega) - gram_matrix_Y_train.slice(k) + (rho_cv(k, 0) * Ip));    
      }
    }
    else if(selected_tau1 == 0 && max(tau2) == 0) {
      mat svd_U2, Phi_oldc;
      vec singular_value2;
      svd_econ(svd_U2, singular_value2, Phi_oldc, Y, "right");
      estimated_Phi = Phi_oldc.cols(0,K - 1); 
    }
    else {
      mat Phigg, Cgg, Lambda2gg;
      mat Y_validation, Y_train;
      mat svd_U, Phi_oldg;
      vec singular_value;

      for(unsigned k = 0; k < M; ++k) {
        Y_train = Y.rows(find(nk != (k + 1)));
        gram_matrix_Y_train.slice(k) = Y_train.t() * Y_train;
        svd_econ(svd_U, singular_value, Phi_oldg, Y_train, "right");
        Phigg = Cgg = Phi_oldg.cols(0, K - 1);
        rho_cv(k, 0) = 10 * pow(singular_value[0], 2.0);
        Lambda2gg = Phigg * (diagmat(rho_cv(k, 0) - 1 / (rho_cv(k, 0) - 2 * pow(singular_value.subvec(0, K - 1), 2))));
        if(selected_tau1 != 0)
          spatpcaCore2(gram_matrix_Y_train.slice(k), Phigg, Cgg, Lambda2gg, Omega, selected_tau1, rho_cv(k, 0), maxit, tol);
        Phi_cv.slice(k) = Phigg;
        Lambd2_cv.slice(k) = Lambda2gg;    
      }
    }
    cv_score_tau1.zeros(1);
  }
  
  if(tau2.n_elem > 1) {
    mat cv2(M, tau2.n_elem);
    spatpcaCVPhi2 spatpcaCVPhi2(Y, gram_matrix_Y_train, Phi_cv, Lambd2_cv, rho_cv, K, selected_tau1, Omega, tau2, nk, maxit, tol, cv2, tempinv_cv); 
    RcppParallel::parallelFor(0, M, spatpcaCVPhi2);

    (sum(cv2, 0)).min(index2);
    selected_tau2 = tau2[index2];

    mat tempinv = inv_sympd((selected_tau1 * Omega) - gram_matrix_Y + (estimated_rho * Ip));
    mat estimated_R = estimated_Phi;
    mat estimated_Lambda1 = 0 * estimated_Phi;

    for(uword i = 0; i <= index2; i++)
      spatpcaCore3(tempinv, estimated_Phi, estimated_R, estimated_C, estimated_Lambda1, estimated_Lambda2, tau2[i], estimated_rho, maxit, tol);

    cv_score_tau2 = sum(cv2, 0) / M;
  }
  else {  
    selected_tau2 = max(tau2);
    if(selected_tau2 > 0) {
      mat tempinv = inv_sympd((selected_tau1 * Omega) - gram_matrix_Y + (estimated_rho * Ip));
      mat estimated_R= estimated_Phi;
      mat estimated_Lambda1 = 0 * estimated_Phi;
      for(uword i = 0; i < l2.n_elem; i++)
        spatpcaCore3(tempinv, estimated_Phi, estimated_R, estimated_C, estimated_Lambda1, estimated_Lambda2, l2[i], estimated_rho, maxit, tol);
    }
    cv_score_tau2.zeros(1);
  }
  
  spatpcaCVPhi3 spatpcaCVPhi3(Y, gram_matrix_Y_train, Phi_cv, Lambd2_cv, rho_cv, tempinv_cv, index2, K, Omega, selected_tau1, tau2, gamma, nk, maxit, tol, cv3);
  RcppParallel::parallelFor(0, M, spatpcaCVPhi3);
  if(gamma.n_elem > 1) {
    (sum(cv3, 0)).min(index3);
    selected_gamma = gamma[index3];
  }
  else {
    selected_gamma = max(gamma);
  }
  cv_score_gamma = sum(cv3, 0) / M;
  
  return List::create(Named("cv_score_tau1") = cv_score_tau1,
                      Named("cv_score_tau2") = cv_score_tau2,
                      Named("cv_score_gamma") = cv_score_gamma,
                      Named("estimated_eigenfn") = estimated_Phi,
                      Named("selected_tau1") = selected_tau1,
                      Named("selected_tau2") = selected_tau2,
                      Named("selected_gamma") = selected_gamma);
}

//' Internal function: Spatial prediction
//' @keywords internal
//' @param phir A matrix of estimated eigenfunctions based on original locations
//' @param Yr A data matrix
//' @param gammar A gamma value
//' @param predicted_eignefunction A vector of values of an eigenfunction on new locations
//' @return A list of objects
//' \item{prediction}{A vector of spatial predictions}
//' \item{estimated_covariance}{An estimated covariance matrix.}
//' \item{eigenvalue}{A vector of estimated eigenvalues.}
//' \item{error}{Error rate for the ADMM algorithm}
// [[Rcpp::export]]
List spatialPrediction(NumericMatrix phir, NumericMatrix Yr, double gamma, NumericMatrix predicted_eignefunction) {
  int n = Yr.nrow(), p = phir.nrow(), K = phir.ncol(), p2 = predicted_eignefunction.nrow() ;
  mat phi(phir.begin(), p, K, false);
  mat predicted_phi(predicted_eignefunction.begin(), p2, K, false);
  mat Y(Yr.begin(), n, p, false);
  mat Vc, Vc2;
  vec Sc, Sc2;
  
  mat cov = Y.t() * Y / n;
  eig_sym(Sc, Vc, phi.t() * cov * phi);
  int tempL = K;
  double total_variance = trace(cov);
  double tempSc2 = accu(Sc);
  double temp_v = total_variance - tempSc2;
  double error = (temp_v + K * gamma) / (p - tempL);
  
  Sc2 = sort(Sc, "descend");
  Vc2 = Vc.cols(sort_index(Sc, "descend"));
  double tempSc = tempSc2, temp;
  mat Sct, Sctz;
  Sct.ones(Sc2.n_elem);
  Sctz.zeros(Sc2.n_elem);
  
  if(Sc2[0] > gamma) {
    error = (total_variance - tempSc + K * gamma) / (p - tempL);
    temp = Sc2[tempL - 1];  
    while(temp - gamma < error) {
      if(tempL == 1) {
        error = (total_variance - Sc2[0] + gamma) / (p - 1);
        break;
      }
      tempSc -= Sc2[tempL - 1];
      tempL --;
      error = (total_variance - tempSc + tempL * gamma) / (p - tempL);
      temp = Sc2[tempL - 1];
    }
    if(Sc2[0] - gamma < error)
      error = total_variance / p;
  }
  else
    error = total_variance / p;
  
  vec eigenvalue = max(Sc2 - (error + gamma) * Sct, Sctz);
  vec eigenvalue2 = eigenvalue + error;
  mat estimated_covariance = phi * Vc2 * diagmat(eigenvalue / eigenvalue2) * trans(predicted_phi * Vc2);
  mat prediction = Y * estimated_covariance;
  return List::create(Named("prediction") = prediction,
                      Named("estimated_covariance") = estimated_covariance,
                      Named("eigenvalue") = eigenvalue,
                      Named("error") = error);
}

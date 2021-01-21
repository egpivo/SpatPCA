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

arma::mat cubicMatrix(const vec z1) {
  vec h(z1.n_elem - 1);
  for(unsigned i = 0; i < z1.n_elem - 1; i++)
    h[i] = z1[i + 1] - z1[i];
  arma::mat Q, R;
  int p = z1.n_elem;
  Q.zeros(p - 2, p);
  R.zeros(p - 2, p - 2);
  
  for(unsigned j = 0; j < p - 2; ++j) {
    Q(j, j) = 1 / h[j];
    Q(j, j + 1) = -1 / h[j] - 1 / h[j + 1];
    Q(j, j + 2) = 1 / h[j + 1];
    R(j, j) = (h[j] + h[j + 1]) / 3;
  }
  for(unsigned j = 0; j < p - 3; ++j)
    R(j, j + 1) = R(j + 1, j) = h[j + 1] / 6;
  
  return(Q.t() * solve(R, Q));
}

struct tpm: public RcppParallel::Worker {
  const mat& P;
  mat& L;  
  int p;
  int d;
  tpm(const mat &P, mat& L, int p, int d): P(P), L(L), p(p), d(d) {}

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
            L(i, j) = r * r * log(r) / (8.0 * arma::datum::pi);
          }
          else if(d == 3) {
            double r = sqrt(pow(P(i, 0) - P(j, 0), 2) +
                            pow(P(i, 1) - P(j, 1), 2) +
                            pow(P(i, 2) - P(j, 2), 2));
            L(i, j) = -r / (8.0 * arma::datum::pi);
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
//' @description Produce a thin-plane spline matrix based on a given location matric
//' 
//' @param location A location matrix
//' @return A thin-plane spline matrix
//' @examples
//' pesudo_sequence <- seq(-5, 5, length = 5)
//' two_dim_location <- as.matrix(expand.grid(x = pesudo_sequence, y = pesudo_sequence))
//' thin_plate_matrix <- thinPlateMatrix(two_dim_location)
// [[Rcpp::export]]
arma::mat thinPlateMatrix(const arma::mat location) {
  int p = location.n_rows, d = location.n_cols;
  int total_size = p + d;
  arma::mat L, Lp, Ip;
  
  L.zeros(total_size + 1, total_size + 1);
  Ip.eye(total_size + 1, total_size + 1);
  tpm tpm(location, L, p, d);
  parallelFor(0, p, tpm);
  L = symmatu(L);
  Lp = inv(L + 1e-8 * Ip);   
  Lp.shed_cols(p, total_size);
  Lp.shed_rows(p, total_size);
  L.shed_cols(p, total_size);
  L.shed_rows(p, total_size);
  arma::mat result = Lp.t() * (L * Lp);
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
//' @return A predictive estimte matrix
//' @examples
//' pesudo_sequence <- seq(-5, 5, length = 2)
//' original_location <- as.matrix(expand.grid(x = pesudo_sequence, y = pesudo_sequence))
//' new_location <- matrix(c(0.1, 0.2), nrow = 1, ncol = 2)
//' Phi <- matrix(c(1, 0, 0, 0), nrow = 4, ncol = 1)
//' thin_plate_matrix <- eigenFunction(new_location, original_location, Phi)
// [[Rcpp::export]]
arma::mat eigenFunction(const arma::mat new_location, const arma::mat original_location, const arma::mat Phi) {
  arma::mat L;
  int p = original_location.n_rows, d = original_location.n_cols, K = Phi.n_cols;
  int total_size = p + d;
  L.zeros(total_size + 1, total_size + 1);
  tpm tpm(original_location, L, p, d);
  parallelFor(0, p, tpm);
  L = L + L.t();
  
  arma::mat Phi_star, para(total_size + 1, K);
  Phi_star.zeros(total_size + 1, K);
  Phi_star.rows(0, p - 1) = Phi;
  para = solve(L, Phi_star);
  int pnew = new_location.n_rows;
  arma::mat eigen_fn(pnew, K);
  double psum, r;

  for(unsigned newi = 0; newi < pnew ; newi++) {
    for(unsigned i = 0; i < K; i++) {
      psum = 0;
      for(unsigned j = 0; j < p; j++) {
        if(d == 1) {  
          r = norm(new_location.row(newi) - original_location.row(j), "f");
          if(r != 0)
            psum += para(j, i) * pow(r, 3)/12;
        }
        else if(d == 2) {
          r = sqrt(pow(new_location(newi, 0) - original_location(j, 0), 2) + 
            (pow(new_location(newi, 1) - original_location(j, 1), 2)));
          if(r != 0)
            psum += para(j, i) * r * r * log(r) / (8.0 * arma::datum::pi);
        }
        else if(d == 3) {
          double r = sqrt(pow(new_location(newi, 0) - original_location(j, 0), 2) +
                          pow(new_location(newi, 1) - original_location(j, 1), 2) +
                          pow(new_location(newi, 2) - original_location(j, 2), 2));
          if(r != 0)
            psum -= para(j, i) * r / (8.0 * arma::datum::pi);
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
void spatpcacore2(const mat YY, mat& Phi, mat& C, mat& Lambda2, const mat Omega, const double tau1, const double rho, const int maxit, const double tol) {
  int p = Phi.n_rows;
  int K = Phi.n_cols;
  int iter = 0;
  arma::mat Ip, Sigtau1, temp;
  arma::vec er(2);
  Ip.eye(p,p);
  Sigtau1 = tau1 * Omega - YY;
  
  arma::mat U, V, diff;
  arma::vec S;
  arma::mat Cold = C;
  arma::mat Lambda2old = Lambda2;
  arma::mat tempinv = arma::inv_sympd(2 * Sigtau1 + rho * Ip);
  for(iter = 0; iter < maxit; iter++) {
    Phi = tempinv*((rho * Cold) - Lambda2old);
    temp = Phi + (Lambda2old / rho);
    arma::svd_econ(U, S, V, temp);
    C = U.cols(0, V.n_cols - 1) * V.t();
    diff = Phi - C;
    Lambda2 = Lambda2old + rho * diff;
    
    er[0] = arma::norm(diff, "fro") / sqrt(p * K / 1.0);
    er[1] = arma::norm(C - Cold, "fro") / sqrt(p * K / 1.0);
    
    if(max(er) <= tol)
      break;
    Cold = C;
    Lambda2old = Lambda2;
  }
  
  iter++;
}

arma::mat spatpcacore2p(const arma::mat YY, arma::mat& C, arma::mat& Lambda2, const arma::mat Omega,const double tau1, const double rho, const int maxit, const double tol) {
  int p = C.n_rows;
  int K = C.n_cols;
  int iter = 0;
  arma::mat Ip, Sigtau1, temp;
  arma::vec er(2);
  Ip.eye(p,p);
  Sigtau1 = tau1 * Omega - YY;
  
  arma::mat U, V, diff;
  arma::vec S;
  arma::mat Phi;
  arma::mat Cold = C;
  arma::mat Lambda2old = Lambda2;
  arma::mat tempinv = arma::inv_sympd(2 * Sigtau1 + rho * Ip);
  for(iter = 0; iter < maxit; iter++) {
    Phi = tempinv * ((rho * Cold) - Lambda2old);
    temp = Phi + (Lambda2old / rho);
    arma::svd_econ(U, S, V, temp);
    C = U.cols(0, V.n_cols - 1) * V.t();
    diff = Phi - C;
    Lambda2 = Lambda2old + rho * diff;
    
    er[0] = arma::norm(diff, "fro") / sqrt(p * K / 1.0);
    er[1] = arma::norm(C - Cold, "fro") / sqrt(p * K / 1.0);
    
    if(max(er) <= tol)
      break;
    Cold = C;
    Lambda2old = Lambda2;
  }
  iter++;
  return(Phi);
}

void spatpcacore3(const arma::mat tempinv, arma::mat& Phi, arma::mat& R, arma::mat& C, arma::mat& Lambda1, arma::mat& Lambda2, const double tau2, double rho, const int maxit, const double tol) {
  int p = Phi.n_rows;
  int K = Phi.n_cols;
  int iter = 0;
  arma::mat temp, tau2onerho, zero, one;
  arma::vec er(4);
  
  zero.zeros(p, K);
  one.ones(p, K);
  tau2onerho = tau2 * one / rho;
  
  arma::mat U, V, diffPR, diffPC;
  arma::vec S;
  arma::mat Phiold = Phi;
  arma::mat Rold = R;
  arma::mat Cold = C;
  arma::mat Lambda1old = Lambda1;
  arma::mat Lambda2old = Lambda2;
  
  for(iter = 0; iter < maxit; iter++) {
    Phi = 0.5 * tempinv * (rho * (Rold + Cold) - (Lambda1old + Lambda2old));
    R = arma::sign(((Lambda1old / rho) + Phi)) % arma::max(zero, arma::abs(((Lambda1old / rho) + Phi)) - tau2onerho);
    temp = Phi + Lambda2old / rho;
    arma::svd_econ(U, S, V, temp);
    C = U.cols(0, V.n_cols - 1) * V.t();
    diffPR = Phi - R;
    diffPC = Phi - C;    
    Lambda1 = Lambda1old + rho * diffPR;
    Lambda2 = Lambda2old + rho * diffPC;
    
    er[1] = arma::norm(diffPR, "fro") / sqrt(p / 1.0);
    er[2] = arma::norm(R-Rold, "fro") / sqrt(p / 1.0);
    er[3] = arma::norm(diffPC, "fro") / sqrt(p / 1.0);
    er[4] = arma::norm(C-Cold, "fro") / sqrt(p / 1.0);
    
    if(max(er) <= tol)
      break;
    Phiold = Phi;
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
  cube& YYtrain;
  cube& Phicv;
  cube& Lmbd2cv;
  mat& rho;
  spatpcaCVPhi(const mat& Y, int K, const mat& Omega, const vec& tau1, const vec& nk, int maxit, double tol, mat& output, cube& YYtrain, cube& Phicv, cube& Lmbd2cv, mat& rho): Y(Y), K(K), Omega(Omega), tau1(tau1), nk(nk), maxit(maxit), tol(tol), output(output), YYtrain(YYtrain), Phicv(Phicv), Lmbd2cv(Lmbd2cv), rho(rho) {}

  void operator()(std::size_t begin, std::size_t end) {
    arma::mat Ip;
    Ip.eye(Y.n_cols,Y.n_cols);
    for(std::size_t k = begin; k < end; k++) {
      arma::mat UPhi, Phiold,Phi,C, Lambda2;
      vec SPhi;
      mat Ytrain = Y.rows(arma::find(nk != (k + 1)));
      mat Yvalid = Y.rows(arma::find(nk == (k + 1)));
      
      arma::svd_econ(UPhi, SPhi, Phiold, Ytrain, "right");  
      rho(k, 0) = 10 * pow(SPhi[0], 2.0);
      
      Phicv.slice(k) = Phiold.cols(0, K - 1);
      Phi = C = Phicv.slice(k);
      Lambda2 = Lmbd2cv.slice(k) = Phi*(diagmat(rho(k, 0) - 1 / (rho(k, 0) - 2 * pow(SPhi.subvec(0, K - 1), 2))));
      output(k, 0) = pow(arma::norm(Yvalid * (Ip - (Phicv.slice(k)) * (Phicv.slice(k)).t()), "fro"), 2.0); 
      YYtrain.slice(k) = Ytrain.t() * Ytrain;
      for(uword i = 1; i < tau1.n_elem; i++) {
        spatpcacore2(YYtrain.slice(k), Phi, C, Lambda2, Omega, tau1[i], rho(k, 0), maxit, tol);
        output(k, i) = pow(arma::norm(Yvalid * (Ip - Phi * Phi.t()), "fro"), 2.0); 
      }
    }
  }
};

struct spatpcaCVPhi2: public RcppParallel::Worker {
  const mat& Y;
  const cube& YYtrain;
  cube& Phicv;
  cube& Lmbd2cv;
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
    const cube& YYtrain,
    cube& Phicv,
    cube& Lmbd2cv,
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
    Y(Y), YYtrain(YYtrain), Phicv(Phicv), Lmbd2cv(Lmbd2cv), rho(rho),K(K), tau1(tau1),Omega(Omega), tau2(tau2), nk(nk), maxit(maxit), tol(tol), output(output), tempinv(tempinv) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t k = begin; k < end; k++) {
      arma::mat Phi, R, C, Lambda1,Ip, Lambda2;
      mat Yvalid = Y.rows(arma::find(nk == (k + 1)));
      Ip.eye(Y.n_cols,Y.n_cols);
      Phi = C = Phicv.slice(k);
      Lambda1 = 0 * Phicv.slice(k);
      Lambda2 = Lmbd2cv.slice(k);
      if(tau1 != 0)
        spatpcacore2(YYtrain.slice(k), Phi,C, Lambda2, Omega, tau1, rho(k, 0), maxit, tol);
      R = Phi;
      Phicv.slice(k) = Phi;
      Lmbd2cv.slice(k) = Lambda2;
      tempinv.slice(k) = arma::inv_sympd((tau1 * Omega) - YYtrain.slice(k) + (rho(k, 0) * Ip));
      for(uword i = 0; i < tau2.n_elem; i++) {
        spatpcacore3(tempinv.slice(k), Phi, R, C, Lambda1, Lambda2, tau2[i], rho(k, 0), maxit, tol);
        output(k,i) = pow(arma::norm(Yvalid * (Ip - Phi * Phi.t()), "fro"), 2.0); 
      }
    }
  }
};

struct spatpcaCVPhi3: public RcppParallel::Worker {
  const mat& Y;
  const cube& YYtrain;
  const cube& Phicv;
  const cube& Lmbd2cv;
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
    const cube& YYtrain,
    const cube& Phicv,
    const cube& Lmbd2cv,
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
    Y(Y), YYtrain(YYtrain), Phicv(Phicv), Lmbd2cv(Lmbd2cv), rho(rho), tempinv(tempinv), index(index),K(K), Omega(Omega), tau1(tau1), tau2(tau2), gamma(gamma), nk(nk), maxit(maxit), tol(tol), output(output) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t k = begin; k < end; k++) {
      int p = Y.n_cols, tempL;
      arma::mat  Phi, C, Ip, Lambda2;
      double totalvar, err, temp, tempSc, tempSc2;
      mat Yvalid = Y.rows(arma::find(nk == (k+1)));
      Ip.eye(Y.n_cols,Y.n_cols);
      Phi = C = Phicv.slice(k);
      int K = Phi.n_cols; 
      Lambda2 = Lmbd2cv.slice(k);
      
      if(max(tau2) != 0 || max(tau1) != 0) {
        if(tau2.n_elem == 1) {
          spatpcacore2(YYtrain.slice(k), Phi, C, Lambda2, Omega, tau1, rho(k, 0), maxit, tol);
        }
        else {
          mat R = Phi;
          mat Lambda1 = 0 * Phicv.slice(k);
          for(uword i = 0; i <= index; i++) {
            spatpcacore3(tempinv.slice(k), Phi, R, C, Lambda1, Lambda2, tau2[i], rho(k, 0), maxit, tol);
          }
        }
      }
      else {
        Phi = Phicv.slice(k);
      }
      arma::mat Vc, Vc2, covtrain, covvalid, estimated_covariance, vec;
      arma::vec Sc, Sc2, Sct, Sctz, cv;
      arma::mat eigenvalue;
      covtrain = YYtrain.slice(k) / (Y.n_rows - Yvalid.n_rows);
      covvalid = arma::trans(Yvalid) * Yvalid / Yvalid.n_rows;
      totalvar = arma::trace(covtrain);
      arma::eig_sym(Sc, Vc, trans(Phi) * covtrain * Phi);
      tempSc2 = accu(Sc);
      Sc2 = sort(Sc,"descend");
      Vc2 = Vc.cols(sort_index(Sc,"descend"));
      Sct.ones(Sc2.n_elem);
      Sctz.zeros(Sc2.n_elem);

      for(unsigned int gj = 0; gj < gamma.n_elem; gj++) {
        tempSc = tempSc2;
        tempL = K;
        if(Sc2[0] > gamma[gj]) {
          err = (totalvar - tempSc + K * gamma[gj]) / (p - tempL);
          temp = Sc2[tempL - 1];
          while(temp-gamma[gj] < err) {
            if(tempL == 1) {
              err = (totalvar - Sc2[0] + gamma[gj]) / (p - 1);
              break;
            }
            tempSc -= Sc2[tempL - 1];
            tempL --;
            err = (totalvar - tempSc + tempL * gamma[gj]) / (p - tempL);
            temp = Sc2[tempL - 1];
          }
          if(Sc2[0]-gamma[gj] < err)
            err = totalvar / p;
        }
        else
          err = totalvar / p;
        eigenvalue = arma::max(Sc2 - (err + gamma[gj]) * Sct, Sctz);
        estimated_covariance =  Phi * Vc2 * diagmat(eigenvalue) * trans(Vc2) * trans(Phi);
        output(k,gj) = pow(arma::norm(covvalid - estimated_covariance - err * Ip, "fro"), 2.0);
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
List spatpcaCV(NumericMatrix sxyr, NumericMatrix Yr, int M, int K, NumericVector tau1r, NumericVector tau2r, NumericVector gammar, NumericVector nkr, int maxit, double tol, NumericVector l2r) {
  int n = Yr.nrow(), p = Yr.ncol(), d = sxyr.ncol();
  arma::mat Y(Yr.begin(), n, p, false);
  arma::mat sxy(sxyr.begin(), p, d, false);
  colvec tau1(tau1r.begin(), tau1r.size(),false);
  colvec tau2(tau2r.begin(), tau2r.size(), false);
  colvec gamma(gammar.begin(), gammar.size(), false);
  colvec nk(nkr.begin(), nkr.size(), false);
  colvec l2(l2r.begin(), l2r.size(), false);
  arma::mat cv(M, tau1.n_elem), cv3(M, gamma.n_elem), cv_score_tau1, cv_score_tau2, cv_score_gamma;
  arma::mat Omega;
  double selected_tau1, selected_tau2 = 0, selected_gamma;
  arma::mat YYest = Y.t()*Y;
  arma::mat UPhiest, Phiest2;
  vec SPhiest;
  arma::svd_econ(UPhiest, SPhiest, Phiest2, Y, "right");
  double rhoest = 10 * pow(SPhiest[0], 2.0);
  arma::mat Phiest = Phiest2.cols(0, K - 1); 
  arma::mat Cest = Phiest;
  arma::mat Lambda2est = Phiest * (diagmat(rhoest - 1 / (rhoest-2 * pow(SPhiest.subvec(0, K - 1), 2))));
  arma::mat rhocv(M, 1);
  arma::cube YYtrain(p, p, M), Phicv(p, K, M), Lmbd2cv(p, K, M), tempinvcv(p, p, M);
  uword index1, index2 = 0, index3;
  mat Ip;
  Ip.eye(Y.n_cols,Y.n_cols);
  if(max(tau1) != 0 || max(tau2) != 0) {
    Omega = thinPlateMatrix(sxy) + 1e-8 * Ip;
  }
  else {
    if(gamma.n_elem > 1) {
      Omega = Ip;
      mat Ytrain, Yvalid; 
      arma::mat UPhi, Phioldg;
      vec SPhi;
      for(unsigned k = 0; k < M; ++k) {
        Ytrain = Y.rows(arma::find(nk != (k+1)));
        Yvalid = Y.rows(arma::find(nk == (k+1)));
        arma::svd_econ(UPhi, SPhi, Phioldg, Ytrain, "right");
        Phicv.slice(k) = Phioldg.cols(0, K - 1);
        YYtrain.slice(k) = Ytrain.t() * Ytrain;
      }
    }
  }

  if(tau1.n_elem > 1) {  
    spatpcaCVPhi spatpcaCVPhi(Y, K, Omega, tau1, nk, maxit, tol, cv, YYtrain, Phicv, Lmbd2cv, rhocv);
    RcppParallel::parallelFor(0, M, spatpcaCVPhi);
    (sum(cv, 0)).min(index1);
    selected_tau1 = tau1[index1];  
    if(index1 > 0)
      Phiest = spatpcacore2p(YYest, Cest, Lambda2est, Omega, selected_tau1, rhoest, maxit, tol);
    cv_score_tau1 = sum(cv, 0) / M;
  }
  else {
    selected_tau1 = max(tau1);
    if(selected_tau1 != 0 && max(tau2) == 0) {
      Phiest = spatpcacore2p(YYest, Cest, Lambda2est, Omega, selected_tau1, rhoest, maxit, tol);
      arma::mat Phigg, Cgg, Lambda2gg;
      mat Yvalid, Ytrain;
      arma::mat UPhi, Phioldg;
      vec SPhi;

      for(unsigned k = 0; k < M; ++k) {
        Ytrain = Y.rows(arma::find(nk != (k + 1)));
        YYtrain.slice(k) = Ytrain.t()*Ytrain;
        arma::svd_econ(UPhi, SPhi, Phioldg, Ytrain, "right");
        Phigg = Cgg = Phioldg.cols(0, K - 1);
        rhocv(k, 0) = 10 * pow(SPhi[0], 2.0);
        Lambda2gg = Phigg * (diagmat(rhocv(k, 0) - 1 / (rhocv(k, 0) - 2 * pow(SPhi.subvec(0, K - 1), 2))));
        spatpcacore2(YYtrain.slice(k), Phigg,Cgg, Lambda2gg, Omega, selected_tau1, rhocv(k, 0), maxit, tol);
        Phicv.slice(k) = Phigg;
        Lmbd2cv.slice(k) = Lambda2gg;
        tempinvcv.slice(k) = arma::inv_sympd((selected_tau1 * Omega) - YYtrain.slice(k) + (rhocv(k, 0) * Ip));    
      }
    }
    else if(selected_tau1 == 0 && max(tau2) == 0) {
      arma::mat UPhi2, Phioldc;
      vec SPhi2;
      arma::svd_econ(UPhi2, SPhi2, Phioldc, Y, "right");
      Phiest = Phioldc.cols(0,K - 1); 
    }
    else {
      arma::mat Phigg, Cgg, Lambda2gg;
      mat Yvalid, Ytrain;
      arma::mat UPhi, Phioldg;
      vec SPhi;

      for(unsigned k = 0; k < M; ++k) {
        Ytrain = Y.rows(arma::find(nk != (k + 1)));
        YYtrain.slice(k) = Ytrain.t() * Ytrain;
        arma::svd_econ(UPhi, SPhi, Phioldg, Ytrain, "right");
        Phigg = Cgg = Phioldg.cols(0, K - 1);
        rhocv(k, 0) = 10 * pow(SPhi[0], 2.0);
        Lambda2gg = Phigg * (diagmat(rhocv(k, 0) - 1 / (rhocv(k, 0) - 2 * pow(SPhi.subvec(0, K - 1), 2))));
        if(selected_tau1 != 0)
          spatpcacore2(YYtrain.slice(k), Phigg,Cgg, Lambda2gg, Omega, selected_tau1, rhocv(k, 0), maxit, tol);
        Phicv.slice(k) = Phigg;
        Lmbd2cv.slice(k) = Lambda2gg;    
      }
    }
    cv_score_tau1.zeros(1);
  }
  
  if(tau2.n_elem > 1) {
    arma::mat cv2(M, tau2.n_elem);
    spatpcaCVPhi2 spatpcaCVPhi2(Y, YYtrain, Phicv, Lmbd2cv, rhocv, K, selected_tau1, Omega, tau2, nk, maxit, tol, cv2, tempinvcv); 
    RcppParallel::parallelFor(0, M, spatpcaCVPhi2);

    (sum(cv2, 0)).min(index2);
    selected_tau2 = tau2[index2];

    mat tempinv = arma::inv_sympd((selected_tau1 * Omega) - YYest + (rhoest * Ip));
    mat Rest = Phiest;
    mat Lambda1est = 0 * Phiest;

    for(uword i = 0; i <= index2; i++)
      spatpcacore3(tempinv, Phiest, Rest, Cest, Lambda1est, Lambda2est, tau2[i], rhoest, maxit, tol);

    cv_score_tau2 = sum(cv2, 0) / M;
  }
  else {  
    selected_tau2 = max(tau2);
    if(selected_tau2 > 0) {
      mat tempinv = arma::inv_sympd((selected_tau1 * Omega) - YYest + (rhoest * Ip));
      mat Rest= Phiest;
      mat Lambda1est = 0 * Phiest;
      for(uword i = 0; i < l2.n_elem; i++)
        spatpcacore3(tempinv, Phiest, Rest, Cest, Lambda1est, Lambda2est, l2[i], rhoest, maxit, tol);
    }
    cv_score_tau2.zeros(1);
  }
  
  spatpcaCVPhi3 spatpcaCVPhi3(Y, YYtrain, Phicv, Lmbd2cv, rhocv, tempinvcv, index2, K, Omega, selected_tau1, tau2, gamma, nk, maxit, tol, cv3);
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
                      Named("estimated_eigenfn") = Phiest,
                      Named("selected_tau1") = selected_tau1,
                      Named("selected_tau2") = selected_tau2,
                      Named("selected_gamma") = selected_gamma);
}

//' Internal function: Spatial prediction
//' @keywords internal
//' @param phir A matrix of estimated eigenfunctions based on original locations
//' @param Yr A data matrix
//' @param gammar A gamma value
//' @param phi2r A vector of values of an eigenfunction on new locations
//' @return A list of objects
//' \item{prediction}{A vector of spatial predicitons}
//' \item{estimated_covariance}{An estimated covariance matrix.}
//' \item{eigenvalue}{A vecotor of estimated eigenvalues.}
//' \item{error}{Error rate for the ADMM algorithm}
// [[Rcpp::export]]
List spatialPrediction(NumericMatrix phir, NumericMatrix Yr, double gamma, NumericMatrix phi2r) {
  int n = Yr.nrow(), p = phir.nrow(), K = phir.ncol(), p2 = phi2r.nrow() ;
  arma::mat phi(phir.begin(), p, K, false);
  arma::mat phi2(phi2r.begin(), p2, K, false);
  arma::mat Y(Yr.begin(), n, p, false);
  arma::mat Vc, Vc2;
  arma::vec Sc, Sc2;
  
  arma::mat cov = Y.t() * Y / n;
  arma::eig_sym(Sc, Vc, phi.t() * cov * phi);
  int tempL = K;
  double totalvar = trace(cov);
  double tempSc2 = accu(Sc);
  double temp_v = totalvar - tempSc2;
  double err = (temp_v + K * gamma) / (p - tempL);
  
  Sc2 = sort(Sc, "descend");
  Vc2 = Vc.cols(sort_index(Sc, "descend"));
  double tempSc = tempSc2, temp;
  arma::mat Sct, Sctz;
  Sct.ones(Sc2.n_elem);
  Sctz.zeros(Sc2.n_elem);
  
  if(Sc2[0] > gamma) {
    err = (totalvar - tempSc + K * gamma) / (p - tempL);
    temp = Sc2[tempL - 1];  
    while(temp - gamma < err) {
      if(tempL == 1) {
        err = (totalvar - Sc2[0] + gamma) / (p - 1);
        break;
      }
      tempSc -= Sc2[tempL - 1];
      tempL --;
      err = (totalvar - tempSc + tempL * gamma) / (p - tempL);
      temp = Sc2[tempL - 1];
    }
    if(Sc2[0] - gamma < err)
      err = totalvar / p;
  }
  else
    err = totalvar / p;
  
  arma::vec eigenvalue = arma::max(Sc2 - (err + gamma) * Sct, Sctz);
  arma::vec eigenvalue2 = eigenvalue + err;
  arma::mat estimated_covariance = phi * Vc2 * diagmat(eigenvalue / eigenvalue2) * trans(phi2 * Vc2);
  arma::mat prediction = Y * estimated_covariance;
  return List::create(Named("prediction") = prediction,
                      Named("estimated_covariance") = estimated_covariance,
                      Named("eigenvalue") = eigenvalue,
                      Named("error") = err);
}

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

arma::mat cubicmatrix(const vec z1){
  
  //vec h = diff(z1);
  vec h(z1.n_elem-1);
  for(unsigned i = 0; i < z1.n_elem-1; i++)
    h[i] = z1[i+1]-z1[i];
  arma::mat Q, R;
  int p = z1.n_elem;
  Q.zeros(p-2,p);
  R.zeros(p-2,p-2);
  
  for(unsigned j = 0; j < p-2; ++j){
    Q(j,j) = 1/h[j];
    Q(j,j+1) = -1/h[j]-1/h[j+1];
    Q(j,j+2) = 1/h[j+1];
    R(j,j) = (h[j]+h[j+1])/3;
  }
  for(unsigned j = 0; j < p-3; ++j)
    R(j,j+1) = R(j+1,j) = h[j+1]/6;
  
  return(Q.t()*solve(R,Q));
}

using namespace arma;
using namespace Rcpp;
using namespace std;

struct tpm: public RcppParallel::Worker {
  const mat& P;
  mat& L;  
  int p;
  int d;
  tpm(const mat &P, mat& L, int p, int d) : P(P), L(L), p(p), d(d){}
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++){
      for(unsigned j = 0; j < p; ++j){
        if(j >i){
          if(d==2){  
            double r  = sqrt(pow(P(i,0)-P(j,0),2)+(pow(P(i,1)-P(j,1),2)));
            L(i,j) = r*r*log(r)/(8.0*arma::datum::pi);
          }
          else{
            double r  = sqrt(pow(P(i,0)-P(j,0),2));
            L(i,j) = sqrt(2)/(16*sqrt(arma::datum::pi))*pow(r,3);
          }  
        }
      }  
      
      L(i,p) = 1;
      for(unsigned k = 0; k < d; ++k){
        L(i,p+k+1) = P(i,k);
      }
    }
  }
};

arma::mat tpmatrix(const arma::mat P){
  arma::mat L, Lp;
  int p = P.n_rows, d = P.n_cols;
  
  L.zeros(p+d+1, p+d+1);
  tpm tpm(P,L,p,d);
  parallelFor(0, p,tpm);
  L = L + L.t();
  Lp = inv(L);   
  Lp.shed_cols(p,p+2);
  Lp.shed_rows(p,p+2);
  L.shed_cols(p,p+2);
  L.shed_rows(p,p+2);
  arma::mat result = Lp.t()*(L*Lp);
  return(result);
}

// [[Rcpp::export]]
arma::mat tpm2(const arma::mat z,const arma::mat P, const arma::mat Phi){
  arma::mat L;//, Lp;
  int p = P.n_rows, d = P.n_cols, K = Phi.n_cols;
  L.zeros(p+d+1, p+d+1);
  tpm tpm(P,L,p,d);
  parallelFor(0, p,tpm);
  L = L + L.t();
  // Lp = inv(L);  
  arma::mat Phi_star, para(p+d+1, K);
  Phi_star.zeros(p+d+1, K);
  Phi_star.rows(0,p-1) = Phi;
  para = solve(L, Phi_star);
  int pnew = z.n_rows;
  arma::mat eigen_fn(pnew, K);
  double psum, r;
  for(unsigned newi = 0; newi < pnew ; newi++){
    for(unsigned i = 0; i < K; i++){
      psum = 0;
      for(unsigned j = 0; j < p; j++){
        if(d==2){
          r  = sqrt(pow(z(newi,0)-P(j,0),2)+(pow(z(newi,1)-P(j,1),2)));
          if(r!=0)
            psum += para(j,i)* r*r*log(r)/(8.0*arma::datum::pi);
        }
        else{
          r = norm(z.row(newi)-P.row(j),'f');
          if(r!=0)
            psum += para(j,i)*((sqrt(2)/(16*sqrt(arma::datum::pi)))*pow(r,3));
        }
      }
      if(d==1)
        eigen_fn(newi,i) = psum + para(p+1,i)*z(newi,0) + para(p,i);
      else
        eigen_fn(newi,i) = psum + para(p+1,i)*z(newi,0) + para(p+2,i)*z(newi,1) + para(p,i); 
    }
  }
  return(eigen_fn);
}

// user includes
using namespace arma;
using namespace Rcpp;
using namespace std;
void spatpcacore2(const mat YY, mat& Phi,  mat& C, mat& Lambda2, const mat Omega,const double tau1, const double rho,const int maxit,const double tol){
  
  int p = Phi.n_rows;
  int K = Phi.n_cols;
  int iter = 0;
  arma::mat Ip, Sigtau1, temp;
  arma::vec er(2);
  Ip.eye(p,p);
  Sigtau1 = tau1*Omega - YY;
  
  arma::mat U;
  arma::vec S;
  arma::mat V;
  arma::mat Cold = C;
  arma::mat Lambda2old = Lambda2;
  arma::mat tempinv = arma::inv_sympd(2*Sigtau1 + rho*Ip);
  for (iter = 0; iter < maxit; iter++){
    Phi = tempinv*((rho*Cold)-Lambda2old);
    temp = Phi + (Lambda2old/rho);
    arma::svd_econ(U, S, V,temp);
    C = U.cols(0,V.n_cols-1)*V.t();
    
    Lambda2 = Lambda2old +rho*(Phi-C);
    
    er[0] = arma::norm(Phi-C,"fro")/sqrt(p*K/1.0);
    er[1] = arma::norm((C-Cold),"fro")/sqrt(p*K/1.0);
    
    if(max(er) <= tol)
      break;
    Cold = C;
    Lambda2old = Lambda2;
  }
  
  iter++;
  //  if(iter == maxit)
  //    Rcpp::Rcout<<"Not converge at tau1="<<tau1<<"\n"<<std::endl;
}

using namespace arma;
using namespace Rcpp;
using namespace std;

arma::mat spatpcacore2p(const arma::mat YY, arma::mat& C, arma::mat& Lambda2, const arma::mat Omega,const double tau1, const double rho,const int maxit,const double tol){
  
  int p = C.n_rows;
  int K = C.n_cols;
  int iter = 0;
  arma::mat Ip, Sigtau1, temp;
  arma::vec er(2);
  Ip.eye(p,p);
  Sigtau1 = tau1*Omega - YY;
  
  arma::mat U;
  arma::vec S;
  arma::mat V;
  arma::mat Phi;
  arma::mat Cold = C;
  arma::mat Lambda2old = Lambda2;
  arma::mat tempinv = arma::inv_sympd(2*Sigtau1 + rho*Ip);
  for (iter = 0; iter < maxit; iter++){
    Phi = tempinv*((rho*Cold)-Lambda2old);
    temp = Phi + (Lambda2old/rho);
    arma::svd_econ(U, S, V,temp);
    C = U.cols(0,V.n_cols-1)*V.t();
    
    Lambda2 = Lambda2old +rho*(Phi-C);
    
    er[0] = arma::norm(Phi-C,"fro")/sqrt(p*K/1.0);
    er[1] = arma::norm((C-Cold),"fro")/sqrt(p*K/1.0);
    
    if(max(er) <= tol)
      break;
    Cold = C;
    Lambda2old = Lambda2;
  }
  
  iter++;
  //  if(iter == maxit)
  //    Rcpp::Rcout<<"Not converge at tau1="<<tau1<<"\n"<<std::endl;
  // cout<<"iter"<<iter<<endl;
  return(Phi);
}



void spatpcacore3(const arma::mat tempinv, arma::mat& Phi, arma::mat& R,  arma::mat& C,  arma::mat& Lambda1, arma::mat& Lambda2, const double tau2, double rho,const int maxit,const double tol){
  
  int p = Phi.n_rows;
  int K = Phi.n_cols;
  int iter = 0;
  arma::mat temp, tau2onerho, zero, one;
  arma::vec er(4);
  
  zero.zeros(p,K);
  one.ones(p,K);
  tau2onerho = tau2*one/rho;
  //Sigtau1 = tau1*Omega - arma::trans(Y)*Y;
  
  arma::mat U;
  arma::vec S;
  arma::mat V;
  arma::mat Phiold = Phi;
  arma::mat Rold = R;
  arma::mat Cold = C;
  arma::mat Lambda1old = Lambda1;
  arma::mat Lambda2old = Lambda2;
  //  arma::mat tempinv = arma::inv_sympd(Sigtau1 + (rho*Ip));
  for (iter = 0; iter < maxit; iter++){
    // Phi =0.5*arma::solve(Sigtau1 + rho*Ip,(rho*(Rold+Cold)-Lambda1old-Lambda2old));
    Phi =0.5*tempinv*(rho*(Rold+Cold)-(Lambda1old+Lambda2old));
    R = arma::sign(((Lambda1old/rho)+Phi))%arma::max(zero, arma::abs(((Lambda1old/rho)+Phi)) - tau2onerho);
    temp = Phi+Lambda2old/rho;
    arma::svd_econ(U, S, V,temp);
    C = U.cols(0,V.n_cols-1)*V.t();
    
    Lambda1 = Lambda1old +rho*(Phi-R);
    Lambda2 = Lambda2old +rho*(Phi-C);
    
    
    er[1] = arma::norm(Phi-R,"fro")/sqrt(p/1.0);
    er[2] = arma::norm((R-Rold),"fro")/sqrt(p/1.0);
    er[3] = arma::norm(Phi-C,"fro")/sqrt(p/1.0);
    er[4] = arma::norm((C-Cold),"fro")/sqrt(p/1.0);
    
    if(max(er) <= tol)
      break;
    Phiold = Phi;
    Rold = R;
    Cold = C;
    Lambda1old = Lambda1;
    Lambda2old = Lambda2;
  }
  
  iter++;
  // if(iter == maxit)
  //    Rcpp::Rcout<<"Not converge at tau2="<<tau2<<"\n"<<std::endl;
}

using namespace arma;
using namespace Rcpp;
using namespace std;


struct spatpcacv_p: public RcppParallel::Worker {
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
  spatpcacv_p(const mat& Y, int K, const mat& Omega, const vec& tau1, const vec& nk, int maxit, double tol, mat& output, cube& YYtrain, cube& Phicv, cube& Lmbd2cv, mat& rho) : Y(Y), K(K), Omega(Omega), tau1(tau1), nk(nk),  maxit(maxit),tol(tol),output(output), YYtrain(YYtrain), Phicv(Phicv), Lmbd2cv(Lmbd2cv), rho(rho){}
  void operator()(std::size_t begin, std::size_t end) {
    arma::mat Ip;
    Ip.eye(Y.n_cols,Y.n_cols);
    for(std::size_t k = begin; k < end; k++){
      arma::mat UPhi, Phiold,Phi,C, Lambda2;
      vec SPhi;
      mat Ytrain = Y.rows(arma::find(nk!=(k+1)));
      mat Yvalid = Y.rows(arma::find(nk==(k+1)));
      
      arma::svd_econ(UPhi, SPhi, Phiold, Ytrain, "right");
    
      rho(k,0) = 10*pow(SPhi[0],2.0);
      Phicv.slice(k) = Phiold.cols(0,K-1);
      Phi = C = Phicv.slice(k);
      Lambda2 = Lmbd2cv.slice(k) = Phi*(diagmat(rho(k,0) -1/(rho(k,0)-2*pow(SPhi.subvec(0,K-1),2))));
      output(k,0) = pow(arma::norm(Yvalid*(Ip-(Phicv.slice(k))*(Phicv.slice(k)).t()),"fro"),2.0); 
      YYtrain.slice(k) = Ytrain.t()*Ytrain;
      for(uword  i = 1; i < tau1.n_elem; i++){
        spatpcacore2(YYtrain.slice(k),Phi,C,Lambda2, Omega, tau1[i], rho(k,0), maxit,tol);
        output(k,i) = pow(arma::norm(Yvalid*(Ip-Phi*Phi.t()),"fro"),2.0); 
      }
    }
    
  }
};

struct spatpcacv_p2: public RcppParallel::Worker {
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
  
  spatpcacv_p2(const mat& Y, const cube& YYtrain, cube& Phicv, cube& Lmbd2cv, const mat& rho, int K, double tau1, const mat& Omega, const vec& tau2, const vec& nk, int maxit, double tol, mat& output, cube& tempinv) : Y(Y), YYtrain(YYtrain), Phicv(Phicv), Lmbd2cv(Lmbd2cv), rho(rho),K(K), tau1(tau1),Omega(Omega), tau2(tau2), nk(nk),  maxit(maxit),tol(tol),output(output), tempinv(tempinv) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t k = begin; k < end; k++){
      arma::mat  Phi, R, C, Lambda1,Ip, Lambda2;
      mat Yvalid = Y.rows(arma::find(nk==(k+1)));
      Ip.eye(Y.n_cols,Y.n_cols);
      Phi = C = Phicv.slice(k);
      Lambda1 = 0*Phicv.slice(k);
      Lambda2 = Lmbd2cv.slice(k);
      if(tau1 !=0)
        spatpcacore2(YYtrain.slice(k), Phi,C, Lambda2, Omega, tau1, rho(k,0), maxit,tol);
      R = Phi;
      Phicv.slice(k) = Phi;
      Lmbd2cv.slice(k) = Lambda2;
      // output(k,0) = pow(arma::norm(Yvalid*(Ip-Phi*Phi.t()),"fro"),2.0); 
      tempinv.slice(k) = arma::inv_sympd((tau1*Omega) - YYtrain.slice(k) + (rho(k,0)*Ip));
      for(uword  i = 0; i < tau2.n_elem; i++){
        spatpcacore3(tempinv.slice(k), Phi,R,C,Lambda1,Lambda2, tau2[i], rho(k,0), maxit,tol);
        output(k,i) = pow(arma::norm(Yvalid*(Ip-Phi*Phi.t()),"fro"),2.0); 
      }
    }
    
  }
};

struct spatpcacv_p3: public RcppParallel::Worker {
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
  
  
  spatpcacv_p3(const mat& Y, const cube& YYtrain, const cube& Phicv, const cube& Lmbd2cv, const mat& rho,  const cube& tempinv, const uword index,int K,  const mat& Omega, const double tau1, const vec& tau2,const vec& gamma, const vec& nk, int maxit, double tol, mat& output) : Y(Y), YYtrain(YYtrain), Phicv(Phicv), Lmbd2cv(Lmbd2cv), rho(rho), tempinv(tempinv), index(index),K(K), Omega(Omega), tau1(tau1),tau2(tau2), gamma(gamma), nk(nk),  maxit(maxit),tol(tol),output(output) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t k = begin; k < end; k++){
      int p = Y.n_cols, tempL;
      arma::mat  Phi, C, Ip, Lambda2;
      double totalvar, err, temp, tempSc, tempSc2;
      mat Yvalid = Y.rows(arma::find(nk==(k+1)));
      Ip.eye(Y.n_cols,Y.n_cols);
      Phi = C = Phicv.slice(k);
      int K = Phi.n_cols; 
      Lambda2 = Lmbd2cv.slice(k);
     
      if(max(tau2) != 0 || max(tau1) !=0){
        if(tau2.n_elem == 1){
          spatpcacore2(YYtrain.slice(k), Phi,C, Lambda2, Omega, tau1, rho(k,0), maxit,tol);
        }
        else{
          mat R = Phi;
          mat Lambda1 = 0*Phicv.slice(k);
          for(uword  i = 0; i <= index; i++){
            spatpcacore3(tempinv.slice(k),Phi,R,C,Lambda1,Lambda2, tau2[i], rho(k,0), maxit,tol);
          }
        }
      }
      else{
        Phi = Phicv.slice(k);
      }
      arma::mat  Vc, Vc2, covtrain, covvalid, covest, vec;
      arma::vec Sc, Sc2, Sct, Sctz, cv;
      arma::mat eigenvalue;
      covtrain = YYtrain.slice(k)/(Y.n_rows - Yvalid.n_rows);
      covvalid = arma::trans(Yvalid)*Yvalid/Yvalid.n_rows;
      totalvar = arma::trace(covtrain);  
      arma::eig_sym(Sc, Vc, trans(Phi)*covtrain*Phi);
      tempSc2 = accu(Sc);
      Sc2 = sort(Sc,"descend");
      Vc2 = Vc.cols(sort_index(Sc,"descend"));
      Sct.ones(Sc2.n_elem);
      Sctz.zeros(Sc2.n_elem);
    
      for(unsigned int gj = 0; gj < gamma.n_elem; gj++){
        tempSc = tempSc2;
        tempL = K;
        if(Sc2[0]> gamma[gj]){
          err = (totalvar - tempSc + K*gamma[gj])/(p - tempL);
          temp = Sc2[tempL-1];  
          while( temp-gamma[gj] < err){
            if(tempL == 1){
              err = (totalvar - Sc2[0] + gamma[gj])/(p - 1);
              break;
            }
            tempSc -= Sc2[tempL-1];
            tempL--;
            err = (totalvar - tempSc + tempL*gamma[gj])/(p-tempL);
            temp = Sc2[tempL-1];
          }
          if(Sc2[0]-gamma[gj] < err)
            err = (totalvar)/(p);
        }
        else{
          err = (totalvar)/(p);
        }
        eigenvalue = arma::max(Sc2-(err+gamma[gj])*Sct,Sctz);
        covest =  Phi*Vc2*diagmat(eigenvalue)*trans(Vc2)*trans(Phi);
        output(k,gj) = pow(arma::norm(covvalid - covest - err*Ip,"fro"),2.0);
      }
    }
    
  }
};


using namespace Rcpp;
using namespace arma;
using namespace std;

List spatpcacv_rcpp(NumericMatrix  sxyr, NumericMatrix Yr, int M, int K,  NumericVector  tau1r, NumericVector  tau2r,  NumericVector  nkr, int maxit, double tol, NumericVector  l2r) {
  int n = Yr.nrow(), p = Yr.ncol(), d = sxyr.ncol();
  arma::mat Y(Yr.begin(), n, p, false);
  arma::mat sxy(sxyr.begin(), p, d, false);
  colvec tau1(tau1r.begin(), tau1r.size(),false);
  colvec tau2(tau2r.begin(), tau2r.size(), false);
  colvec nk(nkr.begin(), nkr.size(), false);
  colvec l2(l2r.begin(), l2r.size(), false);
  arma::mat cv(M,tau1.n_elem), out, out2;
  arma::mat Omega;
  double cvtau1, cvtau2;
  arma::mat YYest = Y.t()*Y;
  arma::mat UPhiest, Phiest2;
  vec SPhiest;
  arma::svd_econ(UPhiest, SPhiest, Phiest2, Y, "right");
  double rhoest = 10*pow(SPhiest[0],2.0);
  arma::mat Phiest = Phiest2.cols(0,K-1); 
  arma::mat Cest = Phiest;
  arma::mat  Lambda2est = Phiest*(diagmat(rhoest -1/(rhoest-2*pow(SPhiest.subvec(0,K-1),2))));
  arma::mat rhocv(M,1);
  arma::cube YYtrain(p,p,M), Phicv(p,K,M), Lmbd2cv(p,K,M), tempinvcv(p,p,M);
  if(d == 2)
    Omega = tpmatrix(sxy);
  else
    Omega = cubicmatrix(sxy);
  
  if(tau1.n_elem > 1){  
    spatpcacv_p spatpcacv_p(Y, K, Omega,  tau1, nk,  maxit, tol, cv,YYtrain, Phicv,Lmbd2cv,rhocv);
    RcppParallel::parallelFor(0, M, spatpcacv_p);
    uword  index1;
    (sum(cv,0)).min(index1);
    cvtau1=tau1[index1];  
    if(index1 > 0)
      Phiest = spatpcacore2p(YYest, Cest,  Lambda2est, Omega, cvtau1, rhoest, maxit,tol);
    out = sum(cv,0)/M;
  }
  else{
    cvtau1= max(tau1); 
    Phiest = spatpcacore2p(YYest, Cest, Lambda2est, Omega, cvtau1, rhoest, maxit,tol);
    out.zeros(1);
  }
  
  if(tau2.n_elem > 1){
    
    arma::mat cv2(M,tau2.n_elem);
    uword  index2;

    spatpcacv_p2 spatpcacv_p2(Y, YYtrain, Phicv, Lmbd2cv, rhocv, K, cvtau1, Omega, tau2, nk,  maxit, tol, cv2,tempinvcv); 
    RcppParallel::parallelFor(0, M, spatpcacv_p2);
    //out = join_cols(-sum(cv,0)/M, -sum(cv2,0)/M);
    (sum(cv2,0)).min(index2);
    cvtau2 = tau2[index2];
    mat Ip;
    Ip.eye(Y.n_cols,Y.n_cols);
    mat tempinv = arma::inv_sympd((cvtau1*Omega) - YYest + (rhoest*Ip));
    mat Rest = Phiest;
    mat Lambda1est = 0*Phiest;
     
    for(uword  i = 0; i <= index2; i++){
      spatpcacore3(tempinv, Phiest, Rest, Cest, Lambda1est, Lambda2est, tau2[i], rhoest, maxit,tol);
    }
   
    out2 = sum(cv2,0)/M;
  }
  else{  
    cvtau2 = max(tau2);
    if(cvtau2 > 0){
      mat Ip;
      Ip.eye(Y.n_cols,Y.n_cols);
      mat tempinv = arma::inv_sympd((cvtau1*Omega) - YYest + (rhoest*Ip));
      mat Rest= Phiest;
      mat Lambda1est = 0*Phiest;
      for(uword  i = 0; i < l2.n_elem; i++){
        spatpcacore3(tempinv, Phiest, Rest, Cest, Lambda1est, Lambda2est, l2[i], rhoest, maxit, tol);
      }
    }
    out2.zeros(1);
  }
  
  
  return List::create(Named("cv1") = out, Named("cv2") = out2,Named("est") = Phiest, Named("cvtau1") = cvtau1,Named("cvtau2") = cvtau2);
}

using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::export]]
List spatpcacv2_rcpp(NumericMatrix  sxyr, NumericMatrix Yr, int M, int K,  NumericVector  tau1r, NumericVector  tau2r, NumericVector  gammar,  NumericVector  nkr, int maxit, double tol, NumericVector  l2r) {
  int n = Yr.nrow(), p = Yr.ncol(), d = sxyr.ncol();
  arma::mat Y(Yr.begin(), n, p, false);
  arma::mat sxy(sxyr.begin(), p, d, false);
  colvec tau1(tau1r.begin(), tau1r.size(),false);
  colvec tau2(tau2r.begin(), tau2r.size(), false);
  colvec gamma(gammar.begin(), gammar.size(), false);
  colvec nk(nkr.begin(), nkr.size(), false);
  colvec l2(l2r.begin(), l2r.size(), false);
  arma::mat cv(M,tau1.n_elem), cv3(M,gamma.n_elem), out, out2, out3;
  arma::mat Omega;
  double cvtau1, cvtau2 = 0, cvgamma;
  arma::mat YYest = Y.t()*Y;
  arma::mat UPhiest, Phiest2;
  vec SPhiest;
  arma::svd_econ(UPhiest, SPhiest, Phiest2, Y, "right");
  double rhoest = 10*pow(SPhiest[0],2.0);
  arma::mat Phiest = Phiest2.cols(0,K-1); 
  arma::mat Cest = Phiest;
  arma::mat  Lambda2est = Phiest*(diagmat(rhoest -1/(rhoest-2*pow(SPhiest.subvec(0,K-1),2))));
  arma::mat rhocv(M,1);
  arma::cube YYtrain(p,p,M), Phicv(p,K,M), Lmbd2cv(p,K,M), tempinvcv(p,p,M);
  uword  index1, index2 = 0, index3;
  mat Ip;
  Ip.eye(Y.n_cols,Y.n_cols);
  if(max(tau1) !=0 || max(tau2)!=0){
    if(d == 2)
      Omega = tpmatrix(sxy);
    else
      Omega = cubicmatrix(sxy);
  }
  else{
    if(gamma.n_elem >1){
      Omega = Ip;
      mat Ytrain, Yvalid; 
      arma::mat UPhi, Phioldg;
      vec SPhi;
      for(unsigned k = 0; k < M; ++k){
        Ytrain = Y.rows(arma::find(nk!=(k+1)));
        Yvalid = Y.rows(arma::find(nk==(k+1)));
        arma::svd_econ(UPhi, SPhi, Phioldg, Ytrain, "right");
        Phicv.slice(k) = Phioldg.cols(0,K-1);
        YYtrain.slice(k) = Ytrain.t()*Ytrain;
      }
    }
  }
 
  if(tau1.n_elem > 1){  
    spatpcacv_p spatpcacv_p(Y, K, Omega,  tau1, nk,  maxit, tol, cv,YYtrain, Phicv,Lmbd2cv,rhocv);
    RcppParallel::parallelFor(0, M, spatpcacv_p);
    (sum(cv,0)).min(index1);
    cvtau1=tau1[index1];  
    if(index1 > 0)
      Phiest = spatpcacore2p(YYest, Cest,  Lambda2est, Omega, cvtau1, rhoest, maxit,tol);
    out = sum(cv,0)/M;
  }
  else{
    cvtau1= max(tau1);
    if(cvtau1 != 0 && max(tau2) == 0){
      Phiest = spatpcacore2p(YYest, Cest, Lambda2est, Omega, cvtau1, rhoest, maxit,tol);
      arma::mat  Phigg, Cgg, Lambda2gg;
      mat Yvalid, Ytrain;
      arma::mat UPhi, Phioldg;
      vec SPhi;
      for(unsigned k = 0; k < M; ++k){
        Ytrain = Y.rows(arma::find(nk!=(k+1)));
        YYtrain.slice(k) = Ytrain.t()*Ytrain;
        arma::svd_econ(UPhi, SPhi, Phioldg, Ytrain, "right");
        Phigg = Cgg = Phioldg.cols(0,K-1);
        rhocv(k,0) = 10*pow(SPhi[0],2.0);
        Lambda2gg = (Phigg)*(diagmat(rhocv(k,0) -1/(rhocv(k,0)-2*pow(SPhi.subvec(0,K-1),2))));
        spatpcacore2(YYtrain.slice(k), Phigg,Cgg, Lambda2gg, Omega, cvtau1, rhocv(k,0), maxit,tol);
        Phicv.slice(k) = Phigg;
        Lmbd2cv.slice(k) = Lambda2gg;
        tempinvcv.slice(k) = arma::inv_sympd((cvtau1*Omega) - YYtrain.slice(k) + (rhocv(k,0)*Ip));    
      }
    }
    else if(cvtau1 == 0 && max(tau2) == 0){
      arma::mat UPhi2, Phioldc;
      vec SPhi2;
      arma::svd_econ(UPhi2, SPhi2, Phioldc, Y, "right");
      Phiest = Phioldc.cols(0,K-1); 
    }
    else{
      arma::mat  Phigg, Cgg, Lambda2gg;
      mat Yvalid, Ytrain;
      arma::mat UPhi, Phioldg;
      vec SPhi;
      for(unsigned k = 0; k < M; ++k){
        Ytrain = Y.rows(arma::find(nk!=(k+1)));
        YYtrain.slice(k) = Ytrain.t()*Ytrain;
        arma::svd_econ(UPhi, SPhi, Phioldg, Ytrain, "right");
        Phigg = Cgg = Phioldg.cols(0,K-1);
        rhocv(k,0) = 10*pow(SPhi[0],2.0);
        Lambda2gg = (Phigg)*(diagmat(rhocv(k,0) -1/(rhocv(k,0)-2*pow(SPhi.subvec(0,K-1),2))));
        if(cvtau1 != 0)
          spatpcacore2(YYtrain.slice(k), Phigg,Cgg, Lambda2gg, Omega, cvtau1, rhocv(k,0), maxit,tol);
        Phicv.slice(k) = Phigg;
        Lmbd2cv.slice(k) = Lambda2gg;
       // tempinvcv.slice(k) = arma::inv_sympd((cvtau1*Omega) - YYtrain.slice(k) + (rhocv(k,0)*Ip));    
      }
    }
      
    
    out.zeros(1);
  }
  if(tau2.n_elem > 1){
    
    arma::mat cv2(M,tau2.n_elem);

    spatpcacv_p2 spatpcacv_p2(Y, YYtrain, Phicv, Lmbd2cv, rhocv, K, cvtau1, Omega, tau2, nk,  maxit, tol, cv2,tempinvcv); 
    RcppParallel::parallelFor(0, M, spatpcacv_p2);
  
    //out = join_cols(-sum(cv,0)/M, -sum(cv2,0)/M);
    (sum(cv2,0)).min(index2);
    cvtau2 = tau2[index2];
    

    mat tempinv = arma::inv_sympd((cvtau1*Omega) - YYest + (rhoest*Ip));
    mat Rest = Phiest;
    mat Lambda1est = 0*Phiest;
     
    for(uword  i = 0; i <= index2; i++){
      spatpcacore3(tempinv, Phiest, Rest, Cest, Lambda1est, Lambda2est, tau2[i], rhoest, maxit,tol);
    }
   
    out2 = sum(cv2,0)/M;
  }
  else{  
    cvtau2 = max(tau2);
    if(cvtau2 > 0){
      mat tempinv = arma::inv_sympd((cvtau1*Omega) - YYest + (rhoest*Ip));
      mat Rest= Phiest;
      mat Lambda1est = 0*Phiest;
      for(uword  i = 0; i < l2.n_elem; i++){
        spatpcacore3(tempinv, Phiest, Rest, Cest, Lambda1est, Lambda2est, l2[i], rhoest, maxit, tol);
      }
    }
    out2.zeros(1);
  }
 
  spatpcacv_p3 spatpcacv_p3(Y,YYtrain, Phicv, Lmbd2cv, rhocv, tempinvcv, index2,K, Omega,  cvtau1, tau2, gamma, nk,  maxit, tol, cv3);
  RcppParallel::parallelFor(0, M, spatpcacv_p3);
  if(gamma.n_elem > 1){
    (sum(cv3,0)).min(index3);
    cvgamma = gamma[index3];
  }
  else{
    cvgamma = max(gamma);
  }
  out3 = sum(cv3,0)/M;

  return List::create(Named("cv1") = out, Named("cv2") = out2, Named("cv3") = out3, Named("est") = Phiest, Named("cvtau1") = cvtau1, Named("cvtau2") = cvtau2,Named("cvgamma") = cvgamma);
}



// [[Rcpp::export]]
List eigenest_rcpp(NumericMatrix  phir, NumericMatrix Yr, double gamma, NumericMatrix  phi2r) {
  int n = Yr.nrow(), p = phir.nrow(), K = phir.ncol(),p2 = phi2r.nrow() ;
  arma::mat phi(phir.begin(), p, K, false);
  arma::mat phi2(phi2r.begin(), p2, K, false);
  arma::mat Y(Yr.begin(), n, p, false);
  arma::mat Vc, Vc2;
  arma::vec Sc, Sc2;
  
  arma::mat cov = Y.t()*Y/n;
  arma::eig_sym(Sc,Vc, phi.t()*cov*phi);
  int tempL = K;
  double totalvar = trace(cov);
  double tempSc2 = accu(Sc);
  double temp_v = totalvar - tempSc2;
  double err = (temp_v + K*gamma)/(p-tempL);

  Sc2 = sort(Sc,"descend");
  Vc2 = Vc.cols(sort_index(Sc,"descend"));
  double tempSc = tempSc2, temp;
  arma::mat Sct, Sctz, Phi;
  Sct.ones(Sc2.n_elem);
  Sctz.zeros(Sc2.n_elem);
  
  if(Sc2[0] > gamma){
    err = (totalvar - tempSc + K*gamma)/(p - tempL);
    temp = Sc2[tempL-1];  
    while(temp-gamma < err){
      if(tempL == 1){
        err = (totalvar - Sc2[0] + gamma)/(p-1);
        break;
      }
      tempSc -= Sc2[tempL-1];
      tempL--;
      err = (totalvar - tempSc + tempL*gamma)/(p-tempL);
      temp = Sc2[tempL-1];
    }
    if(Sc2[0]-gamma < err)
      err = (totalvar)/(p);
  }
  else{
    err = (totalvar)/(p);
  }
  
  arma::vec eigenvalue = arma::max(Sc2-(err+gamma)*Sct,Sctz);
  arma::vec eigenvalue2 = eigenvalue+err;
  arma::mat covest =  phi*Vc2*diagmat(eigenvalue/eigenvalue2)*trans(phi2*Vc2);
  arma::mat predict = Y*covest;
  return List::create(Named("err") = err,Named("Phi") = Phi, Named("eigenvalue") = eigenvalue,Named("covest") = covest, Named("predict") = predict);
}


arma::mat eigenest2(mat phi, mat cov, double gamma) {
  int p = phi.n_rows, K = phi.n_cols;
  arma::mat Vc, Vc2;
  arma::vec Sc, Sc2;
  arma::eig_sym(Sc,Vc, phi.t()*cov*phi);
  int tempL = K;
  double err;
  arma::vec eigenvalue;
  arma::mat covest;
  double totalvar = trace(cov);
  double tempSc2 = accu(Sc);
  double temp_v = totalvar - tempSc2;
  err = (temp_v + K*gamma)/(p-tempL);
  
  Sc2 = sort(Sc,"descend");
  Vc2 = Vc.cols(sort_index(Sc,"descend"));
  double tempSc = tempSc2, temp;
  arma::mat Sct, Sctz, Phi, Ip;
  Sct.ones(Sc2.n_elem);
  Sctz.zeros(Sc2.n_elem);
  Ip.eye(p,p);
  if(Sc2[0] > gamma){
    err = (totalvar - tempSc + K*gamma)/(p - tempL);
    temp = Sc2[tempL-1];  
    while(temp-gamma < err){
      if(tempL == 1){
        err = (totalvar - Sc2[0] + gamma)/(p-1);
        break;
      }
      tempSc -= Sc2[tempL-1];
      tempL--;
      err = (totalvar - tempSc + tempL*gamma)/(p-tempL);
      temp = Sc2[tempL-1];
    }
    if(Sc2[0] - gamma < err)
      err = (totalvar)/(p);
  }
  else{
    err = (totalvar)/(p);
  }
  
  eigenvalue = arma::max(Sc2 - (err+gamma)*Sct,Sctz);
  Phi = phi*Vc2;
  covest =  Phi*diagmat(eigenvalue)*trans(Phi) + err*Ip;
  return(covest);
}


using namespace arma;
using namespace Rcpp;
using namespace std;

arma::mat spatpca_rcpp(const arma::mat  sxy, const arma::mat Y, const int K, const double l1, const arma::vec l2,  const int maxit, const double tol){
  
  int d = sxy.n_cols;
  double rho;
  arma::mat UPhi, Phiold, Phi, R, C, Lambda1, Lambda2;
  arma::vec SPhi;
  arma::mat Omega;
  if(d ==2)
    Omega = tpmatrix(sxy);
  else
    Omega = cubicmatrix(sxy);
  arma::svd_econ(UPhi, SPhi, Phiold, Y,"right");
  Phi = Phiold.cols(0,K-1);
  C = Phi;
  rho = 10*pow(SPhi[0],2.0);
  
  Lambda2= Phi*(diagmat(rho -1/(rho-2*pow(SPhi.subvec(0,K-1),2)))); 
  mat YY = Y.t()*Y;
  
  if(l1 > 0)
    Phi = spatpcacore2p(YY,C,Lambda2, Omega, l1, rho, maxit,tol);
  
  
  if(max(l2)!=0){
    mat Ip;
    Ip.eye(Y.n_cols,Y.n_cols);
    mat tempinv = arma::inv_sympd((l1*Omega) - YY + (rho*Ip));
    Lambda1 = 0*Phi;
    R = Phi;
    for (unsigned int j = 0; j < l2.n_elem; j++){
      spatpcacore3(tempinv,Phi,R,C,Lambda1,Lambda2, l2[j], rho, maxit,tol);
    }
  }
  return Phi;
}

using namespace arma;
using namespace Rcpp;
using namespace std;

arma::mat spatpca2_rcpp(const arma::mat Omega, const arma::mat Y, const int K, const double l1, const arma::vec l2,  const int maxit, const double tol){
  
 
  double rho;
  arma::mat UPhi, Phiold, Phi, R, C, Lambda1, Lambda2;
  arma::vec SPhi;
  
  arma::svd_econ(UPhi, SPhi, Phiold, Y,"right");
  Phi = Phiold.cols(0,K-1);
  C = Phi;
  rho = 10*pow(SPhi[0],2.0);
  
  Lambda2= Phi*(diagmat(rho -1/(rho-2*pow(SPhi.subvec(0,K-1),2)))); 
  mat YY = Y.t()*Y;
  
  if(l1 > 0)
    Phi = spatpcacore2p(YY,C,Lambda2, Omega, l1, rho, maxit,tol);
  
  
  if(max(l2)!=0){
    mat Ip;
    Ip.eye(Y.n_cols,Y.n_cols);
    mat tempinv = arma::inv_sympd((l1*Omega) - YY + (rho*Ip));
    Lambda1 = 0*Phi;
    R = Phi;
    for (unsigned int j = 0; j < l2.n_elem; j++){
      spatpcacore3(tempinv,Phi,R,C,Lambda1,Lambda2, l2[j], rho, maxit,tol);
    }
  }
  return Phi;
}


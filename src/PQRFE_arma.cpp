#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::vec fill_v(int a, int b, arma::vec x){
  int n = b - a + 1;
  arma::vec y(n);
  for(int i=0; i<n; i++){
    y(i) = x(i+a-1);
  }
  return(y);
}

//' Rho Koenker
//'
//' @param x generic vector
//' @param tau percentile
//' 
//' @return y vector, linear transformation by rho 
// [[Rcpp::export(rho_koenker)]]
arma::vec rho_koenker(arma::vec x, double tau){
  int n = x.n_elem;
  arma::vec y(n); 
  for(int i = 0; i < n; ++i){
    if(x(i)<0){
      y(i) = x(i)*(tau-1);
    } else {
      y(i) = x(i)*tau;
    }
  }  
  return(y);
}

//' Rho M-quantile
//'
//' @param x generic vector
//' @param tau percentile
//' @param c tuning
//' 
//' @return y vector, linear transformation by m-rho
// [[Rcpp::export(rho_mq)]]
arma::vec rho_mq(arma::vec x, double tau, double c){
  int n = x.n_elem;
  arma::vec y(n);
  arma::vec a(n);
  arma::vec y1(n);
  arma::vec y2(n);
  arma::vec b1(n);
  arma::vec b2(n);
  a = abs(x);
  for(int i = 0; i < n; ++i){
    if(a(i) > c){
      b1(i) = 1;
      b2(i) = 0;
    } else {
      b1(i) = 0;
      b2(i) = 1;
    }
    if(x(i) < 0){
      y1(i) = 1-tau;
      y2(i) = 1-tau;
    } else {
      y1(i) = tau;
      y2(i) = tau;
    }
    y1(i) = y1(i) * (c*a(i) - 0.5*pow(c,2));
    y2(i) = y2(i) * (0.5*pow(a(i),2));
    y(i) = (b1(i) * y1(i)) + (b2(i) * y2(i));
  }
  return(y);
}

//' Psi M-quantile
//'
//' @param x generic vector
//' @param tau percentile
//' @param c tuning
//' 
//' @return y vector, linear transformation by m-rho derivative
// [[Rcpp::export(psi_mq)]]
arma::vec psi_mq(arma::vec x, double tau, double c){
  int n = x.n_elem;
  arma::vec y(n);
  arma::vec a(n);
  arma::vec y1(n);
  arma::vec y2(n);
  arma::vec b1(n);
  arma::vec b2(n);
  a = abs(x);
  for(int i = 0; i < n; ++i){
    if(a(i) > c){
      b1(i) = 1;
      b2(i) = 0;
    } else {
      b1(i) = 0;
      b2(i) = 1;
    }
    if(x(i) < 0){
      y1(i) = 1-tau;
      y2(i) = 1-tau;
    } else {
      y1(i) = tau;
      y2(i) = tau;
    }
    y1(i) = y1(i) * c * x(i)/a(i);
    y2(i) = y2(i) * x(i);
    y(i) = (b1(i) * y1(i)) + (b2(i) * y2(i));
  }
  return(y);
}

//' D Psi M-quantile
//'
//' @description Derivative of psi M-quantile
//'
//' @param x generic vector
//' @param tau percentile
//' @param c tuning
//' 
//' @return y vector, linear transformation by second derivative m-rho
// [[Rcpp::export(d_psi_mq)]]
arma::vec d_psi_mq(arma::vec x, double tau, double c){
  int n = x.n_elem;
  arma::vec y(n);
  arma::vec a(n);
  arma::vec y1(n);
  arma::vec y2(n);
  arma::vec b1(n);
  arma::vec b2(n);
  a = abs(x);
  for(int i = 0; i < n; ++i){
    if(a(i) > c){
      b1(i) = 1;
      b2(i) = 0;
    } else {
      b1(i) = 0;
      b2(i) = 1;
    }
    if(x(i) < 0){
      y1(i) = 1-tau;
      y2(i) = 1-tau;
    } else {
      y1(i) = tau;
      y2(i) = tau;
    }
    y1(i) = y1(i) * 0;
    y2(i) = y2(i) * 1;
    y(i) = (b1(i) * y1(i)) + (b2(i) * y2(i));
  }
  return(y);
}

arma::vec rho_als(arma::vec x, double tau){
  int n = x.n_elem;
  arma::vec y(n); 
  for(int i = 0; i < n; ++i){
    if(x(i)<0){
      y(i) = pow(x(i),2)*(1-tau);
    } else {
      y(i) = pow(x(i),2)*tau;
    }
  }  
  return(y);
}

//' Psi ALS
//'
//' @description Psi asymetric least square
//'
//' @param x generic vector
//' @param tau percentile
//' 
//' @return y vector, linear transformation by ALS psi
// [[Rcpp::export(psi_als)]]
arma::vec psi_als(arma::vec x, double tau){
  int n = x.n_elem;
  arma::vec y(n); 
  for(int i = 0; i < n; ++i){
    if(x(i)<0){
      y(i) = 2*x(i)*(1-tau);
    } else {
      y(i) = 2*x(i)*tau;
    }
  }  
  return(y);
}

//' D Psi ALS
//'
//' @description Derivative of Psi asymetric least square
//'
//' @param x generic vector
//' @param tau percentile
//' 
//' @return y vector, linear transformation by derivative ALS psi
// [[Rcpp::export(d_psi_als)]]
arma::vec d_psi_als(arma::vec x, double tau){
  int n = x.n_elem;
  arma::vec y(n); 
  for(int i = 0; i < n; ++i){
    if(x(i)<0){
      y(i) = 2*(1-tau);
    } else {
      y(i) = 2*tau;
    }
  }  
  return(y);
}

//' Loss quantile regression
//' 
//' @description This function returns the core of quantile regression to be minimized
//'
//' @param beta initial values
//' @param x design matrix
//' @param y vector output
//' @param tau percentile
//' @param N sample size
//' @param d columns of x  
//' 
//' @return eta Numeric, sum of quantile regression
// [[Rcpp::export(loss_qr)]]
double loss_qr(arma::vec beta, arma::mat x, arma::vec y, double tau, int N, int d){
  double eta = 0;
  arma::vec res(N);
  arma::vec rho(N);
  res = y - (x * beta);
  rho = rho_koenker(res,tau);
  eta = accu(rho);
  return(log(eta));
}

//' Loss quantile regression with fixed effects
//' 
//' @description This function returns the core of quantile regression with fixed effects to be minimized
//'
//' @param theta initial values
//' @param x design matrix
//' @param y vector output
//' @param z incident matrix
//' @param tau percentile
//' @param n N sample size
//' @param d columns of x
//' @param mm n columns of z   
//' 
//' @return eta Numeric, sum of quantile regression with fixed effects
// [[Rcpp::export(loss_qrfe)]]
double loss_qrfe(arma::vec theta,arma::mat x,arma::vec y,arma::mat z,double tau,int n,int d,int mm){
  double eta;
  arma::vec beta(d);
  arma::vec alpha(mm);
  arma::vec res(n);
  arma::vec rho(n);
  beta = fill_v(1,d, theta); 
  alpha = fill_v((d+1),d+mm, theta); 
  res = y -  (z * alpha) -  (x * beta);
  rho = rho_koenker(res,tau);
  eta = accu(rho);
  return(log(eta));
}

//' Loss lasso quantile regression with fixed effects
//'
//' @description This function returns the core of lasso quantile regression with fixed effects to be minimized
//'
//' @param theta initial values
//' @param x design matrix
//' @param y vector output
//' @param z incident matrix
//' @param tau percentile
//' @param n N sample size
//' @param d columns of x
//' @param mm n columns of z  
//' @param lambda constriction parameter
//' 
//' @return eta Numeric, sum of lasso quantile regression with fixed effects
// [[Rcpp::export(loss_qrlasso)]]
double loss_qrlasso(arma::vec theta,arma::mat x,arma::vec y,arma::mat z,double tau,int n,int d,int mm, double lambda){
  double eta;
  arma::vec beta(d);
  arma::vec alpha(mm);
  arma::vec res(n);
  arma::vec rho(n);
  beta = fill_v(1,d, theta); 
  alpha = fill_v(d+1,d+mm, theta); 
  res = y -  (z * alpha) -  (x * beta);   
  rho = rho_koenker(res,tau);
  eta = accu(rho)/n + lambda * accu(abs(alpha));
  return(log(eta));
}

//' Loss M-quantile regression
//'
//' @description This function returns the core of M-quantile regression to be minimized
//'
//' @param beta initial values
//' @param x design matrix
//' @param y vector output
//' @param tau percentile
//' @param N sample size
//' @param d columns of x  
//' @param c tuning
//' 
//' @return eta Numeric, sum of M-quantile regression
// [[Rcpp::export(loss_mqr)]]
double loss_mqr(arma::vec beta, arma::mat x, arma::vec y, double tau, int N, int d, double c){
  double eta = 0;
  arma::vec res(N);
  arma::vec rho(N);
  res = y -  (x * beta);
  rho = rho_mq(res,tau,c);
  eta = accu(rho);
  return(log(eta));
}

//' Loss M-quantile regression with fixed effects
//'
//' @description This function returns the core of M-quantile regression with fixed effects to be minimized
//'
//' @param theta initial values
//' @param x design matrix
//' @param y vector output
//' @param z incident matrix
//' @param tau percentile
//' @param n N sample size
//' @param d columns of x
//' @param mm n columns of z 
//' @param c tuning
//' 
//' @return eta Numeric, sum of M-quantile regression with fixed effects
// [[Rcpp::export(loss_mqrfe)]]
double loss_mqrfe(arma::vec theta,arma::mat x,arma::vec y,arma::mat z,double tau,int n,int d,int mm, double c){
  double eta;
  arma::vec beta(d);
  arma::vec alpha(mm);
  arma::vec res(n);
  arma::vec rho(n);
  beta = fill_v(1,d, theta); 
  alpha = fill_v((d+1),d+mm, theta); 
  res = y -  (z * alpha) -  (x * beta);
  rho = rho_mq(res,tau,c);
  eta = accu(rho);
  return(log(eta));
}

//' Loss lasso M-quantile regression with fixed effects
//' 
//' @description This function returns the core of lasso M-quantile regression with fixed effects to be minimized   
//'
//' @param theta initial values
//' @param x design matrix
//' @param y vector output
//' @param z incident matrix
//' @param tau percentile
//' @param n N sample size
//' @param d columns of x
//' @param mm n columns of z  
//' @param c tuning
//' @param lambda constriction parameter
//' 
//' @return eta Numeric, sum of lasso M-quantile regression with fixed effects
// [[Rcpp::export(loss_mqrlasso)]]
double loss_mqrlasso(arma::vec theta,arma::mat x,arma::vec y,arma::mat z,double tau,int n,int d,int mm, double c, double lambda){
  double eta;
  arma::vec beta(d);
  arma::vec alpha(mm);
  arma::vec res(n);
  arma::vec rho(n);
  beta = fill_v(1,d, theta); 
  alpha = fill_v(d+1,d+mm, theta); 
  res = y -  (z * alpha) -  (x * beta);   
  rho = rho_mq(res,tau,c);
  eta = accu(rho)/n + lambda * accu(abs(alpha));
  return(log(eta));
}

//' Loss expectile regression
//'
//' @description This function returns the core of expectile regression to be minimized   
//' 
//' @param beta initial values
//' @param x design matrix
//' @param y vector output
//' @param tau percentile
//' @param N sample size
//' @param d columns of x  
//' 
//' @return eta Numeric, sum of expectile regression
// [[Rcpp::export(loss_er)]]
double loss_er(arma::vec beta, arma::mat x, arma::vec y, double tau, int N, int d){
  double eta = 0;
  arma::vec res(N);
  arma::vec rho(N);
  res = y -  (x * beta);
  rho = rho_als(res,tau);
  eta = accu(rho);
  return(log(eta));
}

//' Loss expectile regression with fixed effects
//'
//' @description This function returns the core of expectile regression with fixed effects to be minimized   
//'
//' @param theta initial values
//' @param x design matrix
//' @param y vector output
//' @param z incident matrix
//' @param tau percentile
//' @param n N sample size
//' @param d columns of x
//' @param mm n columns of z  
//' 
//' @return eta Numeric, sum of expectile regression with fixed effects
// [[Rcpp::export(loss_erfe)]]
double loss_erfe(arma::vec theta,arma::mat x,arma::vec y,arma::mat z,double tau,int n,int d,int mm){
  double eta;
  arma::vec beta(d);
  arma::vec alpha(mm);
  arma::vec res(n);
  arma::vec rho(n);
  beta = fill_v(1,d, theta); 
  alpha = fill_v((d+1),d+mm, theta); 
  res = y -  (z * alpha) -  (x * beta);
  rho = rho_als(res,tau);
  eta = accu(rho);
  return(log(eta));
}

//' Loss lasso expectile regression with fixed effects
//' 
//' @description This function returns the core of lasso expectile regression with fixed effects to be minimized   
//'  
//' @param theta initial values
//' @param x design matrix
//' @param y vector output
//' @param z incident matrix
//' @param tau percentile
//' @param n N sample size
//' @param d columns of x
//' @param mm n columns of z  
//' @param lambda constriction parameter
//' 
//' @return eta Numeric, sum of lasso expectile regression with fixed effects
// [[Rcpp::export(loss_erlasso)]]
double loss_erlasso(arma::vec theta,arma::mat x,arma::vec y,arma::mat z,double tau,int n,int d,int mm, double lambda){
  double eta;
  arma::vec beta(d);
  arma::vec alpha(mm);
  arma::vec res(n);
  arma::vec rho(n);
  beta = fill_v(1,d, theta); 
  alpha = fill_v(d+1,d+mm, theta); 
  res = y -  (z * alpha) -  (x * beta);   
  rho = rho_als(res,tau);
  eta = accu(rho)/n + lambda * accu(abs(alpha));
  return(log(eta));
}


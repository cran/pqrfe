#' @title Penalized quantile regression with fixed effects 
#'
#' @description Estimate parameters and tuning parameter.
#'
#' @param y      Numeric vector, outcome.
#' @param x      Numeric matrix, covariates
#' @param subj   Numeric vector, identifies the unit to which the observation belongs.
#' @param tau    Numeric scalar between zero and one, identifies the percentile.
#' @param effect Factor, "simple" simple regression, "fixed" regression with fixed effects, "lasso" penalized regression with fixed effects.
#' @param c  Numeric, 0 is quantile, Inf is expectile, any number between zero and infinite is M-quantile.
#  
#' @return alpha       Numeric vector, intercepts' coefficients.
#' @return beta        Numeric vector, exploratory variables' coefficients.
#' @return lambda      Numeric, estimated lambda.
#' @return res         Numeric vector, percentile residuals.
#' @return tau         Numeric scalar, the percentile.
#' @return penalty     Numeric scalar, indicate the chosen effect.
#' @return c           Numeric scalar, indicate the chosen c.
#' @return sig2_alpha  Numeric vector, intercepts' standard errors.
#' @return sig2_beta   Numeric vector, exploratory variables' standard errors.
#' @return Tab_alpha   Data.frame, intercepts' summary.
#' @return Tab_beta    Data.frame, exploratory variables' summary.
#' @return Mat_alpha   Numeric matrix, intercepts' summary.
#' @return Mat_beta    Numeric matrix, exploratory variables' summary.
#' 
#' @importFrom MASS stats graphics
#' @import Rcpp RcppArmadillo
#' 
#' @examples 
#' n = 10
#' m = 5
#' d = 4
#' N = n*m
#' x = matrix(rnorm(d*N), ncol=d, nrow=N)
#' subj = rep(1:n, each=m)
#' alpha = rnorm(n)
#' beta = rnorm(d)
#' eps = rnorm(N)
#' y = as.vector(x %*% beta + rep(alpha, each=m) + eps)
#' m1 = pqr(x=x, y=y, subj=subj, tau=0.75, effect="lasso", c = 0)
#' m1$Tab_beta
#' 
#' @references
#' Koenker, R. (2004) "Quantile regression for longitudinal data", J. Multivar. Anal., 91(1): 74-89, <doi:10.1016/j.jmva.2004.05.006>
#'      
#' @export
pqr =  function(x,y,subj,tau=0.5, effect="simple", c=1){
  penalty = choice_p(effect) 
  d  = ifelse(is.null(dim(x)[2]), 1, dim(x)[2])
  if(d==1) x = as.matrix(x)
  dclean = clean_data(y, x, id=subj)
  subj = as.vector(dclean$id)
  y = as.vector(dclean$y)
  N  = length(y)
  x = as.matrix(dclean$x)
  beta = as.vector(solve(t(x)%*%x)%*%t(x)%*%y)  
  if(c==0)         beta = optim_qr(beta,x,y,tau,N,d)$beta 
  if(c>0 && c<Inf) beta = optim_mqr(beta,x,y,tau,N,d,c)$beta          
  if(c==Inf)       beta = optim_er(beta,x,y,tau,N,d)$beta
  n = max(subj)
  alpha = rep(0,n)
  z = matrix(0, nrow = N, ncol = n)
  j = 0 
  for(i in 1:N){
    j = subj[i] 
    z[i,j] = 1
  }
  if(penalty==2){
    if(c==0){
      opt = optim_qrfe(beta,alpha,x,y,z,tau,N,d,n)    
    }
    if(c>0 && c<Inf){
      opt = optim_mqrfe(beta,alpha,x,y,z,tau,N,d,n,c)          
    }
    if(c==Inf){
      opt = optim_erfe(beta,alpha,x,y,z,tau,N,d,n)          
    }
  }
  if(penalty==3){
    if(c==0){
      opt = optim_qrlasso(beta,alpha,x,y,z,tau,N,d,n)    
    }
    if(c>0 && c<Inf){
      opt = optim_mqrlasso(beta,alpha,x,y,z,tau,N,d,n,c)          
    }
    if(c==Inf){
      opt = optim_erlasso(beta,alpha,x,y,z,tau,N,d,n)          
    }
  }  
  if(penalty==1){
    if(c==0){
      opt = optim_qr(beta,x,y,tau,N,d)    
    }
    if(c>0 && c<Inf){
      opt = optim_mqr(beta,x,y,tau,N,d,c)          
    }
    if(c==Inf){
      opt = optim_er(beta,x,y,tau,N,d)          
    }
  }
  alpha  = opt$alpha
  beta   = opt$beta
  lambda = opt$lambda
  res    = opt$res
  Sigma = q_cov(n, N, d, Z=z, X=x, tau, res, penalty,c)
  sig2_alpha = Sigma$sig2_alpha
  sig2_beta = Sigma$sig2_beta
  Tab_alpha = f_tab(N, n, d, alpha, sig2_alpha, 1)$Core
  Tab_beta  = f_tab(N, n, d, beta,  sig2_beta,  2)$Core
  Mat_beta  = f_tab(N, n, d, beta,  sig2_beta,  2)$Matx
  obj = list(alpha=alpha, beta=beta, lambda=lambda, res=res,tau=tau, penalty=penalty, c=c, sig2_alpha=sig2_alpha,sig2_beta=sig2_beta, Tab_alpha=Tab_alpha, Tab_beta=Tab_beta, Mat_beta=Mat_beta)
  class(obj) = "PQR"
  return(obj)
}

#' @title optim quantile regression
#'
#' @description This function solves a quantile regression
#' 
#' @param beta Numeric vector, initials values.
#' @param x Numeric matrix, covariates.
#' @param y Numeric vector, output.
#' @param tau Numeric scalar, the percentile.
#' @param N Numeric integer, sample size.
#' @param d Numeric integer, X number of columns.
#' 
#' @return parametric vector and residuals.
optim_qr = function(beta,x,y,tau,N,d){
  Opt    = stats::optim(par=beta, fn = loss_qr, method = "L-BFGS-B",  
                 N=N, y=y, x=x, tau=tau, d=d)
  beta   = Opt$par
  if(d==1) res = y -  (x * beta)
  else     res = y -  (x %*% beta)   
  return(list(alpha=0, beta=beta, lambda=0, res=res))
}

#' @title optim quantile regression with fixed effects
#'
#' @description This function solves a quantile regression  with fixed effects
#' 
#' @param beta Numeric vector, initials values beta.
#' @param alpha Numeric vector, initials values alpha.
#' @param x Numeric matrix, covariates.
#' @param y Numeric vector, output.
#' @param z Numeric matrix, incidence matrix.
#' @param tau Numeric scalar, the percentile.
#' @param N Numeric integer, sample size.
#' @param d Numeric integer, X number of columns.
#' @param n Numeric integer, length of alpha.
#' 
#' @return parametric vector and residuals.
optim_qrfe = function(beta,alpha,x,y,z,tau,N,d,n){
  M      = n
  theta  = c(beta,alpha)
  Opt    = stats::optim(par=theta, fn = loss_qrfe, method = "L-BFGS-B",  
            n=N, y=y, x=x, z=z, tau=tau, d=d, mm=n)
  beta   = Opt$par[1:d]
  alpha  = Opt$par[(d+1):(d+n)]
  if(d==1) res = y -  (z %*% alpha) -  (x * beta)
  else     res = y -  (z %*% alpha) -  (x %*% beta) 
  sdu    = stats::sd(res, na.rm = TRUE)
  sda    = stats::sd(alpha, na.rm = TRUE)
  lambda = check_lambda(sdu/sda, 1/M, M)
  return(list(alpha=alpha, beta=beta, lambda=lambda, res=res))
}

#' @title optim quantile regression with fixed effects and LASSO
#'
#' @description This function solves a quantile regression  with fixed effects and LASSO
#' 
#' @param beta Numeric vector, initials values beta.
#' @param alpha Numeric vector, initials values alpha.
#' @param x Numeric matrix, covariates.
#' @param y Numeric vector, output.
#' @param z Numeric matrix, incidence matrix.
#' @param tau Numeric scalar, the percentile.
#' @param N Numeric integer, sample size.
#' @param d Numeric integer, X number of columns.
#' @param n Numeric integer, length of alpha.
#' 
#' @return parametric vector and residuals.
optim_qrlasso = function(beta,alpha,x,y,z,tau,N,d,n){
  M       = n
  theta  = c(beta,alpha)
  Opt    = stats::optim(par=theta, fn = loss_qrlasso, method = "L-BFGS-B",   
                 n=N, y=y, x=x, z=z, tau=tau, d=d, mm=n,lambda=1)
  beta   = Opt$par[1:d]
  alpha  = Opt$par[(d+1):(d+n)]
  if(d==1) res = y -  (z %*% alpha) -  (x * beta)
  else     res = y -  (z %*% alpha) -  (x %*% beta)   
  sdu    = stats::sd(res, na.rm = TRUE)
  sda    = stats::sd(alpha, na.rm = TRUE)
  lambda = check_lambda(sdu/sda, 1/M, M)
  return(list(alpha=alpha, beta=beta, lambda=lambda, res=res))
}

#' @title optim M-quantile regression
#'
#' @description This function solves a M-quantile regression
#' 
#' @param beta Numeric vector, initials values beta.
#' @param x Numeric matrix, covariates.
#' @param y Numeric vector, output.
#' @param tau Numeric scalar, the percentile.
#' @param N Numeric integer, sample size.
#' @param d Numeric integer, X number of columns.
#' @param c Numeric, positive real value.
#' 
#' @return parametric vector and residuals.
optim_mqr = function(beta,x,y,tau,N,d,c){
  Opt    = stats::optim(par=beta, fn = loss_mqr, method = "L-BFGS-B",  
                 N=N, y=y, x=x, tau=tau, d=d, c=c)
  beta   = Opt$par[1:d]
  if(d==1) res = y -  (x * beta)
  else     res = y -  (x %*% beta)   
  return(list(alpha=0, beta=beta, lambda=0, res=res))
}

#' @title optim quantile regression with fixed effects 
#'
#' @description This function solves a quantile regression  with fixed effects
#' 
#' @param beta Numeric vector, initials values beta.
#' @param alpha Numeric vector, initials values alpha.
#' @param x Numeric matrix, covariates.
#' @param y Numeric vector, output.
#' @param z Numeric matrix, incidence matrix.
#' @param tau Numeric scalar, the percentile.
#' @param N Numeric integer, sample size.
#' @param d Numeric integer, X number of columns.
#' @param n Numeric integer, length of alpha.
#' @param c Numeric, positive real value.
#' 
#' @return parametric vector and residuals.
optim_mqrfe = function(beta,alpha,x,y,z,tau,N,d,n,c){
  M      = n
  theta  = c(beta,alpha)
  Opt    = stats::optim(par=theta, fn = loss_mqrfe, method = "L-BFGS-B", 
                 n=N, y=y, x=x, z=z, tau=tau, d=d, mm=n, c=c)
  beta   = Opt$par[1:d]
  alpha  = Opt$par[(d+1):(d+n)]
  if(d==1) res = y -  (z %*% alpha) -  (x * beta)
  else     res = y -  (z %*% alpha) -  (x %*% beta) 
  sdu    = stats::sd(res, na.rm = TRUE)
  sda    = stats::sd(alpha, na.rm = TRUE)
  lambda = check_lambda(sdu/sda, 1/M, M)
  return(list(alpha=alpha, beta=beta, lambda=lambda, res=res))
}

#' @title optim M-quantile regression with fixed effects and LASSO
#'
#' @description This function solves a M-quantile regression  with fixed effects and LASSO
#' 
#' @param beta Numeric vector, initials values beta.
#' @param alpha Numeric vector, initials values alpha.
#' @param x Numeric matrix, covariates.
#' @param y Numeric vector, output.
#' @param z Numeric matrix, incidence matrix.
#' @param tau Numeric scalar, the percentile.
#' @param N Numeric integer, sample size.
#' @param d Numeric integer, X number of columns.
#' @param n Numeric integer, length of alpha.
#' @param c Numeric, positive real value.
#' 
#' @return parametric vector and residuals.
optim_mqrlasso = function(beta,alpha,x,y,z,tau,N,d,n,c){
  M      = n
  theta  = c(beta,alpha)
  Opt    = stats::optim(par=theta, fn = loss_mqrlasso, method = "L-BFGS-B",   
                 n=N, y=y, x=x, z=z, tau=tau, d=d, mm=n, c=c,lambda=1)
  beta   = Opt$par[1:d]
  alpha  = Opt$par[(d+1):(d+n)]
  if(d==1) res = y -  (z %*% alpha) -  (x * beta)
  else     res = y -  (z %*% alpha) -  (x %*% beta)   
  sdu    = stats::sd(res, na.rm = TRUE)
  sda    = stats::sd(alpha, na.rm = TRUE)
  lambda = check_lambda(sdu/sda, 1/M, M)
  return(list(alpha=alpha, beta=beta, lambda=lambda, res=res))
}

#' @title optim expectile regression
#'
#' @description This function solves a expectile regression
#' 
#' @param beta Numeric vector, initials values beta.
#' @param x Numeric matrix, covariates.
#' @param y Numeric vector, output.
#' @param tau Numeric scalar, the percentile.
#' @param N Numeric integer, sample size.
#' @param d Numeric integer, X number of columns.
#' 
#' @return parametric vector and residuals.
optim_er = function(beta,x,y,tau,N,d){
  Opt    = stats::optim(par=beta, fn = loss_er, method = "L-BFGS-B",  
                 N=N, y=y, x=x, tau=tau, d=d)
  beta   = Opt$par
  if(d==1) res = y -  (x * beta)
  else     res = y -  (x %*% beta)   
  return(list(alpha=0, beta=beta, lambda=0, res=res))
}

#' @title optim expectile regression with fixed effects
#'
#' @description This function solves a expectile regression  with fixed effects
#' 
#' @param beta Numeric vector, initials values beta.
#' @param alpha Numeric vector, initials values alpha.
#' @param x Numeric matrix, covariates.
#' @param y Numeric vector, output.
#' @param z Numeric matrix, incidence matrix.
#' @param tau Numeric scalar, the percentile.
#' @param N Numeric integer, sample size.
#' @param d Numeric integer, X number of columns.
#' @param n Numeric integer, length of alpha.
#' 
#' @return parametric vector and residuals.
optim_erfe = function(beta,alpha,x,y,z,tau,N,d,n){
  M      = n
  theta  = c(beta,alpha)
  Opt    = stats::optim(par=theta, fn = loss_erfe, method = "L-BFGS-B",  
                 n=N, y=y, x=x, z=z, tau=tau, d=d, mm=n)
  beta   = Opt$par[1:d]
  alpha  = Opt$par[(d+1):(d+n)]
  if(d==1) res = y -  (z %*% alpha) -  (x * beta)
  else     res = y -  (z %*% alpha) -  (x %*% beta) 
  sdu    = stats::sd(res, na.rm = TRUE)
  sda    = stats::sd(alpha, na.rm = TRUE)
  lambda = check_lambda(sdu/sda, 1/M, M)
  return(list(alpha=alpha, beta=beta, lambda=lambda, res=res))
}

#' @title optim expectile regression with fixed effects and LASSO
#'
#' @description This function solves a expectile regression  with fixed effects and LASSO
#' 
#' @param beta Numeric vector, initials values beta.
#' @param alpha Numeric vector, initials values alpha.
#' @param x Numeric matrix, covariates.
#' @param y Numeric vector, output.
#' @param z Numeric matrix, incidence matrix.
#' @param tau Numeric scalar, the percentile.
#' @param N Numeric integer, sample size.
#' @param d Numeric integer, X number of columns.
#' @param n Numeric integer, length of alpha.
#' 
#' @return parametric vector and residuals.
optim_erlasso = function(beta,alpha,x,y,z,tau,N,d,n){
  M      = n
  theta  = c(beta,alpha)
  Opt    = stats::optim(par=theta, fn = loss_erlasso, method = "L-BFGS-B",   
                 n=N, y=y, x=x, z=z, tau=tau, d=d, mm=n,lambda=1)
  beta   = Opt$par[1:d]
  alpha  = Opt$par[(d+1):(d+n)]
  if(d==1) res = y -  (z %*% alpha) -  (x * beta)
  else     res = y -  (z %*% alpha) -  (x %*% beta)   
  sdu    = stats::sd(res, na.rm = TRUE)
  sda    = stats::sd(alpha, na.rm = TRUE)
  lambda = check_lambda(sdu/sda, 1/M, M)
  return(list(alpha=alpha, beta=beta, lambda=lambda, res=res))
}

#' @title Kernel density
#'
#' @param x Numeric vector.
#' 
#' @return y vector, kernel density estimation.
#' 
#' @examples 
#' x = rnorm(10)
#' f_den(x)
#' 
f_den = function(x){
  inf = 1e-08
  n  = length(x)
  dh = stats::density(x, kernel = "epanechnikov")
  dx = dh$x
  dy = dh$y
  M  = length(dx)
  y  = rep(max(min(dy),inf),n)
  for(i in 1:n){
    for(j in 2:(M-2)){
      if(x[i]>=dx[j-1] && x[i]<=dx[j+1]) y[i] = dy[j] 
    }
  }
  return(y)  
}

#' @title Covariance
#'
#' @description Estimate Covariance matrix
#'
#' @param n length of alpha.
#' @param N sample size.
#' @param d length of beta.
#' @param Z Numeric matrix, incident matrix.
#' @param X Numeric matrix, covariates.
#' @param tau Numeric, identifies the percentile.
#' @param res Numeric vector, residuals.
#' @param penalty Numeric, 1 quantile regression, 2 quantile regression with fixed effects, 3 Lasso quantile regression with fixed effects
#' @param c Numeric, tuning
#' 
#' @return a list with two matrices: sig2_alpha (which is the matrix of covariance of estimated alpha) and sig2_beta (which is the matrix of covariance of estimated beta)
q_cov = function(n, N, d, Z, X, tau, res, penalty, c){
  omega = tau*(1-tau)
  inf = 1/(10^8)
  if(c==Inf){
    gres = as.vector(ifelse(res<0, 1-tau, tau)) 
    gsd  = as.vector(gres^2 * res^2) 
    G    = diag(gres)
    H    = diag(gsd)
    zgz = matrix(n * t(Z) %*% G %*% Z, nrow = n, ncol = n)
    zgx = matrix(sqrt(n) * t(Z) %*% G %*% X, nrow = n, ncol = d)
    xgz = matrix(sqrt(n) * t(X) %*% G %*% Z, nrow = d, ncol = n)
    xgx = matrix(t(X) %*% G %*% X, nrow = d, ncol = d) 
    Ga  = (1/N) * matrix(rbind(cbind(zgz,zgx), cbind(xgz,xgx)), nrow = (n+d), ncol = (n+d))
    D1inv = MASS::ginv(Ga)    
    zhz = matrix(n * t(Z) %*% H %*% Z, nrow = n, ncol = n)
    zhx = matrix(sqrt(n) * t(Z) %*% H %*% X, nrow = n, ncol = d)
    xhz = matrix(sqrt(n) * t(X) %*% H %*% Z, nrow = d, ncol = n)
    xhx = matrix(t(X) %*% H %*% X, nrow = d, ncol = d) 
    D0  = (1/N) * matrix(rbind(cbind(zhz,zhx), cbind(xhz,xhx)), nrow = (n+d), ncol = (n+d))
    if(penalty==1){
      sig2_alpha = 0
      if(d==1) sig2_beta = (1/N) * (N^2)*solve(xgx) %*% xhx %*% solve(xgx)
      if(d> 1) sig2_beta = diag((1/N) * (N^2)*solve(xgx) %*% xhx %*% solve(xgx)) 
    }else{
      Sigma = D1inv %*% D0 %*% D1inv
      sig2_alpha = diag(Sigma[(1:n),(1:n)])
      if(d==1) sig2_beta  = Sigma[(n+1),(n+1)]
      if(d> 1) sig2_beta  = diag(Sigma[(n+1):(n+d),(n+1):(n+d)])
    }
  }
  if(c<Inf && c>0){
    gres = as.vector(ifelse(res<0, 1-tau, tau)) 
    gsd  = as.vector(rho_mq(x=res,tau=tau,c=c))
    G    = diag(gres)
    H    = diag(gsd)
    zgz = matrix(n * t(Z) %*% G %*% Z, nrow = n, ncol = n)
    zgx = matrix(sqrt(n) * t(Z) %*% G %*% X, nrow = n, ncol = d)
    xgz = matrix(sqrt(n) * t(X) %*% G %*% Z, nrow = d, ncol = n)
    xgx = matrix(t(X) %*% G %*% X, nrow = d, ncol = d) 
    Ga  = (1/N) * matrix(rbind(cbind(zgz,zgx), cbind(xgz,xgx)), nrow = (n+d), ncol = (n+d))
    D1inv = MASS::ginv(Ga)    
    zhz = matrix(n * t(Z) %*% H %*% Z, nrow = n, ncol = n)
    zhx = matrix(sqrt(n) * t(Z) %*% H %*% X, nrow = n, ncol = d)
    xhz = matrix(sqrt(n) * t(X) %*% H %*% Z, nrow = d, ncol = n)
    xhx = matrix(t(X) %*% H %*% X, nrow = d, ncol = d) 
    D0  = (1/N) * matrix(rbind(cbind(zhz,zhx), cbind(xhz,xhx)), nrow = (n+d), ncol = (n+d))
    if(penalty==1){
      sig2_alpha = 0
      if(d==1) sig2_beta = (1/N) * (N^2)*solve(xgx) %*% xhx %*% solve(xgx)
      if(d> 1) sig2_beta = diag((1/N) * (N^2)*solve(xgx) %*% xhx %*% solve(xgx)) 
    }else{
      Sigma = D1inv %*% D0 %*% D1inv
      sig2_alpha = diag(Sigma[(1:n),(1:n)])
      if(d==1) sig2_beta  = Sigma[(n+1),(n+1)]
      if(d> 1) sig2_beta  = diag(Sigma[(n+1):(n+d),(n+1):(n+d)])
    }
  }  
  if(c==0){
    zz  = matrix(n * t(Z) %*% Z, nrow = n, ncol = n)
    zx  = matrix(sqrt(n) * t(Z) %*% X, nrow = n, ncol = d)
    xz  = matrix(sqrt(n) * t(X) %*% Z, nrow = d, ncol = n)
    xx  = matrix(t(X) %*% X, nrow = d, ncol = d) 
    D0  = drop(omega/N) * matrix(rbind(cbind(zz,zx), cbind(xz,xx)), nrow = (n+d), ncol = (n+d))
    fij = f_den(res) + inf + stats::runif(N,inf,2*inf)
    Phi = diag(fij)
    zpz = matrix(n * t(Z) %*% Phi %*% Z, nrow = n, ncol = n)
    zpx = matrix(sqrt(n) * t(Z) %*% Phi %*% X, nrow = n, ncol = d)
    xpz = matrix(sqrt(n) * t(X) %*% Phi %*% Z, nrow = d, ncol = n)
    xpx = matrix(t(X) %*% Phi %*% X, nrow = d, ncol = d) 
    D1  = (1/N) * matrix(rbind(cbind(zpz,zpx), cbind(xpz,xpx)), nrow = (n+d), ncol = (n+d))
    D1inv = MASS::ginv(D1)
    if(penalty==1){
      sig2_alpha = 0
      if(d==1) sig2_beta = (omega/N) * (N^2)*solve(xpx) %*% xx %*% solve(xpx)
      if(d> 1) sig2_beta = diag((omega/N) * (N^2)*solve(xpx) %*% xx %*% solve(xpx)) 
    }else{
      Sigma = D1inv %*% D0 %*% D1inv
      sig2_alpha = diag(Sigma[(1:n),(1:n)])
      if(d==1) sig2_beta  = Sigma[(n+1),(n+1)]
      if(d> 1) sig2_beta  = diag(Sigma[(n+1):(n+d),(n+1):(n+d)])
    }
  }
  return(list(sig2_alpha=sig2_alpha, sig2_beta=sig2_beta))
}

#' @title Tabular function
#' 
#' @param N sample size.
#' @param n length of alpha.
#' @param d length of beta.
#' @param theta Numeric vector.
#' @param sig2 Numeric vector.
#' @param kind Numeric, 1 means alpha, 2 means beta
#'
#' @return a list with a dataframe Core and a matrix Matx, both display the same information
f_tab = function(N, n, d, theta, sig2, kind){
  inf = 1e-08
  m = N/n
  len   = stats::qnorm(0.975)
  p     = length(theta)
  for(i in 1:p){
    if(is.na(sig2[i])) sig2[i] = inf
    if(sig2[i]<=0)     sig2[i] = inf
  }
  if(kind==1) SE = sqrt(sig2/(m))
  if(kind==2) SE = sqrt(sig2/(N))  
  infb  = theta - (len * SE) 
  supb  = theta + (len * SE)
  zval  = theta/SE
  pval  = 2 * stats::pnorm(abs(theta/SE), lower.tail = F)
  sig0  = sgf(as.vector(pval))
  Matx  = matrix(cbind(theta, SE, infb, supb, zval, pval), ncol=6, nrow=p)
  Core  = data.frame(cbind(theta, SE, infb, supb, zval, pval, sig0))
  colnames(Core) = c("Coef", "Std. Error", "Inf CI95%", "Sup CI95%", "z value", "Pr(>|z|)", "Sig")
  if(kind==1){
    if(p==1) rownames(Core)[1] = "alpha 1"
    if(p> 1) rownames(Core) = paste("alpha", 1:p)    
  } 
  if(kind==2){
    if(d==1) rownames(Core)[1] = "beta 1"
    if(d> 1) rownames(Core) = paste("beta", 1:p)
  }   
  return(list(Core=Core, Matx=Matx))
}

#' @title Identify significance
#'
#' @param x Numeric vector.
#' 
#' @return y vector Factor, symbol flag of significant p-values.
#' 
#' @examples 
#' n = 10
#' pvalue = rgamma(10,1,10)
#' sgf(pvalue) 
#' 
#' @return a vector of Factors, i.e., the symbols to help p-value interpretation
sgf = function(x){
  m = length(x)
  y = rep(" ", m)
  for(i in 1:m){
    if(is.na(x[i]))x[i]=0
    if(x[i]<0.001) y[i] = "***"
    if(x[i]>=0.001 && x[i]<0.01) y[i] = "**"
    if(x[i]>=0.01 && x[i]<0.05) y[i] = "*"
    if(x[i]>=0.05 && x[i]<0.1) y[i] = "."
  }
  return(y)
}

#' @title check lambda
#'
#' @param lambda Numeric, value of lambda.
#' @param infb Numeric, lower bound of lambda.
#' @param supb Numeric, upper bound of lambda.
#' 
#' @return lambda Numeric, valid value of lambda.
check_lambda = function(lambda, infb, supb){
  if(is.na(lambda) || is.null(lambda)){
    lambda = stats::runif(1,infb, supb)
  }else{
    if(lambda <= infb){
      lambda = infb;
    }
    if(lambda >= supb){
      lambda = supb;
    }    
  }
  return(lambda);  
}

#' @title choice model
#'
#' @param effect Factor, simple, fixed or lasso.
#' 
#' @return penalty Numeric, 1, 2 and 3.
choice_p = function(effect){
  penalty = 1
  if(effect == "fixed") penalty = 2
  if(effect == "lasso") penalty = 3
  return(penalty)  
}

#' @title Clean missings
#'
#' @param y Numeric vector, outcome.
#' @param x Numeric matrix, covariates
#' @param id Numeric vector, identifies the unit to which the observation belongs.
#' 
#' @returns list with the same objects y, x, id, but without missings. 
#' 
#' @examples
#' n = 10
#' m = 4
#' d = 3
#' N = n*m
#' L = N*d
#' x = matrix(rnorm(L), ncol=d, nrow=N)
#' subj = rep(1:n, each=m)
#' alpha = rnorm(n)
#' beta = rnorm(d)
#' eps = rnorm(N)
#' y = x %*% beta  + matrix(rep(alpha, each=m) + eps)
#' y = as.vector(y)
#' x[1,3] = NA
#' clean_data(y=y, x=x, id=subj)  
#'  
clean_data = function(y, x, id){
  n = max(id)
  N = length(id)
  d = dim(x)[2]
  ynew = NULL
  xnew = NULL
  nj   = NULL
  inew = NULL
  for(j in 1:n){
    n0 = length(ynew)
    for(i in 1:N){
      if(id[i] == j){
        if( any(c(is.null(y[i]), is.null(x[i,]), is.na(y[i]), is.na(x[i,])))){
        } else {
            ynew = c(ynew, y[i])
            xnew = matrix(rbind(xnew, x[i,]), ncol=d)
            inew = c(inew, id[i])
        }  
      }
    }
    n1 = length(ynew) - n0
    if(n1>0) nj = rbind(nj, n1)  
  }
  y = as.vector(ynew)
  x = as.matrix(xnew)
  nj = as.vector(nj)
  id = as.vector(inew)
  return(list(y=y,x=x,nj=nj,id=id))
}

#' @title Print an PQR
#'
#' @description Define the visible part of the object class PQR
#'
#' @param x An object of class "PQR"
#' @param ... further arguments passed to or from other methods.
#' 
#' @return None
print.PQR = function(x,...){
  if(x$penalty == 1){
    if(x$c==0){
      cat("\n Quantile regression \n \n")
      base::print(x$Tab_beta)
      invisible(x)    
    }
    if(x$c>0 && x$c<Inf){
      cat("\n M-Quantile regression \n \n")
      base::print(x$Tab_beta)
      invisible(x)    
    }
    if(x$c==Inf){
      cat("\n Expectile regression \n \n")
      base::print(x$Tab_beta)
      invisible(x)    
    }
  }  
  if(x$penalty == 2){
    if(x$c==0){
      cat("\n Quantile regression with fixed effects \n \n")
      base::print(x$Tab_beta)
      invisible(x)    
    }
    if(x$c>0 && x$c<Inf){
      cat("\n M-Quantile regression with fixed effects \n \n")
      base::print(x$Tab_beta)
      invisible(x)    
    }
    if(x$c==Inf){
      cat("\n Expectile regression with fixed effects \n \n")
      base::print(x$Tab_beta)
      invisible(x)    
    }
  }
  if(x$penalty == 3){
    if(x$c==0){
      cat("\n Penalized Quantile regression with fixed effects \n \n")
      base::print(x$Tab_beta)
      invisible(x)    
    }
    if(x$c>0 && x$c<Inf){
      cat("\n Penalized M-Quantile regression with fixed effects \n \n")
      base::print(x$Tab_beta)
      invisible(x)   
    }
    if(x$c==Inf){
      cat("\n Penalized Expectile regression with fixed effects \n \n")
      base::print(x$Tab_beta)
      invisible(x)    
    }
  }
}

#' @title Multiple penalized quantile regression
#'
#' @description Estimate penalized quantile regression for several taus
#'
#' @param y Numeric vector, outcome.
#' @param x Numeric matrix, covariates
#' @param subj Numeric vector, identifies the unit to which the observation belongs.
#' @param tau Numeric vector, identifies the percentiles.
#' @param effect Factor, "simple" simple regression, "fixed" regression with fixed effects, "lasso" penalized regression with fixed effects.
#' @param c  Numeric, 0 is quantile, Inf is expectile, any number between zero and infinite is M-quantile.
#'
#' @return Beta Numeric array, with three dimmensions: 1) tau, 2) coef., lower bound, upper bound, 3) exploratory variables.
#'
#' @examples
#' n = 10
#' m = 5
#' d = 4
#' N = n*m
#' L = N*d
#' x = matrix(rnorm(L), ncol=d, nrow=N)
#' subj = rep(1:n, each=m)
#' alpha = rnorm(n)
#' beta = rnorm(d)
#' eps = rnorm(N)
#' y = as.vector(x %*% beta + rep(alpha, each=m) + eps)
#'
#' Beta = mpqr(x,y,subj,tau=1:9/10, effect="fixed", c = 1.2)
#' Beta
#' 
#' @returns Beta array with dimension (ntau, 3, d), where Beta[i,1,k] is the i-th tau estimation
#' of beta_k, Beta[i,2,k] is the i-th tau lower bound 95\% confidence of beta_k, and Beta[i,3,k] 
#' is the i-th tau lower bound 95\% confidence of beta_k.   
#'
#' @export
mpqr = function(x,y,subj,tau=1:9/10, effect="simple", c=0){
  ntau = length(tau)
  d    = ifelse(is.null(dim(x)[2]), 1, dim(x)[2])
  Beta = array(dim = c(ntau, 3, d))
  for(i in 1:ntau){
    Est = pqr(x,y,subj,tau[i], effect, c)
    Beta[i,1,] = Est$Mat_beta[,1] 
    Beta[i,2,] = Est$Mat_beta[,3] 
    Beta[i,3,] = Est$Mat_beta[,4] 
  }
  return(Beta)
}

#' @title Plot multiple penalized quantile regression
#'
#' @description plot penalized quantile regression for several taus
#'
#' @param Beta Numeric array, with three dimmensions: 1) tau, 2) coef., lower bound, upper bound, 3) exploratory variables.
#' @param tau Numeric vector, identifies the percentiles.
#' @param D covariate's number.
#' @param col color.
#' @param lwd line width.
#' @param lty line type.
#' @param pch point character.
#' @param cex.axis cex axis length.
#' @param cex.lab cex axis length.
#' @param main title. 
#' @param shadow color of the Confidence Interval 95\%
#'
#' @return None
#' 
#' @examples
#' n = 10
#' m = 5
#' d = 4
#' N = n*m
#' L = N*d
#' x = matrix(rnorm(L), ncol=d, nrow=N)
#' subj = rep(1:n, each=m)
#' alpha = rnorm(n)
#' beta = rnorm(d)
#' eps = rnorm(N)
#' y = as.vector(x %*% beta + rep(alpha, each=m) + eps)
#'
#' Beta = mpqr(x,y,subj,tau=1:9/10, effect="lasso", c = Inf)
#' plot_taus(Beta,tau=1:9/10,D=1)
#'
#' @export
plot_taus = function(Beta, tau=1:9/10, D, col=2, lwd=1, lty=2, pch=16, cex.axis=1, cex.lab=1, main="", shadow="gray90"){
  ntau  = dim(Beta)[1]
  d     = dim(Beta)[3]
  Beta  = matrix(Beta[,,D], ncol = 3, nrow = ntau)
  Mbeta = max(Beta) 
  mbeta = min(Beta)
  Mtau  = max(tau)
  mtau  = min(tau)
  graphics::plot(c(mtau,Mtau),c(mbeta, Mbeta), xlab=expression(tau), ylab=expression(paste(beta,"(",tau,")")), main=main, type="n", cex.axis=cex.axis, cex.lab=cex.lab)
  graphics::polygon(c(tau,tau[ntau:1]), c(Beta[,2],Beta[ntau:1,3]), col=shadow, border = NA)
  graphics::lines(tau, Beta[,1], col=col, lty=lty, lwd=lwd)
  graphics::lines(tau, Beta[,1], col=col, type = "p", pch=pch)
  graphics::lines(tau, rep(0,ntau), lty=3, lwd=lwd)
}  

  
  
  
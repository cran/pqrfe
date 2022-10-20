// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// rho_koenker
arma::vec rho_koenker(arma::vec x, double tau);
RcppExport SEXP _pqrfe_rho_koenker(SEXP xSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(rho_koenker(x, tau));
    return rcpp_result_gen;
END_RCPP
}
// rho_mq
arma::vec rho_mq(arma::vec x, double tau, double c);
RcppExport SEXP _pqrfe_rho_mq(SEXP xSEXP, SEXP tauSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(rho_mq(x, tau, c));
    return rcpp_result_gen;
END_RCPP
}
// psi_mq
arma::vec psi_mq(arma::vec x, double tau, double c);
RcppExport SEXP _pqrfe_psi_mq(SEXP xSEXP, SEXP tauSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(psi_mq(x, tau, c));
    return rcpp_result_gen;
END_RCPP
}
// d_psi_mq
arma::vec d_psi_mq(arma::vec x, double tau, double c);
RcppExport SEXP _pqrfe_d_psi_mq(SEXP xSEXP, SEXP tauSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(d_psi_mq(x, tau, c));
    return rcpp_result_gen;
END_RCPP
}
// psi_als
arma::vec psi_als(arma::vec x, double tau);
RcppExport SEXP _pqrfe_psi_als(SEXP xSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(psi_als(x, tau));
    return rcpp_result_gen;
END_RCPP
}
// d_psi_als
arma::vec d_psi_als(arma::vec x, double tau);
RcppExport SEXP _pqrfe_d_psi_als(SEXP xSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(d_psi_als(x, tau));
    return rcpp_result_gen;
END_RCPP
}
// loss_qr
double loss_qr(arma::vec beta, arma::mat x, arma::vec y, double tau, int N, int d);
RcppExport SEXP _pqrfe_loss_qr(SEXP betaSEXP, SEXP xSEXP, SEXP ySEXP, SEXP tauSEXP, SEXP NSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(loss_qr(beta, x, y, tau, N, d));
    return rcpp_result_gen;
END_RCPP
}
// loss_qrfe
double loss_qrfe(arma::vec theta, arma::mat x, arma::vec y, arma::mat z, double tau, int n, int d, int mm);
RcppExport SEXP _pqrfe_loss_qrfe(SEXP thetaSEXP, SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP tauSEXP, SEXP nSEXP, SEXP dSEXP, SEXP mmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    rcpp_result_gen = Rcpp::wrap(loss_qrfe(theta, x, y, z, tau, n, d, mm));
    return rcpp_result_gen;
END_RCPP
}
// loss_qrlasso
double loss_qrlasso(arma::vec theta, arma::mat x, arma::vec y, arma::mat z, double tau, int n, int d, int mm, double lambda);
RcppExport SEXP _pqrfe_loss_qrlasso(SEXP thetaSEXP, SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP tauSEXP, SEXP nSEXP, SEXP dSEXP, SEXP mmSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(loss_qrlasso(theta, x, y, z, tau, n, d, mm, lambda));
    return rcpp_result_gen;
END_RCPP
}
// loss_mqr
double loss_mqr(arma::vec beta, arma::mat x, arma::vec y, double tau, int N, int d, double c);
RcppExport SEXP _pqrfe_loss_mqr(SEXP betaSEXP, SEXP xSEXP, SEXP ySEXP, SEXP tauSEXP, SEXP NSEXP, SEXP dSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(loss_mqr(beta, x, y, tau, N, d, c));
    return rcpp_result_gen;
END_RCPP
}
// loss_mqrfe
double loss_mqrfe(arma::vec theta, arma::mat x, arma::vec y, arma::mat z, double tau, int n, int d, int mm, double c);
RcppExport SEXP _pqrfe_loss_mqrfe(SEXP thetaSEXP, SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP tauSEXP, SEXP nSEXP, SEXP dSEXP, SEXP mmSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(loss_mqrfe(theta, x, y, z, tau, n, d, mm, c));
    return rcpp_result_gen;
END_RCPP
}
// loss_mqrlasso
double loss_mqrlasso(arma::vec theta, arma::mat x, arma::vec y, arma::mat z, double tau, int n, int d, int mm, double c, double lambda);
RcppExport SEXP _pqrfe_loss_mqrlasso(SEXP thetaSEXP, SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP tauSEXP, SEXP nSEXP, SEXP dSEXP, SEXP mmSEXP, SEXP cSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(loss_mqrlasso(theta, x, y, z, tau, n, d, mm, c, lambda));
    return rcpp_result_gen;
END_RCPP
}
// loss_er
double loss_er(arma::vec beta, arma::mat x, arma::vec y, double tau, int N, int d);
RcppExport SEXP _pqrfe_loss_er(SEXP betaSEXP, SEXP xSEXP, SEXP ySEXP, SEXP tauSEXP, SEXP NSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(loss_er(beta, x, y, tau, N, d));
    return rcpp_result_gen;
END_RCPP
}
// loss_erfe
double loss_erfe(arma::vec theta, arma::mat x, arma::vec y, arma::mat z, double tau, int n, int d, int mm);
RcppExport SEXP _pqrfe_loss_erfe(SEXP thetaSEXP, SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP tauSEXP, SEXP nSEXP, SEXP dSEXP, SEXP mmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    rcpp_result_gen = Rcpp::wrap(loss_erfe(theta, x, y, z, tau, n, d, mm));
    return rcpp_result_gen;
END_RCPP
}
// loss_erlasso
double loss_erlasso(arma::vec theta, arma::mat x, arma::vec y, arma::mat z, double tau, int n, int d, int mm, double lambda);
RcppExport SEXP _pqrfe_loss_erlasso(SEXP thetaSEXP, SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP tauSEXP, SEXP nSEXP, SEXP dSEXP, SEXP mmSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(loss_erlasso(theta, x, y, z, tau, n, d, mm, lambda));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pqrfe_rho_koenker", (DL_FUNC) &_pqrfe_rho_koenker, 2},
    {"_pqrfe_rho_mq", (DL_FUNC) &_pqrfe_rho_mq, 3},
    {"_pqrfe_psi_mq", (DL_FUNC) &_pqrfe_psi_mq, 3},
    {"_pqrfe_d_psi_mq", (DL_FUNC) &_pqrfe_d_psi_mq, 3},
    {"_pqrfe_psi_als", (DL_FUNC) &_pqrfe_psi_als, 2},
    {"_pqrfe_d_psi_als", (DL_FUNC) &_pqrfe_d_psi_als, 2},
    {"_pqrfe_loss_qr", (DL_FUNC) &_pqrfe_loss_qr, 6},
    {"_pqrfe_loss_qrfe", (DL_FUNC) &_pqrfe_loss_qrfe, 8},
    {"_pqrfe_loss_qrlasso", (DL_FUNC) &_pqrfe_loss_qrlasso, 9},
    {"_pqrfe_loss_mqr", (DL_FUNC) &_pqrfe_loss_mqr, 7},
    {"_pqrfe_loss_mqrfe", (DL_FUNC) &_pqrfe_loss_mqrfe, 9},
    {"_pqrfe_loss_mqrlasso", (DL_FUNC) &_pqrfe_loss_mqrlasso, 10},
    {"_pqrfe_loss_er", (DL_FUNC) &_pqrfe_loss_er, 6},
    {"_pqrfe_loss_erfe", (DL_FUNC) &_pqrfe_loss_erfe, 8},
    {"_pqrfe_loss_erlasso", (DL_FUNC) &_pqrfe_loss_erlasso, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_pqrfe(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
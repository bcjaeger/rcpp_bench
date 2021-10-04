
#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec subvec_mult_v1(arma::vec& u,
                         arma::vec& v,
                         arma::uword n){

 arma::vec out = u % v;
 return(out);

}

// [[Rcpp::export]]
arma::vec subvec_mult_v2(arma::vec& u,
                         arma::vec& v,
                         arma::uword n){

 u = u % v;
 return(u);

}

// [[Rcpp::export]]
arma::vec subvec_mult_v3(arma::vec& u,
                         arma::vec& v,
                         arma::uword n){

 u.subvec(0, n-1) = u % v;
 return(u);

}

// [[Rcpp::export]]
arma::vec subvec_mult_v4(arma::vec& u,
                         arma::vec& v,
                         arma::uword n){

 u(arma::span(0, n-1)) = u % v;
 return(u);

}

// [[Rcpp::export]]
arma::vec subvec_mult_v5(arma::vec& u,
                         arma::vec& v,
                         arma::uword n){

 for(arma::uword i = 0; i < n; i++) u[i] *= v[i];
 return(u);

}

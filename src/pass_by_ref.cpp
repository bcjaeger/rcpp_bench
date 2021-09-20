#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::mat x_node;

arma::mat* x_ptr;

// [[Rcpp::export]]
void global_bind(arma::mat& x){

 x_node = x;

}

// [[Rcpp::export]]
void advanced_constructor(arma::mat& x){

 arma::mat xx(x.begin(), x.n_rows, x.n_cols, false);

}


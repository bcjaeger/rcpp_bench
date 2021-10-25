#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::umat A;
arma::imat B;
arma::vec temp;

arma::uword i;

// [[Rcpp::export]]
void foo(arma::imat& x,
         arma::mat& y){

   A = arma::conv_to<arma::umat>::from(x);

   for(i = 0; i < y.n_cols; i++){
      temp = y.col(A(0,i));
   }

}

// [[Rcpp::export]]
void bar(IntegerMatrix& x,
         arma::mat& y){

   A = arma::conv_to<arma::umat>::from(
      arma::imat(x.begin(), x.nrow(), x.ncol(), false)
   );

   for(i = 0; i < y.n_cols; i++){
      temp = y.col(A(0,i));
   }

}

// [[Rcpp::export]]
void baz(IntegerMatrix& x,
         arma::mat& y){

   B = arma::imat(x.begin(), x.nrow(), x.ncol(), false);

   for(i = 0; i < y.n_cols; i++){
      temp = y.col(B(0,i));
   }

}

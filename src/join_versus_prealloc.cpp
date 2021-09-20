
#include<RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat arma_join(uword n_row, uword n_col){

 mat big;

 for(uword i = 0; i < n_row; i++){

  rowvec a(n_col);
  big = join_vert(big, a);

 }

 return(big);

}

// [[Rcpp::export]]
mat arma_prealloc(uword n_row, uword n_col, uword max_row){

 mat big(n_row, n_col);

 return(big.rows(span(0,max_row)));

}

// [[Rcpp::export]]
List arma_join_list(uword n_row, uword n_col){

 List big(n_row);

 for(uword i = 0; i < n_row; i++){

  big(i) = rowvec(n_col);

 }

 return(big);

}

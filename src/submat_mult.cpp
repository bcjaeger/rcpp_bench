#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec iter_2loop(arma::mat& x,
                     arma::uvec& x_rows,
                     arma::uvec& x_cols,
                     arma::vec& beta){

  arma::vec out (x_rows.size());
  arma::uvec::iterator iit;
  arma::uvec::iterator jit;
  arma::uword i = 0;
  arma::uword j = 0;

  for(iit = x_rows.begin(); iit < x_rows.end(); ++iit){
    j=0;
    for(jit = x_cols.begin(); jit < x_cols.end(); ++jit){
      out[i] += x.at(*iit, *jit) * beta[j];
      j++;
    }
    i++;
  }

  return(out);

}

// [[Rcpp::export]]
arma::vec submat(arma::mat& x,
                 arma::uvec& x_rows,
                 arma::uvec& x_cols,
                 arma::vec& beta){

  arma::mat x_sub = x(x_rows, x_cols);
  arma::vec out = x_sub * beta;
  return(out);

}

// [[Rcpp::export]]
arma::vec submat2(arma::mat& x,
                 arma::uvec& x_rows,
                 arma::uvec& x_cols,
                 arma::vec& beta){

  arma::vec out = x(x_rows, x_cols) * beta;
  return(out);

}

// [[Rcpp::export]]
arma::vec submat_cheat(arma::mat& x,
                       arma::vec& beta){

  arma::vec out = x * beta;
  return(out);

}



// [[Rcpp::export]]
arma::vec iter_bigout(arma::mat& x,
                      arma::uvec& x_rows,
                      arma::uvec& x_cols,
                      arma::vec& beta){

  arma::vec out (x.n_rows);
  arma::uvec::iterator iit;
  arma::uvec::iterator jit;
  arma::uword j = 0;

  for(iit = x_rows.begin(); iit < x_rows.end(); ++iit){
    j=0;
    for(jit = x_cols.begin(); jit < x_cols.end(); ++jit){
      out[*iit] += x.at(*iit, *jit) * beta[j];
      j++;
    }
  }

  return(out);

}

// [[Rcpp::export]]
arma::vec iter_bigbeta(arma::mat& x,
                       arma::uvec& x_rows,
                       arma::uvec& x_cols,
                       arma::vec& beta){

  arma::vec out (x_rows.size());
  arma::uvec::iterator iit;
  arma::uword i;

  arma::vec beta_big(x.n_cols);

  i = 0;

  for(iit = x_cols.begin(); iit < x_cols.end(); ++iit){
    beta_big[*iit] = beta[i];
    i++;
  }

  i=0;

  for(iit = x_rows.begin(); iit < x_rows.end(); ++iit){

    out[i] = arma::dot(x.row(*iit), beta_big);
    i++;

  }

  return(out);

}
